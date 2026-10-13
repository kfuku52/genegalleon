"""Filter modeled events using exact Pfam hits in bilateral supported gene pairs.

Uses saved query RPS-BLAST records, including explicit no-hit rows, rather than
best-hit protein annotations. No search, domain architecture inference or
neighbor-gene substitution is performed.
"""

import csv
import hashlib
import io
import math
import re
from collections import defaultdict
from contextlib import nullcontext
from functools import lru_cache
from itertools import product

from focus_hgt_gene_trees import background_supported, number, validate_link_identity
from gene_family_output_store import GeneFamilyOutputStore, read_only_observation

EVENT_FIELDS = ["pfam_filter_status", "pfam_filter_reason", "pfam_compared_pair_count",
                "pfam_passing_pair_count", "pfam_shared_pair_count", "pfam_both_no_hit_pair_count",
                "pfam_shared_accessions", "pfam_min_shared_query_coverage", "pfam_attention_flags",
                "pfam_best_pair_donor_gene", "pfam_best_pair_recipient_gene",
                "pfam_best_pair_donor_query_coverage", "pfam_best_pair_recipient_query_coverage",
                "pfam_filter_enabled", "length_ratio_filter_enabled", "minimum_protein_length_ratio",
                "pfam_coverage_passing_pair_count", "length_ratio_passing_pair_count",
                "pfam_best_pair_protein_length_ratio"]
PAIR_FIELDS = ["event_id", "orthogroup", "branch_id", "node_name", "event_index", "generax_transfer",
               "generax_donor_node", "generax_recipient_node", "donor_gene_id", "recipient_gene_id",
               "donor_annotation_status", "recipient_annotation_status", "donor_pfam_accessions",
               "recipient_pfam_accessions", "shared_pfam_accessions", "pair_status", "passes_pfam_filter",
               "donor_query_length_aa", "recipient_query_length_aa", "donor_shared_pfam_covered_aa",
               "recipient_shared_pfam_covered_aa", "donor_shared_pfam_query_coverage",
               "recipient_shared_pfam_query_coverage", "min_shared_pfam_query_coverage",
               "coverage_status", "pair_attention_flags", "protein_length_ratio", "min_protein_length_ratio",
               "length_ratio_status", "passes_length_ratio_filter", "passes_shared_pfam_filter", "passes_pair_filter"]
GENE_FIELDS = ["event_id", "orthogroup", "side", "gene_id", "query_length_aa",
               "annotation_status", "pfam_accessions", "source_rpsblast", "source_sha256", "gene_attention_flags"]

# Review flags only: these model names never constitute a domain blacklist.
GENERIC_DOMAIN_NAME = re.compile(r"^(?:Ank(?:_|$)|WD(?:40|_|$)|zf(?:_|$)|RRM(?:_|$)|PUF(?:_|$)|Homeobox(?:_|$)|PHD(?:_|$)|RING(?:_|$))", re.IGNORECASE)
SHORT_QUERY_AA = 100


def validate_shared_pfam_coverage(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError("Shared Pfam query coverage must be a finite fraction from 0 to 1") from None
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("Shared Pfam query coverage must be a finite fraction from 0 to 1")
    return value


def validate_length_ratio(value):
    try:
        return validate_shared_pfam_coverage(value)
    except ValueError:
        raise ValueError("Protein length ratio must be a finite fraction from 0 to 1") from None


def covered_aa(record, accessions):
    """Union of 1-based inclusive saved query coordinates, in amino acids."""
    merged = []
    intervals = [interval for accession in accessions for interval in record['intervals'].get(accession, ())]
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    return sum(end - start + 1 for start, end in merged)


def gene_flags(record):
    return 'short_query_protein_lt100aa' if record['length'] != '' and record['length'] < SHORT_QUERY_AA else ''


def query_record(rows):
    if not rows:
        return dict(status="annotation_record_unavailable", length="", pfam=set(), intervals={}, names={})
    lengths = {number(row["qlen"]) for row in rows}
    if len(lengths) != 1 or None in lengths:
        raise ValueError("Missing or conflicting Pfam query lengths")
    length = lengths.pop()
    if not length.is_integer() or length <= 0:
        raise ValueError("Invalid Pfam protein query length")
    accessions, no_hits, intervals, names = set(), 0, defaultdict(list), defaultdict(set)
    for row in rows:
        if not row["sacc"]:
            if any(value for field, value in row.items() if field not in {"qacc", "qlen"}):
                raise ValueError("Malformed Pfam no-hit record")
            no_hits += 1
            continue
        # GeneGalleon's Pfam_LE database titles contain pfamNNNNN, name, description.
        match = re.match(r"pfam(\d{5})\s*,", row["stitle"], re.IGNORECASE)
        if match is None:
            raise ValueError("Unmapped Pfam model: " + row["stitle"])
        evalue = number(row["evalue"])
        qstart, qend = number(row["qstart"]), number(row["qend"])
        if (evalue is None or evalue < 0 or qstart is None or qend is None
                or not qstart.is_integer() or not qend.is_integer() or not 1 <= qstart <= qend <= length):
            raise ValueError("Invalid saved Pfam hit coordinates or E-value")
        accession = "PF" + match[1]
        accessions.add(accession)
        intervals[accession].append((int(qstart), int(qend)))
        names[accession].add(row['stitle'].split(',')[1].strip())
    if no_hits and (accessions or no_hits != 1):
        raise ValueError("Conflicting or duplicate Pfam no-hit records")
    return dict(status="pfam_hit_detected" if accessions else "searched_no_pfam_hit",
                length=int(length), pfam=accessions, intervals=intervals, names=names)


def binary_option(value, name):
    """Normalize explicit API flags without treating the string '0' as true."""
    if isinstance(value, str):
        value = value.strip().lower()
        if value in {'0', 'false'}:
            return False
        if value in {'1', 'true'}:
            return True
    elif isinstance(value, (bool, int)) and value in (0, 1):
        return bool(value)
    raise ValueError(name + ' must be a boolean or an explicit 0/1 flag')


def filter_events(events, links, family_root, allow_both_no_pfam=False, min_shared_pfam_coverage=0.5,
                  require_length_ratio=True, min_length_ratio=0.5, require_shared_pfam=True):
    """Keep an event if ANY exact, retained, scaffold-passing pair qualifies.

Missing query/search records never qualify for the optional bilateral no-hit
exception. All eligible context genes remain linked to surviving events.
Coverage is the union of saved shared-domain query intervals / protein length,
required on both genes of the same pair. Explicit bilateral no-hit opt-in is an
exception with unmeasured coverage, never a fabricated coverage of zero or one.
Protein length ratio is min(lengths)/max(lengths), from the saved protein-query
records of the SAME pair, inclusive. Missing lengths never satisfy an enabled
length rule, including the no-domain exception. Pfam and length rules can be
disabled independently. The legacy passes_pfam_filter column is the combined
decision; passes_shared_pfam_filter records the domain-only decision.
"""
    minimum = validate_shared_pfam_coverage(min_shared_pfam_coverage)
    minimum_length = validate_length_ratio(min_length_ratio)
    allow_both_no_pfam = binary_option(allow_both_no_pfam, 'allow_both_no_pfam')
    require_length_ratio = binary_option(require_length_ratio, 'require_length_ratio')
    require_shared_pfam = binary_option(require_shared_pfam, 'require_shared_pfam')
    if any(set(EVENT_FIELDS) & set(row) for row in events):
        raise ValueError("Reserved Pfam filter columns already exist in event input")
    ids = {row["event_id"]: row for row in events}
    if len(ids) != len(events) or "" in ids:
        raise ValueError("Duplicate or empty Pfam filtering event ID")
    linked, identities = defaultdict(list), set()
    family_genes = defaultdict(set)
    for link in links:
        if link["event_id"] not in ids:
            continue
        event = ids[link["event_id"]]
        identity = link["event_id"], link["side"], link["gene_id"]
        if identity in identities or link["side"] not in {"donor", "recipient"}:
            raise ValueError("Duplicate or invalid Pfam event-gene link")
        identities.add(identity)
        validate_link_identity(event, link)
        if link.get("lineage_status", "retained") == "retained" and background_supported(link):
            linked[link["event_id"], link["side"]].append(link)
            family_genes[event["orthogroup"]].add(link["gene_id"])

    records, sources = {}, {}
    with read_only_observation():
        store = GeneFamilyOutputStore(family_root) if family_root else None
        with store.read_snapshot() if store else nullcontext():
            for family in sorted({row["orthogroup"] for row in events}):
                if not re.fullmatch(r"[A-Za-z0-9_.-]+", family) or family in {".", ".."}:
                    raise ValueError("Unsafe orthogroup identifier")
                # Index qualifying genes once rather than scanning all event
                # links again for every family. Event-specific pairs stay separate.
                genes = family_genes[family]
                if not genes:
                    # The scaffold criterion already failed. Unused domain
                    # inputs cannot affect this event's decision or its audit.
                    continue
                raw, logical = None, "rpsblast/" + family + "_rpsblast.tsv"
                if store:
                    # Both existing public filename conventions are supported.
                    for name in (family + "_rpsblast.tsv", family + ".rpsblast.tsv"):
                        try:
                            with store.open_binary("rpsblast", name) as handle:
                                raw = handle.read()
                            logical = "rpsblast/" + name
                            break
                        except FileNotFoundError:
                            continue
                grouped = defaultdict(list)
                if raw is not None:
                    sources[logical] = hashlib.sha256(raw).hexdigest()
                    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
                    fields = reader.fieldnames or []
                    required = {"qacc", "sacc", "qlen", "stitle", "qstart", "qend", "evalue"}
                    if len(fields) != len(set(fields)) or not required <= set(fields):
                        raise ValueError("Malformed saved Pfam RPS-BLAST columns: " + logical)
                    for row in reader:
                        if None in row or any(value is None for value in row.values()) or not row["qacc"]:
                            raise ValueError("Malformed saved Pfam RPS-BLAST row: " + logical)
                        if row["qacc"] in genes:
                            grouped[row["qacc"]].append(row)
                for gene in genes:
                    record = dict(query_record(grouped[gene]), source=logical, sha256=sources.get(logical, ""))
                    record['accession_text'] = '; '.join(sorted(record['pfam']))
                    record['short_flag'] = gene_flags(record)
                    record['generic_accessions'] = frozenset(accession for accession, names in record['names'].items()
                                                            if all(GENERIC_DOMAIN_NAME.match(name) for name in names))
                    records[family, gene] = record
            for logical, expected in sources.items():
                subdir, name = logical.split("/", 1)
                with store.open_binary(subdir, name) as handle:
                    if hashlib.sha256(handle.read()).hexdigest() != expected:
                        raise ValueError("Pfam filtering input changed during generation")

    # Interval unions also recur across different partner genes. Bound that
    # smaller cache separately so unique-pair cohorts benefit without retaining
    # a template for every combination.
    @lru_cache(maxsize=32768)
    def query_coverage(family, gene, shared):
        return covered_aa(records[family, gene], shared)

    signatures, cacheable_events = {}, set()
    for event in events:
        eid = event['event_id']
        signature = (event['orthogroup'],
                     tuple(sorted(link['gene_id'] for link in linked[eid, 'donor'])),
                     tuple(sorted(link['gene_id'] for link in linked[eid, 'recipient'])))
        if signature in signatures:
            cacheable_events.update((eid, signatures[signature]))
        else:
            signatures[signature] = eid
    del signatures

    # Repeated events may compare the same exact pair. Keep a bounded cache of
    # immutable evidence templates, never event-specific rows or decisions from
    # another family. Unique-pair workloads cannot grow this cache without bound.
    def pair_evidence(family, donor_gene, recipient_gene, result=None):
        dr, rr = records[family, donor_gene], records[family, recipient_gene]
        shared = frozenset(dr['pfam'] & rr['pfam'])
        both_empty = dr['status'] == rr['status'] == 'searched_no_pfam_hit'
        measured = bool(shared)
        dc, rc = (query_coverage(family, gene, shared) for gene in (donor_gene, recipient_gene)) if shared else (0, 0)
        dfrac, rfrac = (dc / dr['length'], rc / rr['length']) if measured else ('', '')
        coverage_passes = measured and dfrac >= minimum and rfrac >= minimum
        domain_passes = coverage_passes or (allow_both_no_pfam and both_empty)
        ratio = min(dr['length'], rr['length']) / max(dr['length'], rr['length']) if dr['length'] and rr['length'] else ''
        length_passes = ratio != '' and ratio >= minimum_length
        passed = (not require_shared_pfam or domain_passes) and (not require_length_ratio or length_passes)
        flags = []
        if measured and shared <= dr['generic_accessions'] and shared <= rr['generic_accessions']:
            flags.append('shared_repeat_or_generic_binding_domain_only')
        if dr['pfam'] and rr['pfam'] and dr['pfam'] != rr['pfam']:
            flags.append('pfam_domain_sets_differ_architecture_review')
        if ratio != '' and ratio < .5:
            flags.append('query_lengths_differ_over2fold')
        for side, record in (('donor', dr), ('recipient', rr)):
            if record['short_flag']:
                flags.append(side + '_' + record['short_flag'])
        status = ('shared_pfam_detected' if shared else
                  'annotation_record_unavailable' if 'annotation_record_unavailable' in (dr['status'], rr['status']) else
                  'both_searched_no_pfam_hit' if both_empty else
                  'one_searched_no_pfam_hit' if not dr['pfam'] or not rr['pfam'] else
                  'detected_pfam_sets_disjoint')
        result = {} if result is None else result
        result.update(donor_gene_id=donor_gene, recipient_gene_id=recipient_gene,
                    donor_annotation_status=dr['status'], recipient_annotation_status=rr['status'],
                    donor_pfam_accessions=dr['accession_text'], recipient_pfam_accessions=rr['accession_text'],
                    shared_pfam_accessions='; '.join(sorted(shared)), pair_status=status,
                    passes_pfam_filter=str(passed), passes_pair_filter=str(passed),
                    passes_shared_pfam_filter=str(domain_passes), passes_length_ratio_filter=str(length_passes),
                    donor_query_length_aa=dr['length'], recipient_query_length_aa=rr['length'],
                    donor_shared_pfam_covered_aa=dc if measured else '',
                    recipient_shared_pfam_covered_aa=rc if measured else '',
                    donor_shared_pfam_query_coverage=dfrac, recipient_shared_pfam_query_coverage=rfrac,
                    min_shared_pfam_query_coverage=minimum,
                    coverage_status='passed' if coverage_passes else 'below_minimum' if measured else
                    'explicit_bilateral_no_hit_exception' if allow_both_no_pfam and both_empty else 'unmeasured',
                    pair_attention_flags='; '.join(flags), protein_length_ratio=ratio,
                    min_protein_length_ratio=minimum_length,
                    length_ratio_status='passed' if length_passes else 'below_minimum' if ratio != '' else 'unmeasured')
        return result

    cached_pair_evidence = lru_cache(maxsize=8192)(pair_evidence)

    selected, event_audit, pairs, genes = [], [], [], []
    for event in events:
        sides = {side: linked[event["event_id"], side] for side in ("donor", "recipient")}
        comparisons = []
        for donor, recipient in product(sides["donor"], sides["recipient"]):
            row = {field: event.get(field, "") for field in PAIR_FIELDS[:8]}
            row.update(branch_id=event.get("branch_id", event.get("gene_tree_branch_id", "")),
                       node_name=event.get("node_name", event.get("gene_tree_node", "")))
            if event['event_id'] in cacheable_events:
                row.update(cached_pair_evidence(event['orthogroup'], donor['gene_id'], recipient['gene_id']))
            else:
                pair_evidence(event['orthogroup'], donor['gene_id'], recipient['gene_id'], row)
            pairs.append(row)
            comparisons.append(row)
        passing = [row for row in comparisons if row["passes_pfam_filter"] == "True"]
        shared_pairs = [row for row in comparisons if row['shared_pfam_accessions']]
        # A surviving event's representative must be one of its passing pairs.
        # No-hit exceptions have unmeasured coverage and rank below measured hits.
        candidates = passing or shared_pairs
        best = max(candidates, key=lambda row: min(row['donor_shared_pfam_query_coverage'],
                   row['recipient_shared_pfam_query_coverage']) if row['shared_pfam_accessions'] else -1) if candidates else {}
        reason = ("no_bilateral_scaffold_supported_gene_pair" if not comparisons else
                  "protein_length_ratio_below_minimum_or_unmeasured" if not passing and require_length_ratio and
                  any(not require_shared_pfam or row['passes_shared_pfam_filter'] == 'True' for row in comparisons) else
                  "shared_pfam_below_minimum_query_coverage" if require_shared_pfam and shared_pairs and not passing else
                  "no_qualifying_pfam_pair" if not passing else "")
        annotated = dict(event, pfam_filter_status="passed" if passing else "withheld", pfam_filter_reason=reason,
                         pfam_compared_pair_count=len(comparisons), pfam_passing_pair_count=len(passing),
                         pfam_shared_pair_count=sum(row["pair_status"] == "shared_pfam_detected" for row in comparisons),
                         pfam_both_no_hit_pair_count=sum(row["pair_status"] == "both_searched_no_pfam_hit" for row in comparisons),
                         pfam_shared_accessions="; ".join(sorted({p for row in passing
                             for p in row["shared_pfam_accessions"].split("; ") if p})),
                         pfam_min_shared_query_coverage=minimum,
                         pfam_attention_flags='; '.join(sorted({flag for row in passing
                             for flag in row['pair_attention_flags'].split('; ') if flag})),
                         pfam_best_pair_donor_gene=best.get('donor_gene_id', ''),
                         pfam_best_pair_recipient_gene=best.get('recipient_gene_id', ''),
                         pfam_best_pair_donor_query_coverage=best.get('donor_shared_pfam_query_coverage', ''),
                         pfam_best_pair_recipient_query_coverage=best.get('recipient_shared_pfam_query_coverage', ''),
                         pfam_filter_enabled=str(require_shared_pfam), length_ratio_filter_enabled=str(require_length_ratio),
                         minimum_protein_length_ratio=minimum_length,
                         pfam_coverage_passing_pair_count=sum(row['coverage_status'] == 'passed' for row in comparisons),
                         length_ratio_passing_pair_count=sum(row['passes_length_ratio_filter'] == 'True' for row in comparisons),
                         pfam_best_pair_protein_length_ratio=best.get('protein_length_ratio', ''))
        event_audit.append(annotated)
        if passing:
            selected.append(annotated)
        for side, side_links in sides.items():
            for link in side_links:
                record = records[event["orthogroup"], link["gene_id"]]
                genes.append(dict(event_id=event["event_id"], orthogroup=event["orthogroup"], side=side,
                                  gene_id=link["gene_id"], query_length_aa=record["length"],
                                  annotation_status=record["status"], pfam_accessions="; ".join(sorted(record["pfam"])),
                                  source_rpsblast=record["source"], source_sha256=record["sha256"],
                                  gene_attention_flags=gene_flags(record)))
    return selected, event_audit, pairs, genes, sources
