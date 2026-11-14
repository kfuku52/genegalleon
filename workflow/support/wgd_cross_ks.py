"""Relative-age diagnostics for exact, paired cross-synteny copy matches."""

import itertools
import random
from collections import Counter, defaultdict
from dataclasses import asdict
from math import ceil
from pathlib import Path
from statistics import median

import numpy as np
from Bio.Data import CodonTable
from cdskit.translate import translate_sequence_string

try:
    from fasta_sequence_store import fasta_records
    from pairwise_synteny import source_file
    from wgd_evidence import number, read_table
    from wgd_ks_outgroups import (
        OUTGROUP_FIELDS,
        SISTER_FIELDS,
        common_outgroup_requests,
        common_outgroup_rows,
        select_common_outgroups,
        select_sister_panels,
        sibling_panel_requests,
        sibling_panel_rows,
    )
    from wgd_quartet_ks import QUARTET_FIELDS, quartet_requests, quartet_results
    from wgd_rooted_copy_ks import (
        RATE_BALANCE_FIELDS,
        ROOTED_COPY_FIELDS,
        ROOTED_COPY_SUMMARY,
        rooted_copy_summary,
        rooted_rate_balance,
        run_rooted_copy_ks,
    )
    from wgd_terminal_ks import TERMINAL_FIELDS, TERMINAL_PAIR_FIELDS
except ImportError:
    from .fasta_sequence_store import fasta_records
    from .pairwise_synteny import source_file
    from .wgd_evidence import number, read_table
    from .wgd_ks_outgroups import (
        OUTGROUP_FIELDS,
        SISTER_FIELDS,
        common_outgroup_requests,
        common_outgroup_rows,
        select_common_outgroups,
        select_sister_panels,
        sibling_panel_requests,
        sibling_panel_rows,
    )
    from .wgd_quartet_ks import QUARTET_FIELDS, quartet_requests, quartet_results
    from .wgd_rooted_copy_ks import (
        RATE_BALANCE_FIELDS,
        ROOTED_COPY_FIELDS,
        ROOTED_COPY_SUMMARY,
        rooted_copy_summary,
        rooted_rate_balance,
        run_rooted_copy_ks,
    )
    from .wgd_terminal_ks import TERMINAL_FIELDS, TERMINAL_PAIR_FIELDS


PAIR_FIELDS = ("species_event_id", "focal", "family_id", "gene_a", "gene_b",
               "contig_a", "window_a", "contig_b", "window_b", "status", "dS")
FOCAL_FIELDS = ("species_event_id", "branch_id", "focal", "status", "placed_event_id",
                "num_eligible_pairs", "num_estimated_pairs", "num_window_components", "median_ds",
                "ci_lower", "ci_upper", "younger_ds", "younger_lower", "younger_upper",
                "older_ds", "older_lower", "older_upper", "older_boundary_status", "younger_boundary_status",
                "copy_interval_method")
MATCHED_FIELDS = ("species_event_id", "branch_id", "focal", "younger_sister", "younger_outgroup",
                  "older_sister", "older_outgroup", "num_complete_families", "younger_ds", "older_ds",
                  "older_minus_younger", "gap_lower", "gap_upper", "status")
MATCHED_TRIO_FIELDS = ("species_event_id", "branch_id", "focal", "sister", "outgroup", "num_complete_families",
                       "num_negative_family_contrasts", "matched_median_ds", "matched_lower", "matched_upper", "status")


def trio_pairs(trio):
    return [tuple(sorted(taxa)) for taxa in ((trio.focal, trio.sister), (trio.focal, trio.outgroup),
                                             (trio.sister, trio.outgroup))]


def complete_distances(keys, contrasts):
    families = sorted(set.intersection(*(set(contrasts.get(key, {})) for key in keys)))
    return np.array([[contrasts[key][family] for key in keys] for family in families], dtype=float).reshape(-1, len(keys))


def signed_correction(first, second, third):
    """Subtract from the larger distance first to preserve small residuals."""
    return (np.maximum(first, second) - third) + np.minimum(first, second)


def matched_boundary_gap(trio_young, trio_old, contrasts, draws, rng):
    """Estimate boundary ordering while keeping all six measurements within each family."""
    from nwkit.ksrate_model import corrected_ks

    values = complete_distances(trio_pairs(trio_young) + trio_pairs(trio_old), contrasts)
    n = len(values)
    if not n:
        return n, None, None, None, None, None, "no_complete_families"
    med = np.median(values, axis=0)
    young, old = corrected_ks(*med[:3]), corrected_ks(*med[3:])
    gap = old - young
    if min(young, old) < 0:
        return n, young, old, gap, None, None, "negative_correction"
    if n < 3 or draws < 40:
        return n, young, old, gap, None, None, "insufficient_family_bootstrap"
    gaps = np.empty(draws)
    for start in range(0, draws, 128):
        end = min(start + 128, draws)
        sample = rng.integers(n, size=(end - start, n))
        first, second, third, fourth, fifth, sixth = np.median(values[sample], axis=1).T
        gaps[start:end] = signed_correction(fourth, fifth, sixth) - signed_correction(first, second, third)
    lo, hi = map(float, np.quantile(gaps, (.025, .975)))
    return n, young, old, gap, lo, hi, "ordered" if lo > 0 else "reversed" if hi < 0 else "order_unresolved"


def ortholog_distances(contrasts):
    by_pair = defaultdict(dict)
    for row in contrasts:
        key = tuple(sorted((row["species_a"], row["species_b"])))
        by_pair[key][row["family_id"]] = float(row["ks"])
    return by_pair


def matched_trio_rows(trios, by_pair):
    """Describe the median signed correction across complete families."""
    from nwkit.ksrate_model import median_order_statistic_interval

    for trio in trios:
        first, second, third = complete_distances(trio_pairs(trio), by_pair).T
        # Keep negative contrasts to avoid upward truncation bias.
        corrected = signed_correction(first, second, third)
        center = lo = hi = None
        status = "no_complete_families"
        if len(corrected):
            center = float(np.median(corrected))
            lo, hi = median_order_statistic_interval(corrected, .05 / len(trios))
            if not np.isfinite(lo):
                lo = hi = None
            status = ("negative_median" if center < 0 else "interval_unavailable" if lo is None
                      else "interval_includes_negative" if lo < 0 else "estimated")
        values = (trio.species_event_id, trio.branch_id, trio.focal, trio.sister, trio.outgroup,
                  len(corrected), int(np.count_nonzero(corrected < 0)), center, lo, hi, status)
        yield dict(zip(MATCHED_TRIO_FIELDS, values, strict=True))


def matched_boundary_rows(tree, focal_rows, trios, by_pair, parameters):
    from nwkit.clade_index import CladeIndex

    index = CladeIndex(tree)
    parent_ids = {index.clade_id_for_node(node): index.clade_id_for_node(node.up) for node in tree.descendants()}
    by_trio = {(t.focal, t.species_event_id): t for t in trios}
    rng = np.random.default_rng(parameters["seed"])
    for row in focal_rows:
        event, focal = row["species_event_id"], row["focal"]
        younger = by_trio.get((focal, event))
        older = by_trio.get((focal, parent_ids[event]))
        trio_fields = (younger.sister, younger.outgroup) if younger else (None, None)
        trio_fields += (older.sister, older.outgroup) if older else (None, None)
        estimates = (matched_boundary_gap(younger, older, by_pair, parameters["ks_bootstrap"], rng)
                     if younger and older else (0, None, None, None, None, None, "missing_trio"))
        yield dict(zip(MATCHED_FIELDS, (event, row["branch_id"], focal, *trio_fields, *estimates), strict=True))


def cross_ks_cds_sources(path, workspace, species, codes, default_code):
    rows = read_table(path)
    if not rows or set(rows[0]) != {"species", "cds"}:
        raise ValueError("Cross-block Ks CDS manifest must contain species and cds columns and at least one row")
    result = {}
    for row in rows:
        name = row["species"]
        if name not in species or name in result:
            raise ValueError(f"Cross-block Ks CDS species must be unique species-tree tips: {name}")
        code = codes.get(name, default_code)
        if code not in CodonTable.generic_by_id:
            raise ValueError(f"Unknown NCBI genetic code for {name}: {code}")
        result[name] = {"cds": source_file(workspace, row["cds"], "species_cds", name,
                                           (".fa", ".fas", ".fasta", ".fna")), "genetic_code": code}
    return result


def selected_cds(source, species, genes):
    """Read only requested CDS; exact and unprefixed aliases must not compete."""
    aliases = {}
    for gene in genes:
        for alias in (gene, gene.removeprefix(species + "_")):
            if alias in aliases and aliases[alias] != gene:
                raise ValueError(f"Ambiguous requested CDS alias for {species}: {alias}")
            aliases[alias] = gene
    records, rejected, seen = {}, {}, set()
    for identifier, _, dna in fasta_records(Path(source["cds"]), identifiers=aliases):
        gene = aliases[identifier]
        if gene in seen:
            raise ValueError(f"Multiple CDS records match {gene}")
        seen.add(gene)
        dna = dna.upper().replace("U", "T")
        if not dna or len(dna) % 3 or set(dna) - set("ACGTRYSWKMBDHVN"):
            rejected[gene] = "invalid_cds"
            continue
        protein = translate_sequence_string(dna, source["genetic_code"], False)
        if protein.endswith("*"):
            dna, protein = dna[:-3], protein[:-1]
        if not protein or "*" in protein:
            rejected[gene] = "invalid_cds"
            continue
        records[gene] = dna, protein
    return records, rejected


def quality_status(pair, species, sequences, report):
    if report["status"] != "ok":
        return report["status"]
    retained, sites = number(report["retained_codons"]), number(report["S"])
    longest = max(len(sequences[name][gene][0]) // 3 for name, gene in zip(species, pair, strict=True))
    if retained is None or sites is None:
        raise ValueError("Estimable CDSKIT dS lacks retained codons or synonymous sites")
    return "ok" if retained >= max(60, ceil(7 * longest / 10)) and sites >= 20 else "low_codon_information"


def component_medians(hits):
    """Merge observations that share either genomic window or a gene family."""
    roots = list(range(len(hits)))

    def root(i):
        while roots[i] != i:
            roots[i] = roots[roots[i]]
            i = roots[i]
        return i

    seen = {}
    for i, (family, first, second, _) in enumerate(hits):
        for key in (("family", family), ("window", first), ("window", second)):
            roots[root(i)] = root(seen.setdefault(key, i))
    groups = defaultdict(list)
    for i, (_, _, _, value) in enumerate(hits):
        groups[root(i)].append(value)
    return sorted(median(group) for group in groups.values())


COPY_HOLDOUT_FIELDS = ("species_event_id", "branch_id", "focal", "num_estimated_pairs", "num_window_components",
                       "point_ds", "drop_component_min", "drop_component_max", "status")


def component_median_holdout_bounds(values):
    """Return extrema of deleting one item from sorted component medians."""
    if len(values) < 2:
        return None, None
    middle = len(values) // 2
    if len(values) % 2 == 0:
        return values[middle - 1], values[middle]
    return median(values[middle - 1:middle + 1]), median(values[middle:middle + 2])


def copy_median_interval(values):
    from nwkit.ksrate_model import median_order_statistic_interval

    if not values:
        return None
    lo, hi = median_order_statistic_interval(values, .05)
    return (lo, hi) if np.isfinite(lo) else None


def boundary(row):
    lo, hi = number(row.get("ci_lower")), number(row.get("ci_upper"))
    if (row.get("status") == "ok" and row.get("interval_status") == "ok"
            and row.get("monotone_from_younger_node") in {"yes", "not_comparable"}
            and lo is not None and hi is not None and lo <= hi):
        return lo, hi
    return None


def boundary_status(row):
    if not row:
        return "not_estimated"
    if boundary(row) is not None:
        return "ok"
    if row["status"] != "ok":
        return row["status"]
    if row["interval_status"] != "ok":
        return row["interval_status"]
    return "nonmonotone_or_invalid_interval"


class CrossKs:
    @staticmethod
    def age_pairs_from_trace(hits, age_species):
        pairs = set()
        for event, a, b, _, family, genes, _, positions in hits:
            for i, sp in enumerate((a, b)):
                if sp in age_species:
                    pair = tuple(sorted((genes[i], genes[i + 2])))
                    first, second = sorted((positions[i][1:], positions[i + 2][1:]))
                    pairs.add((event, sp, family, pair, first, second))
        return pairs

    def __init__(self, tree, sources, parameters, family_tables):
        from nwkit.clade_index import CladeIndex

        self.tree = tree
        self.sources = sources
        self.codes = {name: source["genetic_code"] for name, source in sources.items()}
        self.parameters = parameters
        self.index = CladeIndex(tree)
        self.nodes = {self.index.clade_id_for_node(node): node for node in tree.traverse()}
        self.leaves = {leaf.name: leaf for leaf in tree.leaves()}
        members, key, counts = family_tables
        self.single_copy = {
            row[key]: {name: row[name].strip() for name in sources if counts[row[key]][name] == 1}
            for row in members if all(value in {0, 1} for value in counts[row[key]].values())
        }
        self.archived = {}
        self.pairs = []
        self.focals = ()
        self.calibration_pairs = set()
        self.trios = []
        self.outgroup_selections = ()
        self.outgroup_rows = ()
        self.sister_panels = ()
        self.terminal_rows = ()
        self.copy_holdout_rows = ()

    def observe_family(self, family, leaves):
        if family in self.single_copy:
            present = {sp: gene for sp, gene in self.single_copy[family].items() if gene in leaves}
            if len(present) >= 2:
                self.archived[family] = present

    def prepare(self, pairs, candidates, evidence):
        from nwkit.ksrate import make_trios

        eligible = {row["species_event_id"] for row in candidates
                    if row["tree_support"] == "tree_score_supported"
                    and evidence.get(row["species_event_id"], {}).get("supported")}
        self.pairs = sorted(row for row in pairs if row[0] in eligible)
        self.focals = sorted({row[1] for row in self.pairs})
        target_boundaries = defaultdict(set)
        for event, focal, *_ in self.pairs:
            node = self.nodes[event]
            # The candidate's boundaries and its immediately neighboring branches.
            target_boundaries[focal].update(self.index.clade_id_for_node(n) for n in
                                            (node, node.up, node.up.up, *node.children)
                                            if n is not None and not n.is_leaf)
        trios, _ = make_trios(self.tree, self.focals)
        available = Counter(pair for present in self.archived.values()
                            for pair in itertools.combinations(sorted(present), 2))
        grouped = defaultdict(list)
        for trio in trios:
            taxa = (trio.focal, trio.sister, trio.outgroup)
            if (trio.species_event_id in target_boundaries[trio.focal]
                    and all(name in self.sources for name in taxa) and len({self.codes[name] for name in taxa}) == 1):
                pairs = tuple(itertools.combinations(sorted(taxa), 2))
                coverage = min(available[pair] for pair in pairs)
                if coverage:
                    grouped[trio.focal, trio.species_event_id].append((coverage, trio, pairs))
        for options in grouped.values():
            # Select before looking at dS so a convenient age cannot pick its own calibration.
            _, trio, pairs = min(options, key=lambda item: (-item[0], item[1].sister, item[1].outgroup))
            self.trios.append(trio)
            self.calibration_pairs.update(pairs)
        if limit := self.parameters.get("max_ks_outgroup_families", 0):
            self.outgroup_selections = list(select_common_outgroups(
                ((event, focal) for event, focal, *_ in self.pairs), self.trios, self.nodes, self.index,
                self.archived, self.codes, limit, self.parameters["seed"]))
            if self.parameters.get("ks_sister_panels"):
                self.sister_panels = list(select_sister_panels(
                    self.outgroup_selections, self.nodes, self.archived, self.codes, limit, self.parameters["seed"]))

    def requests(self):
        calibration = defaultdict(list)
        rng = random.Random(self.parameters["seed"])
        ordered = sorted(self.archived)
        rng.shuffle(ordered)
        limit = self.parameters["max_ks_families"]
        for family in ordered:
            present = self.archived[family]
            for a, b in itertools.combinations(sorted(present), 2):
                if (a, b) in self.calibration_pairs and len(calibration[a, b]) < limit:
                    calibration[a, b].append(((present[a], present[b]), (a, b, family, "ortholog")))
        requested = dict(itertools.chain.from_iterable(calibration.values()))
        by_focal = defaultdict(list)
        for branch, focal, _, pair, first, _ in self.pairs:
            by_focal[branch, focal].append((first, pair))
        for (_, focal), pairs in sorted(by_focal.items()):
            unique = sorted({pair for _, pair in pairs})
            limit = self.parameters["max_pairs"]
            if len(unique) > limit:
                # Equal-window round robin prevents large blocks from filling a cap.
                by_window = defaultdict(list)
                for first, pair in pairs:
                    by_window[first].append(pair)
                queues = [iter(sorted(set(by_window[key]))) for key in sorted(by_window)]
                sampled = {}
                for stripe in itertools.zip_longest(*queues):
                    sampled.update((pair, None) for pair in stripe if pair is not None)
                    if len(sampled) >= limit:
                        break
                unique = list(itertools.islice(sampled, limit))
            for pair in unique:
                requested.setdefault(pair, (focal, focal, "", "paralog"))
        return requested

    def run(self, output, cpus, candidates, estimate, write_table, traced_hits=(), *, terminal=None, rooted_copy=None,
            run_tool=None, method_label=None):
        requested = self.requests()
        outgroup_pairs = dict(common_outgroup_requests(self.outgroup_selections, self.archived))
        for pair, info in outgroup_pairs.items():
            requested.setdefault(pair, info)
        sister_pairs = dict(sibling_panel_requests(self.sister_panels, self.archived))
        for pair, info in sister_pairs.items():
            requested.setdefault(pair, info)
        quartet_limit = self.parameters.get("max_ks_quartets", 0)
        selected = ()
        if quartet_limit:
            selected, pairs = quartet_requests(traced_hits, candidates, self.codes, quartet_limit, self.parameters["seed"])
            for pair, info in pairs.items():
                requested.setdefault(pair, info)
        if terminal is not None:
            for pair, info in terminal.requested.items():
                requested.setdefault(pair, info)
        genes = defaultdict(set)
        for pair, info in requested.items():
            for name, gene in zip(info[:2], pair, strict=True):
                genes[name].add(gene)
        rooted_selection = list(rooted_copy.selection(selected, genes)) if rooted_copy is not None else ()
        sequences, rejected = {}, {}
        for name, targets in genes.items():
            sequences[name], rejected[name] = selected_cds(self.sources[name], name, targets)
        ready = {pair: info for pair, info in requested.items()
                 if all(gene in sequences[sp] for sp, gene in zip(info[:2], pair, strict=True))}
        print(f"Cross-block dS: {len(ready):,} unique CDS pairs ready ({len(requested):,} requested)", flush=True)
        reports = estimate(output, sequences, self.codes, ready, cpus, prefix="cross_ks_")
        results = {}
        for pair, info in requested.items():
            row = reports.get("|".join(pair))
            if row:
                status = quality_status(pair, info[:2], sequences, row)
                results[pair] = status, row["dS"] if status == "ok" else "NA"
            else:
                status = ("invalid_cds" if any(gene in rejected[sp] for sp, gene in zip(info[:2], pair, strict=True))
                          else "cds_unavailable")
                results[pair] = status, "NA"
        selection = [{"species_a": a, "species_b": b, "family_id": family, "gene_a": pair[0], "gene_b": pair[1],
                      "status": results[pair][0], "dS": results[pair][1]}
                     for pair, (a, b, family, role) in sorted(requested.items()) if role == "ortholog"]
        write_table(output / "cross_ks_orthologs.tsv", selection,
                    ("species_a", "species_b", "family_id", "gene_a", "gene_b", "status", "dS"))
        contrasts = [{**{field: row[field] for field in ("species_a", "species_b", "family_id")}, "ks": row["dS"]}
                     for row in selection if row["status"] == "ok"]
        write_table(output / "cross_ks_contrasts.tsv", contrasts, ("species_a", "species_b", "family_id", "ks"))
        bounds = self.divergence_bounds(contrasts, output, write_table)
        audit, hits, totals = [], defaultdict(list), defaultdict(int)
        for event, focal, family, pair, first, second in self.pairs:
            status, ds = results.get(pair, ("not_sampled", "NA"))
            key = event, focal
            totals[key] += 1
            if status == "ok":
                hits[key].append((family, first, second, float(ds)))
            audit.append(dict(zip(PAIR_FIELDS, (event, focal, family, *pair, *first, *second, status, ds), strict=True)))
        write_table(output / "cross_ks_paralogs.tsv", audit, PAIR_FIELDS)
        holdout = [] if self.parameters.get("ks_component_holdout") else None
        rows = self.summarize(candidates, bounds, hits, totals, holdout=holdout)
        write_table(output / "cross_ks_focals.tsv", rows, FOCAL_FIELDS)
        if holdout is not None:
            self.copy_holdout_rows = holdout
            write_table(output / "cross_ks_component_holdout.tsv", holdout, COPY_HOLDOUT_FIELDS)
        by_pair = ortholog_distances(contrasts)
        write_table(output / "cross_ks_matched_boundaries.tsv",
                    matched_boundary_rows(self.tree, rows, self.trios, by_pair, self.parameters), MATCHED_FIELDS)
        write_table(output / "cross_ks_matched_trios.tsv", matched_trio_rows(self.trios, by_pair), MATCHED_TRIO_FIELDS)
        if self.parameters.get("max_ks_outgroup_families"):
            sampled = defaultdict(dict)
            audit = []
            pair_fields = ("species_a", "species_b", "family_id", "gene_a", "gene_b", "status", "dS")
            for pair, (a, b, family, _) in sorted(outgroup_pairs.items()):
                status, value = results[pair]
                if status == "ok":
                    sampled[a, b][family] = float(value)
                audit.append(dict(zip(pair_fields, (a, b, family, *pair, status, value), strict=True)))
            write_table(output / "cross_ks_common_outgroup_pairs.tsv", audit, pair_fields)
            mapping = {r["species_event_id"]: r["branch_id"] for r in candidates}
            self.outgroup_rows = list(common_outgroup_rows(self.outgroup_selections, sampled, mapping))
            write_table(output / "cross_ks_common_outgroups.tsv", self.outgroup_rows, OUTGROUP_FIELDS)
            if self.parameters.get("ks_sister_panels"):
                audit = []
                for pair, (a, b, family, _) in sorted(sister_pairs.items()):
                    status, value = results[pair]
                    if status == "ok":
                        sampled[a, b][family] = float(value)
                    audit.append(dict(zip(pair_fields, (a, b, family, *pair, status, value), strict=True)))
                write_table(output / "cross_ks_sister_panel_pairs.tsv", audit, pair_fields)
                write_table(output / "cross_ks_sister_panels.tsv",
                            sibling_panel_rows(self.sister_panels, sampled, mapping), SISTER_FIELDS)
        if quartet_limit:
            write_table(output / "cross_ks_quartets.tsv", quartet_results(selected, results), QUARTET_FIELDS)
        if rooted_copy is not None:
            rooted_rows = run_rooted_copy_ks(rooted_selection, sequences, self.codes, cpus, output,
                                            self.parameters.get("cross_ks_method", "yn00"), method_label, run_tool,
                                            write_table, read_table)
            write_table(output / "cross_ks_rooted_copies.tsv", rooted_rows, ROOTED_COPY_FIELDS)
            write_table(output / "cross_ks_rooted_copy_summary.tsv", rooted_copy_summary(rooted_rows), ROOTED_COPY_SUMMARY)
            write_table(output / "cross_ks_rooted_rate_balance.tsv", rooted_rate_balance(rooted_rows), RATE_BALANCE_FIELDS)
        if terminal is not None:
            pairs, summaries = terminal.summarize(results)
            self.terminal_rows = summaries
            write_table(output / "terminal_ks_quartets.tsv", pairs, TERMINAL_PAIR_FIELDS)
            write_table(output / "terminal_ks_summary.tsv", summaries, TERMINAL_FIELDS)
        return rows

    def divergence_bounds(self, contrasts, output, write_table):
        from nwkit.ksrate_model import KsObservation, PairKs, simultaneous_median_corrections, trio_values

        trio_fields = ("focal", "sister", "outgroup", "species_event_id", "branch_id",
                       "num_families_focal_sister", "num_families_focal_outgroup", "num_families_sister_outgroup")
        fields = (*trio_fields, "status", "corrected_ks", "ci_lower", "ci_upper", "interval_status",
                  "monotone_from_younger_node", "ci_method", "num_simultaneous_pairs")
        pairs = PairKs([KsObservation(row["species_a"], row["species_b"], row["family_id"], float(row["ks"]))
                        for row in contrasts], self.sources) if contrasts else None
        values = trio_values(self.trios, pairs.estimates()) if pairs else [np.nan] * len(self.trios)
        intervals, n = simultaneous_median_corrections(pairs, self.trios, .95) if pairs else ({}, 0)
        bounds, previous = {}, {}
        for trio, value in zip(self.trios, values, strict=True):
            key = trio.focal, trio.species_event_id
            lo, hi = (number(v) for v in intervals.get(key, (None, None)))
            current = number(value)
            prior = previous.get(trio.focal)
            monotone = "not_comparable" if prior is None or current is None else ("yes" if current >= prior else "no")
            if current is not None:
                previous[trio.focal] = current
            bounds[key] = {**asdict(trio), "status": ("ok" if current is not None else
                           "negative_correction" if np.isfinite(value) else "missing_comparisons"),
                           "num_families_focal_sister": pairs.count(trio.focal, trio.sister) if pairs else 0,
                           "num_families_focal_outgroup": pairs.count(trio.focal, trio.outgroup) if pairs else 0,
                           "num_families_sister_outgroup": pairs.count(trio.sister, trio.outgroup) if pairs else 0,
                           "corrected_ks": current, "ci_lower": lo, "ci_upper": hi,
                           "interval_status": "ok" if lo is not None and hi is not None else "unavailable_interval",
                           "monotone_from_younger_node": monotone, "ci_method": "pair-median-bonferroni",
                           "num_simultaneous_pairs": n}
        write_table(output / "cross_ks_trios.tsv", ({k: row[k] for k in trio_fields} for row in bounds.values()), trio_fields)
        write_table(output / "cross_ks_boundaries.tsv", bounds.values(), fields)
        return bounds

    def summarize(self, candidates, bounds, hits, totals, *, holdout=None):
        by_event = defaultdict(list)
        minimum = max(6, self.parameters["min_cross_windows"])
        ids = {row["species_event_id"]: row["branch_id"] for row in candidates} if holdout is not None else None
        for (event, focal), total in sorted(totals.items()):
            records = hits[event, focal]
            values = component_medians(records)
            center = median(values) if values else None
            if holdout is not None:
                lo, hi = component_median_holdout_bounds(values)
                status = "ok" if lo is not None else "single_component" if values else "no_measured_pairs"
                result = (event, ids[event], focal, len(records), len(values), center, lo, hi, status)
                holdout.append(dict(zip(COPY_HOLDOUT_FIELDS, result, strict=True)))
            node = self.nodes[event]
            younger = bounds.get((focal, event), {})
            older = bounds.get((focal, self.index.clade_id_for_node(node.up)), {})
            ci = copy_median_interval(values) if len(values) >= minimum else None
            status, placed = "insufficient_window_components", ""
            if ci is not None:
                status, placed = self.place(focal, event, ci, bounds, younger, older)
            row = {"species_event_id": event, "focal": focal, "status": status, "placed_event_id": placed,
                   "num_eligible_pairs": total, "num_estimated_pairs": len(records),
                   "num_window_components": len(values), "median_ds": center if center is not None else "NA",
                   "ci_lower": ci[0] if ci else "NA", "ci_upper": ci[1] if ci else "NA",
                   "copy_interval_method": "component-median-order-statistic"}
            for name, data in (("younger", younger), ("older", older)):
                row.update({name + "_ds": data.get("corrected_ks", "NA"), name + "_lower": data.get("ci_lower", "NA"),
                            name + "_upper": data.get("ci_upper", "NA"), name + "_boundary_status": boundary_status(data)})
            by_event[event].append(row)
        for row in candidates:
            event = row["species_event_id"]
            focal_rows = by_event[event]
            supported = {r["focal"] for r in focal_rows if r["status"] == "supports_branch"}
            conflicts = {r["focal"] for r in focal_rows if r["status"] in {"older_than_branch", "younger_than_branch"}}
            if conflicts and supported:
                status = "discordant_focals"
            elif conflicts:
                status = "age_conflict"
            elif supported and all(supported.intersection(child.leaf_names()) for child in self.nodes[event].children):
                status = "supports_branch"
            elif supported:
                status = "partial_support"
            elif focal_rows:
                status = "inconclusive"
            else:
                status = "not_assessed"
            row.update(cross_ks_status=status, cross_ks_supporting_species=",".join(sorted(supported)),
                       cross_ks_conflicting_species=",".join(sorted(conflicts)), cross_ks_num_focals=len(focal_rows))
            for r in focal_rows:
                r["branch_id"] = row["branch_id"]
        return [r for candidate in candidates for r in by_event[candidate["species_event_id"]]]

    def place(self, focal, event, ci, bounds, younger, older):
        low, high = boundary(younger), boundary(older)
        if low and high and low[1] < high[0] and low[1] < ci[0] <= ci[1] < high[0]:
            return "supports_branch", event
        younger_conflict = low is not None and ci[1] < low[0]
        older_conflict = high is not None and ci[0] > high[1]
        status = ("younger_than_branch" if younger_conflict else "older_than_branch"
                  ) if younger_conflict != older_conflict else ""
        if not status:
            if low and high and low[1] >= high[0]:
                return "overlapping_divergence_boundaries", ""
            return ("overlaps_divergence_boundary" if low and high else "missing_divergence_boundary"), ""
        alternate = []
        node = self.leaves[focal]
        while node.up is not None:
            branch = self.index.clade_id_for_node(node)
            parent = self.index.clade_id_for_node(node.up)
            a = (0, 0) if node.is_leaf else boundary(bounds.get((focal, branch), {}))
            b = boundary(bounds.get((focal, parent), {}))
            if a and b and a[1] < ci[0] <= ci[1] < b[0]:
                alternate.append(branch)
            node = node.up
        return status, alternate[0] if len(alternate) == 1 else ""


def plot_cross_ks(rows, output):
    if not rows:
        return []
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    shown = rows[:30]
    fig, ax = plt.subplots(figsize=(11, max(3.5, 2 + len(shown) * 0.36)))
    for i, row in enumerate(reversed(shown)):
        for name, color, shift in (("younger", "#237a75", -0.13), ("older", "#9b5937", 0.13)):
            lo, hi = (number(row[name + "_" + field]) for field in ("lower", "upper"))
            if row[name + "_boundary_status"] == "ok" and lo is not None and hi is not None:
                ax.plot((lo, hi), (i + shift, i + shift), lw=5, alpha=.65, color=color, solid_capstyle="butt")
        center, lo, hi = (number(row[field]) for field in ("median_ds", "ci_lower", "ci_upper"))
        if center is not None:
            if lo is not None and hi is not None:
                ax.plot((lo, hi), (i, i), lw=2, color="#262e3d")
            ax.scatter(center, i, color="#262e3d", s=18, zorder=3)
        ax.axhline(i - 0.5, color="#e6e8e9", lw=.6)
        ax.text(1.02, i, row["status"].replace("_", " "), transform=ax.get_yaxis_transform(),
                va="center", fontsize=8)
    ax.set_yticks(range(len(shown)), [f"B{r['branch_id']}  {r['focal']}" for r in reversed(shown)], fontsize=8)
    ax.set(xlabel="Synonymous substitutions per synonymous site (focal scale)", xlim=(0, None))
    title = "WGD branch age checks" + (f" (first {len(shown)} of {len(rows)} focal estimates)" if len(rows) > len(shown) else "")
    ax.set_title(title, loc="left", pad=25)
    ax.legend(handles=[Line2D([], [], color=color, lw=width, label=label) for color, width, label in
                      (("#262e3d", 2, "Copy-pair median and order-statistic interval"),
                       ("#237a75", 4, "Younger divergence interval"), ("#9b5937", 4, "Older divergence interval"))],
              frameon=False, fontsize=8, ncols=3, loc="lower left", bbox_to_anchor=(-.02, 1))
    fig.subplots_adjust(left=.26, right=.77, bottom=max(.10, .34 / fig.get_figheight()), top=.86)
    paths = ["wgd_cross_ks." + suffix for suffix in ("pdf", "png")]
    for path in paths:
        fig.savefig(output / path, dpi=175, bbox_inches="tight")
    plt.close(fig)
    return paths
