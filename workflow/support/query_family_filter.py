#!/usr/bin/env python3
"""Classify query2family candidates using overlapping Swiss-Prot hits and InterPro Family IDs."""

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory

if __package__:
    from .fasta_sequence_store import open_text
else:
    from fasta_sequence_store import open_text

MAX_EVALUE = 1e-5
NEAR_BEST_FRACTION = 0.90
MIN_DISCOVERY_OVERLAP = 0.50
RESULT_FIELDS = ("query", "target", "bits", "evalue", "qstart", "qend")


def reachable(start, edges):
    remaining = [start]
    visited = {start}
    while remaining:
        related = edges.get(remaining.pop(), set()) - visited
        visited.update(related)
        remaining.extend(related)
    return visited


def read_family_context(entries, tree, mapping, query):
    with open_text(entries) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        types = {row["ENTRY_AC"]: row["ENTRY_TYPE"].lower() for row in rows}
    with open_text(mapping) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        expected = [row for row in rows if row["query"] == query]
    if len(expected) != 1:
        raise ValueError(f"{query} needs exactly one InterPro mapping row")
    target = expected[0]["interpro_family_id"]
    broader = {value.strip() for value in expected[0].get("broader_interpro_family_ids", "").split(";")
               if value.strip()}
    if target == "unmapped" and not broader:
        return None, set(), set()
    if types.get(target) != "family":
        raise ValueError(f"{query} needs an InterPro Family ID or an explicit unmapped entry")
    if any(types.get(family) != "family" for family in broader):
        raise ValueError(f"{query} has a broader ID that is not an InterPro Family")
    parents, children = defaultdict(set), defaultdict(set)
    stack = []
    with open_text(tree) as handle:
        for line in handle:
            match = re.match(r"^((?:--)*)((?:IPR)[0-9]+)::", line)
            if not match:
                raise ValueError(f"Invalid InterPro hierarchy row: {line.rstrip()}")
            indentation, current = match.groups()
            depth = len(indentation) // 2
            del stack[depth:]
            if depth != len(stack) or current in stack:
                raise ValueError(f"Invalid InterPro hierarchy for {current}")
            if stack:
                parents[current].add(stack[-1])
                children[stack[-1]].add(current)
            stack.append(current)
    families = {key for key, kind in types.items() if kind == "family"}
    compatible = reachable(target, children) & families
    if broader & compatible:
        raise ValueError(f"{query} lists a target or descendant as a broader Family")
    ancestors = reachable(target, parents) | broader
    return target, compatible, families - compatible - ancestors


def accession(identifier):
    fields = identifier.split("|")
    return fields[1] if len(fields) >= 3 and fields[0] == "sp" else identifier


def read_reference_membership(metadata, accessions, compatible, outside):
    membership = {}
    if not accessions:
        return membership
    with open_text(metadata) as handle:
        rows = csv.reader(handle, delimiter="\t")
        fields = next(rows, [])
        if not {"accession", "interpro_ids"} <= set(fields):
            raise ValueError("Swiss-Prot metadata requires accession and interpro_ids columns")
        accession_index, interpro_index = fields.index("accession"), fields.index("interpro_ids")
        for row in rows:
            ref = row[accession_index]
            if ref in accessions:
                ids = set(row[interpro_index].split(";"))
                matching, alternate = ids & compatible, ids & outside
                membership[ref] = (matching, alternate)
    return membership


def discovery_intervals(path, genes, min_coverage):
    intervals = defaultdict(set)
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            subject = row["sacc"]
            if subject in genes and float(row["qjointcov"]) >= min_coverage:
                starts = row["sstart"].split(";")
                ends = row["send"].split(";")
                for start, end in zip(starts, ends, strict=True):
                    lo, hi = sorted((int(start), int(end)))
                    if lo < 1:
                        raise ValueError(f"Invalid discovery alignment coordinates for {subject}")
                    intervals[subject].add((lo, hi))
    if genes - intervals.keys():
        raise ValueError("Missing discovery alignment coordinates for ranked candidates")
    return intervals


def overlapping_hits(hits, intervals):
    matches = defaultdict(list)
    with open_text(hits) as handle:
        for hit in csv.DictReader(handle, fieldnames=RESULT_FIELDS, delimiter="\t"):
            gene = hit["query"]
            score, evalue = float(hit["bits"]), float(hit["evalue"])
            start, end = sorted((int(hit["qstart"]), int(hit["qend"])))
            if gene not in intervals or start < 1 or not isfinite(score) or score <= 0 or evalue < 0:
                raise ValueError(f"Invalid Swiss-Prot alignment for {gene}")
            if evalue <= MAX_EVALUE and any(
                min(end, hi) - max(start, lo) + 1 >= MIN_DISCOVERY_OVERLAP * (hi - lo + 1)
                for lo, hi in intervals[gene]
            ):
                matches[gene].append((score, evalue, accession(hit["target"]), start, end))
    return matches


def classify(genes, hits, metadata, compatible, outside):
    accessions = {hit[2] for matches in hits.values() for hit in matches}
    membership = read_reference_membership(metadata, accessions, compatible, outside)
    for gene in genes:
        matches = sorted(hits.get(gene, ()), key=lambda hit: (-hit[0], hit[1], hit[2], hit[3], hit[4]))
        best_target = best_other = None
        target_ids = other_ids = ()
        focal = next((hit for hit in matches if membership.get(hit[2], ((), ()))[0]), None)
        if focal:
            start, end = focal[3:]
            # A fusion partner on a different part of the candidate does not
            # compete with the target-family alignment.
            matches = [hit for hit in matches
                       if min(end, hit[4]) - max(start, hit[3]) + 1
                       >= MIN_DISCOVERY_OVERLAP * min(end - start + 1, hit[4] - hit[3] + 1)]
        score_cutoff = matches[0][0] * NEAR_BEST_FRACTION if matches else 0
        for hit in matches:
            if hit[0] < score_cutoff:
                break
            matching, alternate = membership.get(hit[2], ((), ()))
            if matching and best_target is None:
                best_target, target_ids = hit, sorted(matching)
            if alternate and best_other is None:
                best_other, other_ids = hit, sorted(alternate)
        if best_target and best_other:
            status = "competing_families"
        elif best_target:
            status = "target_supported"
        elif best_other:
            status = "other_family_supported"
        else:
            status = "unresolved"
        yield gene, status, best_target, target_ids, best_other, other_ids


def describe(hit, ids):
    if hit is None:
        return [""] * 6
    score, evalue, ref, start, end = hit
    return [ref, ";".join(ids), score, evalue, start, end]


def ranked_genes(path):
    genes = path.read_text().splitlines()
    if len(genes) != len(set(genes)) or not all(genes):
        raise ValueError("Expected unique ranked candidate IDs")
    return genes


def screen(args, genes, hits, context):
    target, compatible, outside = context
    if target and genes:
        bounds = discovery_intervals(args.discovery, set(genes), args.min_coverage)
        aligned_hits = overlapping_hits(hits, bounds)
        decisions = list(classify(genes, aligned_hits, args.metadata, compatible, outside))
    else:
        decisions = [(gene, "unmapped_query", None, (), None, ()) for gene in genes]

    supported = [gene for gene, status, *_ in decisions if status == "target_supported"]
    provisional = [gene for gene, status, *_ in decisions
                   if status in {"unresolved", "competing_families", "unmapped_query"}]
    retained = set((supported + provisional)[:args.cap])
    with args.audit.open("w", newline="") as audit, args.selected.open("w") as selected:
        writer = csv.writer(audit, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "target_family", "status", "retained",
                         "target_accession", "target_family_ids", "target_bits", "target_evalue",
                         "target_start", "target_end",
                         "other_accession", "other_family_ids", "other_bits", "other_evalue",
                         "other_start", "other_end"])
        for gene, status, t, t_ids, o, o_ids in decisions:
            keep = gene in retained
            writer.writerow([gene, target, status, int(keep), *describe(t, t_ids), *describe(o, o_ids)])
            if keep:
                selected.write(gene + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("entries", "tree", "mapping", "metadata", "genes", "discovery", "sequence-store",
                 "swissprot-db", "audit", "selected"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--min-coverage", required=True, type=float)
    parser.add_argument("--cap", required=True, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--memory-limit", default="0", help="MMseqs2 prefilter memory limit.")
    args = parser.parse_args()
    if not 0 <= args.min_coverage <= 1 or min(args.cap, args.threads) < 1:
        parser.error("Invalid candidate coverage, cap or CPU count")
    context = read_family_context(args.entries, args.tree, args.mapping, args.query)
    genes = ranked_genes(args.genes)
    if context[0] is None or not genes:
        screen(args, genes, None, context)
        return
    target_db = args.swissprot_db.resolve()
    with TemporaryDirectory(prefix="query-family-swissprot-", dir=".") as directory:
        prefix = Path(directory)
        query_db, result_db, hits = (prefix / name for name in ("query", "result", "hits.tsv"))
        fasta = prefix / "candidates.fa"
        subprocess.run([sys.executable, str(Path(__file__).with_name("fasta_sequence_store.py")),
                        "extract", "--database", str(args.sequence_store), "--pattern-file", str(args.genes),
                        "--output", str(fasta), "--require-all"], check=True)
        subprocess.run(["mmseqs", "createdb", str(fasta), str(query_db)], check=True)
        subprocess.run(["mmseqs", "search", str(query_db), str(target_db), str(result_db), str(prefix / "tmp"),
                        "--threads", str(args.threads), "--max-seqs", "200", "-s", "7.5",
                        "--split-memory-limit", args.memory_limit,
                        "-e", str(MAX_EVALUE)], check=True)
        subprocess.run(["mmseqs", "convertalis", str(query_db), str(target_db), str(result_db), str(hits),
                        "--format-output", ",".join(RESULT_FIELDS), "--threads", str(args.threads)], check=True)
        screen(args, genes, hits, context)


if __name__ == "__main__":
    main()
