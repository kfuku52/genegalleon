"""Check the rooted copy split with dS from a common five-sequence alignment."""

import io
import itertools
import os
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from functools import lru_cache
from math import isfinite
from tempfile import TemporaryDirectory

import numpy as np
from Bio import SeqIO
from cdskit.backalign import backalign_sequence_strings

try:
    from pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS
    from parallel_io import ordered_bounded_results
except ImportError:
    from .pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS
    from .parallel_io import ordered_bounded_results

SLOTS = tuple(itertools.combinations(range(5), 2))
EXPECTED = (4, (0, 1), (2, 3))
ROOTED_COPY_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b", "contig", "window", "family_id",
                      "gene_a", "gene_b", "gene_a2", "gene_b2", "outgroup_species", "outgroup_gene", "common_codons",
                      "pair_statuses", "ten_ds", "expected_residual", "best_alternative_residual", "residual_gap",
                      "shorter_copy_edge", "status")
ROOTED_COPY_SUMMARY = (*ROOTED_COPY_FIELDS[:4], "sampled_quartets", "complete_quartets", "num_expected_best",
                       "num_other_best", "num_tied", "complete_families", "mixed_families")


@lru_cache(maxsize=1)
def designs():
    result = {}
    for middle in range(5):
        a, b, c, d = (i for i in range(5) if i != middle)
        for left, right in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            splits = (*((i,) for i in range(5)), left, right)
            result[middle, left, right] = np.array([[(x in side) != (y in side) for side in splits]
                                                  for x, y in SLOTS], dtype=float)
    return result


def rooted_fit(values):
    from scipy.optimize import nnls

    distances = np.asarray(values)
    models = designs()
    edges, expected = nnls(models[EXPECTED], distances)
    alternative = min(nnls(matrix, distances)[1] for key, matrix in models.items() if key != EXPECTED)
    gap = alternative - expected
    tolerance = 1e-10 * max(1., np.linalg.norm(distances))
    status = "expected_best" if gap > tolerance else "other_best" if gap < -tolerance else "tied"
    return expected, alternative, gap, min(edges[5:]), status


def common_alignment(records, genes, code, env):
    ids = tuple(map(str, range(5)))
    fasta = "".join(f">{name}\n{protein}\n" for name, (_, protein) in zip(ids, records, strict=True))
    aligned = subprocess.run(["mafft", *MAFFT_ARGS, "-"], input=fasta, capture_output=True,
                             text=True, check=True, env=env)
    rows = list(SeqIO.parse(io.StringIO(aligned.stdout), "fasta"))
    peptides = {row.id: str(row.seq).upper() for row in rows}
    if len(rows) != 5 or peptides.keys() != set(ids) or len({len(s) for s in peptides.values()}) != 1:
        raise ValueError(f"Invalid five-copy MAFFT alignment: {genes}")
    if any(peptides[name].replace("-", "") != protein for name, (_, protein) in zip(ids, records, strict=True)):
        raise ValueError(f"MAFFT changed a protein sequence: {genes}")
    codons = [backalign_sequence_strings(dna, peptides[name], code, gene, False)
              for (dna, _), name, gene in zip(records, ids, genes, strict=True)]
    allowed = frozenset("ACGT")
    keep = [i for i in range(0, len(codons[0]), 3) if all(allowed.issuperset(s[i:i + 3]) for s in codons)]
    n = len(keep)
    floor = max(60, max((7 * len(dna) + 29) // 30 for dna, _ in records))
    return n, tuple("".join(s[i:i + 3] for i in keep) for s in codons) if n >= floor else ()


class RootedCopyKs:
    def __init__(self, tree, codes):
        from nwkit.clade_index import CladeIndex

        index = CladeIndex(tree)
        leaves = {tip.name: tip for tip in tree.leaves()}
        self.sisters = {index.clade_id_for_node(node):
                        {sp: tree.get_distance(node.up, leaves[sp]) for sp in node.up.leaf_names()
                         if sp in codes and sp not in index.names_for_node(node)}
                        for node in tree.descendants() if not node.is_leaf}
        self.codes = codes
        self.outgroups = {}

    def observe(self, family, tree, membership, hits):
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        family_counts = Counter(membership.values())
        for event, a, b, _, _, genes, *_ in hits:
            key = event, family, genes
            if key in self.outgroups or a not in self.codes or self.codes.get(b) != self.codes[a]:
                continue
            duplication = tree.common_ancestor([leaves[gene] for gene in genes])
            if duplication.up is None:
                continue
            sisters = self.sisters[event]
            options = [(sisters[sp], sp, gene) for child in duplication.up.children if child is not duplication
                       for gene in child.leaf_names() if (sp := membership[gene]) in sisters
                       and family_counts[sp] == 1 and self.codes[sp] == self.codes[a]]
            if options:
                _, name, gene = min(options)
                self.outgroups[key] = name, gene

    def selection(self, selected, targets):
        for fields, _ in selected:
            event, _, a, b, _, _, family, *genes = fields
            outside = self.outgroups.get((event, family, tuple(genes)))
            if outside:
                sp, gene = outside
                for name, identifier in zip((a, b, a, b, sp), (*genes, gene), strict=True):
                    targets[name].add(identifier)
            yield fields, outside


def run_rooted_copy_ks(selection, sequences, codes, cpus, output, method, method_label, run_tool, write_table, read_table):
    rows, jobs = [], defaultdict(list)
    with TemporaryDirectory(prefix="gg-wgd-copy-mafft-", dir="/tmp") as scratch:
        env = os.environ if os.environ.get("MAFFT_TMPDIR") else {**os.environ, "MAFFT_TMPDIR": scratch}

        def align(fields, outside):
            if outside is None:
                return "no_single_sister_copy", 0, ()
            a, b = fields[2:4]
            genes = (*fields[7:], outside[1])
            records = [sequences[sp].get(gene) for sp, gene in zip((a, b, a, b, outside[0]), genes, strict=True)]
            if None in records:
                return "cds_unavailable", 0, ()
            count, codons = common_alignment(records, genes, codes[a], env)
            return ("dS_incomplete" if codons else "low_common_codon_coverage"), count, codons

        with ThreadPoolExecutor(max_workers=cpus) as pool:
            aligned = ordered_bounded_results(pool, align, selection, max_pending=cpus * 4)
            with closing(aligned):
                for i, ((fields, outside), (status, n, codons)) in enumerate(zip(selection, aligned, strict=True)):
                    metadata = (*fields, *(outside or (None, None)), n)
                    rows.append(dict(zip(ROOTED_COPY_FIELDS, (*metadata, *(None,) * 6, status), strict=True)))
                    if codons:
                        genes = (*fields[7:], outside[1])
                        jobs[codes[fields[2]]].extend(dict(zip(PAIR_COLUMNS, (f"{i}:{a}:{b}", genes[a], genes[b],
                                                                             codons[a], codons[b]), strict=True))
                                                       for a, b in SLOTS)

    reports = {}
    for code, records in sorted(jobs.items()):
        pairs = output / f"cross_ks_rooted_alignments.code{code}.tsv"
        table = output / f"cross_ks_rooted_pair_ds.code{code}.tsv"
        write_table(pairs, records, PAIR_COLUMNS)
        run_tool(["cdskit", "dnds", "--pairs_file", pairs, "--method", method, "--codon_table", code,
                  "--threads", cpus, "--out_file", table], output, f"cross_ks_rooted_dnds.code{code}")
        wanted = {row["pair_id"] for row in records}
        for row in read_table(table):
            key = row["pair_id"]
            if (key not in wanted or row.get("method") != method_label
                    or int(row["retained_codons"]) != rows[int(key.split(":", 1)[0])]["common_codons"]):
                raise ValueError("Unexpected rooted-copy dS pair, method or retained codon count")
            wanted.remove(key)
            reports[key] = row
        if wanted:
            raise ValueError("Incomplete rooted-copy dS report")
    for i, row in enumerate(rows):
        if row["status"] != "dS_incomplete":
            continue
        measurements = [reports[f"{i}:{a}:{b}"] for a, b in SLOTS]
        values = []
        states = []
        for item in measurements:
            status, ds = item["status"], None
            if status == "ok":
                sites, ds = float(item["S"]), float(item["dS"])
                if not (isfinite(ds) and ds >= 0 and isfinite(sites)):
                    raise ValueError("Invalid rooted-copy dS report value")
                if sites < 20:
                    status, ds = "low_codon_information", None
            states.append(status)
            values.append(ds)
        row["pair_statuses"] = ",".join(states)
        row["ten_ds"] = ",".join("NA" if value is None else str(value) for value in values)
        if None not in values:
            row.update(zip(ROOTED_COPY_FIELDS[-5:], rooted_fit(values), strict=True))
    return rows


def rooted_copy_summary(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[field] for field in ROOTED_COPY_FIELDS[:4])].append(row)
    complete = {"expected_best", "other_best", "tied"}
    for key, group in sorted(groups.items()):
        states = Counter(row["status"] for row in group)
        families = defaultdict(set)
        for row in group:
            if row["status"] in complete:
                families[row["family_id"]].add(row["status"])
        counts = (states[s] for s in ("expected_best", "other_best", "tied"))
        values = (*key, len(group), sum(states[s] for s in complete), *counts, len(families),
                  sum(len(status) > 1 for status in families.values()))
        yield dict(zip(ROOTED_COPY_SUMMARY, values, strict=True))
