#!/usr/bin/env python3
"""Extend an existing test input with real, complete AHA genomic neighborhoods.

All existing biological CDS records and unrelated input files are preserved.
Only explicitly marked gg_dummy annotations/CDS are separated into fixtures.
The output must be a new directory; source and seed inputs are never modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

try:
    from build_minimal_test_dataset import fasta_iter, find_first_file, open_text, write_fasta_record
except ImportError:
    from .build_minimal_test_dataset import fasta_iter, find_first_file, open_text, write_fasta_record


def normalized(value):
    return value.replace("-", "_")


def attributes(text):
    return {k: unquote(v) for field in text.split(";") if "=" in field
            for k, v in [field.split("=", 1)]}


def gff_rows(path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) != 9:
                raise ValueError(f"Invalid GFF row in {path}: {line[:100]}")
            yield row


class GeneCatalog:
    """Resolve exact gene/transcript identifiers and their parent relationships."""

    def __init__(self, path):
        self.genes = {}
        self.parents = defaultdict(set)
        names = defaultdict(set)
        for row in gff_rows(path):
            a = attributes(row[8])
            identifier = a.get("ID", "")
            if row[2] == "gene":
                if not identifier or identifier in self.genes:
                    raise ValueError(f"Missing/duplicate gene ID: {identifier}")
                self.genes[identifier] = (row[0], int(row[3]), int(row[4]), row[6])
            if identifier:
                self.parents[identifier].update(p for p in a.get("Parent", "").split(",")
                                                if p and p != identifier)
            if row[2] in {"gene", "mRNA", "transcript"}:
                for key in ["ID", "Name", "Alias", "gene", "gene_id", "transcript_id"]:
                    for name in a.get(key, "").split(","):
                        if name:
                            names[normalized(name)].add(identifier)
        self.names = {}
        for name, identifiers in names.items():
            roots = set().union(*(self.roots(i) for i in identifiers))
            if roots:
                self.names[name] = roots

    def roots(self, identifier, trail=frozenset()):
        if identifier in self.genes:
            return {identifier}
        if identifier in trail:
            raise ValueError(f"Cyclic GFF parent: {identifier}")
        return set().union(*(self.roots(p, trail | {identifier})
                             for p in self.parents.get(identifier, ())))

    def row_roots(self, row):
        a = attributes(row[8])
        roots = self.roots(a.get("ID", ""))
        for parent in a.get("Parent", "").split(","):
            roots |= self.roots(parent)
        return roots

    def map_cds(self, records, species):
        mapping = {}
        for identifier in records:
            core = identifier.removeprefix(species + "_")
            roots = self.names.get(normalized(core), set())
            if len(roots) > 1:
                raise ValueError(f"Ambiguous CDS locus: {identifier}: {sorted(roots)}")
            if roots:
                mapping[identifier] = next(iter(roots))
        return mapping


def read_records(path):
    result = {}
    for header, identifier, sequence in fasta_iter(path):
        if identifier in result:
            raise ValueError(f"Duplicate FASTA ID: {identifier}")
        result[identifier] = (header, sequence)
    return result


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_file(root, kind, species):
    matches = list((root / kind).glob(species + "*"))
    if len(matches) != 1:
        raise ValueError(f"Expected one {kind} file for {species}, found {matches}")
    return find_first_file(root / kind, species)


def source_cds_order(files, mapping):
    """Use the same CDS bounds and identifier ordering as synteny_neighbors."""
    script = Path(__file__).with_name("gff2genestat.py")
    with tempfile.TemporaryDirectory(prefix="gg_real_gene_order_") as directory:
        output = Path(directory) / "genes.tsv"
        command = [sys.executable, str(script), "--mode", "gene_delim",
                   "--dir_gff", str(files["species_gff"].parent),
                   "--seqfile", str(files["species_cds"]), "--outfile", str(output), "--ncpu", "1"]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            raise ValueError(f"Full-source GFF parsing failed:\n{result.stdout}\n{result.stderr}")
        order = {}
        with output.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                identifier = row["gene_id"]
                if identifier not in mapping:
                    continue
                gene = mapping[identifier]
                start, end = sorted([int(row["start"]), int(row["end"])])
                if gene in order:
                    raise ValueError(f"Source CDS must have one representative per gene: {gene}")
                order[gene] = (row["chromosome"], start, end, identifier)
        if set(order) != set(mapping.values()):
            raise ValueError(f"Source GFF parser missed loci: {set(mapping.values()) - set(order)}")
        return order


def build_species(source, seed, out, species, anchors, neighbors, flank):
    files = {kind: input_file(source, kind, species)
             for kind in ["species_cds", "species_gff", "species_genome"]}
    seed_files = {kind: input_file(seed, kind, species) for kind in files
                  if list((seed / kind).glob(species + "*"))}
    if "species_cds" not in seed_files:
        raise ValueError(f"Seed CDS is required for {species}")
    catalog = GeneCatalog(files["species_gff"])
    source_cds = read_records(files["species_cds"])
    old_cds = read_records(seed_files["species_cds"])
    mapping = catalog.map_cds(source_cds, species)
    by_gene = defaultdict(set)
    for identifier, gene in mapping.items():
        by_gene[gene].add(identifier)

    dummy_rows = [row for row in gff_rows(seed_files["species_gff"])
                  if "coge_fid=gg_dummy_" in row[8]] if "species_gff" in seed_files else []
    dummy_names = {normalized(attributes(row[8])["Alias"]) for row in dummy_rows}
    dummies = {identifier for identifier in old_cds
               if normalized(identifier.removeprefix(species + "_")) in dummy_names}
    retained = set(old_cds) - dummies
    for identifier in retained:
        if identifier not in mapping:
            raise ValueError(f"Existing CDS has no source gene: {identifier}")
        if old_cds[identifier][1].upper() != source_cds[identifier][1].upper():
            raise ValueError(f"Existing/source CDS sequence differs: {identifier}")
    if not anchors <= retained:
        raise ValueError(f"Anchors missing from seed: {sorted(anchors - retained)}")
    if dummies & set(source_cds):
        raise ValueError("Marked dummy IDs collide with real source CDS")

    # Match synteny_neighbors.py: order loci with available CDS, once per gene.
    coding_order = source_cds_order(files, mapping)
    chrom_genes = defaultdict(list)
    for gene in by_gene:
        if coding_order[gene][0] != catalog.genes[gene][0]:
            raise ValueError(f"CDS and gene annotation disagree on chromosome: {gene}")
        chrom_genes[catalog.genes[gene][0]].append(gene)
    for genes in chrom_genes.values():
        genes.sort(key=lambda g: coding_order[g][1:])
    selected = {mapping[i] for i in retained}
    coverage = []
    anchor_windows = defaultdict(list)
    for identifier in sorted(anchors):
        gene = mapping[identifier]
        chrom, start, end, strand = catalog.genes[gene]
        genes = chrom_genes[chrom]
        index = genes.index(gene)
        left = genes[max(0, index - neighbors):index]
        right = genes[index + 1:index + 1 + neighbors]
        selected.update(left + right)
        neighborhood = left + [gene] + right
        anchor_windows[chrom].append([
            max(1, min(catalog.genes[g][1] for g in neighborhood) - flank),
            max(catalog.genes[g][2] for g in neighborhood) + flank,
        ])
        coverage.append(dict(species=species, anchor=identifier, source_seqid=chrom,
                             source_start=start, source_end=end, strand=strand,
                             left_count=len(left), right_count=len(right),
                             left_genes=left, right_genes=right))

    # Include complete gene models intersecting a window, never clip features.
    # Only selected targets receive padding; closure does not recursively pad.
    # Keep the entire interval between neighbors, including long intergenic DNA.
    # Separate per-gene snippets would create artificial contig boundaries and
    # prevent synteny_neighbors from seeing these genes as neighbors.
    windows = anchor_windows
    for gene in selected:
        chrom, start, end, _ = catalog.genes[gene]
        windows[chrom].append([max(1, start - flank), end + flank])
    changed = True
    while changed:
        changed = False
        for chrom, intervals in windows.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            windows[chrom] = merged
        for gene, (chrom, start, end, _) in catalog.genes.items():
            for interval in windows.get(chrom, []):
                if start <= interval[1] and end >= interval[0]:
                    selected.add(gene)
                    new = [min(start, interval[0]), max(end, interval[1])]
                    if interval != new:
                        interval[:] = new
                        changed = True

    paths = {}
    for kind in files:
        filename = seed_files.get(kind, files[kind]).name.removesuffix(".gz")
        if filename.endswith(".fa.masked"):
            filename = filename.removesuffix(".fa.masked") + ".masked.fa"
        paths[kind] = out / kind / filename
    window_rows = []
    seen_chroms = set()
    with paths["species_genome"].open("w") as handle:
        for _, chrom, sequence in fasta_iter(files["species_genome"]):
            if chrom not in windows:
                continue
            seen_chroms.add(chrom)
            for interval in windows[chrom]:
                interval[1] = min(interval[1], len(sequence))
                start, end = interval
                name = f"{chrom}:{start}-{end}"
                subseq = sequence[start - 1:end]
                write_fasta_record(handle, name, subseq)
                window_rows.append(dict(species=species, source_seqid=chrom,
                                        start=start, end=end, window_id=name,
                                        sequence_sha256=hashlib.sha256(subseq.encode()).hexdigest()))
    if set(windows) != seen_chroms:
        raise ValueError(f"Source genome lacks chromosomes: {set(windows) - seen_chroms}")

    with paths["species_gff"].open("w") as handle:
        handle.write("##gff-version 3\n")
        for row in gff_rows(files["species_gff"]):
            if not catalog.row_roots(row) & selected:
                continue
            start, end = int(row[3]), int(row[4])
            matches = [(s, e) for s, e in windows[row[0]] if s <= start <= end <= e]
            if len(matches) != 1:
                raise ValueError(f"Feature not contained in one window: {row}")
            offset, limit = matches[0]
            row[0] = f"{row[0]}:{offset}-{limit}"
            row[3:5] = [str(start - offset + 1), str(end - offset + 1)]
            handle.write("\t".join(row) + "\n")
    selected_cds = retained | {i for g in selected for i in by_gene[g]}
    with paths["species_cds"].open("w") as handle:
        for identifier in sorted(selected_cds):
            # Preserve the original header, sequence case and sequence for seed genes.
            header, sequence = old_cds.get(identifier, source_cds[identifier])
            write_fasta_record(handle, header, sequence)
    if dummies:
        fixture = out / "dataset_manifest" / "synthetic_synteny"
        fixture.mkdir(parents=True, exist_ok=True)
        with (fixture / f"{species}.fa").open("w") as handle:
            for identifier in sorted(dummies):
                write_fasta_record(handle, *old_cds[identifier])
        (fixture / f"{species}.gff").write_text(
            "##gff-version 3\n" + "".join("\t".join(row) + "\n" for row in dummy_rows))
    for kind, seed_file in seed_files.items():
        copied = out / kind / seed_file.name
        if copied != paths[kind]:
            copied.unlink()
    return dict(species=species, retained_ids=sorted(retained), removed_dummy_ids=sorted(dummies),
                retained_cds_sha256={i: hashlib.sha256(old_cds[i][1].encode()).hexdigest()
                                    for i in sorted(retained)},
                added_ids=sorted(selected_cds - retained),
                seed_files={k: dict(filename=p.name, sha256=sha256(p)) for k, p in seed_files.items()},
                source_files={k: dict(filename=p.name, sha256=sha256(p)) for k, p in files.items()},
                output_files={k: dict(filename=p.name, bytes=p.stat().st_size, sha256=sha256(p))
                              for k, p in paths.items()}), coverage, window_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pg", required=True, type=Path, help="Full original input workspace")
    parser.add_argument("--seed-pg", required=True, type=Path, help="Existing test input to preserve")
    parser.add_argument("--out-pg", required=True, type=Path, help="New output directory (must not exist)")
    parser.add_argument("--anchor-ids", required=True, type=Path, help="One AHA CDS ID per line")
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--flank-bp", type=int, default=5000)
    args = parser.parse_args()
    if args.neighbors < 1 or args.flank_bp < 0:
        parser.error("neighbors must be positive and flank-bp nonnegative")
    source, seed, out = (p.resolve() for p in (args.source_pg, args.seed_pg, args.out_pg))
    if out.exists() or seed in out.parents or source in out.parents:
        parser.error("out-pg must be a new directory outside source-pg and seed-pg")
    anchors = defaultdict(set)
    for identifier in args.anchor_ids.read_text().splitlines():
        if identifier.strip():
            anchors["_".join(identifier.split("_")[:2])].add(identifier.strip())
    if not anchors:
        parser.error("anchor-ids is empty")
    shutil.copytree(seed, out, ignore=shutil.ignore_patterns(".DS_Store"))
    manifest = out / "dataset_manifest"
    manifest.mkdir(exist_ok=True)
    # A failed rebuild must not leave the seed's completion manifest looking current.
    for filename in ["real_neighborhoods.json", "aha_coverage.tsv"]:
        (manifest / filename).unlink(missing_ok=True)
    (manifest / "aha_anchors.txt").write_text(
        "".join(identifier + "\n" for identifier in sorted(set().union(*anchors.values()))))
    summaries, coverage, windows = [], [], []
    for species in sorted(anchors):
        print(f"Extending {species}", flush=True)
        summary, rows, regions = build_species(source, seed, out, species,
                                               anchors[species], args.neighbors, args.flank_bp)
        summaries.append(summary)
        coverage.extend(rows)
        windows.extend(regions)
    replaced = {f"{kind}/{details['filename']}"
                for summary in summaries for kind, details in summary["seed_files"].items()}
    preserved = {}
    for path in sorted(seed.rglob("*")):
        relative = str(path.relative_to(seed))
        if (not path.is_file() or path.name == ".DS_Store" or relative in replaced
                or relative.startswith("dataset_manifest/")):
            continue
        original_hash = sha256(path)
        if sha256(out / relative) != original_hash:
            raise ValueError(f"Unrelated seed input changed: {relative}")
        preserved[relative] = original_hash
    (manifest / "real_neighborhoods.json").write_text(json.dumps(
        dict(neighbors=args.neighbors, flank_bp=args.flank_bp, species=summaries,
             ordering="CDS bounds and IDs from gff2genestat.py, matching synteny_neighbors.py",
             gff_parser_sha256=sha256(Path(__file__).with_name("gff2genestat.py")),
             coverage=coverage, windows=windows, preserved_files=preserved), indent=2) + "\n")
    with (manifest / "aha_coverage.tsv").open("w") as handle:
        keys = [k for k in coverage[0] if k not in {"left_genes", "right_genes"}]
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(coverage)
    print(f"Complete: {out}", flush=True)


if __name__ == "__main__":
    main()
