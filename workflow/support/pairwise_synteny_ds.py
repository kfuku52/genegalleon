"""Codon-alignment preparation for synonymous-divergence coloring of JCVI anchors."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import random
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    from fasta_sequence_store import fasta_records
except ImportError:
    from .fasta_sequence_store import fasta_records

MAFFT_ARGS = ("--auto", "--amino", "--thread", "1", "--quiet")
PAIR_COLUMNS = ("pair_id", "target_gene", "query_gene", "sequence_1", "sequence_2")
DS_STATUSES = frozenset(("ok", "saturated", "no_synonymous_sites", "no_aligned_sense_codons", "internal_stop"))


def ds_tool_identity():
    try:
        import cdskit
        import cdskit.atomicio
        import cdskit.backalign
        import cdskit.codonutil
        import cdskit.dnds
        import cdskit.tsvio
        import cdskit.util
        import numpy
    except ImportError as exc:
        raise RuntimeError("dS coloring requires the CDSKIT dnds command; update CDSKIT and the GeneGalleon runtime") from exc
    files = (Path(__file__), Path(__file__).with_name("fasta_sequence_store.py"), *(Path(module.__file__) for module in (
        cdskit.dnds, cdskit.codonutil, cdskit.backalign, cdskit.tsvio, cdskit.util, cdskit.atomicio)))
    mafft = shutil.which("mafft")
    if mafft is None:
        raise RuntimeError("dS coloring requires MAFFT in the GeneGalleon runtime")
    return {"cdskit": cdskit.__version__, "numpy": numpy.__version__, "method": cdskit.dnds.METHOD,
            "mafft_version": subprocess.check_output(["mafft", "--version"], stderr=subprocess.STDOUT, text=True).strip(),
            "mafft_executable_sha256": hashlib.sha256(Path(mafft).read_bytes()).hexdigest(),
            "source_hashes": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}


def load_cds(path, analysis, side, species, code):
    """Match selected synteny IDs and verify their CDS translation exactly."""
    if code not in CodonTable.generic_by_id:
        raise ValueError("Unknown NCBI genetic code")
    records = {}
    for identifier, _, sequence in fasta_records(Path(path)):
        if identifier in records:
            raise ValueError(f"Duplicate CDS identifier: {identifier}")
        if not sequence.isascii() or set(sequence.upper()) - set("ACGTURYSWKMBDHVN"):
            raise ValueError(f"Invalid CDS alphabet: {identifier}")
        sequence = sequence.upper().replace("U", "T")
        if not sequence or len(sequence) % 3:
            raise ValueError(f"Invalid CDS frame: {identifier}")
        protein = str(Seq(sequence).translate(table=code))
        if protein.endswith("*"):
            sequence, protein = sequence[:-3], protein[:-1]
        if not protein or "*" in protein:
            raise ValueError(f"Empty CDS translation or internal stop: {identifier}")
        records[identifier] = sequence, protein
    proteins = {}
    for identifier, _, sequence in fasta_records(analysis / f"{side}.pep"):
        if identifier in proteins or not sequence:
            raise ValueError(f"Duplicate or empty synteny protein: {identifier}")
        proteins[identifier] = sequence
    result = {}
    with (analysis / f"{side}.id_map.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fields = reader.fieldnames or []
        if (any(not field.strip() for field in fields) or len(fields) != len(set(fields))
                or not {"original_id", "jcvi_id", "status"}.issubset(fields)):
            raise ValueError("Invalid synteny ID-map header")
        seen_originals = set()
        for row in reader:
            if None in row or any(value is None for value in row.values()) or not row["original_id"].strip():
                raise ValueError("Invalid synteny ID-map row")
            if row["original_id"] in seen_originals:
                raise ValueError(f"Duplicate original ID in synteny ID-map: {row['original_id']}")
            seen_originals.add(row["original_id"])
            if row["status"] != "selected":
                continue
            original = row["original_id"]
            short = original.removeprefix(species + "_")
            matches = {original, short, species + "_" + short} & records.keys()
            if len(matches) != 1:
                raise ValueError(f"Expected one matching CDS for selected protein: {original}")
            record = records[matches.pop()]
            if row["jcvi_id"] not in proteins or row["jcvi_id"] in result:
                raise ValueError("Duplicate or unknown selected synteny protein")
            if record[1] != proteins[row["jcvi_id"]]:
                raise ValueError(f"CDS translation differs from synteny protein: {original}")
            result[row["jcvi_id"]] = record
    if result.keys() != proteins.keys():
        raise ValueError("ID map does not cover every selected synteny protein")
    return result


def anchor_pairs(path):
    pairs = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f"Invalid anchor: {line.rstrip()}")
                pairs.add(tuple(fields[:2]))
    if not pairs:
        raise ValueError("No anchor pairs to align")
    return sorted(pairs)


def align_pair(pair, sequences, code):
    from cdskit.backalign import backalign_sequence_strings

    records = [sequences[i][gene] for i, gene in enumerate(pair)]
    fasta = "".join(f">{name}\n{record[1]}\n" for name, record in zip(("target", "query"), records, strict=True))
    alignment = subprocess.run(["mafft", *MAFFT_ARGS, "-"], input=fasta,
                               capture_output=True, text=True, check=True)
    rows = list(SeqIO.parse(io.StringIO(alignment.stdout), "fasta"))
    aligned = {record.id: str(record.seq).upper() for record in rows}
    if len(rows) != 2 or set(aligned) != {"target", "query"} or len(aligned["target"]) != len(aligned["query"]):
        raise ValueError(f"Invalid MAFFT pair alignment: {pair}")
    if any(aligned[name].replace("-", "") != record[1] for name, record in zip(("target", "query"), records, strict=True)):
        raise ValueError(f"MAFFT changed a protein sequence: {pair}")
    codons = [backalign_sequence_strings(record[0], aligned[name], code, gene, False)
              for name, gene, record in zip(("target", "query"), pair, records, strict=True)]
    return dict(zip(PAIR_COLUMNS, (pair[0] + "|" + pair[1], *pair, *codons), strict=True))


def prepare_pairs(analysis, cds_paths, species, code, output, cpus=1, limit=0):
    from cdskit.atomicio import atomic_text_writer, validate_distinct_paths

    if cpus < 1 or limit < 0:
        raise ValueError("cpus must be positive; limit must be nonnegative")
    inputs = [analysis, *cds_paths]
    validate_distinct_paths(inputs, [output])
    sequences = [load_cds(path, analysis, side, name, code)
                 for path, side, name in zip(cds_paths, ("target", "query"), species, strict=True)]
    pairs = anchor_pairs(analysis / "target.query.lifted.anchors")
    if limit:
        pairs = sorted(random.Random(1).sample(pairs, min(limit, len(pairs))))
    if any(pair[i] not in sequences[i] for pair in pairs for i in (0, 1)):
        raise ValueError("An anchor gene has no selected CDS")
    with atomic_text_writer(output, newline="") as handle, ThreadPoolExecutor(max_workers=cpus) as workers:
        writer = csv.DictWriter(handle, fieldnames=PAIR_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, row in enumerate(workers.map(lambda pair: align_pair(pair, sequences, code), pairs), 1):
            writer.writerow(row)
            if i % 100 == 0:
                print(f"Aligned {i}/{len(pairs)} anchor pairs", flush=True)
        validate_distinct_paths(inputs, [output])
    return len(pairs)


def read_ds(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fields = reader.fieldnames or []
        if (any(not field.strip() for field in fields) or len(fields) != len(set(fields))
                or not {"pair_id", "dS", "status"}.issubset(fields)):
            raise ValueError("Invalid dS report header")
        try:
            rows = list(reader)
        except csv.Error as exc:
            raise ValueError("Invalid dS report TSV") from exc
    lookup = {}
    for row in rows:
        if (None in row or any(value is None for value in row.values()) or not row["pair_id"].strip()
                or row["status"] not in DS_STATUSES):
            raise ValueError("Invalid dS report row")
        if row["pair_id"] in lookup:
            raise ValueError("Duplicate dS pair ID")
        value = float(row["dS"]) if row["dS"] else math.nan
        if (row["dS"] and (not math.isfinite(value) or value < 0 or row["status"] != "ok")) or (row["status"] == "ok" and not row["dS"]):
            raise ValueError("Invalid estimable dS value or inconsistent status")
        lookup[row["pair_id"]] = value
    return fields, rows, lookup


def estimate_ds(plan, output, cpus):
    from collections import Counter
    from types import SimpleNamespace

    from cdskit.codonutil import CODON_SEMANTICS_VERSION
    from cdskit.dnds import METHOD, REPORT_COLUMNS, dnds_main
    from cdskit.tsvio import TSV_REPORT_SCHEMA_VERSION

    output.mkdir(parents=True)
    for pair in plan["pairs"]:
        directory = output / pair["analysis_id"]
        directory.mkdir()
        analysis = Path(plan["workspace"]) / "output/genome_evolution/synteny/analysis" / pair["analysis_id"]
        sources = pair["ds"]
        code = sources["target"]["genetic_code"]
        if code != sources["query"]["genetic_code"]:
            raise ValueError("dS estimation requires the same genetic code for both species")
        audit = directory / "aligned_pairs.tsv"
        count = prepare_pairs(analysis, [sources[side]["fasta"] for side in ("target", "query")],
                              [pair[f"{side}_species"] for side in ("target", "query")], code, audit, cpus)
        dnds_main(SimpleNamespace(pairs_file=str(audit), outfile=str(directory / "ds.tsv"), codontable=code, threads=cpus))
        fields, rows, lookup = read_ds(directory / "ds.tsv")
        expected = {a + "|" + b for a, b in anchor_pairs(analysis / "target.query.lifted.anchors")}
        if fields != list(REPORT_COLUMNS) or set(lookup) != expected or len(rows) != count:
            raise ValueError("CDSKIT did not retain all unique anchor pairs")
        if any(row["method"] != METHOD or int(row["codon_table"]) != code for row in rows):
            raise ValueError("CDSKIT report method or genetic code differs from the requested estimator")
        if any(row["schema_version"] != str(TSV_REPORT_SCHEMA_VERSION)
               or row["codon_semantics_version"] != str(CODON_SEMANTICS_VERSION) for row in rows):
            raise ValueError("CDSKIT report schema or codon semantics differs from the requested estimator")
        with audit.open("rb") as source, gzip.open(directory / "aligned_pairs.tsv.gz", "wb") as dest:
            shutil.copyfileobj(source, dest)
        audit.unlink()  # Only our scratch uncompressed copy; the audited alignment remains.
        summary = {"unique_anchor_pairs": count, "status_counts": dict(Counter(row["status"] for row in rows)),
                   "genetic_code": code, "tools": plan["ds_tools"], "alignment": list(MAFFT_ARGS),
                   "pairwise_deletion": "nonconcrete codons, gaps and terminal unconditional stops",
                   "saturated_policy": "NA estimate; retain diagnostic separately; never color it as dS"}
        (directory / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def render_ds_dotplot(directory, ds_directory, pair, fmt, color_max, filtered=False):
    if not math.isfinite(color_max) or color_max <= 0:
        raise ValueError("dS color limit must be finite and positive")
    os.environ.setdefault("MPLCONFIGDIR", str(directory / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from jcvi.formats.bed import Bed
    from jcvi.graphics.dotplot import plot_breaks_and_labels
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    _fields, _rows, lookup = read_ds(ds_directory / "ds.tsv")
    prefix = "dotplot." if filtered else ""
    beds = [Bed(str(directory / f"{prefix}{side}.bed"), sorted=False) for side in ("target", "query")]
    orders = [bed.order for bed in beds]
    points = []
    with (directory / ("dotplot.anchors" if filtered else "target.query.lifted.anchors")).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("Invalid dotplot anchor")
            a, b = fields[:2]
            if a not in orders[0] or b not in orders[1] or a + "|" + b not in lookup:
                raise ValueError("An anchor has no BED coordinate or dS result")
            points.append((orders[0][a][0], orders[1][b][0], lookup[a + "|" + b]))
    if not points:
        raise ValueError("No anchors to render")
    data = np.asarray(points)
    fig = plt.figure(figsize=(12, 10))
    root = fig.add_axes((0, 0, 1, 1), frameon=False)
    root.set_axis_off()
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#b0b0b0")
    scatter = ax.scatter(data[:, 0], data[:, 1], c=np.ma.masked_invalid(data[:, 2]),
                         cmap=cmap, norm=Normalize(vmin=0, vmax=color_max, clip=True),
                         s=2, edgecolors="none", linewidths=0, plotnonfinite=True)
    if np.ma.count(scatter.get_offsets()) != 2 * len(points):
        raise ValueError("The renderer discarded an anchor")
    plot_breaks_and_labels(fig, root, ax, pair["target_species"].replace("_", " ") + " (gene rank)",
                          pair["query_species"].replace("_", " ") + " (gene rank)",
                          len(beds[0]), len(beds[1]), beds[0].get_breaks(), beds[1].get_breaks(),
                          chpf=False, usetex=False, sepcolor="#b0b0b0")
    cax = fig.add_axes((0.3, 0.015, 0.4, 0.016))
    fig.colorbar(scatter, cax=cax, orientation="horizontal", extend="max" if np.any(data[:, 2] > color_max) else "neither")
    cax.set_xlabel("dS (CDSKIT YN00; upper limit is display clipping only)", fontsize=10)
    missing = int(np.isnan(data[:, 2]).sum())
    root.legend(handles=[Line2D([], [], marker="o", linestyle="", color="#b0b0b0",
                                label=f"Unestimable / saturated: {missing:,} pairs")],
                loc="upper left", bbox_to_anchor=(0.09, 0.988), frameon=False, fontsize=9)
    root.text(0.5, 1.02, f"Pairwise synteny colored by dS ({len(points):,} gene pairs)", ha="center", fontsize=15)
    if fmt == "pdf":
        try:
            from pairwise_synteny_dotplot import save_pdf_dotplot
        except ImportError:
            from .pairwise_synteny_dotplot import save_pdf_dotplot
        save_pdf_dotplot(fig, root, ax, pair, directory / "dotplot.pdf", colorbar_ax=cax)
    else:
        fig.savefig(directory / f"dotplot.{fmt}", format=fmt, dpi=300, bbox_inches="tight")
    plt.close(fig)
    metadata = {"coordinate_system": "gene_rank", "color": "dS", "cmap": "viridis", "color_range": [0, color_max],
                "anchor_count": len(points), "missing_count": missing, "above_color_max_count": int(np.sum(data[:, 2] > color_max)),
                "downsampled": False, "filtered_by_dS": False, "ds_source": str(ds_directory / "ds.tsv")}
    (directory / "dotplot_ds.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, type=Path)
    parser.add_argument("--target-cds", required=True, type=Path)
    parser.add_argument("--query-cds", required=True, type=Path)
    parser.add_argument("--target-species", required=True)
    parser.add_argument("--query-species", required=True)
    parser.add_argument("--genetic-code", type=int, default=1)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.cpus < 1 or args.limit < 0:
        parser.error("cpus must be positive; limit must be nonnegative")
    count = prepare_pairs(args.analysis, (args.target_cds, args.query_cds),
                          (args.target_species, args.query_species), args.genetic_code, args.output, args.cpus, args.limit)
    print(json.dumps({"prepared_pairs": count, "genetic_code": args.genetic_code, "aligner": "MAFFT --auto --amino --thread 1"}))


if __name__ == "__main__":
    main()
