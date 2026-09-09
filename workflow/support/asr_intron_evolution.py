#!/usr/bin/env python3
"""Adapt GFF intron counts to NWKIT fixed-Q ASR and its probability-tree plot.

NWKIT owns inference, tree parsing, annotations, and drawing. This adapter only
translates the GFF table and joins observed counts back to the native ASR rows.
The input dated tree and its branch lengths are used without modification.
"""

import argparse
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

STATES = ("intron_absent", "intron_present")
MISSING_COUNTS = {"", "NA", "NaN", "nan"}


def nonnegative_rate(value):
    rate = float(value)
    if not math.isfinite(rate) or rate < 0:
        raise argparse.ArgumentTypeError("Rates must be finite and non-negative.")
    return rate


def read_intron_counts(path):
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"gene_id", "num_intron"}
    if not required.issubset(table.columns):
        raise ValueError("GFF table must contain gene_id and num_intron columns.")
    if table.gene_id.eq("").any() or table.gene_id.duplicated().any():
        raise ValueError("GFF gene_id values must be non-empty and unique.")
    counts = table.num_intron.str.strip().replace(dict.fromkeys(MISSING_COUNTS, None))
    counts = pd.to_numeric(counts, errors="raise")
    observed = counts.dropna()
    if any(not math.isfinite(value) or value < 0 or value != int(value) for value in observed):
        raise ValueError("num_intron must be a non-negative integer or missing.")
    return pd.Series(pd.array(counts, dtype="Int64"), index=table.gene_id, name="num_intron")


def run_asr(args):
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction, validate_output_targets

    outputs = {
        source: str(Path(args.output_prefix).absolute()) + suffix
        for source, suffix in (
            ("model.tsv", ".model.tsv"), ("tree.nhx", ".nhx"),
            ("plot.pdf", ".pdf"), ("summary.tsv", ".tsv"),
        )
    }
    validate_outputs_do_not_replace_inputs(
        [("tree", args.tree_file), ("traits", args.trait_file)], list(outputs.items())
    )
    validate_output_targets(outputs.values())
    counts = read_intron_counts(args.trait_file)
    prefix = Path(args.output_prefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="intron-asr-", dir=prefix.parent) as tmp:
        work = Path(tmp)
        traits = pd.DataFrame({"leaf_name": counts.index, "intron": [
            "NA" if pd.isna(value) else STATES[int(value > 0)] for value in counts
        ]})
        traits.to_csv(work / "traits.tsv", sep="\t", index=False)
        gain, loss = args.intron_gain_rate, args.retrotransposition_rate
        q = pd.DataFrame([[-gain, gain], [loss, -loss]], index=STATES, columns=STATES)
        q.to_csv(work / "q.tsv", sep="\t", index_label="state")
        subprocess.run([
            "nwkit", "asr", "--infile", str(Path(args.tree_file).resolve()),
            "--trait", str(work / "traits.tsv"), "--state-column", "intron",
            "--trait-type", "discrete", "--states", ",".join(STATES),
            "--model", "CUSTOM", "--rate-matrix", str(work / "q.tsv"),
            "--root-prior", "equal", "--target", "all", "--output", "probabilities",
            # A family's GFF table may still contain genes removed by tree pruning.
            "--unmatched", "ignore", "--outfile", str(work / "summary.tsv"),
            "--model-out", str(work / "model.tsv"), "--tree-out", str(work / "raw.nhx"),
            "--tree-annotation", "all",
        ], check=True)
        summary = pd.read_csv(work / "summary.tsv", sep="\t", dtype={"name": str}, keep_default_na=False)
        # Only observed leaves carry counts; internal names may equal gene IDs.
        summary["num_intron"] = pd.array(
            summary["name"].map(counts).where(summary.node_class.eq("leaf")), dtype="Int64"
        )
        summary.to_csv(work / "summary.tsv", sep="\t", index=False, na_rep="NA")
        # Human-readable legend keys while retaining all native ASR annotations.
        subprocess.run([
            "nwkit", "transfer", "--infile", str(work / "raw.nhx"),
            "--infile2", str(work / "raw.nhx"),
            "--property-map", "asr_p_intron_absent=Absent",
            "--property-map", "asr_p_intron_present=Present",
            "--outfile", str(work / "tree.nhx"),
        ], check=True)
        subprocess.run([
            "nwkit", "draw", "--infile", str(work / "tree.nhx"),
            "--node-pie-properties", "Absent,Present", "--node-pie-target", "all",
            "--property-color", "Absent=#D9D9D9", "--property-color", "Present=#222222",
            "--node-label-property", "name", "--node-label-target", "root,intnode",
            "--species-overlap-node-plot", "no", "--figure-width", "7.2",
            "--figure-height", str(max(3.0, 0.26 * summary.node_class.eq("leaf").sum() + 1.0)),
            "--tip-label-wrap", "auto", "--tip-spacing", "label-aware",
            "--outfile", str(work / "plot.pdf"),
        ], check=True)
        # NWKIT restores the entire previous set on a handled installation failure.
        with output_transaction(outputs.values()) as staged:
            for source, target in outputs.items():
                shutil.copyfile(work / source, staged[target])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree-file", required=True)
    parser.add_argument("--trait-file", required=True)
    parser.add_argument("--intron-gain-rate", required=True, type=nonnegative_rate)
    parser.add_argument("--retrotransposition-rate", required=True, type=nonnegative_rate)
    parser.add_argument("--output-prefix", required=True)
    run_asr(parser.parse_args())


if __name__ == "__main__":
    main()
