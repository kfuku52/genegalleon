#!/usr/bin/env python3

import argparse
import csv
import shlex
from pathlib import Path
import re


PACBIO_RE = re.compile(r"pacbio|smrt|sequel|revio", re.IGNORECASE)
ONT_RE = re.compile(r"nanopore|minion|gridion|promethion|flongle", re.IGNORECASE)
ILLUMINA_RE = re.compile(r"illumina|hiseq|miseq|nextseq|novaseq|iseq", re.IGNORECASE)
SHORT_OTHER_RE = re.compile(r"ion torrent|dnbseq|bgiseq|capillary|abi|454|ls454", re.IGNORECASE)
DIRECT_RNA_RE = re.compile(r"direct[ _-]*rna|\bdrna\b", re.IGNORECASE)

TEXT_FIELDS = (
    "instrument",
    "platform",
    "instrument_model",
    "exp_title",
    "design",
    "sample_title",
    "study_title",
    "lib_name",
    "library_name",
    "sample_description",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify amalgkit metadata rows into short-read / PacBio / ONT transcriptome inputs."
    )
    parser.add_argument("--metadata", required=True, help="Input amalgkit metadata TSV.")
    parser.add_argument("--classification-out", required=True, help="Output TSV with per-run classifications.")
    parser.add_argument("--summary-sh", required=True, help="Output shell fragment with aggregated classification variables.")
    return parser.parse_args()


def determine_instrument_field(fieldnames):
    for name in ("instrument", "platform", "instrument_model"):
        if name in fieldnames:
            return name
    return "none"


def classify_row(row):
    text_parts = []
    for field in TEXT_FIELDS:
        value = str(row.get(field, "") or "").strip()
        if value:
            text_parts.append(value)
    combined = " ".join(text_parts)

    if PACBIO_RE.search(combined):
        return "pacbio", "pacbio", "0"
    if ONT_RE.search(combined):
        if DIRECT_RNA_RE.search(combined):
            return "ont_direct_rna", "ont", "1"
        return "ont_cdna", "ont", "0"
    if ILLUMINA_RE.search(combined):
        return "short_read", "illumina", "0"
    if SHORT_OTHER_RE.search(combined):
        return "short_read", "short_other", "0"
    return "short_read", "unknown", "0"


def write_summary(path, values):
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        for key, value in values.items():
            handle.write(f"{key}={shlex.quote(str(value))}\n")


def main():
    args = parse_args()
    metadata_path = Path(args.metadata)
    classification_path = Path(args.classification_out)
    summary_path = Path(args.summary_sh)

    classification_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        instrument_field = determine_instrument_field(fieldnames)

        rows = []
        seen_runs = set()
        short_runs = set()
        pacbio_runs = set()
        ont_cdna_runs = set()
        ont_direct_rna_runs = set()

        for row in reader:
            run = str(row.get("run", "") or "").strip()
            if not run or run in seen_runs:
                continue
            seen_runs.add(run)

            read_class, platform_family, ont_direct_rna = classify_row(row)
            instrument_value = ""
            if instrument_field != "none":
                instrument_value = str(row.get(instrument_field, "") or "").strip()

            rows.append(
                {
                    "run": run,
                    "read_class": read_class,
                    "platform_family": platform_family,
                    "long_read": "1" if read_class != "short_read" else "0",
                    "ont_direct_rna": ont_direct_rna,
                    "lib_layout": str(row.get("lib_layout", "") or "").strip(),
                    "instrument_value": instrument_value,
                }
            )

            if read_class == "short_read":
                short_runs.add(run)
            elif read_class == "pacbio":
                pacbio_runs.add(run)
            elif read_class == "ont_cdna":
                ont_cdna_runs.add(run)
            elif read_class == "ont_direct_rna":
                ont_direct_rna_runs.add(run)

    with classification_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "run",
                "read_class",
                "platform_family",
                "long_read",
                "ont_direct_rna",
                "lib_layout",
                "instrument_value",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    has_pacbio = int(bool(pacbio_runs))
    has_ont_cdna = int(bool(ont_cdna_runs))
    has_ont_direct_rna = int(bool(ont_direct_rna_runs))
    has_ont = int(bool(ont_cdna_runs or ont_direct_rna_runs))
    has_short_reads = int(bool(short_runs))
    has_long_reads = int(bool(pacbio_runs or ont_cdna_runs or ont_direct_rna_runs))

    if has_pacbio and has_ont:
        input_class = "mixed_long_platforms"
    elif has_ont_cdna and has_ont_direct_rna:
        input_class = "ont_mixed_protocols"
    elif has_long_reads and has_short_reads:
        input_class = "hybrid_long_short"
    elif has_pacbio:
        input_class = "pacbio"
    elif has_ont_direct_rna:
        input_class = "ont_direct_rna"
    elif has_ont_cdna:
        input_class = "ont_cdna"
    elif has_short_reads:
        input_class = "short_read"
    else:
        input_class = "unknown"

    write_summary(
        summary_path,
        {
            "detected_metadata_run_count": len(seen_runs),
            "detected_short_read_run_count": len(short_runs),
            "detected_pacbio_run_count": len(pacbio_runs),
            "detected_ont_cdna_run_count": len(ont_cdna_runs),
            "detected_ont_direct_rna_run_count": len(ont_direct_rna_runs),
            "detected_has_long_reads": has_long_reads,
            "detected_has_short_reads": has_short_reads,
            "detected_has_pacbio": has_pacbio,
            "detected_has_ont": has_ont,
            "detected_has_ont_cdna": has_ont_cdna,
            "detected_has_ont_direct_rna": has_ont_direct_rna,
            "detected_input_class": input_class,
            "detected_metadata_instrument_field": instrument_field,
        },
    )


if __name__ == "__main__":
    main()
