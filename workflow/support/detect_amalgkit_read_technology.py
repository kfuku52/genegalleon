#!/usr/bin/env python3

import argparse
import csv
import re
import shlex
from pathlib import Path


SHORT_OTHER_RE = re.compile(r"ion torrent|dnbseq|bgiseq|capillary|abi|454|ls454", re.IGNORECASE)
SHORT_READ_SPOT_LENGTH_THRESHOLD = 1000
PLATFORM_REGEX_BY_FAMILY = {
    "ont": (
        r"oxford[\s_/-]*nanopore",
        r"\bnanopore\b",
        r"\bpromethion\b",
        r"\bgridion\b",
        r"\bminion\b",
        r"\bont\b",
    ),
    "pacbio": (
        r"\bpacbio\b",
        r"pacific[\s_/-]*biosciences",
        r"\bsmrt\b",
        r"\bsequel\b",
        r"\brevio\b",
        r"\brs[\s_-]*ii\b",
    ),
    "illumina": (
        r"\billumina\b",
        r"\bnovaseq\b",
        r"\bhiseq\b",
        r"\bnextseq\b",
        r"\bmiseq\b",
        r"\biseq\b",
        r"genome[\s_-]*analyzer",
    ),
}
DIRECT_RNA_PATTERNS = (
    r"\bdirect[\s_-]*rna\b",
    r"\bd[\s_-]*rna\b",
    r"\bdrna\b",
)

TEXT_FIELDS = (
    "platform",
    "instrument",
    "instrument_model",
    "protocol",
    "lib_strategy",
    "lib_selection",
    "lib_source",
    "lib_name",
    "library_name",
    "exp_title",
    "design",
    "sample_title",
    "study_title",
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


def normalize_text(value):
    text = str(value or "").strip()
    if text.lower() in {"", "nan", "none"}:
        return ""
    return re.sub(r"\s+", " ", text)


def normalize_search_text(text):
    return re.sub(r"[^0-9a-z]+", " ", normalize_text(text).lower()).strip()


def collect_text_fragments(row):
    fragments = []
    for field in TEXT_FIELDS:
        value = normalize_text(row.get(field, ""))
        if value:
            fragments.append(value)
    return fragments


def parse_positive_float(value):
    text = normalize_text(value)
    if text == "":
        return None
    try:
        parsed = float(text.replace(",", ""))
    except ValueError:
        return None
    if parsed <= 0:
        return None
    return parsed


def estimate_spot_length(row):
    spot_length = parse_positive_float(row.get("spot_length", ""))
    if spot_length is not None:
        return spot_length
    total_spots = parse_positive_float(row.get("total_spots", ""))
    total_bases = parse_positive_float(row.get("total_bases", ""))
    if total_spots is None or total_bases is None or total_spots <= 0:
        return None
    return total_bases / total_spots


def matches_search_patterns(search_text, patterns):
    if search_text == "":
        return False
    for pattern in patterns:
        if re.search(pattern, search_text):
            return True
    return False


def infer_platform_family(row):
    search_text = normalize_search_text(" ".join(collect_text_fragments(row)))
    for family in ("ont", "pacbio", "illumina"):
        if matches_search_patterns(search_text, PLATFORM_REGEX_BY_FAMILY[family]):
            return family
    if SHORT_OTHER_RE.search(search_text):
        return "short_other"

    lib_layout = normalize_text(row.get("lib_layout", "")).lower()
    if lib_layout == "paired":
        return "illumina"

    spot_length = estimate_spot_length(row)
    if spot_length is None:
        return "unknown"
    if spot_length <= SHORT_READ_SPOT_LENGTH_THRESHOLD:
        return "illumina"
    return "long_read_unknown"


def classify_row(row):
    search_text = normalize_search_text(" ".join(collect_text_fragments(row)))
    platform_family = infer_platform_family(row)

    if platform_family == "pacbio":
        return "pacbio", "pacbio", "0"
    if platform_family == "ont":
        if matches_search_patterns(search_text, DIRECT_RNA_PATTERNS):
            return "ont_direct_rna", "ont", "1"
        return "ont_cdna", "ont", "0"
    if platform_family == "illumina":
        return "short_read", "illumina", "0"
    if platform_family == "short_other":
        return "short_read", "short_other", "0"
    if platform_family == "long_read_unknown":
        return "long_read_unknown", "long_read_unknown", "0"
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
        long_read_unknown_runs = set()

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
            elif read_class == "long_read_unknown":
                long_read_unknown_runs.add(run)

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
    has_long_read_unknown = int(bool(long_read_unknown_runs))
    has_long_reads = int(bool(pacbio_runs or ont_cdna_runs or ont_direct_rna_runs or long_read_unknown_runs))

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
    elif has_long_read_unknown:
        input_class = "long_read_unknown"
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
            "detected_long_read_unknown_run_count": len(long_read_unknown_runs),
            "detected_has_long_reads": has_long_reads,
            "detected_has_short_reads": has_short_reads,
            "detected_has_pacbio": has_pacbio,
            "detected_has_ont": has_ont,
            "detected_has_ont_cdna": has_ont_cdna,
            "detected_has_ont_direct_rna": has_ont_direct_rna,
            "detected_has_long_read_unknown": has_long_read_unknown,
            "detected_input_class": input_class,
            "detected_metadata_instrument_field": instrument_field,
        },
    )


if __name__ == "__main__":
    main()
