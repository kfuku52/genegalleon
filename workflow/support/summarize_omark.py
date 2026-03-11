#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas


HEADER = [
    "species",
    "sum_file",
    "selected_clade",
    "conserved_hogs",
    "single_count",
    "duplicated_count",
    "duplicated_unexpected_count",
    "duplicated_expected_count",
    "missing_count",
    "single_pct",
    "duplicated_pct",
    "duplicated_unexpected_pct",
    "duplicated_expected_pct",
    "missing_pct",
    "protein_count",
    "consistent_count",
    "consistent_partial_count",
    "consistent_fragment_count",
    "inconsistent_count",
    "inconsistent_partial_count",
    "inconsistent_fragment_count",
    "contamination_count",
    "contamination_partial_count",
    "contamination_fragment_count",
    "unknown_count",
    "consistent_pct",
    "consistent_partial_pct",
    "consistent_fragment_pct",
    "inconsistent_pct",
    "inconsistent_partial_pct",
    "inconsistent_fragment_pct",
    "contamination_pct",
    "contamination_partial_pct",
    "contamination_fragment_pct",
    "unknown_pct",
    "top_detected_clade",
    "top_detected_taxid",
    "top_detected_proteins",
    "top_detected_pct",
    "detected_clade_count",
]

COMPLETENESS_COUNTS_RE = re.compile(
    r"^S:(?P<single>\d+),D:(?P<duplicated>\d+)\[U:(?P<duplicated_unexpected>\d+),E:(?P<duplicated_expected>\d+)\],M:(?P<missing>\d+)$"
)
COMPLETENESS_PCTS_RE = re.compile(
    r"^S:(?P<single>[0-9.]+)%,D:(?P<duplicated>[0-9.]+)%\[U:(?P<duplicated_unexpected>[0-9.]+)%,E:(?P<duplicated_expected>[0-9.]+)%\],M:(?P<missing>[0-9.]+)%$"
)
PROTEOME_COUNTS_RE = re.compile(
    r"^A:(?P<consistent>\d+)\[P:(?P<consistent_partial>\d+),F:(?P<consistent_fragment>\d+)\],"
    r"I:(?P<inconsistent>\d+)\[P:(?P<inconsistent_partial>\d+),F:(?P<inconsistent_fragment>\d+)\],"
    r"C:(?P<contamination>\d+)\[P:(?P<contamination_partial>\d+),F:(?P<contamination_fragment>\d+)\],"
    r"U:(?P<unknown>\d+)$"
)
PROTEOME_PCTS_RE = re.compile(
    r"^A:(?P<consistent>[0-9.]+)%\[P:(?P<consistent_partial>[0-9.]+)%,F:(?P<consistent_fragment>[0-9.]+)%\],"
    r"I:(?P<inconsistent>[0-9.]+)%\[P:(?P<inconsistent_partial>[0-9.]+)%,F:(?P<inconsistent_fragment>[0-9.]+)%\],"
    r"C:(?P<contamination>[0-9.]+)%\[P:(?P<contamination_partial>[0-9.]+)%,F:(?P<contamination_fragment>[0-9.]+)%\],"
    r"U:(?P<unknown>[0-9.]+)%$"
)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Summarize per-species OMArk .sum outputs into one TSV table."
    )
    parser.add_argument(
        "--omark_outdir",
        metavar="PATH",
        required=True,
        help="Directory containing per-species OMArk output folders.",
    )
    parser.add_argument(
        "--outfile",
        metavar="PATH",
        default="omark_summary.tsv",
        help="Output TSV path.",
    )
    return parser


def _default_row(sum_path: Path) -> dict:
    species = sum_path.parent.name if sum_path.parent.name else sum_path.stem
    return {
        "species": species,
        "sum_file": str(sum_path),
        "selected_clade": "-",
        "conserved_hogs": pandas.NA,
        "single_count": pandas.NA,
        "duplicated_count": pandas.NA,
        "duplicated_unexpected_count": pandas.NA,
        "duplicated_expected_count": pandas.NA,
        "missing_count": pandas.NA,
        "single_pct": pandas.NA,
        "duplicated_pct": pandas.NA,
        "duplicated_unexpected_pct": pandas.NA,
        "duplicated_expected_pct": pandas.NA,
        "missing_pct": pandas.NA,
        "protein_count": pandas.NA,
        "consistent_count": pandas.NA,
        "consistent_partial_count": pandas.NA,
        "consistent_fragment_count": pandas.NA,
        "inconsistent_count": pandas.NA,
        "inconsistent_partial_count": pandas.NA,
        "inconsistent_fragment_count": pandas.NA,
        "contamination_count": pandas.NA,
        "contamination_partial_count": pandas.NA,
        "contamination_fragment_count": pandas.NA,
        "unknown_count": pandas.NA,
        "consistent_pct": pandas.NA,
        "consistent_partial_pct": pandas.NA,
        "consistent_fragment_pct": pandas.NA,
        "inconsistent_pct": pandas.NA,
        "inconsistent_partial_pct": pandas.NA,
        "inconsistent_fragment_pct": pandas.NA,
        "contamination_pct": pandas.NA,
        "contamination_partial_pct": pandas.NA,
        "contamination_fragment_pct": pandas.NA,
        "unknown_pct": pandas.NA,
        "top_detected_clade": "-",
        "top_detected_taxid": pandas.NA,
        "top_detected_proteins": pandas.NA,
        "top_detected_pct": pandas.NA,
        "detected_clade_count": 0,
    }

def parse_sum_file(sum_path: Path) -> dict:
    row = _default_row(sum_path)
    lines = [line.rstrip("\n") for line in sum_path.read_text(encoding="utf-8").splitlines()]
    detected_rows = []
    in_detected_table = False

    for line in lines:
        if line.startswith("#The selected clade was "):
            row["selected_clade"] = line.replace("#The selected clade was ", "", 1).strip() or "-"
            continue
        if line.startswith("#Number of conserved HOGs is: "):
            row["conserved_hogs"] = int(line.rsplit(": ", 1)[1])
            continue
        if line.startswith("#On the whole proteome, there are "):
            match = re.search(r"(\d+)", line)
            if match:
                row["protein_count"] = int(match.group(1))
            continue
        if line.startswith("#Clade\tNCBI taxid\tNumber of associated proteins\tPercentage of proteome's total"):
            in_detected_table = True
            continue

        if in_detected_table:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 4:
                detected_rows.append(
                    {
                        "clade": fields[0],
                        "taxid": int(fields[1]) if fields[1].isdigit() else pandas.NA,
                        "proteins": int(fields[2]) if fields[2].isdigit() else pandas.NA,
                        "pct": float(fields[3].rstrip("%")) if fields[3].rstrip("%") else pandas.NA,
                    }
                )
            continue

        match = COMPLETENESS_COUNTS_RE.fullmatch(line)
        if match:
            for key, value in match.groupdict().items():
                row[f"{key}_count"] = int(value)
            continue

        match = COMPLETENESS_PCTS_RE.fullmatch(line)
        if match:
            for key, value in match.groupdict().items():
                row[f"{key}_pct"] = float(value)
            continue

        match = PROTEOME_COUNTS_RE.fullmatch(line)
        if match:
            for key, value in match.groupdict().items():
                row[f"{key}_count"] = int(value)
            continue

        match = PROTEOME_PCTS_RE.fullmatch(line)
        if match:
            for key, value in match.groupdict().items():
                row[f"{key}_pct"] = float(value)
            continue

    if detected_rows:
        detected_rows.sort(
            key=lambda entry: (
                entry["proteins"] if entry["proteins"] is not pandas.NA else -1,
                entry["clade"],
            ),
            reverse=True,
        )
        top = detected_rows[0]
        row["top_detected_clade"] = top["clade"]
        row["top_detected_taxid"] = top["taxid"]
        row["top_detected_proteins"] = top["proteins"]
        row["top_detected_pct"] = top["pct"]
        row["detected_clade_count"] = len(detected_rows)

    return row


def find_sum_files(omark_outdir: Path):
    if not omark_outdir.exists():
        return []
    return sorted(
        path
        for path in omark_outdir.rglob("*.sum")
        if path.is_file() and not any(part.startswith(".") for part in path.parts)
    )


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    omark_outdir = Path(args.omark_outdir).expanduser().resolve()
    outfile = Path(args.outfile).expanduser().resolve()
    rows = [parse_sum_file(sum_path) for sum_path in find_sum_files(omark_outdir)]
    df = pandas.DataFrame(rows, columns=HEADER)
    df = df.sort_values(["species", "sum_file"], kind="mergesort").reset_index(drop=True)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
