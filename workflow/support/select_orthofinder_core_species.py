#!/usr/bin/env python3

import argparse
import csv
from functools import cmp_to_key
import gzip
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Sequence, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from species_labeling import extract_species_label, matches_species_label


FASTA_EXTENSIONS = (
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".faa",
    ".fa.gz",
    ".fas.gz",
    ".fasta.gz",
    ".fna.gz",
    ".faa.gz",
)
DEFAULT_COLUMNS = ("leaf_name", "fasta_file", "fasta_basename", "num_seq", "busco_complete_pct")
NUMERIC_COLUMNS = {"num_seq", "busco_complete_pct"}
SUPPORTED_FILTER_OPS = {"ge", "gt", "le", "lt", "eq", "ne"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select size-aware OrthoFinder core species, using nwkit sample when "
            "a species tree is available."
        )
    )
    parser.add_argument("--protein-dir", required=True, type=Path)
    parser.add_argument("--busco-short-dir", required=True, type=Path)
    parser.add_argument("--max-core-species", required=True, type=int)
    parser.add_argument("--species-tree", default="", type=str)
    parser.add_argument("--filters", default="busco_complete_pct:ge:80,num_seq:le:100000")
    parser.add_argument("--rank", default="num_seq:asc,busco_complete_pct:desc")
    parser.add_argument("--method", default="max-pd")
    parser.add_argument("--candidates-table", required=True, type=Path)
    parser.add_argument("--selected-table", required=True, type=Path)
    parser.add_argument("--selected-list", required=True, type=Path)
    parser.add_argument("--core-tree", required=True, type=Path)
    return parser.parse_args()


def fasta_files(protein_dir: Path) -> List[Path]:
    if not protein_dir.is_dir():
        raise SystemExit(f"Protein FASTA directory was not found: {protein_dir}")
    files = [
        path
        for path in protein_dir.iterdir()
        if path.is_file() and not path.name.startswith(".") and path.name.endswith(FASTA_EXTENSIONS)
    ]
    return sorted(files, key=lambda path: path.name)


def leaf_name_from_fasta(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).stem


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def count_fasta_records(path: Path) -> int:
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def parse_busco_complete_pct(path: Path) -> str:
    if not path.is_file():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"(?:^|[,\s])C:([0-9]+(?:\.[0-9]+)?)%", text)
    if match:
        return match.group(1)
    return ""


def find_busco_short_file(busco_short_dir: Path, leaf_name: str) -> Optional[Path]:
    if not busco_short_dir.is_dir():
        return None
    leaf_label = extract_species_label(leaf_name, strip_extension=True) or leaf_name
    candidates = [
        busco_short_dir / f"{leaf_name}.busco.short.txt",
        busco_short_dir / f"{leaf_name}_busco.short.txt",
        busco_short_dir / f"{leaf_label}.busco.short.txt",
        busco_short_dir / f"{leaf_label}_busco.short.txt",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    matches = sorted(
        path
        for path in busco_short_dir.iterdir()
        if path.is_file()
        and not path.name.startswith(".")
        and path.name.endswith("busco.short.txt")
        and matches_species_label(path.name, leaf_label, strip_extension=True)
    )
    if len(matches) > 1:
        print(
            "Multiple BUSCO short files matched species label {}: {}".format(
                leaf_label,
                ", ".join(path.name for path in matches),
            ),
            file=sys.stderr,
        )
        return None
    return matches[0] if matches else None


def build_candidate_rows(protein_dir: Path, busco_short_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for fasta in fasta_files(protein_dir):
        leaf_name = leaf_name_from_fasta(fasta)
        busco_file = find_busco_short_file(busco_short_dir, leaf_name)
        rows.append(
            {
                "leaf_name": leaf_name,
                "fasta_file": str(fasta),
                "fasta_basename": fasta.name,
                "num_seq": str(count_fasta_records(fasta)),
                "busco_complete_pct": parse_busco_complete_pct(busco_file) if busco_file else "",
            }
        )
    return rows


def write_table(path: Path, rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(DEFAULT_COLUMNS)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_table(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def split_specs(specs: str) -> List[str]:
    return [part.strip() for part in specs.split(",") if part.strip()]


def parse_filter_spec(spec: str) -> Tuple[str, str, str]:
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise SystemExit(f"Invalid filter spec: {spec}. Expected COLUMN:OP:VALUE.")
    column, op, value = parts
    if op not in SUPPORTED_FILTER_OPS:
        raise SystemExit(f"Invalid filter operator in {spec}: {op}")
    return column, op, value


def parse_rank_spec(spec: str) -> Tuple[str, str]:
    parts = spec.split(":", 1)
    if len(parts) != 2:
        raise SystemExit(f"Invalid rank spec: {spec}. Expected COLUMN:asc|desc.")
    column, direction = parts[0], parts[1].lower()
    if direction not in {"asc", "desc"}:
        raise SystemExit(f"Invalid rank direction in {spec}: {direction}")
    return column, direction


def as_float(value: str) -> Optional[float]:
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def row_passes_filter(row: Dict[str, str], filter_spec: str) -> bool:
    column, op, value = parse_filter_spec(filter_spec)
    row_value = row.get(column, "")
    if column in NUMERIC_COLUMNS or op in {"ge", "gt", "le", "lt"}:
        left = as_float(row_value)
        right = as_float(value)
        if left is None or right is None:
            return False
        if op == "ge":
            return left >= right
        if op == "gt":
            return left > right
        if op == "le":
            return left <= right
        if op == "lt":
            return left < right
        if op == "eq":
            return left == right
        if op == "ne":
            return left != right
    if op == "eq":
        return row_value == value
    if op == "ne":
        return row_value != value
    return False


def filter_rows(rows: Sequence[Dict[str, str]], filter_specs: Sequence[str]) -> List[Dict[str, str]]:
    return [row for row in rows if all(row_passes_filter(row, spec) for spec in filter_specs)]


def choose_active_filters(rows: Sequence[Dict[str, str]], filter_specs: Sequence[str]) -> List[str]:
    filtered = filter_rows(rows, filter_specs)
    if filtered:
        return list(filter_specs)
    busco_filters = [spec for spec in filter_specs if parse_filter_spec(spec)[0] == "busco_complete_pct"]
    if not busco_filters:
        return list(filter_specs)
    if any(str(row.get("busco_complete_pct", "")).strip() for row in rows):
        return list(filter_specs)
    fallback = [spec for spec in filter_specs if parse_filter_spec(spec)[0] != "busco_complete_pct"]
    fallback_filtered = filter_rows(rows, fallback)
    if fallback_filtered:
        print(
            "BUSCO completeness values were not found for OrthoFinder core selection; "
            "retrying without BUSCO completeness filters.",
            file=sys.stderr,
        )
        return fallback
    return list(filter_specs)


def rank_rows(rows: Sequence[Dict[str, str]], rank_specs: Sequence[str]) -> List[Dict[str, str]]:
    parsed = [parse_rank_spec(spec) for spec in rank_specs]

    def compare(left: Dict[str, str], right: Dict[str, str]) -> int:
        for column, direction in parsed:
            left_value = left.get(column, "")
            right_value = right.get(column, "")
            left_num = as_float(left_value)
            right_num = as_float(right_value)
            left_missing = left_value == "" or (column in NUMERIC_COLUMNS and left_num is None)
            right_missing = right_value == "" or (column in NUMERIC_COLUMNS and right_num is None)
            if left_missing and right_missing:
                continue
            if left_missing:
                return 1
            if right_missing:
                return -1
            if left_num is not None and right_num is not None:
                if left_num == right_num:
                    continue
                result = -1 if left_num < right_num else 1
            else:
                if left_value == right_value:
                    continue
                result = -1 if left_value < right_value else 1
            return result if direction == "asc" else -result
        if left.get("leaf_name", "") == right.get("leaf_name", ""):
            return 0
        return -1 if left.get("leaf_name", "") < right.get("leaf_name", "") else 1

    return sorted(rows, key=cmp_to_key(compare))


def command_for_log(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(arg) for arg in command)


def run_nwkit_sample(
    species_tree: Path,
    candidates_table: Path,
    selected_table: Path,
    core_tree: Path,
    max_core_species: int,
    method: str,
    filter_specs: Sequence[str],
    rank_specs: Sequence[str],
) -> None:
    nwkit = shutil.which("nwkit")
    if nwkit is None:
        raise SystemExit(
            "nwkit was not found in PATH. Install a nwkit version with the sample command "
            "or run without a species tree to use ranked fallback selection."
        )
    command = [
        nwkit,
        "sample",
        "--infile",
        str(species_tree),
        "--trait",
        str(candidates_table),
        "--n",
        str(max_core_species),
        "--method",
        method,
        "--output-table",
        str(selected_table),
        "--outfile",
        str(core_tree),
        "--allow-fewer",
        "yes",
    ]
    for spec in filter_specs:
        command.extend(["--filter", spec])
    for spec in rank_specs:
        command.extend(["--rank", spec])
    print("Command:", command_for_log(command))
    try:
        completed = subprocess.run(command, check=False)
    except OSError as exc:
        raise SystemExit(f"Failed to run nwkit sample: {exc}") from exc
    if completed.returncode != 0:
        raise SystemExit(
            "nwkit sample failed with exit code "
            f"{completed.returncode}. GeneGalleon requires a nwkit version "
            "that provides the sample subcommand for tree-aware core selection."
        )


def write_selected_list(path: Path, selected_rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = [row.get("fasta_basename", "") for row in selected_rows if row.get("fasta_basename", "")]
    if not names:
        raise SystemExit("No OrthoFinder core species were selected.")
    path.write_text("".join(f"{name}\n" for name in names), encoding="utf-8")


def validate_selected_rows(selected_rows: Sequence[Dict[str, str]], candidate_rows: Sequence[Dict[str, str]]) -> None:
    candidate_names = {row["fasta_basename"] for row in candidate_rows}
    missing = [
        row.get("fasta_basename", "")
        for row in selected_rows
        if row.get("fasta_basename", "") not in candidate_names
    ]
    if missing:
        raise SystemExit(f"Selected FASTA basenames were not found among candidates: {', '.join(missing)}")


def main() -> int:
    args = parse_args()
    if args.max_core_species < 1:
        raise SystemExit("--max-core-species must be at least 1.")
    candidate_rows = build_candidate_rows(args.protein_dir, args.busco_short_dir)
    if not candidate_rows:
        raise SystemExit(f"No protein FASTA files were found in: {args.protein_dir}")
    write_table(args.candidates_table, candidate_rows)

    requested_filters = split_specs(args.filters)
    rank_specs = split_specs(args.rank)
    active_filters = choose_active_filters(candidate_rows, requested_filters)
    filtered_rows = filter_rows(candidate_rows, active_filters)
    if not filtered_rows:
        raise SystemExit(
            "No OrthoFinder core species candidates passed filters: "
            + (", ".join(active_filters) if active_filters else "(none)")
        )

    species_tree = Path(args.species_tree) if args.species_tree else None
    if species_tree is not None and species_tree.is_file():
        run_nwkit_sample(
            species_tree=species_tree,
            candidates_table=args.candidates_table,
            selected_table=args.selected_table,
            core_tree=args.core_tree,
            max_core_species=args.max_core_species,
            method=args.method,
            filter_specs=active_filters,
            rank_specs=rank_specs,
        )
        selected_rows = read_table(args.selected_table)
        if not selected_rows:
            raise SystemExit("nwkit sample did not select any OrthoFinder core species.")
    else:
        print("No species tree was provided; using ranked OrthoFinder core-species selection.")
        selected_rows = rank_rows(filtered_rows, rank_specs)[: args.max_core_species]
        write_table(args.selected_table, selected_rows)
        if args.core_tree.exists():
            args.core_tree.unlink()

    validate_selected_rows(selected_rows, candidate_rows)
    write_selected_list(args.selected_list, selected_rows)
    print(f"Selected {len(selected_rows)} OrthoFinder core species.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
