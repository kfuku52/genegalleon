#!/usr/bin/env python3
"""Extract all BUSCO marker FASTAs while streaming species inputs once."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Iterator, TextIO

from busco_table_utils import BUSCO_ID_PATTERN


class ExtractionError(RuntimeError):
    pass


def validate_busco_id(busco_id: str, line_number: int) -> None:
    if busco_id in {".", ".."} or BUSCO_ID_PATTERN.fullmatch(busco_id) is None:
        raise ExtractionError(f"Invalid BUSCO identifier on row {line_number}: {busco_id!r}")


def validate_suffix(suffix: str) -> None:
    if re.fullmatch(r"\.[A-Za-z0-9][A-Za-z0-9._-]*", suffix) is None:
        raise ExtractionError(f"Invalid output suffix: {suffix!r}")


def safe_output_path(output_dir: Path, busco_id: str, suffix: str) -> Path:
    path = output_dir / f"{busco_id}{suffix}"
    if path.parent.resolve() != output_dir.resolve():
        raise ExtractionError(f"BUSCO output escapes its output directory: {path}")
    if path.is_symlink():
        raise ExtractionError(f"Refusing to follow BUSCO output symlink: {path}")
    return path


def open_output_append(path: Path) -> TextIO:
    flags = os.O_WRONLY | os.O_CREAT | os.O_APPEND
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o666)
    return os.fdopen(descriptor, "a", encoding="utf-8")


def open_output_replace(path: Path) -> TextIO:
    flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o666)
    return os.fdopen(descriptor, "w", encoding="utf-8")


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, list[str]]]:
    record: list[str] = []
    identifier = ""
    has_sequence = False
    with open_text(path) as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                if record:
                    if not has_sequence:
                        raise ExtractionError(f"Empty FASTA sequence for {identifier!r}: {path}")
                    yield identifier, record
                fields = raw_line[1:].split(None, 1)
                if not fields:
                    raise ExtractionError(f"Empty FASTA identifier: {path}")
                identifier = fields[0]
                record = [raw_line if raw_line.endswith("\n") else raw_line + "\n"]
                has_sequence = False
            elif record:
                record.append(raw_line if raw_line.endswith("\n") else raw_line + "\n")
                has_sequence = has_sequence or bool(raw_line.strip())
            elif raw_line.strip():
                raise ExtractionError(f"Sequence data precedes the first FASTA header: {path}")
    if record:
        if not has_sequence:
            raise ExtractionError(f"Empty FASTA sequence for {identifier!r}: {path}")
        yield identifier, record


class BestScoreSelection:
    """Opt-in representative selection; HMM ranking does not establish orthology."""

    def __init__(self, full_dir: Path):
        self.full_dir = full_dir
        self.tables = {}
        self.table_paths = []
        self.audit = []

    def choose(self, species: str, busco: str, token: str) -> str:
        if BUSCO_ID_PATTERN.fullmatch(species) is None or species in {".", ".."}:
            raise ExtractionError(f"Invalid species column: {species!r}")
        if species not in self.tables:
            path = self.full_dir / (species + ".busco.full.tsv")
            if path.is_symlink() or not path.is_file():
                raise ExtractionError(f"Missing regular BUSCO full table: {path}")
            table = {}
            with path.open() as stream:
                for line in stream:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) >= 2 and fields[1] == "Missing":
                        continue
                    if len(fields) < 5:
                        raise ExtractionError(f"Incomplete scored BUSCO row in {path}")
                    if fields[1] not in {"Complete", "Duplicated", "Fragmented"}:
                        raise ExtractionError(f"Unsupported scored BUSCO status in {path}: {fields[1]}")
                    try:
                        score, length = float(fields[3]), int(fields[4])
                    except ValueError as exc:
                        raise ExtractionError(f"Invalid BUSCO score/length in {path}") from exc
                    if not math.isfinite(score) or length <= 0:
                        raise ExtractionError(f"Invalid BUSCO score/length in {path}")
                    key = (fields[0], fields[2])
                    if key in table:
                        raise ExtractionError(f"Ambiguous duplicate BUSCO scoring row: {key}")
                    table[key] = (score, length, fields[1])
            self.tables[species] = table
            self.table_paths.append(path)
        candidates = token.split(",")
        if len(candidates) != len(set(candidates)) or any(not gene for gene in candidates):
            raise ExtractionError(f"Invalid BUSCO candidates: {busco}/{species}")
        scored = []
        for gene in candidates:
            row = self.tables[species].get((busco, gene))
            if row is None:
                raise ExtractionError(f"Missing BUSCO score: {busco}/{species}/{gene}")
            score, length, status = row
            scored.append((gene, score, length, status))
        ranked = sorted(scored, key=lambda row: (-row[1], -row[2], row[0]))
        selected = ranked[0][0]
        for gene, score, length, status in ranked:
            self.audit.append((busco, species, gene, score, length, status, int(gene == selected)))
        return selected

    def write(self, path: Path) -> None:
        if path.is_symlink():
            raise ExtractionError(f"Refusing BUSCO selection report symlink: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with open_output_replace(path) as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["busco_id", "species", "gene_id", "hmm_score", "alignment_length", "status", "selected"])
            writer.writerows(self.audit)


def read_rows(summary: Path, mode: str, strict: bool, selection=None) -> tuple[dict[str, list[str]], dict[str, int]]:
    gene_to_buscos: dict[str, list[str]] = defaultdict(list)
    expected: dict[str, int] = {}
    with summary.open(encoding="utf-8") as handle:
        header = next(handle, "").rstrip("\r\n").split("\t")
        species = header[3:]
        if selection is not None and (not species or len(species) != len(set(species))):
            raise ExtractionError("BUSCO summary requires unique species columns")
        for line_number, raw_line in enumerate(handle, start=2):
            fields = raw_line.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                continue
            busco_id = fields[0]
            validate_busco_id(busco_id, line_number)
            tokens = fields[3:] if len(fields) > 3 else []
            if selection is not None and len(tokens) != len(species):
                raise ExtractionError(f"BUSCO summary row width mismatch: {busco_id}")
            if mode == "single-copy":
                if strict and any(token == "-" or "," in token for token in tokens):
                    continue
                genes = ([selection.choose(sp, busco_id, token) for sp, token in zip(species, tokens, strict=True)
                          if token != "-"] if selection is not None else
                         [token for token in tokens if token != "-" and "," not in token])
            else:
                genes = []
                for token in tokens:
                    if token == "-":
                        continue
                    genes.extend(gene for gene in token.split(",") if gene and gene != "-")
            if not genes:
                continue
            if len(genes) != len(set(genes)):
                raise ExtractionError(f"BUSCO row {line_number} contains duplicate gene identifiers: {busco_id}")
            expected[busco_id] = len(genes)
            for gene in genes:
                gene_to_buscos[gene].append(busco_id)
    return gene_to_buscos, expected


class OutputPool:
    def __init__(self, output_dir: Path, suffix: str, max_open: int):
        self.output_dir = output_dir
        self.suffix = suffix
        self.max_open = max(1, max_open)
        self.handles: OrderedDict[str, TextIO] = OrderedDict()
        self.paths: dict[str, Path] = {}

    def write(self, busco_id: str, record: list[str]) -> None:
        handle = self.handles.pop(busco_id, None)
        if handle is None:
            if len(self.handles) >= self.max_open:
                _old_id, old_handle = self.handles.popitem(last=False)
                old_handle.close()
            if busco_id not in self.paths:
                self.paths[busco_id] = safe_output_path(self.output_dir, busco_id, self.suffix)
            path = self.paths[busco_id]
            # O_NOFOLLOW in open_output_append still checks every reopen.
            handle = open_output_append(path)
        self.handles[busco_id] = handle
        handle.writelines(record)

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def validate_output_layout(args, fasta_paths, expected, selection):
    """Check all destinations before clearing outputs or truncating a report."""
    protected = [args.summary, args.fasta_list, *fasta_paths]
    if selection is not None:
        protected += selection.table_paths
        if args.species_fasta_map:
            protected.append(args.species_fasta_map)
        if args.species_genetic_codes:
            protected.append(args.species_genetic_codes)
    sources = {path.resolve() for path in protected}
    source_inodes = {(path.stat().st_dev, path.stat().st_ino) for path in protected}
    root = args.output_dir.resolve()
    if any(root == path or root in path.parents for path in sources):
        raise ExtractionError("Output paths overlap input files: output directory contains inputs")
    if selection is not None and root.is_relative_to(args.busco_full_dir.resolve()):
        raise ExtractionError("Output paths overlap input BUSCO full-table directory")
    reports = [args.report] + ([args.selection_report] if args.selection_report else [])
    output_paths = {safe_output_path(args.output_dir, marker, args.suffix).resolve() for marker in expected}
    output_inodes = {(path.stat().st_dev, path.stat().st_ino) for path in output_paths if path.exists()}
    seen = set()
    seen_inodes = set()
    for report in reports:
        resolved = report.resolve()
        inode = (report.stat().st_dev, report.stat().st_ino) if report.exists() else None
        if (resolved in sources or resolved in output_paths or resolved in seen or
                inode is not None and inode in source_inodes | output_inodes | seen_inodes):
            raise ExtractionError(f"Output paths overlap input or other output files: {report}")
        if selection is not None and resolved.is_relative_to(args.busco_full_dir.resolve()):
            raise ExtractionError(f"Output paths overlap input BUSCO full-table directory: {report}")
        if report.is_symlink():
            raise ExtractionError(f"Refusing to replace BUSCO report symlink: {report}")
        if report.exists() and not report.is_file():
            raise ExtractionError(f"Unexpected non-file BUSCO report path: {report}")
        if any(parent.exists() and not parent.is_dir() for parent in report.parents):
            raise ExtractionError(f"BUSCO report parent is not a directory: {report}")
        seen.add(resolved)
        if inode is not None:
            seen_inodes.add(inode)


def run(args: argparse.Namespace) -> int:
    validate_suffix(args.suffix)
    if args.copy_policy == "best-score" and args.duplicate_policy == "best_score":
        raise ExtractionError("--copy-policy best-score conflicts with --duplicate-policy best_score")
    if args.species_headers and args.duplicate_policy != "best_score":
        raise ExtractionError("--species-headers requires best_score and its explicit species FASTA map")
    if args.duplicate_policy == "best_score":
        if args.mode != "single-copy" or args.strict:
            raise ExtractionError("best_score requires single-copy mode without --strict")
        if not args.busco_full_dir or not args.species_fasta_map or not args.selection_report:
            raise ExtractionError("best_score requires --busco-full-dir, --species-fasta-map and --selection-report")
    fasta_paths = [
        Path(line.strip()) for line in args.fasta_list.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if len(fasta_paths) != len(set(fasta_paths)):
        raise ExtractionError("Repeated FASTA input paths")
    for path in fasta_paths:
        if not path.is_file() or path.is_symlink():
            raise ExtractionError(f"Invalid FASTA input: {path}")
    selection = None
    legacy_selection = None
    if args.duplicate_policy == "best_score":
        from busco_copy_selection import BestScoreSelection as QualitySelection

        selection = QualitySelection(args, fasta_paths)
        expected = selection.expected
        for line_number, marker in enumerate(expected, start=2):
            validate_busco_id(marker, line_number)
    else:
        if args.copy_policy == "best-score":
            if args.mode != "single-copy" or args.strict or not args.busco_full_dir or not args.selection_report:
                raise ExtractionError("best-score requires single-copy mode, full tables and selection report, without --strict")
            legacy_selection = BestScoreSelection(args.busco_full_dir)
        gene_to_buscos, expected = read_rows(args.summary, args.mode, args.strict, legacy_selection)
    if args.output_dir.is_symlink():
        raise ExtractionError(f"Refusing to use a symlink as the BUSCO output directory: {args.output_dir}")
    validate_output_layout(args, fasta_paths, expected, selection or legacy_selection)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path in args.output_dir.glob(f"*{args.suffix}"):
        if path.is_symlink():
            raise ExtractionError(f"Refusing to replace BUSCO output symlink: {path}")
        if not path.is_file():
            raise ExtractionError(f"Unexpected non-file BUSCO output path: {path}")
        path.unlink()
    found = defaultdict(int)
    pool = OutputPool(args.output_dir, args.suffix, args.max_open_files)
    try:
        if selection is not None:
            found = selection.extract(fasta_paths, fasta_records, pool)
        else:
            for fasta_path in fasta_paths:
                for identifier, record in fasta_records(fasta_path):
                    for busco_id in gene_to_buscos.get(identifier, ()):
                        pool.write(busco_id, record)
                        found[busco_id] += 1
    finally:
        pool.close()
    if args.selection_report and legacy_selection is None:
        from busco_copy_selection import write_selection_report

        args.selection_report.parent.mkdir(parents=True, exist_ok=True)
        with open_output_replace(args.selection_report) as handle:
            write_selection_report(handle, selection.report_rows if selection else [])
    missing = [busco_id for busco_id, count in expected.items() if found.get(busco_id, 0) != count]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    if args.report.is_symlink():
        raise ExtractionError(f"Refusing to replace BUSCO report symlink: {args.report}")
    with open_output_replace(args.report) as handle:
        handle.write("busco_id\texpected\tfound\toutput\n")
        for busco_id in expected:
            output = safe_output_path(args.output_dir, busco_id, args.suffix)
            handle.write(f"{busco_id}\t{expected[busco_id]}\t{found.get(busco_id, 0)}\t{output}\n")
    print(
        f"Extracted {len(expected) - len(missing)} complete BUSCO marker FASTAs "
        f"from {len(fasta_paths)} species FASTA file(s) in one pass."
    )
    if args.require_complete and missing:
        print(
            "BUSCO sequence counts did not match the summary: " + ", ".join(missing[:20]),
            file=sys.stderr,
        )
        return 3
    if legacy_selection is not None:
        legacy_selection.write(args.selection_report)
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--summary", required=True, type=Path)
    result.add_argument("--fasta-list", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--suffix", required=True)
    result.add_argument("--report", required=True, type=Path)
    result.add_argument("--mode", choices=("single-copy", "duplicate-aware"), required=True)
    result.add_argument("--copy-policy", choices=("single-copy", "best-score"), default="single-copy",
                        help="Opt-in best-score: one candidate per species/BUSCO, ranked by HMM score, aligned length, then gene ID; not proof of orthology")
    result.add_argument("--strict", action="store_true")
    result.add_argument("--require-complete", action="store_true")
    result.add_argument("--max-open-files", type=int, default=64)
    result.add_argument("--duplicate-policy", choices=("exclude", "best_score"), default="exclude")
    result.add_argument("--busco-full-dir", type=Path)
    result.add_argument("--species-fasta-map", type=Path, help="TSV with species and fasta columns")
    result.add_argument("--selection-report", type=Path, help="Persistent per-candidate duplicate selection TSV")
    result.add_argument("--sequence-mode", choices=("cds", "protein"), default="protein")
    result.add_argument("--genetic-code", type=int, default=1)
    result.add_argument("--species-genetic-codes", type=Path)
    result.add_argument("--species-headers", action="store_true", help="Use explicitly mapped species as output FASTA IDs")
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        return run(args)
    except (OSError, UnicodeError, ExtractionError, ValueError) as exc:
        print(f"BUSCO batch extraction error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
