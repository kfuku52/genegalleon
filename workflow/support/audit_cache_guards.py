#!/usr/bin/env python3
"""Inventory shell cache guards that lack content-and-parameter provenance."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

CONDITION_START = re.compile(r"^\s*(?:if|elif)\s+\[\[")
OUTPUT_MISSING = re.compile(r"!\s+-s\s+\"?\$\{([A-Za-z_][A-Za-z0-9_]*)\}\"?")
RUN_FLAG = re.compile(r"\$\{run_[A-Za-z0-9_]+\}")
TASK_ASSIGNMENT = re.compile(r'^\s*task=(?:"([^"]+)"|\'([^\']+)\')')
PROVENANCE_CALL = "gg_artifact_needs_run"


@dataclass(frozen=True)
class Finding:
    fingerprint: str
    file: str
    line: int
    task: str
    guard_type: str
    output_variables: str
    condition: str


def normalize_condition(condition: str) -> str:
    return " ".join(condition.split())


def finding_fingerprint(relative_file: str, task: str, condition: str) -> str:
    payload = f"{relative_file}\n{task}\n{normalize_condition(condition)}\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def iter_conditions(lines: list[str]):
    index = 0
    while index < len(lines):
        if not CONDITION_START.search(lines[index]):
            index += 1
            continue
        start = index
        parts = [lines[index].strip()]
        while "then" not in lines[index] and index + 1 < len(lines):
            index += 1
            parts.append(lines[index].strip())
        yield start, " ".join(parts)
        index += 1


def scan_script(path: Path, source_root: Path) -> list[Finding]:
    lines = path.read_text(encoding="utf-8").splitlines()
    task_by_line: list[str] = []
    current_task = ""
    provenance_seen_in_task = False
    provenance_by_line: list[bool] = []
    for line in lines:
        task_match = TASK_ASSIGNMENT.match(line)
        if task_match:
            current_task = task_match.group(1) or task_match.group(2) or ""
            provenance_seen_in_task = False
        if PROVENANCE_CALL in line:
            provenance_seen_in_task = True
        task_by_line.append(current_task)
        provenance_by_line.append(provenance_seen_in_task)

    relative_file = path.relative_to(source_root).as_posix()
    findings: list[Finding] = []
    for start, condition in iter_conditions(lines):
        output_variables = sorted(set(OUTPUT_MISSING.findall(condition)))
        nearby_lines = "\n".join(lines[start : min(len(lines), start + 12)])
        sets_manual_update_flag = re.search(
            r"[A-Za-z_][A-Za-z0-9_]*_needs_update=1", nearby_lines
        )
        if not output_variables or (
            RUN_FLAG.search(condition) is None and sets_manual_update_flag is None
        ):
            continue
        if provenance_by_line[start] or "_needs_update" in condition:
            continue
        guard_type = "mtime_only" if " -nt " in condition or " -ot " in condition else "existence_only"
        task = task_by_line[start]
        findings.append(
            Finding(
                fingerprint=finding_fingerprint(relative_file, task, condition),
                file=relative_file,
                line=start + 1,
                task=task,
                guard_type=guard_type,
                output_variables=",".join(output_variables),
                condition=normalize_condition(condition),
            )
        )
    return findings


def scan_paths(paths: list[Path], source_root: Path) -> list[Finding]:
    scripts: set[Path] = set()
    for path in paths:
        if path.is_dir():
            scripts.update(path.rglob("*.sh"))
        elif path.suffix == ".sh":
            scripts.add(path)
    findings = [finding for script in sorted(scripts) for finding in scan_script(script, source_root)]
    return sorted(findings, key=lambda finding: (finding.file, finding.line))


def load_baseline(path: Path) -> set[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }


def write_baseline(path: Path, findings: list[Finding]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "# Known shell cache guards awaiting provenance migration.\n"
    text += "\n".join(finding.fingerprint for finding in findings) + "\n"
    path.write_text(text, encoding="utf-8")


def write_tsv(path: Path | None, findings: list[Finding]) -> None:
    header = "fingerprint\tfile\tline\ttask\tguard_type\toutput_variables\tcondition\n"
    rows = "".join(
        "\t".join(
            str(value).replace("\t", " ").replace("\n", " ")
            for value in (
                finding.fingerprint,
                finding.file,
                finding.line,
                finding.task,
                finding.guard_type,
                finding.output_variables,
                finding.condition,
            )
        )
        + "\n"
        for finding in findings
    )
    if path is None:
        sys.stdout.write(header + rows)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(header + rows, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[1]
    parser.add_argument("paths", nargs="*", type=Path, default=[default_root / "core"])
    parser.add_argument("--source-root", type=Path, default=default_root.parent)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--write-baseline", type=Path)
    parser.add_argument("--fail-on-findings", action="store_true")
    parser.add_argument("--json-summary", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = scan_paths(args.paths, args.source_root.resolve())
    if args.write_baseline is not None:
        write_baseline(args.write_baseline, findings)
    write_tsv(args.output_tsv, findings)

    unapproved = findings
    if args.baseline is not None:
        baseline = load_baseline(args.baseline)
        unapproved = [finding for finding in findings if finding.fingerprint not in baseline]
    summary = {
        "findings": len(findings),
        "unapproved_findings": len(unapproved),
        "existence_only": sum(finding.guard_type == "existence_only" for finding in findings),
        "mtime_only": sum(finding.guard_type == "mtime_only" for finding in findings),
    }
    if args.json_summary:
        print(json.dumps(summary, sort_keys=True), file=sys.stderr)
    else:
        print(
            "Cache guard audit: " + ", ".join(f"{key}={value}" for key, value in summary.items()),
            file=sys.stderr,
        )
    if args.baseline is not None and unapproved:
        for finding in unapproved:
            print(
                f"New unprovenanced cache guard: {finding.file}:{finding.line}: {finding.condition}",
                file=sys.stderr,
            )
        return 1
    if args.fail_on_findings and findings:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
