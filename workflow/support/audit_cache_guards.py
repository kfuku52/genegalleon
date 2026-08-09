#!/usr/bin/env python3
"""Inventory shell and Python cache guards lacking content/parameter provenance."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

CONDITION_START = re.compile(r"^\s*(?:if|elif)\s+\[\[")
OUTPUT_GUARD = re.compile(
    r"!\s+-(?:s|f|d|e)\s+\"?\$\{([A-Za-z_][A-Za-z0-9_]*)\}\"?"
)
ANY_EXISTENCE_GUARD = re.compile(
    r"-(?:s|f|d|e)\s+\"?\$\{([A-Za-z_][A-Za-z0-9_]*)\}\"?"
)
MTIME_VARIABLE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}[^\n]*(?:-nt|-ot)|(?:-nt|-ot)[^\n]*\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
COUNT_COMPARISON = re.compile(
    r"\$\{(num_[A-Za-z0-9_]+)\}\s+-(?:eq|ne|lt|le|gt|ge)\s+\$\{(num_[A-Za-z0-9_]+)\}"
)
READINESS_VARIABLE = re.compile(
    r"\$\{([A-Za-z_][A-Za-z0-9_]*(?:_ready|_complete|_current))\}"
)
MTIME_HELPER = re.compile(r"\bis_output_older_than_inputs\b")
RUN_FLAG = re.compile(
    r"\$\{run_[A-Za-z0-9_]+\}\s+-(?:eq|ne|lt|le|gt|ge)\s+[01](?:\s|$)"
)
TASK_ASSIGNMENT = re.compile(r'^\s*task=(?:"([^"]+)"|\'([^\']+)\')')
EARLY_EXIT = re.compile(r"^\s*(?:return|continue)(?:\s|$)", re.M)
AUDITED_MARKER = "gg-cache-guard: audited"
OUTPUT_LIKE_VARIABLE = re.compile(
    r"(?:^|_)(?:out(?:file)?[0-9]*|output|result|summary|cache|stamp|ready)(?:_|$)"
)


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
    for line in lines:
        task_match = TASK_ASSIGNMENT.match(line)
        if task_match:
            current_task = task_match.group(1) or task_match.group(2) or ""
        task_by_line.append(current_task)

    relative_file = path.relative_to(source_root).as_posix()
    findings: list[Finding] = []
    for start, condition in iter_conditions(lines):
        if AUDITED_MARKER in condition or (
            start > 0 and AUDITED_MARKER in lines[start - 1]
        ):
            continue
        output_variables = set(OUTPUT_GUARD.findall(condition))
        all_existence_variables = set(ANY_EXISTENCE_GUARD.findall(condition))
        positive_existence_variables = all_existence_variables - output_variables
        mtime_matches = MTIME_VARIABLE.findall(condition)
        for left, right in mtime_matches:
            output_variables.add(left or right)
        count_matches = COUNT_COMPARISON.findall(condition)
        for left, right in count_matches:
            output_variables.update((left, right))
        readiness_variables = READINESS_VARIABLE.findall(condition)
        output_variables.update(readiness_variables)
        nearby_lines = "\n".join(lines[start : min(len(lines), start + 12)])
        early_guard_lines = "\n".join(lines[start : min(len(lines), start + 5)])
        sets_manual_update_flag = re.search(
            r"[A-Za-z_][A-Za-z0-9_]*_needs_update=1", nearby_lines
        )
        has_run_flag = RUN_FLAG.search(condition) is not None
        has_existence_guard = OUTPUT_GUARD.search(condition) is not None
        has_mtime_guard = bool(mtime_matches) or MTIME_HELPER.search(condition) is not None
        has_count_guard = bool(count_matches)
        has_readiness_guard = bool(readiness_variables)
        early_reuse_variables = {
            variable
            for variable in positive_existence_variables
            if OUTPUT_LIKE_VARIABLE.search(variable)
        }
        has_early_reuse_guard = (
            bool(task_by_line[start])
            and bool(early_reuse_variables)
            and bool(EARLY_EXIT.search(early_guard_lines))
        )
        if not (
            has_existence_guard
            or has_mtime_guard
            or has_count_guard
            or has_readiness_guard
            or has_early_reuse_guard
        ):
            continue
        # mtime-based freshness is always suspicious, including when a run
        # flag is supplied by an enclosing shell condition.
        if (
            not has_mtime_guard
            and not has_early_reuse_guard
            and not has_run_flag
            and sets_manual_update_flag is None
        ):
            continue
        if has_mtime_guard:
            guard_type = "mtime_only"
        elif has_count_guard:
            guard_type = "count_only"
        elif has_readiness_guard:
            guard_type = "readiness_only"
        elif has_early_reuse_guard:
            guard_type = "early_exit_existence_only"
            output_variables.update(early_reuse_variables)
        else:
            guard_type = "existence_only"
        task = task_by_line[start]
        findings.append(
            Finding(
                fingerprint=finding_fingerprint(relative_file, task, condition),
                file=relative_file,
                line=start + 1,
                task=task,
                guard_type=guard_type,
                output_variables=",".join(sorted(output_variables)),
                condition=normalize_condition(condition),
            )
        )
    return findings


def python_test_has_positive_existence_call(node: ast.AST, negated: bool = False) -> bool:
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return python_test_has_positive_existence_call(node.operand, not negated)
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in {"exists", "is_file", "is_dir"}
        and not negated
    ):
        return True
    return any(
        python_test_has_positive_existence_call(child, negated)
        for child in ast.iter_child_nodes(node)
    )


def python_body_has_early_exit(statements: list[ast.stmt]) -> bool:
    def statement_has_early_exit(statement: ast.stmt) -> bool:
        if isinstance(statement, (ast.Return, ast.Continue)):
            return True
        if isinstance(statement, ast.If):
            return python_body_has_early_exit(statement.body) or python_body_has_early_exit(
                statement.orelse
            )
        return False

    for statement in statements[:3]:
        if statement_has_early_exit(statement):
            return True
    return False


def scan_python(path: Path, source_root: Path) -> list[Finding]:
    source = path.read_text(encoding="utf-8")
    lines = source.splitlines()
    tree = ast.parse(source, filename=str(path))
    relative_file = path.resolve().relative_to(source_root).as_posix()
    functions = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]
    findings: list[Finding] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        line_index = node.lineno - 1
        if AUDITED_MARKER in lines[line_index] or (
            line_index > 0 and AUDITED_MARKER in lines[line_index - 1]
        ):
            continue
        if not python_test_has_positive_existence_call(node.test):
            continue
        if not python_body_has_early_exit(node.body):
            continue
        names = {
            candidate.id
            for candidate in ast.walk(node.test)
            if isinstance(candidate, ast.Name)
            and OUTPUT_LIKE_VARIABLE.search(candidate.id)
        }
        names.update(
            candidate.attr
            for candidate in ast.walk(node.test)
            if isinstance(candidate, ast.Attribute)
            and OUTPUT_LIKE_VARIABLE.search(candidate.attr)
        )
        if not names:
            continue
        enclosing_functions = [
            function
            for function in functions
            if function.lineno <= node.lineno <= (function.end_lineno or function.lineno)
        ]
        task = ""
        if enclosing_functions:
            task = min(
                enclosing_functions,
                key=lambda function: (function.end_lineno or function.lineno) - function.lineno,
            ).name
        condition = normalize_condition(ast.unparse(node.test))
        findings.append(
            Finding(
                fingerprint=finding_fingerprint(relative_file, task, condition),
                file=relative_file,
                line=node.lineno,
                task=task,
                guard_type="python_early_exit_existence_only",
                output_variables=",".join(sorted(names)),
                condition=condition,
            )
        )
    return findings


def scan_paths(paths: list[Path], source_root: Path) -> list[Finding]:
    scripts: set[Path] = set()
    python_files: set[Path] = set()
    for path in paths:
        if path.is_dir():
            scripts.update(path.rglob("*.sh"))
            python_files.update(path.rglob("*.py"))
        elif path.suffix == ".sh":
            scripts.add(path)
        elif path.suffix == ".py":
            python_files.add(path)
    findings = [
        finding
        for script in sorted(scripts)
        for finding in scan_script(script.resolve(), source_root)
    ]
    findings.extend(
        finding
        for python_file in sorted(python_files)
        for finding in scan_python(python_file.resolve(), source_root)
    )
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
    text = "# Known workflow cache guards awaiting provenance migration.\n"
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
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[default_root / "core", default_root / "support"],
    )
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
        "count_only": sum(finding.guard_type == "count_only" for finding in findings),
        "readiness_only": sum(finding.guard_type == "readiness_only" for finding in findings),
        "early_exit_existence_only": sum(
            finding.guard_type == "early_exit_existence_only" for finding in findings
        ),
        "python_early_exit_existence_only": sum(
            finding.guard_type == "python_early_exit_existence_only"
            for finding in findings
        ),
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
