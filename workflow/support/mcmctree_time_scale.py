#!/usr/bin/env python3
"""Scale MCMCTree time values between public and internal units."""

from __future__ import annotations

import argparse
import math
import re
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path

NUMERIC_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
CALIBRATION_RE = re.compile(r"(?P<kind>\b[BUL])\((?P<body>[^)]*)\)")
BRANCH_LENGTH_RE = re.compile(
    rf"(?P<prefix>:\s*)(?P<number>{NUMERIC_RE})(?=\s*(?:\[|,|\)|;))"
)
HPD_COMMENT_RE = re.compile(r"\[&[^\]]*(?:95%|HPD)[^\]]*\]")
BRACE_RE = re.compile(r"\{(?P<body>[^{}]*)\}")
INTERNAL_NUMERIC_LABEL_RE = re.compile(
    rf"(?P<prefix>\))(?P<label>'?\s*{NUMERIC_RE}(?:\s*,\s*{NUMERIC_RE})?\s*'?)(?=\s*(?::|,|\)|;|\[))"
)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def parse_decimal(text: str) -> Decimal | None:
    stripped = text.strip()
    if not stripped or stripped == "-":
        return None
    try:
        return Decimal(stripped)
    except InvalidOperation:
        return None


def format_decimal(value: Decimal) -> str:
    if value.is_zero():
        return "0"
    normalized = value.normalize()
    text = format(normalized, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    if text == "-0":
        return "0"
    return text


def scale_decimal(value: Decimal, scale: Decimal, direction: str) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 28
        if direction == "down":
            return value / scale
        if direction == "up":
            return value * scale
    raise ValueError(f"invalid direction: {direction}")


def scale_number_text(text: str, scale: Decimal, direction: str) -> str:
    value = parse_decimal(text)
    if value is None:
        return text
    return format_decimal(scale_decimal(value, scale, direction))


def calibration_age_field_count(kind: str, field_count: int) -> int:
    if kind == "B":
        return min(2, field_count)
    if kind in {"L", "U"}:
        return min(1, field_count)
    return 0


def split_field_spacing(field: str) -> tuple[str, str, str]:
    leading_len = len(field) - len(field.lstrip())
    trailing_len = len(field) - len(field.rstrip())
    leading = field[:leading_len]
    trailing = field[len(field) - trailing_len :] if trailing_len else ""
    core_end = len(field) - trailing_len if trailing_len else len(field)
    return leading, field[leading_len:core_end], trailing


def calibration_ages(text: str) -> list[Decimal]:
    ages: list[Decimal] = []
    for match in CALIBRATION_RE.finditer(text):
        kind = match.group("kind")
        fields = match.group("body").split(",")
        for idx in range(calibration_age_field_count(kind, len(fields))):
            value = parse_decimal(fields[idx])
            if value is not None:
                ages.append(value)
    return ages


def choose_scale_factor(text: str, target_max: Decimal) -> Decimal:
    ages = calibration_ages(text)
    if not ages:
        return Decimal(1)
    max_age = max(abs(age) for age in ages)
    factor = Decimal(1)
    while max_age / factor > target_max:
        factor *= Decimal(10)
    return factor


def scale_calibration_labels(text: str, scale: Decimal, direction: str) -> str:
    if scale == 1:
        return text

    def replace(match: re.Match[str]) -> str:
        kind = match.group("kind")
        fields = match.group("body").split(",")
        for idx in range(calibration_age_field_count(kind, len(fields))):
            leading, core, trailing = split_field_spacing(fields[idx])
            fields[idx] = (
                leading + scale_number_text(core, scale, direction) + trailing
            )
        return f"{kind}({','.join(fields)})"

    return CALIBRATION_RE.sub(replace, text)


def scale_hpd_comment(comment: str, scale: Decimal, direction: str) -> str:
    def replace_brace(match: re.Match[str]) -> str:
        body = match.group("body")
        scaled = re.sub(
            NUMERIC_RE,
            lambda num_match: scale_number_text(num_match.group(0), scale, direction),
            body,
        )
        return "{" + scaled + "}"

    return BRACE_RE.sub(replace_brace, comment)


def scale_newick_time_values(text: str, scale: Decimal, direction: str) -> str:
    if scale == 1:
        return text

    scaled = BRANCH_LENGTH_RE.sub(
        lambda match: match.group("prefix")
        + scale_number_text(match.group("number"), scale, direction),
        text,
    )
    scaled = HPD_COMMENT_RE.sub(
        lambda match: scale_hpd_comment(match.group(0), scale, direction),
        scaled,
    )

    def replace_internal_label(match: re.Match[str]) -> str:
        label = match.group("label")
        scaled_label = re.sub(
            NUMERIC_RE,
            lambda num_match: scale_number_text(num_match.group(0), scale, direction),
            label,
        )
        return match.group("prefix") + scaled_label

    return INTERNAL_NUMERIC_LABEL_RE.sub(replace_internal_label, scaled)


def looks_like_figtree_tree_line(line: str) -> bool:
    stripped = line.strip()
    if "UTREE" in line:
        return True
    return stripped.startswith("(") and stripped.endswith(";") and ")" in stripped


def scale_figtree_text(text: str, scale: Decimal, direction: str) -> str:
    lines = text.splitlines(keepends=True)
    scaled_lines = []
    for line in lines:
        if looks_like_figtree_tree_line(line) and figtree_tree_kind(line) == "dated":
            scaled_lines.append(scale_newick_time_values(line, scale, direction))
        else:
            scaled_lines.append(line)
    return "".join(scaled_lines)


def scale_rootage_line(line: str, scale: Decimal, direction: str) -> str:
    if not re.match(r"^\s*RootAge\s*=", line):
        return line
    if "#" in line:
        value_part, comment = line.split("#", 1)
        comment = "#" + comment
    else:
        value_part = line
        comment = ""
    if CALIBRATION_RE.search(value_part):
        return scale_calibration_labels(value_part, scale, direction) + comment
    return re.sub(
        NUMERIC_RE,
        lambda match: scale_number_text(match.group(0), scale, direction),
        value_part,
    ) + comment


def scale_ctl_rootage_text(text: str, scale: Decimal, direction: str) -> str:
    if scale == 1:
        return text
    return "".join(
        scale_rootage_line(line, scale, direction)
        for line in text.splitlines(keepends=True)
    )


TREE_ASSIGNMENT_RE = re.compile(r"^\s*(?:UTREE|TREE)\s+(?:\*\s+)?(?:'[^']*'|[^\s=]+)\s*=\s*", re.IGNORECASE)


def uncomment_newick(text: str, comment_spans: list[tuple[int, int]] | None = None) -> str:
    """Strip Newick annotations while preserving quoted tip/internal labels."""
    result = []
    depth = 0
    comment_start = 0
    quoted = False
    index = 0
    while index < len(text):
        char = text[index]
        if depth:
            if char == "[":
                depth += 1
            elif char == "]":
                depth -= 1
                if depth == 0 and comment_spans is not None:
                    comment_spans.append((comment_start, index + 1))
        elif char == "'":
            result.append(char)
            if quoted and index + 1 < len(text) and text[index + 1] == "'":
                result.append("'")
                index += 1
            else:
                quoted = not quoted
        elif char == "[" and not quoted:
            depth = 1
            comment_start = index
        else:
            result.append(char)
        index += 1
    if depth or quoted:
        raise ValueError("Unterminated Newick comment or quoted label")
    return "".join(result)


def figtree_newick(statement: str) -> str:
    newick = TREE_ASSIGNMENT_RE.sub("", statement).strip()
    # Native NEXUS trees may carry an [&R] prefix before the root.
    while newick.startswith("["):
        closing = newick.find("]")
        if closing < 0:
            raise ValueError("Unterminated Newick root annotation")
        newick = newick[closing + 1:].lstrip()
    return newick


def figtree_tree_kind(statement: str) -> str | None:
    from ete4 import Tree
    from ete4.parser.newick import NewickError

    try:
        newick = figtree_newick(statement)
        comment_spans = []
        text = uncomment_newick(newick, comment_spans)
        tree = Tree(text, parser=1)
        leaves = list(tree.leaves())
        names = [leaf.name for leaf in leaves]
        if len(names) < 2 or any(not name for name in names) or len(set(names)) != len(names):
            return None
        distances = [node.dist for node in tree.traverse() if node is not tree]
        if all(distance is None for distance in distances):
            return "topology"
        if any(distance is None or not math.isfinite(distance) or distance < 0 for distance in distances):
            return None
        if tree.dist is not None and (not math.isfinite(tree.dist) or tree.dist < 0):
            return None
        for start, end in comment_spans:
            annotation = newick[start:end]
            if not HPD_COMMENT_RE.fullmatch(annotation):
                continue
            intervals = list(BRACE_RE.finditer(annotation))
            if not intervals:
                return None
            for interval in intervals:
                values = [parse_decimal(value) for value in interval.group("body").split(",")]
                if len(values) != 2 or any(value is None or not value.is_finite() or value < 0 for value in values):
                    return None
                if values[0] > values[1]:
                    return None
        return "dated"
    except (ValueError, NewickError):
        return None


def figtree_statements(text: str) -> tuple[str, list[str]]:
    """Read validated tree statements from a public block or native NEXUS file."""
    header = "Species tree for FigTree"
    in_block = False
    nexus = bool(re.search(r"^\s*#NEXUS\s*$", text, re.IGNORECASE | re.MULTILINE))
    if nexus and re.search(r"^\s*TRANSLATE\b", text, re.IGNORECASE | re.MULTILINE):
        return header, []
    statements = []
    kinds = []
    pending = []
    for line in text.splitlines():
        if "Species tree for FigTree" in line:
            if pending:
                return header, []
            header = line
            in_block = True
            continue
        if nexus and re.match(r"^\s*BEGIN\s+TREES\s*;", line, re.IGNORECASE):
            in_block = True
            continue
        if not in_block:
            continue
        if not pending:
            if not (line.lstrip().startswith("(") or re.match(r"^\s*(?:UTREE|TREE)\b", line, re.IGNORECASE)):
                continue
        pending.append(line)
        statement = "\n".join(pending)
        try:
            complete = uncomment_newick(statement).rstrip().endswith(";")
        except ValueError:
            complete = False
        if not complete:
            continue
        kind = figtree_tree_kind(statement)
        if kind is None:
            return header, []
        statements.append(statement)
        kinds.append(kind)
        pending = []
        if len(statements) == 3:
            break
    if pending or "dated" not in kinds:
        return header, []
    return header, statements


def extract_figtree_text(text: str, scale: Decimal, direction: str) -> str:
    header, statements = figtree_statements(text)
    if not statements:
        return ""
    return "\n".join([header] + [
        scale_newick_time_values(line, scale, direction) if figtree_tree_kind(line) == "dated" else line
        for line in statements
    ]) + "\n"


def has_figtree_tree(text: str) -> bool:
    return bool(figtree_statements(text)[1])


def cmd_conversion_inputs(args: argparse.Namespace) -> int:
    _, statements = figtree_statements(read_text(args.infile))
    if not statements:
        print(f"No valid dated FigTree tree was found in {args.infile}", file=sys.stderr)
        return 1
    newick = figtree_newick([line for line in statements if figtree_tree_kind(line) == "dated"][-1])
    comment_spans = []
    uncomment_newick(newick, comment_spans)
    without_ci = newick
    for start, end in reversed(comment_spans):
        if HPD_COMMENT_RE.fullmatch(newick[start:end]):
            without_ci = without_ci[:start] + without_ci[end:]
    without_ci = re.sub(r":\s+", ":", without_ci)
    write_text(args.outdir / "mcmctree_95CI.nwk", newick + "\n")
    write_text(args.outdir / "mcmctree_no95CI.nwk", without_ci + "\n")
    return 0


def parse_scale(value: str) -> Decimal:
    scale = parse_decimal(value)
    if scale is None or scale <= 0:
        raise argparse.ArgumentTypeError("scale must be a positive number")
    return scale


def cmd_factor(args: argparse.Namespace) -> int:
    target_max = parse_scale(args.target_max)
    factor = choose_scale_factor(read_text(args.infile), target_max)
    print(format_decimal(factor))
    return 0


def cmd_scale_calibrations(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = scale_calibration_labels(read_text(args.infile), scale, args.direction)
    write_text(args.outfile, text)
    return 0


def cmd_scale_figtree(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = scale_figtree_text(read_text(args.infile), scale, args.direction)
    write_text(args.outfile, text)
    return 0


def cmd_scale_ctl_rootage(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = scale_ctl_rootage_text(read_text(args.infile), scale, args.direction)
    write_text(args.outfile, text)
    return 0


def cmd_extract_figtree(args: argparse.Namespace) -> int:
    scale = parse_scale(args.scale)
    text = extract_figtree_text(read_text(args.infile), scale, args.direction)
    if not text:
        print(
            f"No FigTree tree block was found in {args.infile}",
            file=sys.stderr,
        )
        return 1
    write_text(args.outfile, text)
    return 0


def cmd_validate_figtree(args: argparse.Namespace) -> int:
    if has_figtree_tree(read_text(args.infile)):
        return 0
    print(
        f"No FigTree tree block was found in {args.infile}",
        file=sys.stderr,
    )
    return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Scale MCMCTree time values between public and internal units."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    factor = subparsers.add_parser(
        "factor",
        help="Choose a power-of-ten divisor that keeps calibration ages <= target.",
    )
    factor.add_argument("--infile", type=Path, required=True)
    factor.add_argument("--target-max", default="10")
    factor.set_defaults(func=cmd_factor)

    scale_cal = subparsers.add_parser(
        "scale-calibrations",
        help="Scale B/L/U calibration labels in a Newick file.",
    )
    scale_cal.add_argument("--infile", type=Path, required=True)
    scale_cal.add_argument("--outfile", type=Path, required=True)
    scale_cal.add_argument("--scale", required=True)
    scale_cal.add_argument("--direction", choices=["down", "up"], required=True)
    scale_cal.set_defaults(func=cmd_scale_calibrations)

    scale_fig = subparsers.add_parser(
        "scale-figtree",
        help="Scale MCMCTree FigTree/Newick time values.",
    )
    scale_fig.add_argument("--infile", type=Path, required=True)
    scale_fig.add_argument("--outfile", type=Path, required=True)
    scale_fig.add_argument("--scale", required=True)
    scale_fig.add_argument("--direction", choices=["down", "up"], required=True)
    scale_fig.set_defaults(func=cmd_scale_figtree)

    scale_ctl = subparsers.add_parser(
        "scale-ctl-rootage",
        help="Scale RootAge values in a MCMCTree control file.",
    )
    scale_ctl.add_argument("--infile", type=Path, required=True)
    scale_ctl.add_argument("--outfile", type=Path, required=True)
    scale_ctl.add_argument("--scale", required=True)
    scale_ctl.add_argument("--direction", choices=["down", "up"], required=True)
    scale_ctl.set_defaults(func=cmd_scale_ctl_rootage)

    extract = subparsers.add_parser(
        "extract-figtree",
        help="Extract and scale the FigTree block from an MCMCTree output file.",
    )
    extract.add_argument("--infile", type=Path, required=True)
    extract.add_argument("--outfile", type=Path, required=True)
    extract.add_argument("--scale", required=True)
    extract.add_argument("--direction", choices=["down", "up"], required=True)
    extract.set_defaults(func=cmd_extract_figtree)

    validate = subparsers.add_parser(
        "validate-figtree",
        help="Require at least one Newick tree in a FigTree section.",
    )
    validate.add_argument("--infile", type=Path, required=True)
    validate.set_defaults(func=cmd_validate_figtree)

    conversion = subparsers.add_parser("conversion-inputs", help="Extract validated CI sidecars from a public FigTree artifact.")
    conversion.add_argument("--infile", type=Path, required=True)
    conversion.add_argument("--outdir", type=Path, required=True)
    conversion.set_defaults(func=cmd_conversion_inputs)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
