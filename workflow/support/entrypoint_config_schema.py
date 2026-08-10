#!/usr/bin/env python3
"""Generate and validate configuration metadata from editable entrypoint blocks."""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
REGISTRY_PATH = WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh"
ENTRYPOINT_RE = re.compile(r"^\s{4}(gg_[a-z_]+_entrypoint\.sh)\)$")
ASSIGNMENT_RE = re.compile(r"^([a-z][a-z0-9_]*)=(.*?)(?:\s+#\s*(.*))?$")
COMMON_RE = re.compile(r'^: "\$\{(GG_COMMON_[A-Z0-9_]+):=(.*)\}"\s+#\s*(.*)$')
CHOICES_RE = re.compile(r"\{([^{}]+)\}")
PREFIXES = {
    "gg_input_generation_entrypoint.sh": "GG_INPUT_",
    "gg_transcriptome_generation_entrypoint.sh": "GG_TRANSCRIPTOME_",
    "gg_genome_annotation_entrypoint.sh": "GG_GENOME_ANNOTATION_",
    "gg_genome_evolution_entrypoint.sh": "GG_GENOME_EVOLUTION_",
    "gg_gene_evolution_entrypoint.sh": "GG_GENE_EVOLUTION_",
    "gg_gene_summary_entrypoint.sh": "GG_GENE_SUMMARY_",
    "gg_progress_summary_entrypoint.sh": "GG_PROGRESS_SUMMARY_",
}


@dataclass(frozen=True)
class Parameter:
    name: str
    environment: str
    default_expression: str
    default: str
    type: str
    choices: list[str]
    description: str
    line: int


def strip_shell_value(expression: str) -> str:
    value = expression.strip()
    expansion = re.fullmatch(r'"?\$\{[A-Za-z_][A-Za-z0-9_]*:-(.*)\}"?', value)
    if expansion:
        value = expansion.group(1)
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        value = value[1:-1]
    return value


def infer_type(name: str, default: str, choices: list[str]) -> str:
    values = choices or [default]
    if values and all(value in {"0", "1"} for value in values) and name.startswith(("run_", "delete_")):
        return "boolean"
    if default and re.fullmatch(r"-?\d+", default):
        return "integer"
    if default and re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)", default):
        return "number"
    return "string"


def parse_choices(description: str) -> list[str]:
    match = CHOICES_RE.search(description)
    if match is None:
        return []
    return [item.strip().strip('"\'') for item in match.group(1).split(",")]


def parse_registry() -> dict[str, list[str]]:
    registry: dict[str, list[str]] = {}
    current = ""
    collecting = False
    for line in REGISTRY_PATH.read_text(encoding="utf-8").splitlines():
        match = ENTRYPOINT_RE.match(line)
        if match:
            current = match.group(1)
            registry[current] = []
            collecting = False
            continue
        if current and line.strip() == "cat <<'EOF'":
            collecting = True
            continue
        if collecting and line == "EOF":
            collecting = False
            current = ""
            continue
        if collecting:
            registry[current].append(line.strip())
    return registry


def editable_lines(path: Path) -> list[tuple[int, str]]:
    active = False
    result: list[tuple[int, str]] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.startswith("### Start: Modify this block"):
            active = True
            continue
        if line.startswith("### End: Modify this block"):
            break
        if active:
            result.append((number, line))
    return result


def parse_entrypoint(name: str) -> list[Parameter]:
    prefix = PREFIXES[name]
    parameters: list[Parameter] = []
    for line_number, line in editable_lines(WORKFLOW_DIR / name):
        match = ASSIGNMENT_RE.match(line)
        if match is None:
            continue
        variable, expression, description = match.groups()
        description = (description or "").strip()
        default = strip_shell_value(expression)
        choices = parse_choices(description)
        parameters.append(
            Parameter(
                name=variable,
                environment=prefix + variable.upper(),
                default_expression=expression.strip(),
                default=default,
                type=infer_type(variable, default, choices),
                choices=choices,
                description=description,
                line=line_number,
            )
        )
    return parameters


def parse_all_assignments(name: str) -> dict[str, Parameter]:
    prefix = PREFIXES[name]
    result: dict[str, Parameter] = {}
    path = WORKFLOW_DIR / name
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = ASSIGNMENT_RE.match(line)
        if match is None:
            continue
        variable, expression, description = match.groups()
        description = (description or "").strip()
        default = strip_shell_value(expression)
        choices = parse_choices(description)
        result[variable] = Parameter(
            name=variable,
            environment=prefix + variable.upper(),
            default_expression=expression.strip(),
            default=default,
            type=infer_type(variable, default, choices),
            choices=choices,
            description=description,
            line=line_number,
        )
    return result


def parse_common() -> list[Parameter]:
    result: list[Parameter] = []
    for number, line in enumerate((WORKFLOW_DIR / "gg_common_params.sh").read_text(encoding="utf-8").splitlines(), 1):
        match = COMMON_RE.match(line)
        if match is None:
            continue
        name, default, description = match.groups()
        choices = parse_choices(description)
        result.append(
            Parameter(
                name=name,
                environment=name,
                default_expression=f"${{{name}:={default}}}",
                default=strip_shell_value(default),
                type=infer_type(name.lower(), strip_shell_value(default), choices),
                choices=choices,
                description=description.strip(),
                line=number,
            )
        )
    return result


def build_schema() -> dict[str, object]:
    registry = parse_registry()
    common = parse_common()
    common_by_local_name = {
        item.name.removeprefix("GG_COMMON_").lower(): item for item in common
    }
    entrypoints: dict[str, object] = {}
    for name in sorted(PREFIXES):
        editable = {item.name: item for item in parse_entrypoint(name)}
        assigned = parse_all_assignments(name)
        effective: list[Parameter] = []
        for parameter_name in registry.get(name, []):
            parameter = editable.get(parameter_name) or assigned.get(parameter_name)
            if parameter is None and parameter_name in common_by_local_name:
                inherited = common_by_local_name[parameter_name]
                parameter = Parameter(
                    name=parameter_name,
                    environment=PREFIXES[name] + parameter_name.upper(),
                    default_expression=inherited.default_expression,
                    default=inherited.default,
                    type=inherited.type,
                    choices=inherited.choices,
                    description=f"Inherited common setting. {inherited.description}",
                    line=inherited.line,
                )
            if parameter is not None:
                effective.append(parameter)
        entrypoints[name] = {
            "environment_prefix": PREFIXES[name],
            "parameters": [asdict(item) for item in effective],
            "editable_parameters": sorted(editable),
            "registered_parameters": registry.get(name, []),
        }
    return {
        "schema_version": 1,
        "source": "workflow entrypoint editable blocks",
        "common_parameters": [asdict(item) for item in common],
        "entrypoints": entrypoints,
    }


def validate(schema: dict[str, object]) -> list[str]:
    errors: list[str] = []
    registry = parse_registry()
    if set(registry) != set(PREFIXES):
        errors.append("config registry entrypoints do not match the supported entrypoint set")
    entrypoints = schema["entrypoints"]
    assert isinstance(entrypoints, dict)
    for name, entrypoint in entrypoints.items():
        assert isinstance(entrypoint, dict)
        parameters = entrypoint["parameters"]
        assert isinstance(parameters, list)
        parsed_names = [item["name"] for item in parameters]
        registered = registry.get(name, [])
        editable = entrypoint["editable_parameters"]
        missing = sorted(set(editable) - set(registered))
        unresolved = sorted(set(registered) - set(parsed_names))
        if missing:
            errors.append(f"{name}: editable parameters missing from registry: {', '.join(missing)}")
        if unresolved:
            errors.append(f"{name}: registered parameters have no editable, entrypoint, or common default: {', '.join(unresolved)}")
        if len(parsed_names) != len(set(parsed_names)):
            errors.append(f"{name}: duplicate editable parameter assignments")
        for parameter in parameters:
            if not parameter["description"]:
                errors.append(f"{name}:{parameter['line']}: {parameter['name']} has no description")
    return errors


def render_markdown(schema: dict[str, object]) -> str:
    lines = ["# Generated entrypoint configuration reference", "", "Generated from the editable blocks in `workflow/gg_*_entrypoint.sh`.", ""]
    entrypoints = schema["entrypoints"]
    assert isinstance(entrypoints, dict)
    for name, entrypoint in entrypoints.items():
        lines.extend([f"## `{name}`", "", "| Parameter | Environment override | Default | Description |", "|---|---|---|---|"])
        for item in entrypoint["parameters"]:
            default = str(item["default"]).replace("|", "\\|") or "_(empty)_"
            description = str(item["description"]).replace("|", "\\|")
            lines.append(f"| `{item['name']}` | `{item['environment']}` | `{default}` | {description} |")
        lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail if editable blocks and the forwarding registry differ.")
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    schema = build_schema()
    errors = validate(schema)
    if args.check and errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    if args.check:
        print(f"Validated {len(schema['entrypoints'])} entrypoints and {len(schema['common_parameters'])} common parameters.")
        return 0
    if args.format == "markdown":
        print(render_markdown(schema))
    else:
        print(json.dumps(schema, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
