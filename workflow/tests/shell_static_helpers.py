import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
CORE_DIR = WORKFLOW_DIR / "core"
CONTAINER_SCRIPTS_DIR = REPO_ROOT / "container" / "scripts"
GITHUB_WORKFLOWS_DIR = REPO_ROOT / ".github" / "workflows"


def read_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    if path == WORKFLOW_DIR / "support" / "gg_util.sh":
        module_dir = path.parent / "gg_util"
        module_text = "\n".join(module.read_text(encoding="utf-8") for module in sorted(module_dir.glob("*.sh")))
        return f"{text}\n{module_text}"
    if path.parent == CORE_DIR and path.name.endswith("_core.sh"):
        stage_library = CORE_DIR / "stages" / f"{path.stem}_functions.sh"
        execution_stage_dir = CORE_DIR / "stages" / path.stem.removesuffix("_core")
        parts = []
        if stage_library.is_file():
            parts.append(stage_library.read_text(encoding="utf-8"))
        parts.append(text)
        if execution_stage_dir.is_dir():
            parts.extend(
                stage.read_text(encoding="utf-8") for stage in sorted(execution_stage_dir.glob("*.sh"))
            )
        return "\n".join(parts)
    return text


def function_body(text: str, function_name: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(function_name)}\(\)\s*\{{", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"Function not found: {function_name}")
    start = match.start()
    next_match = re.search(r"^\s*[A-Za-z_][A-Za-z0-9_]*\(\)\s*\{", text[match.end() :], re.MULTILINE)
    if next_match is None:
        return text[start:]
    return text[start : match.end() + next_match.start()]


def workflow_shell_scripts():
    return sorted(WORKFLOW_DIR.rglob("*.sh"))


def container_shell_scripts():
    return sorted(CONTAINER_SCRIPTS_DIR.rglob("*.sh"))


def core_and_entrypoint_scripts():
    core_scripts = sorted(CORE_DIR.glob("*.sh"))
    entrypoint_scripts = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    return core_scripts + entrypoint_scripts


def strict_mode_header(script: Path) -> str:
    max_lines = 60 if script.name.endswith("_entrypoint.sh") else 30
    return "\n".join(script.read_text(encoding="utf-8").splitlines()[:max_lines])


def entrypoint_scheduler_header(script: Path) -> str:
    text = read_text(script)
    marker = "set -euo pipefail"
    assert marker in text, f"Missing strict mode marker in {script}"
    return text.split(marker, 1)[0]
