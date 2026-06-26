import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RUNTIME_SOURCE_ROOTS = (
    REPO_ROOT / "workflow" / "core",
    REPO_ROOT / "workflow" / "support",
    REPO_ROOT / "container",
    REPO_ROOT / ".github" / "workflows",
)
RUNTIME_SOURCE_SUFFIXES = {".py", ".sh", ".R", ".r", ".yml", ".yaml", ".md"}
BANNED_RUNTIME_DEBT_PATTERN = re.compile(r"\b(TODO|FIXME|HACK|workaround)\b", re.IGNORECASE)


def runtime_source_files():
    for root in RUNTIME_SOURCE_ROOTS:
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.suffix in RUNTIME_SOURCE_SUFFIXES:
                yield path


def test_runtime_sources_do_not_carry_dependency_workaround_debt_markers():
    offenders = []
    for path in runtime_source_files():
        text = path.read_text(encoding="utf-8", errors="replace")
        for lineno, line in enumerate(text.splitlines(), start=1):
            if BANNED_RUNTIME_DEBT_PATTERN.search(line):
                offenders.append("{}:{}:{}".format(path.relative_to(REPO_ROOT), lineno, line.strip()))

    assert offenders == []
