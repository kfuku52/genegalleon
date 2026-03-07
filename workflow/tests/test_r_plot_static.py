from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
NUMERIC_LEGEND_POSITION_RE = re.compile(r"legend\.position\s*=\s*c\(")


def test_r_scripts_do_not_use_deprecated_numeric_legend_position():
    offending = []
    for path in SUPPORT_DIR.rglob("*"):
        if path.suffix not in {".r", ".R"}:
            continue
        text = path.read_text(encoding="utf-8")
        if NUMERIC_LEGEND_POSITION_RE.search(text):
            offending.append(path.relative_to(REPO_ROOT).as_posix())

    assert offending == [], (
        "Found deprecated numeric legend.position usage in: "
        + ", ".join(offending)
    )
