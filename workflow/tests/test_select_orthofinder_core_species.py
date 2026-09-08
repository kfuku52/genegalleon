import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "select_orthofinder_core_species.py"


def _write_protein(path: Path, count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    records = [f">{path.stem}_gene{i}\nMPEP\n" for i in range(1, count + 1)]
    path.write_text("".join(records), encoding="utf-8")


def _write_busco(path: Path, complete_pct: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"C:{complete_pct}%[S:{complete_pct}%,D:0.0%],F:0.0%,M:0.0%,n:100\n",
        encoding="utf-8",
    )


def _run_selector(tmp_path: Path, extra_args: list[str], env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(SCRIPT),
        "--protein-dir",
        str(tmp_path / "protein"),
        "--busco-short-dir",
        str(tmp_path / "busco"),
        "--max-core-species",
        "1",
        "--filters",
        "busco_complete_pct:ge:80,num_seq:le:100000",
        "--rank",
        "num_seq:asc,busco_complete_pct:desc",
        "--method",
        "max-pd",
        "--candidates-table",
        str(out_dir / "candidates.tsv"),
        "--selected-table",
        str(out_dir / "selected.tsv"),
        "--selected-list",
        str(out_dir / "selected.txt"),
        "--core-tree",
        str(out_dir / "core.nwk"),
        *extra_args,
    ]
    return subprocess.run(command, capture_output=True, text=True, env=env, check=False)


def test_select_orthofinder_core_species_falls_back_when_busco_values_are_missing(tmp_path: Path):
    _write_protein(tmp_path / "protein" / "Large_species.fa", 3)
    _write_protein(tmp_path / "protein" / "Small_species.fa", 1)

    completed = _run_selector(tmp_path, [])

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "retrying without BUSCO completeness filters" in completed.stderr
    assert (tmp_path / "out" / "selected.txt").read_text(encoding="utf-8") == "Small_species.fa\n"
    selected = (tmp_path / "out" / "selected.tsv").read_text(encoding="utf-8")
    assert "Small_species" in selected
    assert "Large_species" not in selected


def test_find_busco_short_file_does_not_match_longer_species_label(tmp_path: Path):
    import importlib.util

    spec = importlib.util.spec_from_file_location("select_orthofinder_core_species", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    _write_busco(tmp_path / "busco" / "Species_a_subsp_x.busco.short.txt", 95.0)

    assert mod.find_busco_short_file(tmp_path / "busco", "Species_a") is None


def test_find_busco_short_file_rejects_ambiguous_species_matches(tmp_path: Path):
    import importlib.util

    spec = importlib.util.spec_from_file_location("select_orthofinder_core_species", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _write_busco(tmp_path / "busco" / "Species_a.first.busco.short.txt", 95.0)
    _write_busco(tmp_path / "busco" / "Species_a.second.busco.short.txt", 96.0)

    with pytest.raises(SystemExit, match="Refusing to treat this ambiguity"):
        mod.find_busco_short_file(tmp_path / "busco", "Species_a")


def test_select_orthofinder_core_species_calls_nwkit_sample_with_default_filters(tmp_path: Path):
    _write_protein(tmp_path / "protein" / "Good_species.fa", 1)
    _write_protein(tmp_path / "protein" / "Lowbusco_species.fa", 1)
    _write_busco(tmp_path / "busco" / "Good_species.busco.short.txt", 90.0)
    _write_busco(tmp_path / "busco" / "Lowbusco_species.busco.short.txt", 70.0)
    tree = tmp_path / "tree.nwk"
    tree.write_text("(Good_species:0.1,Lowbusco_species:0.1);\n", encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    nwkit = bin_dir / "nwkit"
    nwkit.write_text(
        """#!/usr/bin/env python3
import csv
from pathlib import Path
import sys

args = sys.argv[1:]
Path(sys.argv[0]).with_suffix(".args").write_text(" ".join(args) + "\\n", encoding="utf-8")
trait = args[args.index("--trait") + 1]
report = args[args.index("--report") + 1]
outfile = args[args.index("--outfile") + 1]
filters = [args[idx + 1] for idx, arg in enumerate(args) if arg == "--filter"]

def passes(row, spec):
    column, op, value = spec.split(":", 2)
    left = float(row[column])
    right = float(value)
    return (op == "ge" and left >= right) or (op == "le" and left <= right)

with open(trait, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\\t"))
selected = [row for row in rows if all(passes(row, spec) for spec in filters)][:1]
selected = [{"sample_order": index, **row} for index, row in enumerate(selected, start=1)]
with open(report, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sample_order", *rows[0].keys()], delimiter="\\t", lineterminator="\\n")
    writer.writeheader()
    writer.writerows(selected)
Path(outfile).write_text("(" + selected[0]["leaf_name"] + ":0.1);\\n", encoding="utf-8")
""",
        encoding="utf-8",
    )
    nwkit.chmod(nwkit.stat().st_mode | stat.S_IXUSR)
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"

    completed = _run_selector(tmp_path, ["--species-tree", str(tree)], env=env)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    nwkit_args = nwkit.with_suffix(".args").read_text(encoding="utf-8")
    assert "sample" in nwkit_args
    assert "--filter busco_complete_pct:ge:80" in nwkit_args
    assert "--filter num_seq:le:100000" in nwkit_args
    assert "--method max-pd" in nwkit_args
    assert "--report" in nwkit_args
    assert "--output-table" not in nwkit_args
    selected_header = (tmp_path / "out" / "selected.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert selected_header.startswith("sample_order\tleaf_name\t")
    assert (tmp_path / "out" / "selected.txt").read_text(encoding="utf-8") == "Good_species.fa\n"
