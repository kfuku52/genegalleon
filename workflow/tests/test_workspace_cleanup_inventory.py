import importlib.util
import os
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "support/workspace_cleanup_inventory.py"
spec = importlib.util.spec_from_file_location("cleanup_inventory", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_inventory_distinguishes_scratch_from_durable_state_and_keeps_every_byte(tmp_path):
    files = {
        "output/genome_evolution/omark/Species_a/Species_a.query.fa": b">a\nMP\n",
        "output/genome_evolution/omark/Species_a/Species_a.sum": b"complete\n",
        "downloads/tmp/species_genetic_code.resolved.tsv": b"code\n",
        "output/species_tree/tmp/failed/input.fa": b"scratch\n",
        "output/query2family/tmp/1_AHA/failed.txt": b"failed\n",
        "output/input_generation/tmp/task_plan.json": b"immutable plan\n",
        "output/input_generation/tmp/task_plan.json.completed/1.json": b"receipt\n",
        "output/species_tree/mcmctree_main/runs/chain-01/mcmc.txt": b"draws\n",
        "downloads/example/x.gz.corrupt.20261001123456.42": b"bad\n",
        "downloads/example/x.gz.part": b"resumable\n",
        "downloads/example/.archive_cache/a.zip": b"cached\n",
        "output/orthofinder/core/Orthogroups/Orthogroups.tsv": b"scientific result\n",
    }
    for relative, data in files.items():
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.rglob("*") if p.is_file()}
    report = module.inventory(tmp_path)
    records = {row["path"]: row for row in report["records"]}
    assert report["read_only"]
    assert len(records) == 7
    assert records[str(tmp_path / "output/species_tree/tmp")]["files"] == 1
    assert records[str(tmp_path / "output/orthofinder/core")]["policy"] == "retain_unless_explicitly_retired"
    assert all("task_plan" not in path and "mcmctree_main" not in path and not path.endswith(".part")
               for path in records)
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == value for p, value in before.items())


def test_inventory_never_follows_symlinked_scratch(tmp_path):
    workspace, foreign = tmp_path / "workspace", tmp_path / "foreign"
    foreign.mkdir()
    (foreign / "large.fa").write_bytes(b"external data\n")
    root = workspace / "output/species_tree"
    root.mkdir(parents=True)
    (root / "tmp").symlink_to(foreign, target_is_directory=True)
    record = module.inventory(workspace)["records"][0]
    assert record["files"] == 0
    assert record["symlinks"] == 1
    assert record["identity"]["uid"] == os.getuid()
    assert (foreign / "large.fa").read_bytes() == b"external data\n"


def test_legacy_output_layout(tmp_path):
    scratch = tmp_path / "gfe_data/species_tree/tmp/job"
    scratch.mkdir(parents=True)
    (scratch / "keep").write_bytes(b"failed\n")
    report = module.inventory(tmp_path / "gfe_data", legacy=True)
    assert report["records"][0]["files"] == 1


def test_inventory_does_not_list_scratch_through_symlinked_parent(tmp_path, monkeypatch):
    workspace, foreign = tmp_path / "workspace", tmp_path / "foreign"
    (foreign / "tmp/1_task").mkdir(parents=True)
    (foreign / "tmp/1_task/input.fa").write_bytes(b"external data\n")
    (workspace / "output").mkdir(parents=True)
    alias = workspace / "output/orthogroup"
    alias.symlink_to(foreign, target_is_directory=True)
    original = Path.iterdir

    def guarded_iterdir(path):
        assert not path.is_relative_to(alias), "traversed a symlinked parent"
        return original(path)

    monkeypatch.setattr(Path, "iterdir", guarded_iterdir)
    report = module.inventory(workspace)
    assert not report["records"]
    assert any(row["error"] == "symlinked parent: not traversed" for row in report["errors"])
    assert (foreign / "tmp/1_task/input.fa").read_bytes() == b"external data\n"
