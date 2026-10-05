import csv
import gzip
import hashlib
import json
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/plan_input_cohort_export.py"


def fixture(tmp_path):
    rows = []
    for directory in ("cds", "gff", "genome", "short", "full", "provenance"):
        (tmp_path / directory).mkdir()
    for species, complete in (("Plant_a", 64), ("Plant_b", 65), ("Plant_c", 129)):
        files = {}
        for role in ("cds", "gff", "genome"):
            path = tmp_path / role / f"{species}_{role}.gz"
            path.write_bytes(gzip.compress(b"synthetic fixture\n", mtime=0))
            files[role] = path
        full = tmp_path / "full" / f"{species}.busco.full.tsv"
        # Duplicated BUSCO hits count one BUSCO group, not two genes.
        full.write_text("\n".join(f"id{i}\t{'Duplicated' if i < complete else 'Missing'}\tgene{i}"
                                  for i in range(129)) + "\nid0\tDuplicated\tcopy2\n")
        short = tmp_path / "short" / f"{species}.busco.short.txt"
        short.write_text("# BUSCO version is: 6.1.0\n# The lineage dataset is: eukaryota_odb12\n"
                         "# BUSCO was run in mode: transcriptome\n"
                         f"C:{100*complete/129:.1f}%[S:0.0%,D:{100*complete/129:.1f}%],F:0.0%,M:0.0%,n:129\n")
        def entry(label, path):
            return {"label": label, "scope": "workspace", "path": str(path.relative_to(tmp_path)),
                    "size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        provenance = {"schema_version": 1, "step": "input_generation_species_busco", "family_id": species,
                      "parameters": {"busco_mode": "transcriptome", "busco_lineage_resolved": "eukaryota_odb12"},
                      "inputs": [entry("species_cds", files["cds"])],
                      "outputs": [entry("busco_full", full), entry("busco_short", short)]}
        (tmp_path / "provenance" / f"busco.{species}.json").write_text(json.dumps(provenance))
        rows.append({"species_prefix": species, "taxid": "123", **{k + "_output_path": str(v) for k, v in files.items()}})
    with (tmp_path / "summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (tmp_path / "traits.tsv").write_text("species\tcarnivory\nPlant_a\t0\nPlant_b\t1\nPlant_c\t0\n")
    return ["--workspace-root", str(tmp_path), "--species-summary", str(tmp_path / "summary.tsv"),
            "--busco-short-dir", str(tmp_path / "short"), "--busco-full-dir", str(tmp_path / "full"),
            "--busco-provenance-dir", str(tmp_path / "provenance"), "--minimum-complete-busco", "50",
            "--target-part-gb", "0.00000015", "--metadata-table", f"traits={tmp_path / 'traits.tsv'}",
            "--output", str(tmp_path / "export")]


def test_export_uses_exact_busco_groups_and_keeps_triplets(tmp_path):
    args = fixture(tmp_path)
    before = {path: path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}
    result = subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / "export/export_plan.json").read_text())
    assert report["selected_species"] == 2 and report["triplet_file_count"] == 6
    assert report["excluded"] == [{"species": "Plant_a", "complete_percent": 100 * 64 / 129}]
    assert sorted(x for part in report["parts"] for x in part["species"]) == ["Plant_b", "Plant_c"]
    assert len({file["export_path"] for file in report["source_files"]}) == 6
    assert report["gzip_bytes"] == sum(file["size_bytes"] for file in report["source_files"])
    for species in ("Plant_b", "Plant_c"):
        assert len({file["part"] for file in report["source_files"] if file["species"] == species}) == 1
    assert (tmp_path / "export/traits.tsv").read_text() == "species\tcarnivory\nPlant_b\t1\nPlant_c\t0\n"
    assert all(path.read_bytes() == data for path, data in before.items())
    for line in (tmp_path / "export/SHA256SUMS").read_text().splitlines():
        sha, name = line.split("  ")
        assert hashlib.sha256((tmp_path / "export" / name).read_bytes()).hexdigest() == sha
    again = subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)
    assert again.returncode != 0 and "new directory" in again.stderr


@pytest.mark.parametrize("invalid", ["changed_cds", "changed_busco", "short_full_disagree", "mixed_lineage", "missing_triplet", "duplicate_species", "missing_trait", "contradictory_groups"])
def test_export_rejects_inconsistent_inputs_before_publishing(tmp_path, invalid):
    args = fixture(tmp_path)
    if invalid == "changed_cds":
        (tmp_path / "cds/Plant_b_cds.gz").write_bytes(b"changed")
    elif invalid == "changed_busco":
        (tmp_path / "full/Plant_b.busco.full.tsv").write_text("id\tMissing\n")
    elif invalid == "short_full_disagree":
        path = tmp_path / "short/Plant_b.busco.short.txt"
        path.write_text(path.read_text().replace("50.4", "50.0"))
    elif invalid == "mixed_lineage":
        path = tmp_path / "short/Plant_b.busco.short.txt"
        path.write_text(path.read_text().replace("eukaryota_odb12", "other_odb12"))
    elif invalid == "missing_triplet":
        (tmp_path / "gff/Plant_b_gff.gz").unlink()
    elif invalid == "duplicate_species":
        path = tmp_path / "summary.tsv"
        path.write_text(path.read_text() + path.read_text().splitlines()[1] + "\n")
    elif invalid == "missing_trait":
        (tmp_path / "traits.tsv").write_text("species\tcarnivory\nPlant_c\t0\n")
    elif invalid == "contradictory_groups":
        path = tmp_path / "full/Plant_b.busco.full.tsv"
        path.write_text(path.read_text() + "id0\tComplete\tgene\n")
    result = subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)
    assert result.returncode != 0
    assert not (tmp_path / "export").exists()


def test_export_fences_sources_from_size_collection_through_hashing(tmp_path, monkeypatch):
    import argparse

    fixture(tmp_path)
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    spec = spec_from_file_location("export_planner", SCRIPT)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    real_digest = module.digest_paths

    def modify_then_hash(paths):
        path = tmp_path / "genome/Plant_b_genome.gz"
        path.write_bytes(gzip.compress(b"different size and content\n", mtime=0))
        return real_digest(paths)

    monkeypatch.setattr(module, "digest_paths", modify_then_hash)
    args = argparse.Namespace(workspace_root=tmp_path, species_summary=tmp_path / "summary.tsv",
                              busco_short_dir=tmp_path / "short", busco_full_dir=tmp_path / "full",
                              busco_provenance_dir=tmp_path / "provenance", minimum_complete_busco=50,
                              target_part_gb=20, metadata_table=[])
    with pytest.raises(ValueError, match="Input changed while planning"):
        module.build(args)
