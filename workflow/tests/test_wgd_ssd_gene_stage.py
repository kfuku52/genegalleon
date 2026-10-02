"""Execute the gene-evolution origin stage with real provenance and NWKIT."""

import csv
import hashlib
import json
import os
import shlex
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def stage(core, title):
    start = core.index(f'task="{title}"')
    end = core.index('\ntask=', start + 1)
    return core[start:end]


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(results, manifest):
    paths = [results, manifest]
    return {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


def test_gene_origin_stage_forwards_parameters_uses_full_tree_and_caches(tmp_path):
    import nwkit

    from workflow.support.gene_family_output_store import (
        GeneFamilyOutputStore,
        archive_queue_status,
        drain_archive_queue,
        enqueue_family_archive,
        family_inventory_path,
        query_id_extractor,
    )

    workspace = tmp_path / "workspace"
    inputs = workspace / "input"
    inputs.mkdir(parents=True)
    full_tree, pruned_tree = inputs / "species.full.nwk", inputs / "species.pruned.nwk"
    full_tree.write_text("(A:1,(B:1,C:1):1);\n", encoding="utf-8")
    pruned_tree.write_text("(A:1,B:1);\n", encoding="utf-8")
    gene_tree = inputs / "OGstage.rooted.nwk"
    gene_tree.write_text("((sampleA_dup1:0.1,override_dup2:0.1):0.2,"
                         "(sampleA_far1:0.1,sampleA_far2:0.1):0.2);\n", encoding="utf-8")
    species_map = inputs / "species-map.tsv"
    write_tsv(species_map, ["leaf_name", "species_label"], [
        {"leaf_name": "override_dup2", "species_label": "A"}])
    regex = r"^sample([A-Z])_"
    identity = f"nwkit {nwkit.__version__} gene-stage-test"

    genes = {"sampleA_dup1": 1, "override_dup2": 2, "sampleA_far1": 4, "sampleA_far2": 7}
    fasta, gff = inputs / "A.protein.fa", inputs / "A.gff3"
    fasta.write_text("".join(f">{name}\nMAAAA\n" for name in genes), encoding="utf-8")
    loci = ["sampleA_dup1", "override_dup2", "intervening3", "sampleA_far1",
            "intervening5", "intervening6", "sampleA_far2"]
    gff.write_text("##gff-version 3\n" + "".join(
        f"chr1\ttest\tgene\t{rank * 100 + 1}\t{rank * 100 + 90}\t.\t+\t.\tID={name}\n"
        for rank, name in enumerate(loci, 1)), encoding="utf-8")
    evidence = workspace / "output/genome_evolution/wgd_ssd"
    evidence.mkdir(parents=True)
    write_tsv(evidence / "gene_positions.tsv",
              ["gene_id", "species", "locus_id", "seqid", "rank", "start", "end"], [
                  {"gene_id": name, "species": "A", "locus_id": name, "seqid": "chr1", "rank": rank,
                   "start": rank * 100 + 1, "end": rank * 100 + 90} for name, rank in genes.items()])
    write_tsv(evidence / "anchor_evidence.tsv", ["species", "block_id", "gene_a", "gene_b", "ks",
                                                 "ks_status", "species_event_id", "placement_status"], [])
    write_tsv(evidence / "wgd_events.tsv", ["species_event_id", "event_support"], [])
    evidence_outputs = [evidence / name for name in ("gene_positions.tsv", "anchor_evidence.tsv", "wgd_events.tsv")]
    if (evidence / "count_model.json").is_file():
        evidence_outputs.append(evidence / "count_model.json")
    (evidence / "summary.json").write_text(json.dumps({
        "schema_version": 1, "plan": {"species_tree": str(full_tree),
                                     "input_hashes": {str(path): digest(path) for path in (full_tree, fasta, gff)},
                                     "absent_inputs": []},
        "num_positioned_genes": len(genes), "num_raw_anchor_rows": 0,
        "species_synteny": {"A": {"num_genes": len(genes), "num_annotated_loci": len(loci)}},
        "output_hashes": {path.name: digest(path) for path in evidence_outputs},
    }), encoding="utf-8")

    temporary = workspace / "tmp"
    temporary.mkdir()
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text(encoding="utf-8")
    destination = next(line for line in core.splitlines() if line.startswith("file_og_wgd_ssd="))
    register_start = core.index("# Include declared destinations")
    register_end = core.index("# Define intermediate files", register_start)
    script = tmp_path / "gene-origin-stage.sh"
    script.write_text(f'''set -euo pipefail
source {shlex.quote(str(SUPPORT / "gg_util.sh"))}
gg_support_dir={shlex.quote(str(SUPPORT))}
gg_workspace_dir={shlex.quote(str(workspace))}
gg_workspace_output_dir="${{gg_workspace_dir}}/output"
dir_output_active="${{gg_workspace_output_dir}}/query2family"
dir_tmp={shlex.quote(str(temporary))}
file_og_rooted_tree_analysis={shlex.quote(str(gene_tree))}
species_tree={shlex.quote(str(full_tree))}
species_tree_pruned={shlex.quote(str(pruned_tree))}
og_id=OGstage
GG_FAMILY_OUTPUT_INVENTORY=$(python "${{gg_support_dir}}/gene_family_output_store.py" inventory-path --root "${{dir_output_active}}" --family-id "${{og_id}}")
mkdir -p "${{GG_FAMILY_OUTPUT_INVENTORY}}"
export GG_FAMILY_OUTPUT_INVENTORY
export GG_FAMILY_OUTPUT_ROOT="${{dir_output_active}}"
run_wgd_ssd_classification=1
wgd_native_tree_likelihood=0
wgd_evidence_dir="${{WGD_TEST_EVIDENCE_DIR}}"
wgd_proximal_distance="${{WGD_TEST_PROXIMAL_DISTANCE}}"
species_label_parser=legacy
species_label_regex={shlex.quote(regex)}
species_label_map_tsv={shlex.quote(str(species_map))}
gene_nwkit_identity={shlex.quote(identity)}
artifact_stale_policy=rebuild
delete_tmp_dir=1
gene_family_output_storage=zip
{destination}
{core[register_start:register_end]}
gg_step_start() {{ printf 'run:%s\\n' "$1" >> {shlex.quote(str(tmp_path / "ran.txt"))}; }}
gg_step_skip() {{ printf 'skip:%s\\n' "$1" >> {shlex.quote(str(tmp_path / "ran.txt"))}; }}
{stage(core, "Duplication-origin evidence")}
''', encoding="utf-8")

    def run(distance, evidence_dir=""):
        result = subprocess.run(["bash", str(script)], cwd=tmp_path, text=True, capture_output=True,
                                env={**os.environ, "WGD_TEST_PROXIMAL_DISTANCE": str(distance),
                                     "WGD_TEST_EVIDENCE_DIR": evidence_dir}, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
        with GeneFamilyOutputStore(output_root).open_binary("wgd_ssd", result_zip.name) as handle:
            with zipfile.ZipFile(handle) as archive:
                assert archive.testzip() is None
                archive.extractall(tmp_path / "inspection")

    output_root = workspace / "output/query2family"
    result_zip = output_root / "wgd_ssd/OGstage_wgd_ssd.zip"
    legacy = output_root / "wgd_ssd/OGstage/prior.tsv"
    legacy.parent.mkdir(parents=True)
    legacy.write_text("prior results\n")
    results = tmp_path / "inspection/OGstage"
    manifest = workspace / "output/query2family/artifact_provenance/OGstage.wgd_ssd.json"
    run(2)
    origins = read_tsv(results / "duplication_origins.tsv")
    assert len(origins) == 3
    assert {row["family_id"] for row in origins} == {"OGstage"}
    assert sum(row["classification"] == "SSD-supported" for row in origins) == 1
    assert next(row["reason"] for row in origins if row["classification"] == "SSD-supported") == "terminal_tandem_adjacency"
    assert sum(row["classification"] == "unresolved" for row in origins) == 2
    assert {row["event_source"] for row in origins} == {"lca"}
    assert "duplication_origin=SSD-supported" in (results / "classified_gene_tree.nhx").read_text()
    reconciliation = read_tsv(results / "reconciliation.tsv")
    assert {row["gene_name"] for row in reconciliation if row["event_type"] == "leaf"} == set(genes)
    assert {row["mapping_status"] for row in reconciliation} == {"mapped"}
    assert {row["species_name"] for row in reconciliation} == {"A"}
    pairs = {frozenset((row["gene_a"], row["gene_b"])): row for row in read_tsv(results / "pair_evidence.tsv")}
    tandem = frozenset(("sampleA_dup1", "override_dup2"))
    nonadjacent = frozenset(("sampleA_far1", "sampleA_far2"))
    assert pairs[tandem]["position_feature"] == "tandem"
    assert pairs[nonadjacent]["position_feature"] == "distant_same_chromosome"
    assert pairs[nonadjacent]["gene_rank_distance"] == "3"

    provenance = json.loads(manifest.read_text())
    assert provenance["step"] == "wgd_ssd_classification"
    assert provenance["family_id"] == "OGstage"
    assert provenance["parameters"] == {"species_parser": "legacy", "species_regex": regex,
                                        "proximal_distance": "2", "nwkit_identity": identity,
                                        "native_tree_likelihood": "0"}
    sources = {row["label"]: row for row in provenance["inputs"]}
    assert sources["full_species_tree"]["sha256"] == digest(full_tree) != digest(pruned_tree)
    assert sources["full_species_tree"]["path"] == "input/species.full.nwk"
    assert sources["species_map"]["sha256"] == digest(species_map)
    assert sources["genome_evidence"]["path"] == "output/genome_evolution/wgd_ssd"
    assert not any(row["path"].endswith("species.pruned.nwk") for row in provenance["inputs"])
    previous = snapshot(result_zip, manifest)
    # A fresh worker must register cached outputs without publishing them again.
    inventory = family_inventory_path(output_root, "OGstage")
    for journal in inventory.glob("*.paths"):
        journal.unlink()
    run(2)
    assert snapshot(result_zip, manifest) == previous
    assert any(b"wgd_ssd/OGstage_wgd_ssd.zip\0" in journal.read_bytes()
               for journal in inventory.glob("*.paths"))

    relative_evidence = evidence.relative_to(workspace).as_posix()
    run(5, relative_evidence)
    pairs = {frozenset((row["gene_a"], row["gene_b"])): row for row in read_tsv(results / "pair_evidence.tsv")}
    assert pairs[nonadjacent]["position_feature"] == "proximal"
    assert pairs[tandem]["position_feature"] == "tandem"
    assert json.loads(manifest.read_text())["parameters"]["proximal_distance"] == "5"
    assert sum(row["classification"] == "SSD-supported" for row in read_tsv(results / "duplication_origins.tsv")) == 1
    previous = snapshot(result_zip, manifest)
    run(5, relative_evidence)
    assert snapshot(result_zip, manifest) == previous
    assert (tmp_path / "ran.txt").read_text().splitlines() == [
        "run:Duplication-origin evidence", "skip:Duplication-origin evidence (current artifacts)",
        "run:Duplication-origin evidence", "skip:Duplication-origin evidence (current artifacts)",
    ]
    assert not list(temporary.glob("wgd_ssd.*"))
    assert list(family_inventory_path(output_root, "OGstage").glob("*.paths"))
    logical = GeneFamilyOutputStore(output_root)
    logical.mark_family_state("OGstage", "running", "stage-test")
    logical.mark_family_state("OGstage", "complete", "stage-test")
    payload = result_zip.read_bytes()
    manifest_state = (manifest.read_bytes(), manifest.stat().st_mtime_ns)
    enqueue_family_archive(output_root, "query2family", "OGstage", "stage-test")
    drain_archive_queue(output_root, "query2family", query_id_extractor(["OGstage"]))
    assert not result_zip.exists()
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == manifest_state
    assert archive_queue_status(output_root)["pending_families"] == 0
    logical = GeneFamilyOutputStore(output_root)
    assert logical.file_names("wgd_ssd") == [result_zip.name]
    with logical.open_binary("wgd_ssd", result_zip.name) as handle:
        assert handle.read() == payload
    logical.verify()
    # Match the worker's materialize-family prelude before reusing a stage.
    logical.materialize_family("OGstage", query_id_extractor(["OGstage"]))
    run(5, relative_evidence)
    assert result_zip.read_bytes() == payload
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == manifest_state
    assert (tmp_path / "ran.txt").read_text().splitlines()[-1] == "skip:Duplication-origin evidence (current artifacts)"
    assert legacy.read_text() == "prior results\n"
