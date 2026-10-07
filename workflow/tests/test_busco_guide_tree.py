"""BUSCO protein preservation and input-integrity contracts."""
import gzip
import inspect
import json
from pathlib import Path

import pytest
from Bio import Phylo

from workflow.support import busco_guide_tree as guide
from workflow.support import busco_reference_quality as quality
from workflow.support import pairwise_synteny, rescue_gene_models
from workflow.support.rescue_gene_models import patristic_distances


def test_existing_quality_and_identifier_imports_keep_the_same_implementations():
    assert rescue_gene_models.COMPARABLE_QUALITY == guide.COMPARABLE_QUALITY == quality.COMPARABLE_QUALITY
    for name, consumers in (("busco_quality", (rescue_gene_models, guide)),
                            ("patristic_distances", (rescue_gene_models, guide)),
                            ("safe_token", (pairwise_synteny, rescue_gene_models, guide))):
        for consumer in consumers:
            function = getattr(consumer, name)
            assert inspect.getsource(function) == inspect.getsource(getattr(quality, name))
            assert Path(function.__globals__["__file__"]).resolve() == Path(quality.__file__).resolve()


def make_archive(tmp_path, name, sequences, duplicated=()):
    for folder in ("cds", "full", "short", "full/single_copy"):
        (tmp_path / folder).mkdir(parents=True, exist_ok=True)
    cds = tmp_path / "cds" / (name + ".fa")
    cds.write_text(">" + name + "_gene\nATGAAATAA\n")
    full = tmp_path / "full" / (name + ".busco.full.tsv")
    short = tmp_path / "short" / (name + ".busco.short.txt")
    run = tmp_path / "raw" / name / "run_test_odb12"
    proteins = run / "busco_sequences/single_copy_busco_sequences"
    proteins.mkdir(parents=True)
    rows = []
    for marker, sequence in sequences.items():
        rows.append(f"{marker}\tComplete\t{name}_gene:24-100\t1\t100\n")
        (proteins / (marker + ".faa")).write_text(f">{name}_gene:24-100 original\n{sequence}\n")
    for marker in duplicated:
        rows += [f"{marker}\tDuplicated\td1\n", f"{marker}\tDuplicated\td2\n"]
    full.write_text("# fixture\n" + "".join(rows))
    short.write_text("# BUSCO version is: fixture\n"
                     "# The lineage dataset is: test_odb12 (Creation date: 2026-01-01)\n"
                     "# BUSCO was run in mode: transcriptome\n"
                     f"C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:{len(sequences)+len(duplicated)}\n")
    archive = tmp_path / "full/single_copy" / (name + ".json.gz")
    args = guide.parser().parse_args(["preserve", "--run-dir", str(run), "--full", str(full),
                                      "--short", str(short), "--input", str(cds), "--species", name,
                                      "--output", str(archive)])
    guide.preserve(args)
    return args


def test_only_sole_complete_rows_are_single_copy(tmp_path):
    path = tmp_path / "full.tsv"
    path.write_text("# fixture\nA\tComplete\tid\nB\tComplete\tid\nB\tComplete\tid2\n"
                    "C\tDuplicated\tid\nD\tFragmented\tid\nE\tMissing\n")
    rows, universe = guide.single_copy_rows(path)
    assert rows == {"A": "id"} and universe == set("ABCDE")


def test_preserve_exact_busco_protein_and_mapping_deterministically(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK*", "BUSCO2": "MXXDEK"}, ["BUSCO3"])
    before = args.output.read_bytes()
    guide.preserve(args)
    assert args.output.read_bytes() == before
    records, _ = guide.read_archive(args.output, args.species, args.full, args.short, args.input)
    assert set(records) == {"BUSCO1", "BUSCO2"}
    assert records["BUSCO1"]["sequence"] == "MDEK"
    assert records["BUSCO1"]["busco_sequence_id"] == "Plant_example_gene:24-100"
    assert records["BUSCO2"]["sequence"] == "MXXDEK"


@pytest.mark.parametrize("packed", [False, True])
def test_resume_imports_verified_single_copy_archive_and_legacy_tables(tmp_path, packed):
    from workflow.support import input_generation_stage_resume as resume
    name = "Plant_example"
    args = make_archive(tmp_path / "source-data", name, {"BUSCO1": "MDEK"})
    source_root, target_root = tmp_path / "source-output", tmp_path / "target-output"
    for root in (source_root, target_root):
        (root / "tmp").mkdir(parents=True)
        (root / "tmp/busco_lineage.resolved.txt").write_text("test_odb12\n")
    (source_root / "artifact_provenance").mkdir()
    target_cds = tmp_path / "target-data" / args.input.name
    target_cds.parent.mkdir()
    target_cds.write_bytes(args.input.read_bytes())
    settings = dict(run_species_busco="1", busco_lineage="test_odb12",
                    species_busco_full_dir=str(args.full.parent), species_busco_short_dir=str(args.short.parent))
    target_settings = {**settings, "species_busco_full_dir": str(tmp_path / "target-full"),
                       "species_busco_short_dir": str(tmp_path / "target-short")}
    entries = [("busco_full", args.full), ("busco_short", args.short)]
    if packed:
        entries.append(("busco_single_copy", args.output))
    parameters = dict(busco_lineage_request="test_odb12", busco_lineage_resolved="test_odb12",
                      busco_mode="transcriptome", evalue="1e-03", limit="20")
    payload = dict(schema_version=1, step="input_generation_species_busco", family_id=name,
                   inputs=[dict(label="species_cds", sha256=guide.digest(args.input))],
                   outputs=[dict(label=label, sha256=guide.digest(path)) for label, path in entries],
                   parameters=parameters)
    (source_root / "artifact_provenance" / ("busco."+name+".json")).write_text(json.dumps(payload))
    assert resume.import_busco(source_root, target_root, name, args.input, target_cds, settings, target_settings)
    destination = tmp_path / "target-full/single_copy" / args.output.name
    assert destination.exists() == packed
    if packed:
        assert destination.read_bytes() == args.output.read_bytes()
        guide.read_archive(destination, name, tmp_path / "target-full" / args.full.name,
                           tmp_path / "target-short" / args.short.name, target_cds)
        args.output.write_bytes(args.output.read_bytes()+b"changed")
        with pytest.raises(ValueError, match="BUSCO output/input differs"):
            resume.import_busco(source_root, target_root, name, args.input, target_cds, settings, target_settings)


def test_changed_cds_and_full_table_are_rejected(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    args.input.write_text(">changed\nATGCCCTAA\n")
    with pytest.raises(ValueError, match="inputs changed"):
        guide.read_archive(args.output, args.species, args.full, args.short, args.input)


def test_internal_stop_and_duplicate_proteins_fail_preservation(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    protein = args.run_dir / "busco_sequences/single_copy_busco_sequences/BUSCO1.faa"
    protein.write_text(">a\nMD*EK\n")
    with pytest.raises(ValueError, match="Invalid BUSCO protein"):
        guide.preserve(args)
    protein.write_text(">a\nMDEK\n>b\nMDEK\n")
    with pytest.raises(ValueError, match="one record"):
        guide.preserve(args)


def test_flat_busco_layout(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    # The core passes run_LINEAGE; BUSCO 6 may publish directly under its root.
    args.run_dir = args.run_dir / "run_test_odb12"
    guide.preserve(args)
    assert args.output.exists()


def test_patristic_matrix_matches_biopython_including_polytomies(tmp_path):
    tree_path = tmp_path / "tree.nwk"
    tree_path.write_text("((A:0.2,B:0.7):0.3,C:0,D:0.8):10;")
    tree = Phylo.read(tree_path, "newick")
    distances = patristic_distances(tree)
    for a in distances:
        for b in distances:
            assert distances[a][b] == pytest.approx(tree.distance(a,b))


def test_corrupt_cache_metadata_is_a_miss(tmp_path):
    path, meta = tmp_path / "cache.bin", tmp_path / "cache.json"
    path.write_bytes(b"cache")
    meta.write_text("{bad")
    assert not guide.cache_matches(path, meta)
    meta.write_text(json.dumps({"sha256": guide.digest(path)}))
    assert guide.cache_matches(path, meta)
    path.write_bytes(b"changed")
    assert not guide.cache_matches(path, meta)


def test_archive_record_tampering_is_rejected(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    with gzip.open(args.output, "rt") as handle:
        data = json.load(handle)
    data["records"]["BUSCO1"]["busco_sequence_id"] = "different"
    with gzip.open(args.output, "wt") as handle:
        json.dump(data, handle)
    with pytest.raises(ValueError, match="Invalid archived"):
        guide.read_archive(args.output, args.species, args.full, args.short, args.input)


def test_preserve_rejects_protein_from_another_busco_match(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    protein = args.run_dir / "busco_sequences/single_copy_busco_sequences/BUSCO1.faa"
    protein.write_text(">Plant_other_gene:24-100\nMDEK\n")
    with pytest.raises(ValueError, match="protein ID"):
        guide.preserve(args)


def test_preserve_accepts_busco_metaeuk_wrapped_target_id(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    protein = args.run_dir / "busco_sequences/single_copy_busco_sequences/BUSCO1.faa"
    protein.write_text(">BUSCO1_ancestral|Plant_example_gene:24-100|+\nMDEK\n")
    guide.preserve(args)
    records, _ = guide.read_archive(args.output, args.species, args.full, args.short, args.input)
    assert records["BUSCO1"]["protein_id"] == "BUSCO1_ancestral|Plant_example_gene:24-100|+"


def test_preserve_accepts_busco_transcriptome_strand_target_id(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    protein = args.run_dir / "busco_sequences/single_copy_busco_sequences/BUSCO1.faa"
    protein.write_text(">Plant_example_gene:24-100|+ orig_seq_frame_1\nMDEK\n")
    guide.preserve(args)
    records, _ = guide.read_archive(args.output, args.species, args.full, args.short, args.input)
    assert records["BUSCO1"]["protein_id"] == "Plant_example_gene:24-100|+"


def test_preserve_cannot_replace_its_cds_input(tmp_path):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    original = args.input.read_bytes()
    args.output = args.input
    with pytest.raises(ValueError, match="overlap"):
        guide.preserve(args)
    assert args.input.read_bytes() == original


def test_preserve_detects_protein_change_during_read(tmp_path, monkeypatch):
    args = make_archive(tmp_path, "Plant_example", {"BUSCO1": "MDEK"})
    before = args.output.read_bytes()
    original = guide.fasta_records
    def changed_after_read(path):
        yield from original(path)
        path.write_text(path.read_text().replace("MDEK", "MDDD"))
    monkeypatch.setattr(guide, "fasta_records", changed_after_read)
    with pytest.raises(OSError, match="File changed"):
        guide.preserve(args)
    assert args.output.read_bytes() == before
