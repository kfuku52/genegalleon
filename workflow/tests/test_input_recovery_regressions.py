"""Bounded regressions for native input retries; raw scientific inputs stay intact."""

import gzip
import json
import sys
from pathlib import Path

import pandas
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import input_generation_stage_resume as resume
from cds_model_normalisation import CdsModelNormaliser
from format_species_annotation.gff_repair import gene_reference_repairs
from gff2genestat import read_gff_table
from gff_attribute_syntax import normalise_attributes, normalise_line


@pytest.mark.parametrize("note", [
    "Retrovirus-related Pol polyprotein [Includes: Protease (EC 3.4.23.-); Reverse transcriptase (EC 2.7.7.49); Endonuclease]",
    "Alpha-mannosidase [Cleaved into: Alpha-mannosidase, heavy subunit; Alpha-mannosidase, light subunit]",
    "Aquaporin (TIP1;1) (rTIP1)",
    "Cysteine synthase (ARAth-Bsas4;1) (AtcysD1)",
    "Protein [Cleaved into: Subunit alpha; Subunit beta] [Includes: Enzyme A; Enzyme B]",
    "literal%253B (A;B) [Includes: C;D]",
    "Nucleosome assembly protein 1;2 (NtNAP1;2) (Nucleosome assembly protein 1-like 2) (NtNAP1_L2)",
    "Cation transporter HKT1;3 (OsHKT1;3) (Cation transporter HKT6) (OsHKT6)",
])
def test_maker_note_balanced_separators_preserve_description_and_identity(note):
    from urllib.parse import unquote
    raw = "ID=t;Parent=g;Name=t;Note=" + note + ";Dbxref=PFAM:PF00001;transl_table=1"
    expected = raw.replace(note, note.replace(";", "%3B").replace(",", "%2C"))
    fixed = normalise_attributes(raw, "maker", "mRNA")
    assert fixed == expected
    assert unquote(fixed) == unquote(raw)
    assert normalise_attributes(fixed, "maker", "mRNA") == fixed
    changes = []
    normalise_line("chr1\tmaker\tmRNA\t1\t9\t.\t+\t.\t" + raw + "\n", "source", 1, changes)
    assert changes[0]["reason"] == "escaped_maker_note_balanced_semicolon"


@pytest.mark.parametrize("attrs,source,feature", [
    ("ID=t;Parent=g;Note=protein (alpha;Parent:g)", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha;unknown:key)", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha; b;Parent=g)", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha; beta", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha; beta]", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha;; beta)", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein [Includes: alpha; beta];unknown", "maker", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha; beta)", "other", "mRNA"),
    ("ID=t;Parent=g;Note=protein (alpha; beta)", "maker", "gene"),
    ("ID=t;Parent=g (alpha; beta)", "maker", "mRNA"),
    ("Note=Nucleosome assembly protein 1;2 (NtNAP1;3)", "maker", "mRNA"),
    ("Note=Cation transporter HKT1;3 (NtNAP1;3)", "maker", "mRNA"),
    ("Note=unknown protein 1;2 (NtNAP1;2)", "maker", "mRNA"),
])
def test_maker_note_repair_rejects_structural_unknown_and_incomplete_boundaries(attrs, source, feature):
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes(attrs, source, feature)


@pytest.mark.parametrize("mode", ["off", "safe", "strict"])
@pytest.mark.parametrize("paired", [False, True])
def test_maker_note_is_repaired_before_full_gff_formatting_without_model_changes(tmp_path, mode, paired):
    import format_species_discovery as formatting
    from gff_attribute_syntax import validate_gff

    note = "Polyprotein [Includes: Protease; Reverse transcriptase; Endonuclease] (TIP1;1)"
    raw = tmp_path / "raw.gff3"
    text = ("##gff-version 3\nchr1\tmaker\tgene\t1\t12\t.\t+\t.\tID=g1;Name=g1\n"
            "chr1\tmaker\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1;Name=t1;Note=" + note + "\n"
            "chr1\tmaker\tCDS\t1\t12\t.\t+\t0\tParent=t1\n")
    raw.write_text(text)
    cds = tmp_path / "Species_one.cds.fa"
    cds.write_text(">Species_one_g1\nATGAAACCCTAA\n")
    output = tmp_path / "out"
    output.mkdir()
    task = {"provider": "local", "species_key": "Species_one", "species_prefix": "Species_one",
            "gff_path": raw, "gff_repair_mode": mode}
    result = formatting.format_gff(task, output, overwrite=False, dry_run=False,
                                   formatted_cds_path=cds if paired else None)
    validate_gff(result["output_path"])
    with gzip.open(result["output_path"], "rt") as handle:
        assert handle.read() == text.replace(note, note.replace(";", "%3B"))
    audit = json.loads(Path(str(result["output_path"]) + ".repair.json").read_text())["attribute_syntax"]
    assert audit["changed_rows"] == 1 and audit["changes"][0]["source_line"] == 3
    assert audit["changes"][0]["reason"] == "escaped_maker_note_balanced_semicolon"
    assert raw.read_text() == text


@pytest.mark.parametrize("note", [
    "Geneious type: Repeat Region; Geneious type: repeat_region",
    "Geneious type: Repeat_Regiion; IRa; Geneious type: repeat_region",
])
def test_geneious_organelle_repeat_notes_preserve_literal_semicolons(note):
    from urllib.parse import unquote
    raw = "ID=repeat;note=" + note + ";standard_name=IRb Repeat region;transl_table=11"
    fixed = normalise_attributes(raw, ".", "misc_feature")
    assert unquote(fixed) == raw
    assert normalise_attributes(fixed, ".", "misc_feature") == fixed
    changes = []
    normalise_line("organelle\t.\tmisc_feature\t1\t50\t.\t+\t.\t" + raw + "\n", "source", 1, changes)
    assert changes[0]["reason"] == "escaped_geneious_repeat_note_semicolon"


@pytest.mark.parametrize("attrs,source,feature", [
    ("note=Geneious type: Repeat Region;Parent:t", ".", "misc_feature"),
    ("note=Geneious type: Repeat Region; IRa; arbitrary", ".", "misc_feature"),
    ("note=Geneious type: Repeat Region; Geneious type: repeat_region", ".", "gene"),
    ("note=Geneious type: Repeat Region; Geneious type: repeat_region", "other", "misc_feature"),
    ("Note=Geneious type: Repeat Region; Geneious type: repeat_region", ".", "misc_feature"),
])
def test_geneious_repair_rejects_ambiguous_or_structural_fragments(attrs, source, feature):
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes(attrs, source, feature)


def test_description_linkage_separator_is_encoded_and_audited():
    raw = "ID=Sb01g07380;description=PREDICTED: endo-1,3;1,4-beta-D-glucanase-like [Sesamum indicum];transl_table=1"
    expected = raw.replace("1,3;1,4-", "1%2C3%3B1%2C4-")
    assert normalise_attributes(raw, "EVM", "gene") == expected
    assert normalise_attributes(expected, "EVM", "gene") == expected
    changes = []
    line = "chr1\tEVM\tgene\t1\t9\t.\t-\t.\t" + raw + "\n"
    assert normalise_line(line, "source.gff", 23978, changes).endswith(expected + "\n")
    assert changes[0]["reason"] == "escaped_description_linkage_semicolon"


@pytest.mark.parametrize(
    "attrs",
    [
        "ID=g;1,4-beta-D-glucanase",
        "ID=g;Parent=t;1,4-beta-D-glucanase",
        "description=endo-1,3;Parent:t",
        "description=endo-1,3;1,4-beta;extra",
        "description=endo;1,4-beta",
        "description=endo-1,3;ID:t",
    ],
)
def test_description_repair_never_consumes_structural_or_ambiguous_fragments(attrs):
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes(attrs, "EVM", "gene")


@pytest.mark.parametrize("mismatch", ["", "strand", "span", "mixed", "duplicate", "in_bounds", "missing_cds", "child_axis", "child_span", "duplicate_rna"])
def test_out_of_bounds_parent_gene_requires_exact_unanimous_child_span(tmp_path, mismatch):
    (genome, gff) = (tmp_path / "genome.fa", tmp_path / "source.gff")
    genome.write_text(">wrong\n" + "A" * (90 if mismatch == "in_bounds" else 6) + "\n>right\n" + "A" * 90 + "\n")
    gene = "wrong\tEVM\tgene\t10\t30\t.\t-\t.\tID=g\n"
    rna = "right\tEVM\tmRNA\t10\t30\t.\t-\t.\tID=t;Parent=g\n"
    if mismatch == "strand":
        rna = rna.replace("\t-\t", "\t+\t")
    if mismatch == "span":
        rna = rna.replace("\t30\t", "\t29\t")
    text = gene + rna + "right\tEVM\tCDS\t10\t30\t.\t-\t0\tParent=t\n"
    if mismatch == "mixed":
        text += rna.replace("right", "wrong").replace("ID=t", "ID=u")
    if mismatch == "duplicate":
        text += gene
    if mismatch == "duplicate_rna":
        text += rna
    if mismatch == "missing_cds":
        text = gene + rna
    if mismatch == "child_axis":
        text = text.replace("right\tEVM\tCDS", "wrong\tEVM\tCDS")
    if mismatch == "child_span":
        text = text.replace("CDS\t10\t30", "CDS\t10\t31")
    gff.write_text(text)
    repairs = gene_reference_repairs(gff, genome)
    assert bool(repairs) == (mismatch in {"", "in_bounds"})
    if repairs:
        assert repairs["g"]["to_seqid"] == "right"
    assert gff.read_text() == text


@pytest.mark.parametrize("description", [
    "PREDICTED: nucleosome assembly protein 1;4 isoform X2 [Nelumbo nucifera]",
    "RecName: Full=Stamen-specific protein FIL1; Flags: Precursor [Antirrhinum majus]",
    "RecName: Full=Baicalin-beta-D-glucuronidase; AltName: Full=Baicalinase.; Flags: Precursor [Scutellaria baicalensis]",
    "RecName: Full=Elongation factor 1-delta; Short=EF-1-delta; AltName: Full=Elongation factor 1B-beta; AltName: Full=eEF-1B beta [Spuriopimpinella brachycarpa]",
    "RecName: Full=Luminal-binding protein 8; Short=BiP 8; AltName: Full=78 kDa glucose-regulated protein homolog 8; Short=GRP-78-8, partial [Nicotiana tabacum]",
])
def test_publisher_description_metadata_is_losslessly_encoded(description):
    from urllib.parse import unquote

    raw = "ID=g;description=" + description + ";transl_table=1"
    fixed = normalise_attributes(raw, ".", "gene")
    assert unquote(fixed) == raw
    assert normalise_attributes(fixed, ".", "gene") == fixed
    assert fixed.startswith("ID=g;description=") and fixed.endswith(";transl_table=1")
    changes = []
    normalise_line("chr1\t.\tgene\t1\t9\t.\t+\t.\t" + raw + "\n", "source", 1, changes)
    assert changes[0]["reason"] == "escaped_description_metadata_semicolon"


@pytest.mark.parametrize("fragment", ["Parent:t", "AltName: Full=x;Parent:t", " Flags: Arbitrary", "4 isoform X2;extra"])
def test_description_metadata_never_swallows_unknown_or_structural_fields(fragment):
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes("ID=g;description=RecName: Full=x;" + fragment, ".", "gene")


def annotation(alignment_count=0):
    return (
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g;pseudo=true\nchr1\ts\tmatch\t1\t9\t.\t+\t.\tID=bridge;Parent=g\nchr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=bridge\nchr1\ts\texon\t1\t9\t.\t+\t.\tParent=t\nchr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t\n"
        + "".join(
            (f"chr1\ts\tmatch_part\t1\t9\t.\t+\t.\tID=alignment{i};Parent=other{i}\n" for i in range(alignment_count))
        )
    )


def test_coding_normalisation_retains_complete_ancestor_closure_and_same_decision(tmp_path):
    gff = tmp_path / "source.gff"
    gff.write_text(annotation(3000))
    args = dict(
        source=dict(gff=str(gff), genome="unused", species="S", genetic_code=1), directory=tmp_path, side="format"
    )
    (full, bounded) = (CdsModelNormaliser(**args), CdsModelNormaliser(**args, coding_only=True))
    assert len(full.features) == 3005 and len(bounded.features) == 5
    assert full.ancestry(["t"]) == bounded.ancestry(["t"]) == {"t", "bridge", "g"}
    assert full.models == bounded.models and full.exons == bounded.exons
    full.normalise("seq", "ATGAAATAA", {"t"})
    bounded.normalise("seq", "ATGAAATAA", {"t"})
    assert full.rows == bounded.rows and full.corrected == bounded.corrected


@pytest.mark.parametrize("compressed", [False, True])
def test_coding_gff_table_retains_referenced_alignment_ancestors(tmp_path, compressed):
    gff = tmp_path / ("source.gff" + (".gz" if compressed else ""))
    text = annotation(3000)
    if compressed:
        with gzip.open(gff, "wt") as handle:
            handle.write(text)
    else:
        gff.write_text(text)
    full = read_gff_table(str(gff))
    bounded = read_gff_table(str(gff), coding_only=True)
    assert len(bounded) == 5
    pandas.testing.assert_frame_equal(
        full.iloc[:5].reset_index(drop=True), bounded.sort_index().iloc[[0, 4, 1, 2, 3]].reset_index(drop=True)
    )


@pytest.fixture
def busco_pair(tmp_path):
    (donor, target) = (tmp_path / "old/output/input_generation", tmp_path / "new/output/input_generation")
    settings = []
    for root in (donor, target):
        (root / "tmp").mkdir(parents=True)
        (root / "tmp/busco_lineage.resolved.txt").write_text("embryophyta_odb12\n")
        item = dict(run_species_busco="1", busco_lineage="embryophyta_odb12")
        for kind in ("full", "short"):
            directory = root / ("busco_" + kind)
            directory.mkdir()
            item["species_busco_" + kind + "_dir"] = str(directory)
        settings.append(item)
    (source_cds, target_cds) = (donor / "cds.fa", target / "cds.fa")
    for f in (source_cds, target_cds):
        f.write_text(">S_g\nATGAAATAA\n")
    outputs = []
    for kind in ("full", "short"):
        f = Path(settings[0]["species_busco_" + kind + "_dir"]) / (
            "S.busco." + ("full.tsv" if kind == "full" else "short.txt")
        )
        f.write_text("BUSCO validated result " + kind + "\n")
        outputs.append(dict(label="busco_" + kind, path=str(f), sha256=resume.digest(f)))
    payload = dict(
        schema_version=1,
        step="input_generation_species_busco",
        family_id="S",
        parameters=dict(
            busco_lineage_request="embryophyta_odb12",
            busco_lineage_resolved="embryophyta_odb12",
            busco_mode="transcriptome",
            evalue="1e-03",
            limit="20",
        ),
        inputs=[dict(label="species_cds", path=str(source_cds), sha256=resume.digest(source_cds))],
        outputs=outputs,
    )
    (donor / "artifact_provenance").mkdir()
    (donor / "artifact_provenance/busco.S.json").write_text(json.dumps(payload))
    return (donor, target, source_cds, target_cds, *settings)


def test_identical_busco_import_writes_native_target_provenance(busco_pair):
    assert resume.import_busco(*busco_pair[:2], "S", *busco_pair[2:])
    payload = json.loads((busco_pair[1] / "artifact_provenance/busco.S.json").read_text())
    assert payload["inputs"][0]["path"] == "output/input_generation/cds.fa"
    assert all((x["path"].startswith("output/input_generation/busco_") for x in payload["outputs"]))


@pytest.mark.parametrize("change", ["lineage", "resolved", "cds", "output", "parameter"])
def test_busco_import_rejects_changed_scientific_contract(busco_pair, change):
    (donor, target, source_cds, target_cds, source_settings, target_settings) = busco_pair
    if change == "lineage":
        target_settings["busco_lineage"] = "viridiplantae_odb12"
    elif change == "resolved":
        (target / "tmp/busco_lineage.resolved.txt").write_text("viridiplantae_odb12\n")
    elif change == "cds":
        target_cds.write_text("changed")
    elif change == "output":
        next(Path(source_settings["species_busco_full_dir"]).iterdir()).write_text("changed")
    else:
        path = donor / "artifact_provenance/busco.S.json"
        payload = json.loads(path.read_text())
        payload["parameters"]["limit"] = "21"
        path.write_text(json.dumps(payload))
    if change in ("lineage", "resolved"):
        assert not resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    else:
        with pytest.raises(ValueError):
            resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    assert not (target / "artifact_provenance/busco.S.json").exists()


def test_parent_reference_repair_runs_through_standard_formatter_and_strict_validator(tmp_path):
    from format_species_annotation.reference import validate_gff_genome_references
    from format_species_inputs import format_cds, format_gff

    (genome, gff) = (tmp_path / "genome.fa", tmp_path / "source.gff")
    genome.write_text(">wrong\nATGAAA\n>right\nATGAAATAA\n")
    raw = "wrong\tEVM\tgene\t1\t9\t.\t+\t.\tID=g\nright\tEVM\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\nright\tEVM\tCDS\t1\t9\t.\t+\t0\tParent=t\n"
    gff.write_text(raw)
    task = dict(
        provider="local",
        species_prefix="Species_one",
        species_key="Species_one",
        gff_path=gff,
        genome_path=genome,
        gff_repair_mode="safe",
        gene_grouping_mode="rescue_overlap",
        format_strict=True,
    )
    cds = format_cds(task, tmp_path, False, False)
    formatted = format_gff(task, tmp_path, False, False, formatted_cds_path=cds["output_path"])
    assert validate_gff_genome_references(formatted["output_path"], genome) == 3
    with gzip.open(formatted["output_path"], "rt") as h:
        assert h.read() == raw.replace("wrong\t", "right\t")
    audit = json.loads(Path(str(formatted["output_path"]) + ".repair.json").read_text())
    assert audit["gene_reference_repairs"]["g"]["from_seqid"] == "wrong"
    assert gff.read_text() == raw


def test_busco_import_rejects_a_source_change_during_copy(busco_pair, monkeypatch):
    copy = resume.copy_atomic
    (donor, target, source_cds, target_cds, source_settings, target_settings) = busco_pair

    def racing_copy(source, destination, **kwargs):
        copy(source, destination, **kwargs)
        source_cds.write_text("raced source")

    monkeypatch.setattr(resume, "copy_atomic", racing_copy)
    with pytest.raises(OSError, match="File changed"):
        resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    assert not (target / "artifact_provenance/busco.S.json").exists()


@pytest.mark.parametrize("extra", ["\t", "\t\t\t"])
@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_empty_trailing_gff_columns_preserve_features_and_audit(tmp_path, extra, ending):
    from format_species_annotation.organelle import iter_non_organelle_gff_lines
    from gff_attribute_syntax import validated_lines

    feature = "Chr10\tEVM\tCDS\t1029939\t1030655\t.\t-\t0\tID=cds.CcoG0000001.1;Parent=CcoG0000001.1;"
    raw = feature + extra + ending
    changes = []
    fixed = normalise_line(raw, "Cornus.gff", 1, changes)
    assert fixed == feature + ending
    assert changes == [{"source_line": 1, "before": feature + extra,
                        "after": feature, "reason": "removed_empty_trailing_columns"}]
    assert list(validated_lines([fixed], "formatted.gff")) == [fixed]
    with pytest.raises(ValueError, match="nine GFF columns"):
        list(validated_lines([raw], "unformatted.gff"))
    source = tmp_path / "Cornus.gff"
    source.write_bytes(raw.encode())
    assert list(iter_non_organelle_gff_lines(source, seqids=())) == [feature + ("\n" if ending else "")]
    assert source.read_bytes() == raw.encode()


@pytest.mark.parametrize("extra", ["\tvalue", "\t ", "\tvalue\t", "\t\tvalue"])
def test_nonempty_extra_gff_columns_remain_rejected(extra):
    feature = "chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t"
    with pytest.raises(ValueError, match="nine GFF columns"):
        normalise_line(feature + extra + "\n", "source.gff", 1, [])
