"""Bounded real-tool rescue tests; all biological inputs are synthetic and temporary."""
import copy
import gzip
import hashlib
import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

from workflow.support import rescue_gene_models as rescue
from workflow.support.rescue_anchor_admission import prepare_rescue_genome

SCRIPT = Path(rescue.__file__)


def test_load_refuses_plan_replaced_while_validating_tools(tmp_path, monkeypatch):
    plan = {"request": {"schema": rescue.SCHEMA, "tools": {}, "parameters": {"minimum_coverage": .95}}}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    def changing_identities():
        new = copy.deepcopy(plan)
        new["request"]["parameters"]["minimum_coverage"] = .5
        path.write_text(json.dumps(new))
        return {}
    monkeypatch.setattr(rescue, "identities", changing_identities)
    with pytest.raises(ValueError, match="Frozen rescue plan changed"):
        rescue.load(tmp_path)


def test_stale_plan_cannot_be_used_to_stamp_prepared_annotation(hidden_models):
    output, names, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    changed = copy.deepcopy(plan)
    changed["request"]["sources"][names[0]]["genetic_code"] = 4
    rescue.atomic_json(output / "plan.json", changed)
    with pytest.raises(ValueError, match="Frozen rescue plan changed"):
        rescue.prepared(output, plan, names[0])
    assert not (output / "prepared" / names[0] / "receipt.json").exists()


def anchor_fixture(tmp_path, models, code=1):
    """Real exon coordinates on both strands, with untouched source records."""
    fasta, genome, features = [], [], ["##gff-version 3"]
    for number, model in enumerate(models):
        identifier, contig = f"g{number}", f"chr{number}"
        blocks = model.get("blocks", [(model.get("cds", "ATGAAATAA"), 0)])
        utr5, utr3 = model.get("utr5", ""), model.get("utr3", "")
        parts = [sequence for sequence, _ in blocks]
        parts[0] = utr5 + parts[0]
        parts[-1] += utr3
        dna = "NNNNNNN".join(parts)
        strand = model.get("strand", "+")
        if strand == "-":
            dna = str(Seq(dna).reverse_complement())
        genome.append(f">{contig}\n{dna}\n")
        def row(kind, start, end, phase, attr, strand=strand, dna=dna, contig=contig):
            if strand == "-":
                start, end = len(dna) - end, len(dna) - start
            return f"{contig}\tsynthetic\t{kind}\t{start + 1}\t{end}\t.\t{strand}\t{phase}\t{attr}"
        features.extend([row("gene", 0, len(dna), ".", "ID=" + identifier),
                         row("mRNA", 0, len(dna), ".", f"ID={identifier}.t1;Parent={identifier}" + model.get("attributes", ""))])
        offset = 0
        for index, ((sequence, phase), part) in enumerate(zip(blocks, parts, strict=True)):
            left = offset + (len(utr5) if index == 0 else 0)
            features.append(row("CDS", left, left + len(sequence), phase, f"Parent={identifier}.t1"))
            if not model.get("no_exons"):
                features.append(row("exon", offset, offset + len(part), ".", f"Parent={identifier}.t1"))
            offset += len(part) + 7
        supplied = model.get("supplied", "".join(sequence for sequence, _ in blocks))
        fasta.append(f">Plant_example_{identifier} original header\n{supplied}\n")
    paths = {key: tmp_path / (key + suffix) for key, suffix in (("fasta", ".fa"), ("gff", ".gff3"), ("genome", ".fa"))}
    for key, contents in (("fasta", fasta), ("genome", genome), ("gff", ["\n".join(features) + "\n"])):
        paths[key].write_text("".join(contents))
    output = tmp_path / "prepared"
    output.mkdir()
    return {"species": "Plant_example", "mode": "cds", "feature": "gene", "attribute": "ID",
            "genetic_code": code, **{key: str(path) for key, path in paths.items()}}, output


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("convention,phase", [("gff3", 2), ("complementary", 1)])
def test_admission_utr_partial_and_original_preservation(tmp_path, strand, convention, phase):
    coding, partial = "ATGAAACCCTAA", "TAATGAAACCCTAA"
    second = (phase + 4 if convention == "complementary" else phase - 4) % 3
    source, output = anchor_fixture(tmp_path, [
        {"cds": coding},
        {"cds": coding, "utr5": "TAA", "supplied": "TAA" + coding, "strand": strand},
        {"blocks": [(partial[:4], phase), (partial[4:], second)], "strand": strand},
        {"cds": "ATGTGACCCTAA"}], code=1)
    before = {key: Path(source[key]).read_bytes() for key in ("fasta", "gff", "genome")}
    genes, metadata = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 3
    assert metadata["anchor_admission"]["counts"] == {"unchanged": 1, "normalised": 2, "excluded": 1}
    assert {key: Path(source[key]).read_bytes() for key in before} == before
    assert not list(tmp_path.glob("*.fai")) and not list(output.glob(".anchor-genome-*"))
    proteins = {identifier: sequence for identifier, _, sequence in rescue.fasta_records(output / "genes.pep")}
    assert proteins["Plant_example_g1"] == proteins["Plant_example_g2"] == "MKP"
    assert "Plant_example_g3" not in proteins
    mapping = rescue.table(output / "genes.id_map.tsv")
    assert mapping[-1]["status"] == "translation_excluded"
    audit = json.loads((output / "genes.anchor_admission.json").read_text())
    assert audit["records"][1]["selected_evidence"][0]["phase_convention"] == convention


@pytest.mark.parametrize("attributes,reason", [
    (";transl_except=(pos:4..6,aa:OTHER)", "annotated_translation_exception"),
    (";exception=unclassified transcription discrepancy", "annotated_translation_exception"),
    (";pseudo=true", "annotated_pseudogene")])
def test_admission_preserves_and_withholds_annotation_exceptions(tmp_path, attributes, reason):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGTGACCCTAA", "attributes": attributes}])
    genes, _ = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 1
    record = json.loads((output / "genes.anchor_admission.json").read_text())["records"][0]
    assert record["reason"] == reason


@pytest.mark.parametrize("case", ["unannotated_utr", "mismatched_utr", "ambiguous_phase", "inconsistent_phase", "mixed_file_phase"])
def test_admission_requires_independent_model_evidence(tmp_path, case):
    model = {"cds": "ATGAAATAA", "utr5": "TAA", "supplied": "TAAATGAAATAA"}
    models = [{}, model]
    if case == "unannotated_utr":
        model["no_exons"] = True
    elif case == "mismatched_utr":
        model["supplied"] = "TAACCCATGAAATAA"
    elif case == "ambiguous_phase":
        model.clear()
        model["blocks"] = [("TAATGAAATAA", 1)]
    elif case == "inconsistent_phase":
        model.clear()
        model["blocks"] = [("TAAT", 1), ("GAAATAA", 1)]
    else:
        models = [{"blocks": [("ATGA", 0), ("AATAA", 2)]},
                  {"blocks": [("ATGA", 0), ("AATAA", 1)]}, {"blocks": [("TAATGAAATAA", 1)]}]
    genes, metadata = prepare_rescue_genome(*anchor_fixture(tmp_path, models), "genes", 1.0)
    assert len(genes) == len(models) - 1
    assert metadata["anchor_admission"]["counts"]["excluded"] == 1


def test_admission_genetic_code_and_strict_new_models(tmp_path):
    source, output = anchor_fixture(tmp_path, [{"cds": "ATGTGACCCTAA"},
                                               {"cds": "ATGAAATAA", "utr5": "TAA", "supplied": "TAAATGAAATAA"}], code=4)
    genes, _ = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 2
    assert "MWP" in (output / "genes.pep").read_text()
    with pytest.raises(ValueError, match="strict original-translation"):
        prepare_rescue_genome(source, output, "genes", 1.0, required_ids=["Plant_example_g1"])
    with pytest.raises(ValueError, match="internal stop"):
        rescue.prepare_genome(source, output, "strict", 1.0)


def test_admission_all_excluded_keeps_failure_audit(tmp_path):
    source, output = anchor_fixture(tmp_path, [{"cds": "ATGTGACCCTAA"}])
    with pytest.raises(ValueError, match="No usable rescue anchors"):
        prepare_rescue_genome(source, output, "genes", 1.0)
    assert json.loads((output / "genes.anchor_admission.json").read_text())["counts"] == {"excluded": 1}


def test_admission_ambiguous_isoforms_are_withheld(tmp_path):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGAAATAA", "utr5": "TAA", "supplied": "TAAATGAAATAA"}])
    gff = Path(source["gff"])
    lines = ["chr1\tsynthetic\tmRNA\t1\t12\t.\t+\t.\tID=g1.t2;Parent=g1",
             "chr1\tsynthetic\texon\t1\t12\t.\t+\t.\tParent=g1.t2",
             "chr1\tsynthetic\tCDS\t7\t12\t.\t+\t0\tParent=g1.t2"]
    gff.write_text(gff.read_text() + "\n".join(lines) + "\n")
    genes, _ = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 1
    assert json.loads((output / "genes.anchor_admission.json").read_text())["records"][0]["reason"] == "ambiguous_genomic_translations"


def test_admission_does_not_relabel_disrupted_cds_as_another_isoforms_utr(tmp_path):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGTGACCCTAA"}])
    gff = Path(source["gff"])
    gff.write_text(gff.read_text() + "chr1\tsynthetic\tmRNA\t1\t12\t.\t+\t.\tID=g1.t2;Parent=g1\n"
                  "chr1\tsynthetic\texon\t1\t12\t.\t+\t.\tParent=g1.t2\n"
                  "chr1\tsynthetic\tCDS\t7\t12\t.\t+\t0\tParent=g1.t2\n")
    genes, _ = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 1
    record = json.loads((output / "genes.anchor_admission.json").read_text())["records"][0]
    assert record["reason"] == "genomic_internal_stop_or_incomplete_codon" and record["bound_transcripts"] == ["g1.t1"]


def test_admission_removes_only_proven_formatter_padding(tmp_path):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGAAATAA", "supplied": "ATGAAATAANN"},
                                              {"cds": "ATGAAATAA", "supplied": "ATGAAATAANNN"}])
    genes, metadata = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 2
    assert metadata["anchor_admission"]["reasons"]["formatter_terminal_padding"] == 1


@pytest.mark.parametrize("attribute,source_id,formatted_id", [("Dbxref", "GeneID:123", "GeneID123"),
                                                            ("Accession", "GWHGAASQ123", "GWHGAASQ123")])
def test_admission_maps_genegalleon_provider_gene_identifiers(tmp_path, attribute, source_id, formatted_id):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGAAATAA", "utr5": "TAA", "supplied": "TAAATGAAATAA"}])
    fasta, gff = Path(source["fasta"]), Path(source["gff"])
    fasta.write_text(fasta.read_text().replace("_g0 ", "_GeneID999 " if attribute == "Dbxref" else "_GWHGAASQ999 ")
                    .replace("_g1 ", "_" + formatted_id + " "))
    gff.write_text(gff.read_text().replace("ID=g0\n", "ID=g0;" + attribute + "=" + ("GeneID:999" if attribute == "Dbxref" else "GWHGAASQ999") + "\n")
                  .replace("ID=g1\n", "ID=g1;" + attribute + "=" + source_id + "\n"))
    source.update(feature="", attribute="")
    genes, metadata = prepare_rescue_genome(source, output, "genes", 1.0)
    assert len(genes) == 2 and metadata["attribute"] == attribute
    assert genes[1].gene_id == "Plant_example_" + formatted_id
    # Explicit source mapping follows the same canonical-ID contract.
    source.update(feature="gene", attribute=attribute)
    assert len(prepare_rescue_genome(source, output, "explicit", 1.0)[0]) == 2


def test_admission_source_phase_evidence_is_not_changed_by_new_rescued_models(tmp_path):
    source, output = anchor_fixture(tmp_path, [
        {"blocks": [("ATGA", 0), ("AATAA", 1)]},
        {"blocks": [("TAATGAAATAA", 1)]},
        {"blocks": [("ATGA", 0), ("AATAA", 2)]}])
    gff = Path(source["gff"])
    gff.write_text("\n".join(line.replace("\tsynthetic\t", "\tgenegalleon_rescue\t") if "g2" in line else line
                             for line in gff.read_text().splitlines()) + "\n")
    genes, metadata = prepare_rescue_genome(source, output, "genes", 1.0, required_ids=["Plant_example_g2"])
    assert len(genes) == 3 and metadata["anchor_admission"]["file_phase_convention"] == "complementary"
    assert metadata["anchor_admission"]["counts"] == {"unchanged": 2, "normalised": 1}


def test_gemoma_withholds_complementary_phase_reference_without_launching_java(tmp_path, monkeypatch):
    from types import SimpleNamespace

    donor, target = "Donor_species", "Target_species"
    prepared = tmp_path / "prepared" / donor
    prepared.mkdir(parents=True)
    rescue.write_tsv(prepared / "genes.id_map.tsv", ("original_id", "jcvi_id", "locus_id", "status"),
                     [(donor + "_g1", donor + "_g1", "g1", "selected")])
    evidence = [{"phase_convention": "complementary", "offset": 2}]
    (prepared / "genes.anchor_admission.json").write_text(json.dumps({"records": [
        {"original_id": donor + "_g1", "status": "normalised", "selected_evidence": evidence}]}))
    gff, jar, java = tmp_path / "donor.gff3", tmp_path / "test.jar", tmp_path / "java"
    gff.write_text("chr1\tsynthetic\tmRNA\t1\t10\t.\t+\t.\tID=t1;Parent=g1\n")
    jar.write_text("synthetic jar")
    java.write_text("synthetic executable")
    plan = {"request": {"gemoma_jar": str(jar), "gemoma_java": str(java),
                        "files": {str(jar): rescue.digest(jar), str(java): rescue.digest(java)}, "parameters": {},
                        "sources": {donor: {"gff": str(gff), "genome": str(tmp_path / "donor.fa")}}}}
    monkeypatch.setattr(rescue, "verify_sources", lambda *_: None)
    def unexpected(*_):
        pytest.fail("GeMoMa must not launch on complementary source phases")
    monkeypatch.setattr(rescue, "run", unexpected)
    region = {"donor": donor, "query": donor + "_g1", "id": "query"}
    rescue.refine_gemoma(tmp_path, tmp_path, plan, {"species": target}, [region], SimpleNamespace(), [], 1)
    assert json.loads((tmp_path / ("gemoma_" + donor) / "query.skipped.json").read_text())["reason"] == "unsupported_reference_phase_convention"


@pytest.mark.parametrize("invalid", ["coordinates", "duplicate_contig"])
def test_admission_rejects_invalid_genome_and_leaves_no_source_index(tmp_path, invalid):
    source, output = anchor_fixture(tmp_path, [{}, {"cds": "ATGTGACCCTAA"}])
    genome = Path(source["genome"])
    genome.write_text(genome.read_text() + ">chr1\nA\n" if invalid == "duplicate_contig" else ">chr0\nA\n>chr1\nA\n")
    with pytest.raises(ValueError, match="Annotation.*genome|FASTA index warning"):
        prepare_rescue_genome(source, output, "genes", 1.0)
    assert not list(tmp_path.glob("*.fai")) and not list(output.glob(".anchor-genome-*"))


def cli(*args):
    # Existing saved-array fixtures continue to exercise the explicit public
    # legacy mode. Compact/default storage has a full parity test below.
    if args[0] == "plan" and "--model-storage" not in args:
        args = (*args, "--model-storage", "legacy", "--retain-search-inputs")
    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_rescue_cli_help_in_owned_runtime(tmp_path):
    result = subprocess.run([sys.executable, str(SCRIPT), "--help"], cwd=tmp_path,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower()
    assert list(tmp_path.iterdir()) == []


@pytest.fixture
def hidden_models(tmp_path):
    rng = random.Random(919)
    codons = ["GCT", "CGT", "AAC", "GAC", "TGC", "CAA", "GAA", "GGT", "CAC", "ATT",
              "CTG", "AAA", "ATG", "TTC", "CCT", "TCT", "ACT", "TGG", "TAC", "GTT"]
    cds = ["ATG" + "".join(rng.choice(codons) for _ in range(159)) + "TAA" for _ in range(18)]
    species = [f"Plant_species{i}" for i in range(6)]
    for directory in ("cds", "gff", "genome", "busco"):
        (tmp_path / directory).mkdir()
    for name in species:
        genome = ""
        fasta, features = [], ["##gff-version 3", "# Note: ##FASTA is a directive only on its own line"]
        for i, sequence in enumerate(cds):
            start = len(genome)
            encoded = sequence
            if name == species[0] and i == 9:
                encoded = sequence[:240] + "TGA" + sequence[243:]
            if name == species[0] and i == 10:
                encoded = sequence[:240] + "A" + sequence[240:]
            genome += encoded + "N" * 60
            # Hide intact g8, stop-disrupted g9 and frameshift-disrupted g10.
            if name == species[0] and i in {8, 9, 10}:
                continue
            identifier = name + f"_g{i}"
            fasta.append(f">{identifier}\n{sequence}\n")
            end = start + len(sequence)
            features += [f"chr1\tsynthetic\tgene\t{start + 1}\t{end}\t.\t+\t.\tID=g{i}",
                         f"chr1\tsynthetic\tmRNA\t{start + 1}\t{end}\t.\t+\t.\tID=g{i}.t1;Parent=g{i}",
                         f"chr1\tsynthetic\tCDS\t{start + 1}\t{end}\t.\t+\t0\tParent=g{i}.t1"]
        (tmp_path / "cds" / (name + ".cds.fa")).write_text("".join(fasta))
        (tmp_path / "gff" / (name + ".gff3")).write_text("\n".join(features) + "\n")
        (tmp_path / "genome" / (name + ".genome.fa")).write_text(">chr1\n" + genome + "\n")
        (tmp_path / "busco" / (name + ".busco.short.txt")).write_text(
            "# BUSCO version is: 6.0.0\n# The lineage dataset is: embryophyta_odb12\n"
            "# BUSCO was run in mode: transcriptome\nC:95.0%[S:95.0%,D:0.0%],F:0.0%,M:5.0%,n:100\n")
    (tmp_path / "tree.nwk").write_text("((Plant_species0:1,Plant_species1:1):1,"
                                      "(Plant_species2:1,Plant_species3:1):1,(Plant_species4:1,Plant_species5:1):1);\n")
    return tmp_path, species, cds


def make_plan(fixture, *, model_storage="legacy", retain_search_inputs=True, directory_name="rescue"):
    root, species, cds = fixture
    output = root / directory_name
    storage_args = ["--model-storage", model_storage]
    if retain_search_inputs:
        storage_args.append("--retain-search-inputs")
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", output, *storage_args)
    return output, species, cds


def test_compact_pipeline_preserves_every_legacy_model_and_export(hidden_models):
    from workflow.support.rescue_model_store import (
        iter_models,
        iter_partial_models,
        iter_revision_models,
    )
    legacy, names, _ = make_plan(hidden_models, directory_name="legacy_storage")
    compact, _, _ = make_plan(hidden_models, model_storage="compact", retain_search_inputs=False,
                              directory_name="compact_storage")
    for output in (legacy, compact):
        cli("run", "--output", output, "--cpus", "4", "--interval-workers", "4")
    for name in names:
        a, b = legacy / "rescued" / name, compact / "rescued" / name
        assert list(iter_models(b)) == json.loads((a / "models.json").read_text())
        assert list(iter_partial_models(b)) == json.loads((a / "partial_models.json").read_text())
        assert list(iter_revision_models(b)) == json.loads((a / "revision_candidates.json").read_text())
        assert (b / "model_store/manifest.json").is_file()
        assert not any((b / filename).exists() for filename in (
            "models.json", "partial_models.json", "revision_candidates.json", "regions.fa", "queries.fa",
            "genome.gff", "unresolved.fa", "unresolved.unique.fa", "genome.covered.unique.fa"))
        assert not list((b / "intervals").glob("*/region.fa"))
        assert not list((b / "intervals").glob("*/queries.fa"))
        for folder, suffix in (("species_cds", ".rescue.cds.fa"), ("species_gff", ".rescue.gff3")):
            assert (legacy / "augmented" / folder / (name + suffix)).read_bytes() == (
                compact / "augmented" / folder / (name + suffix)).read_bytes()
    target = names[0]
    old_receipt = (compact / "rescued" / target / "receipt.json").read_bytes()
    exported = compact.parent / "legacy_export"
    cli("export-models", "--output", compact, "--task-index", "1", "--destination", exported, "--cpus", "4")
    for filename in ("models.json", "partial_models.json", "revision_candidates.json"):
        assert json.loads((exported / filename).read_text()) == json.loads((legacy / "rescued" / target / filename).read_text())
    diagnostics = compact.parent / "diagnostic_export"
    cli("export-search-inputs", "--output", compact, "--task-index", "1", "--destination", diagnostics,
        "--combined", "--cpus", "4")
    for path in (legacy / "rescued" / target / "intervals").glob("*/*.fa"):
        assert (diagnostics / "intervals" / path.parent.name / path.name).read_bytes() == path.read_bytes()
    for filename in ("regions.fa", "queries.fa", "unresolved.unique.fa", "unresolved.fa", "genome.gff"):
        path = legacy / "rescued" / target / filename
        if path.is_file():
            assert (diagnostics / filename).read_bytes() == path.read_bytes()
    assert (compact / "rescued" / target / "receipt.json").read_bytes() == old_receipt
    # Restart reuses verified publications rather than running a new prediction.
    result = cli("rescue", "--output", compact, "--task-index", "1", "--cpus", "4")
    assert "Reused rescued/" in result.stdout
    assert (compact / "rescued" / target / "receipt.json").read_bytes() == old_receipt


def test_parallel_intervals_preserve_every_real_alignment_and_order(hidden_models):
    import pysam
    root, names, sequences = hidden_models
    genome_path = root / "genome" / (names[0] + ".genome.fa")
    pysam.faidx(str(genome_path))
    donor = names[1]
    proteins = {donor: {donor + f"_g{i}": str(Seq(sequence).translate())[:-1] for i, sequence in enumerate(sequences)}}
    windows = {}
    for i in [8, 9, 10, 11, 12, 13]:
        start = i * (len(sequences[0]) + 60)
        windows[("chr1", start, start + len(sequences[i]) + 10)] = [
            {"id": f"query_{i}", "donor": donor, "query": donor + f"_g{i}"}]
    serial, parallel = root / "serial", root / "parallel"
    serial.mkdir()
    parallel.mkdir()
    with pysam.FastaFile(str(genome_path)) as genome:
        expected = rescue.search_intervals(serial, windows, proteins, genome, 1, 20000, 4, 1)
        actual = rescue.search_intervals(parallel, windows, proteins, genome, 1, 20000, 4)
    assert actual == expected
    assert actual and len({m["query"] for m in actual}) >= 3
    for i in range(1, len(windows) + 1):
        relative = Path("intervals") / str(i)
        assert (serial / relative / "models.gff").read_bytes() == (parallel / relative / "models.gff").read_bytes()
        command = json.loads((parallel / relative / "logs/miniprot.command.json").read_text())
        assert command[command.index("-t") + 1] == "1"


def test_local_input_diagnostics_fetch_each_new_window_once(tmp_path):
    class Genome:
        fetched = []
        def fetch(self, seqid, start, end):
            self.fetched.append((seqid, start, end))
            return "ACGT"[start:end]
    genome = Genome()
    proteins = {"donor": {"g1": "MK", "g2": "MP"}}
    regions = [{"id": "q1", "donor": "donor", "query": "g1"},
               {"id": "q2", "donor": "donor", "query": "g2"}]
    rescue.write_local_search_inputs(tmp_path, {}, proteins, genome)
    assert list(tmp_path.iterdir()) == [] and genome.fetched == []
    rescue.write_local_search_inputs(tmp_path, {("chr1", 1, 4): regions}, proteins, genome)
    assert genome.fetched == [("chr1", 1, 4)]
    assert (tmp_path / "regions.fa").read_text() == ">q1\nCGT\n>q2\nCGT\n"
    assert (tmp_path / "queries.fa").read_text() == ">q1\nMK\n>q2\nMP\n"


def test_partial_local_cache_writes_only_new_queries_and_replays_empty_results(hidden_models, monkeypatch):
    fixture, names, sequences = hidden_models
    original, _, _ = make_plan(hidden_models)
    plan = rescue.load(original)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    plan["request"]["parameters"]["genome_fallback"] = 0
    rescue.atomic_json(original / "plan.json", plan)
    old = []
    for number, query in enumerate((8, 0)):
        start = 8 * (len(sequences[0]) + 60) if number == 0 else len(sequences[0])
        end = start + len(sequences[query]) if number == 0 else start + 60
        old.append({"id": "old_" + str(number), "donor": names[1], "query": names[1] + f"_g{query}",
                    "seqid": "chr1", "start": start, "end": end,
                    "expected_start": start, "expected_end": end})
    start = 11 * (len(sequences[0]) + 60) + 1
    new = {"id": "new_local", "donor": names[1], "query": names[1] + "_g11", "seqid": "chr1",
           "start": start, "end": start + len(sequences[11]),
           "expected_start": start, "expected_end": start + len(sequences[11])}
    monkeypatch.setattr(rescue, "candidates", lambda root, *_: copy.deepcopy(old if root == original else [*old, new]))
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    first = rescue.rescue(original, plan, names[0], 2)
    assert "old_1" not in {m["query"] for m in json.loads((first / "models.json").read_text())}

    def derivative(label, previous=None):
        output = fixture / label
        output.mkdir()
        current = copy.deepcopy(plan)
        if previous:
            current["request"]["prediction_cache"] = rescue.frozen_prediction_cache_key(previous, [names[0]])
        rescue.atomic_json(output / "plan.json", current)
        return output, rescue.load(output)
    fresh_root, fresh_plan = derivative("fresh")
    fresh = rescue.rescue(fresh_root, fresh_plan, names[0], 2)
    partial_root, partial_plan = derivative("partial", original)
    partial = rescue.rescue(partial_root, partial_plan, names[0], 2)
    assert [i for i, _, _ in rescue.fasta_records(fresh / "queries.fa")] == ["old_0", "old_1", "new_local"]
    assert [i for i, _, _ in rescue.fasta_records(partial / "queries.fa")] == ["new_local"]
    assert [i for i, _, _ in rescue.fasta_records(partial / "regions.fa")] == ["new_local"]
    expected = json.loads((fresh / "models.json").read_text())
    fields = ("query", "cds", "sequence", "problems", "status", "evidence", "coverage", "identity",
              "terminal_completion", "model_id", "quality_evidence", "partial_evidence", "raw_prediction",
              "search", "start", "end")
    def scientific_rows(directory):
        rows = []
        for model in json.loads((directory / "models.json").read_text()):
            row = {k: model.get(k) for k in fields}
            row["raw_prediction"] = {k: v for k, v in row["raw_prediction"].items() if k != "id"}
            rows.append(json.dumps(row, sort_keys=True))
        return sorted(rows)
    assert expected and scientific_rows(partial) == scientific_rows(fresh)
    assert json.loads((partial / "prediction_reuse.json").read_text())["searched_local_windows"] == 1
    replay_root, replay_plan = derivative("replay", partial_root)
    def no_prediction(*args, **kwargs):
        raise AssertionError("Verified local alignment or empty result was searched again")
    with monkeypatch.context() as patch:
        patch.setattr(rescue, "run", no_prediction)
        replay = rescue.rescue(replay_root, replay_plan, names[0], 2)
    assert scientific_rows(replay) == scientific_rows(fresh)
    assert not (replay / "queries.fa").exists() and not (replay / "regions.fa").exists()
    assert not (replay / "genome.fa").exists() and not (replay / "genome.mpi").exists()
    reuse = json.loads((replay / "prediction_reuse.json").read_text())
    assert reuse["cached_local_queries"] == 3 and reuse["searched_local_windows"] == 0
    assert [r["id"] for r in json.loads((replay / "candidates.json").read_text())] == ["old_0", "old_1", "new_local"]


def test_genome_only_search_uses_fallback_inputs_without_local_diagnostics(hidden_models, monkeypatch):
    output, names, sequences = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    rescue.atomic_json(output / "plan.json", plan)
    region = {"id": "genome_only", "target": names[0], "donor": names[1], "query": names[1] + "_g8",
              "comparison": "synthetic", "genome_only": True, "placement": "unanchored",
              "nomination": {"reason": "no_target_match"}, "orthology": "unassigned", "expected_copy": "unassigned"}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [])
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([region], {"nominated": 1}))
    directory = rescue.rescue(output, plan, names[0], 2)
    assert not (directory / "regions.fa").exists() and not (directory / "queries.fa").exists()
    assert [i for i, _, _ in rescue.fasta_records(directory / "unresolved.unique.fa")] == [region["id"]]
    assert [row["candidate"] for row in rescue.table(directory / "genome_query_mapping.tsv")] == [region["id"]]
    models = json.loads((directory / "models.json").read_text())
    assert any(m["search"] == "genome_fallback" and m["sequence"] == sequences[8] for m in models)
    assert all(m["status"] != "accepted" and "unanchored_genome_search" in m["problems"] for m in models)


def test_broader_genome_contract_researches_cached_positive_and_empty_coverage(hidden_models, monkeypatch):
    """Real miniprot, strict QC and two cache generations, without old decisions."""
    fixture, names, sequences = hidden_models
    original, _, _ = make_plan(hidden_models)
    plan = rescue.load(original)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    old_contract = rescue.prediction_search_contract()
    old_contract["genome"]["output_score_ratio"] = .99
    plan["request"]["tools"]["prediction_search_contract"] = old_contract
    rescue.atomic_json(original / "plan.json", plan)
    # Local empty search; the hidden intact copy is elsewhere on this genome.
    region = {"id": "negative_local", "donor": names[1], "query": names[1] + "_g8",
              "seqid": "chr1", "start": len(sequences[0]), "end": len(sequences[0]) + 60,
              "expected_start": len(sequences[0]), "expected_end": len(sequences[0]) + 60}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [copy.deepcopy(region)])
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    original_run = rescue.run
    def narrower(command, directory, label, stdout=None):
        if label == "miniprot_genome":
            command = ["--outs=0.99" if str(arg).startswith("--outs=") else arg for arg in command]
        original_run(command, directory, label, stdout)
    with monkeypatch.context() as patch:
        patch.setattr(rescue, "prediction_search_contract", lambda: copy.deepcopy(old_contract))
        patch.setattr(rescue, "run", narrower)
        narrow = rescue.rescue(original, plan, names[0], 2)
    assert json.loads((narrow / "models.json").read_text())  # A raw positive was produced.
    assert (narrow / "genome_query_mapping.tsv").exists()  # The negative result was also fully covered.
    assert not any(m["search"] == "synteny_interval" for m in json.loads((narrow / "models.json").read_text()))

    def derivative(label, previous=None):
        output = fixture / label
        output.mkdir()
        current = copy.deepcopy(plan)
        current["request"]["tools"]["prediction_search_contract"] = rescue.prediction_search_contract()
        if previous:
            current["request"]["prediction_cache"] = rescue.frozen_prediction_cache_key(previous, [names[0]])
        rescue.atomic_json(output / "plan.json", current)
        return output, rescue.load(output)
    fresh_root, fresh_plan = derivative("broader_fresh")
    fresh = rescue.rescue(fresh_root, fresh_plan, names[0], 2)
    current_root, current_plan = derivative("broader_cached", original)
    calls = []
    def counted(command, directory, label, stdout=None):
        calls.append((label, list(map(str, command))))
        original_run(command, directory, label, stdout)
    with monkeypatch.context() as patch:
        patch.setattr(rescue, "run", counted)
        current = rescue.rescue(current_root, current_plan, names[0], 2)
    assert [label for label, _ in calls] == ["miniprot_index", "miniprot_genome"]
    command = calls[-1][1]
    assert "--outs=0.5" in command and command[command.index("-N") + 1] == "30"
    reuse = json.loads((current / "prediction_reuse.json").read_text())
    assert reuse["cached_local_queries"] == 1 and reuse["searched_local_windows"] == 0
    assert reuse["cached_genome_queries"] == 0 and reuse["cached_genome_search_compatible"] is False
    assert not (current / "regions.fa").exists() and not (current / "queries.fa").exists()
    fields = ("query", "cds", "sequence", "problems", "status", "evidence", "coverage", "identity",
              "model_id", "quality_evidence", "search", "start", "end")
    def scientific_rows(directory):
        return sorted(json.dumps({k: m.get(k) for k in fields}, sort_keys=True)
                      for m in json.loads((directory / "models.json").read_text()))
    assert scientific_rows(current) == scientific_rows(fresh)
    replay_root, replay_plan = derivative("broader_replay", current_root)
    def no_search(*args, **kwargs):
        raise AssertionError("The completed .5 search was repeated")
    with monkeypatch.context() as patch:
        patch.setattr(rescue, "run", no_search)
        replay = rescue.rescue(replay_root, replay_plan, names[0], 2)
    assert scientific_rows(replay) == scientific_rows(fresh)
    assert json.loads((replay / "prediction_reuse.json").read_text())["cached_genome_queries"] == 1


def test_genome_query_dedup_restores_exact_real_gff_and_each_candidate(hidden_models):
    root, names, sequences = hidden_models
    genome = root / "genome" / (names[0] + ".genome.fa")
    protein = str(Seq(sequences[8]).translate())[:-1]
    proteins = {names[1]: {"intact": protein, "unmapped": "W" * 160, "other": str(Seq(sequences[9]).translate())[:-1]},
                names[2]: {"same": protein}}
    specs = [("first", names[1], "intact"), ("unmapped", names[1], "unmapped"),
             ("outside_copy", names[2], "same"), ("other", names[1], "other"), ("again", names[1], "intact")]
    regions = [{"id": identifier, "donor": donor, "query": query} for identifier, donor, query in specs]
    full, unique = root / "full.fa", root / "unique.fa"
    full.write_text("".join(f">{r['id']}\n{proteins[r['donor']][r['query']]}\n" for r in regions))
    mapping = rescue.write_unique_queries(regions, proteins, unique)
    assert len(list(rescue.fasta_records(unique))) == 3 and len(mapping) == 5
    baseline, dedup, expanded = root / "full.gff", root / "unique.gff", root / "expanded.gff"
    for query, out, label in [(full, baseline, "full"), (unique, dedup, "unique")]:
        rescue.run(["miniprot", "-u", "-t", 2, "-G", 20000, "--gff", genome, query], root, label, out)
    rescue.expand_miniprot_queries(dedup, expanded, mapping)
    assert expanded.read_bytes() == baseline.read_bytes()
    models = rescue.read_miniprot(expanded)
    assert models == rescue.read_miniprot(baseline)
    assert {m["query"] for m in models} >= {"first", "outside_copy", "again", "other"}
    # Identical proteins retain their different expected intervals and therefore
    # cannot turn an alignment to a paralog elsewhere into an accepted model.
    first = next(m for m in models if m["query"] == "first")
    outside = next(m for m in models if m["query"] == "outside_copy")
    for model in [first, outside]:
        model.update(start=min(e[0] for e in model["cds"]), end=max(e[1] for e in model["cds"]), problems=[])
    first["evidence"] = {"seqid": first["seqid"], "expected_start": first["start"], "expected_end": first["end"]}
    outside["evidence"] = {"seqid": outside["seqid"], "expected_start": 0, "expected_end": 1}
    assert rescue.check_interval(first)["problems"] == []
    assert "outside_expected_synteny_interval" in rescue.check_interval(outside)["problems"]


def test_query_expansion_refuses_missing_or_unknown_representatives(tmp_path):
    source, dest = tmp_path / "raw.gff", tmp_path / "expanded.gff"
    source.write_text("##gff-version 3\n")
    with pytest.raises(ValueError, match="every representative"):
        rescue.expand_miniprot_queries(source, dest, {"original": "known"})
    source.write_text("##PAF\tunknown\t100\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n")
    with pytest.raises(ValueError, match="unknown representative"):
        rescue.expand_miniprot_queries(source, dest, {"original": "known"})
    assert not dest.exists()


def benchmark_evidence(fixture, intervals=2, fallback=True):
    """Generate real miniprot evidence with the producer's frozen file/plan binding."""
    root, names, sequences = fixture
    evidence = root / "benchmark_evidence"
    source = evidence / "rescued" / names[0]
    source.mkdir(parents=True)
    genome = root / "genome" / (names[0] + ".genome.fa")
    plan = {"request": {"sources": {names[0]: {"genome": str(genome), "genetic_code": 1}},
                        "files": {str(genome): rescue.digest(genome)}, "parameters": {"max_intron": 20000},
                        "tools": {"prediction_search_contract": rescue.prediction_search_contract()}}}
    rescue.atomic_json(evidence / "plan.json", plan)
    rescue.atomic_json(source / rescue.SEARCH_CONTRACT_FILE, rescue.prediction_search_contract())
    dna = next(rescue.fasta_records(genome))[2]
    class Genome:
        def fetch(self, _, start, end):
            return dna[start:end]
    proteins = {"donor": {f"query{i}": str(Seq(sequences[i]).translate())[:-1] for i in [8, 9]}}
    windows = {}
    for i in [8, 9][:intervals]:
        start = i * (len(sequences[0]) + 60)
        windows[("chr1", start, start + len(sequences[i]))] = [
            {"id": f"query{i}", "donor": "donor", "query": f"query{i}"}]
    rescue.search_intervals(source, windows, proteins, Genome(), 1, 20000, 2, retain_inputs=True)
    if fallback:
        queries = source / "unresolved.fa"
        queries.write_text(">query9\n" + proteins["donor"]["query9"] + "\n")
        rescue.run(["miniprot", "-u", "-t", 2, "-G", 20000,
                    "-N", rescue.MINIPROT_MAX_SECONDARY, f"--outs={rescue.MINIPROT_OUTPUT_SCORE_RATIO}",
                    "--gff", genome, queries],
                   source, "map", source / "genome.gff")
    files = {str(p.relative_to(source)): rescue.digest(p) for p in source.rglob("*") if p.is_file()}
    rescue.atomic_json(source / "receipt.json", {"key": {"plan": rescue.digest(evidence / "plan.json"),
                                                        "species": names[0]}, "files": files})
    return evidence, source, names[0]


def run_benchmark(evidence, species, output, check_existing=True, optimized_python=False):
    script = SCRIPT.parents[1] / "benchmarks/benchmark_rescue_search.py"
    return subprocess.run([sys.executable, *(["-O"] if optimized_python else []), str(script),
                           "--evidence", str(evidence), "--species", species,
                           "--output", str(output), "--cpus", "2", "--repeats", "1",
                           *(["--check-existing"] if check_existing else [])], capture_output=True, text=True, timeout=120)


@pytest.mark.parametrize("check_existing", [False, True])
def test_benchmark_verifies_both_real_search_phases(hidden_models, check_existing):
    evidence, _, species = benchmark_evidence(hidden_models)
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output, check_existing)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((output / "result.json").read_text())
    assert report["outputs_identical"] and not report["skipped"]
    assert report["original_intervals"] == report["intervals"] == 2
    assert report["genome_queries"] == report["unique_queries"] == 1
    assert report["warmups"] == (0 if check_existing else 1)


@pytest.mark.parametrize("intervals", [0, 2])
@pytest.mark.parametrize("check_existing", [False, True])
def test_benchmark_accepts_absent_search_phases(hidden_models, intervals, check_existing):
    evidence, _, species = benchmark_evidence(hidden_models, intervals=intervals, fallback=False)
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output, check_existing)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((output / "result.json").read_text())
    assert report["outputs_identical"] and report["genome_queries"] == report["unique_queries"] == 0
    assert report["samples"]["genome_unique"] == [] and "genome" in report["skipped"]
    assert not (output / "genome_fixture").exists()
    if not intervals:
        assert report["samples"]["interval_parallel"] == [] and "interval" in report["skipped"]


def test_benchmark_refuses_a_missing_original_interval(hidden_models):
    import shutil
    evidence, source, species = benchmark_evidence(hidden_models)
    shutil.rmtree(source / "intervals/2")
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output)
    assert result.returncode != 0 and "interval" in result.stderr.lower(), result.stdout + result.stderr
    assert not (output / "result.json").exists()


@pytest.mark.parametrize("missing", ["unresolved.fa", "genome.gff"])
def test_benchmark_refuses_partial_fallback_evidence(hidden_models, missing):
    evidence, source, species = benchmark_evidence(hidden_models)
    receipt = json.loads((source / "receipt.json").read_text())
    del receipt["files"][missing]
    rescue.atomic_json(source / "receipt.json", receipt)
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output)
    assert result.returncode != 0 and "fallback evidence" in result.stderr.lower(), result.stdout + result.stderr
    assert not (output / "result.json").exists()


def test_benchmark_equivalence_checks_survive_python_optimization(hidden_models):
    evidence, source, species = benchmark_evidence(hidden_models)
    (source / "genome.gff").write_text("# Different producer alignment\n")
    receipt = json.loads((source / "receipt.json").read_text())
    receipt["files"]["genome.gff"] = rescue.digest(source / "genome.gff")
    rescue.atomic_json(source / "receipt.json", receipt)
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output, optimized_python=True)
    assert result.returncode != 0 and "Full genome evidence differs" in result.stderr, result.stdout + result.stderr
    assert not (output / "result.json").exists()


@pytest.mark.parametrize("field", ["plan", "species"])
def test_benchmark_refuses_a_foreign_producer_receipt(hidden_models, field):
    evidence, source, species = benchmark_evidence(hidden_models)
    receipt = json.loads((source / "receipt.json").read_text())
    receipt["key"][field] = "foreign"
    rescue.atomic_json(source / "receipt.json", receipt)
    output = hidden_models[0] / "bench"
    result = run_benchmark(evidence, species, output)
    assert result.returncode != 0 and "receipt" in result.stderr.lower(), result.stdout + result.stderr
    assert not (output / "result.json").exists()


def test_benchmark_rechecks_fallback_source_after_mapping(hidden_models, monkeypatch):
    import importlib.util
    evidence, source, species = benchmark_evidence(hidden_models)
    script = SCRIPT.parents[1] / "benchmarks/benchmark_rescue_search.py"
    spec = importlib.util.spec_from_file_location("rescue_search_benchmark_audit", script)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    original = benchmark.rescue.run
    def mutate_source(command, directory, label, stdout=None):
        original(command, directory, label, stdout)
        if label == "map":
            with (source / "unresolved.fa").open("a") as handle:
                handle.write(">changed_after_mapping\nMWP\n")
    monkeypatch.setattr(benchmark.rescue, "run", mutate_source)
    output = hidden_models[0] / "bench"
    monkeypatch.setattr(sys, "argv", [str(script), "--evidence", str(evidence), "--species", species,
                                     "--output", str(output), "--cpus", "2", "--check-existing"])
    with pytest.raises(ValueError, match="source changed"):
        benchmark.main()
    assert not (output / "result.json").exists()


def test_comparison_cache_reuses_metadata_changes_and_binds_actual_inputs(hidden_models):
    root, names, _ = hidden_models
    first, _, _ = make_plan(hidden_models)
    plan = rescue.load(first)
    job = next(j for j in plan["synteny_jobs"] if j["a"] == names[0] and j["b"] == names[1])
    output = rescue.synteny(first, plan, job["index"], 1)
    original = (output / "blocks.json").read_bytes()
    cache = root / "gene_model_rescue_comparison_cache"
    receipts = {f: f.stat().st_mtime_ns for f in cache.glob("comparisons/*/receipt.json")}
    assert len(receipts) == 1
    gff = root / "gff" / (names[0] + ".gff3")
    gff.write_text(gff.read_text().replace("ID=g0\n", "ID=g0;Name=metadata_only\n"))
    second = root / "second"
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", second)
    new_plan = rescue.load(second)
    new_job = next(j for j in new_plan["synteny_jobs"] if j["a"] == names[0] and j["b"] == names[1])
    reused = rescue.synteny(second, new_plan, new_job["index"], 2)
    assert rescue.digest(first / "plan.json") != rescue.digest(second / "plan.json")
    assert (reused / "blocks.json").read_bytes() == original
    assert receipts == {f: f.stat().st_mtime_ns for f in cache.glob("comparisons/*/receipt.json")}
    assert "--cpus=1" in json.loads((reused / "logs/pair.command.json").read_text())
    entry = next(cache.glob("comparisons/*"))
    (entry / "blocks.json").write_text("corrupted cache")
    assert (reused / "blocks.json").read_bytes() == original  # No shared writable inode.
    third = root / "third"
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", third)
    third_plan = rescue.load(third)
    repaired = rescue.synteny(third, third_plan, new_job["index"], 2)
    assert (repaired / "blocks.json").read_bytes() == original
    assert rescue.verified(entry, rescue.comparison_cache_key(third, third_plan, new_job))
    key = rescue.comparison_cache_key(second, new_plan, new_job)
    assert key == rescue.comparison_cache_key(second, new_plan, {**new_job, "id": "renumbered", "index": 123})
    changed = copy.deepcopy(new_plan)
    changed["request"]["parameters"]["minimum_coverage"] = 0.99
    assert key == rescue.comparison_cache_key(second, changed, new_job)
    changed["request"]["parameters"]["cscore"] = 0.8
    assert key != rescue.comparison_cache_key(second, changed, new_job)
    for tool in ("diamond_sha256", "lastdb_sha256", "numpy", "python"):
        changed = copy.deepcopy(new_plan)
        changed["request"]["tools"][tool] = "changed_tool"
        assert key != rescue.comparison_cache_key(second, changed, new_job)
    changed = copy.deepcopy(new_plan)
    changed["request"]["tools"]["source_hashes"]["jcvi.files/0/algorithms/lis.py"] = "changed_support"
    assert key != rescue.comparison_cache_key(second, changed, new_job)
    protein = second / "prepared" / names[0] / "genes.pep"
    protein.write_text(protein.read_text() + "\n")
    assert key != rescue.comparison_cache_key(second, new_plan, new_job)
    with pytest.raises(ValueError, match="corrupted"):
        rescue.comparison_key(second, new_job)


def test_concurrent_plans_share_one_verified_comparison(hidden_models, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    root, names, _ = hidden_models
    first, _, _ = make_plan(hidden_models)
    first_plan = rescue.load(first)
    first_job = next(j for j in first_plan["synteny_jobs"] if j["a"] == names[0] and j["b"] == names[1])
    rescue.prepared(first, first_plan, names[0])
    rescue.prepared(first, first_plan, names[1])
    second = root / "second"
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", second, "--minimum-coverage", 0.99)
    second_plan = rescue.load(second)
    second_job = next(j for j in second_plan["synteny_jobs"] if j["a"] == names[0] and j["b"] == names[1])
    rescue.prepared(second, second_plan, names[0])
    rescue.prepared(second, second_plan, names[1])
    builds, original = [], rescue.build_comparison
    def counted(*args):
        builds.append(args[0])
        original(*args)
    monkeypatch.setattr(rescue, "build_comparison", counted)
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = [pool.submit(rescue.synteny, directory, plan, job["index"], 1)
                   for directory, plan, job in [(first, first_plan, first_job), (second, second_plan, second_job)]]
        outputs = [future.result(timeout=120) for future in pending]
    assert len(builds) == 1
    assert (outputs[0] / "blocks.json").read_bytes() == (outputs[1] / "blocks.json").read_bytes()
    for directory, plan, job, output in zip([first, second], [first_plan, second_plan], [first_job, second_job], outputs, strict=True):
        assert rescue.verified(output, rescue.comparison_key(directory, job))


def test_rescue_exports_invalid_originals_and_audits_while_adding_intact_model(hidden_models, monkeypatch):
    root, names, sequences = hidden_models
    name = names[0]
    # Two existing annotations cannot translate directly: one contains a UTR,
    # the other has a genuine internal stop. Neither original is overwritten.
    cds = root / "cds" / (name + ".cds.fa")
    gff = root / "gff" / (name + ".gff3")
    genome = root / "genome" / (name + ".genome.fa")
    intact, disrupted = sequences[0], "ATGTGACCCTAA"
    with cds.open("a") as handle:
        handle.write(f">{name}_utr source UTR\nTAA{intact}\n>{name}_stop source disruption\n{disrupted}\n")
    with genome.open("a") as handle:
        handle.write(f">utr_contig\nTAA{intact}\n>stop_contig\n{disrupted}\n")
    with gff.open("a") as handle:
        for contig, identifier, length, first in (("utr_contig", "utr", len(intact) + 3, 4), ("stop_contig", "stop", len(disrupted), 1)):
            handle.write(f"{contig}\tsynthetic\tgene\t1\t{length}\t.\t+\t.\tID={identifier}\n"
                         f"{contig}\tsynthetic\tmRNA\t1\t{length}\t.\t+\t.\tID={identifier}.t;Parent={identifier}\n"
                         f"{contig}\tsynthetic\texon\t1\t{length}\t.\t+\t.\tParent={identifier}.t\n"
                         f"{contig}\tsynthetic\tCDS\t{first}\t{length}\t.\t+\t0\tParent={identifier}.t\n")
    original_records = list(rescue.fasta_records(cds))
    original_gff = gff.read_bytes()
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][name] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    start = 8 * (len(sequences[0]) + 60)
    monkeypatch.setattr(rescue, "candidates", lambda *_: [{"id": "query", "donor": names[1], "query": names[1] + "_g8",
                       "seqid": "chr1", "start": start, "end": start + len(sequences[8]),
                       "expected_start": start, "expected_end": start + len(sequences[8])}])
    rescue.rescue(output, plan, name, 1)
    exported = rescue.finalize(output, plan, [name])
    records = list(rescue.fasta_records(exported / "species_cds" / (name + ".rescue.cds.fa")))
    assert records[:len(original_records)] == original_records and len(records) == len(original_records) + 1
    assert (exported / "species_gff" / (name + ".rescue.gff3")).read_bytes().startswith(original_gff)
    summary = json.loads((exported / "summary.json").read_text())
    assert summary["rescued_models"] == 1
    assert summary["anchor_admission"][name]["counts"] == {"unchanged": 16, "normalised": 1, "excluded": 1}
    assert (exported / "anchor_admission" / (name + ".tsv")).is_file()
    assert name + "_stop" not in (output / "prepared" / name / "genes.pep").read_text()
    assert not list((root / "genome").glob("*.fai"))


def test_reference_plan_sparse_and_frozen(hidden_models):
    output, species, _ = make_plan(hidden_models)
    plan = json.loads((output / "plan.json").read_text())
    assert len(plan["common_references"]) == 5
    assert len([j for j in plan["synteny_jobs"] if j["kind"] == "self"]) == 6
    pairs = [(j["a"], j["b"]) for j in plan["synteny_jobs"] if j["kind"] == "pair"]
    assert len(pairs) == len(set(pairs))
    assert all(n not in plan["donors"][n] for n in species)
    first = (output / "plan.json").read_bytes()
    make_plan(hidden_models)
    assert (output / "plan.json").read_bytes() == first
    with (hidden_models[0] / "cds" / (species[0] + ".cds.fa")).open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="changed"):
        rescue.verify_sources(plan, [species[0]], ["fasta"])


def test_root_stem_is_not_treated_as_pairwise_branch_length(hidden_models):
    root, species, _ = hidden_models
    (root / "tree.nwk").write_text("((" + ",".join(species[:3]) + "),(" + ",".join(species[3:]) + ")):9;\n")
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    assert plan["request"]["tree_metric"] == "unit_edges"
    assert plan["nearest_references"][species[0]][0] == species[1]


def test_partial_branch_lengths_do_not_silently_mean_zero_distance(hidden_models):
    root, names, _ = hidden_models
    (root / "tree.nwk").write_text(f"(({names[0]}:1,{names[1]}):1,({names[2]}:1,{names[3]}:1):1,"
                                   f"({names[4]}:1,{names[5]}:1):1);\n")
    with pytest.raises(AssertionError, match="Partial branch lengths"):
        make_plan(hidden_models)


@pytest.mark.parametrize("invalid", ["unknown_contig", "outside_genome", "duplicate_contig"])
def test_invalid_genome_annotation_pair_fails_even_without_candidates(hidden_models, monkeypatch, invalid):
    root, names, _ = hidden_models
    path = root / "genome" / (names[0] + ".genome.fa")
    if invalid == "unknown_contig":
        path.write_text(path.read_text().replace(">chr1", ">different_assembly"))
    elif invalid == "outside_genome":
        path.write_text(">chr1\nATGAAATAA\n")
    else:
        with path.open("a") as handle:
            handle.write(">chr1\nATGCCCTAA\n")
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    monkeypatch.setattr(rescue, "candidates", lambda *_: [])
    with pytest.raises(ValueError, match="annotation.*genome|FASTA index warning"):
        rescue.rescue(output, plan, names[0], 1)
    assert not (output / "rescued" / names[0] / "receipt.json").exists()


@pytest.mark.parametrize("contig", ["001", "NA", "NULL", "nan"])
def test_compressed_inputs_and_literal_contigs_recover_and_export(hidden_models, monkeypatch, contig):
    root, names, sequences = hidden_models
    # Paths containing spaces and contigs resembling missing/numeric values
    # must survive discovery, annotation mapping, indexing and export literally.
    relocated = root / "inputs with spaces"
    relocated.mkdir()
    for name in ("cds", "gff", "genome", "busco", "tree.nwk"):
        (root / name).rename(relocated / name)
    root = relocated
    for subdir in ("cds", "gff", "genome"):
        for path in (root / subdir).iterdir():
            contents = path.read_bytes().replace(b"chr1", contig.encode())
            with gzip.open(str(path) + ".gz", "wb") as handle:
                handle.write(contents)
            path.unlink()
    output, _, _ = make_plan((root, names, sequences))
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    start = 8 * (len(sequences[0]) + 60)
    region = {"id": "query", "donor": names[1], "query": names[1] + "_g8", "seqid": contig,
              "start": start, "end": start + len(sequences[8]),
              "expected_start": start, "expected_end": start + len(sequences[8])}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [region])
    directory = rescue.rescue(output, plan, names[0], 1)
    models = [m for m in json.loads((directory / "models.json").read_text()) if m["status"] == "accepted"]
    assert len(models) == 1 and models[0]["seqid"] == contig and models[0]["sequence"] == sequences[8]
    assert (directory / "quality_flags.tsv").is_file()
    assert models[0]["quality_evidence"]["start_codon"] == sequences[8][:3]
    assert models[0]["quality_evidence"]["translation_initiation"] == "not_established"
    assert models[0]["quality_evidence"]["native_terminal_completeness"] == "not_established"
    exported = rescue.finalize(output, plan, [names[0]])
    gff = exported / "species_gff" / (names[0] + ".rescue.gff3")
    assert all(line.split("\t")[0] == contig for line in gff.read_text().splitlines() if not line.startswith("#"))
    assert not list((root / "genome").glob("*.fai"))
    originals = list(rescue.fasta_records(root / "cds" / (names[0] + ".cds.fa.gz")))
    augmented = list(rescue.fasta_records(exported / "species_cds" / (names[0] + ".rescue.cds.fa")))
    assert augmented[:len(originals)] == originals


@pytest.mark.parametrize("fallback", [False, True])
def test_mixed_species_genetic_codes_apply_to_local_and_genome_prediction(hidden_models, monkeypatch, fallback):
    root, names, sequences = hidden_models
    codes = root / "codes.tsv"
    codes.write_text("species\tgenetic_code\n" + names[0] + "\t4\n")
    output = root / "rescue"
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", output, "--genetic-codes", codes)
    plan = rescue.load(output)
    assert plan["request"]["sources"][names[0]]["genetic_code"] == 4
    assert plan["request"]["sources"][names[1]]["genetic_code"] == 1
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    start = 9 * (len(sequences[0]) + 60)
    region = {"id": "recoded_query", "donor": names[1], "query": names[1] + "_g9", "seqid": "chr1",
              "start": start, "end": start + len(sequences[9]),
              "expected_start": start, "expected_end": start + len(sequences[9])}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [region])
    if fallback:
        original_run = rescue.run
        def force_genome_search(command, directory, label, stdout=None):
            if label == "miniprot":
                Path(stdout).write_text("")
            else:
                original_run(command, directory, label, stdout)
        monkeypatch.setattr(rescue, "run", force_genome_search)
    directory = rescue.rescue(output, plan, names[0], 1)
    models = [m for m in json.loads((directory / "models.json").read_text()) if m["status"] == "accepted"]
    # TGA at position 80 is a stop under code 1, but tryptophan under code 4.
    expected = sequences[9][:240] + "TGA" + sequences[9][243:]
    assert len(models) == 1 and models[0]["sequence"] == expected
    assert models[0]["search"] == ("genome_fallback" if fallback else "synteny_interval")
    assert "*" not in str(Seq(expected).translate(table=4))[:-1]


def test_fewer_than_requested_eligible_references_fails(hidden_models):
    root, species, _ = hidden_models
    for name in species[-2:]:
        path = root / "busco" / (name + ".busco.short.txt")
        path.write_text(path.read_text().replace("C:95.0%", "C:85.0%"))
    with pytest.raises(AssertionError, match="Need 5 eligible common references; found 4"):
        make_plan(hidden_models)
    assert not (root / "rescue" / "plan.json").exists()


@pytest.mark.parametrize("tool", ["miniprot", "lastdb"])
def test_changed_binary_with_same_version_invalidates_frozen_plan(hidden_models, monkeypatch, tool):
    import os
    import shutil
    output, _, _ = make_plan(hidden_models)
    root = hidden_models[0]
    tools = root / "tools"
    tools.mkdir()
    executable = shutil.which(tool)
    wrapper = tools / tool
    wrapper.write_text(f"#!/bin/sh\nexec '{executable}' \"$@\"\n")
    wrapper.chmod(0o755)
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])
    with pytest.raises(ValueError, match="schema/tools changed"):
        rescue.load(output)


def test_tool_identity_binds_jcvi_support_code(tmp_path, monkeypatch):
    import importlib.util
    from types import SimpleNamespace
    package = tmp_path / "jcvi"
    (package / "algorithms").mkdir(parents=True)
    (package / "__init__.py").write_text("")
    support = package / "algorithms" / "lis.py"
    support.write_text("# Comparison support implementation A\n")
    original = importlib.util.find_spec
    def find_spec(name):
        if name == "jcvi":
            return SimpleNamespace(origin=str(package / "__init__.py"), submodule_search_locations=[str(package)])
        return original(name)
    monkeypatch.setattr(importlib.util, "find_spec", find_spec)
    before = rescue.identities()
    support.write_text("# Comparison support implementation B\n")
    assert rescue.identities() != before


def test_tool_identity_binds_canonical_busco_source_not_call_wrapper(tmp_path, monkeypatch):
    before = rescue.identities()
    implementation = rescue.busco_reference_implementation
    assert before["busco_quality_implementation"] == rescue.digest(implementation.__file__)
    original_quality = rescue.busco_quality
    monkeypatch.setattr(rescue, "busco_quality", lambda path: original_quality(path))
    assert rescue.identities() == before
    source = tmp_path / "busco_reference_quality.py"
    source.write_bytes(Path(implementation.__file__).read_bytes() + b"\n# changed source\n")
    monkeypatch.setattr(implementation, "__file__", str(source))
    assert rescue.identities()["busco_quality_implementation"] != before["busco_quality_implementation"]


@pytest.mark.parametrize("revision", ["date", "markers"])
def test_same_lineage_name_with_different_dataset_revision_is_not_comparable(hidden_models, revision):
    root, species, _ = hidden_models
    path = root / "busco" / (species[0] + ".busco.short.txt")
    if revision == "date":
        path.write_text(path.read_text().replace("embryophyta_odb12", "embryophyta_odb12 (Creation date: 2026-01-01)"))
    else:
        path.write_text(path.read_text().replace("n:100", "n:200"))
    with pytest.raises(AssertionError, match="dataset date and marker count"):
        make_plan(hidden_models)


def test_500_species_plan_has_balanced_references_and_sparse_comparisons(tmp_path):
    for directory in ("cds", "gff", "genome", "busco"):
        (tmp_path / directory).mkdir()
    species = [f"Plant_species{i:03d}" for i in range(500)]
    sequence = "ATG" + "GCT" * 79 + "TAA"
    for name in species:
        (tmp_path / "cds" / (name + ".cds.fa")).write_text(f">{name}_g0\n{sequence}\n")
        (tmp_path / "gff" / (name + ".gff3")).write_text("##gff-version 3\nchr1\tfixture\tgene\t1\t243\t.\t+\t.\tID=g0\n")
        (tmp_path / "genome" / (name + ".genome.fa")).write_text(f">chr1\n{sequence}\n")
        (tmp_path / "busco" / (name + ".busco.short.txt")).write_text(
            "# BUSCO version is: 6.0.0\n# The lineage dataset is: embryophyta_odb12\n"
            "# BUSCO was run in mode: transcriptome\nC:95.0%[S:95.0%,D:0.0%],F:0.0%,M:5.0%,n:100\n")
    clades = ["(" + ",".join(n + ":1" for n in species[i:i + 100]) + "):100" for i in range(0, 500, 100)]
    (tmp_path / "tree.nwk").write_text("(" + ",".join(clades) + ");\n")
    output, _, _ = make_plan((tmp_path, species, []))
    plan = rescue.load(output)
    assert {int(n.removeprefix("Plant_species")) // 100 for n in plan["common_references"]} == set(range(5))
    pairs = [(j["a"], j["b"]) for j in plan["synteny_jobs"] if j["kind"] == "pair"]
    assert len(pairs) <= 4000 and len(pairs) == len(set(pairs))
    assert len([j for j in plan["synteny_jobs"] if j["kind"] == "self"]) == 500
    assert all(len(plan["donors"][n]) <= 8 and n not in plan["donors"][n] for n in species)


def test_real_tools_recover_hidden_model_reject_disruptions_and_resume(hidden_models):
    root, species, _ = hidden_models
    gff = root / "gff" / (species[0] + ".gff3")
    embedded = "##FASTA\n" + (root / "genome" / (species[0] + ".genome.fa")).read_text()
    with gff.open("a") as handle:
        handle.write(embedded)
    original_annotation = gff.read_text().split("\n##FASTA\n")[0] + "\n"
    output, species, cds = make_plan(hidden_models)
    cli("run", "--output", output, "--cpus", 2)
    summary = json.loads((output / "augmented" / "summary.json").read_text())
    assert summary["rescued_models"] == 1
    assert summary["changed_species"] == [species[0]]
    assert summary["gene_loss_calls"] is False
    models = json.loads((output / "rescued" / species[0] / "models.json").read_text())
    accepted = [m for m in models if m["status"] == "accepted"]
    assert len(accepted) == 1 and accepted[0]["sequence"] == cds[8]
    assert accepted[0]["evidence"]["query"].endswith("_g8")
    rejected = [m for m in models if m["status"] == "unresolved"]
    assert any("internal_stop" in m["problems"] for m in rejected)
    assert any("frameshift" in m["problems"] for m in rejected)
    assert len(list(rescue.fasta_records(output / "augmented" / "species_cds" / (species[0] + ".rescue.cds.fa")))) == 16
    exported_gff = (output / "augmented" / "species_gff" / (species[0] + ".rescue.gff3")).read_text()
    assert exported_gff.startswith(original_annotation) and exported_gff.endswith(embedded)
    assert not list((hidden_models[0] / "genome").glob("*.fai"))
    benchmark = SCRIPT.parents[1] / "benchmarks/benchmark_rescue_search.py"
    for directory, extra in [("bounded_benchmark", []), ("complete_evidence_check", ["--check-existing"])]:
        result = subprocess.run([sys.executable, str(benchmark), "--evidence", "rescue", "--species", species[0],
                                 "--output", directory, "--cpus", "2", "--interval-count", "2",
                                 "--query-count", "2", "--repeats", "1", *extra], cwd=root,
                                capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        assert json.loads((root / directory / "result.json").read_text())["outputs_identical"] is True
    times = {str(p): p.stat().st_mtime_ns for p in output.rglob("receipt.json")}
    cli("run", "--output", output, "--cpus", 2)
    assert times == {str(p): p.stat().st_mtime_ns for p in output.rglob("receipt.json")}
    # Dependency corruption cannot hide behind an unchanged receipt digest.
    plan = rescue.load(output)
    dependency = next(j for j in plan["synteny_jobs"] if j["a"] == species[0] and j["kind"] == "pair")
    blocks = output / "synteny" / dependency["id"] / "blocks.json"
    saved = blocks.read_bytes()
    blocks.write_text("[]")
    with pytest.raises(ValueError, match="Comparison incomplete or corrupted"):
        rescue.finalize(output, plan)
    pending = json.loads(cli("status", "--output", output).stdout)
    assert dependency["index"] in pending["synteny_tasks"]
    blocks.write_bytes(saved)
    # Corruption cannot grant a completed stage or false downstream finalization.
    (output / "rescued" / species[0] / "audit.tsv").write_text("corrupted")
    with pytest.raises(ValueError, match="corrupted"):
        rescue.finalize(output, rescue.load(output))
    cli("rescue", "--output", output, "--task-index", 1, "--cpus", 2)
    cli("finalize", "--output", output)
    cli("qc", "--output", output, "--busco-dir", hidden_models[0] / "busco")
    assert (output / "qc_report" / "before_after.tsv").is_file()
    # Exercise the real core dispatch/QC namespace. BUSCO is stubbed here to
    # test scheduling and selective reruns without downloading a lineage DB.
    from workflow.tests.test_gg_input_generation_end_to_end import (
        _core_env,
        _install_fake_toolchain,
        _write_runtime_busco_dataset,
    )
    root = hidden_models[0]
    workspace = root / "workspace"
    workspace.mkdir()
    fake = _install_fake_toolchain(root / "core_tools")
    calls = root / "busco_calls"
    stub = "#!" + sys.executable + "\n" + '''import sys
from pathlib import Path
if "--version" in sys.argv:
    print("BUSCO 6.0.0")
    raise SystemExit(0)
def arg(flag):
    return sys.argv[sys.argv.index(flag)+1]
lineage = Path(arg("--lineage_dataset")).name
out = Path(arg("--out")) / ("run_" + lineage)
out.mkdir(parents=True, exist_ok=True)
(out / "full_table.tsv").write_text("# BUSCO version is: 6.0.0\\n" + "".join(
    f"BUSCO{i}\\tComplete\\tfixture{i}\\n" for i in range(1, 101)))
proteins = out / "busco_sequences/single_copy_busco_sequences"
proteins.mkdir(parents=True)
for i in range(1, 101):
    (proteins / f"BUSCO{i}.faa").write_text(f">fixture{i}\\nMKAAA\\n")
(out / "short_summary.txt").write_text("# BUSCO version is: 6.0.0\\n# The lineage dataset is: " + lineage + "\\n# BUSCO was run in mode: transcriptome\\nC:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:100\\n")
with Path(CALLS).open("a") as log:
    log.write("run\\n")
'''
    (fake / "busco").write_text(stub.replace("CALLS", repr(str(calls))))
    (fake / "busco").chmod(0o755)
    full = root / "full_busco"
    full.mkdir()
    for n in species:
        (full / (n + ".busco.full.tsv")).write_text("# initial fixture\n" + "".join(
            f"BUSCO{i}\tComplete\tfixture{i}\n" if i <= 95 else f"BUSCO{i}\tMissing\n" for i in range(1, 101)))
    _write_runtime_busco_dataset(workspace, "embryophyta_odb12")
    env = _core_env(workspace, None, fake, "rescue_models", task_id=1)
    env.update(gene_model_rescue_dir=str(output), species_busco_full_dir=str(full),
               run_gene_model_rescue_swissprot="0",  # This fixture verifies selective BUSCO scheduling, not annotation downloads.
               species_busco_short_dir=str(root / "busco"), species_cds_dir=str(root / "cds"),
               species_gff_dir=str(root / "gff"), species_genome_dir=str(root / "genome"), overwrite="0")
    core = SCRIPT.parent.parent / "core" / "gg_input_generation_core.sh"
    result = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert rescue.verified(output / "workers" / species[0], {"plan": rescue.digest(output / "plan.json"), "species": species[0]})
    # Host Slurm helpers can read these files without container-path aliases.
    worker_receipt = json.loads((output / "workers" / species[0] / "receipt.json").read_text())
    assert all(not Path(p).is_absolute() for p in worker_receipt["files"])
    env["input_generation_mode"] = "rescue_finalize"
    result = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert calls.read_text().splitlines() == ["run"]
    # Final aggregation must not repeat completed worker BUSCO even when the
    # surrounding entrypoint's overwrite flag is enabled.
    env["overwrite"] = "1"
    result = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert calls.read_text().splitlines() == ["run"]
    # A damaged effective input must fail before any QC can be silently skipped.
    effective = output / "effective" / species[0] / "inputs.tsv"
    original = effective.read_bytes()
    effective.unlink()
    result = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
    assert result.returncode != 0 and "Effective inputs incomplete or corrupted" in result.stderr
    effective.write_bytes(original)
    # A different workspace using the same custom rescue directory must still
    # respect that directory's active worker/finalizer phase lock.
    import os
    lock_script = SCRIPT.parent / "shared_namespace_lock.py"
    lock_path = output / ".array-phase.lock"
    result = subprocess.run([sys.executable, str(lock_script), "acquire-shared", str(lock_path),
                             "--owner-pid", str(os.getpid()), "--nonblocking"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    token = result.stdout.strip()
    try:
        result = subprocess.run(["bash", str(core)], env=env, capture_output=True, text=True)
        assert result.returncode != 0 and "already in use" in result.stderr
    finally:
        subprocess.run([sys.executable, str(lock_script), "release-shared", str(lock_path), "--token", token], check=True)
    qc = rescue.table(output / "qc_report" / "before_after.tsv")
    assert [r["species"] for r in qc if float(r["delta_pct"]) > 0] == [species[0]]
    assert rescue.busco_quality(root / "busco" / (species[0] + ".busco.short.txt"))["complete_pct"] == 95


def gene(identifier, start, seqid="chr1"):
    return {"gene_id": identifier, "seqid": seqid, "start": start, "end": start + 10, "strand": "+"}


def test_wgd_multiple_blocks_and_self_copy_candidates():
    target = {g["gene_id"]: g for g in [gene("a", 0), gene("b", 100), gene("c", 200), gene("d", 300)]}
    donor = {g["gene_id"]: g for g in [gene("x", 0), gene("missing", 50), gene("y", 100)]}
    params = {"max_interval": 1000, "padding": 0}
    found = rescue.candidate_intervals("Target_species", "Donor_species",
                                       [[["a", "x"], ["b", "y"]], [["c", "x"], ["d", "y"]]], target, donor, params)
    assert len(found) == 2 and {r["start"] for r in found} == {10, 210}
    all_genes = {**target, **{k: {**v, "seqid": "chr2"} for k, v in donor.items()}}
    self_found = rescue.candidate_intervals("Target_species", "Target_species",
                                            [[["a", "x"], ["b", "y"]]], all_genes, all_genes, params)
    assert any(r["query"] == "missing" for r in self_found)
    # A retained copy at a different locus never removes this candidate.
    donor["missing"]["start"] = 50
    assert all(r["query"] == "missing" for r in found)


@pytest.mark.parametrize("inverted", [False, True])
def test_parallel_copies_in_one_block_preserve_each_missing_locus(inverted):
    target = {g["gene_id"]: g for g in [gene("a", 0), gene("b", 100), gene("c", 200), gene("d", 300)]}
    donor = {g["gene_id"]: g for g in [gene("x", 0), gene("missing", 50), gene("y", 100)]}
    block = [["a", "x"], ["b", "y"], ["c", "y" if inverted else "x"], ["d", "x" if inverted else "y"]]
    params = {"max_interval": 1000, "padding": 0}
    # Anchor input order must not decide which WGD copy is examined.
    for anchors in (block, list(reversed(block))):
        found = rescue.candidate_intervals("Target_species", "Donor_species", [anchors], target, donor, params)
        expected = {(10, 100, "+"), (210, 300, "-" if inverted else "+")}
        assert expected <= {(r["start"], r["end"], r["orientation"]) for r in found if r["query"] == "missing"}


def test_many_to_many_copies_keep_local_flanks_in_both_directions():
    target = {g["gene_id"]: g for g in [gene("a", 0), gene("b", 100), gene("c", 200), gene("d", 300)]}
    donor = {g["gene_id"]: g for g in [gene("x", 0), gene("missing1", 50), gene("y", 100),
                                     gene("v", 200), gene("missing2", 250), gene("w", 300)]}
    block = [[a, b] for a, bs in (("a", ["x", "v"]), ("b", ["y", "w"]),
                                 ("c", ["x", "v"]), ("d", ["y", "w"])) for b in bs]
    params = {"max_interval": 1000, "padding": 0}
    found = rescue.candidate_intervals("Target_species", "Donor_species", [block], target, donor, params)
    local = [r for r in found if (r["start"], r["end"]) in {(10, 100), (210, 300)}]
    assert {(r["query"], r["start"]) for r in local} == {(q, s) for q in ("missing1", "missing2") for s in (10, 210)}
    # Neither interval is supported by crossing flanks from different donor copies.
    assert all((r["left_anchor"], r["right_anchor"]) in {("a", "b"), ("c", "d")} for r in local)
    assert all((r["donor_left_anchor"], r["donor_right_anchor"]) in {("x", "y"), ("v", "w")} for r in local)


def test_parallel_donor_path_matching_retains_ambiguous_ties():
    left = [(1, "a"), (3, "b")]
    right = [(0, "x"), (2, "y"), (4, "z")]
    assert set(rescue.local_flanks(left, right)) == {
        (left[0], right[0]), (left[0], right[1]), (left[1], right[1]), (left[1], right[2])}
    assert set(rescue.local_flanks(right, left)) == {(b, a) for a, b in rescue.local_flanks(left, right)}


def test_local_flanks_matches_exhaustive_ordered_assignment_oracle():
    rng = random.Random(631)
    for _ in range(100):
        left = [(n, "l" + str(i)) for i, n in enumerate(sorted(rng.sample(range(20), rng.randint(1, 5))))]
        right = [(n, "r" + str(i)) for i, n in enumerate(sorted(rng.sample(range(20), rng.randint(1, 5))))]
        small, large = (left, right) if len(left) <= len(right) else (right, left)
        assignments = [list(zip(small, choice, strict=True)) for choice in itertools.combinations(large, len(small))]
        minimum = min(sum(abs(a[0] - b[0]) for a, b in pairs) for pairs in assignments)
        expected = {pair for pairs in assignments if sum(abs(a[0] - b[0]) for a, b in pairs) == minimum for pair in pairs}
        if len(left) > len(right):
            expected = {(b, a) for a, b in expected}
        assert set(rescue.local_flanks(left, right)) == expected


@pytest.mark.parametrize("copy_contig", ["chr1", "chr2"])
def test_real_unquota_self_synteny_finds_missing_copy_despite_retained_copy(hidden_models, copy_contig):
    root, species, cds = hidden_models
    name = species[0]
    original = root / "genome" / (name + ".genome.fa")
    offset = len(original.read_text().splitlines()[1]) if copy_contig == "chr1" else 0
    genome = ""
    features, fasta = [], []
    for i, sequence in enumerate(cds):
        start = offset + len(genome)
        genome += sequence + "N" * 60
        features += [f"{copy_contig}\tsynthetic\tgene\t{start + 1}\t{start + len(sequence)}\t.\t+\t.\tID=b{i}",
                     f"{copy_contig}\tsynthetic\tmRNA\t{start + 1}\t{start + len(sequence)}\t.\t+\t.\tID=b{i}.t1;Parent=b{i}",
                     f"{copy_contig}\tsynthetic\tCDS\t{start + 1}\t{start + len(sequence)}\t.\t+\t0\tParent=b{i}.t1"]
        fasta.append(f">{name}_b{i}\n{sequence}\n")
    if copy_contig == "chr1":
        original.write_text(original.read_text().rstrip() + genome + "\n")
    else:
        with original.open("a") as handle:
            handle.write(">chr2\n" + genome + "\n")
    with (root / "gff" / (name + ".gff3")).open("a") as handle:
        handle.write("\n".join(features) + "\n")
    with (root / "cds" / (name + ".cds.fa")).open("a") as handle:
        handle.write("".join(fasta))
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    job = next(j for j in plan["synteny_jobs"] if j["kind"] == "self" and j["a"] == name)
    directory = rescue.synteny(output, plan, job["index"], 2)
    positions = {g["gene_id"]: g for g in json.loads((output / "prepared" / name / "positions.json").read_text())}
    blocks = json.loads((directory / "blocks.json").read_text())
    found = rescue.candidate_intervals(name, name, blocks, positions, positions, plan["request"]["parameters"])
    found += rescue.candidate_intervals(name, name, [[[b, a] for a, b in block] for block in blocks],
                                        positions, positions, plan["request"]["parameters"])
    if copy_contig == "chr2":
        assert any(r["query"] == name + "_b8" and r["seqid"] == "chr1" for r in found)
    pair = next(j for j in plan["synteny_jobs"] if j["a"] == name and j["b"] == species[1])
    directory = rescue.synteny(output, plan, pair["index"], 2)
    donor_positions = {g["gene_id"]: g for g in json.loads((output / "prepared" / species[1] / "positions.json").read_text())}
    blocks = json.loads((directory / "blocks.json").read_text())
    anchors = [a for block in blocks for a, _ in block]
    assert any(a.startswith(name + "_g") for a in anchors)
    assert any(a.startswith(name + "_b") for a in anchors)
    found = rescue.candidate_intervals(name, species[1], blocks, positions, donor_positions, plan["request"]["parameters"])
    assert any(r["query"] == species[1] + "_g8" and r["seqid"] == "chr1" for r in found)


def test_padding_does_not_accept_a_model_outside_its_flanking_anchors(hidden_models, monkeypatch):
    output, names, sequences = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    start = 8 * (len(sequences[0]) + 60)
    region = {"id": "padded_window", "donor": names[1], "query": names[1] + "_g8", "seqid": "chr1",
              "start": start, "end": start + len(sequences[8]) + 60,
              "expected_start": start + len(sequences[8]), "expected_end": start + len(sequences[8]) + 60}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [region])
    directory = rescue.rescue(output, plan, names[0], 1)
    models = json.loads((directory / "models.json").read_text())
    assert any(m["sequence"] == sequences[8] for m in models)
    assert all(m["status"] != "accepted" for m in models)
    assert any("outside_expected_synteny_interval" in m["problems"] for m in models)


@pytest.mark.parametrize("query,fallback", [(8, False), (11, False), (8, True)])
def test_only_conflict_free_models_skip_refinement(hidden_models, monkeypatch, query, fallback):
    output, names, sequences = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    plan["request"]["parameters"]["genome_fallback"] = int(fallback)
    plan["request"]["gemoma_jar"] = "test-refinement.jar"
    rescue.atomic_json(output / "plan.json", plan)
    start = query * (len(sequences[0]) + 60) + (1 if query >= 11 else 0)
    region = {"id": "query", "donor": names[1], "query": names[1] + f"_g{query}", "seqid": "chr1",
              "start": start, "end": start + len(sequences[query]),
              "expected_start": start, "expected_end": start + len(sequences[query])}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [region])
    if fallback:
        original_run = rescue.run
        def fail_local_alignment(command, directory, label, stdout=None):
            if label == "miniprot":
                Path(stdout).write_text("")
            else:
                original_run(command, directory, label, stdout)
        monkeypatch.setattr(rescue, "run", fail_local_alignment)
    refined = []
    monkeypatch.setattr(rescue, "refine_gemoma", lambda tmp, root, plan, source, regions, *args: refined.extend(regions))
    directory = rescue.rescue(output, plan, names[0], 1)
    models = json.loads((directory / "models.json").read_text())
    assert bool(refined) == (query == 11)
    assert any(m["sequence"] == sequences[query] for m in models)
    assert bool([m for m in models if m["status"] == "accepted"]) == (query == 8)
    if fallback:
        assert all(m["search"] == "genome_fallback" for m in models)


def test_rescue_refuses_dependencies_replaced_during_prediction(hidden_models, monkeypatch):
    output, names, sequences = make_plan(hidden_models)
    plan = rescue.load(output)
    plan["donors"][names[0]] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    start = 8 * (len(sequences[0]) + 60)
    region = {"id": "query", "donor": names[1], "query": names[1] + "_g8", "seqid": "chr1",
              "start": start, "end": start + len(sequences[8]),
              "expected_start": start, "expected_end": start + len(sequences[8])}
    monkeypatch.setattr(rescue, "candidates", lambda *_: [region])
    original_run = rescue.run
    def replace_dependency(command, directory, label, stdout=None):
        original_run(command, directory, label, stdout)
        if label == "miniprot":
            directory = output / "prepared" / names[1]
            protein = directory / "genes.pep"
            protein.write_text(protein.read_text() + "\n")
            receipt = json.loads((directory / "receipt.json").read_text())
            receipt["files"]["genes.pep"] = rescue.digest(protein)
            rescue.atomic_json(directory / "receipt.json", receipt)
    monkeypatch.setattr(rescue, "run", replace_dependency)
    with pytest.raises(ValueError, match="dependencies changed"):
        rescue.rescue(output, plan, names[0], 1)
    assert not (output / "rescued" / names[0] / "receipt.json").exists()
    assert (output / "rescued" / (names[0] + ".failed") / "models.json").exists()


def test_negative_strand_phase_and_stop_are_reconstructed(tmp_path):
    import pysam
    coding = "ATGAAAGCCTAA"
    genome_path = tmp_path / "genome.fa"
    genome_path.write_text(">chr1\n" + str(Seq(coding).reverse_complement()) + "\n")
    pysam.faidx(str(genome_path))
    model = {"seqid": "chr1", "strand": "-", "cds": [[3, 12, 0]], "frameshift": False,
             "coverage": 1, "identity": 1}
    with pysam.FastaFile(str(genome_path)) as genome:
        checked = rescue.validate_model(model, genome, 1, {"minimum_coverage": 0.95, "minimum_identity": 0.5})
    assert checked["sequence"] == coding and checked["problems"] == []
    assert checked["cds"] == [[0, 12, 0]]
    assert model["cds"] == [[3, 12, 0]]


def test_split_codon_splice_phase_and_assembly_gap_qc(tmp_path):
    import pysam
    coding = "ATGAAAGCCTAA"
    path = tmp_path / "genome.fa"
    path.write_text(">good\n" + coding[:4] + "GTACACAG" + coding[4:] + "\n>gap\n" +
                    coding[:4] + "GTNNNNAG" + coding[4:] + "\n")
    pysam.faidx(str(path))
    base = {"strand": "+", "frameshift": False, "coverage": 1, "identity": 1}
    params = {"minimum_coverage": 0.95, "minimum_identity": 0.5}
    with pysam.FastaFile(str(path)) as genome:
        good = rescue.validate_model({**base, "seqid": "good", "cds": [[0, 4, 0], [12, 20, 2]]}, genome, 1, params)
        bad = rescue.validate_model({**base, "seqid": "good", "cds": [[0, 4, 0], [12, 20, 0]]}, genome, 1, params)
        gap = rescue.validate_model({**base, "seqid": "gap", "cds": [[0, 4, 0], [12, 20, 2]]}, genome, 1, params)
    assert good["sequence"] == coding and good["problems"] == []
    assert "invalid_phase" in bad["problems"]
    assert "assembly_gap_within_model_span" in gap["problems"]


@pytest.mark.parametrize("donor,acceptor,canonical", [("GT", "AG", True), ("GC", "AG", True),
                                                      ("AT", "AC", True), ("AT", "AG", False)])
def test_search_sensitivity_does_not_accept_noncanonical_at_ag_splice(tmp_path, donor, acceptor, canonical):
    import pysam
    coding = "ATGAAAGCCTAA"
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\n" + coding[:4] + donor + "CCCC" + acceptor + coding[4:] + "\n")
    pysam.faidx(str(path))
    model = {"seqid": "chr1", "strand": "+", "cds": [[0, 4, 0], [12, 20, 2]],
             "frameshift": False, "coverage": 1, "identity": 1}
    with pysam.FastaFile(str(path)) as genome:
        checked = rescue.validate_model(model, genome, 1, {"minimum_coverage": .95, "minimum_identity": .5})
    assert checked["sequence"] == coding
    assert checked["problems"] == ([] if canonical else ["noncanonical_splice"])


def test_competing_models_and_existing_annotations_are_not_added():
    base = {"seqid": "chr1", "strand": "+", "problems": [], "query": "q", "evidence": {}, "cds": [[20, 80, 0]]}
    a = {**base, "problems": []}
    b = {**base, "problems": [], "cds": [[20, 83, 0]]}
    result = rescue.consolidate([a, b], [], "Target_species")
    assert all(m["status"] == "unresolved" and "competing_new_models" in m["problems"] for m in result)
    result = rescue.consolidate([{**base, "problems": []}], [gene("annotated", 70)], "Target_species")
    assert result[0]["status"] == "unresolved"


def test_failed_stage_preserves_previous_output_and_diagnostics(tmp_path):
    root = tmp_path / "rescue"
    root.mkdir()
    rescue.stage(root, Path("job"), "first", lambda d: (d / "data").write_text("complete"))
    def fail(directory):
        (directory / "log").write_text("failure evidence")
        raise RuntimeError("deliberate failure")
    with pytest.raises(RuntimeError, match="deliberate"):
        rescue.stage(root, Path("job"), "second", fail)
    assert (root / "job" / "data").read_text() == "complete"
    assert (root / "job.failed" / "log").read_text() == "failure evidence"


@pytest.mark.parametrize("receipt", [None, [], "broken", {}, {"key": "k", "files": []},
                                   {"key": "k", "files": {"data": None}}])
def test_malformed_receipt_is_pending_without_crashing(tmp_path, receipt):
    (tmp_path / "receipt.json").write_text(json.dumps(receipt))
    assert rescue.verified(tmp_path, "k") is False


def test_failed_publication_restores_previous_complete_directory(tmp_path, monkeypatch):
    rescue.stage(tmp_path, Path("job"), "old", lambda d: (d / "data").write_text("old"))
    rename = Path.rename
    def fail_replacement(path, target):
        if path.name.startswith(".working-") and Path(target).name == "job":
            raise OSError("publication failure")
        return rename(path, target)
    monkeypatch.setattr(Path, "rename", fail_replacement)
    with pytest.raises(OSError, match="publication failure"):
        rescue.stage(tmp_path, Path("job"), "new", lambda d: (d / "data").write_text("new"))
    assert rescue.verified(tmp_path / "job", "old")
    assert (tmp_path / "job.failed" / "data").read_text() == "new"


@pytest.mark.parametrize("moment", ["journal_written", "old_renamed", "new_published"])
def test_sigkill_during_publication_recovers_without_losing_previous_result(tmp_path, moment):
    rescue.stage(tmp_path, Path("job"), "old", lambda d: (d / "data").write_text("old"))
    code = '''import os, signal, sys
from pathlib import Path
from workflow.support import rescue_gene_models as r
root = Path(sys.argv[1])
moment = sys.argv[2]
rename = Path.rename
def kill_at_publication(path, target):
    if moment == "journal_written" and path == root / "job":
        os.kill(os.getpid(), signal.SIGKILL)
    result = rename(path, target)
    if ((moment == "old_renamed" and path == root / "job") or
            (moment == "new_published" and path.name.startswith(".working-") and Path(target) == root / "job")):
        os.kill(os.getpid(), signal.SIGKILL)
    return result
Path.rename = kill_at_publication
r.stage(root, Path("job"), "new", lambda d: (d / "data").write_text("new"))
'''
    result = subprocess.run([sys.executable, "-c", code, str(tmp_path), moment], capture_output=True, text=True)
    assert result.returncode == -9, result.stderr
    if moment != "new_published":
        def fail(directory):
            raise RuntimeError("retry cannot rebuild")
        with pytest.raises(RuntimeError, match="cannot rebuild"):
            rescue.stage(tmp_path, Path("job"), "new", fail)
        assert rescue.verified(tmp_path / "job", "old")
    else:
        rescue.stage(tmp_path, Path("job"), "new", lambda d: pytest.fail("must reuse published result"))
        assert rescue.verified(tmp_path / "job", "new")
        assert not list(tmp_path.glob(".previous-*"))


@pytest.mark.parametrize("corruption", ["payload", "other_job_receipt"])
def test_publication_recovery_restores_backup_when_new_output_is_corrupted(tmp_path, corruption):
    rescue.stage(tmp_path, Path("job"), "old", lambda d: (d / "data").write_text("old"))
    token = hashlib.sha256(b"job").hexdigest()[:12]
    backup, temporary = ".previous-" + token + "-test", ".working-" + token + "-test"
    (tmp_path / "job").rename(tmp_path / backup)
    rescue.stage(tmp_path, Path("job"), "new", lambda d: (d / "data").write_text("new"))
    if corruption == "payload":
        (tmp_path / "job" / "data").write_text("corrupted")
    else:
        receipt = json.loads((tmp_path / "job" / "receipt.json").read_text())
        receipt["key"] = "another-job"
        (tmp_path / "job" / "receipt.json").write_text(json.dumps(receipt))
    journal = tmp_path / ".locks" / "job.publish.json"
    journal.write_text(json.dumps({"destination": "job", "backup": backup, "temporary": temporary, "key": "new"}))
    rescue.stage(tmp_path, Path("job"), "old", lambda d: pytest.fail("must recover old output"))
    assert rescue.verified(tmp_path / "job", "old")
    quarantined = list(tmp_path.glob(".interrupted-*"))
    assert len(quarantined) == 1
    assert (quarantined[0] / "data").read_text() == ("corrupted" if corruption == "payload" else "new")


@pytest.mark.parametrize("bad", ["../unrelated", ".previous-other-job-test", None])
def test_unsafe_publication_journal_leaves_previous_and_foreign_data_intact(tmp_path, bad):
    rescue.stage(tmp_path, Path("job"), "old", lambda d: (d / "data").write_text("old"))
    token = hashlib.sha256(b"job").hexdigest()[:12]
    journal = tmp_path / ".locks" / "job.publish.json"
    journal.write_text(json.dumps({"destination": "job", "backup": bad, "temporary": ".working-" + token + "-test", "key": "new"}))
    foreign = tmp_path / "unrelated"
    foreign.mkdir()
    (foreign / "data").write_text("keep")
    with pytest.raises(ValueError, match="Unsafe publication journal"):
        rescue.stage(tmp_path, Path("job"), "new", lambda d: pytest.fail("must not build"))
    assert rescue.verified(tmp_path / "job", "old")
    assert (foreign / "data").read_text() == "keep"


def test_concurrent_workers_publish_the_same_stage_only_once(tmp_path):
    code = '''import os, sys, time
from pathlib import Path
from workflow.support import rescue_gene_models as r
root = Path(sys.argv[1])
def build(directory):
    (root / ("built-" + str(os.getpid()))).write_text("called")
    time.sleep(.15)
    (directory / "data").write_text("complete")
r.stage(root, Path("job"), "same", build)
'''
    processes = [subprocess.Popen([sys.executable, "-c", code, str(tmp_path)], stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, text=True) for _ in range(4)]
    for process in processes:
        stdout, stderr = process.communicate(timeout=30)
        assert process.returncode == 0, stdout + stderr
    assert len(list(tmp_path.glob("built-*"))) == 1
    assert rescue.verified(tmp_path / "job", "same")
    assert not list(tmp_path.glob(".working-*"))


@pytest.mark.parametrize("changed", ["busco", "inputs"])
def test_external_qc_refuses_inputs_changed_during_report(hidden_models, monkeypatch, changed):
    root, _, _ = hidden_models
    output, names, _ = make_plan(hidden_models)
    augmented = output / "augmented"
    augmented.mkdir()
    (augmented / "receipt.json").write_text("{}\n")
    rescue.write_tsv(augmented / "inputs.tsv", ("species", "rescued_models"), [(n, 0) for n in names])
    monkeypatch.setattr(rescue, "finalize", lambda *_, **__: augmented)
    post = root / "post_busco"
    post.mkdir()
    for source in (root / "busco").iterdir():
        (post / source.name).write_bytes(source.read_bytes())
    original_quality = rescue.busco_quality
    def replace_after_read(path):
        result = original_quality(path)
        if Path(path).parent == post:
            target = Path(path) if changed == "busco" else augmented / "inputs.tsv"
            with target.open("a") as handle:
                handle.write("# replaced while reporting\n")
        return result
    monkeypatch.setattr(rescue, "busco_quality", replace_after_read)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "qc", "--output", str(output), "--busco-dir", str(post)])
    with pytest.raises(ValueError, match="changed during execution"):
        rescue.main()
    assert not (output / "qc_report" / "receipt.json").exists()


def test_worker_completion_refuses_qc_changed_after_comparability_check(hidden_models, monkeypatch):
    root, names, _ = hidden_models
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    name = names[0]
    plan["donors"][name] = [names[1]]
    plan["synteny_jobs"] = []
    monkeypatch.setattr(rescue, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
    rescue.atomic_json(output / "plan.json", plan)
    monkeypatch.setattr(rescue, "candidates", lambda *_: [])
    rescue.rescue(output, plan, name, 1)
    rescue.finalize(output, plan, [name], Path("effective") / name)
    short_dir = output / "qc/species_cds_busco_short"
    full_dir = output / "qc/species_cds_busco_full"
    short_dir.mkdir(parents=True)
    full_dir.mkdir()
    summary = short_dir / (name + ".busco.short.txt")
    summary.write_bytes((root / "busco" / summary.name).read_bytes())
    (full_dir / (name + ".busco.full.tsv")).write_text("# synthetic complete table\n")
    original_quality = rescue.busco_quality
    def replace_after_check(path):
        result = original_quality(path)
        if Path(path) == summary:
            summary.write_text(summary.read_text().replace("n:100", "n:101"))
        return result
    monkeypatch.setattr(rescue, "busco_quality", replace_after_check)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "worker-complete", "--output", str(output), "--task-index", "1"])
    with pytest.raises(ValueError, match="changed during execution"):
        rescue.main()
    assert not (output / "workers" / name / "receipt.json").exists()


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("phase", [0, 1, 2])
def test_real_miniprot_spliced_cds_roundtrip(tmp_path, strand, phase):
    import pysam
    rng = random.Random(919)
    codons = ["GCT", "CGT", "AAC", "GAC", "TGC", "CAA", "GAA", "GGT", "CAC", "ATT",
              "CTG", "AAA", "ATG", "TTC", "CCT", "TCT", "ACT", "TGG", "TAC", "GTT"]
    coding = "ATG" + "".join(rng.choice(codons) for _ in range(159)) + "TAA"
    cut = 240 + phase
    intron = "GTAAGT" + "".join(rng.choice("ACGT") for _ in range(300)) + "TCTTTCTTTCAG"
    sequence = coding[:cut] + intron + coding[cut:]
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\n" + sequence + "\n")
    query = tmp_path / "query.fa"
    query.write_text(">donor\n" + str(Seq(coding).translate()).rstrip("*") + "\n")
    rescue.run(["miniprot", "--gff", path, query], tmp_path, "real_miniprot", tmp_path / "models.gff")
    models = rescue.read_miniprot(tmp_path / "models.gff")
    assert len(models) == 1 and len(models[0]["cds"]) == 2
    assert models[0]["strand"] == strand and models[0]["coverage"] == 1
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as genome:
        checked = rescue.validate_model(models[0], genome, 1, {"minimum_coverage": .95, "minimum_identity": .5})
    assert checked["sequence"] == coding and checked["problems"] == []


def test_real_miniprot_internal_deletion_is_not_full_coverage(hidden_models):
    import pysam
    root, _, cds = hidden_models
    coding = cds[8]
    path = root / "deletion.fa"
    path.write_text(">chr1\n" + coding[:240] + coding[300:] + "\n")
    query = root / "query.fa"
    query.write_text(">donor\n" + str(Seq(coding).translate()).rstrip("*") + "\n")
    rescue.run(["miniprot", "--gff", path, query], root, "deletion", root / "deletion.gff")
    models = rescue.read_miniprot(root / "deletion.gff")
    assert len(models) == 1
    assert models[0]["query_span_coverage"] == 1
    assert models[0]["coverage"] == .875
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as genome:
        checked = rescue.validate_model(models[0], genome, 1, {"minimum_coverage": .95, "minimum_identity": .5})
    assert "low_coverage" in checked["problems"]


def test_real_pair_without_anchors_is_complete_empty_evidence(hidden_models):
    root, species, _ = hidden_models
    # Valid, mapped sequences with no homologous proteins to the other species.
    path = root / "cds" / (species[1] + ".cds.fa")
    records = list(rescue.fasta_records(path))
    path.write_text("".join(f">{i}\nATG{'GCT' * 159}TAA\n" for i, _, _ in records))
    output, _, _ = make_plan(hidden_models)
    plan = rescue.load(output)
    job = next(j for j in plan["synteny_jobs"] if j["a"] == species[0] and j["b"] == species[1])
    directory = rescue.synteny(output, plan, job["index"], 2)
    assert json.loads((directory / "blocks.json").read_text()) == []
    assert rescue.verified(directory, rescue.comparison_key(output, job))


@pytest.mark.parametrize("feature,attribute", [("mRNA", "ID"), ("CDS", "Parent")])
def test_real_spliced_rescue_exports_transcript_level_mapping(hidden_models, feature, attribute):
    import pysam
    root, species, cds = hidden_models
    name = species[0]
    rng = random.Random(321)
    intron = "GTAAGT" + "".join(rng.choice("ACGT") for _ in range(300)) + "TCTTTCTTTCAG"
    genome, features = "", ["##gff-version 3"]
    for i, coding in enumerate(cds):
        start = len(genome)
        if i == 8:
            encoded = str(Seq(coding[:241] + intron + coding[241:]).reverse_complement())
        elif i == 9:
            encoded = coding[:240] + "TGA" + coding[243:]
        elif i == 10:
            encoded = coding[:240] + "A" + coding[240:]
        else:
            encoded = coding
        genome += encoded + "N" * 60
        if i in {8, 9, 10}:
            continue
        end = start + len(coding)
        features += [f"chr1\tsynthetic\tgene\t{start + 1}\t{end}\t.\t+\t.\tID=g{i}",
                     f"chr1\tsynthetic\tmRNA\t{start + 1}\t{end}\t.\t+\t.\tID=g{i}.t1;Parent=g{i}",
                     f"chr1\tsynthetic\tCDS\t{start + 1}\t{end}\t.\t+\t0\tParent=g{i}.t1"]
    (root / "genome" / (name + ".genome.fa")).write_text(">chr1\n" + genome + "\n")
    (root / "gff" / (name + ".gff3")).write_text("\n".join(features) + "\n")
    for n in species:
        path = root / "cds" / (n + ".cds.fa")
        records = list(rescue.fasta_records(path))
        path.write_text("".join(f">{i}.t1\n{s}\n" for i, _, s in records))
    output = root / "rescue"
    cli("plan", "--cds-dir", root / "cds", "--gff-dir", root / "gff", "--genome-dir", root / "genome",
        "--busco-dir", root / "busco", "--tree", root / "tree.nwk", "--output", output,
        "--common-references", 1, "--nearest-references", 1, "--feature", feature, "--attribute", attribute)
    cli("run", "--output", output, "--cpus", 2)
    models = [m for m in json.loads((output / "rescued" / name / "models.json").read_text()) if m["status"] == "accepted"]
    assert len(models) == 1 and models[0]["sequence"] == cds[8] and models[0]["strand"] == "-"
    assert len(models[0]["cds"]) == 2
    rows = [line.split("\t") for line in (output / "augmented" / "species_gff" / (name + ".rescue.gff3")).read_text().splitlines()
            if "genegalleon_rescue" in line]
    assert len(rows) == 4
    exons = [f for f in rows if f[2] == "CDS"]
    pysam.faidx(str(root / "genome" / (name + ".genome.fa")))
    with pysam.FastaFile(str(root / "genome" / (name + ".genome.fa"))) as assembly:
        reconstructed = "".join(str(Seq(assembly.fetch(f[0], int(f[3]) - 1, int(f[4]))).reverse_complement())
                                for f in sorted(exons, key=lambda f: int(f[3]), reverse=True))
    assert reconstructed == cds[8]
    assert {rescue.attributes(f[8])["Parent"] for f in exons} == {models[0]["model_id"]}


@pytest.mark.parametrize("command", ["synteny", "rescue"])
def test_zero_task_index_is_rejected_instead_of_running_all(hidden_models, command):
    output, _, _ = make_plan(hidden_models)
    result = subprocess.run([sys.executable, str(SCRIPT), command, "--output", str(output), "--task-index", "0"],
                            capture_output=True, text=True)
    assert result.returncode != 0 and "task index outside frozen plan" in result.stderr
    assert not (output / "synteny").exists()


def test_interval_index_and_competing_sweep_match_brute_force():
    import copy
    rng = random.Random(782)
    existing = [{"seqid": rng.choice(["chr1", "chr2"]), "start": rng.randrange(0, 1000), "end": 0} for _ in range(100)]
    for old in existing:
        old["end"] = old["start"] + rng.randrange(1, 30)
    base = [{"seqid": rng.choice(["chr1", "chr2"]), "strand": "+", "problems": [],
             "query": str(i), "evidence": {}, "cds": [[s, s + rng.randrange(1, 40), 0]]}
            for i in range(100) for s in [rng.randrange(0, 1500)]]
    original = copy.deepcopy(base)
    actual = rescue.consolidate(base, existing, "Target_species")
    # Independent, slow oracle verifies boundaries, nested spans and chains.
    for model in original:
        model["start"], model["end"] = model["cds"][0][:2]
        if any(rescue.overlaps(model, old) for old in existing):
            model["problems"].append("overlap_existing_annotation")
    eligible = [m for m in original if not m["problems"]]
    competing = {m["query"] for m in eligible for n in eligible
                 if m["cds"] != n["cds"] and rescue.overlaps(m, n)}
    for model, expected in zip(actual, original, strict=True):
        if expected["query"] in competing:
            expected["problems"].append("competing_new_models")
        assert set(model["problems"]) == set(expected["problems"])
    identical = {"seqid": "chr1", "strand": "+", "problems": [], "query": "a", "evidence": {}, "cds": [[0, 60, 0]]}
    records = [copy.deepcopy(identical), copy.deepcopy(identical), {**identical, "cds": [[20, 80, 0]], "problems": []}]
    assert all(m["status"] == "unresolved" for m in rescue.consolidate(records, [], "Target_species"))


def test_rescue_preflight_rejects_legacy_invalid_gff_without_mutation(hidden_models):
    root, species, _ = hidden_models
    gff = root / "gff" / (species[0] + ".gff3")
    gff.write_text(gff.read_text().replace("\tsynthetic\tgene\t", "\tfunannotate\tgene\t", 1)
                   .replace("ID=g0\n", "ID=g0;Name=SULTR4;1;\n", 1))
    before = gff.read_bytes()
    result = subprocess.run([sys.executable, str(SCRIPT), "plan", "--cds-dir", str(root / "cds"),
                             "--gff-dir", str(root / "gff"), "--genome-dir", str(root / "genome"),
                             "--busco-dir", str(root / "busco"), "--tree", str(root / "tree.nwk"),
                             "--output", str(root / "rescue")], capture_output=True, text=True)
    assert result.returncode != 0
    assert str(gff) + ":3:" in result.stderr and "regenerate with gg_input_generation" in result.stderr
    assert not (root / "rescue" / "plan.json").exists()
    assert gff.read_bytes() == before


def test_formatted_funannotate_metadata_reaches_real_anchor_reader(tmp_path):
    from kffractbias.io import parse_attributes

    from workflow.support.gff_attribute_syntax import normalise_line, validate_gff

    source, output = anchor_fixture(tmp_path, [{}, {}])
    gff = Path(source["gff"])
    text = gff.read_text().replace("\tsynthetic\t", "\tfunannotate\t")
    text = text.replace("ID=g0\n", "ID=g0;Name=SULTR4;1_1;\n")
    text = text.replace("ID=g0.t1;Parent=g0\n", "ID=g0.t1;Parent=g0;product=Protein 1;3, variant 2;\n")
    changes = []
    formatted = tmp_path / "formatted.gff3"
    formatted.write_text("".join(normalise_line(line, gff, number, changes)
                                for number, line in enumerate(text.splitlines(keepends=True), 1)))
    assert len(changes) == 2
    validate_gff(formatted)
    attributes = [parse_attributes(line.split("\t")[8]) for line in formatted.read_text().splitlines()
                  if not line.startswith("#")]
    assert attributes[0]["Name"] == ("SULTR4;1_1",)
    assert attributes[1]["product"] == ("Protein 1;3, variant 2",)
    assert parse_attributes("Parent=t1,t2;Note=literal%253B;") == {
        "Parent": ("t1", "t2"), "Note": ("literal%3B",)}
    source["gff"] = str(formatted)
    genes, _ = prepare_rescue_genome(source, output, "genes", 1.0)
    assert [gene.gene_id for gene in genes] == ["Plant_example_g0", "Plant_example_g1"]


def test_gemoma_reads_proteins_once_per_donor_and_keeps_each_query_alignment(tmp_path, monkeypatch):
    donor = "Donor_species"
    prepared = tmp_path / "prepared" / donor
    prepared.mkdir(parents=True)
    identifiers = [donor + "_g1", donor + "_g2"]
    rescue.write_tsv(prepared / "genes.id_map.tsv", ("original_id", "jcvi_id", "locus_id", "status"),
                     [(identifier, identifier, identifier, "selected") for identifier in identifiers])
    (prepared / "genes.anchor_admission.json").write_text('{"records": []}')
    (prepared / "genes.pep").write_text(f">{identifiers[0]}\nMK\n>{identifiers[1]}\nMKP\n")
    gff, jar, java = (tmp_path / name for name in ("donor.gff", "test.jar", "java"))
    gff.write_text("".join(f"chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID={identifier};Parent={identifier}\n" for identifier in identifiers))
    jar.write_text("jar")
    java.write_text("java")
    plan = {"request": {"gemoma_jar": str(jar), "gemoma_java": str(java),
                        "files": {str(path): rescue.digest(path) for path in (jar, java)},
                        "parameters": {"minimum_coverage": 0.9, "minimum_identity": 0.9},
                        "sources": {donor: {"gff": str(gff), "genome": str(tmp_path / "donor.fa")}}}}
    regions = [{"donor": donor, "query": query, "id": f"query_{index}", "seqid": "chr1",
                "expected_strand": "+", "start": 0, "end": 9} for index, query in enumerate(identifiers)]
    monkeypatch.setattr(rescue, "verify_sources", lambda *_: None)
    def fake_java(command, directory, query):
        out = Path(next(value.split("=", 1)[1] for value in command if value.startswith("outdir=")))
        out.mkdir()
        (out / "final_annotation.gff").write_text("".join(
            f"chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=model{index}\n"
            f"chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=model{index}\n" for index in (1, 2, 3)))
    monkeypatch.setattr(rescue, "run", fake_java)
    reads = []
    original = rescue.fasta_records
    def counted(path):
        reads.append(path)
        return original(path)
    monkeypatch.setattr(rescue, "fasta_records", counted)
    def checked(model, *args):
        return {**model, "sequence": "ATGAAATAA" if model["query"] == "query_0" else "ATGAAACCCTAA", "problems": []}
    monkeypatch.setattr(rescue, "validate_model", checked)
    monkeypatch.setattr(rescue, "check_interval", lambda model: model)
    validated = []
    rescue.refine_gemoma(tmp_path, tmp_path, plan, {"genetic_code": 1}, regions, None, validated, 1)
    assert reads == [prepared / "genes.pep"]
    assert len(validated) == 6
    assert all(model["coverage"] == model["identity"] == 1.0 and model["problems"] == [] for model in validated)


def test_new_producer_predictions_can_be_reused_transitively_without_search(hidden_models, monkeypatch):
    """Published pristine fields and query coverage survive two cache generations."""
    fixture, species, _ = hidden_models
    original, _, _ = make_plan(hidden_models)
    cli("run", "--output", original, "--cpus", 2)
    expected = json.loads((original / "rescued" / species[0] / "models.json").read_text())
    assert expected and all(m["raw_prediction"]["evidence"] == m["evidence"] for m in expected)
    previous = original
    fields = ("query", "cds", "sequence", "problems", "status", "raw_prediction", "coverage", "identity",
              "terminal_completion", "model_id", "quality_evidence")
    def checked_rows(models):
        # Exact-AA reuse can change record order and private predictor IDs.
        # Public model IDs, query provenance and all biological checks stay equal.
        rows = []
        for model in models:
            row = {key: model.get(key) for key in fields}
            row["raw_prediction"] = {k: v for k, v in row["raw_prediction"].items() if k != "id"}
            rows.append(json.dumps(row, sort_keys=True))
        return sorted(rows)
    for generation in (2, 3):
        output = fixture / ("generation_" + str(generation))
        cli("plan", "--cds-dir", fixture / "cds", "--gff-dir", fixture / "gff", "--genome-dir", fixture / "genome",
            "--busco-dir", fixture / "busco", "--tree", fixture / "tree.nwk", "--output", output,
            "--prediction-cache", previous)
        plan = rescue.load(output)
        for job in plan["synteny_jobs"]:
            if species[0] in {job["a"], job["b"]}:
                rescue.synteny(output, plan, job["index"], 2)
        def no_prediction(*args, **kwargs):
            raise AssertionError("A proven identical search was repeated")
        original_symlink = Path.symlink_to
        def no_predictor_genome(path, target, *args, original_symlink=original_symlink, **kwargs):
            if path.name == "genome.fa":
                raise AssertionError("A covered search materialized a predictor genome alias")
            return original_symlink(path, target, *args, **kwargs)
        with monkeypatch.context() as patch:
            patch.setattr(rescue, "search_intervals", lambda *args, **kwargs: [] if not args[1] else no_prediction())
            patch.setattr(rescue, "run", no_prediction)
            patch.setattr(Path, "symlink_to", no_predictor_genome)
            directory = rescue.rescue(output, plan, species[0], 2)
        actual = json.loads((directory / "models.json").read_text())
        assert checked_rows(actual) == checked_rows(expected)
        assert (directory / "genome_query_mapping.tsv").is_file()
        assert not (directory / "queries.fa").exists() and not (directory / "regions.fa").exists()
        assert not (directory / "genome.fa").exists() and not (directory / "genome.mpi").exists()
        assert not (directory / "logs").exists()
        assert json.loads((directory / "prediction_reuse.json").read_text())["searched_local_windows"] == 0
        previous = output
