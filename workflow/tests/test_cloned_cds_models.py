import csv
import gzip
import hashlib
import json
import random
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import pytest
from Bio.Seq import Seq

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
cloned = import_module("cloned_cds_models")
fsi = import_module("format_species_inputs")


def fixture(root, strand="+", duplicate=False, split=300, seed=812):
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    codons = [
        "".join((a, b, c)) for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")
    ]
    sequence = "ATG" + "".join(rng.choice(codons) for _ in range(240)) + "TAA"
    first = sequence[:split]
    second = sequence[split:]
    intron = "GT" + "".join(rng.choice("ACGT") for _ in range(80)) + "AG"
    region = first + intron + second
    if strand == "-":
        region = str(Seq(region).reverse_complement())
    genome = "C" * 100 + region + "C" * 100 + "ATGAAACCCTAA" + "C" * 100
    if duplicate:
        genome += "C" * 400 + region + "C" * 100
    paths = {
        role: root / name
        for role, name in dict(cds="source.cds.fa", gff="source.gff3", genome="genome.fa", clone="clones.fa").items()
    }
    paths["genome"].write_text(">chr1\n" + genome + "\n")
    paths["clone"].write_text(">experimental/clone:1\n" + sequence + "\n")
    start1, end1 = 101, 100 + split
    start2, end2 = 101 + split + len(intron), 100 + len(region)
    if strand == "-":
        start1, end1 = 101 + len(second) + len(intron), 100 + len(region)
        start2, end2 = 101, 100 + len(second)
    rows = ["##gff-version 3"]
    for gene, start, end in [
        ("N_fragment", start1, end1),
        ("C_fragment", start2, end2),
        ("neighbour", 201 + len(region), 212 + len(region)),
    ]:
        orientation = "+" if gene == "neighbour" else strand
        rows += [
            f"chr1\tsource\tgene\t{start}\t{end}\t.\t{orientation}\t.\tID={gene}",
            f"chr1\tsource\tmRNA\t{start}\t{end}\t.\t{orientation}\t.\tID={gene}.t1;Parent={gene}",
            f"chr1\tsource\tCDS\t{start}\t{end}\t.\t{orientation}\t0\tID={gene}.c1;Parent={gene}.t1",
        ]
    paths["gff"].write_text("\n".join(rows) + "\n")
    paths["cds"].write_text(
        ">N_fragment.t1\n" + first + "\n>C_fragment.t1\n" + second + "\n>neighbour.t1\nATGAAACCCTAA\n"
    )
    task = dict(
        provider="local",
        species_key="Test_species",
        species_prefix="Test_species",
        cds_path=paths["cds"],
        gff_path=paths["gff"],
        genome_path=paths["genome"],
        gene_grouping_mode="strict",
        gff_repair_mode="safe",
        format_strict=True,
        genetic_code=1,
    )
    row = dict(
        cloned_cds_fasta=str(paths["clone"]), cloned_cds_sha256=hashlib.sha256(paths["clone"].read_bytes()).hexdigest()
    )
    cloned.attach(task, row, root)
    return task, sequence, paths, row


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("split", [298, 299, 300])
def test_native_alignment_replaces_both_split_models_and_publishes_consistent_pair(tmp_path, strand, split):
    task, sequence, sources, _ = fixture(tmp_path / "inputs", strand, split=split)
    before = {p: p.read_bytes() for p in sources.values()}
    out = tmp_path / "formatted"
    out.mkdir()
    cds = fsi.format_cds(task, out, False, False)
    gff = fsi.format_gff(task, out, False, False, formatted_cds_path=cds["output_path"])
    genome = fsi.format_genome(task, out, False, False)
    report = json.loads(Path(task["cloned_model_audit"]).read_text())
    assert report["replacements"][0]["replaced_gene_ids"] == ["C_fragment", "N_fragment"]
    assert report["removed_CDS_ids"] == ["N_fragment.t1", "C_fragment.t1"]
    assert report["replacements"][0]["strand"] == strand
    assert report["original_sources_modified"] is False
    with gzip.open(cds["output_path"], "rt") as stream:
        fasta = stream.read()
    assert sequence in dict(fsi.iter_fasta_records(cds["output_path"])).values() and "fragment" not in fasta
    assert ">Test_species_neighbour\nATGAAACCCTAA\n" in fasta
    with gzip.open(gff["output_path"], "rt") as stream:
        annotation = stream.read()
    assert annotation.count("\tgene\t") == 2
    assert annotation.count("\tCDS\t") == 3
    assert "gg_replaces=C_fragment,N_fragment" in annotation
    from gene_model_catalog import indexed_genome, reconstruct_sequence

    model = report["replacements"][0]
    with indexed_genome(genome["output_path"]) as indexed:
        assert reconstruct_sequence(model["CDS_intervals_0based_halfopen"], "chr1", strand, indexed) == sequence
    from validate_longest_cds_selection import collect_expected_longest_records

    expected, _ = collect_expected_longest_records(task)
    assert len(expected) == 2
    from validate_cds_gff_mapping import validate_single_species

    result = validate_single_species(
        dict(
            index=1,
            species_prefix="Test_species",
            cds_file=cds["output_path"],
            gff_file=gff["output_path"],
            genome_file=genome["output_path"],
            strict=True,
        ),
        10,
    )
    assert result["ok"], result
    assert all(p.read_bytes() == data for p, data in before.items())


@pytest.mark.parametrize('full_model', [False, True])
def test_longer_only_genomic_policy_compares_owned_cds_and_binds_stage_cache(tmp_path, full_model):
    task, sequence, paths, row = fixture(tmp_path / 'inputs')
    if full_model:
        rows = paths['gff'].read_text().splitlines()
        end = next(line.split('\t')[4] for line in rows if line.endswith('\tID=C_fragment'))
        merged = ['##gff-version 3',
                  f'chr1\tsource\tgene\t101\t{end}\t.\t+\t.\tID=full',
                  f'chr1\tsource\tmRNA\t101\t{end}\t.\t+\t.\tID=full.t1;Parent=full']
        for line in rows:
            if '\tCDS\t' in line and 'fragment' in line:
                merged.append(line.replace('N_fragment', 'full').replace('C_fragment.c1', 'full.c2')
                              .replace('C_fragment.t1', 'full.t1'))
            elif 'neighbour' in line:
                merged.append(line)
        paths['gff'].write_text('\n'.join(merged) + '\n')
        paths['cds'].write_text('>full.t1\n' + sequence + '\n>neighbour.t1\nATGAAACCCTAA\n')
    original = dict(task)
    cloned.apply_to_task(task, tmp_path / 'curated')
    ordinary = Path(task['cloned_model_audit'])
    guarded = dict(original)
    cloned.attach(guarded, {**row, 'cloned_cds_longer_only_ids': '["experimental/clone:1"]'}, tmp_path)
    cloned.apply_to_task(guarded, tmp_path / 'curated')
    assert Path(guarded['cloned_model_audit']) != ordinary
    audit = json.loads(Path(guarded['cloned_model_audit']).read_text())
    if full_model:
        assert audit['replacements'] == []
        assert audit['skipped'][0]['status'] == 'not_longer'
        assert list(fsi.iter_fasta_records(guarded['cds_path']))[0] == ('full.t1', sequence)
    else:
        assert audit['skipped'] == []
        assert audit['replacements'][0]['replaced_gene_ids'] == ['C_fragment', 'N_fragment']


@pytest.mark.parametrize(
    "damage",
    [
        "duplicate_locus",
        "synonymous_reference_difference",
        "conflicting_clones",
        "unexplained_overlap",
        "noncoding_overlap",
        "missing_genome",
    ],
)
def test_rejects_unresolved_coordinate_or_model_identity_without_source_mutation(tmp_path, damage):
    task, sequence, paths, row = fixture(tmp_path / "inputs", duplicate=damage == "duplicate_locus")
    if damage == "synonymous_reference_difference":
        seq = str(Seq(sequence[:-3]).translate())
        for i in range(3, len(sequence) - 3, 3):
            changed = sequence[: i + 2] + ("C" if sequence[i + 2] != "C" else "T") + sequence[i + 3 :]
            if str(Seq(changed[:-3]).translate()) == seq:
                paths["clone"].write_text(">experimental/clone:1\n" + changed + "\n")
                break
    if damage == "conflicting_clones":
        paths["clone"].write_text(paths["clone"].read_text() + ">another_clone\n" + sequence + "\n")
    if damage == "unexplained_overlap":
        paths["gff"].write_text(paths["gff"].read_text().replace("\t101\t400\t", "\t98\t400\t"))
    if damage == "noncoding_overlap":
        with paths["gff"].open("a") as handle:
            handle.write("chr1\tsource\tgene\t120\t140\t.\t+\t.\tID=independent_noncoding\n")
    if damage == "missing_genome":
        task["genome_path"] = None
    row["cloned_cds_sha256"] = hashlib.sha256(paths["clone"].read_bytes()).hexdigest()
    cloned.attach(task, row, tmp_path)
    before = {p: p.read_bytes() for p in paths.values()}
    with pytest.raises(ValueError):
        cloned.apply_to_task(task, tmp_path / "curated")
    assert task["cds_path"] == paths["cds"]
    assert all(p.read_bytes() == data for p, data in before.items())


@pytest.mark.parametrize(
    "damage", ["partial", "internal_stop", "ambiguous", "duplicate_id", "unknown_selected", "bad_sha"]
)
def test_cloned_fasta_admission_is_explicit_and_strict(tmp_path, damage):
    task, sequence, paths, row = fixture(tmp_path)
    if damage == "partial":
        sequence = sequence[:-1]
    if damage == "internal_stop":
        sequence = sequence[:15] + "TAA" + sequence[18:]
    if damage == "ambiguous":
        sequence = sequence[:15] + "NNN" + sequence[18:]
    paths["clone"].write_text(">clone\n" + sequence + "\n")
    if damage == "duplicate_id":
        paths["clone"].write_text(paths["clone"].read_text() * 2)
    row["cloned_cds_sha256"] = hashlib.sha256(paths["clone"].read_bytes()).hexdigest()
    if damage == "bad_sha":
        row["cloned_cds_sha256"] = "0" * 64
    if damage == "unknown_selected":
        row["cloned_cds_ids"] = '["missing"]'
    with pytest.raises(ValueError):
        cloned.attach(task, row, tmp_path)
        cloned.cloned_records(paths["clone"], task["cloned_cds_ids"], 1)


def test_qualified_pair_reuse_is_byte_stable_and_new_process_rejects_changed_clone(tmp_path):
    task, _, paths, row = fixture(tmp_path / "inputs")
    cloned.apply_to_task(task, tmp_path / "curated")
    members = [task["cds_path"], task["gff_path"], Path(task["cloned_model_audit"])]
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in members}
    other, *_ = fixture(tmp_path / "inputs")
    cloned.apply_to_task(other, tmp_path / "curated")
    assert before == {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in members}
    paths["clone"].write_text(paths["clone"].read_text().replace("ATG", "ACG", 1))
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        cloned.attach({}, row, tmp_path)


@pytest.mark.parametrize("identity", ["gene_prefix", "locus_tag"])
def test_raw_gene_feature_ids_can_differ_from_native_cds_grouping_tokens(tmp_path, identity):
    task, sequence, paths, _ = fixture(tmp_path / "inputs")
    lines = []
    for line in paths["gff"].read_text().splitlines():
        if not line.startswith("#"):
            fields = line.split("\t")
            for gene in ("N_fragment", "C_fragment", "neighbour"):
                if fields[2] == "gene" and fields[8] == "ID=" + gene:
                    fields[8] = "ID=gene-" + gene
                    if identity == "locus_tag":
                        fields[8] += ";locus_tag=LOC_" + gene
                elif fields[2] == "mRNA":
                    fields[8] = fields[8].replace("Parent=" + gene, "Parent=gene-" + gene)
            line = "\t".join(fields)
        lines.append(line)
    paths["gff"].write_text("\n".join(lines) + "\n")
    cloned.apply_to_task(task, tmp_path / "curated")
    records = dict(fsi.iter_fasta_records(task["cds_path"]))
    assert len(records) == 2 and sequence in records.values()
    assert "N_fragment.t1" not in records and "C_fragment.t1" not in records
    audit = json.loads(Path(task["cloned_model_audit"]).read_text())
    assert audit["replacements"][0]["replaced_gene_ids"] == ["gene-C_fragment", "gene-N_fragment"]
    out = tmp_path / "formatted"
    out.mkdir()
    cds = fsi.format_cds(task, out, False, False)
    gff = fsi.format_gff(task, out, False, False, formatted_cds_path=cds["output_path"])
    genome = fsi.format_genome(task, out, False, False)
    from validate_cds_gff_mapping import validate_single_species

    result = validate_single_species(
        dict(
            index=1,
            species_prefix="Test_species",
            cds_file=cds["output_path"],
            gff_file=gff["output_path"],
            genome_file=genome["output_path"],
            strict=True,
        ),
        10,
    )
    assert result["ok"], result


@pytest.mark.parametrize("suffix,kind", [("", "gene"), (".t1", "mRNA"), (".t1.exon1", "exon"), (".t1.CDS1", "CDS")])
def test_generated_model_ids_cannot_collide_with_retained_annotation(tmp_path, suffix, kind):
    task, _, paths, _ = fixture(tmp_path / "inputs")
    gene = "gg_clone_" + hashlib.sha256(b"experimental/clone:1").hexdigest()[:16]
    with paths["gff"].open("a") as handle:
        phase = "0" if kind == "CDS" else "."
        handle.write(f"chr1\tsource\t{kind}\t1011\t1022\t.\t+\t{phase}\tID={gene}{suffix};Parent=neighbour\n")
    with pytest.raises(ValueError, match="collid"):
        cloned.apply_to_task(task, tmp_path / "curated")


def test_implicit_gtf_models_are_completed_and_replaced_without_source_mutation(tmp_path):
    task, sequence, paths, _ = fixture(tmp_path / "inputs")
    lines = []
    for line in paths["gff"].read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 9 and fields[2] == "CDS":
            transcript = fields[8].split("Parent=")[1]
            gene = transcript.removesuffix(".t1")
            fields[8] = f'gene_id "{gene}"; transcript_id "{transcript}";'
            lines.append("\t".join(fields))
    source = "\n".join(lines) + "\n"
    paths["gff"].write_text(source)
    cloned.apply_to_task(task, tmp_path / "curated")
    records = dict(fsi.iter_fasta_records(task["cds_path"]))
    assert len(records) == 2 and sequence in records.values()
    assert paths["gff"].read_text() == source
    assert task["gff_path"].read_text().count("\tgene\t") == 2


def test_source_graph_cycle_remains_an_error_with_shared_parent_cache(tmp_path):
    path = tmp_path / "cycle.gff3"
    path.write_text(
        "chr1\tsource\tmRNA\t1\t30\t.\t+\t.\tID=a;Parent=b\nchr1\tsource\tmRNA\t1\t30\t.\t+\t.\tID=b;Parent=a\n"
    )
    with pytest.raises(ValueError, match="Cyclic"):
        cloned.source_graph(path)


@pytest.mark.parametrize(
    "kind", ["gff3", "gtf", "gtf_percent", "gtf_noncoding", "missing_gene", "implicit_parent", "shared_parent"]
)
def test_streamed_gff_graph_preserves_native_normalization_contract(tmp_path, kind):
    from gene_model_refinement import _gff_source_graph, _selected_gff_source_rows

    if kind.startswith("gtf"):
        text = 'chr1\tsource\tCDS\t1\t30\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'
        if kind == "gtf_percent":
            text = text.replace('"g"', '"g%3BA"').replace('"t"', '"t%3BA"')
        if kind == "gtf_noncoding":
            text += 'chr1\tsource\texon\t40\t60\t.\t+\t.\tgene_id "nc"; transcript_id "nct";\n'
    else:
        text = (
            "##gff-version 3\n# metadata\n"
            "chr1\tsource\tgene\t1\t30\t.\t+\t.\tID=gene-g;gene_id=g\n"
            "chr1\tsource\tmRNA\t1\t30\t.\t+\t.\tID=t;Parent=gene-g;gene_id=g;transcript_id=t\n"
            "chr1\tsource\tCDS\t1\t30\t.\t+\t0\tParent=t\n"
        )
    if kind == "missing_gene":
        text = "\n".join(line for line in text.splitlines() if "\tgene\t" not in line) + "\n"
    if kind == "implicit_parent":
        text = text.replace("Parent=gene-g;gene_id=g;transcript_id=t", "gene_id=gene-g;transcript_id=t")
    if kind == "shared_parent":
        text = (
            "chr1\tsource\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\tsource\tmRNA\t1\t30\t.\t+\t.\tID=t2;Parent=g2\n"
            "chr1\tsource\tCDS\t1\t30\t.\t+\t0\tParent=t1,t2\n"
        )
    path = tmp_path / "annotation.gff3"
    path.write_text(text)
    expected = _gff_source_graph(_selected_gff_source_rows(_gff_source_graph(text), set(), retain_all=True))[0]
    records, _, _, _ = cloned.source_graph(path)
    assert records == expected


def test_synthetic_parent_export_is_complete_scope_only_and_does_not_edit_graph():
    import copy

    from gene_model_refinement import _gff_source_graph, _selected_gff_source_rows

    graph = _gff_source_graph('chr1\tsource\tCDS\t1\t30\t.\t+\t0\tgene_id "g"; transcript_id "t";\n')
    before = copy.deepcopy(graph)
    parents = _selected_gff_source_rows(graph, set(), retain_all=True, synthetic_parents_only=True)
    assert "\tgene\t" in parents and "\tmRNA\t" in parents and "\tCDS\t" not in parents
    assert graph == before
    with pytest.raises(ValueError, match="complete source graph"):
        _selected_gff_source_rows(graph, {"t"}, synthetic_parents_only=True)


def test_overlap_index_preserves_matches_errors_and_source_order():
    rng = random.Random(271)
    genes, coding, models = {}, {}, {}
    for i in range(200):
        start = rng.randrange(60) * 100
        genes[str(i)] = dict(seqid="chr" + str(i % 3), strand="+" if i % 2 else "-", start=start, end=start + 60)
        # Some source coding bounds deliberately disagree with gene bounds.
        coding[str(i)] = [(start + 400, start + 460)] if i % 4 == 0 else [(start, start + 60)] if i % 5 else []
    for i in range(60):
        gene = genes[str(i)]
        blocks = coding[str(i)] or [(gene["start"], gene["end"])]
        models[i] = dict(seqid=gene["seqid"], strand=gene["strand"], cds=[(s, e, 0) for s, e in blocks])
    indexed = cloned.overlap_gene_candidates(genes, coding, models)
    for key, model in models.items():
        try:
            expected = cloned.overlapping_genes(model, genes, coding)
        except ValueError as error:
            with pytest.raises(ValueError) as actual:
                cloned.overlapping_genes(model, indexed[key], coding)
            assert str(actual.value) == str(error)
        else:
            assert cloned.overlapping_genes(model, indexed[key], coding) == expected
    assert cloned.overlap_gene_candidates(genes, coding, {0: models[0]}) is None


def test_overlap_index_excludes_touching_boundaries_and_empty_strand_scope():
    genes = {
        name: dict(seqid="chr1", strand="+", start=s, end=e)
        for name, s, e in [("left", 40, 100), ("right", 160, 220), ("target", 100, 160)]
    }
    coding = {name: [(g["start"], g["end"])] for name, g in genes.items()}
    models = {strand: dict(seqid="chr1", strand=strand, cds=[(100, 160, 0)]) for strand in ("+", "-")}
    indexed = cloned.overlap_gene_candidates(genes, coding, models)
    assert list(indexed["+"]) == ["target"]
    assert not indexed["-"]
    assert cloned.overlapping_genes(models["+"], indexed["+"], coding) == ["target"]
    with pytest.raises(ValueError, match="No existing gene model"):
        cloned.overlapping_genes(models["-"], indexed["-"], coding)


def test_independent_clones_on_two_contigs_and_strands_replace_both_pairs(tmp_path):
    task, seq1, paths, row = fixture(tmp_path / "one")
    _, seq2, other, _ = fixture(tmp_path / "two", strand="-", seed=813)

    def rename(text):
        for source, target in [
            ("chr1", "chr2"),
            ("N_fragment", "N_second"),
            ("C_fragment", "C_second"),
            ("neighbour", "other_neighbour"),
        ]:
            text = text.replace(source, target)
        return text

    for role in ("genome", "gff", "cds"):
        with paths[role].open("a") as handle:
            handle.write(rename(other[role].read_text()))
    with paths["clone"].open("a") as handle:
        handle.write(">second_clone\n" + seq2 + "\n")
    row["cloned_cds_sha256"] = hashlib.sha256(paths["clone"].read_bytes()).hexdigest()
    cloned.attach(task, row, tmp_path)
    out = tmp_path / "formatted"
    out.mkdir()
    cds = fsi.format_cds(task, out, False, False)
    gff = fsi.format_gff(task, out, False, False, formatted_cds_path=cds["output_path"])
    genome = fsi.format_genome(task, out, False, False)
    sequences = dict(fsi.iter_fasta_records(cds["output_path"]))
    assert len(sequences) == 4 and seq1 in sequences.values() and seq2 in sequences.values()
    from validate_cds_gff_mapping import validate_single_species

    result = validate_single_species(
        dict(
            index=1,
            species_prefix="Test_species",
            cds_file=cds["output_path"],
            gff_file=gff["output_path"],
            genome_file=genome["output_path"],
            strict=True,
        ),
        10,
    )
    assert result["ok"], result
    audit = json.loads(Path(task["cloned_model_audit"]).read_text())
    assert {(r["seqid"], r["strand"]) for r in audit["replacements"]} == {("chr1", "+"), ("chr2", "-")}


@pytest.mark.parametrize("ambiguous", [False, True])
def test_cloned_mapping_precedes_annotation_but_valid_mapping_still_validates_gff(tmp_path, ambiguous):
    task, _, paths, _ = fixture(tmp_path / "inputs", duplicate=ambiguous)
    paths["gff"].write_text("not a usable GFF\n")
    error = "complete genome-identical paths" if ambiguous else "Invalid source GFF row"
    with pytest.raises(ValueError, match=error):
        cloned.apply_to_task(task, tmp_path / "curated")
    audits = list((tmp_path / "curated").rglob("audit.json"))
    assert len(audits) == 1 and json.loads(audits[0].read_text())["candidates"]
    assert paths["gff"].read_text() == "not a usable GFF\n"


def test_selected_clone_from_mixed_fasta_preserves_intronic_and_opposite_strand_genes(tmp_path):
    task, sequence, paths, row = fixture(tmp_path / "inputs")
    paths["clone"].write_text(paths["clone"].read_text() + ">published_partial\nATGA\n")
    row.update(
        cloned_cds_sha256=hashlib.sha256(paths["clone"].read_bytes()).hexdigest(),
        cloned_cds_ids='["experimental/clone:1"]',
    )
    for gene, start, end, strand in [("intronic_gene", 411, 422, "+"), ("antisense_gene", 111, 122, "-")]:
        with paths["gff"].open("a") as handle:
            handle.write(
                f"chr1\tsource\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene}\n"
                f"chr1\tsource\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gene}.t1;Parent={gene}\n"
                f"chr1\tsource\tCDS\t{start}\t{end}\t.\t{strand}\t0\tParent={gene}.t1\n"
            )
        with paths["cds"].open("a") as handle:
            handle.write(f">{gene}.t1\nATGAAACCCTAA\n")
    cloned.attach(task, row, tmp_path)
    cloned.apply_to_task(task, tmp_path / "curated")
    annotation = task["gff_path"].read_text()
    assert "ID=intronic_gene\n" in annotation and "ID=antisense_gene\n" in annotation
    assert set(dict(fsi.iter_fasta_records(task["cds_path"])).values()) == {sequence, "ATGAAACCCTAA"}
    assert json.loads(Path(task["cloned_model_audit"]).read_text())["replacements"][0]["replaced_gene_ids"] == [
        "C_fragment",
        "N_fragment",
    ]


def test_manifest_signature_and_frozen_plan_include_only_explicit_selected_clone_inputs(tmp_path):
    task, _, paths, row = fixture(tmp_path / "inputs")
    row.update(
        provider="local",
        id="Test_species",
        species_key="Test_species",
        bind_local_sources="1",
        cds_url=paths["cds"].as_uri(),
        gff_url=paths["gff"].as_uri(),
        genome_url=paths["genome"].as_uri(),
        cloned_cds_ids='["experimental/clone:1"]',
    )
    manifest = tmp_path / "inputs.tsv"
    with manifest.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    assert len(cloned.verify_manifest(manifest)) == 64
    plan = tmp_path / "plan.json"
    subprocess.run(
        [
            sys.executable,
            str(SUPPORT / "plan_input_generation_tasks.py"),
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "downloads"),
            "--stage-downloads",
            "--strict",
            "--outfile",
            str(plan),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    frozen = json.loads(plan.read_text())["tasks"][0]
    assert frozen["input_sha256"][str(paths["clone"])] == row["cloned_cds_sha256"]
    assert frozen["cloned_cds_ids"] == ["experimental/clone:1"]
    from local_input_manifest_binds import source_files

    assert paths["clone"] in source_files(manifest, tmp_path)
    assert paths["clone"] in source_files(plan, tmp_path, plan=True)
    from stage_input_generation_downloads import stage_downloads

    stage_downloads(plan, require_gff=True, require_genome=True)
    paths["clone"].write_text(paths["clone"].read_text() + "\n")
    with pytest.raises(ValueError, match="changed after planning"):
        stage_downloads(plan)


@pytest.mark.parametrize("mode", ["single", "array"])
@pytest.mark.parametrize("curation_mode", ["genomic", "sequence", "experimental", "experimental_new_locus"])
def test_native_strict_workflow_adopts_coupled_clone_models_with_restartable_receipts(tmp_path, mode, curation_mode):
    from test_gg_input_generation_end_to_end import (
        CORE_PATH,
        REPO_ROOT,
        _core_env,
        _install_fake_toolchain,
        _read_tsv_rows,
        _write_minimal_ete_taxonomy_db,
        _write_runtime_busco_dataset,
    )

    _task, sequence, paths, row = fixture(tmp_path / "sources")
    if curation_mode.startswith('experimental'):
        from test_cloned_cds_projection import change_clone, projected_fixture
        _task, sequence, paths, row, _ = projected_fixture(tmp_path / 'sources')
        sequence = sequence[:150] + 'GCT' + sequence[150:]
        change_clone(_task, paths, row, sequence)
        if curation_mode == 'experimental_new_locus':
            paths['gff'].write_text('\n'.join(line for line in paths['gff'].read_text().splitlines() if 'fragment' not in line) + '\n')
            paths['cds'].write_text('>neighbour.t1\nATGAAACCCTAA\n')
            targets = json.loads(row['cloned_cds_projection_targets'])
            targets['experimental/clone:1'].update(source_gene_ids=[], new_locus=True)
            row['cloned_cds_projection_targets'] = json.dumps(targets)
    row.update(
        provider="local",
        id="Arabidopsis_thaliana",
        species_key="Arabidopsis_thaliana",
        bind_local_sources="1",
        cds_url=paths["cds"].as_uri(),
        gff_url=paths["gff"].as_uri(),
        genome_url=paths["genome"].as_uri(),
        cloned_cds_ids='["experimental/clone:1"]',
    )
    if curation_mode == 'sequence':
        row.pop('gff_url')
        row.pop('genome_url')
        row['cloned_cds_mode'] = 'sequence'
    manifest = tmp_path / "cloned_inputs.tsv"
    with manifest.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    workspace = tmp_path / "workspace"
    _write_minimal_ete_taxonomy_db(workspace)
    _write_runtime_busco_dataset(workspace)
    fake_bin = _install_fake_toolchain(tmp_path)

    def run(phase, index=None, success=True, force=False):
        env = _core_env(
            workspace, None, fake_bin, phase, index, require_cds=True,
            require_gff=curation_mode != 'sequence', require_genome=curation_mode != 'sequence'
        )
        env.update(
            provider="all",
            download_manifest=str(manifest),
            strict="1",
            overwrite="1" if force else "0",
            run_species_taxonomy="1",
        )
        result = subprocess.run(
            ["bash", str(CORE_PATH)], cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=180
        )
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result

    root = workspace / "output/input_generation"
    if mode == "array":
        run("array_prepare")
        run("array_worker", 1)
        receipt = root / "tmp/task_plan.json.completed/1.json"
        assert receipt.is_file()
        files = json.loads(receipt.read_text())["files"]
        metadata = json.loads((root / "tmp/task_meta_shards/1.json").read_text())
        assert all(path in files for path in cloned.evidence_paths(metadata))
        before = {p: (Path(p).read_bytes(), Path(p).stat().st_mtime_ns) for p in files}
        run("array_worker", 1)
        assert before == {p: (Path(p).read_bytes(), Path(p).stat().st_mtime_ns) for p in files}
        run("array_finalize")
    else:
        run("single")
        run("single")
    summary = _read_tsv_rows(root / "gg_input_generation_species.tsv")
    assert len(summary) == 1
    clone_audit = Path(summary[0]["cloned_model_audit_path"])
    assert clone_audit.is_file()
    audit = json.loads(clone_audit.read_text())
    assert summary[0]['cloned_cds_mode'] == ('experimental' if curation_mode.startswith('experimental') else curation_mode)
    if curation_mode == 'sequence':
        assert audit['decisions'][0]['status'] == 'replace'
    else:
        assert len(audit["replacements"]) == 1
    cds = Path(summary[0]["cds_output_path"])
    records = dict(fsi.iter_fasta_records(cds))
    assert len(records) == 2 and sequence in records.values()
    assert not any("fragment" in name for name in records)
    if mode == "single":
        clone_audit.unlink()
        missing = run("single", success=False)
        assert "cloned model evidence" in missing.stderr.lower()
        run("single", force=True)
        assert clone_audit.is_file()
        paths["clone"].write_text(paths["clone"].read_text().replace("ATG", "ACG", 1))
        failed = run("single", success=False)
        assert "Cloned CDS SHA256 mismatch" in failed.stderr


def test_retain_source_policy_skips_only_unresolved_clone_and_keeps_default_strict(tmp_path):
    task, sequence, paths, row = fixture(tmp_path / 'inputs')
    paths['clone'].write_text(paths['clone'].read_text() + '>unmapped\nATG' + 'GCT' * 240 + 'TAA\n')
    row['cloned_cds_sha256'] = hashlib.sha256(paths['clone'].read_bytes()).hexdigest()
    cloned.attach(task, row, tmp_path)
    with pytest.raises(ValueError, match='complete genome-identical paths'):
        cloned.apply_to_task(dict(task), tmp_path / 'strict')
    cloned.attach(task, {**row, 'cloned_cds_unresolved_policy': 'retain_source'}, tmp_path)
    cloned.apply_to_task(task, tmp_path / 'retained')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert [x['clone_id'] for x in audit['replacements']] == ['experimental/clone:1']
    assert audit['skipped'][0]['clone_id'] == 'unmapped'
    assert audit['skipped'][0]['complete_genome_identical_paths'] == 0
    assert audit['genomic_CDS_identity_required'] is True


def test_contig_batching_preserves_coordinates_and_detects_cross_batch_paralogs(tmp_path):
    task, sequence, paths, row = fixture(tmp_path / 'inputs')
    reference = dict(task)
    cloned.apply_to_task(reference, tmp_path / 'whole')
    baseline = json.loads(Path(reference['cloned_model_audit']).read_text())
    original = paths['genome'].read_text().splitlines()[1]
    paths['genome'].write_text(paths['genome'].read_text() + '>chr2\n' + 'C' * 1500 + '\n')
    cloned.attach(task, {**row, 'cloned_cds_index_batch_bases': '1000'}, tmp_path)
    batched_task = dict(task)
    cloned.apply_to_task(task, tmp_path / 'parts')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['replacements'] == baseline['replacements']
    batches = json.loads((Path(task['cloned_model_audit']).parent / 'index_batches.json').read_text())
    assert len(batches['batches']) == 2
    assert batches['contigs_split'] is False
    paths['genome'].write_text('>chr1\n' + original + '\n>chr2\n' + original + '\n')
    with pytest.raises(ValueError, match='2 complete genome-identical paths'):
        cloned.apply_to_task(batched_task, tmp_path / 'duplicated')


@pytest.mark.parametrize('field,value', [('cloned_cds_unresolved_policy', 'ignore'),
                                         ('cloned_cds_index_batch_bases', '-1'),
                                         ('cloned_cds_index_batch_bases', '1.5')])
def test_genomic_search_configuration_rejects_invalid_values(tmp_path, field, value):
    _task, _sequence, _paths, row = fixture(tmp_path / 'inputs')
    with pytest.raises(ValueError, match=field):
        cloned.attach({}, {**row, field: value}, tmp_path)


@pytest.mark.parametrize('size', ['10000000', 'auto'])
def test_native_subset_manifest_explicit_options_preserve_source_rows(tmp_path, size):
    from input_generation_array_state import export_manifest
    row = dict(provider='local', id='Test_species', species_key='Test_species', cloned_cds_fasta='/clones.fa',
               cloned_cds_mode='genomic', cloned_cds_sha256='a' * 64)
    plan = dict(tasks=[dict(manifest_row=row)])
    target = tmp_path / 'manifest.tsv'
    options = dict(cloned_cds_unresolved_policy='retain_source', cloned_cds_index_batch_bases=size)
    export_manifest(plan, target, cloned_cds_options=options)
    result = next(csv.DictReader(target.open(), delimiter='\t'))
    assert result == {**row, **options}
    assert plan['tasks'][0]['manifest_row'] == row
    for bad in ({'cloned_cds_ids': '[]'}, {'cloned_cds_index_batch_bases': '-1'},
                {'cloned_cds_unresolved_policy': 'ignore'}):
        with pytest.raises(ValueError):
            export_manifest(plan, target, cloned_cds_options=bad)
    assert next(csv.DictReader(target.open(), delimiter='\t')) == result


@pytest.mark.parametrize('counts,capped', [([100], True), ([60, 60], False)])
def test_search_cap_applies_per_batch_after_global_coordinate_deduplication(tmp_path, monkeypatch, counts, capped):
    from rescue_gene_models import read_miniprot
    task, _sequence, _paths, _row = fixture(tmp_path / 'inputs')
    prepared = dict(task)
    cloned.apply_to_task(prepared, tmp_path / 'baseline')
    model_file = Path(prepared['cloned_model_audit']).parent / 'models.gff'
    raw = list(read_miniprot(model_file))[0]
    index = iter(counts)
    monkeypatch.setattr('cloned_cds_search.alignments', lambda *args, **kwargs: [model_file] * len(counts))
    monkeypatch.setattr('rescue_gene_models.read_miniprot', lambda _path: iter([raw] * next(index)))
    if capped:
        with pytest.raises(ValueError, match='complete genome-identical paths'):
            cloned.apply_to_task(task, tmp_path / 'tested')
    else:
        cloned.apply_to_task(task, tmp_path / 'tested')
        audit = json.loads(Path(task['cloned_model_audit']).read_text())
        assert len(audit['replacements']) == 1
        assert audit['search_limit_reached'] == []
