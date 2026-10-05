"""Held-out exon/path truth tested with real miniprot on temporary genomes."""

import hashlib
import json
import random
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

refinement = import_module("gene_model_refinement")
gff2genestat = import_module("gff2genestat")
fasta_records = import_module("fasta_sequence_store").fasta_records
write_tsv = import_module("test_gene_model_refinement").write_tsv

SCRIPT = Path(refinement.__file__)
GFF_COLS = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
OUT_COLS = ["gene_id", "feature_size", "num_intron", "gff_transcript_id", "feature_blocks"]


def test_refinement_cli_help_uses_real_runtime_without_outputs(tmp_path):
    result = subprocess.run([sys.executable, str(SCRIPT), "--help"], cwd=tmp_path,
                            text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower() and "verify-inputs" in result.stdout
    assert list(tmp_path.iterdir()) == []


def test_refinement_preserves_missing_models_from_completed_real_rescue(tmp_path):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    anchor_root, names, truth = fixture.make_plan(data)
    before = {path: path.read_bytes() for role in ('cds', 'gff', 'genome') for path in (tmp_path / role).iterdir()}
    fixture.cli('run', '--output', anchor_root, '--cpus', 2)
    models = json.loads((anchor_root / 'rescued' / names[0] / 'models.json').read_text())
    accepted = [model for model in models if model['status'] == 'accepted']
    assert len(accepted) == 1 and accepted[0]['sequence'] == truth[8]
    root = tmp_path / 'refinement'
    cli('plan', '--rescue-output', anchor_root, '--output', root, '--mode', 'off')
    value = refinement.load(root)
    source = value['request']['sources'][names[0]]
    assert Path(source['fasta']).is_relative_to(anchor_root / 'augmented')
    assert Path(source['gff']).is_relative_to(anchor_root / 'augmented')
    cli('run', '--output', root, '--cpus', 2)
    cli('qc', '--output', root)
    records = {identifier: sequence for identifier, _header, sequence in
               fasta_records(root / 'effective/species_cds' / (names[0] + '.fa'))}
    assert truth[8] in records.values()
    assert len(records) == 16  # Fifteen original annotations plus one missing model.
    assert accepted[0]['model_id'] in (root / 'effective/full_annotation' / (names[0] + '.gff3')).read_text()
    assert all(path.read_bytes() == content for path, content in before.items())
    augmented = anchor_root / 'augmented/species_cds' / (names[0] + '.rescue.cds.fa')
    augmented.write_text(augmented.read_text() + '\n')
    with pytest.raises(ValueError, match='Frozen refinement input changed'):
        refinement.load(root)
    with pytest.raises(ValueError, match='Augmented rescue inputs incomplete or corrupt'):
        refinement.plan(tmp_path / 'new-refinement', rescue_output=anchor_root, mode='off')


def coding(seed, residues=280):
    rng = random.Random(seed)
    codons = {}
    for codon, amino_acid in sorted(CodonTable.unambiguous_dna_by_id[1].forward_table.items()):
        codons.setdefault(amino_acid, codon)
    amino_acids = sorted(codons)
    return "ATG" + "".join(codons[rng.choice(amino_acids)] for _ in range(residues - 1)) + "TAA"


def spliced_gene(sequence, strand="+", flank=1100):
    pieces = [sequence[:300], sequence[300:540], sequence[540:]]
    intron = "GT" + "AC" * 38 + "AG"
    genomic = "AC" * (flank // 2)
    blocks = []
    for index, piece in enumerate(pieces):
        start = len(genomic)
        genomic += piece
        blocks.append([start, len(genomic), 0])
        if index + 1 < len(pieces):
            genomic += intron
    genomic += "CA" * (flank // 2)
    if strand == "-":
        size = len(genomic)
        blocks = [[size - end, size - start, phase] for start, end, phase in blocks]
        genomic = str(Seq(genomic).reverse_complement())
    return genomic, blocks


def annotation(contig, gene, transcript, blocks, strand="+", attributes=""):
    start, end = min(block[0] for block in blocks), max(block[1] for block in blocks)
    lines = [f"{contig}\tsynthetic\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\tID={gene}{attributes}",
             f"{contig}\tsynthetic\tmRNA\t{start + 1}\t{end}\t.\t{strand}\t.\tID={transcript};Parent={gene}"]
    for number, (a, b, phase) in enumerate(blocks):
        lines += [f"{contig}\tsynthetic\tCDS\t{a + 1}\t{b}\t.\t{strand}\t{phase}\tID={transcript}.c{number};Parent={transcript}",
                  f"{contig}\tsynthetic\texon\t{a + 1}\t{b}\t.\t{strand}\t.\tID={transcript}.e{number};Parent={transcript}"]
    return "\n".join(lines) + "\n"


def truth_fixture(tmp_path, *, rna=False):
    """Target annotations lose exons/paths; DNA truth and donors stay complete."""
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    kinds = ["repair", "minus", "alternate", "intact", "pseudogene", "stop", "gap"]
    truth, sources = {}, []
    for name in names:
        genome, gff, fasta = [], ["##gff-version 3\n"], []
        for number, kind in enumerate(kinds):
            strand = "-" if kind == "minus" else "+"
            sequence = coding(104 + number)
            if name == names[0] and kind == "stop":
                sequence = sequence[:390] + "TGA" + sequence[393:]
            if name == names[0] and kind == "gap":
                sequence = sequence[:390] + "NNN" + sequence[393:]
            dna, full_blocks = spliced_gene(sequence, strand)
            contig = "chr_" + kind
            genome.append(f">{contig}\n{dna}\n")
            target = name == names[0]
            blocks = full_blocks[:2] if target and kind in {"repair", "minus"} else full_blocks
            # Provider parent IDs and formatter locus tokens are deliberately
            # different; original and predicted paths must share both correctly.
            source_gene = "provider_" + kind
            attributes = ";gene_id=" + kind + (";pseudo=true" if target and kind == "pseudogene" else "")
            gff.append(annotation(contig, source_gene, kind + ".long", blocks, strand, attributes))
            supplied = "".join(sequence[a:b] for a, b in ((0, 300), (300, 540))) if target and kind in {"repair", "minus"} else sequence
            fasta.append(f">{kind}.long\n{supplied}\n")
            if kind == "alternate" and not target:
                short = [full_blocks[0], full_blocks[2]]
                # An exon-skipped coding path exists in both donor annotations.
                gff.append(annotation(contig, source_gene, kind + ".short", short, strand).split("\n", 1)[1])
                fasta.append(f">{kind}.short\n{sequence[:300] + sequence[540:]}\n")
            truth[name, kind] = dict(cds=sequence, blocks=full_blocks, strand=strand,
                                    shorter_cds=sequence[:300] + sequence[540:])
        paths = {role: tmp_path / (name + suffix) for role, suffix in
                 (("cds", ".cds.fa"), ("gff", ".gff3"), ("genome", ".genome.fa"))}
        for role, content in (("cds", fasta), ("gff", gff), ("genome", genome)):
            paths[role].write_text("".join(content))
        sources.append(dict(species=name, genetic_code=1, **{role: str(path) for role, path in paths.items()}))
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    write_tsv(inputs, list(sources[0]), sources)
    links = [dict(species_a=names[i], gene_a=names[i] + "_" + kind,
                  species_b=names[j], gene_b=names[j] + "_" + kind, weight=1)
             for kind in kinds for i, j in ((0, 1), (0, 2), (1, 2))]
    write_tsv(edges, list(links[0]), links)
    rna_path = None
    if rna:
        full = truth[names[0], "alternate"]["blocks"]
        rna_path = tmp_path / "rna.tsv"
        evidence = dict(species=names[0], seqid="chr_alternate", strand="+",
                        cds_blocks=json.dumps([[a, b] for a, b, _ in (full[0], full[2])]),
                        transcript_id="target_short_rna", count=12)
        write_tsv(rna_path, list(evidence), [evidence])
    return inputs, edges, rna_path, sources, truth


def cli(*args, success=True):
    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], text=True, capture_output=True)
    if success:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0, result.stdout + result.stderr
    return result


@pytest.mark.parametrize("rna", [False, True])
def test_real_miniprot_repairs_heldout_exons_adds_paths_and_preserves_biological_negatives(tmp_path, rna):
    inputs, edges, evidence, sources, truth = truth_fixture(tmp_path, rna=rna)
    original = {row[role]: Path(row[role]).read_bytes() for row in sources for role in ("cds", "gff", "genome")}
    root = tmp_path / "run"
    params = ["plan", "--inputs", inputs, "--edges", edges, "--output", root,
              "--padding", 650, "--max-intron", 2000]
    if evidence:
        params += ["--rna", evidence]
    cli(*params)
    for command in ("catalog", "correspondence", "select", "predict", "finalize", "qc"):
        cli(command, "--output", root, "--cpus", 2)
    target = "Species_target"
    predictions = json.loads((root / "predictions" / target / "predictions.json").read_text())
    accepted = {row["gene_id"].removeprefix(target + "_"): row for row in predictions if row["status"] == "accepted"}
    for kind in ("repair", "minus"):
        assert kind in accepted, [(row["gene_id"], row["status"], row["problems"]) for row in predictions]
        model = accepted[kind]
        assert model["candidate"]["cds"] == truth[target, kind]["cds"]
        assert model["candidate"]["blocks"] == truth[target, kind]["blocks"]
        assert model["candidate"]["strand"] == truth[target, kind]["strand"]
        assert model["change_type"] == "model_revision"
        assert model["donors"] == ["Species_donor1", "Species_donor2"]
    assert "alternate" in accepted
    alternative = accepted["alternate"]
    assert alternative["candidate"]["cds"] == truth[target, "alternate"]["shorter_cds"]
    assert alternative["change_type"] == "isoform_addition"
    assert alternative["evidence_class"] == ("rna_path_supported" if rna else "homology_only_predicted")
    assert alternative["candidate"]["quality"]["representative_eligible"] is rna
    assert not ({"intact", "pseudogene", "stop", "gap"} & set(accepted))
    effective = root / "effective"
    summary = json.loads((root / "catalog_index_final/summary.json").read_text())
    final_db = root / "catalog_index_final/loci.sqlite3"
    candidates = sum(len(locus['candidates']) for locus in refinement.iter_loci(final_db))
    assert Path(summary['database']) == final_db
    assert summary['candidate_count'] == candidates
    assert summary['accepted_prediction_count'] >= len(accepted)
    selected = {identifier: sequence for identifier, _header, sequence in
                fasta_records(effective / "species_cds" / (target + ".fa"))}
    proteins = {identifier: sequence for identifier, _header, sequence in
                fasta_records(effective / "species_protein" / (target + ".fa"))}
    assert not ({target + "_" + kind for kind in ("pseudogene", "stop", "gap")} & set(proteins))
    for kind in ("repair", "minus", "intact", "alternate"):
        assert target + "_" + kind in proteins
        assert proteins[target + "_" + kind] == str(Seq(selected[target + "_" + kind]).translate()).removesuffix("*")
    assert not any("*" in sequence for sequence in proteins.values())
    for kind in ("repair", "minus", "intact", "pseudogene", "stop", "gap"):
        assert selected[target + "_" + kind] == truth[target, kind]["cds"]
    if not rna:
        assert selected[target + "_alternate"] == truth[target, "alternate"]["cds"]
    full = (effective / "full_annotation" / (target + ".gff3")).read_text()
    assert all("ID=" + kind + ".long;" in full for kind in ("repair", "minus", "alternate"))
    assert "ID=" + alternative["candidate"]["source_transcript_id"] + ";" in full
    reimported = refinement.build_catalog(target, effective / "all_candidates" / (target + ".fa"),
                                          effective / "full_annotation" / (target + ".gff3"),
                                          effective / "species_genome" / (target + ".fa"))
    rebound = {candidate["source_transcript_id"]: (locus, candidate) for locus in reimported["loci"]
               for candidate in locus["candidates"]}
    for kind in ("repair", "minus", "alternate"):
        for transcript in (kind + ".long", accepted[kind]["candidate"]["source_transcript_id"]):
            locus, candidate = rebound[transcript]
            assert locus["gene_id"] == target + "_" + kind
            assert candidate["gene_token"] == kind
            assert candidate["source_gene_id"] == "provider_" + kind
    assert (effective / "source_annotation" / (target + ".gff3")).read_bytes() == original[sources[0]["gff"]]
    # Downstream reader reports the very transcript/coordinates behind selected CDS.
    choices = {row["gene_id"]: row for row in refinement.read_table(effective / "representative_map.tsv") if row["species"] == target}
    traits = gff2genestat.process_single_gff(target + ".gff3", str(effective / "species_gff"),
                                           list(selected), "CDS", "longest", GFF_COLS, OUT_COLS,
                                           representative_map=effective / "representative_map.tsv")
    for row in traits.itertuples(index=False):
        assert row.gff_transcript_id == choices[row.gene_id]["source_transcript_id"]
        assert row.feature_size == len(selected[row.gene_id])
    assert all(Path(path).read_bytes() == content for path, content in original.items())
    assert not list(tmp_path.glob("*.fai"))
    receipts = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("receipt.json")}
    cli("run", "--output", root, "--cpus", 2)
    assert receipts == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("receipt.json")}
    cli("qc", "--output", root)


def test_real_miniprot_wrong_neighbor_hits_are_proposals_not_repairs(tmp_path):
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    sequence = coding(800, residues=300)
    donor_dna, donor_blocks = spliced_gene(sequence)
    rows = []
    for name in names:
        cds, gff, genome = [tmp_path / (name + suffix) for suffix in (".cds.fa", ".gff3", ".genome.fa")]
        if name == names[0]:
            # The owned locus is physically truncated. A nearby annotated duplicate
            # carries the matching full protein and lies inside the search window.
            partial = sequence[:540]
            clone_start = 1100
            dna = "AC" * 100 + partial + "AC" * ((clone_start - 200 - len(partial)) // 2) + sequence + "AC" * 100
            genome.write_text(">chr1\n" + dna + "\n")
            cds.write_text(">wrong.t\n" + partial + "\n>neighbor.t\n" + sequence + "\n")
            gff.write_text("##gff-version 3\n" + annotation("chr1", "wrong", "wrong.t", [[200, 740, 0]])
                           + annotation("chr1", "neighbor", "neighbor.t", [[clone_start, clone_start + len(sequence), 0]]))
        else:
            genome.write_text(">chr1\n" + donor_dna + "\n")
            cds.write_text(">wrong.t\n" + sequence + "\n")
            gff.write_text("##gff-version 3\n" + annotation("chr1", "wrong", "wrong.t", donor_blocks))
        rows.append(dict(species=name, cds=str(cds), gff=str(gff), genome=str(genome), genetic_code=1))
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    write_tsv(inputs, list(rows[0]), rows)
    links = [dict(species_a=names[0], gene_a=names[0] + "_wrong", species_b=donor, gene_b=donor + "_wrong")
             for donor in names[1:]]
    write_tsv(edges, list(links[0]), links)
    root = tmp_path / "run"
    cli("plan", "--inputs", inputs, "--edges", edges, "--output", root, "--padding", 1600)
    cli("predict", "--output", root, "--task-index", 3 if sorted(names).index(names[0]) == 2 else sorted(names).index(names[0]) + 1)
    predictions = json.loads((root / "predictions" / names[0] / "predictions.json").read_text())
    assert predictions
    assert not any(row["status"] == "accepted" for row in predictions)
    assert any("overlap_other_locus" in row["problems"] for row in predictions)
    assert any("no_overlap_owned_locus" in row["problems"] for row in predictions)


def test_real_synteny_reuse_recovers_correspondence_for_existing_anchor_excluded_model(tmp_path):
    """A masked final CDS base removes an anchor while its two flanks survive."""
    rescue = refinement.rescue
    names = ["Plant_species0", "Plant_species1", "Plant_species2"]
    truth = [coding(1400 + number, residues=240) for number in range(7)]
    for role in ("cds", "gff", "genome", "busco"):
        (tmp_path / role).mkdir()
    for name in names:
        genome, annotations, records = "", ["##gff-version 3\n"], []
        for number, sequence in enumerate(truth):
            start = len(genome)
            genome += sequence + "AC" * 60
            end = start + len(sequence) - (1 if name == names[0] and number == 3 else 0)
            annotations.append(annotation("chr1", "g" + str(number), "g" + str(number) + ".t", [[start, end, 0]]))
            records.append(f">{name}_g{number}\n{sequence[:end - start]}\n")
        (tmp_path / "cds" / (name + ".cds.fa")).write_text("".join(records))
        (tmp_path / "gff" / (name + ".gff3")).write_text("".join(annotations))
        (tmp_path / "genome" / (name + ".genome.fa")).write_text(">chr1\n" + genome + "\n")
        (tmp_path / "busco" / (name + ".busco.short.txt")).write_text(
            "# BUSCO version is: 6.0.0\n# The lineage dataset is: embryophyta_odb12\n"
            "# BUSCO was run in mode: transcriptome\nC:95.0%[S:95.0%,D:0.0%],F:0.0%,M:5.0%,n:100\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(Plant_species0:1,(Plant_species1:1,Plant_species2:1):1);\n")
    anchor_root, root = tmp_path / "anchors", tmp_path / "refinement"
    command = [sys.executable, rescue.__file__, "plan", "--cds-dir", str(tmp_path / "cds"),
               "--gff-dir", str(tmp_path / "gff"), "--genome-dir", str(tmp_path / "genome"),
               "--busco-dir", str(tmp_path / "busco"), "--tree", str(tree), "--output", str(anchor_root),
               "--common-references", "2", "--nearest-references", "1", "--min-anchors", "3"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    anchors = rescue.load(anchor_root)
    jobs = [job for job in anchors["synteny_jobs"] if job["a"] != job["b"]]
    assert len(jobs) <= 3
    for job in jobs:
        rescue.synteny(anchor_root, anchors, job["index"], 2, tmp_path / "comparisons")
    receipts = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (anchor_root / "synteny").rglob("receipt.json")}
    mapping = rescue.table(anchor_root / "prepared" / names[0] / "genes.id_map.tsv")
    excluded = next(row for row in mapping if row["original_id"] == names[0] + "_g3")
    assert excluded["status"] == "translation_excluded"
    cli("plan", "--rescue-output", anchor_root, "--output", root, "--padding", 150)
    cli("correspondence", "--output", root, "--cpus", 2, "--comparison-cache", tmp_path / "comparisons")
    actual = json.loads((root / "correspondence" / "edges.json").read_text())
    restored = [edge for edge in actual if edge["gene_a"].endswith("_g3") and edge["gene_b"].endswith("_g3")
                and names[0] in {edge["species_a"], edge["species_b"]}]
    assert len(restored) == 2, actual
    assert all(edge["evidence"]["kind"] == "two_flanking_anchors" and not edge["ambiguous"] for edge in restored)
    assert receipts == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (anchor_root / "synteny").rglob("receipt.json")}
    cli("predict", "--output", root, "--task-index", 1, "--cpus", 2)
    predictions = json.loads((root / "predictions" / names[0] / "predictions.json").read_text())
    repaired = [row for row in predictions if row["gene_id"] == names[0] + "_g3" and row["status"] == "accepted"]
    assert len(repaired) == 1, predictions
    assert repaired[0]["candidate"]["cds"] == truth[3]
    assert repaired[0]["donors"] == names[1:]
