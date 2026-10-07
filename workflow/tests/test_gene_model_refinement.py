"""Frozen refinement contracts, ownership gates, and effective reader agreement."""

import copy
import csv
import gzip
import hashlib
import json
import shutil
import sqlite3
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

refinement = import_module("gene_model_refinement")
validate_candidate = import_module("gene_model_catalog").validate_candidate


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def tiny_inputs(tmp_path):
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    rows = []
    for name in names:
        cds, gff, genome = [tmp_path / (name + suffix) for suffix in (".cds.fa", ".gff3", ".genome.fa")]
        cds.write_text(">g\nATGAAACCCTAA\n")
        genome.write_text(">chr1\nATGAAACCCTAA\n")
        gff.write_text("##gff-version 3\n"
                       "chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g\n"
                       "chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g\n"
                       "chr1\ts\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=t1\n"
                       "chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t2;Parent=g\n"
                       "chr1\ts\tCDS\t1\t6\t.\t+\t0\tID=c2;Parent=t2\n"
                       "chr1\ts\tCDS\t10\t12\t.\t+\t0\tID=c3;Parent=t2\n")
        rows.append(dict(species=name, cds=str(cds), gff=str(gff), genome=str(genome), genetic_code=1))
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    write_tsv(inputs, list(rows[0]), rows)
    links = [dict(species_a=names[i], gene_a=names[i] + "_g", species_b=names[j], gene_b=names[j] + "_g")
             for i, j in ((0, 1), (0, 2), (1, 2))]
    write_tsv(edges, list(links[0]), links)
    return inputs, edges, rows


def classification_fixture(*, valid_original=False, other_locus=False):
    sequence = "ATG" + "AAA" * 18 + "TAA"
    original = {"candidate_id": "Species_target_t", "source_transcript_id": "t", "source_gene_id": "g",
                "seqid": "chr1", "strand": "+", "blocks": [[100, 130, 0]], "cds": sequence[:30],
                "origin": "original", "quality": {"valid_orf": valid_original}}
    locus = {"gene_id": "Species_target_g", "species": "Species_target", "seqid": "chr1", "strand": "+",
             "candidates": [original]}
    loci = [locus]
    if other_locus:
        loci.append({**copy.deepcopy(locus), "gene_id": "Species_target_neighbor",
                     "candidates": [{**copy.deepcopy(original), "blocks": [[150, 180, 0]]}]})
    catalog = {"species": "Species_target", "genetic_code": 1, "loci": loci}
    models = [dict(gene_id="Species_target_g", seqid="chr1", strand="+", cds=[[100, 160, 0]],
                   sequence=sequence, problems=[], donor_species=donor, donor_candidate=donor + "_t",
                   identity=1.0, coverage=1.0) for donor in ("Species_donor1", "Species_donor2")]
    edges = [dict(species_a="Species_target", gene_a="Species_target_g", species_b=donor,
                  gene_b=donor + "_g", weight=1.0, ambiguous=False)
             for donor in ("Species_donor1", "Species_donor2")]
    return catalog, models, edges


def classify(catalog, models, edges, rna=(), **params):
    return refinement.classify_predictions(models, catalog, edges, {**refinement.DEFAULTS, **params},
                                           list(rna), "assembly-sha256")


@pytest.mark.parametrize("problem", ["low_identity", "low_coverage"])
def test_weak_extra_donor_cannot_veto_two_qualifying_donors(problem):
    catalog, models, edges = classification_fixture()
    weak = {**copy.deepcopy(models[0]), "donor_species": "Species_weak",
            "donor_candidate": "Species_weak_t", "problems": [problem]}
    edges.append({**edges[0], "species_b": "Species_weak"})
    result = classify(catalog, models + [weak], edges)[0]
    assert result["status"] == "accepted"
    assert result["donors"] == ["Species_donor1", "Species_donor2"]
    assert result["alignments"][-1]["supports_path"] is False
    assert result["alignments"][-1]["problems"] == [problem]


def test_weak_donor_cannot_supply_independent_support():
    catalog, models, edges = classification_fixture()
    models[1]["problems"] = ["low_identity"]
    result = classify(catalog, models, edges)[0]
    assert result["status"] == "proposal"
    assert result["donors"] == ["Species_donor1"]
    assert result["problems"] == ["insufficient_independent_support"]


@pytest.mark.parametrize("problem", ["frameshift", "noncanonical_splice", "assembly_gap", "outside_search_window"])
def test_structural_defect_still_vetoes_shared_path(problem):
    catalog, models, edges = classification_fixture()
    models[1]["problems"] = [problem]
    result = classify(catalog, models, edges)[0]
    assert result["status"] == "proposal"
    assert problem in result["problems"]


def test_source_coding_path_phase_inference_reaches_verified_analysis_gff_and_protein(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    Path(target['gff']).write_text('##gff-version 3\n'
                                 'chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g\n'
                                 'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g\n'
                                 'chr1\ts\tCDS\t1\t12\t.\t+\t.\tID=c1;Parent=t1\n'
                                 'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t2;Parent=g\n'
                                 'chr1\ts\tCDS\t1\t12\t.\t+\t.\tID=c2;Parent=t2\n')
    root = tmp_path / 'refinement'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    refinement.finalize(root, value)
    refinement.verify_inputs(root / 'effective/inputs.tsv')
    assert (root / 'effective/species_protein/Species_target.fa').read_text() == '>Species_target_g\nMKP\n'
    assert '\tCDS\t1\t12\t.\t+\t0\t' in (root / 'effective/analysis_gff/Species_target.gff3').read_text()
    assert (root / 'effective/source_annotation/Species_target.gff3').read_text() == Path(target['gff']).read_text()
    locus = json.loads((root / 'catalog/Species_target/loci.jsonl').read_text())
    assert not locus['source_baseline_candidate_id']
    assert all(c['source_fasta_ids'] == [] for c in locus['candidates'])
    assert all(c['quality']['phase_inference_evidence'] == 'complete_genomic_cds_and_unique_source_coding_path' for c in locus['candidates'])
    review = import_module('plot_gene_model_refinement').collect(root, max_loci=2)
    assert review['species']['Species_target']['coding_path_phase_resolved_representatives'] == 1
    assert review['coding_path_phase_loci_available'] == 1
    assert [(r['species'], r['gene_id']) for r in review['details']] == [('Species_target', 'Species_target_g')]
    assert all(c['source_blocks'] == [[0, 12, -1]] for c in review['details'][0]['candidates'])


def test_gene_only_mismatch_is_archived_but_not_exported_as_genomic_representative(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    Path(target["gff"]).write_text("##gff-version 3\n"
                                 "chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g\n"
                                 "chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g\n"
                                 "chr1\ts\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=t1\n")
    original = ">g\nATGCCCCCCTAA\n"
    Path(target["cds"]).write_text(original)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    assert (effective / "source_cds" / (target["species"] + ".fa")).read_text() == original
    assert (effective / "species_cds" / (target["species"] + ".fa")).read_text() == ""
    assert (effective / "species_protein" / (target["species"] + ".fa")).read_text() == ""
    assert "source_cds_sequence_mismatch" in (effective / "effective_exclusions.tsv").read_text()
    assert "ID=t1;Parent=g" in (effective / "source_annotation" / (target["species"] + ".gff3")).read_text()
    assert refinement.verify_inputs(effective / "inputs.tsv")


@pytest.mark.parametrize("kind", ["five_prime_partial", "masked_iupac"])
def test_explained_formatter_difference_retains_dna_with_separate_translation_admission(tmp_path, kind):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    partial = kind == "five_prime_partial"
    genomic = "CCATGAAACCC" if partial else "ATGGRGTAA"
    supplied = "NNNCATGAAACCC" if partial else "ATGGNGTAA"
    coding = genomic[2:] if partial else genomic
    n = len(genomic)
    Path(target["genome"]).write_text(">chr1\n" + genomic + "\n")
    original = ">g\n" + supplied + "\n"
    Path(target["cds"]).write_text(original)
    Path(target["gff"]).write_text("##gff-version 3\n"
                                 f"chr1\ts\tgene\t1\t{n}\t.\t+\t.\tID=g\n"
                                 f"chr1\ts\tmRNA\t1\t{n}\t.\t+\t.\tID=t;Parent=g\n"
                                 f"chr1\ts\texon\t1\t{n}\t.\t+\t.\tID=e;Parent=t\n"
                                 f"chr1\ts\tCDS\t{3 if partial else 1}\t{n}\t.\t+\t0\tID=c;Parent=t\n")
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    species = target["species"]
    assert (effective / "source_cds" / (species + ".fa")).read_text() == original
    assert (effective / "species_cds" / (species + ".fa")).read_text() == f">{species}_g\n{coding}\n"
    assert (effective / "source_annotation" / (species + ".gff3")).read_text() == Path(target["gff"]).read_text()
    assert not [r for r in refinement.read_table(effective / "effective_exclusions.tsv") if r["species"] == species]
    admitted = next(r for r in refinement.read_table(effective / "translation_admission.tsv") if r["species"] == species)
    quality = json.loads(admitted["quality"])
    assert not quality["sequence_mismatch"]
    assert quality["partial"] is partial
    assert quality["ambiguous"] is not partial
    assert admitted["status"] == ("included" if partial else "excluded")
    assert bool((effective / "species_protein" / (species + ".fa")).read_text()) is partial
    assert refinement.verify_inputs(effective / "inputs.tsv")


@pytest.mark.parametrize('code,dual', [(27, 'TGA'), (28, 'TAA'), (31, 'TAA')])
def test_whole_rna_cannot_adopt_predicted_dual_coding_translation_context(code, dual):
    catalog, models, edges = classification_fixture()
    catalog['genetic_code'] = code
    for model in models:
        model['sequence'] = model['sequence'][:-3] + dual
    rna = [dict(species='Species_target', seqid='chr1', strand='+', cds_blocks='[[100,160]]',
                transcript_id='whole-target-RNA', count=5)]
    result = classify(catalog, models, edges, rna= rna)
    assert len(result) == 1
    row = result[0]
    assert row['rna_paths'] == ['whole-target-RNA']
    assert row['candidate']['quality']['translation_uncertain']
    assert not row['candidate']['quality']['valid_orf']
    assert row['status'] == 'proposal'
    assert 'invalid_predicted_coding_path' in row['problems']
    assert not row['candidate']['quality']['representative_eligible']


@pytest.mark.parametrize('code,dual', [(27, 'TGA'), (28, 'TAA'), (31, 'TAA')])
def test_effective_withholds_dual_context_protein_and_analysis_but_preserves_dna_gff_source(tmp_path, code, dual):
    inputs, edges, rows = tiny_inputs(tmp_path)
    sequence = 'ATGAAA' + dual
    before = {}
    for row in rows:
        row['genetic_code'] = code
        Path(row['cds']).write_text('>g\n' + sequence + '\n')
        Path(row['genome']).write_text('>chr1\n' + sequence + '\n')
        Path(row['gff']).write_text('##gff-version 3\nchr1\ts\tgene\t1\t9\t.\t+\t.\tID=g\n'
                                    'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n'
                                    'chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=t\n')
        for role in ('cds', 'genome', 'gff'):
            before[row[role]] = Path(row[role]).read_bytes()
    write_tsv(inputs, list(rows[0]), rows)
    evidence = [dict(species=row['species'], seqid='chr1', strand='+', cds_blocks='[[0,9]]',
                     transcript_id='whole-RNA', count=5) for row in rows]
    rna = tmp_path / 'rna.tsv'
    write_tsv(rna, list(evidence[0]), evidence)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, rna=rna, mode='off')
    effective = refinement.finalize(root, value)
    for row in rows:
        name = row['species']
        assert sequence in (effective / 'species_cds' / (name + '.fa')).read_text()
        assert '\tCDS\t1\t9\t' in (effective / 'species_gff' / (name + '.gff3')).read_text()
        assert (effective / 'species_protein' / (name + '.fa')).read_text() == ''
        assert (effective / 'analysis_cds' / (name + '.fa')).read_text() == ''
        assert '\tCDS\t' not in (effective / 'analysis_gff' / (name + '.gff3')).read_text()
        assert sequence in (effective / 'all_candidates' / (name + '.fa')).read_text()
        assert (effective / 'source_cds' / (name + '.fa')).read_bytes() == before[row['cds']]
        assert (effective / 'source_annotation' / (name + '.gff3')).read_bytes() == before[row['gff']]
        assert '\tCDS\t1\t9\t' in (effective / 'full_annotation' / (name + '.gff3')).read_text()
    audit = refinement.read_table(effective / 'translation_admission.tsv')
    assert all(row['status'] == 'excluded' and json.loads(row['quality'])['translation_uncertain']
               and json.loads(row['quality'])['rna_supported'] for row in audit)
    assert all(Path(path).read_bytes() == content for path, content in before.items())


def analysis_source(sequence, blocks, *, owner="t", strand="+", code=1):
    candidate = dict(candidate_id="Species_" + owner, source_transcript_id=owner, source_gene_id="g",
                     seqid="chr1", strand=strand, cds=sequence, blocks=blocks, origin="original")
    candidate["quality"] = validate_candidate(candidate, code)
    assert candidate["quality"]["usable"]
    return candidate


@pytest.mark.parametrize("strand,blocks,expected", [
    ("+", [[0, 1, 2], [10, 20, 1], [30, 31, 0]], [[11, 20, 0]]),
    ("-", [[40, 41, 2], [20, 30, 1], [0, 1, 0]], [[20, 29, 0]]),
])
def test_analysis_coding_view_clips_partial_bases_across_short_blocks_without_source_edits(strand, blocks, expected):
    candidate = analysis_source("TTATGAAATAAC", blocks, strand=strand)
    original = copy.deepcopy(candidate)
    analysis = refinement.analysis_coding_candidate(candidate)
    assert candidate == original
    assert analysis["cds"] == "ATGAAATAA"
    assert analysis["protein"] == candidate["protein"] == "MK"
    assert analysis["blocks"] == expected
    assert analysis["quality"]["translation_offset"] == 0
    assert analysis["analysis"] == dict(removed_first_bases=2, removed_final_bases=1,
                                        source_cds_length=12, analysis_cds_length=9)


def test_analysis_coding_view_recomputes_split_codon_phases_and_respects_genetic_code():
    candidate = analysis_source("TATGAAATAA", [[0, 5, 1], [10, 15, 2]])
    analysis = refinement.analysis_coding_candidate(candidate)
    assert analysis["blocks"] == [[1, 5, 0], [10, 15, 2]]
    assert analysis["junctions"] == [[5, 10, 1]]
    alternative_code = analysis_source("ATGTGATAA", [[0, 9, 0]], code=4)
    assert refinement.analysis_coding_candidate(alternative_code)["protein"] == "MW"
    invalid = copy.deepcopy(candidate)
    invalid["protein"] = "different"
    with pytest.raises(ValueError, match="reproduce the admitted protein"):
        refinement.analysis_coding_candidate(invalid)
    invalid["quality"]["usable"] = False
    with pytest.raises(ValueError, match="admitted source translation"):
        refinement.analysis_coding_candidate(invalid)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_analysis_gff_reconstructs_same_admitted_protein_and_keeps_source_metadata(tmp_path, strand):
    from Bio.Seq import Seq
    owner = "t,%raw"
    block = [20, 30, 1] if strand == "-" else [0, 10, 1]
    candidate = analysis_source("TATGAAATAA", [block], owner=owner, strand=strand)
    analysis = refinement.analysis_coding_candidate(candidate)
    start, end = block[0] + 1, block[1]
    selected = (f"##gff-version 3\nchr1\ts\tgene\t{start}\t{end}\t5.5\t{strand}\t.\tID=g\n"
                f"chr1\ts\tmRNA\t{start}\t{end}\t7.5\t{strand}\t.\tID=t%2C%25raw;Parent=g;Note=source%3Devidence\n"
                f"chr1\ts\texon\t{start}\t{end}\t9.5\t{strand}\t.\tID=e;Parent=t%2C%25raw;custom=kept\n"
                f"chr1\ts\tCDS\t{start}\t{end}\t.\t{strand}\t1\tID=c;Parent=t%2C%25raw;Target=q 1 2 +;Gap=M2\n"
                "chr1\ts\tgene\t40\t50\t.\t+\t.\tID=bad\n"
                "chr1\ts\tmRNA\t40\t50\t.\t+\t.\tID=not_admitted;Parent=bad\n"
                "chr1\ts\tCDS\t40\t50\t.\t+\t0\tParent=not_admitted\n")
    derived = refinement.analysis_gff_rows(selected, [analysis])
    assert "not_admitted" not in derived
    assert "Target=q 1 2 +;Gap=M2" in derived
    assert "Note=source%3Devidence" in derived
    assert f"\texon\t{start}\t{end}\t9.5\t{strand}\t" in derived
    genome = tmp_path / "genome.fa"
    genomic = "N" * 20 + str(Seq(candidate["cds"]).reverse_complement()) + "N" * 30 if strand == "-" else candidate["cds"] + "N" * 50
    genome.write_text(">chr1\n" + genomic + "\n")
    fasta, gff = tmp_path / "analysis.fa", tmp_path / "analysis.gff3"
    fasta.write_text(">" + owner + "\n" + analysis["cds"] + "\n")
    gff.write_text(derived)
    result = refinement.build_catalog("Species", fasta, gff, genome)
    reconstructed = result["loci"][0]["candidates"][0]
    assert reconstructed["cds"] == analysis["cds"]
    assert reconstructed["protein"] == analysis["protein"]
    assert reconstructed["blocks"] == analysis["blocks"]


def test_analysis_direct_gene_cds_keeps_existing_identity_without_synthetic_transcript_collision():
    candidate = analysis_source("TATGAAATAA", [[0, 10, 1]], owner="g")
    selected = ("##gff-version 3\nchr1\ts\tgene\t1\t10\t.\t+\t.\tID=g\n"
                "chr1\ts\tCDS\t1\t10\t.\t+\t1\tID=c;Parent=g\n")
    derived = refinement.analysis_gff_rows(selected, [refinement.analysis_coding_candidate(candidate)])
    assert "\tmRNA\t" not in derived
    assert "\tCDS\t2\t10\t.\t+\t0\tID=c;Parent=g" in derived
    assert derived.count("ID=g\n") == 1


def test_analysis_shared_cds_splits_per_owner_phase_and_avoids_active_gene_id_collisions():
    a = analysis_source("ATGAATGAAATAAC", [[0, 4, 0], [10, 20, 2]], owner="ta,1")
    b = analysis_source("ATGAAATAAC", [[10, 20, 2]], owner="tb%raw")
    candidates = [refinement.analysis_coding_candidate(a), refinement.analysis_coding_candidate(b)]
    selected = ("##gff-version 3\nchr1\ts\tgene\t1\t20\t.\t+\t.\tID=g1\n"
                "chr1\ts\tgene\t11\t20\t.\t+\t.\tID=g2\n"
                "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=ta%2C1;Parent=g1\n"
                "chr1\ts\tmRNA\t11\t20\t.\t+\t.\tID=tb%25raw;Parent=g2\n"
                "chr1\ts\tCDS\t1\t4\t.\t+\t0\tID=up;Parent=ta%2C1\n"
                "chr1\ts\tCDS\t11\t20\t.\t+\t2\tID=shared;Parent=ta%2C1,tb%25raw;custom=source\n"
                "chr1\ts\tsequence_feature\t13\t18\t.\t+\t.\tID=child;Parent=shared\n")
    derived = refinement.analysis_gff_rows(selected, candidates)
    records = [(fields, refinement.parse_gff_attributes(fields[8])) for line in derived.splitlines()
               if not line.startswith("#") and len(fields := line.split("\t")) == 9]
    shared = [(fields, attrs) for fields, attrs in records if fields[2] == "CDS" and attrs.get("custom") == ("source",)]
    assert len(shared) == 2
    assert {(fields[3], fields[4], fields[7], attrs["Parent"]) for fields, attrs in shared} == {
        ("11", "18", "2", ("ta,1",)), ("13", "18", "0", ("tb%raw",))}
    identifiers = {attrs["ID"][0] for _, attrs in records if attrs.get("ID")}
    assert all(set(attrs.get("Parent", ())) <= identifiers for _, attrs in records)
    assert refinement.analysis_gff_rows(selected, list(reversed(candidates))) == derived
    collision = shared[0][1]["ID"][0]
    conflicting = selected.replace("ID=g2", "ID=" + collision).replace("Parent=g2", "Parent=" + collision)
    alternate = refinement.analysis_gff_rows(conflicting, candidates)
    records = [(fields, refinement.parse_gff_attributes(fields[8])) for line in alternate.splitlines()
               if not line.startswith("#") and len(fields := line.split("\t")) == 9]
    genes = {attrs["ID"][0] for fields, attrs in records if fields[2] == "gene"}
    cds = {attrs["ID"][0] for fields, attrs in records if fields[2] == "CDS" and attrs.get("ID")}
    assert not genes & cds


def test_two_donors_can_repair_partial_original_without_overwriting_it():
    catalog, models, edges = classification_fixture()
    original = copy.deepcopy(catalog)
    rows = classify(catalog, models, edges)
    assert len(rows) == 1
    assert rows[0]["status"] == "accepted"
    assert rows[0]["change_type"] == "model_revision"
    assert rows[0]["evidence_class"] == "homology_only_predicted"
    assert rows[0]["candidate"]["quality"]["representative_eligible"]
    assert catalog == original
    assert rows[0]["candidate"]["quality"]["valid_orf"]
    assert validate_candidate(rows[0]["candidate"])["usable"]


def test_many_isoforms_of_one_donor_are_not_independent_species_support():
    catalog, models, edges = classification_fixture()
    models[1]["donor_species"] = models[0]["donor_species"]
    rows = classify(catalog, models, edges)
    assert rows[0]["status"] == "proposal"
    assert "insufficient_independent_support" in rows[0]["problems"]
    assert not rows[0]["candidate"]["quality"]["representative_eligible"]


def test_prediction_correspondence_is_direction_independent_and_locus_specific():
    catalog, models, edges = classification_fixture()
    reverse = [{**e, 'species_a': e['species_b'], 'gene_a': e['gene_b'],
                'species_b': e['species_a'], 'gene_b': e['gene_a']} for e in edges]
    unrelated = {**edges[0], 'species_a': 'Other_target', 'ambiguous': True}
    accepted = classify(catalog, models, [unrelated, *reverse])
    assert accepted == classify(catalog, models, edges)
    ambiguous = {**reverse[0], 'ambiguous': True}
    row = classify(catalog, models, [*reverse, ambiguous])[0]
    assert row['status'] == 'proposal'
    assert 'ambiguous_locus_correspondence' in row['problems']
    untrusted = copy.deepcopy(models)
    untrusted[0]['donor_species'] = 'Unrelated_donor'
    assert 'untrusted_donor_correspondence' in classify(catalog, untrusted, reverse)[0]['problems']


@pytest.mark.parametrize("rna_supported", [False, True])
def test_source_sequence_contradiction_cannot_be_treated_as_homology_only_incompleteness(rna_supported):
    catalog, models, edges = classification_fixture()
    original = catalog["loci"][0]["candidates"][0]
    original["quality"].update(sequence_mismatch=True, partial=False, has_start=True, has_stop=True,
                               usable=False, valid_orf=False)
    before = copy.deepcopy(catalog)
    rna = [dict(species="Species_target", seqid="chr1", strand="+", cds_blocks="[[100,160]]",
                transcript_id="independent-complete-path", count="1")] if rna_supported else []
    row = classify(catalog, models, edges, rna)[0]
    assert row["status"] == ("accepted" if rna_supported else "proposal")
    assert row["candidate"]["quality"]["representative_eligible"] is rna_supported
    assert ("unresolved_source_sequence_mismatch" in row["problems"]) is not rna_supported
    assert catalog == before


def test_homology_addition_to_intact_gene_requires_target_path_for_representative_adoption():
    catalog, models, edges = classification_fixture(valid_original=True)
    row = classify(catalog, models, edges)[0]
    assert row["status"] == "accepted"
    assert row["change_type"] == "isoform_addition"
    assert not row["candidate"]["quality"]["representative_eligible"]
    rna = [dict(species="Species_target", seqid="chr1", strand="+", cds_blocks="[[100,160]]",
                transcript_id="rna-full-path", count="3")]
    row = classify(catalog, models[:1], edges, rna)[0]
    assert row["status"] == "accepted"
    assert row["evidence_class"] == "rna_path_supported"
    assert row["candidate"]["quality"]["representative_eligible"]
    # Individual junction/partial exon evidence is not a coexisting coding path.
    rna[0]["cds_blocks"] = "[[100,130]]"
    row = classify(catalog, models[:1], edges, rna)[0]
    assert row["status"] == "proposal"
    assert row["evidence_class"] == "homology_only_predicted"


def test_conservation_supported_adoption_changes_only_the_independent_gate():
    catalog, models, edges = classification_fixture(valid_original=True)
    baseline = classify(catalog, models, edges)[0]
    relaxed = classify(catalog, models, edges, isoform_adoption="conservation_supported")[0]
    expected = copy.deepcopy(baseline)
    expected["candidate"]["quality"].update(representative_eligible=True,
                                           representative_admission="conservation_supported")
    assert relaxed == expected
    assert relaxed["evidence_class"] == "homology_only_predicted"
    assert not relaxed["candidate"]["quality"]["rna_supported"]


@pytest.mark.parametrize("problem", ["frameshift", "internal_stop", "invalid_phase"])
def test_relaxing_rna_adoption_does_not_admit_failed_predictions(problem):
    catalog, models, edges = classification_fixture(valid_original=True)
    for model in models:
        model["problems"] = [problem]
    row = classify(catalog, models, edges, isoform_adoption="conservation_supported")[0]
    assert row["status"] == "proposal"
    assert not row["candidate"]["quality"]["representative_eligible"]


@pytest.mark.parametrize("policy,adoption", [("conserved", "unknown"), ("longest", "conservation_supported")])
def test_invalid_isoform_adoption_configuration_is_rejected(tmp_path, policy, adoption):
    inputs, edges, _ = tiny_inputs(tmp_path)
    with pytest.raises(ValueError, match="(?i)(adoption|conservation)"):
        refinement.plan(tmp_path / "run", inputs=inputs, edges=edges, mode="off",
                        policy=policy, isoform_adoption=adoption)


@pytest.mark.parametrize("case,problem", [
    ("other_locus", "overlap_other_locus"),
    ("unowned", "no_overlap_owned_locus"),
    ("wrong_strand", "wrong_locus_strand"),
    ("existing_path", "existing_coding_path"),
    ("pseudogene", "protected_annotation_exception"),
    ("translation_exception", "protected_annotation_exception"),
    ("sequence_exception", "protected_annotation_exception"),
])
def test_prediction_ownership_and_biological_exceptions_are_rejection_gates(case, problem):
    catalog, models, edges = classification_fixture(other_locus=case == "other_locus")
    if case == "unowned":
        for model in models:
            model["cds"] = [[300, 360, 0]]
    elif case == "wrong_strand":
        for model in models:
            model["strand"] = "-"
    elif case == "existing_path":
        catalog["loci"][0]["candidates"][0]["blocks"] = [[100, 160, 0]]
    elif case in {"pseudogene", "translation_exception", "sequence_exception"}:
        flag = "annotated_pseudogene" if case == "pseudogene" else case
        catalog["loci"][0]["candidates"][0]["quality"][flag] = True
    row = classify(catalog, models, edges)[0]
    assert row["status"] == "proposal"
    assert problem in row["problems"]
    assert not row["candidate"]["quality"]["representative_eligible"]


def test_prediction_identity_is_stable_across_order_and_donor_count():
    catalog, models, edges = classification_fixture()
    identifier = classify(catalog, models, edges)[0]["candidate"]["candidate_id"]
    assert classify(catalog, list(reversed(models)), edges)[0]["candidate"]["candidate_id"] == identifier
    extra = {**models[0], "donor_species": "Species_donor3"}
    assert classify(catalog, models + [extra], edges)[0]["candidate"]["candidate_id"] == identifier
    changed = copy.deepcopy(models)
    for model in changed:
        model["sequence"] = "ATG" + "CCC" * 18 + "TAA"
    assert classify(catalog, changed, edges)[0]["candidate"]["candidate_id"] != identifier


def test_noncoding_annotation_ownership_rejects_neighbor_transfer(tmp_path):
    catalog, models, edges = classification_fixture()
    source = tmp_path / "source.gff3"
    source.write_text("##gff-version 3\n"
                      "chr1\tprovider\tgene\t101\t130\t.\t+\t.\tID=g\n"
                      "chr1\tprovider\tmRNA\t101\t130\t.\t+\t.\tID=t;Parent=g\n"
                      "chr1\tprovider\tCDS\t101\t130\t.\t+\t0\tParent=t\n"
                      "chr1\tprovider\tgene\t151\t180\t.\t+\t.\tID=noncoding_gene\n"
                      "chr1\tprovider\tncRNA\t151\t180\t.\t+\t.\tID=noncoding_transcript;Parent=noncoding_gene\n")
    catalog["annotation_spans"] = refinement.annotation_ownership_spans(source, catalog)
    assert any(row["gene_id"].startswith("annotation:") for row in catalog["annotation_spans"])
    row = classify(catalog, models, edges)[0]
    assert row["status"] == "proposal"
    assert "overlap_other_locus" in row["problems"]
    assert not row["candidate"]["quality"]["representative_eligible"]


def test_exon_only_gtf_noncoding_owner_protects_the_entire_declared_locus(tmp_path):
    catalog, models, edges = classification_fixture()
    source = tmp_path / "implicit.gtf"
    source.write_text('chr1\tprovider\texon\t141\t150\t.\t+\t.\tgene_id "nc"; transcript_id "nct1";\n'
                      'chr1\tprovider\texon\t171\t180\t.\t+\t.\tgene_id "nc"; transcript_id "nct1";\n'
                      'chr1\tprovider\texon\t151\t160\t.\t+\t.\tgene_id "nc"; transcript_id "nct2";\n')
    catalog["annotation_spans"] = refinement.annotation_ownership_spans(source, catalog)
    assert catalog["annotation_spans"] == [{"seqid": "chr1", "start": 140, "end": 180, "gene_id": "annotation:nc"}]
    row = classify(catalog, models, edges)[0]
    assert row["status"] == "proposal"
    assert "overlap_other_locus" in row["problems"]
    assert not row["candidate"]["quality"]["representative_eligible"]


def test_competing_predictions_for_separate_loci_cannot_both_claim_new_interval():
    catalog, models, edges = classification_fixture()
    second = copy.deepcopy(catalog["loci"][0])
    second["gene_id"] = "Species_target_second"
    second["candidates"][0].update(candidate_id="Species_target_second_t", source_transcript_id="second_t", source_gene_id="second",
                                   blocks=[[180, 210, 0]])
    catalog["loci"].append(second)
    alternatives = [{**model, "gene_id": second["gene_id"], "cds": [[150, 210, 0]]} for model in models]
    edges += [{**edge, "gene_a": second["gene_id"], "gene_b": edge["gene_b"] + "_second"} for edge in edges]
    result = classify(catalog, models + alternatives, edges)
    assert len(result) == 2
    assert all(row["status"] == "proposal" and "competing_predictions_overlap" in row["problems"] for row in result)
    assert all(not row["candidate"]["quality"]["representative_eligible"] for row in result)


def test_off_mode_stages_restart_and_effective_catalog_agrees_with_source_path(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    sources = {row[role]: Path(row[role]).read_bytes() for row in rows for role in ("cds", "gff", "genome")}
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    assert refinement.plan(root, inputs=inputs, edges=edges, mode="off") == value
    effective = refinement.finalize(root, value)
    receipts = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in root.rglob("receipt.json")}
    assert refinement.finalize(root, value) == effective
    assert receipts == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in root.rglob("receipt.json")}
    manifest = refinement.verify_inputs(effective / "inputs.tsv")
    assert len(manifest) == 3
    for directory in ("catalog_index", "catalog_index_final"):
        summary = json.loads((root / directory / "summary.json").read_text())
        assert Path(summary["database"]) == root / directory / "loci.sqlite3"
        assert Path(summary["database"]).is_file()
        assert summary["candidate_count"] == 6
    for row in manifest:
        selected = refinement.build_catalog(row["species"], row["cds"], row["gff"], row["genome"])
        candidate = selected["loci"][0]["candidates"][0]
        assert candidate["source_transcript_id"] == "t1"
        assert candidate["cds"] == "ATGAAACCCTAA"
        full = (effective / "full_annotation" / (row["species"] + ".gff3")).read_text()
        assert "ID=t1;Parent=g" in full and "ID=t2;Parent=g" in full
    assert all(Path(path).read_bytes() == content for path, content in sources.items())
    expected = [row["source_transcript_id"] for row in refinement.read_table(effective / "representative_map.tsv")]
    assert expected == ["t1"] * 3


def test_frozen_plan_rejects_parameter_or_input_mutation(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    with pytest.raises(ValueError, match="Frozen refinement plan differs"):
        refinement.plan(root, inputs=inputs, edges=edges, mode="off", min_margin=.5)
    original = Path(rows[0]["gff"])
    original.write_text(original.read_text() + "# changed annotation metadata\n")
    with pytest.raises(ValueError, match="Frozen refinement input changed"):
        refinement.load(root)


def test_effective_input_tampering_fails_verification(tmp_path):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    manifest = effective / "inputs.tsv"
    row = refinement.read_table(manifest)[0]
    cds = Path(row["cds"])
    cds.write_text(cds.read_text().replace("AAA", "CCC"))
    with pytest.raises(ValueError, match="Effective input changed"):
        refinement.verify_inputs(manifest)


def test_effective_verification_freshly_reads_every_unique_publication_file_once(tmp_path, monkeypatch):
    from collections import Counter

    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    effective = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off'))
    manifest = effective / 'inputs.tsv'
    receipt = json.loads((effective / 'receipt.json').read_text())
    state = import_module('input_generation_array_state')
    original, reads = state.digest, Counter()

    def counted(path):
        reads[Path(path).resolve()] += 1
        return original(path)

    monkeypatch.setattr(state, 'digest', counted)
    rows = refinement.verify_inputs(manifest)
    expected = {effective / 'receipt.json', *(effective / path for path in receipt['files'])}
    assert reads == Counter({path.resolve(): 1 for path in expected})
    # A fresh invocation must still read all bytes and catch a changed source.
    genome = Path(rows[0]['genome'])
    genome.write_text(genome.read_text() + '\n')
    with pytest.raises(ValueError, match='Effective input changed'):
        refinement.verify_inputs(manifest)
    assert reads[genome.resolve()] == 2


@pytest.mark.parametrize('target', ['genome', 'receipt', 'auxiliary'])
def test_effective_verification_rejects_changes_after_an_earlier_batch_hash(tmp_path, monkeypatch, target):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    effective = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off'))
    manifest = effective / 'inputs.tsv'
    row = refinement.read_table(manifest)[0]
    path = {'genome': Path(row['genome']), 'receipt': effective / 'receipt.json',
            'auxiliary': effective / 'species_genetic_code.tsv'}[target]
    state = import_module('input_generation_array_state')
    original = state.digest

    def changed_after_hash(source):
        result = original(source)
        if Path(source).resolve() == path.resolve():
            path.write_bytes(path.read_bytes() + b'\n')
        return result

    monkeypatch.setattr(state, 'digest', changed_after_hash)
    with pytest.raises(OSError, match='changed while hashing'):
        refinement.verify_inputs(manifest)


def test_effective_verification_binds_parsing_before_the_hash_batch(tmp_path, monkeypatch):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    effective = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off'))
    manifest, receipt_path = effective / 'inputs.tsv', effective / 'receipt.json'
    original = refinement.digest_paths

    def changed_before_batch(paths):
        rows = refinement.read_table(manifest)
        rows[0]['genetic_code'] = '4'
        write_tsv(manifest, list(rows[0]), rows)
        receipt = json.loads(receipt_path.read_text())
        receipt['files']['inputs.tsv'] = refinement.digest(manifest)
        refinement.atomic_json(receipt_path, receipt)
        return original(paths)

    monkeypatch.setattr(refinement, 'digest_paths', changed_before_batch)
    with pytest.raises(OSError, match='manifest or receipt changed while verifying'):
        refinement.verify_inputs(manifest)


def test_effective_verification_retains_external_role_paths_and_checks_their_bytes(tmp_path):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    effective = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off'))
    manifest, receipt_path = effective / 'inputs.tsv', effective / 'receipt.json'
    rows = refinement.read_table(manifest)
    external = tmp_path / 'external-genome.fa'
    shutil.copyfile(rows[0]['genome'], external)
    rows[0]['genome'] = str(external)
    write_tsv(manifest, list(rows[0]), rows)
    receipt = json.loads(receipt_path.read_text())
    receipt['files']['inputs.tsv'] = refinement.digest(manifest)
    refinement.atomic_json(receipt_path, receipt)
    assert refinement.verify_inputs(manifest)[0]['genome'] == str(external)
    external.write_text(external.read_text() + '\n')
    with pytest.raises(ValueError, match='Effective input changed'):
        refinement.verify_inputs(manifest)


def test_effective_verification_checks_role_hashes_and_unselected_receipt_members(tmp_path):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    effective = refinement.finalize(root, refinement.plan(root, inputs=inputs, edges=edges, mode='off'))
    manifest, receipt_path = effective / 'inputs.tsv', effective / 'receipt.json'
    rows = refinement.read_table(manifest)
    rows[0]['genome_sha256'] = '0' * 64
    write_tsv(manifest, list(rows[0]), rows)
    receipt = json.loads(receipt_path.read_text())
    receipt['files']['inputs.tsv'] = refinement.digest(manifest)
    refinement.atomic_json(receipt_path, receipt)
    with pytest.raises(ValueError, match='Effective input changed'):
        refinement.verify_inputs(manifest)
    rows[0]['genome_sha256'] = refinement.digest(rows[0]['genome'])
    write_tsv(manifest, list(rows[0]), rows)
    receipt['files']['inputs.tsv'] = refinement.digest(manifest)
    refinement.atomic_json(receipt_path, receipt)
    auxiliary = effective / 'species_genetic_code.tsv'
    auxiliary.write_bytes(auxiliary.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='Effective view is incomplete or corrupt'):
        refinement.verify_inputs(manifest)


def test_correspondence_rejects_absent_locus(tmp_path):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    links = refinement.read_table(edges)
    links[0]["gene_b"] = "Species_donor1_absent"
    write_tsv(edges, list(links[0]), links)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    with pytest.raises(ValueError, match="absent locus"):
        refinement.correspondence(root, value)


def test_one_to_many_correspondence_abstains_for_ambiguous_duplicate_copy(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    donor = rows[1]
    genome, gff, cds = [Path(donor[role]) for role in ("genome", "gff", "cds")]
    genome.write_text(genome.read_text() + ">chr2\nATGAAATTTTAA\n")
    gff.write_text(gff.read_text() + "chr2\ts\tgene\t1\t12\t.\t+\t.\tID=g2\n"
                   "chr2\ts\tmRNA\t1\t12\t.\t+\t.\tID=t3;Parent=g2\n"
                   "chr2\ts\tCDS\t1\t12\t.\t+\t0\tID=c4;Parent=t3\n")
    cds.write_text(cds.read_text() + ">g2\nATGAAATTTTAA\n")
    links = refinement.read_table(edges)
    links.append({**links[0], "gene_b": donor["species"] + "_g2"})
    write_tsv(edges, list(links[0]), links)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    output = refinement.correspondence(root, value)
    actual = json.loads((output / "edges.json").read_text())
    duplicated = [edge for edge in actual if {edge["species_a"], edge["species_b"]} == {rows[0]["species"], donor["species"]}]
    assert len(duplicated) == 2 and all(edge["ambiguous"] for edge in duplicated)
    choices = json.loads((refinement.select(root, value) / "selection.json").read_text())["selections"]
    target = next(row for row in choices if row["species"] == rows[0]["species"])
    assert target["candidate_id"] == "Species_target_t1"
    assert target["status"] != "selected"


def test_qc_fails_on_corrupt_catalog_receipt_artifact(tmp_path):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    artifact = root / "catalog" / "Species_target" / "catalog.json"
    artifact.write_text(artifact.read_text() + "\n")
    result = subprocess.run([sys.executable, refinement.__file__, "qc", "--output", str(root)],
                            text=True, capture_output=True)
    assert result.returncode != 0
    assert "Corrupt refinement stage" in result.stderr


@pytest.mark.parametrize("damage", ["missing_catalog", "missing_index", "foreign_key", "upstream_key"])
def test_qc_requires_every_stage_and_its_current_dependency_key(tmp_path, damage):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    assert all(refinement.inspect_stages(root, value).values())
    if damage.startswith("missing_"):
        directory = root / ("catalog/Species_target" if damage == "missing_catalog" else "catalog_index")
        shutil.rmtree(directory)
    else:
        path = root / "selection_final/receipt.json"
        receipt = json.loads(path.read_text())
        if damage == "foreign_key":
            receipt["key"]["plan"] = "other-run"
        else:
            receipt["key"]["dependencies"]["correspondence"] = "obsolete-stage"
        path.write_text(json.dumps(receipt))
    result = subprocess.run([sys.executable, refinement.__file__, "qc", "--output", str(root)],
                            text=True, capture_output=True)
    assert result.returncode != 0
    assert "Corrupt refinement stage" in result.stderr


@pytest.mark.parametrize("damage", ["genetic_code", "species", "missing_row"])
def test_copied_effective_manifest_cannot_change_unhashed_metadata(tmp_path, damage):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    published = effective / "inputs.tsv"
    copied = tmp_path / "copied.tsv"
    shutil.copyfile(published, copied)
    assert refinement.verify_inputs(copied) == refinement.verify_inputs(published)
    rows = refinement.read_table(copied)
    if damage == "missing_row":
        rows.pop()
    else:
        rows[0][damage] = "4" if damage == "genetic_code" else "Foreign_species"
    write_tsv(copied, list(rows[0]), rows)
    with pytest.raises(ValueError, match="manifest differs from the published view"):
        refinement.verify_inputs(copied)


def test_inconsistent_supplied_cds_is_archived_and_excluded_without_genomic_substitution(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    cds = Path(target["cds"])
    supplied = b">t1 conflicting supplied genotype\naTgCcCcCcTaA\n"
    cds.write_bytes(supplied)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, policy="longest", mode="off")
    effective = refinement.finalize(root, value)
    archived = effective / "source_cds" / (target["species"] + ".fa")
    assert archived.read_bytes() == supplied
    assert hashlib.sha256(archived.read_bytes()).hexdigest() == hashlib.sha256(cds.read_bytes()).hexdigest()
    assert (effective / "source_annotation" / (target["species"] + ".gff3")).read_bytes() == Path(target["gff"]).read_bytes()
    for role in ("species_cds", "species_protein"):
        assert (effective / role / (target["species"] + ".fa")).read_text() == ""
    assert "ID=t1" not in (effective / "species_gff" / (target["species"] + ".gff3")).read_text()
    assert "ID=t1;Parent=g" in (effective / "full_annotation" / (target["species"] + ".gff3")).read_text()
    excluded = refinement.read_table(effective / "effective_exclusions.tsv")
    assert [(r["gene_id"], r["reason"]) for r in excluded] == [(target["species"] + "_g", "source_cds_sequence_mismatch")]
    catalog = json.loads((root / "catalog" / target["species"] / "catalog.json").read_text())
    candidate = next(c for g in catalog["loci"] for c in g["candidates"] if c["source_transcript_id"] == "t1")
    assert candidate["cds"] == "ATGAAACCCTAA"
    assert candidate["source_cds"][0]["cds"] == "aTgCcCcCcTaA"
    assert not candidate["source_cds"][0]["sequence_agreement"]
    assert not candidate["quality"]["usable"]


def test_pseudogenic_transcript_type_preserves_cds_but_cannot_admit_protein(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    source.write_text(source.read_text().replace("\tmRNA\t1\t12\t.\t+\t.\tID=t1", "\tpseudogenic_transcript\t1\t12\t.\t+\t.\tID=t1"))
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, policy="longest", mode="off")
    effective = refinement.finalize(root, value)
    selected = refinement.build_catalog(target["species"], effective / "species_cds" / (target["species"] + ".fa"),
                                        effective / "species_gff" / (target["species"] + ".gff3"), target["genome"])
    candidate = selected["loci"][0]["candidates"][0]
    assert candidate["source_transcript_id"] == "t1"
    assert candidate["cds"] == "ATGAAACCCTAA"
    assert candidate["quality"]["annotated_pseudogene"]
    assert not candidate["quality"]["usable"]
    assert (effective / "species_protein" / (target["species"] + ".fa")).read_text() == ""
    assert "pseudogenic_transcript" in (effective / "full_annotation" / (target["species"] + ".gff3")).read_text()


def test_ancestor_axis_problem_is_excluded_consistently_from_effective_views(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    original = source.read_text().replace("chr1\ts\tgene\t1\t12\t.\t+", "chr2\ts\tgene\t1\t12\t.\t-")
    source.write_text(original)
    genome = Path(target["genome"])
    genome.write_text(genome.read_text() + ">chr2\nATGAAACCCTAA\n")
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    for role in ("species_cds", "species_protein"):
        assert (effective / role / (target["species"] + ".fa")).read_text() == ""
    assert "\tCDS\t" not in (effective / "species_gff" / (target["species"] + ".gff3")).read_text()
    exclusions = refinement.read_table(effective / "effective_exclusions.tsv")
    assert [(row["species"], row["reason"]) for row in exclusions] == [
        (target["species"], "invalid_annotation_structure:ancestor_contig_or_strand_mismatch")]
    assert (effective / "source_annotation" / (target["species"] + ".gff3")).read_bytes() == source.read_bytes()
    assert "chr2\ts\tgene\t1\t12\t.\t-" in (effective / "full_annotation" / (target["species"] + ".gff3")).read_text()


@pytest.mark.parametrize("attributes,identifier", [
    ('gene_id "g"; gene_name "value=retained";', "g"),
    ('gene_id "g%2C1";', "g%2C1"),
    ('ID=g%2C1;Name=value%3Dretained', "g,1"),
])
def test_gene_bound_expansion_uses_typed_gtf_and_gff_identities(attributes, identifier):
    original = f"chr1\ts\tgene\t101\t130\t.\t+\t.\t{attributes}\n"
    derived, changes = refinement.extend_gene_bounds(original, {identifier: (100, 160)})
    assert derived.split("\t")[3:5] == ["101", "160"]
    assert derived.split("\t")[8].strip() == attributes
    assert changes == [{"source_gene_id": identifier, "before": [101, 130], "after": [101, 160]}]


def test_full_gtf_view_preserves_noncoding_metadata_and_synthesizes_missing_graph_nodes():
    original = ('chr1\tprovider\tCDS\t1\t12\t.\t+\t0\tgene_id "g"; transcript_id "t"; custom "a=b;kept";\n'
                'chr1\tprovider\texon\t31\t42\t.\t-\t.\tgene_id "nc"; transcript_id "nct"; note "noncoding=kept";\n'
                'chr1\tprovider\tgene\t61\t72\t5.5\t+\t.\tgene_id "other"; gene_name "name=kept";\n')
    derived = refinement.selected_gff_rows(original, set(), retain_all=True)
    parsed = [(fields[2], refinement.parse_gff_attributes(fields[8])) for line in derived.splitlines()
              if not line.startswith("#") and len(fields := line.split("\t")) == 9]
    assert derived.startswith("##gff-version 3\n")
    assert {attrs["ID"][0] for kind, attrs in parsed if kind == "gene"} == {"g", "nc", "other"}
    assert {attrs["ID"][0] for kind, attrs in parsed if kind == "mRNA"} == {"t"}
    assert {attrs["ID"][0] for kind, attrs in parsed if kind == "transcript"} == {"nct"}
    assert any(kind == "CDS" and attrs["custom"] == ("a=b;kept",) and attrs["Parent"] == ("t",) for kind, attrs in parsed)
    assert any(kind == "exon" and attrs["note"] == ("noncoding=kept",) and attrs["Parent"] == ("nct",) for kind, attrs in parsed)
    assert any(attrs.get("gene_name") == ("name=kept",) for _, attrs in parsed)


def test_existing_gff3_structured_alignment_attributes_keep_original_syntax():
    original = ("##gff-version 3\nchr1\ts\tgene\t1\t12\t.\t+\t.\tID=g\n"
                "chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\n"
                "chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=t;Target=query 1 4 +;Gap=M4\n")
    for retain_all in (False, True):
        derived = refinement.selected_gff_rows(original, {"t"}, retain_all=retain_all)
        assert "Target=query 1 4 +;Gap=M4" in derived


def test_missing_explicit_gff3_gene_parent_is_synthesized_in_derived_views():
    original = ("##gff-version 3\nchr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g%2C1;custom=provider%3Devidence\n"
                "chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=t\n")
    for retain_all in (False, True):
        derived = refinement.selected_gff_rows(original, {"t"}, retain_all=retain_all)
        records = [(fields, refinement.parse_gff_attributes(fields[8])) for line in derived.splitlines()
                   if not line.startswith("#") and len(fields := line.split("\t")) == 9]
        identifiers = {attrs["ID"][0] for _, attrs in records if attrs.get("ID")}
        assert identifiers == {"t", "g,1"}
        assert all(set(attrs.get("Parent", ())) <= identifiers for _, attrs in records)
        assert "custom=provider%3Devidence" in derived


@pytest.mark.parametrize("compressed", [False, True])
def test_source_annotation_archive_preserves_decoded_crlf_bytes(tmp_path, compressed):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    path = Path(target["gff"])
    content = path.read_bytes().replace(b"\n", b"\r\n")
    if compressed:
        path = path.with_suffix(path.suffix + ".gz")
        path.write_bytes(gzip.compress(content))
        target["gff"] = str(path)
        write_tsv(inputs, list(rows[0]), rows)
    else:
        path.write_bytes(content)
    before = path.read_bytes()
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    assert (effective / "source_annotation" / (target["species"] + ".gff3")).read_bytes() == content
    assert path.read_bytes() == before


@pytest.mark.parametrize("implicit", [False, True])
def test_accepted_addition_expands_canonical_full_gtf_graph_without_changing_source(tmp_path, monkeypatch, implicit):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    genomic = "ATGGTACCCTAA" + "C" * 30 + "AGAAATAA" + "C" * 40
    for row in rows:
        Path(row["genome"]).write_text(">chr1\n" + genomic + "\n")
    Path(target["cds"]).write_text(">t1\nATGGTACCCTAA\n")
    source = ('chr1\tprovider\tgene\t1\t12\t5.5\t+\t.\tgene_id "g"; gene_name "coding=name";\n'
              'chr1\tprovider\ttranscript\t1\t12\t7.5\t+\t.\tgene_id "g"; transcript_id "t1";\n'
              'chr1\tprovider\tCDS\t1\t12\t.\t+\t0\tgene_id "g"; transcript_id "t1";\n'
              'chr1\tprovider\ttranscript\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n'
              'chr1\tprovider\tCDS\t1\t6\t.\t+\t0\tgene_id "g"; transcript_id "t2";\n'
              'chr1\tprovider\tCDS\t10\t12\t.\t+\t0\tgene_id "g"; transcript_id "t2";\n'
              'chr1\tprovider\texon\t61\t72\t9.5\t+\t.\tgene_id "nc"; transcript_id "nct"; note "noncoding=kept";\n')
    if implicit:
        source = "".join(line for line in source.splitlines(keepends=True) if line.split("\t")[2] in {"CDS", "exon"})
    original = source.replace("\n", "\r\n").encode()
    Path(target["gff"]).write_bytes(original)
    for donor in rows[1:]:
        Path(donor["cds"]).write_text(">t1\nATGAAATAA\n")
        Path(donor["gff"]).write_text("##gff-version 3\nchr1\ts\tgene\t1\t50\t.\t+\t.\tID=g\n"
                                      "chr1\ts\tmRNA\t1\t50\t.\t+\t.\tID=t1;Parent=g\n"
                                      "chr1\ts\tCDS\t1\t3\t.\t+\t0\tParent=t1\n"
                                      "chr1\ts\tCDS\t45\t50\t.\t+\t0\tParent=t1\n")
    def predicted_paths(_tmp, windows, _proteins, _genome, _code, _max_intron, _cpus):
        return [dict(query=region["id"], seqid=seqid, strand="+", cds=[[0 - start, 3 - start, 0], [44 - start, 50 - start, 0]],
                     frameshift=False, coverage=1.0, identity=1.0)
                for (seqid, start, _end), regions in windows.items() for region in regions]
    monkeypatch.setattr(refinement.rescue, "search_intervals", predicted_paths)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="conservative")
    effective = refinement.finalize(root, value)
    predictions = json.loads((root / "predictions" / target["species"] / "predictions.json").read_text())
    accepted = [row for row in predictions if row["status"] == "accepted"]
    assert len(accepted) == 1
    assert not accepted[0]["candidate"]["quality"]["representative_eligible"]
    assert (effective / "source_annotation" / (target["species"] + ".gff3")).read_bytes() == original
    full = effective / "full_annotation" / (target["species"] + ".gff3")
    records = [fields for line in full.read_text().splitlines() if not line.startswith("#") and len(fields := line.split("\t")) == 9]
    parsed = [(fields, refinement.parse_gff_attributes(fields[8])) for fields in records]
    identifiers = {attrs["ID"][0] for _, attrs in parsed if attrs.get("ID")}
    assert all(set(attrs.get("Parent", ())) <= identifiers for _, attrs in parsed)
    assert next(fields[3:5] for fields, attrs in parsed if attrs.get("ID") == ("g",)) == ["1", "50"]
    assert any(attrs.get("note") == ("noncoding=kept",) for _, attrs in parsed)
    if not implicit:
        assert any(attrs.get("gene_name") == ("coding=name",) for _, attrs in parsed)
    audit = json.loads((effective / "gene_bounds_changes.json").read_text())
    assert any(row["species"] == target["species"] and row["source_gene_id"] == "g" and row["after"] == [1, 50] for row in audit)
    recatalog = refinement.build_catalog(target["species"], effective / "all_candidates" / (target["species"] + ".fa"),
                                        full, effective / "species_genome" / (target["species"] + ".fa"))
    assert recatalog["summary"]["candidates"] == 3
    assert all(candidate["quality"]["usable"] for gene in recatalog["loci"] for candidate in gene["candidates"])


@pytest.mark.parametrize("implicit", [False, True])
def test_gtf_without_id_parent_reexports_exact_selected_path_and_metadata(tmp_path, implicit):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    source.write_text('chr1\tprovider\tgene\t1\t12\t6.5\t+\t.\tgene_id "g"; gene_name "original name";\n'
                      'chr1\tprovider\ttranscript\t1\t12\t8.5\t+\t.\tgene_id "g"; transcript_id "t1"; custom "selected evidence";\n'
                      'chr1\tprovider\texon\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t1"; exon_number "1"; custom "selected exon evidence";\n'
                      'chr1\tprovider\tCDS\t1\t12\t.\t+\t0\tgene_id "g"; transcript_id "t1";\n'
                      'chr1\tprovider\ttranscript\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n'
                      'chr1\tprovider\tCDS\t1\t6\t.\t+\t0\tgene_id "g"; transcript_id "t2";\n'
                      'chr1\tprovider\tCDS\t10\t12\t.\t+\t0\tgene_id "g"; transcript_id "t2";\n')
    if implicit:
        source.write_text("".join(line for line in source.read_text().splitlines(keepends=True) if line.split("\t")[2] in {"CDS", "exon"}))
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, policy="longest", mode="off")
    effective = refinement.finalize(root, value)
    selected_gff = effective / "species_gff" / (target["species"] + ".gff3")
    selected = refinement.build_catalog(target["species"], effective / "species_cds" / (target["species"] + ".fa"), selected_gff, target["genome"])
    assert selected["summary"]["candidates"] == 1
    assert selected["loci"][0]["candidates"][0]["source_transcript_id"] == "t1"
    assert selected["loci"][0]["candidates"][0]["cds"] == "ATGAAACCCTAA"
    parsed = [(fields[2], refinement.parse_gff_attributes(fields[8])) for line in selected_gff.read_text().splitlines()
              if not line.startswith("#") and len(fields := line.split("\t")) == 9]
    assert all("t2" not in attributes.get("transcript_id", ()) for _feature, attributes in parsed)
    assert any(feature == "exon" and attributes.get("custom") == ("selected exon evidence",) for feature, attributes in parsed)
    if not implicit:
        assert any(attributes.get("custom") == ("selected evidence",) for _feature, attributes in parsed)
        assert '\tprovider\tgene\t1\t12\t6.5\t' in selected_gff.read_text()
    assert (effective / "source_annotation" / (target["species"] + ".gff3")).read_bytes() == source.read_bytes()


def test_embedded_fasta_directive_is_exact_and_comment_mentions_are_preserved(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    source.write_text(source.read_text().replace("##gff-version 3\n", "##gff-version 3\n# provider comment mentions ##FASTA without starting a sequence section\n")
                      + "##FASTA\n>chr1\nATGAAACCCTAA\n")
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    selected = (effective / "species_gff" / (target["species"] + ".gff3")).read_text()
    assert "# provider comment mentions ##FASTA" in selected
    assert "ID=t1;Parent=g" in selected
    assert "\n##FASTA\n" not in selected
    full = (effective / "full_annotation" / (target["species"] + ".gff3")).read_text()
    assert "ID=t1;Parent=g" in full and "ID=t2;Parent=g" in full
    assert full.endswith("##FASTA\n>chr1\nATGAAACCCTAA\n")


def test_named_stage_guard_also_freezes_global_rna_evidence(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    rna = tmp_path / "rna.tsv"
    evidence = dict(species=rows[0]["species"], seqid="chr1", strand="+", cds_blocks="[[0,12]]", transcript_id="rna", count=1)
    write_tsv(rna, list(evidence), [evidence])
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, rna=rna, mode="off")
    refinement.catalog_species(root, value, rows[0]["species"])
    write_tsv(rna, list(evidence), [{**evidence, "count": 2}])
    with pytest.raises(ValueError, match="Frozen refinement input changed"):
        refinement.catalog_species(root, value, rows[0]["species"])


def test_saved_plan_rejects_changed_library_identity(tmp_path, monkeypatch):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    assert value["request"]["dependencies"]["pairwise_aligner_sha256"]
    original = refinement.importlib.metadata.version
    monkeypatch.setattr(refinement.importlib.metadata, "version", lambda name: "changed-library-version" if name == "biopython" else original(name))
    with pytest.raises(ValueError, match="Refinement implementation changed"):
        refinement.load(root)


@pytest.mark.parametrize("helper", ["format_species_constants.py", "format_species_writers.py"])
def test_saved_plan_binds_shared_formatter_identity_helpers(tmp_path, monkeypatch, helper):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    assert helper in value["request"]["implementation"]
    original = refinement.digest
    monkeypatch.setattr(refinement, "digest", lambda path: "changed-helper" if Path(path).name == helper else original(path))
    with pytest.raises(ValueError, match="Refinement implementation changed"):
        refinement.load(root)


@pytest.mark.parametrize('damage', ['receipt', 'content'])
def test_selection_does_not_publish_when_dependency_changes_during_build(tmp_path, monkeypatch, damage):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    refinement.correspondence(root, value)
    original = refinement.select_from_store

    def changed_dependency(*args, **kwargs):
        result = original(*args, **kwargs)
        path = root / 'correspondence' / ('receipt.json' if damage == 'receipt' else 'edges.json')
        path.write_text(path.read_text() + '\n')
        return result

    monkeypatch.setattr(refinement, 'select_from_store', changed_dependency)
    with pytest.raises(ValueError, match='Stage dependenc'):
        refinement.select(root, value)
    assert not (root / 'selection_initial/receipt.json').exists()
    assert (root / 'selection_initial.failed/selection.json').exists()


@pytest.mark.parametrize('consumer', ['correspondence', 'selection_initial', 'selection_final', 'prediction', 'effective'])
@pytest.mark.parametrize('damage', ['receipt', 'content'])
def test_all_index_consumers_refuse_changes_during_cached_invocation(tmp_path, monkeypatch, consumer, damage):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    monkeypatch.setattr(refinement, '_INVOCATION_CACHE', {})
    final_index = consumer in {'selection_final', 'effective'}
    db = refinement.catalog_index(root, value, predictions=final_index)
    if consumer in {'selection_initial', 'selection_final'}:
        refinement.correspondence(root, value)
    elif consumer == 'prediction':
        refinement.select(root, value)
    elif consumer == 'effective':
        refinement.select(root, value, predictions=True)
    changed = False

    def change_index():
        nonlocal changed
        if changed:
            return
        if damage == 'receipt':
            receipt = db.parent / 'receipt.json'
            receipt.write_text(receipt.read_text() + '\n')
        else:
            with sqlite3.connect(db) as connection:
                gene, serialized = connection.execute('SELECT gene_id,json FROM loci WHERE species=?', ('Species_target',)).fetchone()
                locus = json.loads(serialized)
                for candidate in locus['candidates']:
                    candidate['quality']['audit_marker'] = 'changed_after_dependency_verification'
                connection.execute('UPDATE loci SET json=? WHERE species=? AND gene_id=?',
                                   (json.dumps(locus), 'Species_target', gene))
        changed = True

    if consumer == 'correspondence':
        original = refinement.iter_locus_keys

        def changing_keys(*args, **kwargs):
            yield from original(*args, **kwargs)
            change_index()

        monkeypatch.setattr(refinement, 'iter_locus_keys', changing_keys)
        def operation():
            return refinement.correspondence(root, value)
        publication = root / 'correspondence'
    elif consumer.startswith('selection_'):
        original = refinement.select_from_store

        def changing_selection(*args, **kwargs):
            result = original(*args, **kwargs)
            change_index()
            return result

        monkeypatch.setattr(refinement, 'select_from_store', changing_selection)
        def operation():
            return refinement.select(root, value, predictions=final_index)
        publication = root / consumer
    elif consumer == 'prediction':
        original = refinement.rescue.stage

        def changing_prediction(run_root, relative, key, builder, guard=None):
            def changed_builder(path):
                builder(path)
                change_index()
            return original(run_root, relative, key, changed_builder, guard)

        monkeypatch.setattr(refinement.rescue, 'stage', changing_prediction)
        def operation():
            return refinement.predict_species(root, value, 'Species_target')
        publication = root / 'predictions' / 'Species_target'
    else:
        original = refinement.iter_loci

        def changing_loci(path, *args, **kwargs):
            if Path(path) == db:
                change_index()
            yield from original(path, *args, **kwargs)

        monkeypatch.setattr(refinement, 'iter_loci', changing_loci)
        def operation():
            return refinement.finalize(root, value)
        publication = root / 'effective'
    with pytest.raises(ValueError, match='Stage dependenc'):
        operation()
    assert changed
    assert not (publication / 'receipt.json').exists()
    assert publication.with_name(publication.name + '.failed').exists()


def test_effective_binds_directly_read_catalog_metadata_content(tmp_path, monkeypatch):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    monkeypatch.setattr(refinement, '_INVOCATION_CACHE', {})
    refinement.select(root, value, predictions=True)
    original = refinement.analysis_coding_candidate
    metadata = root / 'catalog' / 'Species_target' / 'catalog_metadata.json'

    def changing_metadata(*args, **kwargs):
        result = original(*args, **kwargs)
        metadata.write_text(metadata.read_text() + '\n')
        return result

    monkeypatch.setattr(refinement, 'analysis_coding_candidate', changing_metadata)
    with pytest.raises(ValueError, match='Stage dependency content changed'):
        refinement.finalize(root, value)
    assert not (root / 'effective' / 'receipt.json').exists()


def test_correspondence_binds_prepared_id_map_content_while_reading(tmp_path, monkeypatch):
    inputs, edges, _rows = tiny_inputs(tmp_path)
    root, anchor_root = tmp_path / 'run', tmp_path / 'anchors'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    anchor_root.mkdir()
    (anchor_root / 'plan.json').write_text('{}\n')
    value['request'].update(rescue_output=str(anchor_root), edges=None)
    refinement.atomic_json(root / 'plan.json', value)
    for name in value['species']:
        directory = anchor_root / 'prepared' / name
        directory.mkdir(parents=True)
        mapping = directory / 'genes.id_map.tsv'
        write_tsv(mapping, ['jcvi_id', 'locus_id', 'original_id', 'status'],
                  [dict(jcvi_id='g', locus_id='g', original_id='g', status='selected')])
        refinement.atomic_json(directory / 'receipt.json', {'key': {'plan': refinement.digest(anchor_root / 'plan.json'), 'species': name},
                                                           'files': {'genes.id_map.tsv': refinement.digest(mapping)}})
    monkeypatch.setattr(refinement.rescue, 'load', lambda *_args: {'synteny_jobs': []})
    monkeypatch.setattr(refinement.rescue, 'prepared', lambda _root, _value, name: anchor_root / 'prepared' / name)
    original = refinement.read_table
    damaged = anchor_root / 'prepared' / value['species'][-1] / 'genes.id_map.tsv'

    def changing_id_map(path):
        result = original(path)
        if Path(path) == damaged:
            damaged.write_text(damaged.read_text() + '\n')
        return result

    monkeypatch.setattr(refinement, 'read_table', changing_id_map)
    with pytest.raises(ValueError, match='Stage dependency content changed'):
        refinement.correspondence(root, value)
    assert not (root / 'correspondence' / 'receipt.json').exists()


def test_analysis_view_clips_only_known_partial_bases_and_retains_raw_cds(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    Path(target['cds']).write_text('>g\nTATGAAATAA\n')
    Path(target['genome']).write_text('>chr1\nTATGAAATAA\n')
    Path(target['gff']).write_text('##gff-version 3\n'
                                  'chr1\ts\tgene\t1\t10\t.\t+\t.\tID=g\n'
                                  'chr1\ts\tmRNA\t1\t10\t.\t+\t.\tID=t1;Parent=g\n'
                                  'chr1\ts\texon\t1\t10\t.\t+\t.\tParent=t1;Note=original\n'
                                  'chr1\ts\tCDS\t1\t10\t.\t+\t1\tParent=t1\n')
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off', policy='longest')
    effective = refinement.finalize(root, value)
    name = target['species']
    assert (effective / 'species_cds' / (name + '.fa')).read_text() == f'>{name}_g\nTATGAAATAA\n'
    assert (effective / 'analysis_cds' / (name + '.fa')).read_text() == f'>{name}_g\nATGAAATAA\n'
    assert (effective / 'species_protein' / (name + '.fa')).read_text() == f'>{name}_g\nMK\n'
    coding = refinement.build_catalog(name, effective / 'analysis_cds' / (name + '.fa'),
                                      effective / 'analysis_gff' / (name + '.gff3'), effective / 'species_genome' / (name + '.fa'))
    candidate = coding['loci'][0]['candidates'][0]
    assert candidate['blocks'] == [[1, 10, 0]]
    assert candidate['cds'] == 'ATGAAATAA' and candidate['protein'] == 'MK'
    assert candidate['quality']['usable']
    assert 'Note=original' in (effective / 'analysis_gff' / (name + '.gff3')).read_text()
    admission = next(row for row in refinement.read_table(effective / 'coding_admission.tsv') if row['species'] == name)
    assert (admission['head_bases_removed'], admission['tail_bases_removed']) == ('1', '0')
    layout = refinement.verify_inputs(effective / 'inputs.tsv', 'coding_layout').split('\t')
    assert Path(layout[0]) == effective / 'analysis_cds'
    assert Path(layout[2]) == effective / 'analysis_gff'
    assert layout[-1] == '1'


def test_coding_analysis_binds_common_code_and_rejects_mixed_code_models(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    for row in rows:
        row['genetic_code'] = 4
    write_tsv(inputs, list(rows[0]), rows)
    root = tmp_path / 'uniform'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    effective = refinement.finalize(root, value)
    assert refinement.verify_inputs(effective / 'inputs.tsv', 'coding_layout').split('\t')[-1] == '4'
    rows[0]['genetic_code'] = 1
    write_tsv(inputs, list(rows[0]), rows)
    mixed = tmp_path / 'mixed'
    value = refinement.plan(mixed, inputs=inputs, edges=edges, mode='off')
    effective = refinement.finalize(mixed, value)
    assert len(refinement.verify_inputs(effective / 'inputs.tsv', 'layout').split('\t')) == 6
    with pytest.raises(ValueError, match='one common genetic code'):
        refinement.verify_inputs(effective / 'inputs.tsv', 'coding_layout')


@pytest.mark.parametrize("encoded,decoded", [("t%2C1", "t,1"), ("t%252C1", "t%2C1"),
                                             ("t%3B1", "t;1"), ("t%3D1", "t=1")])
def test_effective_gff_preserves_escaped_identity_through_selected_readers(tmp_path, encoded, decoded):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    source.write_text(source.read_text().replace("ID=t1;", "ID=" + encoded + ";").replace("Parent=t1", "Parent=" + encoded))
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, policy="longest", mode="off")
    effective = refinement.finalize(root, value)
    selected_path = effective / "species_gff" / (target["species"] + ".gff3")
    assert "ID=" + encoded + ";Parent=g" in selected_path.read_text()
    assert "Parent=" + encoded in selected_path.read_text()
    selected = refinement.build_catalog(target["species"], effective / "species_cds" / (target["species"] + ".fa"), selected_path, target["genome"])
    assert selected["summary"]["candidates"] == 1
    assert selected["loci"][0]["candidates"][0]["source_transcript_id"] == decoded
    assert selected["loci"][0]["candidates"][0]["cds"] == "ATGAAACCCTAA"
    reader = import_module("gff2genestat")
    columns = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    traits = reader.process_single_gff(selected_path.name, str(selected_path.parent), [target["species"] + "_g"],
                                       "CDS", "longest", columns, ["gene_id", "feature_size", "gff_transcript_id"],
                                       representative_map=effective / "representative_map.tsv")
    assert traits.iloc[0].gff_transcript_id == decoded
    assert traits.iloc[0].feature_size == 12
    resolver = import_module("cds_resolution")
    resolved, traits_path, report_path = resolver.resolve(effective / "species_cds" / (target["species"] + ".fa"),
                                                         selected_path, target["genome"], tmp_path / "resolved",
                                                         representative_map=effective / "representative_map.tsv")
    assert resolved.read_bytes() == (effective / "species_cds" / (target["species"] + ".fa")).read_bytes()
    assert refinement.read_table(traits_path)[0]["gff_transcript_id"] == decoded
    assert json.loads(report_path.read_text())["accepted_count"] == 1


def test_literal_gtf_percent_identity_survives_catalog_export_and_cds_resolution(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    target = rows[0]
    source = Path(target["gff"])
    source.write_text('chr1\tprovider\tgene\t1\t12\t.\t+\t.\tgene_id "g";\n'
                      'chr1\tprovider\ttranscript\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t%2C1";\n'
                      'chr1\tprovider\tCDS\t1\t12\t.\t+\t0\tgene_id "g"; transcript_id "t%2C1"; Note "source,metadata;literal%25";\n')
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    choices = refinement.read_table(effective / "representative_map.tsv")
    assert next(row for row in choices if row["species"] == target["species"])["source_transcript_id"] == "t%2C1"
    selected_path = effective / "species_gff" / (target["species"] + ".gff3")
    assert "ID=t%252C1" in selected_path.read_text()
    assert "Parent=t%252C1" in selected_path.read_text()
    selected = refinement.build_catalog(target["species"], effective / "species_cds" / (target["species"] + ".fa"), selected_path, target["genome"])
    assert selected["loci"][0]["candidates"][0]["source_transcript_id"] == "t%2C1"
    attrs = [refinement.parse_gff_attributes(line.split("\t")[8]) for line in selected_path.read_text().splitlines() if "\tCDS\t" in line]
    assert attrs[0]["Note"] == ("source,metadata;literal%25",)
    resolver = import_module("cds_resolution")
    output, traits, report = resolver.resolve(effective / "species_cds" / (target["species"] + ".fa"), selected_path,
                                              target["genome"], tmp_path / "resolved", representative_map=effective / "representative_map.tsv")
    assert output.read_text() == ">" + target["species"] + "_g\nATGAAACCCTAA\n"
    assert refinement.read_table(traits)[0]["gff_transcript_id"] == "t%2C1"
    assert json.loads(report.read_text())["accepted_count"] == 1


@pytest.mark.parametrize("rna_case", ["absent", "whole_path", "partial", "disconnected"])
def test_original_short_isoform_requires_whole_target_rna_for_major_span_loss(tmp_path, rna_case):
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    first, middle, last = "ATG" + "AAA" * 12, "CCC" * 18, "GGG" * 12 + "TAA"
    genomic = first + "N" * 20 + middle + "N" * 20 + last
    long_cds, short_cds = first + middle + last, first + last
    assert len(short_cds) < .70 * len(long_cds)
    all_blocks = [[0, 39], [59, 113], [133, 172]]
    rows = []
    for name in names:
        cds, gff, genome = [tmp_path / (name + suffix) for suffix in (".fa", ".gff3", ".genome.fa")]
        genome.write_text(">chr1\n" + genomic + "\n")
        annotation = "##gff-version 3\nchr1\tsource\tgene\t1\t172\t.\t+\t.\tID=g\n"
        paths = {"long": all_blocks, "short": [all_blocks[0], all_blocks[2]]} if name == names[0] else {"short": [all_blocks[0], all_blocks[2]]}
        for transcript, blocks in paths.items():
            annotation += f"chr1\tsource\tmRNA\t1\t172\t.\t+\t.\tID={transcript};Parent=g\n"
            annotation += "".join(f"chr1\tsource\tCDS\t{start + 1}\t{end}\t.\t+\t0\tParent={transcript}\n" for start, end in blocks)
        gff.write_text(annotation)
        cds.write_text(">" + ("long" if name == names[0] else "short") + "\n" + (long_cds if name == names[0] else short_cds) + "\n")
        rows.append(dict(species=name, cds=str(cds), gff=str(gff), genome=str(genome), genetic_code=1))
    inputs, edges, rna = tmp_path / "inputs.tsv", tmp_path / "edges.tsv", tmp_path / "rna.tsv"
    write_tsv(inputs, list(rows[0]), rows)
    links = [dict(species_a=names[0], gene_a=names[0] + "_g", species_b=donor, gene_b=donor + "_g") for donor in names[1:]]
    write_tsv(edges, list(links[0]), links)
    evidence_path = None
    if rna_case != "absent":
        blocks = [all_blocks[0], all_blocks[2]] if rna_case == "whole_path" else [all_blocks[0]] if rna_case == "partial" else [all_blocks[0], [130, 169]]
        evidence = dict(species=names[0], seqid="chr1", strand="+", cds_blocks=json.dumps(blocks), transcript_id="target-short-RNA", count=4)
        write_tsv(rna, list(evidence), [evidence])
        evidence_path = rna
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, rna=evidence_path, mode="off")
    output = refinement.select(root, value)
    decisions = json.loads((output / "selection.json").read_text())["selections"]
    target = next(row for row in decisions if row["species"] == names[0])
    catalog = json.loads((root / "catalog" / names[0] / "catalog.json").read_text())
    original_short = next(candidate for candidate in catalog["loci"][0]["candidates"] if candidate["source_transcript_id"] == "short")
    assert original_short["origin"] == "original"
    assert original_short["quality"].get("rna_supported", False) is (rna_case == "whole_path")
    assert original_short.get("support", {}).get("rna_paths", []) == (["target-short-RNA"] if rna_case == "whole_path" else [])
    assert catalog["loci"][0]["source_baseline_candidate_id"] == names[0] + "_long"
    if rna_case == "whole_path":
        assert target["source_transcript_id"] == "short"
        assert target["status"] == "conserved"
    else:
        assert target["source_transcript_id"] == "long"
        assert target["reason"] == "shorter_than_source_major_coding_span"


def test_fractional_rna_coordinates_cannot_claim_a_coding_path():
    catalog, models, edges = classification_fixture()
    rna = [dict(species="Species_target", seqid="chr1", strand="+", cds_blocks="[[100.5,160]]",
                transcript_id="fractional-coordinate", count="3")]
    with pytest.raises(ValueError, match="(?i)(RNA|coordinate|integer)"):
        classify(catalog, models, edges, rna)


def test_selected_gff_preserves_source_utr_exon_metadata_and_filters_shared_parents():
    original = ("##gff-version 3\n# source attribution\n"
                "chr1\tprovider\tgene\t1\t30\t7.5\t+\t.\tID=g;Name=gene-name\n"
                "chr1\tprovider\tmRNA\t1\t30\t4.5\t+\t.\tID=t1;Parent=g;Note=unselected\n"
                "chr1\tprovider\tmRNA\t1\t30\t8.5\t+\t.\tID=t2;Parent=g;Note=selected\n"
                "chr1\tprovider\tfive_prime_UTR\t1\t3\t.\t+\t.\tID=u5;Parent=t1,t2;Note=shared-UTR\n"
                "chr1\tprovider\texon\t1\t15\t.\t+\t.\tID=e1;Parent=t1,t2;custom=exon-evidence\n"
                "chr1\tprovider\tCDS\t4\t15\t.\t+\t0\tID=c1;Parent=t1,t2;protein_id=protein\n"
                "chr1\tprovider\tCDS\t22\t30\t.\t+\t0\tID=c2;Parent=t1\n"
                "chr1\tprovider\tthree_prime_UTR\t19\t30\t.\t+\t.\tID=u3;Parent=t2;Note=own-UTR\n"
                "##FASTA\n>chr1\nACGT\n")
    selected = refinement.selected_gff_rows(original, {"t2"})
    assert "ID=t1;" not in selected and "ID=c2;" not in selected
    assert "Parent=t1" not in selected
    assert "ID=u5;Parent=t2;Note=shared-UTR" in selected
    assert "ID=e1;Parent=t2;custom=exon-evidence" in selected
    assert "ID=c1;Parent=t2;protein_id=protein" in selected
    assert "ID=u3;Parent=t2;Note=own-UTR" in selected
    assert "chr1\tprovider\tgene\t1\t30\t7.5\t+\t.\tID=g;Name=gene-name\n" in selected
    assert "chr1\tprovider\tmRNA\t1\t30\t8.5\t+\t.\tID=t2;Parent=g;Note=selected\n" in selected
    assert "# source attribution" in selected
    assert "##FASTA" not in selected and ">chr1" not in selected


def test_refinement_keeps_full_organelle_annotation_but_effective_representatives_stay_nuclear(tmp_path):
    inputs, edges, rows = tiny_inputs(tmp_path)
    for row in rows:
        genome, gff = Path(row["genome"]), Path(row["gff"])
        genome.write_text(genome.read_text() + ">cp\nATGTGATAA\n")
        gff.write_text(gff.read_text() + "cp\ts\tregion\t1\t9\t.\t+\t.\tID=cp_region;genome=chloroplast\n"
                       "cp\ts\tgene\t1\t9\t.\t+\t.\tID=cp_gene\n"
                       "cp\ts\tmRNA\t1\t9\t.\t+\t.\tID=cp_transcript;Parent=cp_gene\n"
                       "cp\ts\tCDS\t1\t9\t.\t+\t0\tID=cp_cds;Parent=cp_transcript\n")
    root = tmp_path / "run"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    effective = refinement.finalize(root, value)
    for row in rows:
        full = (effective / "full_annotation" / (row["species"] + ".gff3")).read_text()
        assert "ID=cp_transcript;Parent=cp_gene" in full
        selected = (effective / "species_gff" / (row["species"] + ".gff3")).read_text()
        assert "cp_transcript" not in selected
        assert "cp_transcript" not in (effective / "representative_map.tsv").read_text()


@pytest.mark.parametrize("case,expected", [("matched", 1), ("tandem", 0), ("low_identity", 0), ("contig_break", 0)])
def test_flanking_anchors_nominate_only_unambiguous_existing_broken_locus(tmp_path, case, expected):
    catalog_module = import_module("gene_model_catalog")
    store = import_module("gene_model_store")
    directories, position_index = [], {}
    species = ["Species_target", "Species_donor"]
    for name in species:
        positions, loci = {}, []
        genes = ["left", "broken", "right"]
        if case == "tandem" and name == species[0]:
            genes.insert(2, "extra_copy")
        for number, gene in enumerate(genes):
            seqid = "chr2" if case == "contig_break" and name == species[0] and gene == "right" else "chr1"
            start = number * 1000
            gene_id = name + "_" + gene
            protein = "M" * 200 if not (case == "low_identity" and name == species[0] and gene == "broken") else "W" * 200
            positions[gene_id] = dict(gene_id=gene_id, seqid=seqid, start=start, end=start + 603)
            candidate = dict(candidate_id=gene_id + ".t", source_transcript_id=gene + ".t", source_gene_id=gene,
                             seqid=seqid, strand="+", cds="ATG" * 200 + "TAA", protein=protein,
                             blocks=[[start, start + 603, 0]], quality={"usable": True}, origin="original")
            loci.append(dict(species=name, gene_id=gene_id, seqid=seqid, strand="+", candidates=[candidate]))
        tracks = {}
        for seqid in {position["seqid"] for position in positions.values()}:
            ordered = sorted((position for position in positions.values() if position["seqid"] == seqid), key=lambda position: position["start"])
            tracks[seqid] = [position["start"] for position in ordered], ordered
        position_index[name] = positions, tracks
        directory = tmp_path / name
        catalog_module.write_catalog(dict(schema_version=1, species=name, loci=loci), directory)
        directories.append(directory)
    db = tmp_path / "loci.sqlite3"
    store.build_store(directories, db)
    anchors = [(name + "_left" for name in species), (name + "_right" for name in species)]
    anchors = [tuple(pair) for pair in anchors]
    actual = refinement.infer_flanked_loci(db, *species, anchors, position_index, refinement.DEFAULTS)
    assert len(actual) == expected
    if actual:
        assert actual[0]["gene_a"] == "Species_target_broken"
        assert actual[0]["gene_b"] == "Species_donor_broken"
        assert actual[0]["evidence"]["kind"] == "two_flanking_anchors"
        assert not actual[0]["ambiguous"]
