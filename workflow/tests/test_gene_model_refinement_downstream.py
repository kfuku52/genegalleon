"""Selected coding paths must agree across exported views and core consumers."""

import argparse
import csv
import gzip
import hashlib
import importlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote

import pytest
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
sys.path.insert(0, str(SUPPORT))

refinement = importlib.import_module("gene_model_refinement")
process_single_gff = importlib.import_module("gff2genestat").process_single_gff
prepare_genome = importlib.import_module("pairwise_synteny").prepare_genome
load_effective_inputs = importlib.import_module("representative_selection").load_effective_inputs
ensure_species_gene_cache = importlib.import_module("synteny_neighbors").ensure_species_gene_cache


def tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def selected_view(tmp_path, request):
    settings = getattr(request, 'param', {})
    short_id, long_id = settings.get('transcript_ids', ('short', 'long'))
    codons = dict(zip("ACDEFGHIKLMNPQRSTVWY", ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT",
                                             "AAA", "CTG", "ATG", "AAT", "CCT", "CAA", "CGT", "TCT",
                                             "ACT", "GTT", "TGG", "TAT"], strict=True))
    common = "M" + ("ACDEFGHIKLNPQRSTVWY" * 5)[:89]
    extension = "M" + "W" * 34
    coding = "".join(codons[residue] for residue in common) + "TAA"
    extra = "".join(codons[residue] for residue in extension)
    intron = "GTCCCCCCCCAG"
    names = ["Target_species", "Donor_one", "Donor_two"]
    rows = []
    for name in names:
        before = extra if name == names[0] else ""
        start = len(before)
        donor = start + 120
        acceptor = donor + len(intron)
        sequence = before + coding[:120] + intron + coding[120:]
        paths = {role: tmp_path / f"{name}.{role}.fa" for role in ["cds", "genome"]}
        paths["gff"] = tmp_path / f"{name}.gff3"
        paths["genome"].write_text(">chr1\n" + sequence + "\n")
        paths["cds"].write_text('>' + short_id + '\n' + coding + '\n' +
                                 ('>' + long_id + '\n' + extra + coding + '\n' if before else ''))
        lines = [f"chr1\ts\tgene\t1\t{len(sequence)}\t.\t+\t.\tID=g\n"]
        for source_transcript, first in [(short_id, start)] + ([(long_id, 0)] if before else []):
            transcript = quote(source_transcript, safe='._-:')
            lines += [f"chr1\ts\tmRNA\t{first + 1}\t{len(sequence)}\t.\t+\t.\tID={transcript};Parent=g\n",
                      f"chr1\ts\tCDS\t{first + 1}\t{donor}\t.\t+\t0\tID={transcript}.c1;Parent={transcript}\n",
                      f"chr1\ts\tCDS\t{acceptor + 1}\t{len(sequence)}\t.\t+\t0\tID={transcript}.c2;Parent={transcript}\n"]
        paths["gff"].write_text("##gff-version 3\n" + "".join(lines))
        rows.append({"species": name, "genetic_code": 1, **{role: str(path) for role, path in paths.items()}})
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    tsv(inputs, list(rows[0]), rows)
    links = [dict(species_a=names[0], gene_a=names[0] + "_g", species_b=name, gene_b=name + "_g")
             for name in names[1:]]
    tsv(edges, list(links[0]), links)
    output = tmp_path / "refinement"
    value = refinement.plan(output, inputs=inputs, edges=edges, mode="off", policy=settings.get('policy', 'conserved'))
    effective = refinement.finalize(output, value)
    return effective, coding, common


def test_exported_short_model_agrees_with_protein_cds_gff_and_synteny(selected_view, tmp_path):
    effective, expected_cds, expected_protein = selected_view
    manifest = effective / "inputs.tsv"
    rows = load_effective_inputs(manifest)
    target = rows["Target_species"]
    selected = refinement.read_table(effective / "representative_map.tsv")
    target_choice = next(row for row in selected if row["species"] == "Target_species")
    assert target_choice["source_transcript_id"] == "short"
    assert target_choice["status"] == "conserved"
    assert Path(target["cds"]).read_text() == ">Target_species_g\n" + expected_cds + "\n"
    assert Path(target["protein"]).read_text() == ">Target_species_g\n" + expected_protein + "\n"
    assert str(Seq(expected_cds).translate()).removesuffix("*") == expected_protein
    columns = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    output = ["gene_id", "feature_size", "num_intron", "intron_positions", "chromosome", "start", "end", "strand",
              "feature_blocks", "feature_type", "gff_transcript_id", "cds_first_phase", "phase_status"]
    path = Path(target["gff"])
    traits = process_single_gff(path.name, str(path.parent), ["Target_species_g"], "CDS", "longest",
                                columns, output, representative_map=target["representative_map"])
    row = traits.iloc[0]
    assert row.gff_transcript_id == "short" and row.feature_size == len(expected_cds)
    assert row.num_intron == 1 and row.intron_positions == "120"
    assert row.start == 106 and row.cds_first_phase == 0 and row.phase_status == "consistent"
    source = dict(species="Target_species", fasta=target["cds"], gff=target["gff"], mode="cds",
                  feature="", attribute="", genetic_code=1, representative_map=target["representative_map"])
    genes, _metadata = prepare_genome(source, tmp_path, "target", 1)
    assert genes[0].start == 105 and genes[0].end == row.end
    assert (tmp_path / "target.pep").read_text() == Path(target["protein"]).read_text()
    full = (effective / "full_annotation/Target_species.gff3").read_text()
    assert "ID=long;Parent=g" in full and "ID=short;Parent=g" in full


@pytest.mark.parametrize('selected_view,expected_transcript,expected_start', [
    ({'transcript_ids': ('t,1', 't%2C1'), 'policy': 'conserved'}, 't,1', 106),
    ({'transcript_ids': ('t,1', 't%2C1'), 'policy': 'longest'}, 't%2C1', 1),
], indirect=['selected_view'])
def test_percent_escaped_source_identities_remain_distinct_across_selected_readers(
        selected_view, tmp_path, expected_transcript, expected_start):
    effective, _cds, _protein = selected_view
    row = load_effective_inputs(effective / 'inputs.tsv')['Target_species']
    selection = next(choice for choice in refinement.read_table(row['representative_map'])
                     if choice['species'] == 'Target_species')
    assert selection['source_transcript_id'] == expected_transcript
    full = effective / 'full_annotation/Target_species.gff3'
    assert 'ID=t%2C1;' in full.read_text() and 'ID=t%252C1;' in full.read_text()
    columns = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    output = ['gene_id', 'feature_size', 'num_intron', 'chromosome', 'start', 'end', 'strand',
              'feature_blocks', 'feature_type', 'gff_transcript_id']
    # Both source models coexist here; the literal-percent ID must never be
    # treated as an alias for the comma ID or as two Parent list members.
    traits = process_single_gff(full.name, str(full.parent), ['Target_species_g'], 'CDS', 'longest',
                               columns, output, representative_map=row['representative_map'])
    trait = traits.iloc[0]
    assert trait.gff_transcript_id == expected_transcript
    assert trait.start == expected_start
    source = dict(species='Target_species', fasta=row['cds'], gff=row['gff'], mode='cds',
                  feature='', attribute='', genetic_code=1, representative_map=row['representative_map'])
    genes, _metadata = prepare_genome(source, tmp_path, 'target', 1)
    assert (genes[0].start, genes[0].end) == (expected_start - 1, trait.end)
    assert (tmp_path / 'target.pep').read_text() == Path(row['protein']).read_text()
    cache = ensure_species_gene_cache('Target_species', row['cds'], str(Path(row['gff']).parent),
                                     str(tmp_path / 'cache'), str(tmp_path / 'locks'),
                                     str(SUPPORT / 'gff2genestat.py'), 1,
                                     representative_map=row['representative_map'])
    cached = refinement.read_table(cache)[0]
    assert cached['gff_transcript_id'] == expected_transcript
    assert int(cached['start']) == expected_start and int(cached['end']) == trait.end


@pytest.mark.parametrize("name", ["gg_genome_annotation", "gg_fractionation_bias", "gg_transcriptome_generation"])
def test_optional_view_config_is_registered_for_container_forwarding(name):
    entrypoint = name + "_entrypoint.sh"
    result = subprocess.run(["bash", "-c", 'source "$1"\ngg_print_entrypoint_config_vars "$2"',
                             "config", str(SUPPORT / "gg_entrypoint_config_vars.sh"), entrypoint],
                            capture_output=True, text=True, check=True, timeout=10)
    assert "representative_inputs" in result.stdout.splitlines()
    source = (ROOT / "workflow" / entrypoint).read_text()
    assert 'representative_inputs="' in source
    assert 'forward_config_vars_to_container_env "${gg_entrypoint_name}"' in source


@pytest.mark.parametrize("consumer", ["annotation", "fractionation", "quantification"])
def test_actual_core_input_routing_uses_verified_view(selected_view, tmp_path, consumer):
    effective, _cds, _protein = selected_view
    environment = dict(os.environ, gg_support_dir=str(SUPPORT), representative_inputs=str(effective / "inputs.tsv"),
                       gg_workspace_input_dir=str(tmp_path / "raw/input"), gg_workspace_output_dir=str(tmp_path / "raw/output"),
                       sp_ub="Target_species", target_species="Target_species", query_species="Donor_one",
                       analysis_id="pair", analysis_mode="compare", kallisto_reference="species_cds", run_amalgkit_quant="1")
    if consumer == "annotation":
        text = (ROOT / "workflow/core/gg_genome_annotation_core.sh").read_text()
        script = text[text.index('dir_sp_cds="'):text.index('dir_sp_dnaseq="')]
        result = 'printf "%s\\t%s\\t%s\\n" "$dir_sp_cds" "$dir_sp_gff" "$representative_map"\n'
        expected = [str(effective / "analysis_cds"), str(effective / "analysis_gff"), str(effective / "representative_map.tsv")]
    elif consumer == "fractionation":
        text = (ROOT / "workflow/core/gg_fractionation_bias_core.sh").read_text()
        function = text[text.index("resolve_species_file() {"):text.index("\n}\n", text.index("resolve_species_file() {")) + 3]
        script = function + "\n" + text[text.index('dir_sp_cds="'):text.index('file_genes="')]
        result = 'printf "%s\\t%s\\t%s\\n" "$target_cds" "$target_gff" "$query_cds"\n'
        expected = [str(effective / "species_cds/Target_species.fa"), None,
                    str(effective / "species_cds/Donor_one.fa")]
    else:
        text = (ROOT / "workflow/core/gg_transcriptome_generation_core.sh").read_text()
        start = text.index('file_kallisto_reference_fasta=""')
        end = text.index('if [[ ${run_amalgkit_quant} -eq 1 ]] &&', start)
        script = text[start:end]
        result = 'printf "%s\\n" "$file_kallisto_reference_fasta"\n'
        expected = [str(effective / "species_cds/Target_species.fa")]
    completed = subprocess.run(["bash", "-c", 'set -euo pipefail\nsource "$gg_support_dir/gg_util.sh"\n' + script + result],
                                env=environment, text=True, capture_output=True, timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    actual = completed.stdout.strip().split("\t")
    if consumer == 'fractionation':
        reader = Path(actual[1])
        assert reader.name == 'reader.gff3' and reader.is_relative_to(tmp_path / 'raw/output/.gg_cache/representative')
        receipt = json.loads((reader.parent / 'receipt.json').read_text())
        assert receipt['identity']['inputs']['gff']['path'] == str(effective / 'species_gff/Target_species.gff3')
        assert receipt['identity']['inputs']['cds']['path'] == expected[0]
        assert receipt['reader_feature'] == 'gene' and receipt['reader_attribute'] == 'ID'
        assert [(row['gene_id'], row['source_transcript_id']) for row in receipt['mapping']] == [('Target_species_g', 'short')]
        assert 'ID=Target_species_g;' in reader.read_text()
        expected[1] = str(reader)
    assert actual == expected
    # Tampered selected files must stop before a raw-input fallback is consumed.
    selected_cds = effective / "species_cds/Target_species.fa"
    selected_cds.write_text(selected_cds.read_text() + "A\n")
    failed = subprocess.run(["bash", "-c", 'set -euo pipefail\nsource "$gg_support_dir/gg_util.sh"\n' + script + result],
                            env=environment, text=True, capture_output=True, timeout=30)
    assert failed.returncode != 0
    assert "Effective input changed" in failed.stderr


@pytest.mark.parametrize('selected', [False, True])
def test_gene_evolution_only_applies_original_resolution_to_legacy_inputs(tmp_path, selected):
    text = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    start = text.index('  gff_cds_validation_args=()')
    script = text[start:text.index('  python "${gg_support_dir}/gff2genestat.py"', start)]
    environment = dict(os.environ, input_sequence_mode='cds', gg_workspace_output_dir=str(tmp_path),
                       representative_inputs='selected-inputs.tsv' if selected else '')
    result = subprocess.run(['bash', '-c', 'set -euo pipefail\n' + script +
                             '\nprintf "%s\\n" "${gff_cds_validation_args[@]}"'],
                            env=environment, text=True, capture_output=True, check=True, timeout=10)
    arguments = result.stdout.splitlines()
    assert '--validate-cds-length' in arguments
    if selected:
        assert '--cds-resolution-dir' not in arguments
    else:
        index = arguments.index('--cds-resolution-dir')
        assert arguments[index + 1] == str(tmp_path / 'species_cds_resolved')


def test_gene_evolution_gff_trait_cache_binds_reader_and_selection_implementations():
    text = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    stage = text[text.index('task="Gene trait extraction from gff files"'):]
    stage = stage[:stage.index('gg_artifact_prepare_stage gff_info_needs_update')]
    for helper in ('gff2genestat.py', 'representative_selection.py', 'gff_feature_structure.py'):
        assert any('--input ' in line and '${gg_support_dir}/' + helper in line for line in stage.splitlines())


@pytest.mark.parametrize('selected', [False, True])
def test_genome_evolution_resolves_codes_from_the_routed_input_table(tmp_path, selected):
    core = (ROOT / 'workflow/core/gg_genome_evolution_core.sh').read_text()
    functions = core[core.index('species_genetic_code_table_path() {'):core.index('lookup_species_genetic_code() {')]
    # Execute the actual production call as well as its function: a correctly
    # parameterized helper alone would not catch an omitted core argument.
    call = next(line for line in core.splitlines()
                if line.startswith('  prepare_species_genetic_code_table "${dir_sp_cds}"'))
    cds_dir = tmp_path / 'effective/species_cds'
    cds_dir.mkdir(parents=True)
    (cds_dir / 'Target_species.fa').write_text('>Target_species_g\nATGTGATAA\n')
    raw_dir = tmp_path / 'raw/species_genetic_code'
    raw_dir.mkdir(parents=True)
    raw_table = raw_dir / 'species_genetic_code.tsv'
    raw_table.write_text('species\tgenetic_code\nTarget_species\t1\n')
    selected_table = tmp_path / 'effective/species_genetic_code.tsv'
    selected_table.write_text('species\tgenetic_code\nTarget_species\t4\n')
    resolved = tmp_path / 'resolved.tsv'
    environment = dict(os.environ, gg_workspace_input_dir=str(raw_dir.parent), dir_sp_cds=str(cds_dir),
                       genetic_code='1', file_species_genetic_code_resolved=str(resolved),
                       file_species_genetic_code=str(selected_table if selected else raw_table))
    result = subprocess.run(['bash', '-c', 'set -euo pipefail\n' + functions + '\n' + call],
                            env=environment, text=True, capture_output=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
    with resolved.open(newline='') as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['species'] == 'Target_species'
    assert row['genetic_code'] == ('4' if selected else '1')
    if selected:
        translated = subprocess.run(['seqkit', 'translate', '--transl-table', row['genetic_code'],
                                     str(cds_dir / 'Target_species.fa')],
                                    capture_output=True, text=True, check=True, timeout=10)
        assert translated.stdout == '>Target_species_g\nMW*\n'


@pytest.mark.parametrize('change', ['genetic_code', 'species', 'subset', 'corrupt_archive', 'missing_receipt'])
def test_shared_effective_reader_binds_copied_manifests_to_publication(selected_view, tmp_path, change):
    effective, _cds, _protein = selected_view
    manifest = effective / 'inputs.tsv'
    copied = tmp_path / 'copied.tsv'
    copied.write_bytes(manifest.read_bytes())
    assert set(load_effective_inputs(copied)) == {'Target_species', 'Donor_one', 'Donor_two'}
    rows = refinement.read_table(copied)
    if change in {'genetic_code', 'species'}:
        rows[0][change] = '4' if change == 'genetic_code' else 'Different_species'
    elif change == 'subset':
        rows = rows[:1]
    elif change == 'corrupt_archive':
        path = effective / 'source_cds/Target_species.fa'
        path.write_text(path.read_text() + 'A\n')
    else:
        (effective / 'receipt.json').unlink()
    tsv(copied, list(rows[0]), rows)
    with pytest.raises(ValueError):
        load_effective_inputs(copied)


def test_effective_cds_pairwise_uses_admitted_partial_proteins_and_excludes_negatives(tmp_path, monkeypatch):
    names = ['Target_species', 'Donor_one', 'Donor_two']
    cases = [('partial', 'TATGAAATAA', 1, ''), ('trailing', 'ATGAAATAAA', 0, ''),
             ('stop', 'ATGTGATAA', 0, ''), ('pseudo', 'ATGAAATAA', 0, ';pseudo=true')]
    sources = []
    for name in names:
        paths = {role: tmp_path / (name + suffix) for role, suffix in
                 [('cds', '.fa'), ('gff', '.gff3'), ('genome', '.genome.fa')]}
        coding, annotation, genome = [], ['##gff-version 3\n'], []
        for gene, dna, phase, extra in cases:
            coding.append(f'>{gene}\n{dna}\n')
            genome.append(f'>chr_{gene}\n{dna}\n')
            annotation += [f'chr_{gene}\ts\tgene\t1\t{len(dna)}\t.\t+\t.\tID={gene}{extra}\n',
                           f'chr_{gene}\ts\tmRNA\t1\t{len(dna)}\t.\t+\t.\tID=t{gene};Parent={gene}\n',
                           f'chr_{gene}\ts\tCDS\t1\t{len(dna)}\t.\t+\t{phase}\tParent=t{gene}\n']
        for role, lines in [('cds', coding), ('gff', annotation), ('genome', genome)]:
            paths[role].write_text(''.join(lines))
        sources.append(dict(species=name, genetic_code=1, **{role: str(path) for role, path in paths.items()}))
    inputs, edges = tmp_path / 'inputs.tsv', tmp_path / 'edges.tsv'
    tsv(inputs, list(sources[0]), sources)
    tsv(edges, ['species_a', 'gene_a', 'species_b', 'gene_b'], [])
    output = tmp_path / 'refinement'
    value = refinement.plan(output, inputs=inputs, edges=edges, policy='longest', mode='off')
    effective = refinement.finalize(output, value)
    pairs = tmp_path / 'pairs.tsv'
    tsv(pairs, ['analysis_id', 'target_species', 'query_species'],
        [dict(analysis_id='pair', target_species=names[0], query_species=names[1])])
    module = importlib.import_module('pairwise_synteny')
    monkeypatch.setattr(module, 'tool_identity', lambda: {})
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode='cds', genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                              formats='svg', karyotype_sort='both_length', representative_inputs=effective / 'inputs.tsv')
    plan = module.build_plan(args)
    source = plan['pairs'][0]['target']
    assert source['admitted_protein'] == str(effective / 'species_protein/Target_species.fa')
    result = tmp_path / 'comparison'
    result.mkdir()
    genes, metadata = prepare_genome(source, result, 'target', 1)
    assert {gene.gene_id for gene in genes} == {'Target_species_partial', 'Target_species_trailing'}
    assert (result / 'target.pep').read_bytes() == Path(source['admitted_protein']).read_bytes()
    assert metadata['admitted_protein_sha256'] == module.digest(source['admitted_protein'])
    admitted_input = 'pair.target.admitted_protein=' + source['admitted_protein']
    assert admitted_input in module.contract_args(plan, 'analysis')
    # Keep every biological source base and archive the excluded CDS paths.
    target_row = load_effective_inputs(effective / 'inputs.tsv')['Target_species']
    raw = Path(target_row['cds']).read_text()
    assert 'TATGAAATAA' in raw and '>Target_species_stop\nATGTGATAA' in raw
    id_map = refinement.read_table(result / 'target.id_map.tsv')
    if 'analysis_cds' in target_row:
        assert source['fasta'] == target_row['analysis_cds'] and source['gff'] == target_row['analysis_gff']
        assert {row['original_id'] for row in id_map} == {'Target_species_partial', 'Target_species_trailing'}
    else:
        assert {row['original_id'] for row in id_map if row['status'] == 'translation_excluded'} == {
            'Target_species_stop', 'Target_species_pseudo'}


def test_synteny_selected_view_uses_admitted_proteins_and_binds_genetic_codes():
    core = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    stage = core[core.index('task="Synteny neighborhood grouping"'):core.index('task="summary statistics"', core.index('task="Synteny neighborhood grouping"'))]
    assert 'if [[ -n "${representative_inputs}" ]]; then' in stage
    assert 'synteny_source_dir="${dir_sp_protein_input}"' in stage
    assert '--genetic-codes "${file_species_genetic_code}"' in stage
    assert '--input "species_genetic_code=${file_species_genetic_code}"' in stage
    assert '--input "representative_map=${representative_map}"' in stage
    for helper in ('synteny_neighbors.py', 'gff2genestat.py', 'representative_selection.py', 'gff_feature_structure.py'):
        assert any('--input ' in line and '${gg_support_dir}/' + helper in line for line in stage.splitlines())


@pytest.mark.parametrize('mode', ['cds', 'protein'])
@pytest.mark.parametrize('selected', [False, True])
def test_gene_core_cache_paths_bind_selected_manifest_and_sequence_view(tmp_path, mode, selected):
    core = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    script = core[core.index('dir_sp_blastdb="'):core.index('file_species_genetic_code="')]
    manifest = tmp_path / 'inputs.tsv'
    manifest.write_text('species\tgenetic_code\nTarget_species\t1\n')
    output = tmp_path / 'output'
    environment = dict(os.environ, gg_workspace_output_dir=str(output), gg_support_dir=str(SUPPORT),
                       representative_inputs=str(manifest) if selected else '', input_sequence_mode=mode)
    result = subprocess.run(['bash', '-c', 'set -euo pipefail\n' + script +
                             '\nprintf "%s\\n" "$file_species_cds_store_db" "$file_species_cds_store_manifest" '
                             '"$file_species_protein_store_db" "$file_species_protein_store_manifest" "$dir_synteny_gene_cache" "$dir_sp_blastdb"'],
                            env=environment, text=True, capture_output=True, check=True, timeout=10)
    base = output / '.gg_cache/fasta_sequence_store'
    if selected:
        base = base / 'representative' / hashlib.sha256(manifest.read_bytes()).hexdigest() / mode
    expected = [str(base / name) for name in ['species_cds.sqlite3', 'species_cds.json',
                                             'species_protein.sqlite3', 'species_protein.json']]
    expected.append(str(base / 'synteny_gene_info' if selected else output / 'species_gff_info'))
    expected.append(str(base / 'species_cds_blastdb' if selected else output / 'species_cds_blastdb'))
    assert result.stdout.splitlines() == expected
    assert '--cache_dir "${dir_synteny_gene_cache}"' in core


@pytest.mark.parametrize('selected', [False, True])
def test_annotation_core_isolates_selected_resolution_and_gene_info_caches(tmp_path, selected):
    core = (ROOT / 'workflow/core/gg_genome_annotation_core.sh').read_text()
    namespace = core[core.index('representative_annotation_cache="'):core.index('dir_sp_dnaseq="')]
    gene_info = core[core.index('file_sp_gff_info="'):core.index('file_sp_cds_busco_full="')]
    resolution = core[core.index('cds_resolution_dir="'):core.index('cds_resolution_args=')]
    receipt = core[core.index('annotation_provenance_dir="'):core.index('ensure_dir "${dir_sp_tmp}"')]
    manifest = tmp_path / 'inputs.tsv'
    manifest.write_text('species\tgenetic_code\nTarget_species\t1\n')
    output = tmp_path / 'output'
    environment = dict(os.environ, gg_workspace_output_dir=str(output), gg_support_dir=str(SUPPORT),
                       representative_inputs=str(manifest) if selected else '', sp_ub='Target_species')
    result = subprocess.run(['bash', '-c', 'set -euo pipefail\n' + namespace + gene_info + resolution + receipt +
                             '\nprintf "%s\\n" "$file_sp_gff_info" "$cds_resolution_dir" "$gff_info_provenance_file"'],
                            env=environment, text=True, capture_output=True, check=True, timeout=10)
    base = (output / '.gg_cache/genome_annotation/representative'
            / hashlib.sha256(manifest.read_bytes()).hexdigest() / 'cds') if selected else output
    assert result.stdout.splitlines() == [str(base / 'species_gff_info/Target_species_gff_info.tsv'),
                                         str(base / 'species_cds_resolved'),
                                         str(base / 'artifact_provenance/Target_species.gff_info.json' if selected else
                                             output / 'artifact_provenance/genome_annotation/Target_species.gff_info.json')]
    # The namespace must be chosen after verified coding-view routing.
    assert core.index('--field analysis_layout') < core.index('representative_annotation_cache="')


def test_selected_synteny_caches_survive_interleaved_bundle_preparation(tmp_path):
    namespace = importlib.import_module('fasta_sequence_store').selected_input_namespace
    outputs = []
    for label, start in [('a', 1), ('b', 21)]:
        source_dir = tmp_path / label
        source_dir.mkdir()
        fasta = source_dir / 'Target_species.fa'
        fasta.write_text('>Target_species_g\nMK\n')
        gff = source_dir / 'Target_species.gff3'
        gff.write_text(f'chr1\ts\tgene\t{start}\t{start + 8}\t.\t+\t.\tID=g\n'
                       f'chr1\ts\tmRNA\t{start}\t{start + 8}\t.\t+\t.\tID=t;Parent=g\n'
                       f'chr1\ts\tCDS\t{start}\t{start + 8}\t.\t+\t0\tParent=t\n')
        selection = source_dir / 'representative_map.tsv'
        selection.write_text('species\tgene_id\tcandidate_id\tsource_transcript_id\tstatus\tscore\tmargin\treason\n'
                             'Target_species\tTarget_species_g\tc\tt\tselected\t1\t0.2\tconserved\n')
        manifest = source_dir / 'inputs.tsv'
        manifest.write_text('species\tgff\nTarget_species\t' + str(gff) + '\n')
        cache_dir = namespace(tmp_path / 'cache', manifest, 'protein') / 'synteny_gene_info'
        cached = ensure_species_gene_cache('Target_species', str(fasta), str(source_dir),
                                           str(cache_dir), str(cache_dir / '.locks'),
                                           str(SUPPORT / 'gff2genestat.py'), 1,
                                           representative_map=str(selection))
        outputs.append(Path(cached))
    assert outputs[0] != outputs[1]
    # Read A only after B has completed its cache publication.
    assert [int(refinement.read_table(path)[0]['start']) for path in outputs] == [1, 21]


@pytest.mark.parametrize('method', ['diamond', 'tblastn'])
def test_selected_query_and_native_search_share_admitted_mixed_code_proteins(tmp_path, method):
    """Never retranslate partial CDS or put excluded biological models in a DB."""
    codons = dict(zip('ACDEFGHIKLMNPQRSTVWY', ['GCT', 'TGT', 'GAT', 'GAA', 'TTT', 'GGT', 'CAT', 'ATT',
                                             'AAA', 'CTG', 'ATG', 'AAT', 'CCT', 'CAA', 'CGT', 'TCT',
                                             'ACT', 'GTT', 'TGG', 'TAT'], strict=True))
    rng = random.Random(811)
    protein = 'M' + ''.join(rng.choice(list(codons)) for _ in range(119))
    protein = protein[:25] + 'WWWW' + protein[29:]
    sources = []
    for species, code in [('Target_species', 4), ('Donor_one', 1)]:
        dna = ''.join(('TGA' if residue == 'W' and code == 4 else codons[residue])
                      for residue in protein) + 'TAA'
        cases = [('partial', 'T' + dna, 1, ''),
                 ('stop', dna[:150] + 'TAA' + dna[153:], 0, ''),
                 ('pseudo', dna, 0, ';pseudo=true')]
        paths = {role: tmp_path / (species + suffix) for role, suffix in
                 [('cds', '.fa'), ('gff', '.gff3'), ('genome', '.genome.fa')]}
        coding, genome, annotation = [], [], ['##gff-version 3\n']
        for gene, sequence, phase, extra in cases:
            coding.append(f'>{gene}\n{sequence}\n')
            genome.append(f'>chr_{gene}\n{sequence}\n')
            annotation += [f'chr_{gene}\ts\tgene\t1\t{len(sequence)}\t.\t+\t.\tID={gene}{extra}\n',
                           f'chr_{gene}\ts\tmRNA\t1\t{len(sequence)}\t.\t+\t.\tID=t{gene};Parent={gene}\n',
                           f'chr_{gene}\ts\tCDS\t1\t{len(sequence)}\t.\t+\t{phase}\tParent=t{gene}\n']
        for role, lines in [('cds', coding), ('genome', genome), ('gff', annotation)]:
            paths[role].write_text(''.join(lines))
        sources.append(dict(species=species, genetic_code=code, **{role: str(path) for role, path in paths.items()}))
    inputs, edges = tmp_path / 'inputs.tsv', tmp_path / 'edges.tsv'
    tsv(inputs, list(sources[0]), sources)
    tsv(edges, ['species_a', 'gene_a', 'species_b', 'gene_b'], [])
    output = tmp_path / 'refinement'
    value = refinement.plan(output, inputs=inputs, edges=edges, policy='longest', mode='off')
    effective = refinement.finalize(output, value)
    rows = load_effective_inputs(effective / 'inputs.tsv')
    records = importlib.import_module('fasta_sequence_store').fasta_records
    for species, row in rows.items():
        assert {key: sequence for key, _header, sequence in records(Path(row['protein']))} == {
            species + '_partial': protein}
        assert 'pseudo' in Path(row['cds']).read_text() and 'stop' in Path(row['cds']).read_text()

    core = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    stores = core[core.index('dir_sp_blastdb="'):core.index('file_species_genetic_code="')]
    ensure_store = core[core.index('ensure_species_fasta_sequence_store() {'):
                        core.index('# shellcheck shell=bash', core.index('ensure_species_fasta_sequence_store() {'))]
    resolve_evalue = core[core.index('resolve_query_blast_evalue() {'):core.index('prepare_synteny_evalue_fasta() {')]
    stages = core[core.index('task="Query fasta generation"'):core.index('task="Fasta generation"')]
    query = tmp_path / 'query.txt'
    query.write_text('Target_species_partial\n')
    work = tmp_path / 'work'
    work.mkdir()
    active = tmp_path / 'output'
    active.mkdir()
    environment = dict(os.environ, gg_support_dir=str(SUPPORT), gg_workspace_dir=str(tmp_path),
                       gg_workspace_output_dir=str(active), dir_output_active=str(active), dir_tmp=str(work),
                       representative_inputs=str(effective / 'inputs.tsv'), input_sequence_mode='protein',
                       dir_sp_cds=str(effective / 'species_cds'), dir_sp_protein_input=str(effective / 'species_protein'),
                       file_species_genetic_code=str(effective / 'species_genetic_code.tsv'),
                       file_query_gene=str(query), file_og_query_aa_fasta=str(active / 'query.fa.gz'),
                       file_og_query_blast=str(active / 'blast.tsv'), og_id='OGtest', genetic_code='1',
                       mode_gene_evolution='query2family', query_blast_method=method, query_blast_evalue='1e-3',
                       query_blast_auto_evalue_maxlen_cutoffs='inf:1e-3', GG_TASK_CPUS='2',
                       run_extract_query_fasta='1', run_query_blast='1')
    script = ('set -euo pipefail\nsource "${gg_support_dir}/gg_util.sh"\n' + stores + ensure_store +
              resolve_evalue + r'''
gg_step_start() { :; }
gg_step_skip() { :; }
gg_artifact_prepare_stage() { printf -v "$1" '%s' 1; }
gg_artifact_record() { printf '%s\n' "$@" >> "${dir_output_active}/provenance.txt"; }
# Selected ID queries must never invoke the CDS phase/padding translator.
gg_prepare_cds_fasta_stream() { echo 'Unexpected CDS retranslation' >&2; return 99; }
# Save the real builder's input; every DIAMOND operation still runs natively.
diamond() {
  if [[ "$1" == makedb ]]; then
    local reference='' database='' previous='' value
    for value in "$@"; do
      [[ "${previous}" != --in ]] || reference="${value}"
      [[ "${previous}" != --db ]] || database="${value}"
      previous="${value}"
    done
    cp -- "${reference}" "${database}.reference.fa"
  fi
  command diamond "$@"
}
makeblastdb() {
  local reference='' database='' previous='' value
  for value in "$@"; do
    [[ "${previous}" != -in ]] || reference="${value}"
    [[ "${previous}" != -out ]] || database="${value}"
    previous="${value}"
  done
  cp -- "${reference}" "${database}.reference.fa"
  command makeblastdb "$@"
}
tblastn() {
  local database='' code='' threads='' dbsize='' previous='' value
  for value in "$@"; do
    [[ "${previous}" != -db ]] || database="${value}"
    [[ "${previous}" != -db_gencode ]] || code="${value}"
    [[ "${previous}" != -num_threads ]] || threads="${value}"
    [[ "${previous}" != -dbsize ]] || dbsize="${value}"
    previous="${value}"
  done
  printf '%s\t%s\t%s\t%s\n' "$(basename "${database}")" "${code}" "${threads}" "${dbsize}" >> "${dir_output_active}/tblastn.calls.tsv"
  command tblastn "$@"
}
cd "${dir_tmp}"
''' + stages)
    result = subprocess.run(['bash', '-c', script], env=environment, text=True,
                            capture_output=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
    with gzip.open(active / 'query.fa.gz', 'rt') as handle:
        header, *lines = handle.read().splitlines()
        assert header == '>Target_species_partial' and ''.join(lines) == protein
    hits = (active / 'blast.tsv').read_text()
    assert 'Target_species_partial' in hits and 'Donor_one_partial' in hits
    assert 'pseudo' not in hits and 'stop' not in hits
    with (active / 'blast.tsv').open() as handle:
        matches = list(csv.DictReader(handle, delimiter='\t'))
    assert len(matches) == 2
    assert all(float(match['pident']) == 100 and int(match['length']) == len(protein) for match in matches)
    db_root = (active / '.gg_cache/fasta_sequence_store/representative'
               / hashlib.sha256((effective / 'inputs.tsv').read_bytes()).hexdigest() / 'protein/species_cds_blastdb')
    databases = sorted(path for path in db_root.glob('*.dmnd' if method == 'diamond' else '*.n*')
                       if path.is_file() and not path.name.startswith('.'))
    assert (len(databases) == 2) if method == 'diamond' else (len(databases) >= 8)
    before = {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in databases}
    for species, row in rows.items():
        database = db_root / (species + '.fa')
        reference = Path(str(database) + '.reference.fa')
        expected_source = Path(row['protein'] if method == 'diamond' else row['analysis_cds'])
        assert {key: sequence for key, _header, sequence in records(reference)} == {
            key: sequence for key, _header, sequence in records(expected_source)}
        signature = Path(str(database) + '.' + method + '.build.signature').read_text()
        assert 'source_sha256=' + hashlib.sha256(expected_source.read_bytes()).hexdigest() in signature
        reference_kind = 'admitted_protein' if method == 'diamond' else 'admitted_analysis_cds'
        assert ';reference=' + reference_kind in signature and 'translation_filter' not in signature
    if method == 'tblastn':
        total_nt = sum(len(sequence) for row in rows.values() for _key, _header, sequence in records(Path(row['analysis_cds'])))
        assert set((active / 'tblastn.calls.tsv').read_text().splitlines()) == {
            f'Target_species.fa\t4\t1\t{total_nt}', f'Donor_one.fa\t1\t1\t{total_nt}'}
        assert 'species_genetic_code=' in (active / 'provenance.txt').read_text()
        assert f'tblastn_database_size={total_nt}' in (active / 'provenance.txt').read_text()
    assert 'species_protein_index=' in (active / 'provenance.txt').read_text()
    assert 'fasta_sequence_store_schema=2' in (active / 'provenance.txt').read_text()
    # Identical selected inputs must reuse the same real DB after build locks release.
    again = subprocess.run(['bash', '-c', script], env=environment, text=True,
                           capture_output=True, timeout=45)
    assert again.returncode == 0, again.stdout + again.stderr
    assert {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in databases} == before
    # A biological negative queried by ID must fail instead of falling back to raw CDS.
    query.write_text('Target_species_pseudo\n')
    negative = subprocess.run(['bash', '-c', script], env=environment, text=True,
                              capture_output=True, timeout=20)
    assert negative.returncode != 0
    assert 'Query gene not found in the protein inventory: Target_species_pseudo' in negative.stdout
