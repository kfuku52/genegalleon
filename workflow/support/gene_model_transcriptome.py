"""Import assembly coding isoforms without assigning genomic coordinates."""

import csv
import json
import re
import shutil
import subprocess
from pathlib import Path

from Bio.Data import CodonTable

if __package__:
    from . import rescue_gene_models as rescue
    from .fasta_sequence_store import fasta_records, open_text
    from .gene_model_catalog import validate_candidate
    from .input_generation_array_state import FreshDigestBatch, digest
    from .rename_rnabloom_transcripts import cluster_transcripts, load_clusters
else:
    import rescue_gene_models as rescue
    from fasta_sequence_store import fasta_records, open_text
    from gene_model_catalog import validate_candidate
    from input_generation_array_state import FreshDigestBatch, digest
    from rename_rnabloom_transcripts import cluster_transcripts, load_clusters


def _rows(path, required, optional=()):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not set(required) <= set(columns) <= set(required) | set(optional):
            raise ValueError(f'{path} requires {", ".join(required)} columns')
        for row in reader:
            if None in row or None in row.values() or any(not row[key].strip() for key in required):
                raise ValueError(f'Incomplete TSV row at {path}:{reader.line_num}')
            yield {key: value.strip() for key, value in row.items()}


def _genetic_code(value):
    code = int(value)
    if code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f'Unsupported genetic code: {code}')
    return code


def read_sources(path):
    path = Path(path).resolve()
    sources = {}
    for row in _rows(path, ('species', 'cds', 'transcript_gene_map'), ('genetic_code',)):
        species = rescue.safe_token(row['species'], 'Transcriptome species')
        if species in sources:
            raise ValueError(f'Duplicate transcriptome species: {species}')
        sources[species] = {'fasta': str((path.parent / row['cds']).resolve()),
                            'transcript_gene_map': str((path.parent / row['transcript_gene_map']).resolve()),
                            'genetic_code': _genetic_code(row.get('genetic_code') or 1)}
    if not sources:
        raise ValueError(f'No transcriptome species in {path}')
    return sources


def _gene_mapping(path, species):
    mapping, baselines = {}, set()
    for row in _rows(path, ('transcript_id', 'gene_id'), ('baseline',)):
        transcript, gene = row['transcript_id'], rescue.safe_token(row['gene_id'], 'Assembly gene ID')
        if not gene.startswith(species + '_'):
            raise ValueError(f'Assembly gene ID must start with {species}_: {gene}')
        if transcript in mapping or any(char.isspace() for char in transcript):
            raise ValueError(f'Duplicate or invalid assembly transcript ID: {transcript}')
        baseline = row.get('baseline', '')
        if baseline not in ('', '0', '1'):
            raise ValueError(f'Baseline must be 0 or 1 for {transcript}')
        if baseline == '1':
            if gene in baselines:
                raise ValueError(f'Multiple baseline transcripts for {gene}')
            baselines.add(gene)
        mapping[transcript] = gene, baseline == '1'
    return mapping


def _mapped_cds(species, source):
    mapping = _gene_mapping(source['transcript_gene_map'], species)
    saw_record = False
    for transcript, _, sequence in fasta_records(Path(source['fasta'])):
        if transcript not in mapping:
            raise ValueError(f'Unmapped or duplicate transcriptome CDS ID: {transcript}')
        gene, baseline = mapping.pop(transcript)
        if baseline and not sequence:
            raise ValueError(f'Baseline transcript has no CDS: {transcript}')
        saw_record = True
        yield transcript, gene, baseline, sequence.upper()
    if mapping:
        raise ValueError(f'Transcriptome mapping references an absent CDS: {min(mapping)}')
    if not saw_record:
        raise ValueError(f'No transcriptome CDS records for {species}')


def build_catalog(species, source):
    code = source['genetic_code']
    loci = {}
    for transcript, gene, baseline, sequence in _mapped_cds(species, source):
        candidate = dict(candidate_id=transcript, source_transcript_id=transcript, source_gene_id=gene,
                         gene_id=gene, species=species, cds=sequence, blocks=[],
                         origin='original', sequence_only=True)
        candidate['quality'] = validate_candidate(candidate, code)
        locus = loci.setdefault(gene, {'species': species, 'gene_id': gene, 'candidates': []})
        locus['candidates'].append(candidate)
        if baseline:
            locus['baseline_candidate_id'] = transcript
    return {'species': species, 'genetic_code': code, 'sequence_only': True,
            'loci': [loci[gene] for gene in sorted(loci)],
            'summary': {'loci': len(loci), 'candidates': sum(len(g['candidates']) for g in loci.values())}}


def write_initial_cds(species, source, output):
    """Retain and translate only the same baseline each catalog locus would use."""
    representatives = {}
    for transcript, gene, baseline, sequence in _mapped_cds(species, source):
        # Explicit baseline first, then longest CDS and stable transcript ID.
        rank = not baseline, -len(sequence), transcript
        previous = representatives.get(gene)
        if previous is None or rank < previous[0]:
            representatives[gene] = rank, sequence
    with output.open('w') as stream:
        for gene in sorted(representatives):
            _, sequence = representatives[gene]
            quality = validate_candidate({'cds': sequence, 'sequence_only': True}, source['genetic_code'])
            if quality['usable']:
                stream.write(f'>{gene}\n{sequence[:len(sequence) // 3 * 3]}\n')


def _assembly_mapping(assembly, species, assembler, corset, gene_map):
    records = ((transcript, sequence) for transcript, _, sequence in fasta_records(assembly))
    if assembler == 'rna-bloom2':
        clusters = load_clusters(corset, strict=True)
        assigned = set(clusters)
        for transcript, gene, _, _ in cluster_transcripts(records, clusters, species):
            assigned.discard(transcript)
            yield transcript, gene, 0
        if assigned:
            raise ValueError(f'Corset cluster references an absent assembly transcript: {min(assigned)}')
    elif gene_map:
        mapping = _gene_mapping(gene_map, species)
        for transcript, _ in records:
            if transcript not in mapping:
                raise ValueError(f'Unmapped or duplicate assembly transcript: {transcript}')
            gene, baseline = mapping.pop(transcript)
            yield transcript, gene, int(baseline)
        if mapping:
            raise ValueError(f'Gene map references an absent assembly transcript: {min(mapping)}')
    else:
        pattern = (r'(?:TRINITY_)?(DN\d+_c\d+_g\d+)_i\d+' if assembler == 'trinity'
                   else r'.+_(g\d+)_i\d+')
        for transcript, _ in records:
            match = re.fullmatch(pattern, transcript)
            if not match:
                raise ValueError(f'Unrecognized {assembler} transcript ID: {transcript}. Supply --gene-map.')
            yield transcript, species + '_' + match[1].replace('_', '-'), 0


def prepare(assembly, species, output, code, cpus, assembler=None, corset=None, gene_map=None, *, input_proof=None):
    if not assembler and not gene_map:
        raise ValueError('Provide --assembler or --gene-map')
    if (assembler == 'rna-bloom2') != bool(corset):
        raise ValueError('--corset-clusters is required with RNA-Bloom2 and applies only to that assembler')
    if cpus < 1:
        raise ValueError('--cpus must be positive')
    code = _genetic_code(code)
    rescue.safe_token(species, 'Transcriptome species')
    assembly, output = Path(assembly).resolve(), Path(output).resolve()
    paths = [assembly, *(Path(path).resolve() for path in (corset, gene_map) if path)]
    if any(output == path or output in path.parents for path in paths):
        raise ValueError('Keep the assembly and mapping sources outside the output directory')
    executable = shutil.which('cdskit')
    if not executable:
        raise ValueError('cdskit is required to extract assembly ORFs')
    proof = input_proof if input_proof is not None else FreshDigestBatch()
    key = {'inputs': proof.read(paths), 'species': species, 'assembler': assembler, 'genetic_code': code,
           'implementation': {name: digest(Path(__file__).with_name(name)) for name in (
               Path(__file__).name, 'rename_rnabloom_transcripts.py', 'gene_model_catalog.py')},
           'cdskit': {**rescue.package_source_hashes('cdskit'), executable: digest(executable)}}

    def guard():
        proof.check()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            try:
                previous = json.loads((output / 'receipt.json').read_text())
            except (OSError, ValueError):
                previous = {}
            if not isinstance(previous, dict) or previous.get('key') != key:
                raise ValueError(f'Assembly candidate inputs differ or {output} belongs to another run. Use a new --output.')

    def build(tmp):
        # Use the same ORF caller as transcriptome generation, before sequence
        # deduplication, padding or collapsing transcripts by gene.
        local = tmp / 'assembly.fa'
        with open_text(assembly) as source, local.open('w') as out:
            shutil.copyfileobj(source, out)
        cds, mapping = tmp / 'all_isoforms.cds.fa', tmp / 'transcript_gene_map.tsv'
        with (tmp / 'longestorf.log').open('w') as log:
            subprocess.run([executable, 'longestorf', '--seq_file', str(local), '--out_file', str(cds),
                            '--codon_table', str(code), '--annotate_seq_name', 'no', '--threads', str(cpus)],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
        rescue.write_tsv(mapping, ('transcript_id', 'gene_id', 'baseline'),
                         _assembly_mapping(local, species, assembler, corset, gene_map))
        source = {'fasta': str(cds), 'transcript_gene_map': str(mapping), 'genetic_code': code}
        write_initial_cds(species, source, tmp / f'{species}_initial_longestCDS.fa')
        local.unlink()
        rescue.write_tsv(tmp / 'inputs.tsv', ('species', 'cds', 'transcript_gene_map', 'genetic_code'),
                         [(species, cds.name, mapping.name, code)])

    return rescue.stage(output.parent, output.name, key, build, guard=guard, publication_guard=guard)


def discover_assemblies(directory, excluded, proof):
    """Find standard GeneGalleon assembly outputs for species without genomic inputs."""
    directory = Path(directory).resolve()
    if not directory.is_dir():
        raise ValueError(f'Transcriptome assembly directory does not exist: {directory}')
    result = {}
    for path in sorted((directory / 'assembled_transcripts_with_isoforms').glob('*_isoform.fa*')):
        name = path.name.removesuffix('.gz').removesuffix('_isoform.fa')
        if name == path.name.removesuffix('.gz') or name in excluded or not path.is_file():
            continue
        rescue.safe_token(name, 'Transcriptome species')
        path = path.resolve()
        if name in result:
            raise ValueError(f'Multiple assembly files for {name}')
        corset = directory / 'corset_clusters' / f'{name}_corset.clusters.tsv'
        receipt = directory / 'artifact_provenance' / f'{name}.assembly.json'
        receipt = receipt if receipt.is_file() else None
        frozen = proof.read([path, receipt] if receipt else [path])
        if receipt:
            metadata = json.loads(receipt.read_text())
            output, = (row for row in metadata['outputs'] if row['label'] == 'isoforms')
            if metadata['family_id'] != name or output['sha256'] != frozen[str(path)]:
                raise ValueError(f'Assembly output differs from the assembly receipt: {name}')
            assembler = metadata['parameters']['assembly_method'].lower()
        else:
            first = next(fasta_records(path), ('', '', ''))[0]
            assembler = ('trinity' if re.fullmatch(r'(?:TRINITY_)?DN\d+_c\d+_g\d+_i\d+', first)
                         else 'rnaspades' if re.fullmatch(r'.+_g\d+_i\d+', first)
                         else 'rna-bloom2' if corset.is_file() else '')
        if assembler not in {'trinity', 'rnaspades', 'rna-bloom2'}:
            raise ValueError(f'Cannot determine assembly gene identifiers for {name}; supply prepared transcriptome inputs')
        result[name] = (path, assembler, corset if assembler == 'rna-bloom2' else None, receipt)
    proof.check()
    return result


def prepare_assemblies(directory, excluded, output, cpus, genetic_codes, proof):
    """Prepare and checksum all assembly-only species using the existing ORF workflow."""
    assemblies = discover_assemblies(directory, excluded, proof)
    if not assemblies:
        return {}, []
    codes = {}
    if genetic_codes:
        genetic_codes = Path(genetic_codes).resolve()
        proof.read([genetic_codes])
        for row in _rows(genetic_codes, ('species', 'genetic_code')):
            if row['species'] in codes:
                raise ValueError('Duplicate species genetic code: ' + row['species'])
            codes[row['species']] = _genetic_code(row['genetic_code'])
    files, sources = [], {}
    if genetic_codes:
        files.append(str(Path(genetic_codes).resolve()))
    for name, (path, assembler, corset, receipt) in assemblies.items():
        stage = prepare(path, name, output / name, codes.get(name, 1), cpus, assembler, corset, input_proof=proof)
        manifest, stage_receipt = stage / 'inputs.tsv', stage / 'receipt.json'
        proof.read([stage_receipt])
        publication = json.loads(stage_receipt.read_text())
        expected = {str(stage / member): publication['files'][member] for member in
                    ('inputs.tsv', 'all_isoforms.cds.fa', 'transcript_gene_map.tsv')}
        proof.read(expected, expected=expected)
        imported = read_sources(manifest)
        expected_source = {'fasta': str(stage / 'all_isoforms.cds.fa'),
                           'transcript_gene_map': str(stage / 'transcript_gene_map.tsv'), 'genetic_code': codes.get(name, 1)}
        if imported != {name: expected_source}:
            raise ValueError(f'Prepared assembly manifest differs from its verified outputs: {name}')
        files.extend(str(p.resolve()) for p in (path, corset, receipt, manifest, stage_receipt) if p)
        sources.update(imported)
    proof.check()
    return sources, files
