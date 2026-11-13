"""Infer conservative sequence-only correspondences for assembly loci."""
import hashlib
import json
import tempfile
from itertools import combinations
from pathlib import Path

try:
    from gene_model_selection import _baseline
    from gene_model_store import iter_loci
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest
    from model_protein_correspondence import reciprocal_matches
    from rescue_gene_models import run, stage
except ImportError:
    from .gene_model_selection import _baseline
    from .gene_model_store import iter_loci
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest
    from .model_protein_correspondence import reciprocal_matches
    from .rescue_gene_models import run, stage


def infer_edges(db, species, transcriptomes, executable, params, directory, cpus):
    pairs = [(a, b) for a, b in combinations(species, 2) if a in transcriptomes or b in transcriptomes]
    if not pairs:
        return []
    lengths, fasta, databases = {}, {}, {}
    for index, name in enumerate(species):
        path = directory / f'baseline_{index}.fa'
        lengths[name] = _write_proteins(db, name, path)
        fasta[name] = path
    edges = []
    identity = params['minimum_identity']
    for i, (a, b) in enumerate(pairs):
        if not lengths[a] or not lengths[b]:
            continue
        if b not in databases:
            database = directory / (fasta[b].stem + '.dmnd')
            run([executable, 'makedb', '--in', fasta[b], '--db', database, '--threads', cpus],
                directory, fasta[b].stem)
            databases[b] = database
        alignments = directory / 'protein_matches.tsv'
        run([executable, 'blastp', '--query', fasta[a], '--db', databases[b], '--threads', cpus,
             '--sensitive', '--evalue', '1e-5', '--max-target-seqs', '0', '--max-hsps', '1',
             '--outfmt', '6'], directory, f'protein_pair_{i}', alignments)
        for ga, gb, metrics in reciprocal_matches(alignments, lengths[a], lengths[b],
                                                  min_identity=identity, coverage_a=.80, coverage_b=.80):
            edges.append(_edge(a, b, ga, gb, metrics))
        alignments.unlink()
    for path in (*fasta.values(), *databases.values()):
        path.unlink()
    return edges


def fixed_panel_edges(db, target, donors, executable, params, directory, cpus, cache, catalog_keys, guard):
    """Reuse verified protein DBs and comparisons without cohort-wide searches."""
    target_fasta = directory / 'target.fa'
    target_lengths = _write_proteins(db, target, target_fasta)
    if not target_lengths:
        return []
    target_sha = digest(target_fasta)
    tools = {'diamond': digest(executable), 'extractor': digest(Path(__file__)),
             'reciprocal_matches': digest(Path(reciprocal_matches.__code__.co_filename))}
    edges = []
    for donor in donors:
        donor_key = {'catalog': catalog_keys[donor], 'tools': tools}

        def build_database(tmp, species=donor):
            fasta = tmp / 'proteins.fa'
            lengths = _write_proteins(db, species, fasta)
            atomic_json(tmp / 'lengths.json', lengths)
            if lengths:
                run([executable, 'makedb', '--in', fasta, '--db', tmp / 'proteins.dmnd',
                     '--threads', cpus], tmp, 'makedb')

        database = stage(cache, Path('diamond') / _key(donor_key), donor_key, build_database, guard,
                         publication_guard=guard)
        database_hash, database_files, database_proof = _publication(database, donor_key)
        lengths = _read_json(database / 'lengths.json', database_files['lengths.json'])
        database_proof.check()
        if not lengths:
            continue
        pair_key = {'target': target, 'donor': donor, 'target_protein': target_sha,
                    'database_receipt': database_hash,
                    'minimum_identity': params['minimum_identity'], 'tools': tools}

        def pair_guard(database_proof=database_proof):
            guard()
            database_proof.check()

        def build_pair(tmp, species=donor, database=database, lengths=lengths):
            with tempfile.TemporaryDirectory(prefix='diamond-', dir=directory) as scratch:
                scratch = Path(scratch)
                alignments = scratch / 'matches.tsv'
                run([executable, 'blastp', '--query', target_fasta, '--db', database / 'proteins.dmnd',
                     '--threads', cpus, '--sensitive', '--evalue', '1e-5',
                     '--max-target-seqs', '0', '--max-hsps', '1', '--outfmt', '6'],
                    scratch, 'blastp', alignments)
                rows = [_edge(target, species, ga, gb, metrics)
                        for ga, gb, metrics in reciprocal_matches(alignments, target_lengths, lengths,
                            min_identity=params['minimum_identity'], coverage_a=.80, coverage_b=.80)]
            atomic_json(tmp / 'edges.json', rows)

        pair = stage(cache, Path('protein_pairs') / _key(pair_key), pair_key, build_pair, pair_guard,
                     publication_guard=pair_guard)
        _, pair_files, pair_proof = _publication(pair, pair_key)
        edges.extend(_read_json(pair / 'edges.json', pair_files['edges.json']))
        pair_proof.check()
        pair_guard()
    return edges


def _publication(directory, key):
    receipt = directory / 'receipt.json'
    raw = receipt.read_bytes()
    expected = hashlib.sha256(raw).hexdigest()
    publication = json.loads(raw)
    if publication['key'] != key:
        raise ValueError(f'Fixed donor comparison cache changed: {receipt}')
    proof = FreshDigestBatch()
    members = {str(directory / name): sha for name, sha in publication['files'].items()}
    members[str(receipt)] = expected
    proof.read(members, expected=members)
    return expected, publication['files'], proof


def _read_json(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f'Fixed donor comparison cache changed: {path}')
    return json.loads(raw)


def _key(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _write_proteins(db, species, path):
    lengths = {}
    with path.open('w') as out:
        for locus in iter_loci(db, species):
            if locus.get('ambiguous_coordinates'):
                continue
            candidate = _baseline(locus, locus['candidates'])
            protein = candidate['protein']
            if protein and candidate['quality']['usable']:
                gene = locus['gene_id']
                out.write(f'>{gene}\n{protein}\n')
                lengths[gene] = len(protein)
    return lengths


def _edge(a, b, ga, gb, metrics):
    return {'species_a': a, 'gene_a': ga, 'species_b': b, 'gene_b': gb,
            'weight': 1., 'ambiguous': False, 'comparison_basis': 'protein_alignment',
            'evidence': {'kind': 'unique_reciprocal_baseline_protein_match', **metrics}}
