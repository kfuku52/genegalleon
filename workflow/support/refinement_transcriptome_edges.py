"""Infer conservative sequence-only correspondences for assembly loci."""
from itertools import combinations

try:
    from gene_model_selection import _baseline
    from gene_model_store import iter_loci
    from model_protein_correspondence import reciprocal_matches
    from rescue_gene_models import run
except ImportError:
    from .gene_model_selection import _baseline
    from .gene_model_store import iter_loci
    from .model_protein_correspondence import reciprocal_matches
    from .rescue_gene_models import run


def infer_edges(db, species, transcriptomes, executable, params, directory, cpus):
    pairs = [(a, b) for a, b in combinations(species, 2) if a in transcriptomes or b in transcriptomes]
    if not pairs:
        return []
    lengths, fasta, databases = {}, {}, {}
    for index, name in enumerate(species):
        lengths[name] = {}
        path = directory / f'baseline_{index}.fa'
        with path.open('w') as out:
            for locus in iter_loci(db, name):
                if locus.get('ambiguous_coordinates'):
                    continue
                candidate = _baseline(locus, locus['candidates'])
                protein = candidate['protein']
                if protein and candidate['quality']['usable']:
                    out.write(f">{locus['gene_id']}\n{protein}\n")
                    lengths[name][locus['gene_id']] = len(protein)
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
             '--query-cover', '80', '--subject-cover', '80', '--id', str(100 * identity),
             '--outfmt', '6'], directory, f'protein_pair_{i}', alignments)
        for ga, gb, metrics in reciprocal_matches(alignments, lengths[a], lengths[b],
                                                  min_identity=identity, coverage_a=.80, coverage_b=.80):
            edges.append({'species_a': a, 'gene_a': ga, 'species_b': b, 'gene_b': gb,
                          'weight': 1., 'ambiguous': False, 'comparison_basis': 'protein_alignment',
                          'evidence': {'kind': 'unique_reciprocal_baseline_protein_match', **metrics}})
        alignments.unlink()
    for path in (*fasta.values(), *databases.values()):
        path.unlink()
    return edges
