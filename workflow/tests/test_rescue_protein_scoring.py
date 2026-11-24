"""Check exact native score ordering and interrupted score checkpoint recovery."""
import random
import sqlite3
import sys
from itertools import product
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))

import gene_model_selection as selection
import rescue_alignment_kernel as alignment
import rescue_protein_scoring as scoring


def write_fastas(directory, queries, donors):
    for name, proteins in [('native.pep', queries), ('donors.pep', donors)]:
        (directory / name).write_text(''.join(f'>{token}\n{protein}\n' for token, protein in proteins.items()))


@pytest.mark.parametrize('use_kernel', [False, True])
def test_packed_worker_scores_match_pair_score(tmp_path, monkeypatch, use_kernel):
    if not use_kernel:
        monkeypatch.setattr(alignment, 'compiled_kernel', lambda: None)
    elif alignment.compiled_kernel() is None:
        pytest.skip('GeneGalleon SIMD alignment library is unavailable')
    prefix = 'ANQH' * 7
    queries = {'n0': 'ACDJBZXY', 'n1': 'M' * 5001, 'n2': prefix + 'MALLLLLLLLLADTW*'}
    donors = {'d0': 'acdjbzxy*', 'd1': 'ACDDATY', 'd2': prefix + 'MATW', 'd3': 'M' * 5000, 'd4': '*'}
    write_fastas(tmp_path, queries, donors)
    # Exercise both worker pools, two groups and restoration of input order.
    pairs = list(product(queries, donors)) * 69
    pairs = [pair for pair in pairs if pair not in {('n1', 'd2'), ('n2', 'd3')}]
    reference = {(q, d): selection.pair_score({'protein': queries[q]}, {'protein': donors[d]})
                 for q, d in set(pairs)}
    monkeypatch.setattr(scoring, 'CHECKPOINT_PAIRS', 512)
    checkpoint = tmp_path / 'scores.sqlite3'
    fresh = list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 2, checkpoint=checkpoint))
    serial = list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 1, checkpoint=checkpoint))
    assert fresh == serial == [(q, d, reference[q, d]) for q, d in pairs]


def test_checkpoint_recovers_committed_groups_and_rejects_changed_inputs(tmp_path, monkeypatch):
    queries = {'n0': 'MACDW', 'n1': 'MACW', 'n2': 'MRAW'}
    donors = {'d0': 'MACDW', 'd1': 'MYRAW'}
    write_fastas(tmp_path, queries, donors)
    pairs = list(product(queries, donors))
    monkeypatch.setattr(scoring, 'CHECKPOINT_PAIRS', 2)
    checkpoint = tmp_path / 'scores.sqlite3'
    expected = list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 1))
    with monkeypatch.context() as interrupted:
        original = scoring._score_batch
        calls = []

        def stop_second(ids, context):
            calls.append(ids)
            if len(calls) == 2:
                raise RuntimeError('interrupted')
            return original(ids, context)

        interrupted.setattr(scoring, '_score_batch', stop_second)
        with pytest.raises(RuntimeError, match='interrupted'):
            list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 1, checkpoint=checkpoint))
    with sqlite3.connect(checkpoint) as db:
        assert db.execute('SELECT position FROM scores').fetchall() == [(0,)]
    assert list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 3, checkpoint=checkpoint)) == expected
    assert sorted(path.name for path in tmp_path.iterdir()) == ['donors.pep', 'native.pep', 'scores.sqlite3']
    with pytest.raises(ValueError, match='checkpoint inputs changed'):
        list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 1, checkpoint=checkpoint, checkpoint_key='new'))
    with pytest.raises(ValueError, match='checkpoint pair IDs changed'):
        list(scoring.scored_pairs(reversed(pairs), queries, donors, tmp_path, 1, checkpoint=checkpoint))
    with sqlite3.connect(checkpoint) as db:
        db.execute("UPDATE scores SET payload=x'00' WHERE position=0")
    with pytest.raises(ValueError, match='Damaged .* checkpoint'):
        list(scoring.scored_pairs(pairs, queries, donors, tmp_path, 1, checkpoint=checkpoint))


def test_simd_kernel_matches_biopython_ties_and_gap_flags():
    kernel = alignment.compiled_kernel()
    if kernel is None:
        pytest.skip('GeneGalleon SIMD alignment library is unavailable')
    rng = random.Random(3042)
    alphabet = 'ACDEFGHIKLMNPQRSTVWYXBZJ'
    varied = [(''.join(rng.choices(alphabet, k=rng.randrange(1, 200))),
               ''.join(rng.choices(alphabet, k=rng.randrange(1, 200)))) for _ in range(512)]
    cases = [('AAAAA', 'AAAA'), ('AB', 'BA'), ('AXBZJC', 'AXBZJC'),
             ('CEHQA' + 'P' * 9 + 'CMWQV', 'CEHQACMWQV'),
             ('CMWQV', 'P' * 9 + 'CMWQV'), ('Q' * 81, 'Q' * 31),
             ('AC' * 120, 'CAC' * 80), ('TEFGHA', 'EFGHAT')]
    stats = alignment.buffer()
    for first, second in [*cases, *[pair[::-1] for pair in cases], *varied]:
        n, m = len(first), len(second)
        assert kernel(first.encode(), n, second.encode(), m, stats) == 0
        assert alignment.score_counts(n, m, *stats) == selection.pair_score(
            {'protein': first}, {'protein': second}), (first, second)
