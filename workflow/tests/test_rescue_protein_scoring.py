"""Check exact native score ordering and interrupted score checkpoint recovery."""
import sqlite3
import sys
from itertools import product
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))

import gene_model_selection as selection
import rescue_protein_scoring as scoring


def write_fastas(directory, queries, donors):
    for name, proteins in [('native.pep', queries), ('donors.pep', donors)]:
        (directory / name).write_text(''.join(f'>{token}\n{protein}\n' for token, protein in proteins.items()))


def test_packed_worker_scores_match_pair_score(tmp_path, monkeypatch):
    queries = {'n0': 'ACDJBZXY', 'n1': 'M' * 5001, 'n2': 'MALLLLLLLLLADTW*'}
    donors = {'d0': 'acdjbzxy*', 'd1': 'ACDDATY', 'd2': 'MATW', 'd3': 'M' * 5000, 'd4': '*'}
    write_fastas(tmp_path, queries, donors)
    # Force the spawned worker path with several long but bounded comparisons.
    pairs = list(product(queries, donors)) * 35
    reference = {(q, d): selection.pair_score({'protein': queries[q]}, {'protein': donors[d]}) for q, d in pairs[:15]}
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
