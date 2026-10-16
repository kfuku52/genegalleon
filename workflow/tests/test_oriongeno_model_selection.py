"""Taxonomy-backed OrionGeno checkpoints, including remainder-group precedence."""
import argparse
import json
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'support'))
import annevo_model_selection as selection
import input_generation_array_state as state
import oriongeno_model_profiles as profiles
from annevo_runtime import selection_record
from annevo_runtime_cache import RuntimeVerification

EXAMPLES = [('plants', 3702, 'Arabidopsis thaliana'), ('mammals', 9606, 'Homo sapiens'),
            ('birds', 9031, 'Gallus gallus'), ('fish', 7955, 'Danio rerio'),
            ('other_vertebrates', 8364, 'Xenopus tropicalis'), ('arthropods', 7227, 'Drosophila melanogaster'),
            ('other_invertebrates', 6239, 'Caenorhabditis elegans'), ('fungi', 4932, 'Saccharomyces cerevisiae')]


def taxonomy(path):
    with sqlite3.connect(path) as conn:
        conn.executescript('CREATE TABLE species(taxid INTEGER PRIMARY KEY,spname TEXT,rank TEXT,track TEXT);'
                           'CREATE TABLE synonym(taxid INTEGER,spname TEXT);'
                           'CREATE TABLE merged(taxid_old INTEGER,taxid_new INTEGER);')
        conn.execute("INSERT INTO species VALUES(1,'root','no rank','1')")
        tracks = {}
        # Explicit fixture tracks test the nested broad models independently
        # of the selector's own ancestor-profile rules.
        chains = {'plants': [33090, 1], 'mammals': [40674, 7742, 33208, 1],
                  'birds': [8782, 7742, 33208, 1], 'fish': [7898, 7742, 33208, 1],
                  'other_vertebrates': [7742, 33208, 1], 'arthropods': [6656, 33208, 1],
                  'other_invertebrates': [33208, 1], 'fungi': [4751, 1]}
        for profile, (taxid, name) in profiles.ANCHORS.items():
            tracks[profile] = ','.join(map(str, chains[profile]))
            conn.execute('INSERT INTO species VALUES(?,?,?,?)', (taxid, name, 'clade', tracks[profile]))
        for profile, taxid, name in EXAMPLES:
            conn.execute('INSERT INTO species VALUES(?,?,?,?)', (taxid, name, 'species', f'{taxid},'+tracks[profile]))
        for taxid, name, track in [(3218, 'Physcomitrium patens', '3218,33090,1'),
                                   (6945, 'Ixodes scapularis', '6945,6656,33208,1'),
                                   (7868, 'Callorhinchus milii', '7868,7742,33208,1'),
                                   (5833, 'Plasmodium falciparum', '5833,1')]:
            conn.execute('INSERT INTO species VALUES(?,?,?,?)', (taxid, name, 'species', track))
        conn.execute('INSERT INTO merged VALUES(999999,9606)')


def fixture(tmp_path, taxid=9606):
    db, genome, mapping = [tmp_path/name for name in ('taxa.sqlite', 'genome.fa', 'taxids.tsv')]
    taxonomy(db)
    genome.write_text('>chr\nATGAAATAA\n')
    mapping.write_text(f'species\ttaxid\nTest_species\t{taxid}\n')
    return argparse.Namespace(engine='oriongeno', species='Test_species', genome=str(genome), model_profile='auto',
                              taxonomy_db=str(db), taxid_map=str(mapping), genome_metadata='',
                              output=str(tmp_path/'selection.json'))


@pytest.mark.parametrize('profile,taxid', [(p, t) for p, t, _ in EXAMPLES]+
                         [('plants', 3218), ('arthropods', 6945), ('other_vertebrates', 7868)])
def test_most_specific_supported_ancestor_selects_each_checkpoint(tmp_path, profile, taxid):
    args = fixture(tmp_path, taxid)
    record = selection.select(args)
    assert record['kind'] == 'genegalleon-oriongeno-model-selection'
    assert record['model_profile'] == profile
    assert record['anchor_taxid'] == profiles.ANCHORS[profile][0]
    assert record['profiles_sha256'] == selection.file_sha(profiles.__file__)
    assert record['taxonomy']['sha256'] == selection.file_sha(args.taxonomy_db)
    assert record['evidence'][0]['sha256'] == selection.file_sha(args.taxid_map)
    before = Path(args.output).read_bytes(), Path(args.output).stat().st_mtime_ns
    selection.select(args)
    assert (Path(args.output).read_bytes(), Path(args.output).stat().st_mtime_ns) == before
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    assert selection_record(native_args, args.genome, engine='oriongeno', model_profile=profile)[0] == record


@pytest.mark.parametrize('profile', sorted(profiles.PROFILES))
def test_explicit_profile_requires_no_taxonomy(tmp_path, profile):
    genome = tmp_path/'genome.fa'
    genome.write_text('>s\nATGAAATAA\n')
    args = argparse.Namespace(engine='oriongeno', species='Test_species', genome=str(genome), model_profile=profile,
                              taxonomy_db='', taxid_map='', genome_metadata='', output=str(tmp_path/'selected.json'))
    assert selection.select(args)['reason'] == 'explicit_model_profile'


@pytest.mark.parametrize('taxid,message', [(5833, 'No supported'), (123456, 'Unresolved'),
                                          (40674, 'does not identify a species')])
def test_unsupported_taxonomy_retains_previous_selection(tmp_path, taxid, message):
    args = fixture(tmp_path)
    selection.select(args)
    previous = Path(args.output).read_bytes()
    Path(args.taxid_map).write_text(f'species\ttaxid\nTest_species\t{taxid}\n')
    with pytest.raises(ValueError, match=message):
        selection.select(args)
    assert Path(args.output).read_bytes() == previous
    assert json.loads(Path(args.output+'.failed.json').read_text())['status'] == 'error'


def test_merged_taxid_and_hash_bound_metadata_agree(tmp_path):
    args = fixture(tmp_path, 999999)
    metadata = tmp_path/'metadata.tsv'
    metadata.write_text('species\ttaxid\tgenome_sha256\tscientific_name\nTest_species\t9606\t'+
                        selection.file_sha(args.genome)+'\tHomo sapiens\n')
    args.genome_metadata = str(metadata)
    record = selection.select(args)
    assert record['taxid'] == 9606 and record['merged_taxids'] == [999999]
    Path(args.taxid_map).unlink()
    args.taxid_map = ''
    assert selection.select(args)['model_profile'] == 'mammals'
    metadata.write_text(metadata.read_text().replace('9606', '9031'))
    with pytest.raises(ValueError, match='scientific_name'):
        selection.select(args)


@pytest.mark.parametrize('changed', ['species', 'model', 'genome', 'mapping', 'taxonomy'])
def test_native_boundary_rechecks_selection_inputs(tmp_path, changed):
    args = fixture(tmp_path)
    selection.select(args)
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    profile = 'mammals'
    if changed == 'species':
        native_args.species = 'Other_species'
    elif changed == 'model':
        profile = 'birds'
    elif changed == 'genome':
        Path(args.genome).write_text('>changed\nATGAAATAA\n')
    elif changed == 'mapping':
        Path(args.taxid_map).write_text('species\ttaxid\nTest_species\t9031\n')
    else:
        with sqlite3.connect(args.taxonomy_db) as conn:
            conn.execute("INSERT INTO synonym VALUES(9606,'late synonym')")
    with pytest.raises(ValueError, match='differs|changed'):
        selection_record(native_args, args.genome, engine='oriongeno', model_profile=profile)


def test_unrelated_clades_cannot_be_treated_as_nested_remainder_groups(tmp_path):
    args = fixture(tmp_path)
    with sqlite3.connect(args.taxonomy_db) as conn:
        conn.execute("UPDATE species SET track='9606,40674,4751,1' WHERE taxid=9606")
        conn.execute("UPDATE species SET track='40674,4751,1' WHERE taxid=40674")
    with pytest.raises(ValueError, match='Conflicting OrionGeno'):
        selection.select(args)


def test_shared_selection_hashes_once_and_standalone_calls_remain_fresh(tmp_path, monkeypatch):
    args = fixture(tmp_path)
    selection.select(args)
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    options = dict(engine='oriongeno', model_profile='mammals')
    reads = []
    real_digest = state.digest
    def digest(path):
        reads.append(Path(path).resolve())
        return real_digest(path)
    monkeypatch.setattr(state, 'digest', digest)
    proof = RuntimeVerification()
    first = selection_record(native_args, args.genome, verification=proof, **options)
    assert selection_record(native_args, args.genome, verification=proof, **options) == first
    expected = {Path(p).resolve() for p in (args.output, args.genome, args.taxid_map, args.taxonomy_db)}
    assert set(reads) == expected and len(reads) == len(expected)
    assert selection_record(native_args, args.genome, **options) == first
    assert len(reads) == 2*len(expected)


@pytest.mark.parametrize('changed', ['mapping', 'genome', 'taxonomy', 'selection'])
def test_shared_selection_rejects_same_byte_replacements_and_keeps_failure(tmp_path, changed):
    args = fixture(tmp_path)
    selection.select(args)
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    proof = RuntimeVerification()
    options = dict(engine='oriongeno', model_profile='mammals', verification=proof)
    selection_record(native_args, args.genome, **options)
    target = Path({'mapping': args.taxid_map, 'genome': args.genome,
                   'taxonomy': args.taxonomy_db, 'selection': args.output}[changed])
    replacement = target.with_name(target.name+'.replacement')
    replacement.write_bytes(target.read_bytes())
    replacement.replace(target)
    with pytest.raises(OSError, match='changed'):
        selection_record(native_args, args.genome, **options)
    with pytest.raises(RuntimeError, match='previously failed'):
        proof.check()


def test_coarse_timestamp_selection_checks_fresh_bytes_at_every_boundary(tmp_path, monkeypatch):
    args = fixture(tmp_path)
    selection.select(args)
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    real_identity, real_signature = state._stat_identity, selection.stable_signature
    monkeypatch.setattr(state, '_stat_identity', lambda info: (*real_identity(info)[:3], 0, 0))
    monkeypatch.setattr(selection, 'stable_signature', lambda path: [*real_signature(path)[:3], 0, 0])
    proof = RuntimeVerification()
    options = dict(engine='oriongeno', model_profile='mammals', verification=proof)
    selection_record(native_args, args.genome, **options)
    mapping = Path(args.taxid_map)
    mapping.write_text(mapping.read_text().replace('9606', '9607'))
    with pytest.raises(ValueError, match='input changed'):
        selection_record(native_args, args.genome, **options)
    with pytest.raises(RuntimeError, match='previously failed'):
        proof.check()


def test_selection_closing_proof_rechecks_sqlite_sidecars(tmp_path, monkeypatch):
    args = fixture(tmp_path)
    selection.select(args)
    native_args = argparse.Namespace(selection_receipt=args.output, species=args.species)
    proof = RuntimeVerification()
    real_check = proof.check
    def check_then_begin_writer():
        real_check()
        Path(args.taxonomy_db+'-wal').write_bytes(b'late writer')
    monkeypatch.setattr(proof, 'check', check_then_begin_writer)
    with pytest.raises(ValueError, match='active transaction'):
        selection_record(native_args, args.genome, engine='oriongeno', model_profile='mammals', verification=proof)
