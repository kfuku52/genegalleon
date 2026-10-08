"""Page checkpoints retain exact input checks without rereading older pages."""

import hashlib
import re
import sqlite3
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))

import focus_hgt_context_annotations as annotation_module  # noqa: E402
from focus_hgt_context import GenomeCoordinates, render_context  # noqa: E402
from focus_hgt_context_annotations import RANKS, ContextAnnotations  # noqa: E402
from focus_hgt_gene_trees import write  # noqa: E402


class RecordedFamilyStore:
    def __init__(self, root):
        self.root, self.reads = root, []

    def open_binary(self, subdir, name):
        self.reads.append(subdir + '/' + name)
        return (self.root / subdir / name).open('rb')


def annotation_input(tmp_path, species=('A', 'B')):
    family_root = tmp_path / 'families'
    (family_root / 'stat_branch').mkdir(parents=True)
    db = tmp_path / 'taxa.sqlite'
    with sqlite3.connect(db) as conn:
        conn.execute('CREATE TABLE species(taxid INTEGER PRIMARY KEY, spname TEXT, rank TEXT, track TEXT)')
        conn.execute('CREATE TABLE merged(taxid_old INTEGER, taxid_new INTEGER)')
        conn.executemany('INSERT INTO species VALUES(?,?,?,?)', [
            (1, 'root', 'no rank', '1'), (2, 'Own kingdom', 'kingdom', '2,1'),
            (3, 'Own species', 'species', '3,2,1')])
    records = []
    for name in species:
        gene, family = name + '_gene', 'OG' + name
        records.append(dict(gene_id=gene, orthogroup=family, besthit_accession='', besthit_organism='',
                            besthit_source='original upstream provenance', **{f'besthit_{rank}': '' for rank in RANKS}))
        leaf = dict(node_name=gene, child1='-999', child2='-999', sprot_best='')
        write(family_root / 'stat_branch' / (family + '_stat.branch.tsv'), list(leaf), [leaf])
        (tmp_path / (name + '_mmseqs2taxonomy.tsv')).write_text(
            gene + '\t3\tspecies\tOwn species\t1\t1\t1\t1\t1;2;3\n')
    path = tmp_path / 'annotations.tsv'
    write(path, list(records[0]), records)
    store = RecordedFamilyStore(family_root)
    annotations = ContextAnnotations(path, store=store, mmseqs2_taxonomy_dir=tmp_path, taxonomy_dbfile=db)
    return annotations, store, path, db


def test_page_verification_does_not_grow_with_unused_previous_sources(tmp_path, monkeypatch):
    annotations, store, path, db = annotation_input(tmp_path, tuple('ABCDEFGHIJ'))
    for name in 'ABCDEFGHIJ':
        current = annotations.get(name + '_gene', species=name)
    read_paths = []
    digest = annotation_module.source_digest

    def recorded_digest(source):
        read_paths.append(source)
        return digest(source)

    monkeypatch.setattr(annotation_module, 'source_digest', recorded_digest)
    store.reads.clear()
    annotations.verify(records=[current])
    assert set(read_paths) == {str(path.resolve()), str(db.resolve()),
                               str((tmp_path / 'J_mmseqs2taxonomy.tsv').resolve())}
    assert len(read_paths) == 3
    assert store.reads == ['stat_branch/OGJ_stat.branch.tsv']
    read_paths.clear()
    store.reads.clear()
    annotations.verify()
    assert len(read_paths) == 12 and len(store.reads) == 10


@pytest.mark.parametrize('changed', ['canonical', 'classification', 'taxonomy', 'family'])
def test_scoped_checkpoint_catches_mutations_even_after_cached_annotation_and_lineage_hits(tmp_path, changed):
    annotations, store, path, db = annotation_input(tmp_path)
    original = annotations.get('A_gene', species='A')
    # Both genes use the identical LCA/lineage cache key.
    annotations.get('B_gene', species='B')
    if changed == 'canonical':
        path.write_bytes(path.read_bytes() + b'\n')
    elif changed == 'classification':
        target = tmp_path / 'A_mmseqs2taxonomy.tsv'
        target.write_bytes(target.read_bytes().replace(b'Own species', b'Changed species'))
    elif changed == 'taxonomy':
        with sqlite3.connect(db) as conn:
            conn.execute('UPDATE species SET spname="Changed kingdom" WHERE taxid=2')
    else:
        target = store.root / 'stat_branch' / 'OGA_stat.branch.tsv'
        target.write_bytes(target.read_bytes() + b'\n')
    cached = annotations.get('A_gene', species='A')
    assert cached == original
    with pytest.raises(ValueError, match='changed during rendering'):
        annotations.verify(records=[cached])


@pytest.mark.parametrize('changed', ['classification', 'family'])
def test_final_full_checkpoint_catches_mutation_of_an_earlier_unused_page(tmp_path, changed):
    annotations, store, _, _ = annotation_input(tmp_path)
    first = annotations.get('A_gene', species='A')
    annotations.verify(records=[first])
    second = annotations.get('B_gene', species='B')
    target = (tmp_path / 'A_mmseqs2taxonomy.tsv' if changed == 'classification'
              else store.root / 'stat_branch' / 'OGA_stat.branch.tsv')
    target.write_bytes(target.read_bytes() + b'\n')
    annotations.verify(records=[second])
    with pytest.raises(ValueError, match='changed during rendering'):
        annotations.verify()


def test_scoped_page_checkpoint_preserves_pdf_and_complete_audits(tmp_path, monkeypatch):
    stat = [dict(branch_id='2', node_name='root', child1='0', child2='1', support_generax_ufboot='45'),
            dict(branch_id='0', node_name='D_gene', child1='-999', child2='-999', support_generax_ufboot=''),
            dict(branch_id='1', node_name='A_gene', child1='-999', child2='-999', support_generax_ufboot='')]
    events = [dict(event_id='OG1:2:1', orthogroup='OG1', gene_tree_branch_id='2', gene_tree_node='root')]
    links = [dict(event_id='OG1:2:1', orthogroup='OG1', gene_id=species + '_gene', gene_species=species,
                  side=side, eligible_for_context='True', host_scaffold_status='measured', host_scaffold_id='s1',
                  host_scaffold_background_class_total_count='20', host_scaffold_background_class_compatible_count='9',
                  host_scaffold_background_class_incompatible_count='1', host_scaffold_background_class_unresolved_count='10',
                  host_scaffold_background_class_classified_fraction='0.5',
                  host_scaffold_background_class_compatible_fraction='0.9')
             for species, side in [('D', 'donor'), ('A', 'recipient')]]
    canonical = tmp_path / 'annotations.tsv'
    records = [dict(gene_id=species + '_gene', orthogroup='OG1', besthit_accession='', besthit_organism='',
                    **{f'besthit_{rank}': '' for rank in RANKS}) for species in ['A', 'D']]
    write(canonical, list(records[0]), records)
    earlier = tmp_path / 'previous-page.tsv'
    earlier.write_text('Earlier page input\n')

    def render(path):
        annotations = ContextAnnotations(canonical)
        annotations.sources[str(earlier)] = hashlib.sha256(earlier.read_bytes()).hexdigest()
        audit = render_context(path, stat, events, links, GenomeCoordinates(''), gene_tree_panel=False,
                               max_genes_per_side=3, annotations=annotations)
        annotations.verify()
        return audit, annotations.display_audit, annotations.sources, annotations.family_sources

    scoped = tmp_path / 'scoped.pdf'
    actual = render(scoped)
    verify = ContextAnnotations.verify
    monkeypatch.setattr(ContextAnnotations, 'verify', lambda self, records=None: verify(self))
    full = tmp_path / 'full.pdf'
    expected = render(full)
    assert actual == expected
    def normalize(path):
        return re.sub(rb'/(CreationDate|ModDate) \([^)]*\)', b'', path.read_bytes())

    assert normalize(scoped) == normalize(full)
