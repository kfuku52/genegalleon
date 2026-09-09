"""Exercise neighborhood extraction and real DIAMOND searches together."""
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'


@pytest.mark.parametrize('mode', ['cds', 'protein'])
def test_search_windows_strands_boundaries_and_cache(tmp_path, mode):
    for tool in ['diamond', 'seqkit']:
        assert shutil.which(tool), f'{tool} is required; run in the GeneGalleon container'
    rng = random.Random(314)
    codons = dict(zip('ACDEFGHIKLMNPQRSTVWY',
                      ['GCT','TGT','GAT','GAA','TTT','GGT','CAT','ATT','AAA','CTG',
                       'ATG','AAT','CCT','CAA','CGT','TCT','ACT','GTT','TGG','TAT']))
    proteins = [''.join(rng.choices(list(codons), k=100)) for _ in range(45)]
    sequences = proteins if mode == 'protein' else [''.join(codons[a] for a in p) for p in proteins]
    focal = tmp_path / 'focal.fa'
    focal.write_text(''.join(f'>{sp}_g{i}\n{sequences[i]}\n'
                            for sp in ['Species_a', 'Species_b'] for i in [0, 22, 23]))
    for sp, strand in [('Species_a', '+'), ('Species_b', '-')]:
        (tmp_path / f'{sp}.fasta').write_text(''.join(
            f'>{sp}_g{i}\n{seq}\n' for i, seq in enumerate(sequences)))
        rows = []
        for i in range(45):
            # Last two genes lie on a separate chromosome: never cross its boundary.
            chrom = 'chr1' if i < 43 else 'chr2'
            start, end = i * 400 + 1, i * 400 + 300
            for feature, attrs in [('gene', f'ID=g{i}'),
                                   ('mRNA', f'ID=t{i};Parent=g{i}'),
                                   ('CDS', f'ID=c{i};Parent=t{i}')]:
                rows.append(f'{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n')
        (tmp_path / f'{sp}.gff').write_text(''.join(rows))
    cache = tmp_path / 'cache'
    for window in [5, 20]:
        out = tmp_path / f'window{window}.tsv'
        subprocess.run([sys.executable, str(SUPPORT / 'synteny_neighbors.py'),
                        '--focal_cds_fasta', str(focal), '--dir_sp_cds', str(tmp_path),
                        '--dir_sp_gff', str(tmp_path), '--cache_dir', str(cache),
                        '--lock_dir', str(tmp_path / 'locks'),
                        '--gff2genestat_script', str(SUPPORT / 'gff2genestat.py'),
                        '--input_sequence_mode', mode, '--window', str(window),
                        '--evalue', '1e-10', '--threads', '2', '--outfile', str(out)],
                       check=True, capture_output=True, text=True)
        df = pd.read_csv(out, sep='\t')
        expected = set()
        for sp, sign in [('Species_a', 1), ('Species_b', -1)]:
            for center in [0, 22, 23]:
                for offset in range(-window, window + 1):
                    i = center + offset * sign
                    if offset and 0 <= i < 43:
                        expected.add((f'{sp}_g{center}', offset, f'{sp}_g{i}'))
        assert set(df[['node_name', 'offset', 'neighbor_gene']].itertuples(index=False, name=None)) == expected
        # Overlapping windows must reuse a gene's group; identical proteins from
        # the two species must share a group, including the outermost neighbors.
        assert df.groupby('neighbor_gene').group_id.nunique().eq(1).all()
        df['index'] = df.neighbor_gene.str.extract(r'_g(\d+)$').astype(int)
        assert df.groupby('index').group_id.nunique().eq(1).all()
        assert df.group_id.nunique() == df['index'].nunique()
        stamps = {p.name: p.stat().st_mtime_ns for p in cache.glob('*.tsv')}
        if window == 5:
            initial_stamps = stamps
        else:
            assert stamps == initial_stamps  # Search changes reuse the GFF cache.
