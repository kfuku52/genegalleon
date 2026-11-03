"""Compare native single/batch retirement on synthetic family and metadata stores."""
import argparse
import hashlib
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--support', type=Path, required=True)
p.add_argument('--families', type=int, default=200)
p.add_argument('--method', choices=['single', 'batch'], required=True)
p.add_argument('--scratch', type=Path)
a = p.parse_args()
sys.path.insert(0, str(a.support.resolve()))
from gene_family_output_store import GeneFamilyOutputStore  # noqa: E402
from metadata_store import pack  # noqa: E402

with tempfile.TemporaryDirectory(dir=a.scratch) as tmp:
    root = Path(tmp) / 'orthogroup'
    paths, families = [], {}
    for subdir in ['stat_tree', 'stat_branch', 'mafft', 'artifact_provenance']:
        (root / subdir).mkdir(parents=True)
    for index in range(a.families):
        family = f'OG{index:07}'
        for subdir, suffix in [('stat_tree', '_stat.tree.tsv'), ('stat_branch', '_stat.branch.tsv')]:
            name = f'{subdir}/{family}{suffix}'
            (root / name).write_text('synthetic\t1\n')
            paths.append(name)
            families[name] = family
        (root / 'mafft' / f'{family}_cds.aln.fa').write_text('>synthetic\nATG\n')
        (root / 'artifact_provenance' / f'{family}.summary.json').write_text(
            json.dumps({'synthetic': True, 'payload': 'x' * 4096}))
    pack(root / 'artifact_provenance', scratch=tmp, writers_stopped=True, delete_originals=True)
    store = GeneFamilyOutputStore(root)
    started = time.perf_counter()
    if a.method == 'single':
        for name in paths:
            store.delete(name, family_id=families[name])
    else:
        store.delete_many(paths, family_ids=families)
    seconds = time.perf_counter() - started
    store = GeneFamilyOutputStore(root)
    assert not store.file_names('stat_tree') and not store.file_names('stat_branch')
    retained = {}
    for subdir in ['mafft', 'artifact_provenance']:
        for name in store.file_names(subdir):
            with store.open_binary(subdir, name) as handle:
                retained[subdir + '/' + name] = hashlib.sha256(handle.read()).hexdigest()
    print(json.dumps({'method': a.method, 'families': a.families, 'deleted': len(paths),
                      'seconds': seconds, 'max_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
                      'retained_files': len(retained), 'retained_manifest_sha256': hashlib.sha256(
                          json.dumps(retained, sort_keys=True, separators=(',', ':')).encode()).hexdigest()}))
