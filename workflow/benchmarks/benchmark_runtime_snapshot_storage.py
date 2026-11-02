#!/usr/bin/env python3
"""Compare complete raw/runtime-SQLite inventories over identical captured sources."""
import argparse
import hashlib
import json
import os
import resource
import statistics
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
import metadata_store as metadata  # noqa: E402
import runtime_snapshot_storage as runtime  # noqa: E402


def fixture(root, count=2000, versions=24):
    pool = root / '.shared-code-v1'
    pool.mkdir(parents=True)
    bodies = []
    for version in range(versions):
        selected = {}
        for name in ['gg_gene_evolution_core.sh','gg_gene_evolution_entrypoint.sh','gg_common_params.sh','VERSION']:
            prefix = (str(version)+':'+name+'\n').encode()
            size = 262144 if name.endswith('_core.sh') else 16384 if name.endswith('.sh') else 10
            raw = (prefix*(size//len(prefix)+1))[:size]
            mode = 0o444 if name=='VERSION' else 0o555
            path = pool/(hashlib.sha256(raw).hexdigest()+'.'+format(mode,'o'))
            path.write_bytes(raw)
            path.chmod(mode)
            os.utime(path,ns=(1600000000123456789,1600000000123456789))
            selected[name] = path
        bodies.append(selected)
    for index in range(count):
        directory = root/'gg_gene_evolution_entrypoint'/('130000000_'+str(index+1))
        directory.mkdir(parents=True)
        for name,path in bodies[index%versions].items():
            os.link(path,directory/name)


def raw_rows(root):
    rows, seen = [], {}
    for parent,_dirs,names in os.walk(root,followlinks=False):
        for name in sorted(names):
            path = Path(parent)/name
            info = path.lstat()
            key = info.st_dev,info.st_ino
            if key not in seen:
                seen[key] = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append([path.relative_to(root).as_posix(),info.st_mode&0o777,info.st_size,info.st_mtime_ns,seen[key]])
    return rows


def archived_rows(root):
    database = root/metadata.DATABASE
    before = metadata._file_sha(database)
    connection = metadata._connection(root)[0]
    runtime.verify_database(connection)
    rows = [list(r) for r in connection.execute('SELECT path,mode,bytes,mtime_ns,sha256 FROM files')]
    rows.extend([[str(Path('.shared-code-v1')/p.name),p.stat().st_mode&0o777,p.stat().st_size,
                  p.stat().st_mtime_ns,hashlib.sha256(p.read_bytes()).hexdigest()]
                 for p in (root/'.shared-code-v1').iterdir()])
    if metadata._file_sha(database) != before:
        raise RuntimeError('Runtime database changed during inventory')
    return rows


def inventory(root):
    files = directories = 0
    for _parent,dirs,names in os.walk(root):
        directories += len(dirs)
        files += len(names)
    return dict(files=files,directories=directories)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshots',type=int,default=2000)
    parser.add_argument('--versions',type=int,default=24)
    parser.add_argument('--repetitions',type=int,default=3)
    args = parser.parse_args()
    if not 1<=args.snapshots<=25000 or not 1<=args.versions<=128 or not 1<=args.repetitions<=10:
        parser.error('bounded workload required')
    with tempfile.TemporaryDirectory(prefix='runtime-storage-benchmark-')as folder:
        base = Path(folder)
        before,after = base/'before/runtime',base/'after/runtime'
        fixture(before,args.snapshots,args.versions)
        fixture(after,args.snapshots,args.versions)
        collected = metadata.pack(after,writers_stopped=True,delete_originals=True,operation_id='benchmark')
        removed = runtime.prune(after,writers_stopped=True)
        samples = dict(before=[],after=[])
        expected = None
        for repetition in range(args.repetitions+1):
            for key in (['before','after']if repetition%2==0 else ['after','before']):
                started = time.perf_counter()
                rows = raw_rows(before)if key=='before'else archived_rows(after)
                proof = hashlib.sha256(json.dumps(sorted(rows),separators=(',',':')).encode()).hexdigest()
                elapsed = time.perf_counter()-started
                if expected is None:
                    expected = proof
                if proof!=expected:
                    raise RuntimeError('Complete captured bytes, names or metadata differ')
                if repetition:
                    samples[key].append(elapsed)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin'else 1024)
        print(json.dumps(dict(python=sys.version,platform=sys.platform,machine=os.uname().machine,
            snapshots=args.snapshots,versions=args.versions,warmups=1,repetitions=args.repetitions,
            complete_output_identical=True,proof=expected,collection=collected,removed_directories=removed,
            physical_before=inventory(before),physical_after=inventory(after),
            database_bytes=(after/metadata.DATABASE).stat().st_size,peak_rss_bytes=rss,
            seconds=samples,median_seconds={k:statistics.median(v)for k,v in samples.items()})))


if __name__=='__main__':
    main()
