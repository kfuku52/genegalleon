"""Hash-verified whole-contig checkpoints outside final model publication."""
import hashlib
import json
import re
import shutil
from pathlib import Path

from input_generation_array_state import atomic_json, digest


def token(key):
    return hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()


def chunk_key(identity, group):
    return dict(identity=identity, contigs=group)


def search_plan(root, identity, contigs, resources):
    from cloned_cds_resources import groups_for
    from rescue_gene_models import verified
    path = root / 'plan.json'
    if path.exists():
        old = json.loads(path.read_text())
        if old.get('identity') != identity or [x for g in old['groups'] for x in g] != contigs:
            raise ValueError('Search checkpoint plan does not match the complete reference')
        groups = []
        for group in old['groups']:
            key = chunk_key(identity, group)
            # A smaller allocation refines unfinished groups, retaining verified
            # results even when their old group exceeds the new memory budget.
            if verified(root / 'parts' / token(key), key):
                groups.append(group)
            else:
                groups.extend(groups_for(group, resources))
    else:
        groups = groups_for(contigs, resources)
    state = dict(schema_version=1, identity=identity, groups=groups, resources=resources)
    atomic_json(root / 'plans' / (token(state) + '.json'), state, immutable=True)
    atomic_json(path, state)
    return groups


def run_chunk(genome, queries, root, identity, group, code, cpus, guard):
    from rescue_gene_models import run, stage
    key = chunk_key(identity, group)

    def build(folder):
        source = folder / 'genome.fa'
        try:
            with source.open('w') as handle:
                for seqid, length in group:
                    handle.write('>' + seqid + '\n')
                    for start in range(0, length, 1_000_000):
                        handle.write(genome.fetch(seqid, start, min(length, start + 1_000_000)) + '\n')
            run(['miniprot', '-u', '-T', code, '-t', cpus, '-N', 100, '--outs=0', '--gff', source, queries],
                folder, 'miniprot', folder / 'models.gff')
            log = (folder / 'logs/miniprot.log').read_text()
            peak = re.search(r'Peak RSS: ([0-9.]+) GB', log)
            atomic_json(folder / 'metrics.json', dict(bases=sum(x[1] for x in group),
                        peak_rss_bytes=int(float(peak[1]) * 1024 ** 3) if peak else None))
        finally:
            source.unlink(missing_ok=True)
    return stage(root, Path('parts') / token(key), key, build, guard=guard, publication_guard=guard)


def copy_result(part, target):
    """Final publication has its own verified copy; checkpoint files stay reusable."""
    source = part / 'models.gff'
    expected = digest(source)
    shutil.copyfile(source, target)
    if digest(target) != expected or digest(source) != expected:
        raise ValueError('Search checkpoint changed while copying results')
