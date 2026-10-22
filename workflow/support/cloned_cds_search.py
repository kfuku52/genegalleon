"""Memory-aware, resumable miniprot searches without cutting genomic contigs."""
import json
import shutil
from pathlib import Path

from input_generation_array_state import FreshDigestBatch, atomic_json, digest

DEPENDENCIES = ('cloned_cds_search.py', 'cloned_cds_resources.py', 'cloned_cds_checkpoints.py')


def implementation():
    from cloned_cds_checkpoints import token
    return token({name: digest(Path(__file__).with_name(name)) for name in DEPENDENCIES})


def alignments(genome, queries, directory, code, cpus, batch_bases='auto', *, checkpoint_root=None,
               identity=None, guard=None, verification=None):
    from cloned_cds_checkpoints import copy_result, run_chunk, search_plan, token
    from cloned_cds_resources import memory_plan
    from fasta_sequence_store import exclusive_lock
    directory = Path(directory)
    resources = memory_plan(cpus, batch_bases)
    contigs = [[name, length] for name, length in zip(genome.references, genome.lengths, strict=True)]
    tool = shutil.which('miniprot')
    if not tool:
        raise ValueError('miniprot is required for cloned CDS search')
    inputs = verification if verification is not None else FreshDigestBatch()
    paths = dict(genome=Path(genome.filename.decode()).resolve(), queries=Path(queries), tool=Path(tool))
    hashes = inputs.read(paths.values())
    identity = dict(source=identity, inputs={role: hashes[str(path)] for role, path in paths.items()},
                    implementation=implementation(), code=int(code), cpus=int(cpus), batch_bases=str(batch_bases),
                    parameters=['-u', '-N', '100', '--outs=0', '--gff'])
    root = Path(checkpoint_root) if checkpoint_root is not None else directory / 'search_batches' / token(identity)

    def check():
        inputs.check()
        if guard:
            guard()
    with exclusive_lock(root / '.plan.lock'):
        groups = search_plan(root, identity, contigs, resources)
        batches = []
        for index, group in enumerate(groups, 1):
            part = run_chunk(genome, queries, root, identity, group, code, cpus, check)
            target = directory / ('models.gff' if len(groups) == 1 else f'models.{index}.gff')
            copy_result(part, target)
            batches.append(dict(contigs=[r[0] for r in group], bases=sum(r[1] for r in group),
                                checkpoint=str(part), receipt_sha256=digest(part / 'receipt.json'),
                                metrics=json.loads((part / 'metrics.json').read_text())))
            atomic_json(directory / 'index_batches.json', dict(schema_version=2, batch_bases=resources['target_bases'],
                        contigs_split=False, resources=resources, complete=False, batches=batches))
            yield target
        check()
        atomic_json(directory / 'index_batches.json', dict(schema_version=2, batch_bases=resources['target_bases'],
                    contigs_split=False, resources=resources, complete=True, batches=batches))
