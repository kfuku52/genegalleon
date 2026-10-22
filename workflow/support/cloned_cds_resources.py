"""Conservative miniprot index planning from job/container limits, never host RAM."""
import os
from decimal import Decimal, InvalidOperation
from pathlib import Path

GIB = 1024 ** 3
BYTES_PER_BASE = 16
FALLBACK_BASES = 250_000_000


def cgroup_limits(root=Path('/sys/fs/cgroup'), membership=Path('/proc/self/cgroup')):
    """Read effective v1/v2 limits, including constrained ancestors and bind roots."""
    try:
        rows = membership.read_text().splitlines()
    except OSError:
        return {}
    limits = {}
    for row in rows:
        fields = row.split(':', 2)
        if len(fields) != 3:
            continue
        if fields[1] == '':
            base, filename = root, 'memory.max'
        elif 'memory' in fields[1].split(','):
            base, filename = root / 'memory', 'memory.limit_in_bytes'
        else:
            continue
        relative = Path(fields[2].lstrip('/'))
        # Namespace-relative paths may refer outside the visible mount. Its
        # visible root still constrains us; never traverse outside that root.
        parts = [] if '..' in relative.parts else list(relative.parts)
        for n in range(len(parts) + 1):
            path = base.joinpath(*parts[:n], filename)
            try:
                raw = path.read_text().strip()
            except FileNotFoundError:
                continue
            if raw == 'max':
                continue
            if not raw.isdecimal() or int(raw) < 1:
                raise ValueError('Invalid cgroup memory limit: ' + str(path))
            value = int(raw)
            if value < 2 ** 60:
                limits[str(path)] = value
    return limits


def memory_plan(cpus, requested='auto', *, environ=None, limits=None):
    env = os.environ if environ is None else environ
    sources = dict(cgroup_limits() if limits is None else limits)
    cpus = int(cpus)
    if cpus < 1:
        raise ValueError('Search CPU count must be positive')
    for name, unit, multiplier in (
        ('GG_MEM_TOTAL_GB', GIB, 1), ('GG_MEM_TOOL_GB', GIB, 1),
        ('SLURM_MEM_PER_NODE', 1024 ** 2, 1), ('SLURM_MEM_PER_CPU', 1024 ** 2, cpus),
    ):
        if env.get(name) in (None, ''):
            continue
        try:
            value = Decimal(env[name])
            if not value.is_finite() or value < 0 or value == 0 and name.startswith('GG_'):
                raise InvalidOperation
            if value:
                sources[name] = int(value * unit * multiplier)
        except (InvalidOperation, ValueError):
            raise ValueError('Invalid allocated memory: ' + name) from None
    raw = str(requested)
    if raw != 'auto' and (not raw.isdecimal() or int(raw) < 1):
        raise ValueError('cloned_cds_index_batch_bases must be auto or a positive integer')
    allocation = min(sources.values()) if sources else None
    # Reserve 25% for Python/GFF/other processes, plus fixed/per-thread space.
    budget = allocation * 3 // 4 if allocation is not None else None
    overhead = (512 + 64 * cpus) * 1024 ** 2
    possible = (budget - overhead) // BYTES_PER_BASE if budget is not None else FALLBACK_BASES
    if possible < 1:
        raise ValueError('Allocated memory is too small for miniprot index overhead')
    target = possible if raw == 'auto' else min(int(raw), possible)
    return dict(schema_version=1, requested=raw, allocation_bytes=allocation, budget_bytes=budget,
                sources=sources, overhead_bytes=overhead, bytes_per_base=BYTES_PER_BASE,
                target_bases=target, estimate='conservative-heuristic', cpus=cpus)


def groups_for(contigs, plan):
    groups, current, bases = [], [], 0
    for seqid, length in contigs:
        estimate = plan['overhead_bytes'] + length * plan['bytes_per_base']
        if plan['budget_bytes'] is not None and estimate > plan['budget_bytes']:
            raise ValueError(f'Contig {seqid} ({length} bases) exceeds the estimated index budget; '
                             'use a larger allocation. Contigs are never cut across gene models.')
        if current and bases + length > plan['target_bases']:
            groups.append(current)
            current, bases = [], 0
        current.append([seqid, length])
        bases += length
    if current:
        groups.append(current)
    return groups
