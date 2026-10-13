"""Read frozen species plans and completion records without certifying content hashes."""
import csv
import hashlib
import io
from pathlib import Path

from workflow_observation import read_regular_bytes, read_regular_json, strict_json_loads

MAX_BYTES = 8 * 1024 * 1024


def inspect(request):
    if (not isinstance(request, dict) or set(request) != {'cohort', 'sources', 'lineage'}
            or not isinstance(request['cohort'], list) or not isinstance(request['sources'], list)
            or not isinstance(request['lineage'], str) or not request['lineage']
            or not 1 <= len(request['cohort']) <= 32 or not 1 <= len(request['sources']) <= 64):
        raise ValueError('Invalid input-generation inventory request')
    plans = {}
    for item in request['cohort'] + request['sources']:
        if (not isinstance(item, dict) or set(item) != {'path', 'sha256'}
                or not isinstance(item['path'], str) or not Path(item['path']).is_absolute()):
            raise ValueError('Invalid frozen plan binding')
        path = Path(item['path'])
        raw = read_regular_bytes(path, MAX_BYTES)
        if hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise ValueError('Frozen input-generation plan changed')
        plan = strict_json_loads(raw)
        tasks = plan['tasks']
        if not isinstance(tasks, list) or not 1 <= len(tasks) <= 75000 or plan['task_count'] != len(tasks):
            raise ValueError('Invalid input-generation task count')
        names = [task['species_prefix'] for task in tasks]
        if not all(isinstance(n, str) and n for n in names) or len(set(names)) != len(names):
            raise ValueError('Duplicate or invalid species identity')
        if str(path) in plans and plans[str(path)] != (item['sha256'], names):
            raise ValueError('Conflicting frozen plan binding')
        plans[str(path)] = (item['sha256'], names)
    cohort = set().union(*(set(plans[item['path']][1]) for item in request['cohort']))
    if len(cohort) > 75000:
        raise ValueError('Species cohort exceeds inventory limit')
    recorded, invalid, finalizers = set(), set(), []
    for item in request['sources']:
        path = Path(item['path'])
        sha, names = plans[str(path)]
        settings = read_regular_json(str(path) + '.settings.json', MAX_BYTES)
        if settings.get('busco_lineage') != request['lineage']:
            continue
        local = set()
        for index, name in enumerate(names, 1):
            if name not in cohort:
                continue
            receipt_path = Path(str(path) + '.completed') / (str(index) + '.json')
            if not receipt_path.exists():
                continue
            receipt = read_regular_json(receipt_path, MAX_BYTES)
            if (receipt.get('plan_sha256') != sha or receipt.get('task_index') != index
                    or receipt.get('species_prefix') != name
                    or not isinstance(receipt.get('files'), dict) or not receipt['files']):
                invalid.add(name)
                continue
            recorded.add(name)
            local.add(name)
        root = path.parent.parent
        summary = root / 'gg_input_generation_runs.tsv'
        if set(names) == cohort and local == cohort and summary.is_file():
            rows = list(csv.DictReader(io.StringIO(read_regular_bytes(summary, MAX_BYTES).decode()), delimiter='\t'))
            if any(row.get('input_generation_mode') == 'array_finalize'
                   and row.get('exit_code') == '0'
                   and row.get('stage_multispecies_summary_status') in {'ok', 'skipped'}
                   for row in rows):
                finalizers.append(str(path))
    # Publication records are useful progress evidence, but this fast read does
    # not hash raw genomes or outputs and must never become workflow verification.
    return {'schema': 'input-generation-receipt-inventory-v1', 'expected': len(cohort),
            'recorded': len(recorded), 'remaining': len(cohort - recorded),
            'invalid': len(invalid - recorded), 'finalize_recorded': bool(finalizers),
            'checksum_verified': False, 'workflow_verified': False,
            'species': [{'species': name, 'receipt_recorded': name in recorded}
                        for name in sorted(cohort)]}
