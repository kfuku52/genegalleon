"""Publish split candidates from an immutable rescue without replaying Miniprot."""
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

try:
    import annevo_gene_models as native
    import refinement_native_splits as splits
    import rescue_gene_models as rescue
    from gene_model_catalog import build_catalog, indexed_genome
    from gene_model_species_profiles import parameters_for
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest
except ImportError:
    from . import annevo_gene_models as native
    from . import refinement_native_splits as splits
    from . import rescue_gene_models as rescue
    from .gene_model_catalog import build_catalog, indexed_genome
    from .gene_model_species_profiles import parameters_for
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest

DIRECTORY = 'native_split_candidates'
CANDIDATES = 'candidates.json'
POLICY_FIELDS = ('minimum_coverage', 'minimum_identity', 'minimum_broad_coverage',
                 'minimum_lower_coverage_identity', 'max_interval', 'max_intron')


def policy(params):
    return {field: params[field] for field in POLICY_FIELDS if field in params}


def eligible_species(plan):
    predictions = plan['request'].get('oriongeno', {}).get('predictions', {})
    return {name for name in plan['species'] if name in predictions and len(plan['donors'][name]) >= 2}


def _read_receipt(directory, proof):
    path = str(directory / 'receipt.json')
    raw = Path(path).read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    proof.read([path], expected={path: sha})
    return json.loads(raw), sha


def _member(directory, receipt, filename, proof):
    path = str(directory / filename)
    sha = receipt['files'].get(filename)
    if not sha:
        raise ValueError('Incomplete rescue publication: ' + path)
    proof.read([path], expected={path: sha})
    return path, sha


def _generated_inputs(root, plan, name, plan_sha, proof):
    worker, worker_sha = _read_receipt(root / 'rescued' / name, proof)
    if worker['key'].get('plan') != plan_sha or worker['key'].get('species') != name:
        raise ValueError('Rescue worker belongs to a different plan: ' + name)
    effective = root / 'effective' / name
    if (root / 'augmented' / 'receipt.json').is_file():
        effective = root / 'augmented'
    publication, _ = _read_receipt(effective, proof)
    key = publication['key']
    if key.get('plan') != plan_sha or key.get('rescue_receipts', {}).get(name) != worker_sha:
        raise ValueError('Rescue annotation belongs to a different worker: ' + name)
    cds, cds_sha = _member(effective, publication, f'species_cds/{name}.rescue.cds.fa', proof)
    gff, gff_sha = _member(effective, publication, f'species_gff/{name}.rescue.gff3', proof)
    source = plan['request']['sources'][name]
    genome = source['genome']
    genome_sha = plan['request']['files'][genome]
    proof.read([genome], expected={genome: genome_sha})
    return {'cds': cds, 'gff': gff, 'genome': genome}, {'cds': cds_sha, 'gff': gff_sha, 'genome': genome_sha}, worker_sha


def _originals(catalog, targets):
    originals, duplicates = {}, set()
    for locus in catalog['loci']:
        candidates = [c for c in locus['candidates'] if c['origin'] == 'original']
        if locus.get('ambiguous_coordinates') or len(candidates) != 1:
            continue
        model, = candidates
        owner = model['source_gene_id']
        quality = model['quality']
        if owner not in targets or not quality['usable'] or any(quality.get(flag) for flag in (
                'annotated_exception', 'structure_problem', 'internal_stop', 'ambiguous', 'invalid_base',
                'translation_uncertain', 'phase_unresolved', 'sequence_mismatch')):
            continue
        if owner in originals:
            duplicates.add(owner)
        originals[owner] = model
    return {owner: model for owner, model in originals.items() if owner not in duplicates}


def publish(root, plan, names, cpus):
    """Read only the annotation, native paths and proteins used by the new step."""
    names = sorted(set(names) & eligible_species(plan))
    if not names:
        return
    diamond = shutil.which('diamond')
    if not diamond:
        raise ValueError('DIAMOND is required for rescue split candidates')
    proof = FreshDigestBatch()
    executable = {'path': diamond, 'sha256': proof.read([diamond])[diamond]}
    plan_sha = rescue.plan_digest(root, plan)
    plan_path = str(root / 'plan.json')
    proof.read([plan_path], expected={plan_path: plan_sha})
    cohorts = {name: plan['donors'][name] for name in names}
    proteins, protein_hashes = {}, {}
    for donor in sorted({donor for donors in cohorts.values() for donor in donors}):
        directory = root / 'prepared' / donor
        receipt, _ = _read_receipt(directory, proof)
        if receipt['key'] != {'plan': plan_sha, 'species': donor}:
            raise ValueError('Rescue donor preparation belongs to a different plan: ' + donor)
        proteins[donor], protein_hashes[donor] = _member(directory, receipt, 'genes.pep', proof)
    oriongeno = plan['request']['oriongeno']
    panel = {'donors': cohorts, 'proteins': proteins, 'executable': executable}
    implementation = {Path(file).name: digest(file) for file in (__file__, splits.__file__)}
    with tempfile.TemporaryDirectory(prefix='.split-replay-', dir=root) as scratch:
        search = splits.DonorSearch(panel, None, Path(scratch) / 'search', cpus)
        for name in names:
            paths, source_sha, worker_sha = _generated_inputs(root, plan, name, plan_sha, proof)
            native.verify_inputs(oriongeno, [name], verification=proof)
            prediction = oriongeno['predictions'][name]
            prediction_sha = oriongeno['files'][prediction['gff']]
            params = {**parameters_for(plan['request'], name), 'min_support': 2}
            code = plan['request']['sources'][name]['genetic_code']
            key = {'schema': 1, 'plan': plan_sha, 'species': name, 'worker': worker_sha, 'annotation': source_sha,
                   'oriongeno_gff': prediction_sha, 'policy': policy(params), 'implementation': implementation,
                   'donors': {d: protein_hashes[d] for d in cohorts[name]}, 'executable': executable}

            def build(tmp, paths=paths, name=name, prediction=prediction, code=code, params=params):
                with rescue.open_text(Path(paths['gff'])) as handle:
                    ownership = rescue.OwnershipIndex(rescue.original_ownership(handle, rescue.attributes))
                with indexed_genome(paths['genome']) as genome:
                    predicted = native.read_paths(prediction['gff'], genome, engine='oriongeno')
                    native_paths = [path for path in predicted if len(ownership.overlapping(path)) >= 2]
                    targets = {row['gene_id'] for path in native_paths for row in ownership.overlapping(path)}
                    candidates = []
                    if targets:
                        catalog = build_catalog(name, paths['cds'], paths['gff'], paths['genome'], code)
                        nominations = splits.nominate(native_paths, _originals(catalog, targets), ownership, genome,
                                                      code, params, name, set())
                        for shape, support in search.support(name, nominations, params).items():
                            candidates.append({**nominations[shape]['review'], 'support': support})
                atomic_json(tmp / CANDIDATES, sorted(candidates, key=lambda row: row['candidate_id']))

            rescue.stage(root, Path(DIRECTORY) / name, key, build, guard=proof.check, publication_guard=proof.check)
    proof.check()


def read_publication(root, name, plan_sha, worker_sha, source, prediction_sha, params, proof):
    """Bind reusable candidates to the exact annotation and native prediction used by refine."""
    directory = root / DIRECTORY / name
    if not (directory / 'receipt.json').is_file():
        return None
    receipt, _ = _read_receipt(directory, proof)
    key = receipt['key']
    if key.get('schema') != 1 or key.get('plan') != plan_sha or key.get('species') != name or key.get('worker') != worker_sha:
        raise ValueError('Rescue split publication belongs to a different plan or worker: ' + name)
    if key.get('oriongeno_gff') != prediction_sha or key.get('policy') != policy(params):
        return None
    source_paths = {role: str(source['fasta' if role == 'cds' else role]) for role in ('cds', 'gff', 'genome')}
    observed = proof.read(source_paths.values())
    if {role: observed[path] for role, path in source_paths.items()} != key['annotation']:
        raise ValueError('Rescue split publication was built from a different annotation: ' + name)
    path, _ = _member(directory, receipt, CANDIDATES, proof)
    rows = json.loads(Path(path).read_text())
    proof.check()
    return rows, [str(directory / 'receipt.json'), path]
