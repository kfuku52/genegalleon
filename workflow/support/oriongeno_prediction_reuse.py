"""Verify an existing native prediction before CPU-only formatter resumption."""
import argparse
import json
from pathlib import Path

from artifact_provenance import audit_manifest, load_manifest, resolve_reference, runtime_support_root
from gene_family_output_store import GeneFamilyOutputStore


def verify(workspace, species, genome, gff, receipt, support_root):
    workspace = Path(workspace).resolve(strict=True)
    prediction = workspace / 'output/genome_prediction/oriongeno' / species
    manifest = prediction / 'prediction.provenance.json'
    payload = load_manifest(manifest)
    if (payload.get('family_id') != species
            or payload.get('step') != 'genome_annotation_oriongeno_prediction'):
        raise ValueError('Cached prediction belongs to another species or stage')
    logical = workspace / 'output/.gg_global_artifacts'
    for collection, expected in [('inputs', {'genome': genome}),
                                 ('outputs', {'gff': gff, 'native_receipt': receipt})]:
        entries = {item['label']: item for item in payload[collection]}
        for label, path in expected.items():
            if label not in entries or resolve_reference(entries[label], logical, workspace).resolve() != Path(path).resolve():
                raise ValueError('Cached prediction path differs: ' + label)
    with runtime_support_root(Path(support_root)):
        status, reason = audit_manifest(payload, GeneFamilyOutputStore(logical), logical, workspace)
    if status != 'current':
        raise ValueError('Cached prediction is not current: ' + reason)
    native = json.loads(Path(receipt).read_text())
    if native.get('native_returncode') != 0 or native.get('kind') != 'genegalleon-oriongeno-native-prediction':
        raise ValueError('Cached native prediction did not succeed')
    return dict(verified=True, species=species, prediction_reused=True,
                original_threads=payload['parameters'].get('threads'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('workspace', 'species', 'genome', 'gff', 'receipt', 'support-root'):
        parser.add_argument('--' + option, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.workspace, args.species, args.genome, args.gff,
                            args.receipt, args.support_root), sort_keys=True))
