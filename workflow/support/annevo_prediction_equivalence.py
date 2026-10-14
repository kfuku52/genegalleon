"""Strict native-feature equivalence, retaining independent artifact provenance."""
import hashlib
from pathlib import Path

from annevo_runtime_download import file_sha

ARTIFACT_FIELDS = {'receipt_sha256', 'GFF_sha256', 'GFF_bytes'}


def feature_digest(path):
    path = Path(path)
    before = file_sha(path)
    with path.open() as stream:
        rows = []
        for line in stream:
            if line.strip() and not line.startswith('#'):
                row = line.rstrip('\n')
                if len(row.split('\t')) != 9:
                    raise ValueError('GFF comparison requires all nine annotation columns')
                rows.append(row)
    value = hashlib.sha256(('\n'.join(sorted(rows))+'\n').encode()).hexdigest()
    if file_sha(path) != before:
        raise ValueError('GFF changed during feature comparison')
    return value


def matching_features(left, right):
    first, second = feature_digest(left), feature_digest(right)
    if first != second:
        raise ValueError('Native GFF coordinates, phase, score, identifiers or attributes changed')
    return first


def decision_fields(value):
    """Call only after native GFF features and model/genome inputs agree.

    Raw GFF comments and independent run receipts can change their artifact
    hashes. These three provenance fields are not biological decision inputs.
    All other fields, including model/runtime hashes and alignments, remain.
    """
    if isinstance(value, dict):
        return {k: decision_fields(v) for k, v in value.items() if k not in ARTIFACT_FIELDS}
    if isinstance(value, list):
        return [decision_fields(v) for v in value]
    return value
