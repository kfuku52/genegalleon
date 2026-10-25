"""Independently rederive historical transcript ORFs without rewriting audits."""

import json
import re

from transcriptome_orfs import audit_path, contract, input_type, signature


def historical_contract(task, output):
    """Verify every source/tool field, retaining the recorded grouping policy.

    This proof is only used by independent longest-CDS rederivation. It never
    authorizes formatter reuse and is insufficient without exact sequence/ID
    comparisons against the rederived records.
    """
    if input_type(task) != 'transcriptome':
        raise ValueError('Historical ORF validation requires transcript input')
    audit = json.loads(audit_path(output).read_text())
    original, current = audit['contract'], contract(task)
    comparison = dict(current)
    if 'gene_grouping' not in original:
        comparison.pop('gene_grouping')
    old_implementation = original.get('implementation', '')
    if not isinstance(old_implementation, str) or not re.fullmatch('[0-9a-f]{64}', old_implementation):
        raise ValueError('Invalid historical ORF implementation identity')
    comparison['implementation'] = old_implementation
    if original != comparison or audit['output'] != signature(output):
        raise ValueError('Historical ORF source, tool, policy or output changed')
    return dict(original_implementation=old_implementation, current_implementation=current['implementation'],
                gene_grouping=original.get('gene_grouping', 'historical-header-groups-v1'),
                audit_rewritten=False, requires_independent_sequence_validation=True)
