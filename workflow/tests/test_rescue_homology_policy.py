"""Coverage and identity trade-offs preserve explicit hard floors and old plans."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from rescue_homology_policy import homology_problems, homology_supported

POLICY = {'minimum_coverage': .70, 'minimum_identity': .50,
          'minimum_broad_coverage': .90, 'minimum_lower_coverage_identity': .70}


@pytest.mark.parametrize('coverage,identity,accepted', [
    (.90, .50, True), (.899, .699, False), (.899, .70, True),
    (.70, .70, True), (.699, 1., False), (.95, .499, False), (.808, .963, True),
])
def test_tiered_support(coverage, identity, accepted):
    assert homology_supported(coverage, identity, POLICY) is accepted


def test_explicit_stricter_floors_are_not_overridden_by_strong_identity():
    assert homology_problems(.90, 1., {**POLICY, 'minimum_coverage': .95}) == ['low_coverage']
    assert homology_problems(.80, .75, {**POLICY, 'minimum_identity': .80}) == ['low_identity']


def test_legacy_request_and_explicit_fixed_policy():
    legacy = {'minimum_coverage': .80, 'minimum_identity': .50}
    assert homology_supported(.80, .50, legacy)
    assert homology_supported(.70, .50, {**POLICY, 'minimum_broad_coverage': 0})


def test_alignment_numeric_scalars_are_supported():
    import numpy as np
    assert homology_supported(np.float64(.9), np.float64(.5), POLICY)


@pytest.mark.parametrize('coverage,identity', [(True, .8), (.8, float('nan')), (1.1, .8), (.8, -1)])
def test_invalid_alignment_evidence_cannot_support_presence(coverage, identity):
    with pytest.raises(ValueError, match='Invalid homology'):
        homology_supported(coverage, identity, POLICY)


@pytest.mark.parametrize('problems,boundary,supported', [
    ([], {}, True), (['low_coverage'], {}, False), (['frameshift'], {}, False),
    (['missing_stop'], {'supports_intact_homolog': True}, True),
])
def test_presence_evidence_survives_placement_but_not_failed_cds(problems, boundary, supported):
    from rescue_gene_models import check_interval
    m = {'sequence': 'ATGAAATAA', 'problems': list(problems), 'coverage': .8, 'identity': .8,
         'boundary_revalidation': boundary,
         'evidence': {'donor': 'External', 'target': 'Target', 'genome_only': True}}
    check_interval(m)
    assert m['presence_evidence']['supports_intact_homolog'] is supported
    assert 'unanchored_genome_search' in m['problems']
    assert m['presence_evidence']['orthology'] == m['presence_evidence']['expected_copy'] == 'unassigned'
