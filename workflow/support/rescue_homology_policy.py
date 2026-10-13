"""Coverage/identity evidence, independent of genomic structure or copy assignment."""
import math
from numbers import Real

DEFAULT_MINIMUM_COVERAGE = .70
DEFAULT_BROAD_COVERAGE = .90
DEFAULT_LOWER_COVERAGE_IDENTITY = .70


def homology_problems(coverage, identity, params):
    """Legacy requests without the broad-coverage field retain fixed cutoffs."""
    if not all(isinstance(v, Real) and not isinstance(v, bool)
               and math.isfinite(v) and 0 <= v <= 1 for v in (coverage, identity)):
        raise ValueError('Invalid homology coverage or identity')
    problems = []
    if coverage < params['minimum_coverage']:
        problems.append('low_coverage')
    required_identity = params['minimum_identity']
    if coverage < params.get('minimum_broad_coverage', 0):
        required_identity = max(required_identity, params.get('minimum_lower_coverage_identity', required_identity))
    if identity < required_identity:
        problems.append('low_identity')
    return problems


def homology_supported(coverage, identity, params):
    return not homology_problems(coverage, identity, params)
