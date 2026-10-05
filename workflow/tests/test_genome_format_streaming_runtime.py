"""Verify streamed genome output against the real seqkit compressor."""

import pytest

from workflow.tests.test_format_species_inputs import (
    test_streamed_genome_matches_record_writer_across_transport_and_edge_cases as compare_genome,
)
from workflow.tests.test_format_species_inputs import (
    test_streamed_genome_read_failure_preserves_published_output as check_atomic_failure,
)


@pytest.mark.parametrize("transport", ["plain", "gz", "bz2", "tar", "tar-single"])
def test_real_seqkit_streamed_genome_equivalence(tmp_path, monkeypatch, transport):
    compare_genome(tmp_path, monkeypatch, transport, "seqkit")


def test_real_seqkit_streamed_genome_failure_is_atomic(tmp_path, monkeypatch):
    check_atomic_failure(tmp_path, monkeypatch, "seqkit")
