"""Predictor inputs remain reproducible without repeated permanent FASTAs."""
import gzip
from pathlib import Path

import pytest

from workflow.support import rescue_search_inputs as inputs


class Genome:
    def __init__(self, sequence="aaATGAAATAAcc"):
        self.sequence = sequence

    def fetch(self, seqid, start, end):
        assert seqid == "chr1"
        return self.sequence[start:end]

    def get_reference_length(self, seqid):
        assert seqid == "chr1"
        return len(self.sequence)


def fixture(directory):
    proteins = {"Self": {"g": "MK"}, "Relative": {"g": "ML"}}
    regions = [{"id": "q1", "donor": "Self", "query": "g"},
               {"id": "q2", "donor": "Relative", "query": "g"}]
    with gzip.open(directory / inputs.LOCAL_RECORDS, "wt") as out:
        inputs.write_record(out, inputs.local_record(1, ("chr1", 2, 11), regions, proteins, "ATGAAATAA"))
    inputs.write_genome_records(directory, regions, {"q1": "q1", "q2": "q2"}, proteins)
    (directory / "genome_prediction_query_mapping.tsv").write_text(
        "candidate\trepresentative\nq1\tq1\nq2\tq2\n")
    return proteins, regions


def test_exact_input_regeneration_and_optional_combined_diagnostics(tmp_path):
    proteins, _ = fixture(tmp_path)
    result = inputs.export_inputs(tmp_path, tmp_path / "export", Genome(), proteins, combined=True)
    assert result == {"local_windows": 1, "combined_diagnostics": True}
    export = tmp_path / "export"
    assert (export / "intervals/1/region.fa").read_bytes() == b">interval\nATGAAATAA\n"
    assert (export / "intervals/1/queries.fa").read_bytes() == b">q1\nMK\n>q2\nML\n"
    assert (export / "queries.fa").read_bytes() == b">q1\nMK\n>q2\nML\n"
    assert (export / "regions.fa").read_bytes() == b">q1\nATGAAATAA\n>q2\nATGAAATAA\n"
    assert (export / "unresolved.unique.fa").read_bytes() == b">q1\nMK\n>q2\nML\n"
    assert (export / "unresolved.fa").read_bytes() == b">q1\nMK\n>q2\nML\n"
    inputs.export_inputs(tmp_path, tmp_path / "minimal", Genome(), proteins)
    assert not (tmp_path / "minimal/regions.fa").exists()
    assert not (tmp_path / "minimal/queries.fa").exists()
    with pytest.raises(FileExistsError):
        inputs.export_inputs(tmp_path, export, Genome(), proteins)


@pytest.mark.parametrize("change", ["genome", "protein", "coordinate", "index", "schema"])
def test_changed_or_invalid_inputs_cannot_be_regenerated(tmp_path, change):
    proteins, _ = fixture(tmp_path)
    genome = Genome()
    if change == "genome":
        genome.sequence = "aaATGACATAAcc"
    elif change == "protein":
        proteins["Self"]["g"] = "MM"
    else:
        row = next(inputs.records(tmp_path / inputs.LOCAL_RECORDS))
        if change == "coordinate":
            row["end"] = 100
        elif change == "index":
            row["index"] = "../../other"
        else:
            row["schema"] = 999
        with gzip.open(tmp_path / inputs.LOCAL_RECORDS, "wt") as out:
            inputs.write_record(out, row)
    with pytest.raises(ValueError):
        inputs.export_inputs(tmp_path, tmp_path / "out", genome, proteins)


def test_genome_aliases_are_not_repeated_as_query_bodies(tmp_path):
    proteins = {"Self": {"g": "MK"}, "Relative": {"g": "MK"}}
    regions = [{"id": "q1", "donor": "Self", "query": "g"},
               {"id": "q2", "donor": "Relative", "query": "g"}]
    inputs.write_genome_records(tmp_path, regions, {"q1": "q1", "q2": "q1"}, proteins)
    rows = list(inputs.records(tmp_path / inputs.GENOME_RECORDS))
    assert len(rows) == 1 and rows[0]["candidate"] == "q1"
    (tmp_path / "genome_prediction_query_mapping.tsv").write_text(
        "candidate\trepresentative\nq1\tq1\nq2\tq1\n")
    inputs.export_inputs(tmp_path, tmp_path / "export", Genome(), proteins, combined=True)
    assert (tmp_path / "export/unresolved.fa").read_text() == ">q1\nMK\n>q2\nMK\n"


@pytest.mark.parametrize("rows", ["q1\tq2\n", "q1\tq1\nq1\tq1\n", "q1\tq1\textra\n"])
def test_invalid_aliases_cannot_recreate_diagnostic_queries(tmp_path, rows):
    proteins, _ = fixture(tmp_path)
    (tmp_path / "genome_prediction_query_mapping.tsv").write_text("candidate\trepresentative\n" + rows)
    with pytest.raises(ValueError, match="alias"):
        inputs.export_inputs(tmp_path, tmp_path / "export", Genome(), proteins, combined=True)


@pytest.mark.parametrize("retain", [False, True])
def test_successful_interval_inputs_retained_only_on_request(tmp_path, monkeypatch, retain):
    from workflow.support import rescue_gene_models as rescue
    proteins = {"Self": {"g": "MK"}}
    regions = [{"id": "q1", "donor": "Self", "query": "g"}]
    windows = {("chr1", 2, 11): regions}
    def predict(command, directory, label, stdout):
        assert (directory / "region.fa").read_text() == ">interval\nATGAAATAA\n"
        assert (directory / "queries.fa").read_text() == ">q1\nMK\n"
        Path(stdout).write_text("##PAF\tq1\t2\t0\t2\t+\tinterval\t9\t0\t6\t2\t2\t60\tcg:Z:2M\n"
                               "interval\tminiprot\tmRNA\t1\t6\t.\t+\t.\tID=MP000001;Target=q1 1 2;Identity=1\n"
                               "interval\tminiprot\tCDS\t1\t6\t.\t+\t0\tParent=MP000001\n")
    monkeypatch.setattr(rescue, "run", predict)
    result = rescue.search_intervals(tmp_path, windows, proteins, Genome(), 1, 20000, 2,
                                     retain_inputs=retain)
    assert result[0]["cds"] == [[0, 6, 0]]
    assert (tmp_path / "intervals/1/region.fa").exists() is retain
    assert (tmp_path / "intervals/1/queries.fa").exists() is retain
    assert (tmp_path / "intervals/1/models.gff").exists()
    row = next(inputs.records(tmp_path / inputs.LOCAL_RECORDS))
    assert row["region_sha256"] == inputs.sequence_sha("ATGAAATAA")


def test_failed_interval_retains_actual_inputs_and_reconstruction_record(tmp_path, monkeypatch):
    from workflow.support import rescue_gene_models as rescue
    proteins = {"Self": {"g": "MK"}}
    regions = [{"id": "q1", "donor": "Self", "query": "g"}]
    def fail(*args):
        raise RuntimeError("predictor failed")
    monkeypatch.setattr(rescue, "run", fail)
    with pytest.raises(RuntimeError, match="predictor failed"):
        rescue.search_intervals(tmp_path, {("chr1", 2, 11): regions}, proteins, Genome(), 1, 20000, 1)
    assert (tmp_path / "intervals/1/region.fa").read_text() == ">interval\nATGAAATAA\n"
    assert len(list(inputs.records(tmp_path / inputs.LOCAL_RECORDS))) == 1


def test_streamed_alias_expansion_preserves_full_legacy_gff_and_model_order(tmp_path):
    from workflow.support import rescue_gene_models as rescue
    source = tmp_path / "unique.gff"
    source.write_text("##gff-version 3\n##PAF\tq1\t2\t0\t2\t+\tchr1\t9\t0\t6\t2\t2\t60\tcg:Z:2M\n"
                      "chr1\tminiprot\tmRNA\t1\t6\t.\t+\t.\tID=MP000001;Target=q1 1 2;Identity=1\n"
                      "chr1\tminiprot\tCDS\t1\t6\t.\t+\t0\tParent=MP000001\n")
    aliases = {"self_alias": "q1", "relative_alias": "q1"}
    destination = tmp_path / "legacy.gff"
    rescue.expand_miniprot_queries(source, destination, aliases)
    assert "".join(rescue.expand_miniprot_lines(source, aliases)).encode() == destination.read_bytes()
    expected = rescue.read_miniprot(destination)
    actual = list(rescue.iter_miniprot(rescue.expand_miniprot_lines(source, aliases)))
    assert actual == expected
    assert [m["query"] for m in actual] == list(aliases)
    assert [m["id"] for m in actual] == ["MP000001", "MP000002"]


@pytest.mark.parametrize("value", [None, [], {}, {"format": "future", "retain_search_inputs": False},
                                  {"format": [], "retain_search_inputs": False},
                                  {"format": "compact", "retain_search_inputs": 0},
                                  {"format": "compact", "retain_search_inputs": False, "unknown": True}])
def test_unknown_storage_contract_is_rejected_before_any_preparation(value, tmp_path, monkeypatch):
    from workflow.support import rescue_gene_models as rescue
    def unexpected(*args):
        raise AssertionError("Malformed storage must not prepare or publish any model")
    monkeypatch.setattr(rescue, "prepared", unexpected)
    plan = {"species": ["Self"], "request": {"output_storage": value}}
    with pytest.raises(ValueError, match="output-storage"):
        rescue.rescue(tmp_path, plan, "Self", 1)
    assert not list(tmp_path.iterdir())


def test_legacy_plans_keep_their_original_retention_contract():
    from workflow.support import rescue_gene_models as rescue
    assert rescue.output_storage({"request": {}}) == {"format": "legacy", "retain_search_inputs": True}


@pytest.mark.parametrize("command", ["export-search-inputs", "export-models"])
@pytest.mark.parametrize("location", ["worker", "other_worker", "prepared", "future_stage", "symlink"])
def test_review_export_cannot_add_unreceipted_files_to_frozen_rescue_output(tmp_path, monkeypatch, command, location):
    import sys

    from workflow.support import rescue_gene_models as rescue
    root = tmp_path / "producer"
    worker = root / "rescued/Self"
    worker.mkdir(parents=True)
    (worker / "receipt.json").write_text('{"key":"frozen","files":{"evidence.txt":"untouched"}}')
    (worker / "evidence.txt").write_text("original genomic evidence")
    locations = {"worker": worker / "export", "other_worker": root / "rescued/Other/export",
                 "prepared": root / "prepared/Self/export", "future_stage": root / "augmented/review"}
    if location == "symlink":
        alias = tmp_path / "producer_alias"
        alias.symlink_to(root, target_is_directory=True)
        destination = alias / "rescued/Self/export"
    else:
        destination = locations[location]
    before = {str(p.relative_to(root)): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    before_paths = {str(p.relative_to(root)) for p in root.rglob("*")}
    monkeypatch.setattr(rescue, "load", lambda *args, **kwargs: {"species": ["Self"]})
    def unexpected(*args, **kwargs):
        raise AssertionError("Unsafe export must stop before hashing or writing producer members")
    monkeypatch.setattr(rescue, "verified", unexpected)
    monkeypatch.setattr(rescue, "rescue_key", unexpected)
    monkeypatch.setattr(sys, "argv", ["rescue_gene_models.py", command, "--output", str(root),
                                    "--task-index", "1", "--destination", str(destination)])
    with pytest.raises(ValueError, match="outside the frozen rescue output"):
        rescue.main()
    assert not destination.exists()
    assert {str(p.relative_to(root)) for p in root.rglob("*")} == before_paths
    assert {str(p.relative_to(root)): p.read_bytes() for p in root.rglob("*") if p.is_file()} == before
