"""Independent adversarial regressions for acquisition and publication."""

import json
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlparse

import pandas as pd
import pytest
from test_generate_species_trait import load_script_module, run_script

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from gift_retrieval import GiftRetrieval, load_reviewed_mappings  # noqa: E402

pytestmark = [pytest.mark.runtime, pytest.mark.integration]
URL = "https://example.org/index3.2.php?query=traits&traitid=1.1.1&startat=0&limit=2"


@pytest.fixture
def g():
    return load_script_module()


def test_refresh_resume_never_combines_new_and_old_pages(tmp_path):
    tail = URL.replace("startat=0", "startat=2")

    def old(url, timeout):
        return [{"work_ID": "1", "trait_value": "old"}]

    with GiftRetrieval(tmp_path, old) as client:
        client.fetch(URL, 1)
        client.fetch(tail, 1)

    def interrupted(url, timeout):
        if url == tail:
            raise URLError("refresh interrupted")
        return [{"work_ID": "1", "trait_value": "new"}]

    with pytest.raises(URLError):
        with GiftRetrieval(tmp_path, interrupted, mode="refresh", retries=0) as client:
            client.fetch(URL, 1)
            client.fetch(tail, 1)
    calls = []

    def finish(url, timeout):
        calls.append(url)
        return [{"work_ID": "2", "trait_value": "new tail"}]

    with GiftRetrieval(tmp_path, finish) as client:
        assert client.fetch(URL, 1)[0]["trait_value"] == "new"
        assert client.fetch(tail, 1)[0]["trait_value"] == "new tail"
    assert calls == [tail]


def test_existing_reader_pins_its_cache_generation_during_refresh(tmp_path):
    tail = URL.replace("startat=0", "startat=2")

    def old(**kwargs):
        return [{"work_ID": "1", "trait_value": "old"}]

    with GiftRetrieval(tmp_path, old) as writer:
        writer.fetch(URL, 1)
        writer.fetch(tail, 1)
    with GiftRetrieval(tmp_path, old) as reader:
        reader.fetch(URL, 1)
        with GiftRetrieval(
            tmp_path, lambda **kwargs: [{"work_ID": "1", "trait_value": "new"}], mode="refresh"
        ) as writer:
            writer.fetch(URL, 1)
            writer.fetch(tail, 1)
        assert reader.fetch(tail, 1)[0]["trait_value"] == "old"


@pytest.mark.parametrize("bad", [[], None, "not a record", {"schema": "bad"}])
def test_bad_cache_container_is_a_recoverable_cache_miss(tmp_path, bad):
    def fetch(**kwargs):
        return [{"work_ID": "1", "trait_value": "good"}]

    with GiftRetrieval(tmp_path, fetch) as client:
        client.fetch(URL, 1)
    page = next((tmp_path / "pages").rglob("*.json"))
    page.write_text(json.dumps(bad))
    with GiftRetrieval(tmp_path, fetch) as client:
        assert client.fetch(URL, 1)[0]["trait_value"] == "good"


@pytest.mark.parametrize(
    "row",
    [
        {"work_ID": "1", "error": "failed"},
        {"work_ID": "wrong", "trait_value": "x"},
        {"work_ID": "1", "trait_value": "x", "agreement": "broken"},
        {"work_ID": "1", "trait_value": "x", "agreement": "2"},
    ],
)
def test_invalid_api_rows_never_become_cached_data(tmp_path, row):
    with pytest.raises(ValueError):
        with GiftRetrieval(tmp_path, lambda **kwargs: [row]) as client:
            client.fetch(URL, 1)
    assert not list((tmp_path / "pages").rglob("*.json"))


def test_http_date_retry_after_is_respected(tmp_path, monkeypatch):
    monkeypatch.setattr("gift_retrieval.time.sleep", lambda value: None)
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        if len(calls) == 1:
            raise HTTPError(url, 429, "busy", {"Retry-After": "Wed, 01 Jan 2098 00:00:00 GMT"}, None)
        return []

    with pytest.raises(HTTPError):
        with GiftRetrieval(tmp_path, fetch) as client:
            client.fetch(URL, 1)
    assert len(calls) == 1


def test_ambiguous_trait_group_cannot_choose_the_most_common_trait(g):
    meta = [
        {"Lvl3": "1.6.1", "Trait1": "Plant_height", "Trait2": "Plant_height_min", "count": 2},
        {"Lvl3": "1.6.2", "Trait1": "Plant_height", "Trait2": "Plant_height_max", "count": 999},
    ]
    with pytest.raises(ValueError, match="ambiguous"):
        g.resolve_gift_trait_token_map(
            "https://example.org/index3.2.php", ["Plant_height"], 1, fetcher=lambda **kwargs: meta
        )


def test_multiple_tokens_for_one_trait_are_not_double_counted(g, tmp_path, monkeypatch):
    def fetch(url, timeout):
        query = parse_qs(urlparse(url).query)["query"][0]
        if query == "traits_meta":
            return [{"Lvl3": "1.1.1", "Trait1": "Woodiness", "Trait2": "Woodiness_1", "count": 1}]
        if query == "names_matched_unique":
            return [{"work_ID": "1", "work_species": "Arabidopsis thaliana"}]
        return [{"work_ID": "1", "trait_value": "woody"}]

    monkeypatch.setattr(g, "fetch_json_payload", fetch)
    plans = [
        g.TraitPlanRow("gift", "trait_value", "trait", "binary", "sum", {"woody"}, token, "trait_ID")
        for token in ["1.1.1", "Woodiness_1"]
    ]
    table = g.fetch_gift_api_table(
        "gift",
        {"uri": "https://example.org/", "gift_version": "3.2"},
        plans,
        ["Arabidopsis_thaliana"],
        1,
        False,
        tmp_path,
    )
    assert len(table) == 1
    assert table.attrs["gift_trait_token_map"] == {"1.1.1": "1.1.1", "Woodiness_1": "1.1.1"}


def test_explicitly_unresolved_name_cannot_become_an_exact_match(g):
    report = []
    table = g.resolve_gift_species_map(
        "https://example.org/index3.2.php",
        ["Arabidopsis_thaliana"],
        1,
        fetcher=lambda **kwargs: [{"work_ID": "1", "work_species": "Arabidopsis thaliana", "resolved": 0}],
        report=report,
    )
    assert table.empty


@pytest.mark.parametrize(
    "row",
    [
        "gift\tvalue\theight\tnumric\tmean",
        "gift\tvalue\theight\tnumeric\tmen",
        "gift\tvalue\tspecies\tnumeric\tmean",
        "gift\tvalue\theight\tnumeric\tmean\textra",
    ],
)
def test_malformed_plans_fail_before_retrieval(g, tmp_path, row):
    p = tmp_path / "plan.tsv"
    p.write_text("database\tsource_column\toutput_trait\tvalue_type\taggregation\n" + row + "\n")
    with pytest.raises(ValueError):
        g.read_trait_plan(p)


@pytest.mark.parametrize("value", ["inf", "-inf", "broken"])
def test_invalid_numeric_values_cannot_silently_disappear(g, value):
    frame = pd.DataFrame({"__species_norm": ["a"], "value": [value]})
    plan = g.TraitPlanRow("gift", "value", "trait", "numeric", "median", set(), "", "")
    with pytest.raises(ValueError):
        g.aggregate_trait_column(frame, plan)


def test_text_traits_keep_all_distinct_values(g):
    frame = pd.DataFrame({"__species_norm": ["a"] * 4, "value": ["a | b", "c", "c", None]})
    plan = g.TraitPlanRow("gift", "value", "trait", "text", "unique", set(), "", "")
    result = g.aggregate_trait_column(frame, plan)
    assert json.loads(result["a"]) == ["a | b", "c"]


def test_numeric_output_preserves_double_precision(g):
    value = 1.2345678901234567
    assert float(g.format_output_value(value)) == value


def test_output_and_stats_cannot_overwrite_an_input_or_each_other(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("species_key\nArabidopsis_thaliana\n")
    raw = tmp_path / "raw.tsv"
    raw.write_text("species\tvalue\nArabidopsis_thaliana\t1\n")
    plan = tmp_path / "plan.tsv"
    plan.write_text("database\tsource_column\toutput_trait\naustraits\tvalue\theight\n")
    sources = tmp_path / "sources.tsv"
    sources.write_text(f"database\turi\tspecies_column\naustraits\t{raw}\tspecies\n")
    output = tmp_path / "result.tsv"
    output.write_text("previous\n")
    args = [
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(plan),
        "--database-sources",
        str(sources),
        "--downloads-dir",
        str(tmp_path / "downloads"),
        "--output",
        str(output),
    ]
    for destination in [output, manifest, raw]:
        original = destination.read_bytes()
        result = run_script(*args, "--stats-output", str(destination))
        assert result.returncode != 0
        assert destination.read_bytes() == original
        assert output.read_text() == "previous\n"


def test_mapping_duplicate_headers_are_rejected(tmp_path):
    p = tmp_path / "mapping.tsv"
    p.write_text(
        "species\tspecies\tgift_version\twork_ID\twork_species\tdecision\tevidence_url\n"
        "Ignored_label\tA_b\t3.2\t1\tC d\tmap\thttps://example.org/evidence\n"
    )
    with pytest.raises(ValueError):
        load_reviewed_mappings("3.2", p, include_bundled=False)


def test_publication_failure_restores_both_previous_outputs(g, tmp_path, monkeypatch):
    first, second = tmp_path / "first.tsv", tmp_path / "stats.json"
    first.write_text("previous table")
    second.write_text("previous stats")
    replace = g.os.replace

    def fail_second(source, target):
        if target == second:
            raise OSError("injected second installation failure")
        return replace(source, target)

    monkeypatch.setattr(g.os, "replace", fail_second)
    with pytest.raises(OSError, match="injected"):
        g.publish_trait_outputs({first: b"new table", second: b"new stats"})
    assert first.read_text() == "previous table"
    assert second.read_text() == "previous stats"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["first.tsv", "stats.json"]


def test_qualified_taxa_are_not_silently_broadened(g):
    def no_request(**kwargs):
        raise AssertionError("No binomial lookup should run for a qualified target")

    decisions = []
    result = g.resolve_gift_species_map(
        "https://example.org/index3.2.php", ["A_b_subsp._c"], 1, fetcher=no_request, report=decisions
    )
    assert result.empty and decisions[0]["status"] == "qualified_taxon_requires_review"


def test_numeric_aggregation_overflow_is_rejected(g):
    frame = pd.DataFrame({"__species_norm": ["a", "a"], "value": [1e308, 1e308]})
    plan = g.TraitPlanRow("gift", "value", "trait", "numeric", "mean", set(), "", "")
    with pytest.raises(ValueError, match="overflowed"):
        g.aggregate_trait_column(frame, plan)


def test_cli_named_and_id_aliases_share_observations(g, tmp_path, monkeypatch):
    def fetch(url, timeout):
        if "traits_meta" in url:
            return [{"Lvl3": "1.1.1", "Trait1": "Woodiness", "Trait2": "Woodiness_1"}]
        if "names_matched" in url:
            return [{"work_ID": "1", "work_species": "Arabidopsis thaliana"}]
        return [{"work_ID": "1", "trait_value": "woody"}]

    monkeypatch.setattr(g, "fetch_json_payload", fetch)
    manifest, plan, sources, output = [
        tmp_path / name for name in ["manifest.tsv", "plan.tsv", "sources.tsv", "out.tsv"]
    ]
    manifest.write_text("species_key\nArabidopsis_thaliana\n")
    plan.write_text(
        "database\tsource_column\toutput_trait\tvalue_type\taggregation\tpositive_values\ttrait_key\ttrait_key_column\n"
        "gift\ttrait_value\tby_id\tbinary\tsum\twoody\t1.1.1\ttrait_ID\n"
        "gift\ttrait_value\tby_name\tbinary\tsum\twoody\tWoodiness_1\ttrait_ID\n"
    )
    sources.write_text("database\turi\tgift_version\ngift\thttps://example.org/\t3.2\n")
    args = [
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(plan),
        "--database-sources",
        str(sources),
        "--downloads-dir",
        str(tmp_path / "downloads"),
        "--output",
        str(output),
        "--strict",
    ]
    assert g.main(args) == 0
    table = pd.read_csv(output, sep="\t")
    assert table.by_id.tolist() == table.by_name.tolist() == [1]
    original = output.read_bytes()
    assert g.main(args + ["--dry-run"]) == 0
    assert output.read_bytes() == original


def test_truncated_http_body_is_retried(tmp_path, monkeypatch):
    from http.client import IncompleteRead

    monkeypatch.setattr("gift_retrieval.time.sleep", lambda seconds: None)
    calls = []

    def fetch(**kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise IncompleteRead(b"[", 10)
        return []

    with GiftRetrieval(tmp_path, fetch) as client:
        assert client.fetch(URL, 1) == []
    assert len(calls) == 2


def test_text_csv_reader_preserves_literal_missing_words(g, tmp_path):
    path = tmp_path / "text.tsv"
    path.write_text("species\tvalue\nA_b\tNA\nA_b\tnone\nA_b\t\n")
    frame = g.read_table(path, "tsv")
    frame["__species_norm"] = frame.species
    plan = g.TraitPlanRow("austraits", "value", "trait", "text", "unique", set(), "", "")
    assert json.loads(g.aggregate_trait_column(frame, plan)["A_b"]) == ["NA", "none"]
    frame2 = g.read_table_from_text(path.read_text(), "tsv")
    pd.testing.assert_frame_equal(frame.drop(columns="__species_norm"), frame2)
