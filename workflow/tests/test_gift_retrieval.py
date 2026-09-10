"""GIFT retrieval contracts: incomplete requests, cache identity and taxonomy."""

import importlib.util
import json
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlparse

import pandas as pd
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
from gift_retrieval import GiftRetrieval, load_reviewed_mappings  # noqa: E402

pytestmark = [pytest.mark.runtime, pytest.mark.integration]


@pytest.fixture
def g():
    spec = importlib.util.spec_from_file_location(
        "generate_species_trait_gift_test", SUPPORT / "generate_species_trait.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


URL = "https://example.org/index3.2.php?query=traits&traitid=1.1.1&startat=0&limit=2"


def test_cache_reuses_only_validated_exact_requests(tmp_path):
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        return [{"work_ID": "1", "trait_value": "woody"}]

    with GiftRetrieval(tmp_path, fetch) as client:
        original = client.fetch(URL, 10)
    with GiftRetrieval(tmp_path, fetch) as client:
        assert client.fetch(URL, 10) == original
        client.fetch(URL + "&biasref=0", 10)
    assert calls == [URL, URL + "&biasref=0"]
    assert client.report["requests"][0]["source"] == "cache"
    page = next(path for path in (tmp_path / "pages").rglob("*.json") if json.loads(path.read_text())["url"] == URL)
    record = json.loads(page.read_text())
    record["payload"][0]["trait_value"] = "tampered"
    page.write_text(json.dumps(record))
    with GiftRetrieval(tmp_path, fetch, mode="offline") as client:
        with pytest.raises(ValueError, match="Invalid offline"):
            client.fetch(URL, 10)
    with GiftRetrieval(tmp_path, fetch) as client:
        assert client.fetch(URL, 10) == original
    assert calls.count(URL) == 2


@pytest.mark.parametrize("payload", [{"error": "database unavailable"}, [{"error": "bad query"}], [None], ["error"]])
def test_error_responses_are_never_cached_as_empty(tmp_path, payload):
    with pytest.raises(ValueError):
        with GiftRetrieval(tmp_path, lambda **kwargs: payload) as client:
            client.fetch(URL, 10)
    assert not (tmp_path / "pages").exists()
    assert json.loads(client.report_path.read_text())["status"] == "failed"


def test_refresh_offline_and_version_discovery(tmp_path):
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        return [{"version": "3.2"}] if "versions" in url else [{"work_ID": "1", "trait_value": "x"}]

    with GiftRetrieval(tmp_path, fetch) as client:
        client.fetch(URL, 10)
    with GiftRetrieval(tmp_path, fetch, mode="refresh") as client:
        client.fetch(URL, 10)
        client.fetch(URL, 10)
    assert len(calls) == 2
    with GiftRetrieval(tmp_path, fetch, mode="offline") as client:
        client.fetch(URL, 10)
        with pytest.raises(FileNotFoundError):
            client.fetch(URL + "&biasderiv=0", 10)
        with pytest.raises(ValueError, match="explicit stable"):
            client.fetch("https://example.org/index.php?query=versions", 10)
    for _ in range(2):
        with GiftRetrieval(tmp_path, fetch) as client:
            client.fetch("https://example.org/index.php?query=versions", 10)
    assert len(calls) == 4


def test_only_transient_errors_are_retried(tmp_path, monkeypatch):
    monkeypatch.setattr("gift_retrieval.time.sleep", lambda delay: None)
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        if len(calls) < 3:
            raise URLError("temporarily unavailable")
        return []

    with GiftRetrieval(tmp_path, fetch) as client:
        assert client.fetch(URL, 10) == []
    assert len(calls) == 3
    calls.clear()

    def bad_request(url, timeout):
        calls.append(url)
        raise HTTPError(url, 400, "bad request", {}, None)

    with GiftRetrieval(tmp_path, bad_request, mode="refresh") as client:
        with pytest.raises(HTTPError):
            client.fetch(URL, 10)
    assert len(calls) == 1


def config_and_plan(g):
    return {"uri": "https://example.org/", "gift_version": "3.2", "gift_page_size": "2", "gift_retries": "0"}, [
        g.TraitPlanRow("gift", "trait_value", "woodiness", "categorical", "mode", set(), "1.1.1", "trait_ID")
    ]


def test_interrupted_pages_resume_and_match_a_clean_run(g, tmp_path, monkeypatch):
    config, plans = config_and_plan(g)
    requests, fail = [], [True]

    def fetch(url, timeout):
        q = parse_qs(urlparse(url).query)
        requests.append(q)
        if q["query"] == ["names_matched_unique"]:
            return [{"work_ID": "1", "work_species": "Arabidopsis thaliana"}]
        assert q["limit"] == ["2"]
        if q["startat"] == ["0"]:
            return [{"work_ID": "1", "trait_value": "woody"}, {"work_ID": "99", "trait_value": "non-woody"}]
        if fail[0]:
            raise URLError("interrupted page 2")
        return []

    monkeypatch.setattr(g, "fetch_json_payload", fetch)
    with pytest.raises(URLError):
        g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)
    fail[0] = False
    requests.clear()
    resumed = g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)
    assert len(requests) == 1 and requests[0]["startat"] == ["2"]
    config["gift_cache_mode"] = "refresh"
    fresh = g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)
    pd.testing.assert_frame_equal(resumed, fresh)
    config["gift_cache_mode"] = "offline"
    requests.clear()
    offline = g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)
    pd.testing.assert_frame_equal(resumed, offline)
    assert requests == []
    reports = [json.loads(p.read_text()) for p in (tmp_path / "gift/runs").glob("*.json")]
    assert sum(r["status"] == "failed" for r in reports) == 1


def test_page_cap_cannot_publish_an_incomplete_trait(g, tmp_path, monkeypatch):
    config, plans = config_and_plan(g)
    config["gift_max_pages_per_trait"] = "1"

    def fetch(url, timeout):
        if "names_matched" in url:
            return [{"work_ID": "1", "work_species": "Arabidopsis thaliana"}]
        return [{"work_ID": "1", "trait_value": "woody"}, {"work_ID": "2", "trait_value": "woody"}]

    monkeypatch.setattr(g, "fetch_json_payload", fetch)
    with pytest.raises(ValueError, match="acquisition incomplete"):
        g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)
    config.pop("gift_max_pages_per_trait")
    with pytest.raises(ValueError, match="repeated page"):
        g.fetch_gift_api_table("gift", config, plans, ["Arabidopsis_thaliana"], 10, False, tmp_path)


def test_reviewed_taxa_override_old_synonyms_without_merging_broader_taxa(g):
    mappings, receipts = load_reviewed_mappings("3.2")
    assert receipts and len(mappings) == 5
    species = list(mappings)

    def fetch(url, timeout):
        q = parse_qs(urlparse(url).query)
        name = q["genus"][0] + " " + q["epithet"][0]
        expected = [row for row in mappings.values() if row["decision"] == "map" and row["work_species"] == name]
        assert len(expected) == 1
        row = expected[0]
        return [{"work_ID": row["work_ID"], "work_species": name}]

    report = []
    table = g.resolve_gift_species_map(
        "https://gift.uni-goettingen.de/api/extended/index3.2.php",
        species,
        10,
        fetcher=fetch,
        reviewed_mappings=mappings,
        report=report,
    )
    assert len(table) == 4
    assert dict(zip(table.species, table.work_ID, strict=True))["Dendrobium_catenatum"] == "439741"
    assert next(row for row in report if row["species"] == "Citrus_sinensis")["status"] == "excluded_taxonomic_scope"
    report = []
    table = g.resolve_gift_species_map(
        "https://example.org/index3.2.php",
        ["Dendrobium_catenatum"],
        10,
        fetcher=lambda **kwargs: [{"work_ID": "13683", "work_species": "Dendrobium moniliforme"}],
        reviewed_mappings=mappings,
        report=report,
    )
    assert table.empty and report[0]["status"] == "reviewed_mapping_mismatch"
    assert load_reviewed_mappings("4.0")[0] == {}


def test_ambiguous_hits_are_not_ranked_and_duplicate_labels_are_retained(g):
    report = []
    table = g.resolve_gift_species_map(
        "https://example.org/index3.2.php",
        ["Mimulus_guttatus"],
        10,
        fetcher=lambda **kwargs: [
            {"work_ID": "1", "work_species": "Erythranthe guttata", "overallscore": "1"},
            {"work_ID": "2", "work_species": "Erythranthe nasuta", "overallscore": "0.9"},
        ],
        report=report,
    )
    assert table.empty and "requires_review" in report[0]["status"]
    table = g.resolve_gift_species_map(
        "https://example.org/index3.2.php",
        ["Mimulus_guttatus", "Mimulus_guttatus_genome2"],
        10,
        fetcher=lambda **kwargs: [{"work_ID": "1", "work_species": "Mimulus guttatus"}],
    )
    # Match the repository's recognized accession suffix rather than discarding a work_ID collision.
    assert table.species.tolist() == ["Mimulus_guttatus", "Mimulus_guttatus_genome2"]


def test_mapping_file_validation_and_relative_resolution(g, tmp_path):
    path = tmp_path / "reviewed.tsv"
    path.write_text(
        "species\tgift_version\twork_ID\twork_species\tdecision\tevidence_url\nA_b\t3.2\t1\tC d\tmap\thttps://example.org/evidence\n"
    )
    assert load_reviewed_mappings("3.2", path, include_bundled=False)[0]["A_b"]["work_ID"] == "1"
    source = tmp_path / "sources.tsv"
    source.write_text("database\tacquisition_mode\tgift_species_mapping_file\ngift\tgift_api\treviewed.tsv\n")
    assert g.read_database_sources(source)["gift"]["gift_species_mapping_file"] == str(path)
    path.write_text(path.read_text() + "A_b\t3.2\t2\tE f\tmap\thttps://example.org/evidence\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_reviewed_mappings("3.2", path, include_bundled=False)


def test_missing_binary_and_categorical_sentinels_remain_missing(g):
    frame = pd.DataFrame({"__species_norm": ["a", "b", "c"], "value": ["unknown", None, "NA"]})
    for kind, aggregation in [("binary", "any"), ("categorical", "mode")]:
        plan = g.TraitPlanRow(
            "gift", "value", "trait", kind, aggregation, {"yes"} if kind == "binary" else set(), "", ""
        )
        values = g.aggregate_trait_column(frame, plan).reindex(["a", "b", "c"])
        assert values.isna().all()


def test_workflow_output_cache_tracks_custom_mapping_contents(tmp_path):
    """An edit to a referenced mapping must rerun the real trait stage."""
    import subprocess

    import test_gg_input_generation_end_to_end as fixture

    workspace = tmp_path / "workspace"
    input_dir = fixture._write_direct_species_fixture(tmp_path)
    fixture._write_minimal_ete_taxonomy_db(workspace)
    manifest = fixture._write_tsv_download_manifest(workspace, input_dir)
    fake_bin = fixture._install_fake_toolchain(tmp_path)
    downloads = workspace / "downloads" / "trait_datasets"
    config_dir = workspace / "input" / "input_generation"
    mappings = config_dir / "reviewed.tsv"
    header = "species\tgift_version\twork_ID\twork_species\tdecision\tevidence_url\n"
    mappings.write_text(
        header + "Arabidopsis_thaliana\t3.2\t1\tArabidopsis thaliana\tmap\thttps://example.org/test-only\n"
    )
    sources = config_dir / "sources.tsv"
    sources.write_text(
        "database\tacquisition_mode\turi\tgift_version\tgift_cache_mode\tgift_species_mapping_file\n"
        "gift\tgift_api\thttps://example.org/\t3.2\toffline\treviewed.tsv\n"
    )
    plan = config_dir / "plan.tsv"
    plan.write_text(
        "database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\n"
        "gift\ttrait_value\twoodiness\tcategorical\tmode\t1.1.1\ttrait_ID\n"
    )
    species = {"Arabidopsis thaliana": "1", "Oryza sativa": "2", "Synthetic synonym": "3"}

    def fetch(url, timeout):
        query = parse_qs(urlparse(url).query)
        if query["query"] == ["names_matched_unique"]:
            name = query["genus"][0] + " " + query["epithet"][0]
            return [{"work_ID": species[name], "work_species": name}]
        return [
            {"work_ID": "1", "trait_value": "woody"},
            {"work_ID": "2", "trait_value": "non-woody"},
            {"work_ID": "3", "trait_value": "non-woody"},
        ]

    with GiftRetrieval(downloads / "gift", fetch) as client:
        for name in species:
            genus, epithet = name.split()
            client.fetch(
                f"https://example.org/index3.2.php?query=names_matched_unique&genus={genus}&epithet={epithet}", 10
            )
        client.fetch(
            "https://example.org/index3.2.php?query=traits&traitid=1.1.1&biasref=1&biasderiv=1&startat=0&limit=10000",
            10,
        )
    output = workspace / "input/species_trait/species_trait.tsv"
    env = fixture._core_env(workspace, input_dir, fake_bin, "single")
    env.update({key: "0" for key in env if key.startswith("run_")})
    env.update(
        run_generate_species_trait="1",
        overwrite="0",
        download_manifest=str(manifest),
        trait_plan=str(plan),
        trait_database_sources=str(sources),
        trait_download_dir=str(downloads),
        species_trait_output=str(output),
        strict="1",
    )

    def run():
        result = subprocess.run(
            ["bash", str(fixture.CORE_PATH)],
            cwd=fixture.REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return pd.read_csv(output, sep="\t").set_index("species").woodiness.to_dict()

    assert run()["Arabidopsis_thaliana"] == "woody"
    reports = list((downloads / "gift/runs").glob("*.json"))
    assert run()["Arabidopsis_thaliana"] == "woody"
    assert len(list((downloads / "gift/runs").glob("*.json"))) == len(reports)
    mappings.write_text(
        header + "Arabidopsis_thaliana\t3.2\t3\tSynthetic synonym\tmap\thttps://example.org/test-only\n"
    )
    stale = subprocess.run(
        ["bash", str(fixture.CORE_PATH)], cwd=fixture.REPO_ROOT, env=env, capture_output=True, text=True, timeout=120
    )
    assert stale.returncode == 3 and "gift_custom_mapping_0" in stale.stderr
    assert pd.read_csv(output, sep="\t").set_index("species").loc["Arabidopsis_thaliana", "woodiness"] == "woody"
    env["artifact_stale_policy"] = "rebuild"
    assert run()["Arabidopsis_thaliana"] == "non-woody"
    assert len(list((downloads / "gift/runs").glob("*.json"))) == len(reports) + 1
