import sys
from pathlib import Path
from urllib.parse import urlparse

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from format_species_provider_urls import (  # noqa: E402
    ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES,
    NCBI_ASSEMBLY_ACCESSION_PATTERN,
    append_cookie_header,
    extract_ncbi_accession_from_source_id,
    render_id_url_template,
    strip_provider_prefix,
)
from format_species_providers.common import infer_provider_from_id  # noqa: E402
from format_species_providers.ncbi import normalize_ncbi_ftp_path  # noqa: E402

PROVIDER_URLS = [
    ("https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/", "ncbi"),
    ("ftp://ftp.ensembl.org/pub/release-110/", "ensembl"),
    ("https://www.ensembl.org/pub/current_fasta/", "ensembl"),
    ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/", "ensemblplants"),
    ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/", "ensemblplants"),
    ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/", "ensemblmetazoa"),
    ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/gff3/", "ensemblprotists"),
    ("https://genomevolution.org/coge/GenomeInfo.pl?gid=24739", "coge"),
    ("https://ngdc.cncb.ac.cn/gwh/Assembly/123/show", "gwh"),
    ("https://www.citrusgenomedb.org/organism/5799", "citrusgenomedb"),
    ("https://figshare.com/articles/dataset/test/12345", "figshare"),
    ("https://api.figshare.com/v2/articles/12345", "figshare"),
    ("https://plantgarden.jp/en/list/t64480/genome/t64480.G001", "plantgarden"),
    ("https://plantaedb.com/taxa/example", "plantaedb"),
    ("https://db.cngb.org/cnsa/assembly/CNA000001", "cngb"),
    ("https://download.cncb.ac.cn/other/", "cngb"),
    ("https://ddbj.nig.ac.jp/search/entry/bioproject/PRJDB15739", "ddbj"),
    ("https://getentry.ddbj.nig.ac.jp/getentry/na/BAAHMP000000000", "ddbj"),
    ("https://ftp.flybase.org/genomes/", "flybase"),
    ("https://downloads.wormbase.org/releases/", "wormbase"),
    ("https://vectorbase.org/common/downloads/", "vectorbase"),
    ("ftp://ftp.fernbase.org/genomes/", "fernbase"),
    ("https://veupathdb.org/veupathdb/service/", "veupathdb"),
    ("https://dictybase.org/downloads/", "dictybase"),
    ("https://www.insect-genome.com/genome/IBG_00001", "insectbase"),
]


@pytest.mark.parametrize("url,provider", PROVIDER_URLS)
def test_infer_provider_preserves_official_urls_and_dns_spelling(url, provider):
    assert infer_provider_from_id(url) == provider
    parsed = urlparse(url)
    alternate = parsed._replace(netloc=parsed.netloc.upper() + ".:8443").geturl()
    assert infer_provider_from_id(alternate) == provider


@pytest.mark.parametrize("url,provider", PROVIDER_URLS)
@pytest.mark.parametrize("placement", ["host_suffix", "userinfo", "path", "query", "fragment"])
def test_infer_provider_does_not_trust_provider_names_outside_hostname(url, provider, placement):
    parsed = urlparse(url)
    if placement == "host_suffix":
        candidate = parsed._replace(netloc=parsed.netloc + ".example.invalid").geturl()
    elif placement == "userinfo":
        candidate = parsed._replace(netloc=parsed.netloc + "@example.invalid").geturl()
    elif placement == "path":
        candidate = "https://example.invalid/" + parsed.netloc + parsed.path
    elif placement == "query":
        candidate = "https://example.invalid/?next=" + url
    else:
        candidate = "https://example.invalid/#" + url
    assert infer_provider_from_id(candidate) == ""


@pytest.mark.parametrize("domain", [
    "ncbi.nlm.nih.gov", "ensembl.org", "genomevolution.org", "cncb.ac.cn", "cngb.org",
    "citrusgenomedb.org", "figshare.com", "plantgarden.jp", "plantaedb.com", "ddbj.nig.ac.jp",
    "flybase.org", "wormbase.org", "vectorbase.org", "fernbase.org", "veupathdb.org",
    "dictybase.org", "insect-genome.com",
])
def test_infer_provider_requires_domain_label_boundary(domain):
    assert infer_provider_from_id("https://not" + domain + "/pub/current_fasta/") == ""


@pytest.mark.parametrize("source_id,provider", [
    ("GCF_000001405.40", "ncbi"), ("GCA_000001635.9", "ncbi"),
    ("ddbj:BAAHMP000000000", "ddbj"), ("PRJDB15739", "ddbj"), ("BAAHMP000000000", "ddbj"),
    ("ensembl:homo_sapiens", "ensembl"), ("ensemblmetazoa:anopheles_gambiae", "ensemblmetazoa"),
    ("ensemblprotists:phytophthora_parasitica", "ensemblprotists"),
    ("fernbase:example", "fernbase"), ("veupathdb:example", "veupathdb"),
    ("dictybase:example", "dictybase"), ("insectbase:IBG_00001", "insectbase"),
    ("IBG_00001", "insectbase"), ("oryza_minuta:gramene_tetraploids", "oryza_minuta"),
    ("coge:24739", "coge"), ("gwh:GWHIGRM00000000.1", "gwh"),
    ("GWHIGRM00000000.1", "gwh"), ("citrusgenomedb:5799", "citrusgenomedb"),
    ("figshare:12345", "figshare"), ("plantgarden:t64480", "plantgarden"),
    ("plantaedb:example", "plantaedb"), ("cngb:CNA000001", "cngb"), ("CNA000001", "cngb"),
    ("", ""), ("unknown", ""), ("file:///tmp/flybase.org/data.fa", ""),
    ("https://example.invalid/?project=PRJDB15739", ""), ("https://example.invalid/gwh/test", ""),
    ("https://[malformed", ""), ("https://example.invalid/BAAHMP000000000", ""),
])
def test_infer_provider_preserves_ids_without_guessing_from_unrelated_urls(source_id, provider):
    assert infer_provider_from_id(source_id) == provider


@pytest.mark.parametrize("scheme,port", [("ftp", ""), ("http", ""), ("https", ""),
                                         ("ftp", ":21"), ("http", ":80"), ("https", ":443")])
def test_ncbi_normalization_rebases_only_official_origin_and_preserves_url_components(monkeypatch, scheme, port):
    monkeypatch.setenv("GG_NCBI_FTP_BASE_URL", "http://127.0.0.1:8123/mirror")
    tail = "/genomes/all/GCA_123.1_ASM%20test;param?download=1#part"
    assert normalize_ncbi_ftp_path(scheme + "://FTP.NCBI.NLM.NIH.GOV." + port + tail) == (
        "http://127.0.0.1:8123/mirror" + tail
    )


@pytest.mark.parametrize("url,expected", [
    ("https://ftp.ncbi.nlm.nih.gov.example.invalid/genomes/all/test", None),
    ("http://ftp.ncbi.nlm.nih.gov.example.invalid/genomes/all/test", None),
    ("https://ftp.ncbi.nlm.nih.gov@example.invalid/genomes/all/test", None),
    ("https://ftp.ncbi.nlm.nih.gov:8443/genomes/all/test", None),
    ("https://user:password@ftp.ncbi.nlm.nih.gov/genomes/all/test", None),
    ("https://example.invalid/ftp://ftp.ncbi.nlm.nih.gov/genomes/all/test", None),
    ("ftp://ftp.ncbi.nlm.nih.gov.example.invalid/genomes/all/test",
     "https://ftp.ncbi.nlm.nih.gov.example.invalid/genomes/all/test"),
    ("ftp://ftp.ncbi.nlm.nih.gov@example.invalid/genomes/all/test",
     "https://ftp.ncbi.nlm.nih.gov@example.invalid/genomes/all/test"),
    ("ftp://example.invalid/genomes/test", "https://example.invalid/genomes/test"),
    ("file:///tmp/genomes/test", None), ("", None),
])
def test_ncbi_normalization_keeps_custom_origins_separate_from_mirror(monkeypatch, url, expected):
    monkeypatch.setenv("GG_NCBI_FTP_BASE_URL", "http://127.0.0.1:8123/mirror")
    assert normalize_ncbi_ftp_path(url) == (url if expected is None else expected)


def test_ncbi_normalization_preserves_default_and_cannot_replace_mirror_with_resource_host(monkeypatch):
    monkeypatch.delenv("GG_NCBI_FTP_BASE_URL", raising=False)
    assert normalize_ncbi_ftp_path("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/test") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/test"
    )
    monkeypatch.setenv("GG_NCBI_FTP_BASE_URL", "http://127.0.0.1:8123/mirror")
    assert normalize_ncbi_ftp_path("https://ftp.ncbi.nlm.nih.gov//example.invalid/test") == (
        "http://127.0.0.1:8123/mirror//example.invalid/test"
    )


def test_provider_url_helpers_preserve_source_id_normalization():
    assert strip_provider_prefix("figshare:12345", "figshare") == "12345"
    assert strip_provider_prefix("FigShare:12345", "figshare") == "12345"
    assert strip_provider_prefix("other:12345", "figshare") == "other:12345"

    assert extract_ncbi_accession_from_source_id("GCF_000001405.40") == "GCF_000001405.40"
    assert extract_ncbi_accession_from_source_id("https://example.org/GCA_000001635.9/") == "GCA_000001635.9"
    assert NCBI_ASSEMBLY_ACCESSION_PATTERN.match("GCF_000001405.40")


def test_provider_url_template_and_cookie_helpers():
    assert (
        render_id_url_template(
            "https://example.org/{provider}/{id_lower}/{species_key}",
            "ensemblplants",
            "Ostreococcus_lucimarinus",
            "Ostreococcus_lucimarinus",
        )
        == "https://example.org/ensemblplants/ostreococcus_lucimarinus/Ostreococcus_lucimarinus"
    )
    assert append_cookie_header("a=1", "b=2") == "a=1; b=2"
    assert "ensemblplants" in ENSEMBLGENOMES_DEFAULT_ID_URL_TEMPLATES
