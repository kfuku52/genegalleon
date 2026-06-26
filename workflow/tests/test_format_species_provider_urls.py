import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
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


def test_format_species_inputs_delegates_provider_url_helpers():
    script_text = (REPO_ROOT / "workflow" / "support" / "format_species_inputs.py").read_text(encoding="utf-8")

    assert "from format_species_provider_urls import (" in script_text
    assert "def fetch_text_with_headers(" not in script_text
    assert "def strip_provider_prefix(" not in script_text
    assert "ENSEMBL_DEFAULT_ID_URL_TEMPLATES = {" not in script_text
    assert "DEFAULT_NCBI_EUTILS_BASE_URL =" not in script_text
