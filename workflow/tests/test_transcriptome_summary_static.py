from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
R_SCRIPT = REPO_ROOT / "workflow" / "support" / "multispecies_transcriptome_summary.r"
BASE_R_REQUIRED = REPO_ROOT / "container" / "env" / "base.r.required.txt"


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_base_r_required_includes_rtsne():
    lines = {
        line.strip()
        for line in _read_text(BASE_R_REQUIRED).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    assert "r-rtsne" in lines


def test_multispecies_transcriptome_summary_filters_before_imputation_and_skips_optional_methods():
    text = _read_text(R_SCRIPT)

    assert 'has_rtsne = requireNamespace("Rtsne", quietly=TRUE)' in text
    assert "filter_expression_for_dimensional_reduction = function(df_in, min_non_missing = 2)" in text
    assert "Converting %d Inf/-Inf value(s) to NA before dimensional reduction." in text
    assert "Removing %d all-NA gene row(s) before dimensional reduction." in text
    assert "Removing %d all-NA species column(s) before dimensional reduction." in text
    assert "Removing %d zero-variance gene row(s) before dimensional reduction." in text
    assert "Removing %d zero-variance species column(s) before dimensional reduction." in text
    assert "Removing %d gene row(s) with <%d observed species before dimensional reduction." in text
    assert "Removing %d species column(s) with <%d observed genes before dimensional reduction." in text
    assert "df_exp_filtered = filter_expression_for_dimensional_reduction(df_exp, min_non_missing=min_species)" in text
    assert "missMDA::estim_ncpPCA(df_exp_filtered" in text
    assert "missMDA::imputePCA(df_exp_filtered" in text
    assert "Skipping dimensional reduction because too few rows/columns remain after filtering" in text
    assert "Skipping PCA plot because at least 2 gene rows and 3 species columns are required after filtering" in text
    assert "Skipping MDS plot because at least 3 species columns are required after filtering" in text
    assert "Skipping tSNE plot because Rtsne is unavailable." in text
    assert "perplexity = min(30, floor((ncol(input_data) - 1) / 3))" in text
