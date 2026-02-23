from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import types


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_statistics.py"


def load_module():
    if "kftools" not in sys.modules:
        pkg = types.ModuleType("kftools")
        pkg.__path__ = []
        sys.modules["kftools"] = pkg
    if "kftools.kfog" not in sys.modules:
        sys.modules["kftools.kfog"] = types.ModuleType("kftools.kfog")
    if "kftools.kfphylo" not in sys.modules:
        sys.modules["kftools.kfphylo"] = types.ModuleType("kftools.kfphylo")
    spec = spec_from_file_location("orthogroup_statistics", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_fimo_hits_parses_modern_fimo_tsv(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.tsv"
    infile.write_text(
        (
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
            "M1\tALT1\tgeneA\t10\t15\t+\t11.5\t1e-4\t1e-3\tACGTAC\n"
            "M2\t.\tgeneB\t20\t24\t-\t9.3\t1e-3\t2e-2\tTTGGCC\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert df.columns.tolist() == [
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "q-value",
    ]
    assert df.loc[0, "motif_alt_id"] == "ALT1"
    # "." should fall back to motif_id for consistent downstream grouping.
    assert df.loc[1, "motif_alt_id"] == "M2"


def test_load_fimo_hits_parses_legacy_fimo_txt_header(tmp_path):
    mod = load_module()
    infile = tmp_path / "fimo.txt"
    infile.write_text(
        (
            "#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence\n"
            "MA0123.1\tgeneC\t30\t39\t+\t7.0\t1e-3\t4e-2\tACGTTGCAAA\n"
            "MA0567.1\tgeneD\t40\t47\t-\t6.1\t2e-3\t5e-2\tTTAACCGG\n"
        ),
        encoding="utf-8",
    )

    df = mod.load_fimo_hits(str(infile))
    assert df.shape[0] == 2
    assert set(df["sequence_name"].tolist()) == {"geneC", "geneD"}
    # legacy format has no alt id; should be mirrored from motif_id.
    assert df.loc[df["sequence_name"] == "geneC", "motif_alt_id"].iloc[0] == "MA0123.1"
