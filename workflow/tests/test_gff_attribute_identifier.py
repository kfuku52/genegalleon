import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
SCRIPT_PATH = SUPPORT_DIR / "gff_attribute_identifier.py"


def load_module():
    sys.path.insert(0, str(SUPPORT_DIR))
    try:
        spec = spec_from_file_location("gff_attribute_identifier", SCRIPT_PATH)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        try:
            sys.path.remove(str(SUPPORT_DIR))
        except ValueError:
            pass


@pytest.mark.parametrize(
    "species",
    [
        "Oryza_sativa_subsp_japonica",
        "Cenchrus_americanus_x_Cenchrus_purpureus",
    ],
)
def test_identify_attribute_key_strips_species_prefix(tmp_path, species):
    mod = load_module()
    fasta = tmp_path / f"{species}.cds.fa"
    gff = tmp_path / f"{species}.gff3"
    fasta.write_text(f">{species}_gene1\nATG\n", encoding="utf-8")
    gff.write_text("chr1\tsrc\tmRNA\t1\t3\t.\t+\t.\tID=gene1\n", encoding="utf-8")

    args = type("Args", (), {"fasta_file": str(fasta), "gff_file": str(gff), "gff_feature": "mRNA"})()

    assert mod.identify_attribute_key(args) == "ID"
