"""Selected decoding preserves existing GFF/GTF attribute and error contracts."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from format_species_annotation.common import parse_gff_attributes


@pytest.mark.parametrize("text", [
    "ID=x;Parent=a,b;Note=ID%3ABtu.g00001t01%2CCDS;unused=one,two;pseudo=%74rue",
    'gene_id "gene;with,commas%20"; Note "ID:Btu.g00001t01"; unused "literal,value";',
    "ID=gene-TIL';Name=gene-TIL%27;Alias=\"balanced value\";gene_id='quoted-gene'",
    "ID=a;ID=b;;Parent=;Note=%3B%2C; unknown = literal;",
    "Lavan.01G000100.1", "simple%20identifier", ".", "", None,
])
@pytest.mark.parametrize("keys", [frozenset(), frozenset({"ID", "Parent", "Note", "pseudo"}),
                                 frozenset({"genome", "gene_id", "absent"})])
def test_projection_is_exact_subset(text, keys):
    full = parse_gff_attributes(text)
    assert parse_gff_attributes(text, keys=keys) == {k: v for k, v in full.items() if k in keys}
    assert parse_gff_attributes(text, keys=None) == full


@pytest.mark.parametrize("text", [
    'Note "good"; unused "unterminated',
    'Note "good"; unused ;',
    'Note "good"; unused "value" extra "value";',
])
@pytest.mark.parametrize("keys", [frozenset(), frozenset({"Note"})])
def test_unused_malformed_gtf_fields_still_fail(text, keys):
    with pytest.raises(ValueError):
        parse_gff_attributes(text)
    with pytest.raises(ValueError):
        parse_gff_attributes(text, keys=keys)


@pytest.mark.parametrize("feature,text,expected", [
    ("CDS", "ID=c;genome_type=%6ditochondrion", True),
    ("region", "ID=c;chromosome=Mt", True),
    ("region", "ID=c;molecule=chloroplast", True),
    ("gene", "ID=g;organelle=apicoplast", True),
    ("CDS", "ID=c;chromosome=Mt", False),
    ("region", "ID=c;Note=genome%3Dmitochondrion;genome=nuclear", False),
    ("region", 'genome "mitochondrion"; unused "literal;percent%20,comma";', True),
])
def test_organelle_scan_projection_preserves_strong_evidence(feature, text, expected):
    from format_species_annotation.organelle import ORGANELLE_SCAN_KEYS, is_organelle_annotation
    assert is_organelle_annotation(feature, parse_gff_attributes(text)) is expected
    assert is_organelle_annotation(feature, parse_gff_attributes(text, keys=ORGANELLE_SCAN_KEYS)) is expected
