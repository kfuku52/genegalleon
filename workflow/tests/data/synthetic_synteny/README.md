# Historical artificial synteny fixtures

These are the 40 Cephalotus and 100 Nepenthes dummy CDS records and GFF rows
formerly appended to the biological development dataset for YABBY plot examples.
They are retained here for isolated synthetic tests and reproducing those examples.

The sequences are artificial. The GFF coordinates refer to the old extracted
windows and do **not** have matching genomic sequence. These files must never be
copied into `workspace/input/species_cds` or `species_gff`, used to validate promoter
extraction, or treated as real evidence of synteny.

The runtime neighborhood tests in `test_synteny_search_integration.py` generate
their own synthetic complete neighborhoods, and additionally test the biological
AHA dataset against the source-derived neighborhood manifest.
