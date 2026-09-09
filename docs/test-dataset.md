# Test dataset provenance and regeneration

The biological inputs under `workspace/input` retain the original AHA, YABBY,
BUSCO and other test CDS. AHA neighborhoods additionally contain real adjacent
genes, their CDS and complete GFF models, and continuous genome sequence between
them. Existing query lists, expression inputs, RNA-seq inputs and species without
an AHA genomic neighborhood are preserved. Expression values are not invented
for newly added genes.

## Extraction rules

- The 60 genomic AHA IDs are recorded in
  `workspace/input/dataset_manifest/aha_anchors.txt`. They were taken from the
  existing AHA family FASTA; the query lists themselves are unchanged.
- Select up to 20 CDS-bearing gene loci on each genomic side of an AHA locus,
  on the same original chromosome/scaffold. Count loci once, not transcripts.
  Use the CDS bounds and ID ordering from `gff2genestat.py`, as the synteny stage
  does; UTR boundaries can otherwise change the order of overlapping genes.
  Source CDS must contain one representative sequence per gene locus.
- Keep each complete neighborhood as a continuous interval, including long
  intergenic regions. Retain other existing test loci with their own windows.
- Add 5,000 bp of flanking sequence, merge overlapping windows, and include
  complete gene models intersecting the resulting boundaries. Include available
  source CDS for all retained models. Window coordinates are shifted without
  clipping GFF features or altering strand/phase/attributes.
- Never join different scaffolds or fabricate neighbors. A short neighborhood
  reflects an original assembly boundary and is recorded as such.
- Every pre-existing biological CDS must exist in the full source with the same
  sequence (case-insensitive comparison). Its original header, sequence and case
  are retained. Missing or changed source records stop the build.

The complete source is the original seven-species material in the
`20221011_gfe_pipeline/gfe_data` workspace. Exact source filenames and SHA-256
checksums, retained/added IDs, extraction windows, sequence hashes, and unchanged
input hashes are recorded in `dataset_manifest/real_neighborhoods.json`.
`dataset_manifest/aha_coverage.tsv` reports the per-anchor coverage; left/right
refer to genomic coordinate order, independently of transcription strand.

## Current coverage and size

All 134 original biological CDS across the seven test species are retained;
2,026 real CDS were added. The 140 artificial records were moved to the separate
fixture directory. All other input files were preserved byte for byte.

| Species | AHA loci | Both sides reach 20 genes |
| --- | ---: | ---: |
| Amborella trichopoda | 5 | 4 |
| Arabidopsis thaliana | 12 | 12 |
| Cephalotus follicularis | 7 | 0 |
| Dionaea muscipula | 8 | 1 |
| Nepenthes gracilis | 9 | 9 |
| Oryza sativa | 10 | 9 |
| Spinacia oleracea | 9 | 5 |
| Total | 60 | 40 |

The remaining 20 loci are limited by the original scaffolds. All remain in the
dataset. Genome FASTA, GFF and CDS together occupy 34.96 MB uncompressed.

The missing Amborella and Spinacia GFF/genome inputs were added from the matching
original assemblies. The CoGe export had nested mRNAs representing individual
CDS fragments. Ninety-seven chains in the included Cephalotus genes were
normalized only after reconstructed genomic CDS matched the nucleotide input
length and every non-N base. The original CDS sequences and unknown phases were
preserved. `dataset_manifest/coge_source_repair.json` records the initial 92
repairs against the full source; `coge_added_repair.json` records five additional
repairs after window closure. These are explicit data repairs, not runtime
fallbacks for alternative transcript models.

Docker validation found exact agreement with source annotations for all 1,915
AHA synteny neighbor/offset records and source sequence for all 132 genome
windows. All 60 AHA promoters were recovered at 2,000 bp and matched full-source
strand-aware extraction. All 60 selected CDS structures agree in length and
non-N bases with the input CDS. Missing GFF and CDS-length mismatches are zero.
Results are recorded in `dataset_manifest/validation.json`. This validates
these affected stages, not a complete phylogenetic workflow or SIF execution.

## Rebuild without dropping existing test inputs

Run in a GeneGalleon container, mounting the full source read-only if necessary:

```bash
python workflow/support/extend_real_test_dataset.py \
  --source-pg /path/to/full/gfe_data \
  --seed-pg workspace/input \
  --out-pg /path/to/new/test-input \
  --anchor-ids workspace/input/dataset_manifest/aha_anchors.txt \
  --neighbors 20 --flank-bp 5000
```

The output directory must not exist and must be outside the source and seed
directories. The seed is copied before extending the seven selected species.
The original source and seed are never modified. Seed species need CDS but may
lack GFF/genome inputs; those can be added from a matching full source. For the
legacy CoGe material, first use `repair_coge_transcripts.py` on a separate source
copy with the exact CDS/genome and requested CDS IDs (`--gene-ids`). Review its
sequence verification audit; mismatches stop repair without writing outputs. Failed builds have no new
completion manifest and must not be published. Compare all retained IDs and
hashes before copying the three regenerated files per species into
`workspace/input`. Preserve the other inputs and update the dataset manifest
with any intentional future input changes.

The older `build_minimal_test_dataset.py` is useful for a new compact AHA/YABBY/
BUSCO subset, but does not preserve additional seed species or promise 20 real
neighbors. It no longer adds synthetic CDS/GFF annotations and excludes incomplete
boundary models instead of clipping them.

## Synthetic fixtures

The former 40 Cephalotus and 100 Nepenthes artificial CDS/GFF entries are isolated
under `workflow/tests/data/synthetic_synteny`. They are not biological input and
have no matching genome sequence. Their original IDs are listed in the manifest
as `removed_dummy_ids`; only records explicitly marked `gg_dummy_` are removed
from the biological inputs. The extension builder saves these removed records
under the new output's `dataset_manifest/synthetic_synteny` for review and transfer
to that test fixture directory.

## Validation

In the GeneGalleon runtime:

```bash
python -m pytest -q \
  workflow/tests/test_extend_real_test_dataset.py \
  workflow/tests/test_repair_coge_transcripts.py \
  workflow/tests/test_real_aha_dataset.py \
  workflow/tests/test_synteny_search_integration.py \
  workflow/tests/test_get_promoter_fasta_runtime.py
```

The builder tests check preservation, coordinate round trips, continuous windows,
complete boundary models and true scaffold edges. The checked-in dataset tests
also verify the manifest, retained sequence hashes and unchanged inputs.
