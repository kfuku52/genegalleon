# Advisory evidence for rescued models

Ordinary rescue automatically writes `quality_evidence` in `models.json` and
`quality_flags.tsv`. These report the genomic start triplet, donor N/C terminal
alignment, donor species and strand conflicts. They do not alter ORF/conflict
admission. The standard genetic-code table allows TTG/CTG initiation, but this
does not establish translation initiation in a particular target gene; see
[NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Alignment to both donor ends likewise does not establish native completeness,
particularly if the donor annotation itself is partial.

For existing frozen runs, a separate CLI regenerates the evidence report without
resuming an old plan under new tool identities or repeating synteny/prediction:

```bash
python workflow/support/rescue_model_evidence.py \
  --rescue-output /absolute/path/to/frozen/gene_model_rescue \
  --species Genus_species \
  --evidence-manifest /absolute/path/to/evidence.json \
  --output /absolute/path/to/separate/model_evidence/Genus_species
```

`--evidence-manifest` is optional. Omitted evidence has `not_provided` status
and null measurements, distinct from supplied evidence with no detected support.
The CLI binds its implementation, model and plan receipts, exact genome, tracks,
BAM index and package versions to a separate completion receipt. Changed inputs
invalidate this report; altered frozen models/genome fail. All accepted CDSs are
reconstructed from the genome and checked for phase, sequence and ORF agreement.
The large `models.json` array is streamed. Original sequence/annotation, accepted
IDs, admission status and BUSCO results are never changed by this report.

## Optional manifest

Use absolute paths and the SHA-256 of the exact genome file frozen in the rescue
plan. This is an explicit assertion of track reference provenance; matching BAM
contig lengths alone cannot prove that the underlying bases are identical.

```json
{
  "schema_version": 1,
  "species": {
    "Genus_species": {
      "reference_genome_sha256": "SHA256_OF_FROZEN_GENOME_FILE",
      "rna_junctions": [
        {"path": "/data/hintsfile.gff", "format": "braker_hints",
         "independence_group": "individual1_leaf"}
      ],
      "rna_transcripts": [
        {"path": "/data/stringtie.gtf", "format": "exon_gff_gtf",
         "independence_group": "individual1_leaf"}
      ],
      "repeats": [
        {"path": "/data/genome.fa.out", "format": "repeatmasker_out"}
      ],
      "dna": {
        "path": "/data/hifi.bam", "index": "/data/hifi.bam.bai",
        "format": "bam", "min_mapq": 20, "min_baseq": 20
      }
    }
  }
}
```

All track sections are optional. Technical RNA runs from the same sample should
share an `independence_group`. RNA exon paths are grouped by transcript identity;
support from disjoint transcripts is never joined into a whole path. Unknown
strand paths remain separate from stranded support. BRAKER hints use only
`intron` rows with `src=E`; protein hints are not counted as RNA.

An RNA exon chain compatible with the CDS supports transcriptional structure.
It does not prove that the CDS start is translated. It is also distinct from the
reviewed exact RNA coding paths consumed by gene-model refinement. RNA absence
in a sampled tissue is not evidence of gene absence or a reason for automatic
rejection.

Repeat overlap is computed on the union of masked CDS intervals; overlapping
mask records are not counted twice. TE-labelled LINE/SINE/LTR/DNA/RC/Retroposon/PLE
classes are reported separately from other repeats. High TE overlap is a review
flag, not a rejection rule. Different loci encoding identical proteins remain
different models. No protein-sequence deduplication is applied to this report.

DNA coverage excludes unmapped, secondary, supplementary, QC-failed and duplicate
reads, and applies the declared mapping/base quality thresholds. It records CDS
minimum/median depth and per-base support for the first spliced codon, including
minus-strand/split codons. DNA support confirms sequence evidence, not expression,
translation initiation, secretion or enzyme function.

A BAM mapped before contig filtering may contain extra references. Provide
`dna.reference_genome` with the original mapping FASTA to use such a BAM. Its
contig names/lengths must exactly match the BAM header, and every retained target
contig must match the frozen genome base for base (case-insensitive). Missing,
resized or changed target contigs fail. The reference FASTA is hashed in the
receipt; the summary records excluded BAM contigs. Extra references without this
verification fail. Mapping quality still reflects the original, larger reference.

Outputs are `evidence.json`, `evidence.tsv`, `summary.json` and `receipt.json`.
They have independent evidence axes and advisory flags, with no calibrated
confidence probability or automatic change to the representative gene set.

## Reassessment and cache scope

Use start/terminal/copy conflicts to nominate bounded locus reassessment. Compare
target RNA structure and donor reliability as well as protein alignments; retain
alternative models when the data cannot decide. Split/merge, a non-ATG start,
repeat overlap and missing tissue RNA are not sufficient automatic adoption or
exclusion criteria.

An evidence-only refresh needs no BUSCO or new synteny. A changed coding sequence
or model membership requires a new effective CDS/GFF view and QC for changed
species, with initial BUSCO lineage/version/date/marker count preserved.
New rescue plans can reuse comparisons only through the verified comparison
cache when BED/PEP content, parameters and comparison tool identities agree.
Never edit an old plan's implementation hash or relabel its receipts.
