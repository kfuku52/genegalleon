# Advisory evidence for rescued models

## Candidate-only Swiss-Prot support

Input generation searches finalized missing-gene rescue CDSs against its existing
Swiss-Prot MMseqs2 database when `run_gene_model_rescue_swissprot=1` (default).
The scoped override is `GG_INPUT_RUN_GENE_MODEL_RESCUE_SWISSPROT`. This runs only
with a completed rescue publication, including one reused by refinement. Set it
to `0` to omit this advisory search. No additional TE database is downloaded.
The normal shared Swiss-Prot FASTA/index/metadata preparation runs on first use.

The audit is separate from the frozen rescue tree, by default its sibling
`gene_model_rescue.swissprot`; `gene_model_rescue_swissprot_dir` overrides it.
Its CLI can annotate an old frozen run without repeating synteny, prediction or BUSCO:

```bash
python workflow/support/rescue_swissprot_evidence.py \
  --rescue-output /data/gene_model_rescue --output /data/gene_model_rescue.swissprot \
  --db-prefix /data/workspace/downloads/uniprot_sprot/uniprot_sprot \
  --metadata /data/workspace/downloads/uniprot_sprot/uniprot_sprot.meta.tsv.gz \
  --cache /data/workspace/downloads/rescue_swissprot_cache --cpus 4 --memory-gb 8
```

Only published rescued coding sequences are translated, using the frozen
species genetic code and ordinary initiator residues. Context-dependent genetic
codes are marked unassessed. Identical proteins share a search across species;
distinct loci and coding-sequence identities remain in the report. The verified
sequence cache queries only new proteins for the same DB, parameters, MMseqs2
binary and annotator. DB preparation and the cache use shared locks.

Search retains up to 50 hits, sensitivity 7.5 and E-value <= 1e-5. Support
requires at least 50 paired residues and 50% paired-residue coverage of both
query and target, and a bit score within 90% of the best qualifying hit. Gap
spans cannot inflate coverage. These are conservative, configurable screening
heuristics, not calibrated probabilities; coverage-sensitive partial matches
remain inspectable. Explicit transposon/transposase/retrotransposon names,
the exact `Transposable element` keyword or `transposase activity` GO term supply
TE-related support. TE silencing/regulation GO terms and the `Transposition`
process keyword alone do not establish TE origin; host methyltransferases,
helicases and other TE-silencing factors remain other proteins. Other informative Swiss-Prot entries supply other-protein
support. Uncharacterized entries without function annotations do not supply
informative support. A generic polymerase/RNase H match is not itself a TE call.

Outputs are `evidence.json`, `loci.tsv`, `hits.tsv`, `summary.json` and a verified
`receipt.json`. The hit table retains accession, annotation, score, coverage and
positions, including hits that did not meet the support rules. An execution
JSON beside the audit records elapsed time and cached/searched query counts.
Each locus counts once as TE-only, other-only, both, no informative support or
not assessed, taking the union of support across its published coding sequences.
Other-protein support does not prove host function. No support does not exclude
a divergent or unrepresented TE. TE-related homology also occurs in domesticated
host genes. The audit never excludes or changes gene models automatically.

The BUSCO comparison accepts `--rescue-swissprot-dir /data/gene_model_rescue.swissprot`.
This changes the lower rescue bar to **Swiss-Prot protein support**, keeping its
S/R/P donor bar and the original BUSCO palette. Existing DNA repeat-overlap
measurements remain a separate evidence axis; protein homology is not labelled
as genomic repeat overlap. Excluded species remain explicitly not analysed.
The figure legend records the actual search/support thresholds from the audit,
the TE/other annotation rules and locus-counting method. Thresholds absent from a
historical summary are labelled unavailable instead of inferred from defaults.

## Genomic, RNA and repeat evidence

Ordinary rescue automatically writes `quality_evidence` in `models.json` and
`quality_flags.tsv`. These report the genomic start triplet, donor N/C terminal
alignment, donor species and strand conflicts. They do not alter ORF/conflict
admission. The standard genetic-code table allows TTG/CTG initiation, but this
does not establish translation initiation in a particular target gene; see
[NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Alignment to both donor ends likewise does not establish native completeness,
particularly if the donor annotation itself is partial.
`query_alignment` reports aligned donor residues separately from the bounding
query span and the fraction unaligned inside that span. A global alignment of a
short fragment can pair both donor end residues while leaving a large internal
deletion. Read N/C flags with these fractions; they do not mean a complete gene.
The fractions are also appended to both TSV outputs. The additional
`donor_internal_unaligned_query` flag is advisory; admission thresholds are unchanged.

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
Parsed plan, receipt and manifest metadata are bound to the bytes actually read;
an update between parsing, verification and publication fails instead of
rebinding old parsed values to a newer file hash. Duplicate evidence attribute
keys, truncated escaped GTF identifiers, and empty FASTA records are refused or
parsed without losing identity.
Ordinary rescue likewise compares its in-memory plan with the frozen plan before
stamping preparation, rescue, export, worker completion or QC receipts. A changed
plan fails before results can be labelled with the replacement plan's hash.

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

For repeat annotation, reuse an existing RepeatMasker `.out` only when it belongs
to the exact frozen genome. For unannotated plant genomes, a species-specific
library from [EDTA](https://github.com/oushujun/EDTA) followed by
[RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) is a practical
route. EDTA accepts trusted coding sequences to remove gene contamination from
the repeat library. Retain simple/low-complexity repeats as well as TE classes
when the report needs both axes. GeneGalleon currently consumes RepeatMasker
output; it does not launch EDTA or RepeatMasker, or treat softmasking alone as
a classified TE annotation.

To assess whether an overlapping coding model itself encodes a TE protein,
inspect TE protein domains separately, for example with
[TEsorter/REXdb](https://github.com/zhangrengang/TEsorter), along with domain
coverage, RNA structure, synteny and conserved non-TE gene homology. An isolated
domain or repeat overlap is advisory; a domain-negative result cannot exclude
non-autonomous or divergent TEs. The repeat bar in the
[BUSCO/model-change comparison](gene-model-refinement.md) counts any CDS overlap,
with TE hits taking priority over other/unclassified repeats. The existing
`te_overlap_ge_50pct` flag remains a separate high-overlap review flag.

DNA coverage excludes unmapped, secondary, supplementary, QC-failed and duplicate
reads, excludes MAPQ 255 (mapping quality unavailable, as defined by the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)), and applies
the declared mapping/base quality thresholds. It records CDS
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

Swiss-Prot raw alignment and target-annotation caches are independent. Raw
alignment identity binds the exact query protein, MMseqs binary/database,
alignment parser, search E-value, sensitivity and maximum hits. Metadata or
classification changes do not discard raw searches. Target annotations bind
the FASTA/metadata release and annotation classifier. Checksummed SQLite records
and batched cursors avoid per-accession network lock operations.
Coverage, length and competing-score changes therefore only reassess hits.
When omitted, `--search-evalue` uses the larger of `1e-5` and the support
E-value, preserving older calls that only set a looser support cutoff.
An explicitly set search bound must include the support E-value. A tighter
support cutoff alone keeps the default search bound and cache. Actual search parameter,
sequence or database changes correctly require new searches.

Primary TE/other categories remain unchanged. A separate
`partial_te_homology` flag marks any returned TE-labelled hit passing support
E-value, paired-residue and query-coverage thresholds but failing target
coverage, without a competing-score filter. It does not prove TE origin.
`rescue_swissprot_diagnostics.png/svg` distinguishes whole-protein TE support,
partial flags alone, no TE support/flag and unassessed translation, and displays
excluded species. Every locus counts once.

`no_support_reason` partitions only no-informative-support loci: no returned
hit; all E-values fail; E/coverage pass but paired length fails; coverage fails;
or qualifying best-score annotation unknown. Across coding sequences, priority
is annotation unknown, partial, short, weak, no hit. Thresholds and priority
appear in the diagnostic legend. Short proteins retain the same minimum
alignment length; separate counts permit calibration without silent relaxation.
