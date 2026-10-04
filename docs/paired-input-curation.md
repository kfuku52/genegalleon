# Explicit curation of a formatted CDS/GFF/genome pair

Input generation rejects GFF references absent from the paired genome. When a
source review chooses to omit those annotations, use the native curation helper
before a new input-generation run. It operates on an already formatted,
representative CDS set and does not select isoforms again or rewrite the genome.

```bash
python workflow/support/curate_paired_species_inputs.py audit \
  --species Genus_species --cds Genus_species_cds.fa.gz \
  --gff Genus_species_annotations.gff.gz --genome Genus_species_genome.fa.gz \
  --report source-review.json
```

Review the reported IDs and references. Supply a JSON decision manifest containing
the **exact three input SHA-256 values** and every individually approved decision:

```json
{
  "schema_version": 1,
  "species": "Genus_species",
  "decision_basis": "Reason and source of the scientific decision",
  "input_sha256": {"cds": "SHA256", "gff": "SHA256", "genome": "SHA256"},
  "exclude_missing_reference_annotations": true,
  "records": [
    {"cds_id": "Genus_species_gene1", "action": "exclude", "reason": "missing_genome_reference"},
    {"cds_id": "Genus_species_gene2", "action": "retain_and_flag", "reason": "missing_gff_counterpart"},
    {"cds_id": "Genus_species_gene3", "action": "retain_and_flag", "reason": "coding_span_conflict",
     "cds_length": 537, "gff_coding_span_length": 769}
  ]
}
```

Examples illustrate the schema; list only decisions observed in the actual pair.
All CDS on absent references must match the approved exclusion set exactly.
CDS without a coding-feature counterpart require an explicit retention flag.
Coding-span flags must match an observed transcript span and the supplied CDS
length. Flags preserve the original CDS and record the unresolved source conflict;
they do not certify nucleotide identity to the genome or repair source coordinates.

```bash
python workflow/support/curate_paired_species_inputs.py curate \
  --species Genus_species --cds Genus_species_cds.fa.gz \
  --gff Genus_species_annotations.gff.gz --genome Genus_species_genome.fa.gz \
  --decision-manifest decisions.json --output-dir new-curated-source
```

The destination must not exist. The helper creates compressed CDS/GFF files and
`curation.json`, binding the decisions, input/output hashes, excluded IDs and
feature counts, retained exceptions and verified remaining reference bounds.
Retained CDS lines and GFF coordinates/phases are preserved. Unique supported
reference aliases are canonicalized; ambiguous aliases, mixed reference owners,
cross-boundary Parents, cycles, unexpected IDs, and changed inputs fail before
publication. Noncoding annotations on an approved missing reference are removed
with that reference. No organelle inference is made from product descriptions.

For subsequent native input generation, put the resulting CDS and GFF paths and
the unchanged genome path in the existing `direct` local-download manifest. Keep
the decision manifest and `curation.json` alongside that manifest. Use a new
prepare → species worker array → finalize chain; changed CDS need fresh BUSCO.
The ordinary formatter and completion checks still apply. Other CDS/GFF conflicts
may remain source warnings; this helper is not a general scientific validator.
