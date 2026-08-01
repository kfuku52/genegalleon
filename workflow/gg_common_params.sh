#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual workflow/core/gg_*_core.sh scripts may override these values.
# Host-side entrypoint path bootstrap lives in support/gg_entrypoint_bootstrap.sh
# because this file is also sourced inside the container by core scripts.

: "${GG_COMMON_GENETIC_CODE:=1}" # NCBI genetic code table ID used for translation/ORF-related steps.
: "${GG_COMMON_BUSCO_LINEAGE:=eukaryota_odb12}" # Shared BUSCO lineage dataset default; override per workflow when a narrower lineage or auto inference is needed.
: "${GG_COMMON_REFERENCE_SPECIES:=auto}" # Reference species, or "auto" to detect a model species from the dataset.
: "${GG_COMMON_INPUT_SEQUENCE_MODE:=cds}" # Shared default input sequence mode for workflows that consume CDS/protein inputs.
: "${GG_COMMON_CSUBST_NONSYN_RECODE:=no}" # Shared CSUBST nonsynonymous-state recoding scheme for search and site-level summaries.
: "${GG_COMMON_SPECIES_LABEL_PARSER:=taxonomic}" # Shared species-label parser used when downstream tools need to recover species names from gene/tree labels.
: "${GG_COMMON_SPECIES_LABEL_REGEX:=}" # Optional regex for downstream tools that support parser-driven species extraction from nonstandard labels.
: "${GG_COMMON_SPECIES_LABEL_MAP_TSV:=}" # Optional mapping table for downstream tools that support parser-driven species extraction.
: "${GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE:=zip}" # zip|files|raw; raw is an alias for the historical files layout, while ZIP mode uses transparent per-subdirectory shards.
: "${GG_COMMON_GENE_FAMILY_ZIP_MIN_BATCH_FILES:=100}" # Minimum live files per output subdirectory before an array task opportunistically creates a ZIP shard; progress summary flushes smaller completed batches.
: "${GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS:=7}" # Failed gene-family task directories below output/{query2family,orthogroup}/tmp are retained for this many days; 0 disables the age limit.
: "${GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS:=100}" # Maximum failed gene-family task directories retained per query2family/orthogroup root; oldest excess directories are removed and 0 disables the count limit.
: "${GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES:=107374182400}" # Maximum total bytes retained across inactive failed gene-family tmp directories (100 GiB by default); 0 disables the byte limit.
: "${GG_COMMON_GENE_FAMILY_TMP_MAX_FILES:=100000}" # Maximum total files retained across inactive failed gene-family tmp directories; 0 disables the file-count limit.
