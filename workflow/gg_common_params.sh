#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual workflow/core/gg_*_core.sh scripts may override these values.
# Host-side entrypoint path bootstrap lives in support/gg_entrypoint_bootstrap.sh
# because this file is also sourced inside the container by core scripts.

: "${GG_COMMON_GENETIC_CODE:=1}" # NCBI genetic code table ID used for translation/ORF-related steps.
: "${GG_COMMON_BUSCO_LINEAGE:=eukaryota_odb12}" # Shared BUSCO lineage dataset default; override per workflow when a narrower lineage or auto inference is needed.
: "${GG_COMMON_REFERENCE_SPECIES:=auto}" # Reference species, or "auto" to detect a model species from the dataset.
: "${GG_COMMON_INPUT_SEQUENCE_MODE:=cds}" # Shared default input sequence mode for workflows that consume CDS/protein inputs.
: "${GG_COMMON_SPECIES_LABEL_PARSER:=taxonomic}" # Shared species-label parser used when downstream tools need to recover species names from gene/tree labels.
: "${GG_COMMON_SPECIES_LABEL_REGEX:=}" # Optional regex for downstream tools that support parser-driven species extraction from nonstandard labels.
: "${GG_COMMON_SPECIES_LABEL_MAP_TSV:=}" # Optional mapping table for downstream tools that support parser-driven species extraction.
