#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual workflow/core/gg_*_core.sh scripts may override these values.
# Host-side entrypoint path bootstrap lives in support/gg_entrypoint_bootstrap.sh
# because this file is also sourced inside the container by core scripts.

: "${GG_COMMON_GENETIC_CODE:=1}" # NCBI genetic code table ID used for translation/ORF-related steps.
: "${GG_COMMON_BUSCO_LINEAGE:=auto}" # BUSCO lineage dataset name, or "auto" to infer a shared dataset from species names.
: "${GG_COMMON_REFERENCE_SPECIES:=auto}" # Reference species, or "auto" to detect a model species from the dataset.
: "${GG_COMMON_CONTAMINATION_REMOVAL_RANK:=domain}" # Default taxonomic rank for contamination-removal workflows.
: "${GG_COMMON_CONTAMINATION_REMOVAL_TARGET_TAXON:=}" # Optional lineage anchor shared by contamination-removal workflows.
