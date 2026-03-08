#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual workflow/core/gg_*_core.sh scripts may override these values.
# Host-side entrypoint path bootstrap lives in support/gg_entrypoint_bootstrap.sh
# because this file is also sourced inside the container by core scripts.

: "${GG_COMMON_GENETIC_CODE:=1}" # NCBI genetic code table ID used for translation/ORF-related steps.
: "${GG_COMMON_BUSCO_LINEAGE:=auto}" # BUSCO lineage dataset name, or "auto" to infer a shared dataset from species names.
: "${GG_COMMON_CONTAMINATION_REMOVAL_RANK:=domain}" # Taxonomic rank used when deciding contamination removal.
: "${GG_COMMON_OUTGROUP_LABELS:=Oryza_sativa}" # Comma-separated outgroup species labels for species-tree rooting.
: "${GG_COMMON_ANNOTATION_SPECIES:=auto}" # Annotation reference species, or "auto" to detect a model species from the dataset.
