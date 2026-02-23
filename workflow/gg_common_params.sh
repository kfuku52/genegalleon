#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual workflow/core/gg_*_core.sh scripts may override these values.

: "${GG_COMMON_GENETIC_CODE:=1}"
: "${GG_COMMON_BUSCO_LINEAGE:=embryophyta_odb12}"
: "${GG_COMMON_CONTAMINATION_REMOVAL_RANK:=phylum}"
: "${GG_COMMON_OUTGROUP_LABELS:=Oryza_sativa}"
: "${GG_COMMON_ANNOTATION_REPRESENTATIVE_SPECIES:=Arabidopsis_thaliana}"
: "${GG_COMMON_MCMCTREE_DIVERGENCE_TIME_CONSTRAINTS_STR:=Arabidopsis_thaliana,Oryza_sativa,130,-}"
: "${GG_COMMON_TREEVIS_CLADE_ORTHOLOG_PREFIX:=Arabidopsis_thaliana_}"
: "${GG_COMMON_GRAMPA_H1:=}"
: "${GG_COMMON_TARGET_BRANCH_GO:=<1>}"
