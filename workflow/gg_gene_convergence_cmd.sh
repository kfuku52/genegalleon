#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

arity_range='2-10'
trait='all'
skip_lower_order='yes'
min_fg_stem_ratio='0.5'
min_OCNany2spe='1.8'
min_omegaCany2spe='3.0'
min_OCNCoD='0'
max_per_K='100'
file_trait='/workspace/input/species_trait/species_trait.tsv'
dir_orthogroup='/workspace/output/orthogroup'
dir_orthofinder='/workspace/output/orthofinder'
dir_out='/workspace/output/csubst_site'

### End: Modify this block to tailor your analysis ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

# Ensure pymol/csubst plotting works in headless scheduler environments.
export PYMOL_HEADLESS="${PYMOL_HEADLESS:-1}"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

python "${dir_myscript}/csubst_site_wrapper.py" \
--arity_range "${arity_range}" \
--trait "${trait}" \
--skip_lower_order "${skip_lower_order}" \
--min_fg_stem_ratio "${min_fg_stem_ratio}" \
--min_OCNany2spe "${min_OCNany2spe}" \
--min_omegaCany2spe "${min_omegaCany2spe}" \
--min_OCNCoD "${min_OCNCoD}" \
--max_per_K "${max_per_K}" \
--ncpu "${NSLOTS:-1}" \
--file_trait "${file_trait}" \
--dir_orthogroup "${dir_orthogroup}" \
--dir_orthofinder "${dir_orthofinder}" \
--dir_out "${dir_out}"

echo "$(date): Exiting Singularity environment"
