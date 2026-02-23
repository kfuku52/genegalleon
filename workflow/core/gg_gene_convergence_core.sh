#!/usr/bin/env bash
set -euo pipefail

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_convergence_entrypoint.sh.
### End: Job-supplied configuration ###

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

# Ensure pymol/csubst plotting works in headless scheduler environments.
export PYMOL_HEADLESS="${PYMOL_HEADLESS:-1}"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

python "${dir_script}/csubst_site_wrapper.py" \
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
