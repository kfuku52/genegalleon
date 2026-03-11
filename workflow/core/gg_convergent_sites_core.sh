#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_convergent_sites_entrypoint.sh.
### End: Job-supplied configuration ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

# Ensure pymol/csubst plotting works in headless scheduler environments.
export PYMOL_HEADLESS="${PYMOL_HEADLESS:-1}"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

if [[ "${file_trait}" == "auto" ]]; then
  file_trait="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
fi
if [[ "${dir_orthogroup}" == "auto" ]]; then
  dir_orthogroup="${gg_workspace_output_dir}/orthogroup"
fi
if [[ "${dir_orthofinder}" == "auto" ]]; then
  dir_orthofinder="${gg_workspace_output_dir}/orthofinder"
fi
if [[ "${dir_out}" == "auto" ]]; then
  dir_out="${gg_workspace_output_dir}/csubst_site"
fi

python "${gg_support_dir}/csubst_site_wrapper.py" \
--arity_range "${arity_range}" \
--trait "${trait}" \
--skip_lower_order "${skip_lower_order}" \
--min_fg_stem_ratio "${min_fg_stem_ratio}" \
--min_OCNany2spe "${min_OCNany2spe}" \
--min_omegaCany2spe "${min_omegaCany2spe}" \
--min_OCNCoD "${min_OCNCoD}" \
--max_per_K "${max_per_K}" \
--ncpu "${GG_TASK_CPUS:-1}" \
--file_trait "${file_trait}" \
--dir_orthogroup "${dir_orthogroup}" \
--dir_orthofinder "${dir_orthofinder}" \
--dir_out "${dir_out}"

echo "$(date): Exiting Singularity environment"
