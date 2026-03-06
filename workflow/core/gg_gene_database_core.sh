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
# Configuration variables are provided by gg_gene_database_entrypoint.sh.
### End: Job-supplied configuration ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

file_orthogroup_db="${gg_workspace_output_dir}/orthogroup/gg_orthogroup.db"

enable_all_run_flags_for_debug_mode

if [[ ${run_database_prep} -eq 1 ]]; then
    missing_input=0
    for required_dir in \
      "$(dirname "${file_orthogroup_db}")/stat_tree" \
      "$(dirname "${file_orthogroup_db}")/stat_branch"
    do
      if [[ ! -d "${required_dir}" ]]; then
        echo "Skipping database prep because required directory is missing: ${required_dir}"
        missing_input=1
      fi
    done
    if [[ ${missing_input} -eq 0 ]]; then
	      python "${gg_support_dir}/generate_orthogroup_database.py" \
	      --overwrite 1 \
	      --dbpath "${file_orthogroup_db}" \
	      --dir_stat_tree "$(dirname "${file_orthogroup_db}")/stat_tree" \
	      --dir_stat_branch "$(dirname "${file_orthogroup_db}")/stat_branch" \
	      --dir_csubst_cb_prefix "$(dirname "${file_orthogroup_db}")/csubst_cb_" \
	      --row_threshold 8000 \
	      --cutoff_stat "OCNany2spe,0.8" \
	      --ncpu "${NSLOTS}"
	    fi
	fi

echo "$(date): Exiting Singularity environment"
