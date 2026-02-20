#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

run_database_prep=1

### End: Modify this block to tailor your analysis ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

dir_og="${dir_pg_output}/orthogroup"
dir_og_stat_branch="${dir_og}/stat.branch"
dir_og_stat_tree="${dir_og}/stat.tree"
dir_og_csubst_cb_base="${dir_og}/csubst.cb_"

if [[ ${run_database_prep} -eq 1 ]]; then
    missing_input=0
    for required_dir in "${dir_og_stat_tree}" "${dir_og_stat_branch}"; do
      if [[ ! -d "${required_dir}" ]]; then
        echo "Skipping database prep because required directory is missing: ${required_dir}"
        missing_input=1
      fi
    done
    if [[ ${missing_input} -eq 0 ]]; then
      python ${dir_myscript}/generate_orthogroup_database.py \
      --overwrite 1 \
      --dbpath ${dir_og}/gg_orthogroup.db \
      --dir_stat_tree ${dir_og_stat_tree} \
      --dir_stat_branch ${dir_og_stat_branch} \
      --dir_csubst_cb_prefix ${dir_og_csubst_cb_base} \
      --row_threshold 8000 \
      --cutoff_stat "OCNany2spe,0.8" \
      --ncpu ${NSLOTS}
    fi
fi

echo "$(date): Exiting Singularity environment"
