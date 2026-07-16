#!/usr/bin/env bash

# Compatibility entry point for the modular GeneGalleon shell helper library.
gg_util_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd 2>/dev/null || true)"
GG_UTIL_SUPPORT_DIR="${gg_util_dir}"
if [[ -n "${gg_util_dir}" && -s "${gg_util_dir}/gg_site_runtime.sh" ]]; then
  # shellcheck disable=SC1090
  source "${gg_util_dir}/gg_site_runtime.sh"
fi
if [[ -n "${gg_util_dir}" && -s "${gg_util_dir}/gg_entrypoint_config_vars.sh" ]]; then
  # shellcheck disable=SC1090
  source "${gg_util_dir}/gg_entrypoint_config_vars.sh"
fi
if [[ -n "${gg_util_dir}" && -s "${gg_util_dir}/gg_shared_lock.sh" ]]; then
  # shellcheck disable=SC1090
  source "${gg_util_dir}/gg_shared_lock.sh"
fi
for gg_util_module in \
  01_runtime_config.sh \
  02_container_scheduler.sh \
  03_species_helpers.sh \
  04_busco_runtime.sh \
  05_workspace_io.sh \
  06_workspace_validation.sh \
  07_execution_runtime.sh \
  08_execution_reporting.sh \
  09_sequence_databases.sh \
  10_reference_databases.sh \
  11_fasta.sh; do
  # shellcheck disable=SC1090
  source "${gg_util_dir}/gg_util/${gg_util_module}"
done
unset gg_util_module gg_util_dir
