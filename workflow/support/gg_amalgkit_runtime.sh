#!/usr/bin/env bash
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
  set -euo pipefail
fi
# Check the executable in the same environment that will run the workflow.
declare -A gg_amalgkit_help_by_command=()

require_amalgkit_supported_options() {
  local subcommand=$1
  shift
  local option_name="" help_text=""
  if [[ -n ${gg_amalgkit_help_by_command[${subcommand}]+present} ]]; then
    help_text=${gg_amalgkit_help_by_command[${subcommand}]}
  elif help_text=$(amalgkit "${subcommand}" --help 2>&1); then
    gg_amalgkit_help_by_command[${subcommand}]=${help_text}
  else
    echo "Could not read installed amalgkit ${subcommand} help. Check the runtime executable." >&2
    return 1
  fi
  for option_name in "$@"; do
    if [[ ! ${help_text} =~ (^|[^[:alnum:]_-])"${option_name}"([^[:alnum:]_-]|$) ]]; then
      echo "Installed amalgkit ${subcommand} does not support ${option_name}. Update amalgkit to a build that provides this option. Exiting." >&2
      return 1
    fi
  done
}

gg_amalgkit_getfastq_options() {
  local -n options=$1
  options=(
    --download_lock_dir --ncbi_metadata_max_concurrency --ncbi_download_max_concurrency
    --aws_download_max_concurrency --gcp_download_max_concurrency --ena
    --rrna_filter_jobs --rrna_filter_chunk_spots --rrna_filter_memory_limit
    --read_source_format --allow_read_format_fallback
  )
  if [[ ${amalgkit_allow_read_format_fallback:-auto} == auto ]]; then
    options+=('--allow_read_format_fallback auto|yes|no')
  fi
  if [[ -n ${amalgkit_gsa_download_max_concurrency:-} ]]; then
    options+=(--gsa_download_max_concurrency)
  fi
  if [[ ${amalgkit_use_prefetched_reads:-0} == 1 ]]; then
    options+=(--read_stage --prefetch_dir)
  fi
}

gg_amalgkit_require_quant_options() {
  local options=(--quant_backend --oarfish_seq_tech --verified_fastq_proof)
  if [[ -n ${amalgkit_oarfish_options:-} ]]; then
    options+=(--oarfish_options)
  fi
  require_amalgkit_supported_options quant "${options[@]}"
}

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
  if [[ ${1:-} != quantification || $# != 1 ]]; then
    echo "Usage: bash gg_amalgkit_runtime.sh quantification" >&2
    exit 2
  fi
  amalgkit_gsa_download_max_concurrency=${amalgkit_gsa_download_max_concurrency:-2}
  required_options=()
  gg_amalgkit_getfastq_options required_options
  require_amalgkit_supported_options getfastq "${required_options[@]}"
  gg_amalgkit_require_quant_options
  for program in kallisto oarfish fastp; do
    command -v "${program}" >/dev/null || { echo "Quantification runtime is missing ${program}." >&2; exit 1; }
  done
  printf '%s\n' "GeneGalleon quantification runtime verified"
fi
