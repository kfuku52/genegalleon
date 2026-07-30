#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../.." && pwd -P)"
workflow_dir="${repo_root}/workflow"

validate_age_task_range() {
  local task_range=$1
  local task_spec=""
  local task_start=""
  local task_end=""
  local -a task_specs=()

  if [[ ! "${task_range}" =~ ^[1-9][0-9]*(-[1-9][0-9]*(:[1-9][0-9]*)?)?(,[1-9][0-9]*(-[1-9][0-9]*(:[1-9][0-9]*)?)?)*$ ]]; then
    return 1
  fi

  IFS=',' read -r -a task_specs <<< "${task_range}"
  for task_spec in "${task_specs[@]}"; do
    if [[ "${task_spec}" =~ ^([1-9][0-9]*)-([1-9][0-9]*)(:[1-9][0-9]*)?$ ]]; then
      task_start=${BASH_REMATCH[1]}
      task_end=${BASH_REMATCH[2]}
      if (( task_start > task_end )); then
        return 1
      fi
    fi
  done
}

usage() {
  cat <<'EOF'
Usage:
  bash workflow/sites/shirokane/gg_shirokane_submit.sh \
    --entrypoint <gg_*_entrypoint.sh> [options]

Required:
  --entrypoint PATH          GeneGalleon entrypoint under workflow/.

Options:
  --tasks RANGE              AGE array range, for example 1-20 or 1-100:5.
  --slots N                  Override the entrypoint def_slot request.
  --mem-per-slot NG          Override the entrypoint s_vmem request.
  --resource NAME            ljob (default), mjob, lmem, or exclusive.
  --job-name NAME            Override the AGE job name.
  --verify                   Validate with qsub -verify without submitting.
  --dry-run                  Print the qsub command without running it.
  -h, --help                 Show this help.

Examples:
  bash workflow/sites/shirokane/gg_shirokane_submit.sh \
    --entrypoint gg_gene_evolution_entrypoint.sh --tasks 1-20

  bash workflow/sites/shirokane/gg_shirokane_submit.sh \
    --entrypoint gg_transcriptome_generation_entrypoint.sh \
    --tasks 1-8 --slots 4 --mem-per-slot 32G --verify
EOF
}

entrypoint=""
tasks=""
slots=""
mem_per_slot=""
resource="ljob"
job_name=""
verify=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --entrypoint)
      [[ $# -ge 2 ]] || { echo "--entrypoint requires a value." >&2; exit 2; }
      entrypoint=$2
      shift 2
      ;;
    --tasks)
      [[ $# -ge 2 ]] || { echo "--tasks requires a value." >&2; exit 2; }
      tasks=$2
      shift 2
      ;;
    --slots)
      [[ $# -ge 2 ]] || { echo "--slots requires a value." >&2; exit 2; }
      slots=$2
      shift 2
      ;;
    --mem-per-slot)
      [[ $# -ge 2 ]] || { echo "--mem-per-slot requires a value." >&2; exit 2; }
      mem_per_slot=$2
      shift 2
      ;;
    --resource)
      [[ $# -ge 2 ]] || { echo "--resource requires a value." >&2; exit 2; }
      resource=$2
      shift 2
      ;;
    --job-name)
      [[ $# -ge 2 ]] || { echo "--job-name requires a value." >&2; exit 2; }
      job_name=$2
      shift 2
      ;;
    --verify)
      verify=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "${entrypoint}" ]]; then
  echo "--entrypoint is required." >&2
  usage >&2
  exit 2
fi

case "${entrypoint}" in
  */*)
    if [[ "${entrypoint}" == /* ]]; then
      entrypoint_path="${entrypoint}"
    else
      entrypoint_path="${PWD}/${entrypoint}"
    fi
    ;;
  *)
    entrypoint_path="${workflow_dir}/${entrypoint}"
    ;;
esac
if [[ ! -f "${entrypoint_path}" ]]; then
  echo "Entrypoint was not found: ${entrypoint_path}" >&2
  exit 2
fi
entrypoint_path="$(cd "$(dirname "${entrypoint_path}")" && pwd -P)/$(basename "${entrypoint_path}")"
case "${entrypoint_path}" in
  "${workflow_dir}"/gg_*_entrypoint.sh)
    ;;
  *)
    echo "Entrypoint must be a workflow/gg_*_entrypoint.sh file: ${entrypoint_path}" >&2
    exit 2
    ;;
esac

if [[ -n "${tasks}" ]] && ! validate_age_task_range "${tasks}"; then
  echo "Invalid AGE task range: ${tasks}" >&2
  exit 2
fi
if [[ -n "${slots}" ]] && [[ ! "${slots}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--slots must be a positive integer: ${slots}" >&2
  exit 2
fi
if [[ -n "${mem_per_slot}" ]] && [[ ! "${mem_per_slot}" =~ ^[1-9][0-9]*G$ ]]; then
  echo "--mem-per-slot must be a positive whole number of GiB such as 8G: ${mem_per_slot}" >&2
  exit 2
fi
case "${resource}" in
  ljob|mjob|lmem|exclusive)
    ;;
  *)
    echo "--resource must be one of: ljob, mjob, lmem, exclusive" >&2
    exit 2
    ;;
esac
if [[ -n "${job_name}" ]] && [[ ! "${job_name}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
  echo "Invalid AGE job name: ${job_name}" >&2
  exit 2
fi

qsub_command=( qsub -terse )
if [[ "${verify}" -eq 1 ]]; then
  qsub_command+=( -verify )
fi
qsub_command+=( -S /bin/bash -cwd )
if [[ -n "${slots}" ]]; then
  qsub_command+=( -pe def_slot "${slots}" )
fi
resource_request="${resource}"
if [[ -n "${mem_per_slot}" ]]; then
  resource_request="s_vmem=${mem_per_slot},${resource_request}"
fi
qsub_command+=( -l "${resource_request}" )
if [[ -n "${tasks}" ]]; then
  qsub_command+=( -t "${tasks}" )
fi
if [[ -n "${job_name}" ]]; then
  qsub_command+=( -N "${job_name}" )
fi
qsub_command+=( "${entrypoint_path}" )

printf 'SHIROKANE qsub command:'
printf ' %q' "${qsub_command[@]}"
printf '\n'

if [[ "${dry_run}" -eq 1 ]]; then
  exit 0
fi
if ! command -v qsub >/dev/null 2>&1; then
  echo "qsub was not found. Run this helper on a SHIROKANE login node." >&2
  exit 1
fi

cd "${repo_root}"
exec "${qsub_command[@]}"
