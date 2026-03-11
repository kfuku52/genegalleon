#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "" 1 1

dir_db="/usr/local/db"
dir_rpsblastdb="${dir_db}/Pfam_LE"
dir_jaspardb="${dir_db}/jaspar"
if [[ -d /opt/conda/bin ]]; then
  export PATH="/opt/conda/bin:${PATH}"
fi

gg_print_section "genegalleon version"
echo "${GG_VERSION:-unknown}"
gg_print_spacer

if gg_initialize_conda_shell && command -v conda >/dev/null 2>&1; then
  mapfile -t conda_envs < <(
    conda env list 2>/dev/null \
    | awk 'NR>2 && $1 !~ /^#/ {print $1}' \
    | sed '/^$/d'
  )
  for conda_env in "${conda_envs[@]}"; do
    gg_print_section "Installed program versions in the conda environment: ${conda_env}"
    if gg_activate_conda_env "${conda_env}" >/dev/null 2>&1; then
      conda list || true
      gg_deactivate_conda_env
    else
      echo "Failed to activate conda environment: ${conda_env}"
    fi
    gg_print_spacer
  done
else
  gg_print_section "Conda environments"
  echo "conda command was not found."
  gg_print_spacer
fi

# conda deactivate can drop /opt/conda/bin from PATH.
if [[ -d /opt/conda/bin ]] && [[ ":${PATH}:" != *":/opt/conda/bin:"* ]]; then
  export PATH="/opt/conda/bin:${PATH}"
fi

gg_print_section "Versions of programs installed by apt"
if command -v apt >/dev/null 2>&1; then
  apt list --installed
else
  echo "apt command was not found."
fi
gg_print_spacer

gg_print_section "Versions of programs installed without conda"

gg_print_labelled_path_listing "Pfam database for RPS-BLAST search" "${dir_rpsblastdb}"

gg_print_labelled_path_listing "UNIPROT database for DIAMOND annotation" "${dir_db}/uniprot_sprot.pep"

gg_print_labelled_path_listing "JASPAR database for cis element annotation" "${dir_jaspardb}"

echo 'Map NH used in mapdNdS'
if command -v mapnh >/dev/null 2>&1; then
  mapnh
else
  echo "mapnh: not found"
fi
gg_print_spacer

echo 'CDSKIT'
if command -v cdskit >/dev/null 2>&1; then
  if cdskit --version >/dev/null 2>&1; then
    cdskit --version
  else
    cdskit -h | sed -n '1,2p'
  fi
  echo 'CDSKIT maxalign'
  cdskit maxalign -h | sed -n '1p'
else
  echo "cdskit: not found"
fi
gg_print_spacer

echo 'seqkit'
if command -v seqkit >/dev/null 2>&1; then
	seqkit version
else
	echo "seqkit: not found"
fi
gg_print_spacer

echo 'OMAmer'
if command -v omamer >/dev/null 2>&1; then
  omamer --version || omamer --help | sed -n '1,2p'
else
  echo "omamer: not found"
fi
gg_print_spacer

echo 'OMArk'
if command -v omark >/dev/null 2>&1; then
  omark --help | sed -n '1,2p'
else
  echo "omark: not found"
fi
gg_print_spacer

gg_print_section "In-house scripts in: ${gg_support_dir}"
ls -la "${gg_support_dir}"
gg_print_spacer

gg_print_section "Container OS info"
if [[ -f /etc/os-release ]]; then
	cat /etc/os-release
else
	echo "/etc/os-release was not found in the container."
	uname -a
fi
gg_print_spacer

echo "$(date): Exiting Singularity environment"
