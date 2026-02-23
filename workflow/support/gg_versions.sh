#!/usr/bin/env bash
set -euo pipefail

dir_pg="/workspace"
dir_myscript="/script/support"
dir_db="/usr/local/db"
dir_rpsblastdb="${dir_db}/Pfam_LE"
dir_jaspardb="${dir_db}/jaspar"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "" 1 1
if [[ -d /opt/conda/bin ]]; then
  export PATH="/opt/conda/bin:${PATH}"
fi

gg_print_section "genegalleon version"
echo "${GG_VERSION:-unknown}"
gg_print_spacer

if command -v conda >/dev/null 2>&1; then
  mapfile -t conda_envs < <(
    conda env list 2>/dev/null \
    | awk 'NR>2 && $1 !~ /^#/ {print $1}' \
    | sed '/^$/d'
  )
  for conda_env in "${conda_envs[@]}"; do
    gg_print_section "Installed program versions in the conda environment: ${conda_env}"
    if conda activate "${conda_env}" >/dev/null 2>&1; then
      conda list || true
      conda deactivate >/dev/null 2>&1 || true
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

gg_print_section "In-house scripts in: ${dir_myscript}"
ls -la "${dir_myscript}"
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
