#!/usr/bin/env bash

dir_pg="/workspace"
dir_myscript="/script/script"
dir_db="/usr/local/db"
dir_rpsblastdb="${dir_db}/Pfam_LE"
dir_jaspardb="${dir_db}/jaspar"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "" 1 1

mapfile -t conda_envs < <(conda env list | awk '/^[A-Za-z0-9]/ {print $1}')
for conda_env in "${conda_envs[@]}"; do
  gg_print_section "Installed program versions in the conda environment: ${conda_env}"
  conda activate "${conda_env}"
  conda list
  conda deactivate
  gg_print_spacer
done

gg_print_section "Versions of programs installed by apt"
apt list --installed
gg_print_spacer

gg_print_section "Versions of programs installed without conda"

gg_print_labelled_path_listing "Pfam database for RPS-BLAST search" "${dir_rpsblastdb}"

gg_print_labelled_path_listing "UNIPROT database for DIAMOND annotation" "${dir_db}/uniprot_sprot.pep"

gg_print_labelled_path_listing "JASPAR database for cis element annotation" "${dir_jaspardb}"

echo 'Map NH used in mapdNdS'
mapnh
gg_print_spacer

echo 'MaxAlign'
maxalign.pl -V
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
