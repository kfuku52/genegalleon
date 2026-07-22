# shellcheck shell=bash
# Genetic-code and FASTA collection helpers.
# This file is sourced by workflow/support/gg_util.sh.

get_hyphy_genetic_code() {
  local ncbi_genetic_code=$1
  if [[ ${ncbi_genetic_code} -eq 1 ]]; then
    echo "Universal"
  elif [[ ${ncbi_genetic_code} -eq 2 ]]; then
    echo "Vertebrate-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 3 ]]; then
    echo "Yeast-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 4 ]]; then
    echo "Mold-Protozoan-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 5 ]]; then
    echo "Invertebrate-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 6 ]]; then
    echo "Ciliate-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 9 ]]; then
    echo "Echinoderm-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 10 ]]; then
    echo "Euplotid-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 12 ]]; then
    echo "Alt-Yeast-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 13 ]]; then
    echo "Ascidian-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 14 ]]; then
    echo "Flatworm-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 15 ]]; then
    echo "Blepharisma-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 16 ]]; then
    echo "Chlorophycean-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 21 ]]; then
    echo "Trematode-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 22 ]]; then
    echo "Scenedesmus-obliquus-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 23 ]]; then
    echo "Thraustochytrium-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 24 ]]; then
    echo "Pterobranchia-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 25 ]]; then
    echo "SR1-and-Gracilibacteria"
  elif [[ ${ncbi_genetic_code} -eq 26 ]]; then
    echo "Pachysolen-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 29 ]]; then
    echo "Mesodinium-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 30 ]]; then
    echo "Peritrich-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 33 ]]; then
    echo "Cephalodiscidae-mtDNA"
  else
    echo "Unknown"
    >&2 echo "This NCBI genetic code cannot be used in HyPhy: ${ncbi_genetic_code}"
  fi
}

get_fasta_extensions_for_grep() {
  printf -- "-e %s " ".fa$" ".fa.gz$" ".fas$" ".fas.gz$" ".fasta$" ".fasta.gz$" ".fna$" ".fna.gz$"
}

gg_find_fasta_files() {
  local search_dir=$1
  local maxdepth=${2:-1}
  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi
  find -H "${search_dir}" -maxdepth "${maxdepth}" -type f ! -name '.*' \
  \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) \
  | sort
}

gg_find_file_basenames() {
  local search_dir=$1
  local name_pattern=${2:-*}
  local maxdepth=${3:-1}
  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi
  find -H "${search_dir}" -maxdepth "${maxdepth}" -type f ! -name '.*' -name "${name_pattern}" \
  | awk -F'/' '{print $NF}' \
  | sort
}

is_species_set_identical() {
  local return_flag=0
  local dir1=$1
  local dir2=$2
  if [[ ! -d "${dir1}" || ! -d "${dir2}" ]]; then
    echo "Directory not found for species-set comparison: ${dir1} / ${dir2}"
    return 1
  fi
  local files1=()
  local files2=()
  local f
  shopt -s nullglob dotglob
  for f in "${dir1}"/*; do
    [[ -f "${f}" ]] || continue
    [[ "$(basename "${f}")" == .* ]] && continue
    files1+=( "${f}" )
  done
  for f in "${dir2}"/*; do
    [[ -f "${f}" ]] || continue
    [[ "$(basename "${f}")" == .* ]] && continue
    files2+=( "${f}" )
  done
  shopt -u nullglob dotglob
  local sp1=()
  local sp2=()
  for f in "${files1[@]}"; do
    sp1+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  for f in "${files2[@]}"; do
    sp2+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  local sp1_unique=()
  local sp2_unique=()
  local sp_name
  while IFS= read -r sp_name; do
    sp1_unique+=( "${sp_name}" )
  done < <(printf '%s\n' "${sp1[@]}" | sort -u)
  while IFS= read -r sp_name; do
    sp2_unique+=( "${sp_name}" )
  done < <(printf '%s\n' "${sp2[@]}" | sort -u)
  sp1=( "${sp1_unique[@]}" )
  sp2=( "${sp2_unique[@]}" )
  local num_sp1=${#sp1[@]}
  local num_sp2=${#sp2[@]}
  if [[ ${num_sp1} -ne ${num_sp2} ]]; then
    echo "Number of unique species in ${dir1} and ${dir2} are different."
    return_flag=1
  fi
  if [[ "${sp1[*]}" != "${sp2[*]}" ]]; then
    echo "Species names are different between ${dir1} and ${dir2}"
    echo "dir1 (${dir1}): ${sp1[*]}"
    echo "dir2 (${dir2}): ${sp2[*]}"
    return_flag=1
  fi
  return "${return_flag}"
}
