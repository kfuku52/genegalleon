# shellcheck shell=bash
# Species-label, annotation, and BUSCO helpers.
# This file is sourced by workflow/support/gg_util.sh.

gg_tmp_root() {
	local candidate
	for candidate in "${TMPDIR:-}" "/tmp" "."; do
		if [[ -n "${candidate}" && -d "${candidate}" && -w "${candidate}" ]]; then
			printf '%s\n' "${candidate}"
			return 0
		fi
	done
	printf '.\n'
}

# Optional template arguments are used by helpers outside the ShellCheck driver set.
# shellcheck disable=SC2120
gg_mktemp() {
	local make_dir=0
	local template=${1:-}
	local tmp_root=""

	if [[ "${template}" == "-d" ]]; then
		make_dir=1
		template=${2:-}
	fi

	if [[ -n "${template}" ]]; then
		tmp_root=$(dirname "${template}")
		if [[ ! -d "${tmp_root}" || ! -w "${tmp_root}" ]]; then
			tmp_root=$(gg_tmp_root)
			template="${tmp_root%/}/$(basename "${template}")"
		fi
	else
		tmp_root=$(gg_tmp_root)
		template="${tmp_root%/}/gg.tmp.XXXXXX"
	fi

	if [[ ${make_dir} -eq 1 ]]; then
		command mktemp -d "${template}"
	else
		command mktemp "${template}"
	fi
}

ensure_dir() {
	local d=$1
	if [[ -z "${d}" ]]; then
		return 0
	fi
	if [[ -d "${d}" ]]; then
		return 0
	fi

	local attempt
	local mkdir_err
	mkdir_err=$(gg_mktemp)
	for attempt in 1 2 3 4 5; do
		if mkdir -p "${d}" 2>"${mkdir_err}"; then
			rm -f -- "${mkdir_err}"
			return 0
		fi
		if grep -qi "Resource deadlock avoided" "${mkdir_err}"; then
			sleep 1
			continue
		fi
		cat "${mkdir_err}" >&2
		rm -f -- "${mkdir_err}"
		return 1
	done

	cat "${mkdir_err}" >&2
	rm -f -- "${mkdir_err}"
	return 1
}

ensure_parent_dir() {
	local p=$1
	if [[ -z "${p}" ]]; then
		return 0
	fi
	ensure_dir "$(dirname -- "${p}")"
}

gg_species_name_from_path() {
  local path=$1
  gg_species_name_from_path_or_dot "${path}"
}

_gg_species_rank_token_key() {
  local token=${1:-}
  local multiplication_sign=$'\303\227'
  token=${token%.}
  token=$(printf '%s' "${token}" | tr '[:upper:]' '[:lower:]')
  case "${token}" in
    x|hybrid|"${multiplication_sign}")
      printf '%s\n' "x"
      ;;
    spp)
      printf '%s\n' "sp"
      ;;
    ssp|subspecies)
      printf '%s\n' "subsp"
      ;;
    variety)
      printf '%s\n' "var"
      ;;
    form)
      printf '%s\n' "forma"
      ;;
    *)
      printf '%s\n' "${token}"
      ;;
  esac
}

_gg_species_prefix_token_count() {
  local -a parts=("$@")
  local second=""
  local third=""

  if [[ ${#parts[@]} -lt 2 ]]; then
    printf '0\n'
    return 0
  fi

  second=$(_gg_species_rank_token_key "${parts[1]}")
  if [[ ${#parts[@]} -ge 3 ]]; then
    third=$(_gg_species_rank_token_key "${parts[2]}")
  fi

  if [[ "${second}" == "sp" ]]; then
    if [[ ${#parts[@]} -ge 3 ]]; then
      printf '3\n'
    else
      printf '2\n'
    fi
    return 0
  fi
  if [[ "${second}" == "x" ]]; then
    if [[ ${#parts[@]} -ge 3 ]]; then
      printf '3\n'
    else
      printf '2\n'
    fi
    return 0
  fi
  if [[ "${second}" == "cf" || "${second}" == "aff" || "${second}" == "nr" ]]; then
    if [[ ${#parts[@]} -ge 3 ]]; then
      printf '3\n'
    else
      printf '2\n'
    fi
    return 0
  fi
  if [[ "${third}" == "x" && ${#parts[@]} -ge 5 && "${parts[3]}" =~ ^[[:upper:]] ]]; then
    printf '5\n'
    return 0
  fi
  if [[ "${third}" == "cf" || "${third}" == "aff" || "${third}" == "nr" ]]; then
    printf '3\n'
    return 0
  fi
  case "${third}" in
    subsp|var|forma|f|strain|substrain|serovar|serotype|serogroup|pathovar|pv|biovar|biotype|chemovar|morphovar|cultivar|cv|isolate|group|subgroup|complex|clade|lineage|section|series|ecotype|breed)
    if [[ ${#parts[@]} -ge 4 ]]; then
      printf '4\n'
    else
      printf '3\n'
    fi
    return 0
      ;;
  esac
  printf '2\n'
}

_gg_species_is_rank_or_qualifier_token() {
  local key
  key=$(_gg_species_rank_token_key "${1:-}")
  case "${key}" in
    x|sp|cf|aff|nr|subsp|var|forma|f|strain|substrain|serovar|serotype|serogroup|pathovar|pv|biovar|biotype|chemovar|morphovar|cultivar|cv|isolate|group|subgroup|complex|clade|lineage|section|series|ecotype|breed)
      return 0
      ;;
  esac
  return 1
}

_gg_species_label_prefix_part() {
  local part=${1:-}
  local key
  if _gg_species_is_rank_or_qualifier_token "${part}"; then
    key=$(_gg_species_rank_token_key "${part}")
    if [[ "${key}" == "x" ]]; then
      printf '%s\n' "x"
      return 0
    fi
    printf '%s\n' "${part}"
    return 0
  fi
  printf '%s\n' "${part%%.*}"
}

_gg_strip_species_terminal_suffixes() {
  local stem=${1:-}
  local stem_lc
  local suffix
  local -a suffixes=(
    ".fasta.busco.full.tsv" ".fa.busco.full.tsv" ".faa.busco.full.tsv" ".fna.busco.full.tsv" ".ffn.busco.full.tsv"
    ".fasta.busco.short.txt" ".fa.busco.short.txt" ".faa.busco.short.txt" ".fna.busco.short.txt" ".ffn.busco.short.txt"
    ".busco.full.tsv" "_busco.full.tsv" ".busco.short.txt" "_busco.short.txt"
    ".busco.full" "_busco.full" ".busco.short" "_busco.short" ".busco" "_busco"
    ".derived.cds.fa.gz" "_derived.cds.fa.gz" ".derived.genome.fa.gz" "_derived.genome.fa.gz" ".derived.gff.gz" "_derived.gff.gz"
    ".cds.all.fa.gz" "_cds.all.fa.gz" ".cds.fa.gz" "_cds.fa.gz" ".cds.fna.gz" "_cds.fna.gz"
    ".genome.fa.gz" "_genome.fa.gz" ".genomic.fna.gz" "_genomic.fna.gz" ".dna.primary_assembly.fa.gz" ".dna.toplevel.fa.gz"
    ".pep.fa.gz" "_pep.fa.gz" ".protein.fa.gz" "_protein.fa.gz"
    ".fastq.gz" ".fq.gz" ".fasta.gz" ".fa.gz" ".faa.gz" ".fna.gz" ".ffn.gz"
    ".gff3.gz" ".gff.gz" ".gtf.gz" ".tsv.gz" ".txt.gz" ".csv.gz"
    ".fastq" ".fq" ".fasta" ".fa" ".faa" ".fna" ".ffn"
    ".gff3" ".gff" ".gtf" ".tsv" ".txt" ".csv"
  )

  stem_lc=$(printf '%s' "${stem}" | tr '[:upper:]' '[:lower:]')
  for suffix in "${suffixes[@]}"; do
    if [[ "${stem_lc}" == *"${suffix}" ]]; then
      stem="${stem:0:${#stem}-${#suffix}}"
      break
    fi
  done
  printf '%s\n' "${stem}"
}

gg_species_name_from_path_or_dot() {
  local path=$1
  local basename_path
  local stem
  local prefix_count
  local -a parts=()
  local -a selected=()
  local species_name=""
  local part=""

  basename_path=$(basename "${path}")
  stem=$(_gg_strip_species_terminal_suffixes "${basename_path}")
  IFS='_' read -r -a parts <<< "${stem}"
  prefix_count=$(_gg_species_prefix_token_count "${parts[@]}")
  if [[ "${prefix_count}" -eq 0 ]]; then
    printf '%s\n' ""
    return 0
  fi

  selected=("${parts[@]:0:${prefix_count}}")
  for part in "${selected[@]}"; do
    [[ -n "${part}" ]] || continue
    part=$(_gg_species_label_prefix_part "${part}")
    [[ -n "${part}" ]] || continue
    if [[ -n "${species_name}" ]]; then
      species_name+="_"
    fi
    species_name+="${part}"
  done
  printf '%s\n' "${species_name}"
}

gg_species_file_matches_label() {
  local file_path=$1
  local species_name=$2
  local parsed_species=""

  parsed_species=$(gg_species_name_from_path_or_dot "$(basename "${file_path}")")
  [[ -n "${parsed_species}" && "${parsed_species}" == "${species_name}" ]]
}

gg_find_species_files_by_label() {
  local search_dir=$1
  local species_name=$2
  local file

  [[ -d "${search_dir}" ]] || return 0
  while IFS= read -r file; do
    if gg_species_file_matches_label "${file}" "${species_name}"; then
      printf '%s\n' "${file}"
    fi
  done < <(find "${search_dir}" -maxdepth 1 -type f ! -name '.*' | sort)
}

gg_orthogroup_file_matches_id() {
  local file_name=$1
  local orthogroup_id=$2

  [[ "${file_name}" == "${orthogroup_id}" \
    || "${file_name}" == "${orthogroup_id}".* \
    || "${file_name}" == "${orthogroup_id}"_* \
    || "${file_name}" == "${orthogroup_id}"-* ]]
}

gg_annotation_species_priority() {
  cat <<'EOF'
Arabidopsis_thaliana
Oryza_sativa
Homo_sapiens
Mus_musculus
Danio_rerio
Drosophila_melanogaster
Caenorhabditis_elegans
Saccharomyces_cerevisiae
Schizosaccharomyces_pombe
Escherichia_coli
EOF
}

gg_normalize_annotation_species() {
  local species=${1:-}
  species=$(printf '%s' "${species}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e 's/[[:space:]]\+/_/g' -e 's/_$//')
  printf '%s\n' "${species}"
}

gg_normalize_busco_lineage_request() {
  local requested=${1:-}
  requested=$(printf '%s' "${requested}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  printf '%s\n' "${requested}"
}

gg_finalize_auto_busco_lineage_name() {
  local lineage=${1:-}
  local odb_version=${2:-12}
  if [[ -z "${lineage}" ]]; then
    return 1
  fi
  if [[ "${lineage}" == *_odb* ]]; then
    printf '%s\n' "${lineage}"
    return 0
  fi
  printf '%s_odb%s\n' "${lineage}" "${odb_version}"
}

gg_cdskit_localize_default_model() {
  printf '%s\n' "targeting5-perox-deeploc21-et-v1"
}

gg_normalize_cdskit_localize_organism_group() {
  local group=${1:-}
  group=$(printf '%s' "${group}" | tr '[:upper:]' '[:lower:]' | sed -e 's/[[:space:]-]\+/_/g')
  case "${group}" in
    ""|auto)
      printf '%s\n' "auto"
      ;;
    unknown)
      printf '%s\n' "unknown"
      ;;
    plant|plants|viridiplantae)
      printf '%s\n' "plant"
      ;;
    nonplant|non_plant|non_plants|other|metazoa|fungi|animal|animals)
      printf '%s\n' "non_plant"
      ;;
    *)
      echo "Invalid cdskit_localize_organism_group: ${group}" >&2
      echo "cdskit_localize_organism_group must be auto, unknown, plant, or non_plant." >&2
      return 1
      ;;
  esac
}

gg_cdskit_localize_organism_group_from_busco_lineage() {
  local lineage=${1:-}
  local lineage_lc=""
  lineage_lc=$(printf '%s' "${lineage}" | tr '[:upper:]' '[:lower:]')
  lineage_lc=${lineage_lc%%_odb*}

  case "${lineage_lc}" in
    *viridiplantae*|*embryophyta*|*streptophyta*|*tracheophyta*|*chlorophyta*|*rhodophyta*|*glaucophyta*|*brassicales*|*fabales*|*poales*|*solanales*|*asterales*)
      printf '%s\n' "plant"
      ;;
    *metazoa*|*vertebrata*|*mammalia*|*aves*|*actinopterygii*|*arthropoda*|*insecta*|*diptera*|*hymenoptera*|*nematoda*|*mollusca*|*annelida*|*cnidaria*|*fungi*|*ascomycota*|*basidiomycota*|*bacteria*|*archaea*)
      printf '%s\n' "non_plant"
      ;;
    *)
      printf '%s\n' "unknown"
      ;;
  esac
}

gg_resolve_cdskit_localize_organism_group() {
  local requested_group=${1:-auto}
  local gg_workspace_dir=$2
  local busco_lineage_request=${3:-}
  shift 3 || true
  local normalized_group=""
  local resolved_lineage=""
  local inferred_group=""

  normalized_group=$(gg_normalize_cdskit_localize_organism_group "${requested_group}") || return 1
  if [[ "${normalized_group}" != "auto" ]]; then
    printf '%s\n' "${normalized_group}"
    return 0
  fi

  if [[ -n "${busco_lineage_request}" ]]; then
    if resolved_lineage=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage_request}" "$@" 2>/dev/null); then
      inferred_group=$(gg_cdskit_localize_organism_group_from_busco_lineage "${resolved_lineage}")
      echo "cdskit localize organism_group=${inferred_group} inferred from BUSCO lineage: ${resolved_lineage}" >&2
      printf '%s\n' "${inferred_group}"
      return 0
    fi
    inferred_group=$(gg_cdskit_localize_organism_group_from_busco_lineage "${busco_lineage_request}")
    if [[ "${inferred_group}" != "unknown" ]]; then
      echo "cdskit localize organism_group=${inferred_group} inferred from BUSCO lineage request: ${busco_lineage_request}" >&2
      printf '%s\n' "${inferred_group}"
      return 0
    fi
  fi

  echo "cdskit localize organism_group=unknown; BUSCO lineage did not resolve to plant or non_plant." >&2
  printf '%s\n' "unknown"
}

gg_species_names_from_fasta_dir() {
  local search_dir=${1:-}
  local file
  local species_name=""

  while IFS= read -r file; do
    species_name=$(gg_species_name_from_path_or_dot "${file}")
    [[ -n "${species_name}" ]] || continue
    printf '%s\n' "${species_name}"
  done < <(gg_find_fasta_files "${search_dir}" 1) | sort -u
}

gg_species_names_from_annotation_dir() {
  local search_dir=${1:-}
  local file_base
  local species_name=""

  while IFS= read -r file_base; do
    species_name=$(gg_species_name_from_path_or_dot "${file_base}")
    [[ -n "${species_name}" ]] || continue
    printf '%s\n' "${species_name}"
  done < <(gg_find_file_basenames "${search_dir}" "*_annotation.tsv" 1) | sort -u
}

gg_resolve_annotation_species() {
  local requested=${1:-auto}
  shift || true
  local normalized_requested=""
  local requested_lc=""
  local candidates=()
  local candidate=""
  local available=()
  local preferred=""

  normalized_requested=$(gg_normalize_annotation_species "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi

  for candidate in "$@"; do
    candidate=$(gg_normalize_annotation_species "${candidate}")
    [[ -n "${candidate}" ]] || continue
    candidates+=( "${candidate}" )
  done
  if [[ ${#candidates[@]} -eq 0 ]]; then
    return 1
  fi

  while IFS= read -r candidate; do
    [[ -n "${candidate}" ]] || continue
    available+=( "${candidate}" )
  done < <(printf '%s\n' "${candidates[@]}" | sed -e '/^[[:space:]]*$/d' | sort -u)
  while IFS= read -r preferred; do
    [[ -n "${preferred}" ]] || continue
    for candidate in "${available[@]}"; do
      if [[ "${candidate}" == "${preferred}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
      fi
    done
  done < <(gg_annotation_species_priority)

  printf '%s\n' "${available[0]}"
  return 0
}

gg_find_python_exec() {
  local candidate

  for candidate in python python3 /opt/conda/bin/python /usr/bin/python3; do
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
    if command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    fi
  done
  return 1
}

workspace_busco_placement_root() {
  local gg_workspace_dir=$1
  local dir_db
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "${dir_db}/busco_downloads/placement_files"
}

gg_latest_busco_mapping_odb_version_from_dir() {
  local mapping_dir=${1:-}
  local py_exec=""

  if [[ -z "${mapping_dir}" || ! -d "${mapping_dir}" ]]; then
    return 1
  fi
  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    return 1
  fi

  GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
  "${py_exec}" - <<'PY'
import os
import re
import sys
from pathlib import Path

mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
if not mapping_dir:
    raise SystemExit(1)

pattern = re.compile(r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\..*\.txt$")
versions = set()
for mapping_file in mapping_dir.glob("mapping_taxids-busco_dataset_name.*_odb*.txt"):
    match = pattern.fullmatch(mapping_file.name)
    if not match:
        continue
    _, version = match.groups()
    versions.add(int(version))

if not versions:
    raise SystemExit(1)
print(max(versions))
PY
}

gg_fetch_latest_busco_mapping_odb_version() {
  local py_exec=""

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    return 1
  fi

  "${py_exec}" - <<'PY'
import re
import urllib.request

base_url = "https://busco-data.ezlab.org/v5/data/placement_files/"
html = urllib.request.urlopen(base_url, timeout=120).read().decode("utf-8", "replace")
pattern = re.compile(
    r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\.\d{4}-\d{2}-\d{2}\.txt\.tar\.gz"
)
required_domains = {"archaea", "bacteria", "eukaryota"}
versions = {}
for domain, version in pattern.findall(html):
    versions.setdefault(int(version), set()).add(domain)

eligible = [version for version, domains in versions.items() if required_domains.issubset(domains)]
if not eligible:
    raise SystemExit("No common BUSCO ODB placement mapping version was found across archaea, bacteria, and eukaryota.")
print(max(eligible))
PY
}

gg_resolve_omark_database_url() {
  local py_exec=""
  local current_release_url="${GG_OMARK_CURRENT_RELEASE_URL:-https://omabrowser.org/oma/current/}"
  local fallback_url="https://omabrowser.org/All/LUCA.h5"

  if [[ -n "${GG_OMARK_DB_URL:-}" ]]; then
    printf '%s\n' "${GG_OMARK_DB_URL}"
    return 0
  fi

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    printf '%s\n' "${fallback_url}"
    return 0
  fi

  GG_OMARK_CURRENT_RELEASE_URL="${current_release_url}" \
  GG_OMARK_FALLBACK_URL="${fallback_url}" \
  "${py_exec}" - <<'PY'
import re
import os
import urllib.parse
import urllib.request

current_release_url = os.environ.get("GG_OMARK_CURRENT_RELEASE_URL", "").strip()
fallback_url = os.environ.get("GG_OMARK_FALLBACK_URL", "").strip()
if not current_release_url:
    print(fallback_url or "https://omabrowser.org/All/LUCA.h5")
    raise SystemExit(0)
if not fallback_url:
    fallback_url = "https://omabrowser.org/All/LUCA.h5"

request = urllib.request.Request(
    current_release_url,
    headers={"User-Agent": "Mozilla/5.0 (compatible; GeneGalleon OMArk downloader)"},
)
try:
    html = urllib.request.urlopen(request, timeout=120).read().decode("utf-8", "replace")
except Exception:
    print(fallback_url)
    raise SystemExit(0)

matches = re.findall(r'href="([^"]*LUCA\.h5)"', html, flags=re.IGNORECASE)
if not matches:
    print(fallback_url)
    raise SystemExit(0)

print(urllib.parse.urljoin(current_release_url, matches[-1]))
PY
}
