# shellcheck shell=bash
# Sourced by workflow/support/gg_util.sh.

_download_busco_dataset_mapping_files_locked() {
  local mapping_dir=$1
  local stamp_file=$2
  local odb_version=${3:-}
  local py_exec=""

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot prepare BUSCO placement mappings." >&2
    return 1
  fi

  if ! GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
    GG_BUSCO_MAPPING_STAMP="${stamp_file}" \
    GG_BUSCO_MAPPING_ODB_VERSION="${odb_version}" \
    "${py_exec}" - <<'PY'
import io
import os
import re
import tarfile
import urllib.request
from pathlib import Path

base_url = "https://busco-data.ezlab.org/v5/data/placement_files/"
mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
stamp_file = Path(os.environ.get("GG_BUSCO_MAPPING_STAMP", "").strip())
odb_version = os.environ.get("GG_BUSCO_MAPPING_ODB_VERSION", "").strip()
if not mapping_dir:
    raise SystemExit("GG_BUSCO_MAPPING_DIR is empty.")
if not stamp_file:
    raise SystemExit("GG_BUSCO_MAPPING_STAMP is empty.")
if not odb_version:
    raise SystemExit("GG_BUSCO_MAPPING_ODB_VERSION is empty.")

mapping_dir.mkdir(parents=True, exist_ok=True)
html = urllib.request.urlopen(base_url, timeout=120).read().decode("utf-8", "replace")
selected = []
for domain in ("archaea", "bacteria", "eukaryota"):
    pattern = re.compile(
        rf'mapping_taxids-busco_dataset_name\.{domain}_odb{re.escape(odb_version)}\.[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}}\.txt\.tar\.gz'
    )
    matches = sorted(set(pattern.findall(html)))
    if not matches:
        raise SystemExit(f"BUSCO placement mapping for domain/version not found: {domain}, odb{odb_version}")
    archive_name = matches[-1]
    archive_url = base_url + archive_name
    archive_bytes = urllib.request.urlopen(archive_url, timeout=120).read()
    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:gz") as tar:
        members = [member for member in tar.getmembers() if member.isfile() and member.name.endswith(".txt")]
        if not members:
            raise SystemExit(f"No mapping text file found inside archive: {archive_name}")
        member = members[0]
        extracted = tar.extractfile(member)
        if extracted is None:
            raise SystemExit(f"Failed to extract mapping text from archive: {archive_name}")
        out_path = mapping_dir / Path(member.name).name
        out_path.write_bytes(extracted.read())
        selected.append((domain, archive_name, out_path.name))

stamp_lines = [f"odb_version\t{odb_version}", "domain\tarchive\tmapping_file"]
stamp_lines.extend("\t".join(item) for item in selected)
stamp_file.write_text("\n".join(stamp_lines) + "\n", encoding="utf-8")
PY
  then
    echo "Failed to prepare BUSCO placement mappings in: ${mapping_dir}" >&2
    return 1
  fi

  if [[ ! -s "${stamp_file}" ]]; then
    echo "BUSCO placement mapping stamp file is missing: ${stamp_file}" >&2
    return 1
  fi
}

ensure_busco_dataset_mapping_files() {
  local gg_workspace_dir=$1
  local dir_db
  local mapping_dir
  local lock_file
  local stamp_file
  local remote_odb_version=""
  local odb_version=""

  dir_db="$(workspace_downloads_root "${gg_workspace_dir}")/busco_downloads"
  mapping_dir=$(workspace_busco_placement_root "${gg_workspace_dir}")
  lock_file="${dir_db}/locks/busco_dataset_mappings.lock"

  ensure_dir "${dir_db}"
  ensure_dir "${dir_db}/locks"
  ensure_dir "${mapping_dir}"

  if remote_odb_version=$(gg_fetch_latest_busco_mapping_odb_version 2>/dev/null); then
    odb_version="${remote_odb_version}"
  else
    odb_version=$(gg_latest_busco_mapping_odb_version_from_dir "${mapping_dir}" || true)
  fi
  if [[ -z "${odb_version}" ]]; then
    echo "Failed to determine a BUSCO placement mapping ODB version." >&2
    return 1
  fi
  stamp_file="${mapping_dir}/mapping_taxids-busco_dataset_name.odb${odb_version}.ready.tsv"

  if [[ -n "${remote_odb_version}" ]]; then
    gg_array_download_once "${lock_file}" "${stamp_file}" "BUSCO taxid-to-dataset mapping files (odb${odb_version})" \
      _download_busco_dataset_mapping_files_locked "${mapping_dir}" "${stamp_file}" "${odb_version}" || return 1
  elif [[ ! -s "${stamp_file}" ]]; then
    printf 'odb_version\t%s\n' "${odb_version}" > "${stamp_file}"
  fi

  if [[ ! -s "${stamp_file}" ]]; then
    echo "Failed to prepare BUSCO placement mapping stamp file: ${stamp_file}" >&2
    return 1
  fi

  echo "${mapping_dir}"
}

gg_resolve_busco_lineage_from_lineages() {
  local requested=${1:-auto}
  local mapping_dir=${2:-}
  shift 2 || true
  local normalized_requested=""
  local requested_lc=""
  local py_exec=""

  normalized_requested=$(gg_normalize_busco_lineage_request "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi
  if [[ -z "${mapping_dir}" ]]; then
    echo "gg_resolve_busco_lineage_from_lineages: mapping_dir is empty." >&2
    return 1
  fi
  if [[ $# -eq 0 ]]; then
    echo "gg_resolve_busco_lineage_from_lineages: no lineage taxid strings were provided." >&2
    return 1
  fi

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot resolve BUSCO lineage from taxid lineages." >&2
    return 1
  fi

  GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
  "${py_exec}" - "$@" <<'PY'
import os
import re
import sys
from pathlib import Path

mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
if not mapping_dir:
    raise SystemExit("GG_BUSCO_MAPPING_DIR is empty.")

mapping = {}
pattern = re.compile(r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\..*\.txt$")
mapping_files_by_version = {}
for mapping_file in mapping_dir.glob("mapping_taxids-busco_dataset_name.*_odb*.txt"):
    match = pattern.fullmatch(mapping_file.name)
    if not match:
        continue
    _, version = match.groups()
    mapping_files_by_version.setdefault(int(version), []).append(mapping_file)
eligible_versions = sorted(mapping_files_by_version.keys(), reverse=True)
if not eligible_versions:
    raise SystemExit(f"No BUSCO mapping files were found in {mapping_dir}")

lineages = []
for raw_lineage in sys.argv[1:]:
    tokens = [token.strip() for token in str(raw_lineage).split(",") if token.strip()]
    if not tokens:
        raise SystemExit("Encountered an empty lineage taxid string.")
    lineages.append([int(token) for token in tokens])

common_taxids = set(lineages[0])
for lineage in lineages[1:]:
    common_taxids &= set(lineage)

for selected_version in eligible_versions:
    mapping = {}
    for mapping_file in sorted(mapping_files_by_version[selected_version]):
        with mapping_file.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) != 2:
                    raise SystemExit(f"Malformed BUSCO mapping line in {mapping_file}: {raw_line.rstrip()}")
                taxid_str, dataset = parts
                mapping[int(taxid_str)] = dataset.strip()
    for taxid in reversed(lineages[0]):
        if taxid in common_taxids and taxid in mapping:
            print(mapping[taxid])
            raise SystemExit(0)

raise SystemExit(
    "No BUSCO dataset mapping matched the deepest common taxid across the provided species lineages in any available ODB version."
)
PY
}

gg_busco_taxid_lineages_for_species() {
  local gg_workspace_dir=$1
  shift || true
  local py_exec=""
  local db_file=""
  local species_name=""
  local normalized_species=()

  if [[ $# -eq 0 ]]; then
    echo "gg_busco_taxid_lineages_for_species: no species names were provided." >&2
    return 1
  fi

  if ! ensure_ete_taxonomy_db "${gg_workspace_dir}" >&2; then
    echo "Failed to prepare ETE taxonomy DB for BUSCO lineage resolution." >&2
    return 1
  fi
  db_file=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot resolve species taxid lineages." >&2
    return 1
  fi

  for species_name in "$@"; do
    species_name=$(gg_normalize_annotation_species "${species_name}")
    [[ -n "${species_name}" ]] || continue
    normalized_species+=( "${species_name}" )
  done
  if [[ ${#normalized_species[@]} -eq 0 ]]; then
    echo "gg_busco_taxid_lineages_for_species: all provided species names were empty after normalization." >&2
    return 1
  fi

  GG_TAXONOMY_DBFILE="${db_file}" \
  "${py_exec}" - "${normalized_species[@]}" <<'PY'
import os
import sys

db_file = os.environ.get("GG_TAXONOMY_DBFILE", "").strip()
if not db_file:
    raise SystemExit("GG_TAXONOMY_DBFILE is empty.")

from ete4 import NCBITaxa  # pragma: no cover - runtime dependency

ncbi = NCBITaxa(dbfile=db_file)

def resolve_species_taxid(species_name: str) -> int:
    candidates = []
    normalized = species_name.replace("_", " ").strip()
    tokens = normalized.split()
    if normalized:
        candidates.append(normalized)
    if len(tokens) >= 3 and tokens[2].lower() in {"cf", "aff", "nr"}:
        candidates.extend(
            [
                f"{tokens[0]} {tokens[2]}. {tokens[1]}",
                f"{tokens[0]} {tokens[2]} {tokens[1]}",
                " ".join(tokens[:2]),
            ]
        )
    elif len(tokens) >= 3 and tokens[1].lower() in {"cf", "aff", "nr"}:
        candidates.extend(
            [
                f"{tokens[0]} {tokens[1]}. {tokens[2]}",
                f"{tokens[0]} {tokens[1]} {tokens[2]}",
                f"{tokens[0]} {tokens[2]}",
            ]
        )
    elif len(tokens) >= 4 and tokens[2].lower() in {"subsp", "ssp", "var", "forma", "f"}:
        rank = tokens[2].lower()
        if rank in {"subsp", "ssp"}:
            candidates.extend(
                [
                    f"{tokens[0]} {tokens[1]} subsp. {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} subsp {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} ssp. {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} ssp {tokens[3]}",
                ]
            )
        elif rank == "var":
            candidates.extend(
                [
                    f"{tokens[0]} {tokens[1]} var. {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} var {tokens[3]}",
                ]
            )
        else:
            candidates.extend(
                [
                    f"{tokens[0]} {tokens[1]} forma {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} forma. {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} f. {tokens[3]}",
                    f"{tokens[0]} {tokens[1]} f {tokens[3]}",
                ]
            )
        candidates.append(" ".join(tokens[:2]))
    elif len(tokens) >= 3 and tokens[1].lower() == "sp":
        candidates.extend(
            [
                f"{tokens[0]} sp. {tokens[2]}",
                f"{tokens[0]} sp {tokens[2]}",
                tokens[0],
            ]
        )
    elif len(tokens) >= 2:
        genus_species = " ".join(tokens[:2])
        if genus_species not in candidates:
            candidates.append(genus_species)
    for candidate in candidates:
        translated = ncbi.get_name_translator([candidate])
        if translated:
            return next(iter(translated.values()))[0]
    raise SystemExit(f"Species name was not found in the NCBI taxonomy database: {species_name}")

for species_name in sys.argv[1:]:
    taxid = resolve_species_taxid(species_name)
    lineage = ncbi.get_lineage(taxid)
    if not lineage:
        raise SystemExit(f"NCBI taxonomy lineage was empty for species: {species_name}")
    print(",".join(str(taxid_item) for taxid_item in lineage))
PY
}

gg_resolve_busco_lineage() {
  local gg_workspace_dir=$1
  local requested=${2:-auto}
  shift 2 || true
  local normalized_requested=""
  local requested_lc=""
  local mapping_dir=""
  local lineage_output=""
  local lineage_item=""
  local lineages=()

  normalized_requested=$(gg_normalize_busco_lineage_request "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi
  if [[ $# -eq 0 ]]; then
    echo "gg_resolve_busco_lineage: no species names were provided for auto resolution." >&2
    return 1
  fi

  if ! mapping_dir=$(ensure_busco_dataset_mapping_files "${gg_workspace_dir}"); then
    echo "Failed to prepare BUSCO placement mapping files for auto lineage resolution." >&2
    return 1
  fi
  if ! lineage_output=$(gg_busco_taxid_lineages_for_species "${gg_workspace_dir}" "$@"); then
    echo "Failed to resolve NCBI taxid lineages for BUSCO auto resolution." >&2
    return 1
  fi
  while IFS= read -r lineage_item; do
    [[ -n "${lineage_item}" ]] || continue
    lineages+=( "${lineage_item}" )
  done <<< "${lineage_output}"
  if [[ ${#lineages[@]} -eq 0 ]]; then
    echo "No lineage taxid strings were produced for BUSCO auto resolution." >&2
    return 1
  fi

  local auto_lineage=""
  if ! auto_lineage=$(gg_resolve_busco_lineage_from_lineages "auto" "${mapping_dir}" "${lineages[@]}"); then
    return 1
  fi
  gg_finalize_auto_busco_lineage_name "${auto_lineage}"
}

gg_fasta_relabel_headers_to_species() {
  local line=""
  local header=""
  local species_name=""
  while IFS= read -r line; do
    if [[ "${line}" == ">"* ]]; then
      header=${line#>}
      header=${header##*/}
      header=${header%%[[:space:]]*}
      species_name=$(gg_species_name_from_path_or_dot "${header}")
      if [[ -z "${species_name}" ]]; then
        species_name="${header}"
      fi
      printf '>%s\n' "${species_name}"
      continue
    fi
    printf '%s\n' "${line}"
  done
}

gg_busco_gene_tokens() {
  local mode=$1
  shift
  local token
  local split_token

  case "${mode}" in
    exclude_duplicated)
      for token in "$@"; do
        if [[ "${token}" == "-" || "${token}" == *","* ]]; then
          continue
        fi
        echo "${token}"
      done
      ;;
    split_duplicated)
      for token in "$@"; do
        if [[ "${token}" == "-" ]]; then
          continue
        fi
        for split_token in ${token//,/ }; do
          if [[ -n "${split_token}" ]]; then
            echo "${split_token}"
          fi
        done
      done
      ;;
    *)
      echo "gg_busco_gene_tokens: invalid mode: ${mode}" >&2
      return 1
      ;;
  esac
}

gg_seqkit_grep_by_patterns_from_infile_list() {
  local threads=$1
  local infile_list=$2
  shift 2
  local patterns=( "$@" )
  seqkit grep --threads "${threads}" --pattern-file <(printf '%s\n' "${patterns[@]}") --infile-list "${infile_list}"
}

gg_prepare_cds_fasta_stream() {
  local threads=${1:-1}
  local codontable=${2:-}
  if [[ -n "${codontable}" ]]; then
    seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" \
    | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${threads}" \
    | cdskit pad --codon_table "${codontable}"
  else
    seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" \
    | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${threads}" \
    | cdskit pad
  fi
}

gg_prepare_cdskit_localize_cds_input() {
  local infile=$1
  local outfile=$2
  local threads=${3:-1}
  local codontable=${4:-1}

  seqkit seq --remove-gaps --upper-case --only-id --threads "${threads}" "${infile}" |
    seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" |
    cdskit pad --codon_table "${codontable}" |
    cdskit mask --codon_table "${codontable}" --stop_codon no --ambiguous_codon yes --mask_char 'N' \
      > "${outfile}"
}

gg_run_cdskit_localize() {
  local seqfile=$1
  local seqtype=$2
  local report=$3
  local model=$4
  local organism_group=$5
  local include_features=${6:-0}
  local no_model_download=${7:-0}
  local threads=${8:-1}
  local codontable=${9:-1}
  local include_features_arg="no"
  local model_download_arg="yes"

  if [[ "${include_features}" == "1" || "${include_features}" == "yes" || "${include_features}" == "true" ]]; then
    include_features_arg="yes"
  fi
  if [[ "${no_model_download}" == "1" || "${no_model_download}" == "yes" || "${no_model_download}" == "true" ]]; then
    model_download_arg="no"
  fi

  export CDSKIT_MODEL_DIR="${CDSKIT_MODEL_DIR:-${gg_workspace_downloads_dir}/cdskit_models}"
  ensure_dir "${CDSKIT_MODEL_DIR}"

  cdskit localize \
    --seq_file "${seqfile}" \
    --in_seq_format fasta \
    --codon_table "${codontable}" \
    --threads "${threads}" \
    --model "${model}" \
    --model_download "${model_download_arg}" \
    --report "${report}" \
    --include_features "${include_features_arg}" \
    --seq_type "${seqtype}" \
    --organism_group "${organism_group}"
}

_download_busco_lineage_to_runtime() {
  local busco_lineage=$1
  local runtime_busco_db=$2
  local runtime_busco_lineage=$3
  local runtime_ready_marker=$4
  if gg_busco_lineage_is_ready "${runtime_busco_lineage}"; then
    gg_write_ready_marker "${runtime_ready_marker}"
    return 0
  fi
  echo "Starting BUSCO dataset download: ${busco_lineage}" >&2
  if ! busco --download "${busco_lineage}" >&2; then
    echo "BUSCO dataset download failed: ${busco_lineage}" >&2
    return 1
  fi
  if [[ -d busco_downloads ]]; then
    gg_merge_directory_contents "busco_downloads" "${runtime_busco_db}" || return 1
    rm -rf -- busco_downloads
  fi
  if ! gg_busco_lineage_is_ready "${runtime_busco_lineage}"; then
    echo "BUSCO lineage dataset is still missing after download: ${runtime_busco_lineage}" >&2
    return 1
  fi
  gg_write_ready_marker "${runtime_ready_marker}"
  echo "BUSCO dataset download has been finished: ${busco_lineage}" >&2
}

gg_merge_directory_contents() {
  local staged_dir=$1
  local runtime_dir=$2
  local staged_entries=()
  local staged_entry=""
  local runtime_entry=""

  [[ -d "${staged_dir}" ]] || return 0
  ensure_dir "${runtime_dir}"

  shopt -s nullglob dotglob
  staged_entries=("${staged_dir}"/*)
  shopt -u nullglob dotglob

  if [[ ${#staged_entries[@]} -eq 0 ]]; then
    return 0
  fi

  for staged_entry in "${staged_entries[@]}"; do
    runtime_entry="${runtime_dir}/$(basename "${staged_entry}")"
    if [[ -d "${staged_entry}" && -d "${runtime_entry}" ]]; then
      gg_merge_directory_contents "${staged_entry}" "${runtime_entry}" || return 1
      rm -rf -- "${staged_entry}"
    else
      if [[ -e "${runtime_entry}" ]]; then
        rm -rf -- "${runtime_entry}" || return 1
      fi
      mv -- "${staged_entry}" "${runtime_entry}" || return 1
    fi
  done
}

gg_busco_lineage_is_ready() {
  local lineage_dir=$1
  [[ -d "${lineage_dir}" ]] || return 1
  [[ -s "${lineage_dir}/dataset.cfg" || -s "${lineage_dir}/info/dataset.cfg" ]]
}

ensure_busco_download_path() {
  local gg_workspace_dir=$1
  local busco_lineage=$2
  local sys_busco_db="/usr/local/db/busco_downloads"
  local sys_busco_lineage="${sys_busco_db}/lineages/${busco_lineage}"
  local runtime_busco_db
  local runtime_busco_lineage
  local runtime_ready_marker
  local lock_file

  if [[ -z "${busco_lineage}" ]]; then
    echo "ensure_busco_download_path: busco_lineage is empty." >&2
    return 1
  fi

  if gg_busco_lineage_is_ready "${sys_busco_lineage}"; then
    echo "${sys_busco_db}"
    return 0
  fi

  runtime_busco_db="$(workspace_downloads_root "${gg_workspace_dir}")/busco_downloads"
  runtime_busco_lineage="${runtime_busco_db}/lineages/${busco_lineage}"
  runtime_ready_marker="${runtime_busco_lineage}/.download.ready"
  lock_file="${runtime_busco_db}/lineages/.busco_${busco_lineage}.download.lock"
  ensure_dir "${runtime_busco_db}"
  ensure_dir "${runtime_busco_db}/lineages"

  if gg_busco_lineage_is_ready "${runtime_busco_lineage}"; then
    gg_write_ready_marker "${runtime_ready_marker}"
    echo "${runtime_busco_db}"
    return 0
  fi

  gg_array_download_once "${lock_file}" "${runtime_ready_marker}" "BUSCO dataset download (${busco_lineage})" \
    _download_busco_lineage_to_runtime "${busco_lineage}" "${runtime_busco_db}" "${runtime_busco_lineage}" "${runtime_ready_marker}" || return 1

  if ! gg_busco_lineage_is_ready "${runtime_busco_lineage}"; then
    echo "Failed to prepare BUSCO lineage dataset: ${busco_lineage}" >&2
    return 1
  fi

  echo "${runtime_busco_db}"
}
