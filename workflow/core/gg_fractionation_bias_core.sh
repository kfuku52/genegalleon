#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 0 1

run_kffractbias="${run_kffractbias:-1}"
delete_tmp_dir="${delete_tmp_dir:-1}"
kffractbias_pairs_file="${kffractbias_pairs_file:-}"
if [[ -z "${kffractbias_pairs_file}" ]]; then
  kffractbias_pairs_file="${gg_workspace_input_dir}/fractionation_bias_pairs.tsv"
fi

for binary_flag in run_kffractbias delete_tmp_dir; do
  if [[ ! "${!binary_flag}" =~ ^[01]$ ]]; then
    echo "Invalid ${binary_flag}=${!binary_flag}; expected 0 or 1." >&2
    exit 1
  fi
done

if [[ ! -s "${kffractbias_pairs_file}" ]]; then
  if [[ ${run_kffractbias} -eq 0 ]]; then
    gg_step_skip "kfFractBias pairwise comparison (pair table not configured)"
    exit 0
  fi
  echo "kfFractBias pair table not found or empty: ${kffractbias_pairs_file}" >&2
  exit 1
fi
if [[ ! "${GG_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-}; expected a positive integer." >&2
  exit 1
fi

resolve_species_file() {
  local search_dir=$1
  local species_label=$2
  local data_label=$3
  local matches=()

  if [[ ! -d "${search_dir}" ]]; then
    echo "${data_label} directory not found: ${search_dir}" >&2
    return 1
  fi
  mapfile -t matches < <(gg_find_species_files_by_label "${search_dir}" "${species_label}")
  if [[ ${#matches[@]} -eq 0 ]]; then
    echo "No ${data_label} file matched species label '${species_label}' in ${search_dir}." >&2
    return 1
  fi
  if [[ ${#matches[@]} -gt 1 ]]; then
    echo "Multiple ${data_label} files matched species label '${species_label}' in ${search_dir}:" >&2
    printf '  %s\n' "${matches[@]}" >&2
    return 1
  fi
  printf '%s\n' "${matches[0]}"
}

dir_tmp_root="${gg_workspace_output_dir}/tmp/kffractbias"
ensure_dir "${dir_tmp_root}"
pair_values_file="${dir_tmp_root}/.pair.${GG_ARRAY_TASK_ID}.$$.values"
trap 'rm -f -- "${pair_values_file}"' EXIT HUP INT TERM

python - "${kffractbias_pairs_file}" "${GG_ARRAY_TASK_ID}" > "${pair_values_file}" <<'PY'
import csv
import re
import sys
from collections import Counter
from pathlib import Path

path = Path(sys.argv[1])
task_id = int(sys.argv[2])
required = ("analysis_id", "target_species", "query_species", "quota")
defaults = {
    "mode": "compare",
    "window_size": "100",
    "step_size": "1",
    "denominator": "all",
    "target_seqids": "",
    "query_seqids": "",
    "exclude_seqid_regex": "",
    "cscore": "0.7",
    "aligner": "last",
    "target_feature": "",
    "target_attribute": "",
    "query_feature": "",
    "query_attribute": "",
    "minimum_mapping_fraction": "0.5",
    "self_hit_percent": "98",
    "diagonal_bound": "300",
}
fields = required + tuple(defaults)

with path.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"Pair table has no header: {path}")
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise SystemExit("Pair table contains duplicate column names")
    missing = [field for field in required if field not in reader.fieldnames]
    if missing:
        raise SystemExit(f"Pair table is missing required columns: {', '.join(missing)}")
    unknown = [field for field in reader.fieldnames if field not in fields]
    if unknown:
        raise SystemExit(f"Pair table contains unsupported columns: {', '.join(unknown)}")
    rows = []
    for line_number, candidate in enumerate(reader, start=2):
        if None in candidate:
            raise SystemExit(f"Pair table row {line_number} contains more fields than the header")
        if any((value or "").strip() for value in candidate.values()):
            rows.append(candidate)

analysis_ids = [(candidate.get("analysis_id") or "").strip() for candidate in rows]
duplicates = sorted(value for value, count in Counter(analysis_ids).items() if value and count > 1)
if duplicates:
    raise SystemExit(f"analysis_id values must be unique; duplicates: {', '.join(duplicates)}")

if task_id > len(rows):
    raise SystemExit(f"GG_ARRAY_TASK_ID={task_id} is out of range for {len(rows)} pair rows")
row = rows[task_id - 1]
values = {field: (row.get(field) or defaults.get(field, "")).strip() for field in fields}

for field in required:
    if not values[field]:
        raise SystemExit(f"Pair row {task_id} has an empty required value: {field}")
if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", values["analysis_id"]):
    raise SystemExit("analysis_id must be one safe filename component using letters, numbers, dot, underscore, or hyphen")
if not re.fullmatch(r"[1-9][0-9]*:[1-9][0-9]*", values["quota"]):
    raise SystemExit("quota must have the form positive_integer:positive_integer, for example 1:2")
if values["mode"] not in {"compare", "self"}:
    raise SystemExit("mode must be compare or self")
quota_left, quota_right = (int(value) for value in values["quota"].split(":"))
for field in ("window_size", "step_size", "diagonal_bound"):
    try:
        valid = int(values[field]) >= 1
    except ValueError:
        valid = False
    if not valid:
        raise SystemExit(f"{field} must be a positive integer")
if values["denominator"] not in {"all", "syntenic"}:
    raise SystemExit("denominator must be all or syntenic")
if values["aligner"] not in {"last", "blast"}:
    raise SystemExit("aligner must be last or blast")
for field in ("cscore", "minimum_mapping_fraction"):
    try:
        valid = 0 < float(values[field]) <= 1
    except ValueError:
        valid = False
    if not valid:
        raise SystemExit(f"{field} must be greater than 0 and at most 1")
try:
    valid_self_hit_percent = 0 < float(values["self_hit_percent"]) <= 100
except ValueError:
    valid_self_hit_percent = False
if not valid_self_hit_percent:
    raise SystemExit("self_hit_percent must be greater than 0 and at most 100")
for prefix in ("target", "query"):
    if bool(values[f"{prefix}_feature"]) != bool(values[f"{prefix}_attribute"]):
        raise SystemExit(f"{prefix}_feature and {prefix}_attribute must be specified together")
if values["mode"] == "self":
    if values["target_species"] != values["query_species"]:
        raise SystemExit("self mode requires identical target_species and query_species")
    if quota_left != quota_right:
        raise SystemExit("self mode requires a symmetric quota such as 1:1 or 2:2")
    if values["query_seqids"]:
        raise SystemExit("self mode uses target_seqids for the shared sequence filter; query_seqids must be empty")
    if values["query_feature"] or values["query_attribute"]:
        raise SystemExit("self mode uses target_feature and target_attribute; query overrides must be empty")

for field in fields:
    value = values[field]
    if "\n" in value or "\r" in value:
        raise SystemExit(f"Pair row values cannot contain newlines: {field}")
    print(value)
PY

mapfile -t pair_values < "${pair_values_file}"
rm -f -- "${pair_values_file}"
trap - EXIT HUP INT TERM
if [[ ${#pair_values[@]} -ne 20 ]]; then
  echo "Internal error: expected 20 normalized pair values, found ${#pair_values[@]}." >&2
  exit 1
fi

analysis_id=${pair_values[0]}
target_species=${pair_values[1]}
query_species=${pair_values[2]}
quota=${pair_values[3]}
analysis_mode=${pair_values[4]}
window_size=${pair_values[5]}
step_size=${pair_values[6]}
denominator=${pair_values[7]}
target_seqids=${pair_values[8]}
query_seqids=${pair_values[9]}
exclude_seqid_regex=${pair_values[10]}
cscore=${pair_values[11]}
aligner=${pair_values[12]}
target_feature=${pair_values[13]}
target_attribute=${pair_values[14]}
query_feature=${pair_values[15]}
query_attribute=${pair_values[16]}
minimum_mapping_fraction=${pair_values[17]}
self_hit_percent=${pair_values[18]}
diagonal_bound=${pair_values[19]}

target_cds=$(resolve_species_file "${gg_workspace_input_dir}/species_cds" "${target_species}" "target CDS") || exit 1
target_gff=$(resolve_species_file "${gg_workspace_input_dir}/species_gff" "${target_species}" "target GFF") || exit 1
query_cds=""
query_gff=""
if [[ "${analysis_mode}" == "compare" ]]; then
  query_cds=$(resolve_species_file "${gg_workspace_input_dir}/species_cds" "${query_species}" "query CDS") || exit 1
  query_gff=$(resolve_species_file "${gg_workspace_input_dir}/species_gff" "${query_species}" "query GFF") || exit 1
  dir_result="${gg_workspace_output_dir}/kffractbias/${analysis_id}"
  dir_provenance="${gg_workspace_output_dir}/artifact_provenance/kffractbias"
  artifact_stage="fractionation_bias_compare"
else
  dir_result="${gg_workspace_output_dir}/genome_evolution/self_fractionation_bias/${analysis_id}"
  dir_provenance="${gg_workspace_output_dir}/artifact_provenance/genome_evolution/self_fractionation_bias"
  artifact_stage="genome_evolution_self_fractionation_bias"
fi
file_genes="${dir_result}/${analysis_id}.genes.tsv"
file_windows="${dir_result}/${analysis_id}.windows.tsv"
file_summary="${dir_result}/${analysis_id}.summary.json"
file_pdf="${dir_result}/${analysis_id}.plot.pdf"
file_png="${dir_result}/${analysis_id}.plot.png"
file_synteny_archive="${dir_result}/${analysis_id}.synteny.zip"
file_manifest="${dir_provenance}/${analysis_id}.json"

needs_update=0
gg_artifact_contract_init kffractbias_provenance_args \
  "${artifact_stage}" "${analysis_id}" "${file_manifest}"
kffractbias_provenance_args+=(
  --input "target_cds=${target_cds}"
  --input "target_gff=${target_gff}"
  --output "genes=${file_genes}"
  --output "windows=${file_windows}"
  --output "summary=${file_summary}"
  --output "plot_pdf=${file_pdf}"
  --output "plot_png=${file_png}"
  --output "synteny_archive=${file_synteny_archive}"
  --parameter "pair_table=${kffractbias_pairs_file}"
  --parameter "array_task_id=${GG_ARRAY_TASK_ID}"
  --parameter "mode=${analysis_mode}"
  --parameter "target_species=${target_species}"
  --parameter "query_species=${query_species}"
  --parameter "quota=${quota}"
  --parameter "window_size=${window_size}"
  --parameter "step_size=${step_size}"
  --parameter "denominator=${denominator}"
  --parameter "target_seqids=${target_seqids}"
  --parameter "query_seqids=${query_seqids}"
  --parameter "exclude_seqid_regex=${exclude_seqid_regex}"
  --parameter "cscore=${cscore}"
  --parameter "aligner=${aligner}"
  --parameter "target_feature=${target_feature}"
  --parameter "target_attribute=${target_attribute}"
  --parameter "query_feature=${query_feature}"
  --parameter "query_attribute=${query_attribute}"
  --parameter "minimum_mapping_fraction=${minimum_mapping_fraction}"
  --parameter "self_hit_percent=${self_hit_percent}"
  --parameter "diagonal_bound=${diagonal_bound}"
)
if [[ "${analysis_mode}" == "compare" ]]; then
  kffractbias_provenance_args+=(
    --input "query_cds=${query_cds}"
    --input "query_gff=${query_gff}"
  )
fi
gg_artifact_prepare_stage needs_update run_kffractbias "${kffractbias_provenance_args[@]}" || exit $?

if [[ "${analysis_mode}" == "compare" ]]; then
  task="kfFractBias ${target_species} vs ${query_species} (${quota})"
else
  task="kfFractBias self-synteny retention for ${target_species} (${quota})"
fi
if [[ ${needs_update} -eq 1 && ${run_kffractbias} -eq 1 ]]; then
  gg_step_start "${task}"
  dir_task_tmp=$(mktemp -d "${dir_tmp_root}/${GG_ARRAY_TASK_ID}_${analysis_id}.XXXXXX")
  ensure_dir "${dir_task_tmp}/result"

  if [[ "${analysis_mode}" == "compare" ]]; then
    kffractbias_args=(
      compare
      --target-cds "${target_cds}"
      --target-gff "${target_gff}"
      --query-cds "${query_cds}"
      --query-gff "${query_gff}"
      --target-name "${target_species}"
      --query-name "${query_species}"
      --quota "${quota}"
    )
  else
    self_depth=${quota%%:*}
    kffractbias_args=(
      selfcompare
      --cds "${target_cds}"
      --gff "${target_gff}"
      --name "${target_species}"
      --depth "${self_depth}"
      --self-hit-percent "${self_hit_percent}"
      --diagonal-bound "${diagonal_bound}"
    )
  fi
  kffractbias_args+=(
    --window-size "${window_size}"
    --step-size "${step_size}"
    --denominator "${denominator}"
    --cscore "${cscore}"
    --aligner "${aligner}"
    --minimum-mapping-fraction "${minimum_mapping_fraction}"
    --cpus "${GG_TASK_CPUS:-1}"
    --output-dir "${dir_task_tmp}/result"
    --prefix "${analysis_id}"
    --force
  )
  if [[ -n "${target_seqids}" ]]; then
    if [[ "${analysis_mode}" == "compare" ]]; then
      kffractbias_args+=(--target-seqids "${target_seqids}")
    else
      kffractbias_args+=(--seqids "${target_seqids}")
    fi
  fi
  if [[ "${analysis_mode}" == "compare" && -n "${query_seqids}" ]]; then
    kffractbias_args+=(--query-seqids "${query_seqids}")
  fi
  if [[ -n "${exclude_seqid_regex}" ]]; then
    kffractbias_args+=(--exclude-seqid-regex "${exclude_seqid_regex}")
  fi
  if [[ -n "${target_feature}" ]]; then
    if [[ "${analysis_mode}" == "compare" ]]; then
      kffractbias_args+=(--target-feature "${target_feature}" --target-attribute "${target_attribute}")
    else
      kffractbias_args+=(--feature "${target_feature}" --attribute "${target_attribute}")
    fi
  fi
  if [[ "${analysis_mode}" == "compare" && -n "${query_feature}" ]]; then
    kffractbias_args+=(--query-feature "${query_feature}" --query-attribute "${query_attribute}")
  fi

  kffractbias "${kffractbias_args[@]}"

  source_genes="${dir_task_tmp}/result/${analysis_id}.genes.tsv"
  source_windows="${dir_task_tmp}/result/${analysis_id}.windows.tsv"
  source_summary="${dir_task_tmp}/result/${analysis_id}.summary.json"
  source_pdf="${dir_task_tmp}/result/${analysis_id}.plot.pdf"
  source_png="${dir_task_tmp}/result/${analysis_id}.plot.png"
  source_synteny_dir="${dir_task_tmp}/result/${analysis_id}.synteny"
  source_synteny_archive="${dir_task_tmp}/${analysis_id}.synteny.zip"
  for expected_file in "${source_genes}" "${source_windows}" "${source_summary}" "${source_pdf}" "${source_png}"; do
    if [[ ! -s "${expected_file}" ]]; then
      echo "kfFractBias did not create the expected output: ${expected_file}" >&2
      exit 1
    fi
  done
  if [[ ! -d "${source_synteny_dir}" ]]; then
    echo "kfFractBias did not create the expected synteny directory: ${source_synteny_dir}" >&2
    exit 1
  fi
  (
    cd "${dir_task_tmp}/result"
    zip -rq "${source_synteny_archive}" "${analysis_id}.synteny"
  )

  python - "${source_summary}" "${file_genes}" "${file_windows}" "${file_pdf}" "${file_png}" "${file_synteny_archive}" "${analysis_mode}" "${target_species}" "${quota}" <<'PY'
import json
import os
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
genes, windows, pdf, png, archive, analysis_mode, species, quota = sys.argv[2:]
payload = json.loads(summary_path.read_text(encoding="utf-8"))
payload["outputs"] = {
    "genes": genes,
    "windows": windows,
    "plot_pdf": pdf,
    "plot_png": png,
    "synteny_archive": archive,
}
for label in ("synteny", "target_bed", "query_bed"):
    item = payload.get("inputs", {}).get(label)
    if item is not None:
        source_path = Path(item.get("path", ""))
        item["path"] = f"{archive}::{summary_path.stem.removesuffix('.summary')}.synteny/{source_path.name}"
payload["genegalleon"] = {
    "array_task_id": os.environ.get("GG_ARRAY_TASK_ID", ""),
    "analysis_mode": analysis_mode,
    "species": species,
    "quota": quota,
}
temporary = summary_path.with_suffix(summary_path.suffix + ".tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
temporary.replace(summary_path)
PY

  mv_out_bundle \
    "${source_genes}" "${file_genes}" \
    "${source_windows}" "${file_windows}" \
    "${source_summary}" "${file_summary}" \
    "${source_pdf}" "${file_pdf}" \
    "${source_png}" "${file_png}" \
    "${source_synteny_archive}" "${file_synteny_archive}"
  gg_artifact_record "${kffractbias_provenance_args[@]}"

  if [[ ${delete_tmp_dir} -eq 1 ]]; then
    rm -rf -- "${dir_task_tmp}"
  fi
else
  gg_step_skip "${task}"
fi

echo "kfFractBias outputs: ${dir_result}"
