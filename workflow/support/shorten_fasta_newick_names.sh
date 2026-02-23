#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <input.fasta> <output.fasta> <input.newick> <output.newick> <max_length>" >&2
  exit 1
fi

input_fasta=$1
output_fasta=$2
input_newick=$3
output_newick=$4
max_length=$5

if [[ ! -f "${input_fasta}" ]]; then
  echo "Input FASTA was not found: ${input_fasta}" >&2
  exit 1
fi
if [[ ! -f "${input_newick}" ]]; then
  echo "Input Newick was not found: ${input_newick}" >&2
  exit 1
fi
if ! [[ "${max_length}" =~ ^[0-9]+$ ]] || [[ "${max_length}" -lt 2 ]]; then
  echo "max_length must be an integer >= 2: ${max_length}" >&2
  exit 1
fi

tmp_dir=$(mktemp -d)
long_names_file="${tmp_dir}/long_names.txt"
name_pairs_file="${tmp_dir}/name_pairs.tsv"
trap 'rm -rf -- "${tmp_dir}"' EXIT

awk -v max_length="${max_length}" '
  /^>/ {
    seq_name = substr($0, 2)
    if (length(seq_name) > max_length) {
      print seq_name
      print "Detected long sequence name: " seq_name > "/dev/stderr"
    }
  }
' "${input_fasta}" > "${long_names_file}"

if [[ ! -s "${long_names_file}" ]]; then
  cp -- "${input_fasta}" "${output_fasta}"
  cp -- "${input_newick}" "${output_newick}"
  exit 0
fi

python3 - "${input_fasta}" "${long_names_file}" "${name_pairs_file}" "${max_length}" <<'PY'
import pathlib
import sys

input_fasta, long_names_file, name_pairs_file, max_length_raw = sys.argv[1:]
max_length = int(max_length_raw)


def load_fasta_headers(path):
    headers = []
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                headers.append(raw_line.rstrip("\n")[1:])
    return headers


all_headers = load_fasta_headers(input_fasta)
reserved_names = set(all_headers)

long_names = []
seen_long = set()
with open(long_names_file, encoding="utf-8") as handle:
    for raw_line in handle:
        old_name = raw_line.rstrip("\n")
        if not old_name or old_name in seen_long:
            continue
        long_names.append(old_name)
        seen_long.add(old_name)

pairs = []
for index, old_name in enumerate(long_names, start=1):
    counter = index
    while True:
        suffix = f"_{counter}"
        prefix_len = max(1, max_length - len(suffix))
        candidate = f"{old_name[:prefix_len]}{suffix}"
        # Avoid candidate collisions with all original headers and previously assigned names.
        if candidate not in reserved_names:
            break
        counter += 1
        if counter > 10_000_000:
            raise SystemExit(f"Could not find a non-colliding short name for: {old_name}")

    pairs.append((old_name, candidate))
    reserved_names.add(candidate)

pathlib.Path(name_pairs_file).write_text(
    "".join(f"{old}\t{new}\n" for old, new in pairs),
    encoding="utf-8",
)
PY

replace_names() {
  local input_file=$1
  local output_file=$2
  local mode=$3
  python3 - "${input_file}" "${output_file}" "${name_pairs_file}" "${mode}" <<'PY'
import pathlib
import re
import sys

input_file, output_file, name_pairs_file, mode = sys.argv[1:]
text = pathlib.Path(input_file).read_text(encoding="utf-8")
pairs = []
with open(name_pairs_file, encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.rstrip("\n")
        if not line:
            continue
        old_name, new_name = line.split("\t", 1)
        pairs.append((old_name, new_name))

if mode == "newick":
    for old_name, new_name in pairs:
        text = re.sub(
            rf"(?<=\(|,){re.escape(old_name)}(?=:|\)|,)",
            new_name,
            text,
        )
elif mode == "fasta":
    for old_name, new_name in pairs:
        text = re.sub(
            rf"^>{re.escape(old_name)}(?=$| )",
            f">{new_name}",
            text,
            flags=re.MULTILINE,
        )
else:
    raise SystemExit(f"Unsupported replacement mode: {mode}")

pathlib.Path(output_file).write_text(text, encoding="utf-8")
PY
}

replace_names "${input_fasta}" "${output_fasta}" "fasta"
replace_names "${input_newick}" "${output_newick}" "newick"
