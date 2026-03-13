#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
download_source_url="${1:-https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip}"

if [[ -z "${download_source_url}" ]]; then
  echo "Notung download source URL is empty" >&2
  exit 1
fi

if [[ "${download_source_url}" == *.zip || "${download_source_url}" == *.zip\?* ]]; then
  printf '%s\n' "${download_source_url}"
  exit 0
fi

base_url="${download_source_url%/*}/"

html="$(bash "${script_dir}/download_url_with_fallback.sh" "${download_source_url}")"

latest_zip="$(
  printf '%s\n' "${html}" \
    | tr '\r' '\n' \
    | grep -Eo 'Notung-2\.9(\.[0-9]+)+\.zip' \
    | sort -u \
    | sed 's/^Notung-//' \
    | sort -V \
    | tail -n 1 \
    | sed 's/^/Notung-/'
)"

if [[ -z "${latest_zip}" ]]; then
  echo "failed to resolve stable Notung 2.9 zip from ${download_source_url}" >&2
  exit 1
fi

printf '%s%s\n' "${base_url}" "${latest_zip}"
