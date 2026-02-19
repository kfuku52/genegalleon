#!/usr/bin/env bash
set -euo pipefail

download_page_url="${1:-https://amberjack.compbio.cs.cmu.edu/Notung/download29.html}"

if [[ -z "${download_page_url}" ]]; then
  echo "download page URL is empty" >&2
  exit 1
fi

base_url="${download_page_url%/*}/"

html="$(curl -fsSL --retry 5 --retry-all-errors --retry-delay 2 "${download_page_url}")"

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
  echo "failed to resolve stable Notung 2.9 zip from ${download_page_url}" >&2
  exit 1
fi

printf '%s%s\n' "${base_url}" "${latest_zip}"
