#!/usr/bin/env bash
set -euo pipefail

url="${1:-}"
dest="${2:-}"

if [[ -z "${url}" ]]; then
  echo "usage: $0 <url> [dest]" >&2
  exit 1
fi

amberjack_host="amberjack.compbio.cs.cmu.edu"
amberjack_host_ip="${NOTUNG_DOWNLOAD_HOST_IP:-128.2.205.60}"

tmp_output=""
output_path="${dest}"
if [[ -z "${output_path}" ]]; then
  tmp_output="$(mktemp)"
  output_path="${tmp_output}"
fi

cleanup() {
  if [[ -n "${tmp_output}" ]]; then
    rm -f -- "${tmp_output}"
  fi
}
trap cleanup EXIT

url_no_scheme="${url#*://}"
host_port="${url_no_scheme%%/*}"
host="${host_port%%:*}"

curl_common=(
  -fsSL
  --retry 5
  --retry-all-errors
  --retry-delay 2
  --connect-timeout 30
  --max-time 600
)

run_curl() {
  rm -f -- "${output_path}"
  curl "${curl_common[@]}" "$@" "${url}" -o "${output_path}"
}

if ! run_curl; then
  if [[ "${host}" == "${amberjack_host}" && -n "${amberjack_host_ip}" ]]; then
    echo "[download_url_with_fallback] retrying ${url} via ${amberjack_host_ip}" >&2
    run_curl --resolve "${amberjack_host}:443:${amberjack_host_ip}"
  else
    echo "[download_url_with_fallback] failed to download ${url}" >&2
    exit 1
  fi
fi

if [[ -z "${dest}" ]]; then
  cat "${output_path}"
fi
