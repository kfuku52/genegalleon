#!/usr/bin/env bash
set -euo pipefail

url="${1:-}"
dest="${2:-}"

if [[ -z "${url}" ]]; then
  echo "usage: $0 <url> [dest]" >&2
  exit 1
fi

max_attempts=${DOWNLOAD_URL_MAX_ATTEMPTS:-8}
retry_delay_sec=${DOWNLOAD_URL_RETRY_DELAY_SEC:-5}
connect_timeout_sec=${DOWNLOAD_URL_CONNECT_TIMEOUT_SEC:-30}
max_time_sec=${DOWNLOAD_URL_MAX_TIME_SEC:-600}

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

download_with_curl() {
  rm -f -- "${output_path}"
  curl \
    -fL \
    --retry "${max_attempts}" \
    --retry-all-errors \
    --retry-delay "${retry_delay_sec}" \
    --connect-timeout "${connect_timeout_sec}" \
    --max-time "${max_time_sec}" \
    "${url}" \
    -o "${output_path}"
}

download_with_wget() {
  rm -f -- "${output_path}"
  wget \
    --tries="${max_attempts}" \
    --waitretry="${retry_delay_sec}" \
    --timeout="${connect_timeout_sec}" \
    -O "${output_path}" \
    "${url}"
}

if command -v curl >/dev/null 2>&1; then
  if download_with_curl; then
    :
  elif command -v wget >/dev/null 2>&1; then
    download_with_wget
  else
    exit 1
  fi
elif command -v wget >/dev/null 2>&1; then
  download_with_wget
else
  echo "download_url.sh requires curl or wget" >&2
  exit 1
fi

if [[ -z "${dest}" ]]; then
  cat "${output_path}"
fi
