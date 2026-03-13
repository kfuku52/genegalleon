#!/usr/bin/env bash
set -euo pipefail

path="${1:-}"
expected_sha256="${2:-}"

if [[ -z "${path}" || -z "${expected_sha256}" ]]; then
  echo "usage: $0 <path> <expected-sha256>" >&2
  exit 1
fi

if [[ ! -f "${path}" ]]; then
  echo "file not found: ${path}" >&2
  exit 1
fi

compute_sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${path}" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${path}" | awk '{print $1}'
  else
    echo "verify_sha256.sh requires sha256sum or shasum" >&2
    exit 1
  fi
}

actual_sha256="$(compute_sha256 | tr '[:upper:]' '[:lower:]')"
expected_sha256="$(printf '%s\n' "${expected_sha256}" | tr '[:upper:]' '[:lower:]')"

if [[ "${actual_sha256}" != "${expected_sha256}" ]]; then
  echo "sha256 mismatch for ${path}" >&2
  echo "expected: ${expected_sha256}" >&2
  echo "actual:   ${actual_sha256}" >&2
  exit 1
fi
