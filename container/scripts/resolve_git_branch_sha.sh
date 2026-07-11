#!/usr/bin/env bash
set -euo pipefail

repo_url="${1:-}"
branch="${2:-}"

if [[ -z "${repo_url}" || -z "${branch}" ]]; then
  echo "usage: $0 <repo-url> <branch>" >&2
  exit 1
fi

sha="$(git ls-remote --exit-code "${repo_url}" "refs/heads/${branch}" | awk 'NR == 1 { print $1 }')"
if [[ ! "${sha}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "Could not resolve ${repo_url} branch ${branch} to a commit SHA." >&2
  exit 1
fi

printf '%s\n' "${sha}"
