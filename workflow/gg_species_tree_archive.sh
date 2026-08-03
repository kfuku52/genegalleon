#!/usr/bin/env bash
set -euo pipefail

gg_archive_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
exec python3 "${gg_archive_script_dir}/support/species_tree_output_store.py" "$@"
