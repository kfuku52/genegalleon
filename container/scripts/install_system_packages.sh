#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
profile="${1:-}"
runtime_packages=()
development_packages=()
while IFS= read -r package; do
  [[ -z "${package}" || "${package}" == \#* ]] && continue
  runtime_packages+=("${package}")
done < "${script_dir}/../apt/runtime.txt"
while IFS= read -r package; do
  [[ -z "${package}" || "${package}" == \#* ]] && continue
  development_packages+=("${package}")
done < "${script_dir}/../apt/development.txt"

case "${profile}" in
  runtime|development)
    packages=("${runtime_packages[@]}")
    if [[ "${profile}" == "development" ]]; then
      packages+=("${development_packages[@]}")
    fi
    apt-get update
    apt-get upgrade -y
    apt-get install -y --no-install-recommends "${packages[@]}"
    sed -i 's/^# *en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen
    locale-gen en_US.UTF-8
    ;;
  prune-development)
    # Native SIF builds have no retained Docker layers. Keep runtime packages
    # explicit before removing build tools and their unused dependencies.
    apt-mark manual "${runtime_packages[@]}"
    apt-get purge --auto-remove -y "${development_packages[@]}"
    ;;
  *)
    echo "Usage: $0 runtime|development|prune-development" >&2
    exit 2
    ;;
esac
apt-get clean
rm -rf -- /var/lib/apt/lists/* /var/cache/apt/archives/*
