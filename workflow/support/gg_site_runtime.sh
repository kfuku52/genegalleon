#!/usr/bin/env bash
set -euo pipefail

gg_detect_site_profile() {
  local hostname_text

  if [[ -n "${GG_SITE_PROFILE:-}" ]]; then
    echo "${GG_SITE_PROFILE}"
    return 0
  fi

  hostname_text="$(hostname 2>/dev/null || true)"
  case "${hostname_text}" in
    at*|m*|igt*|it*)
      echo "nig"
      return 0
      ;;
    *.nhr.fau.de)
      echo "nhr-fau"
      return 0
      ;;
    *)
      echo "default"
      return 0
      ;;
  esac
}

gg_site_scheduler_prelude() {
  case "$(gg_detect_site_profile)" in
    nig)
      if [[ -n "${PBS_O_WORKDIR:-}" ]]; then
        cd "${PBS_O_WORKDIR}" || return 1
        if [[ -d "/bio/package/singularity/singularity_3.0/bin" ]]; then
          export PATH="${PATH}:/bio/package/singularity/singularity_3.0/bin"
        fi
      fi
      ;;
    *)
      if [[ -n "${PBS_O_WORKDIR:-}" ]]; then
        cd "${PBS_O_WORKDIR}" || return 1
      fi
      ;;
  esac
}

gg_site_container_shell_command() {
  local runtime_bin=$1
  local echo_header="set_singularity_command: "
  local site_profile

  site_profile="$(gg_detect_site_profile)"
  case "${site_profile}" in
    nig)
      echo "${echo_header}site profile = nig" >&2
      if [[ -e /var/spool/uge ]]; then
        gg_add_container_bind_mount "/var/spool/uge:/var/spool/uge"
      fi
      if [[ -e /var/spool/age ]]; then
        gg_add_container_bind_mount "/var/spool/age:/var/spool/age"
      fi
      if [[ -e /opt/pkg ]]; then
        gg_add_container_bind_mount "/opt/pkg:/opt/pkg"
      fi
      if [[ -e /home/geadmin/UGER/uger/spool ]]; then
        gg_add_container_bind_mount "/home/geadmin/UGER/uger/spool:/home/geadmin/UGER/uger/spool"
      fi
      echo "${runtime_bin} shell"
      return 0
      ;;
    nhr-fau)
      echo "${echo_header}site profile = nhr-fau" >&2
      echo "${runtime_bin} shell --contain"
      return 0
      ;;
    *)
      echo "${echo_header}site profile = default" >&2
      echo "${runtime_bin} shell"
      return 0
      ;;
  esac
}
