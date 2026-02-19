#!/usr/bin/env bash

# Shared defaults for commonly reused workflow parameters.
# Individual gg_*_cmd.sh scripts may override these values.

: "${GG_COMMON_GENETIC_CODE:=1}"
: "${GG_COMMON_BUSCO_LINEAGE:=embryophyta_odb12}"
