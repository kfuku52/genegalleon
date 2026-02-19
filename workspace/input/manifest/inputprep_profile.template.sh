# Example profile for workflow/gg_inputPrep_cmd.sh
# Copy this file to:
#   workspace/input/manifest/inputprep_profile.sh
# and adjust values. The job script auto-loads that path.

# 0/1 switches
run_build_manifest=1
run_format_inputs=1
run_validate_inputs=1

# all|ensemblplants|phycocosm|phytozome
provider="all"

# strict checks and overwrite policy
strict=0
overwrite=0
download_only=0
dry_run=0

# download request options
download_timeout=120
auth_bearer_token_env="" # e.g., GFE_DOWNLOAD_BEARER_TOKEN
http_header="" # e.g., "User-Agent: genegalleon-inputprep"

# input source selection
# If using host bind via GG_INPUT_DATASET_HOST_PATH in gg_inputPrep_job.sh,
# this can be "/external/gfe_dataset".
dataset_root=""
input_dir=""
download_manifest=""
download_dir=""
manifest_output=""
summary_output=""
