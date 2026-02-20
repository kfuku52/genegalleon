# Runtime profile for workflow/gg_input_generation_cmd.sh
# Loaded automatically by gg_input_generation_cmd.sh if present.

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

# Input source selection:
# For gg_input_generation_job.sh + GG_INPUT_DATASET_HOST_PATH, use /external/gfe_dataset.
dataset_root="/external/gfe_dataset"
input_dir=""
download_manifest=""
download_dir="/workspace/output/input_download_cache"
manifest_output="/workspace/input/manifest/download_manifest.tsv"
summary_output="/workspace/output/inputprep/inputprep_runs.tsv"
