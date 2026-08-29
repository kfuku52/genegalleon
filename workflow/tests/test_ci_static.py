import re
from pathlib import Path

import yaml
from shell_static_helpers import GITHUB_WORKFLOWS_DIR

REVIEWED_ACTIONS = {
    "actions/cache/restore",
    "actions/cache/save",
    "actions/checkout",
    "actions/download-artifact",
    "actions/setup-python",
    "actions/upload-artifact",
    "crazy-max/ghaction-github-runtime",
    "docker/build-push-action",
    "docker/login-action",
    "docker/setup-buildx-action",
    "r-lib/actions/setup-r",
    "r-lib/actions/setup-r-dependencies",
}

SINGULARITY_CE_VERSION = "4.5.0"
SINGULARITY_CE_DEB_SHA256 = "85e6f7af5e7aad5b1bf28183ce333998bd37eb8f4769af352c47a5153f3373fb"


def load_workflow(name: str) -> dict:
    path = GITHUB_WORKFLOWS_DIR / name
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def workflow_steps(workflow: dict):
    for job in workflow.get("jobs", {}).values():
        yield from job.get("steps", [])


def named_step(job: dict, name: str) -> dict:
    matches = [step for step in job.get("steps", []) if step.get("name") == name]
    assert len(matches) == 1, f"Expected exactly one CI step named {name!r}"
    return matches[0]


def step_run(job: dict, name: str) -> str:
    return str(named_step(job, name).get("run", ""))


def assert_pinned_singularity_runtime(step: dict) -> None:
    run = str(step.get("run", ""))
    env = step.get("env", {})

    assert env["SINGULARITY_CE_VERSION"] == SINGULARITY_CE_VERSION
    assert env["SINGULARITY_CE_DEB_SHA256"] == SINGULARITY_CE_DEB_SHA256
    assert "releases/download/v${SINGULARITY_CE_VERSION}" in run
    assert 'deb_name="singularity-ce_${SINGULARITY_CE_VERSION}-noble_amd64.deb"' in run
    assert "sha256sum --check --strict" in run
    assert 'apt-get install -y "${deb_path}"' in run
    assert "apt-get install -y singularity-container" not in run
    assert 'installed_version="$(singularity version)"' in run
    assert 'expected_version="${SINGULARITY_CE_VERSION}-noble"' in run
    assert '"${installed_version}" != "${expected_version}"' in run


def test_workflows_pin_all_actions_to_reviewed_full_commit_shas():
    uses = {
        step["uses"]
        for path in sorted(GITHUB_WORKFLOWS_DIR.glob("*.yml"))
        for step in workflow_steps(load_workflow(path.name))
        if "uses" in step
    }

    assert uses
    for value in uses:
        action, separator, revision = value.rpartition("@")
        assert separator == "@", f"Action reference lacks a revision: {value}"
        assert action in REVIEWED_ACTIONS, f"Unreviewed action reference: {value}"
        assert re.fullmatch(r"[0-9a-f]{40}", revision), f"Action is not pinned to a full SHA: {value}"
    for action in REVIEWED_ACTIONS:
        assert any(value.startswith(f"{action}@") for value in uses)


def test_workflows_pin_x64_jobs_to_ubuntu_2404():
    workflow_text = "\n".join(path.read_text(encoding="utf-8") for path in sorted(GITHUB_WORKFLOWS_DIR.glob("*.yml")))

    assert "ubuntu-latest" not in workflow_text
    assert "runs-on: ubuntu-24.04" in workflow_text


def test_workflows_limit_github_token_permissions_by_default():
    tests = load_workflow("tests.yml")
    container = load_workflow("container-ghcr.yml")
    release = load_workflow("release-sif.yml")

    assert tests["permissions"] == {"contents": "read"}
    assert container["permissions"] == {"contents": "read"}
    assert release["permissions"] == {"contents": "read"}
    assert container["jobs"]["check-trigger"]["permissions"] == {
        "contents": "read",
        "packages": "read",
    }
    assert container["jobs"]["build-and-push"]["permissions"] == {
        "contents": "read",
        "packages": "write",
    }
    assert container["jobs"]["publish-manifest"]["permissions"] == {
        "contents": "read",
        "packages": "write",
    }
    assert release["jobs"]["build-platform"]["permissions"] == {
        "contents": "read",
        "packages": "write",
    }
    assert release["jobs"]["publish-manifest-and-sif"]["permissions"] == {
        "contents": "write",
        "packages": "write",
    }


def test_tests_workflow_avoids_duplicate_branch_push_and_pull_request_runs():
    triggers = load_workflow("tests.yml")["on"]

    assert triggers["push"]["branches"] == ["main"]
    assert "pull_request" in triggers
    assert "workflow_dispatch" in triggers


def test_workflow_run_blocks_do_not_interpolate_untrusted_event_data():
    unsafe_expression = re.compile(r"\$\{\{\s*(?:github\.event\.|github\.head_ref\b|inputs\.)")
    offenders = []
    for path in sorted(GITHUB_WORKFLOWS_DIR.glob("*.yml")):
        for step in workflow_steps(load_workflow(path.name)):
            run = str(step.get("run", ""))
            if unsafe_expression.search(run):
                offenders.append(f"{path.name}: {step.get('name', '<unnamed>')}")

    assert offenders == []


def test_dependabot_updates_pinned_github_actions():
    dependabot_path = GITHUB_WORKFLOWS_DIR.parent / "dependabot.yml"
    config = yaml.load(dependabot_path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)

    github_actions = [update for update in config["updates"] if update.get("package-ecosystem") == "github-actions"]
    assert len(github_actions) == 1
    assert github_actions[0]["directory"] == "/"
    assert github_actions[0]["schedule"]["interval"] == "weekly"

    pip_updates = [update for update in config["updates"] if update.get("package-ecosystem") == "pip"]
    assert len(pip_updates) == 1
    assert pip_updates[0]["directory"] == "/workflow/tests"
    assert pip_updates[0]["ignore"] == [
        {"dependency-name": "numpy", "versions": [">=2"]},
        {"dependency-name": "matplotlib", "versions": [">=3.11"]},
    ]


def test_immutable_image_tags_fingerprint_resolved_upstream_sources():
    workflow_steps_by_name = {
        "container-ghcr.yml": (
            "prepare-build",
            "Resolve image tags",
            "publish-manifest",
            "Create and push multi-arch manifest",
        ),
        "release-sif.yml": (
            "prepare-release",
            "Resolve release, image tags, and source revisions",
            "publish-manifest-and-sif",
            "Create and push multi-arch release manifest",
        ),
    }
    for workflow_name, (metadata_job, metadata_step, publish_job, publish_step) in workflow_steps_by_name.items():
        jobs = load_workflow(workflow_name)["jobs"]
        metadata_run = step_run(jobs[metadata_job], metadata_step)
        publish_run = step_run(jobs[publish_job], publish_step)

        assert 'repo_name="$(basename "${GITHUB_REPOSITORY}" | tr \'[:upper:]\' \'[:lower:]\')"' in metadata_run
        assert 'source_revision_hash="$(' in metadata_run
        assert "| sha256sum | cut -c1-12" in metadata_run
        assert 'immutable_tag="${date_tag}-${short_sha}-${source_revision_hash}"' in metadata_run
        assert 'imagetools inspect --raw "${IMAGE}:${IMMUTABLE_TAG}"' in publish_run
        assert "imagetools create --dry-run --progress=none" in publish_run
        assert "planned_fingerprint" in publish_run
        assert "existing_fingerprint" in publish_run
        assert "already matches the planned manifest; reusing it" in publish_run
        assert "Refusing to overwrite immutable image tag with different content" in publish_run


def test_sif_runtime_validation_reuses_only_an_exact_runtime_image_or_builds_current():
    jobs = load_workflow("tests.yml")["jobs"]
    sif_job = jobs["sif-runtime-validation"]
    runtime_input = step_run(sif_job, "Resolve exact validation runtime input")
    restore_step = named_step(sif_job, "Restore exact validation SIF")
    selection_step = named_step(sif_job, "Select published or local validation image")
    build_run = step_run(sif_job, "Build current GeneGalleon image")
    build_step = named_step(sif_job, "Build current GeneGalleon image")
    conversion_run = step_run(sif_job, "Build validation SIF from current image")
    identity_run = step_run(sif_job, "Verify exact validation SIF identity")
    save_step = named_step(sif_job, "Save validated SIF by exact runtime input")

    assert "resolve_source_revisions.sh" in runtime_input
    assert "--scope all" in runtime_input
    assert "compute_build_input_hash.sh" in runtime_input
    assert "--runtime" in runtime_input
    assert restore_step["uses"].startswith("actions/cache/restore@")
    assert "steps.runtime-input.outputs.value" in restore_step["with"]["key"]
    assert "io.genegalleon.runtime-input" in selection_step["run"]
    assert "check_published_build_input.py" in selection_step["run"]
    assert "IMAGE_SOURCE=local" in build_run
    assert "IMAGE=local/genegalleon" in build_run
    assert "TAG=ci" in build_run
    assert "FORCE_LOAD_IN_CI=1" in build_run
    assert "bash ./gg_container_build_entrypoint.sh" in build_run
    assert "validation-image.outputs.source == 'local'" in build_step["if"]
    assert build_step["env"]["CACHE_FROM"] == "type=gha,scope=container-linux-amd64"
    assert "type=gha,mode=max,scope=container-linux-amd64" in build_step["env"]["CACHE_TO"]
    assert "event_name != 'pull_request'" in build_step["env"]["CACHE_TO"]
    assert named_step(sif_job, "Expose GitHub runtime for BuildKit cache")["uses"].startswith(
        "crazy-max/ghaction-github-runtime@"
    )
    assert "source_type=docker-daemon" in conversion_run
    assert "source_type=docker" in conversion_run
    assert 'SOURCE="${source_type}"' in conversion_run
    assert "steps.validation-image.outputs.source" in conversion_run
    assert "check_runtime_freshness.sh" in identity_run
    assert "--expected-hash" in identity_run
    assert "steps.runtime-input.outputs.value" in identity_run
    assert save_step["uses"].startswith("actions/cache/save@")
    assert save_step["with"]["key"] == "${{ steps.sif-cache.outputs.cache-primary-key }}"
    assert "github.event_name != 'pull_request'" in save_step["if"]


def test_sif_runtime_validation_starts_after_fast_preflight_guards():
    sif_job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]

    assert set(sif_job["needs"]) == {
        "python-smoke",
        "runtime-change-filter",
        "shell-static",
    }
    assert "runtime-change-filter.outputs.should_run == '1'" in sif_job["if"]


def test_sif_runtime_validation_preserves_disk_headroom_for_conversion():
    sif_job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    runner_cleanup = step_run(sif_job, "Reclaim runner disk space for SIF conversion")
    buildkit_cleanup = step_run(sif_job, "Reclaim BuildKit cache before SIF conversion")
    runtime_install = named_step(sif_job, "Install Singularity runtime")
    conversion = named_step(sif_job, "Build validation SIF from current image")
    disk_report = named_step(sif_job, "Report SIF conversion disk usage")

    assert "docker system prune --all --force --volumes" in runner_cleanup
    for path in ("/opt/ghc", "/usr/local/lib/android", "/usr/share/dotnet"):
        assert path in runner_cleanup
    assert "docker buildx prune --all --force" in buildkit_cleanup
    assert "docker image inspect --format '{{.Size}}' local/genegalleon:ci" in buildkit_cleanup
    assert "Insufficient disk space for SIF conversion" in buildkit_cleanup
    assert_pinned_singularity_runtime(runtime_install)
    assert "apt-get install -y apptainer" not in runtime_install["run"]
    assert conversion["env"]["APPTAINER_DISABLE_CACHE"] == "1"
    assert conversion["env"]["SINGULARITY_DISABLE_CACHE"] == "1"
    assert "runner.temp" in conversion["env"]["APPTAINER_TMPDIR"]
    assert "runner.temp" in conversion["env"]["SINGULARITY_TMPDIR"]
    assert "runner.temp" in conversion["env"]["TMPDIR"]
    assert "2 * image_bytes + reserve_bytes" in buildkit_cleanup
    assert "always()" in disk_report["if"]
    assert "sif-cache.outputs.cache-hit != 'true'" in disk_report["if"]
    assert "df -h /" in disk_report["run"]
    assert "singularity-transport-tmp" in disk_report["run"]


def test_release_sif_conversion_has_matching_disk_safeguards():
    release_job = load_workflow("release-sif.yml")["jobs"]["publish-manifest-and-sif"]
    runner_cleanup = step_run(release_job, "Reclaim runner disk space for SIF conversion")
    runtime_install = named_step(release_job, "Install Singularity runtime")
    conversion = named_step(release_job, "Build SIF from immutable GHCR tag")
    disk_report = named_step(release_job, "Report SIF conversion disk usage")

    assert "docker system prune --all --force --volumes" in runner_cleanup
    for path in ("/opt/ghc", "/usr/local/lib/android", "/usr/share/dotnet"):
        assert path in runner_cleanup
    assert conversion["env"]["APPTAINER_DISABLE_CACHE"] == "1"
    assert conversion["env"]["SINGULARITY_DISABLE_CACHE"] == "1"
    assert "required_bytes" in conversion["run"]
    assert "Insufficient disk space for SIF conversion" in conversion["run"]
    assert_pinned_singularity_runtime(runtime_install)
    assert "apt-get install -y apptainer" not in runtime_install["run"]
    assert "runner.temp" in conversion["env"]["TMPDIR"]
    assert disk_report["if"] == "always()"
    assert "singularity-transport-tmp" in disk_report["run"]


def test_release_images_and_sif_are_validated_before_publication():
    jobs = load_workflow("release-sif.yml")["jobs"]
    release_resolution = step_run(
        jobs["prepare-release"],
        "Resolve release, image tags, and source revisions",
    )
    image_validation = named_step(jobs["build-platform"], "Validate release runtime contracts")
    conversion = named_step(jobs["publish-manifest-and-sif"], "Build SIF from immutable GHCR tag")
    identity = named_step(jobs["publish-manifest-and-sif"], "Verify exact release SIF identity")
    sif_validation = named_step(jobs["publish-manifest-and-sif"], "Validate release SIF runtime contracts")

    assert '"${IMAGE}@${DIGEST}"' in image_validation["run"]
    assert '--user "$(id -u):$(id -g)"' in image_validation["run"]
    assert 'git show-ref --verify --quiet "refs/tags/${release_tag}"' in release_resolution
    assert 'git rev-list -n 1 "refs/tags/${release_tag}"' in release_resolution
    assert "--gg-suite runtime workflow/tests" in image_validation["run"]
    assert conversion["env"]["TAG"] == "${{ needs.prepare-release.outputs.immutable_tag }}"
    assert "release_tag" not in conversion["env"]["TAG"]
    assert "check_runtime_freshness.sh" in identity["run"]
    assert "--expected-hash" in identity["run"]
    assert "steps.runtime-input.outputs.value" in identity["env"]["EXPECTED_RUNTIME_INPUT"]
    assert "run_in_sif.sh" in sif_validation["run"]
    assert "--gg-suite runtime workflow/tests" in sif_validation["run"]


def test_toolchain_dependent_r_integration_test_runs_in_sif_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_job = jobs["r-script-parse"]
    sif_job = jobs["sif-runtime-validation"]
    r_commands = "\n".join(str(step.get("run", "")) for step in r_job["steps"])
    sif_commands = "\n".join(str(step.get("run", "")) for step in sif_job["steps"])
    integration_test = "workflow/tests/test_orthogroup_copy_number_trait_pgls.R"

    assert not any(str(step.get("uses", "")).startswith("actions/setup-python@") for step in r_job["steps"])
    assert integration_test not in r_commands
    assert f"bash workflow/tests/run_in_sif.sh Rscript {integration_test}" in sif_commands


def test_runtime_python_suite_runs_in_authoritative_sif_job():
    sif_job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    commands = "\n".join(str(step.get("run", "")) for step in sif_job["steps"])

    assert "run_in_sif.sh python -m pytest -q --gg-suite runtime workflow/tests" in commands


def test_runtime_suite_contains_real_owned_upstream_contracts():
    conftest = (Path(__file__).parent / "conftest.py").read_text(encoding="utf-8")
    contracts = (Path(__file__).parent / "test_owned_runtime_contracts.py").read_text(encoding="utf-8")

    assert '"test_owned_runtime_contracts.py"' in conftest
    assert "select_orthofinder_core_species.py" in contracts
    assert "SOURCE_REVISIONS" in contracts


def test_treevis_package_validation_runs_only_in_the_container_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_commands = "\n".join(str(step.get("run", "")) for step in jobs["r-script-parse"]["steps"])
    sif_commands = "\n".join(str(step.get("run", "")) for step in jobs["sif-runtime-validation"]["steps"])

    assert "test_treevis_main.R" not in r_commands
    assert "workflow/tests/check_treevis_package.sh" in sif_commands
    assert "workflow/tests/test_treevis_main.R" in sif_commands


def test_workflow_files_parse_to_job_mappings():
    for path in sorted(Path(GITHUB_WORKFLOWS_DIR).glob("*.yml")):
        workflow = load_workflow(path.name)
        assert isinstance(workflow.get("jobs"), dict)
        assert workflow["jobs"], f"No jobs defined in {path.name}"


def test_shell_job_uses_structural_workflow_and_shell_linters():
    shell_job = load_workflow("tests.yml")["jobs"]["shell-static"]
    install_run = step_run(shell_job, "Install shell and workflow linters")
    actionlint_run = step_run(shell_job, "Validate GitHub Actions workflows")
    shellcheck_run = step_run(shell_job, "ShellCheck all tracked shell scripts")

    assert "github.com/rhysd/actionlint/cmd/actionlint@v1.7.12" in install_run
    assert "actionlint" in actionlint_run
    assert "shellcheck -S warning" in shellcheck_run
    assert "git ls-files -z '*.sh'" in shellcheck_run
    assert "scripts+=(dev)" in shellcheck_run


def test_container_workflows_embed_version_and_shared_build_input_hash():
    cases = (
        ("container-ghcr.yml", "build-and-push", "prepare-build"),
        ("release-sif.yml", "build-platform", "prepare-release"),
    )
    for workflow_name, build_job_name, prepare_job_name in cases:
        jobs = load_workflow(workflow_name)["jobs"]
        build_job = jobs[build_job_name]
        hash_step = named_step(build_job, "Compute build input hash")
        build_step = next(
            step for step in build_job["steps"] if str(step.get("uses", "")).startswith("docker/build-push-action@")
        )
        build_args = str(build_step["with"]["build-args"])

        assert "container/scripts/compute_build_input_hash.sh" in hash_step["run"]
        assert "GG_BUILD_VERSION" in hash_step["env"]
        assert "GG_BUILD_SECURITY_REFRESH_EPOCH" in hash_step["env"]
        assert 'build_input_hash="$(bash container/scripts/compute_build_input_hash.sh ' in hash_step["run"]
        assert 'runtime_input_hash="$(bash container/scripts/compute_build_input_hash.sh --runtime ' in hash_step["run"]
        assert "^[0-9a-f]{64}$" in hash_step["run"]
        assert 'echo "value=$(' not in hash_step["run"]
        assert 'echo "runtime_value=$(' not in hash_step["run"]
        assert "GG_VERSION=${{ needs." + prepare_job_name + ".outputs.gg_version }}" in build_args
        assert "BUILD_INPUT_HASH=${{ steps.build-input.outputs.value }}" in build_args
        assert "RUNTIME_INPUT_HASH=${{ steps.build-input.outputs.runtime_value }}" in build_args
        assert (
            "SECURITY_REFRESH_EPOCH=${{ needs."
            + prepare_job_name
            + ".outputs.security_refresh_epoch }}"
        ) in build_args


def test_scheduled_container_build_compares_resolved_inputs_with_published_image():
    workflow = load_workflow("container-ghcr.yml")
    check_job = workflow["jobs"]["check-trigger"]
    decision = step_run(check_job, "Decide whether a scheduled build is needed")
    runtime_contracts = step_run(
        workflow["jobs"]["build-and-push"],
        "Validate repository-owned runtime contracts",
    )

    assert set(workflow["on"]) == {"schedule", "workflow_dispatch"}
    assert workflow["on"]["schedule"] == [{"cron": "0 19 * * *"}]
    assert named_step(check_job, "Checkout")["uses"].startswith("actions/checkout@")
    assert named_step(check_job, "Set up Docker Buildx")["uses"].startswith(
        "docker/setup-buildx-action@"
    )
    assert named_step(check_job, "Login to GHCR")["uses"].startswith("docker/login-action@")
    assert decision.count("resolve_source_sha https://") == 11
    assert decision.count("resolve_source_sha https://gitlab.com/") == 1
    assert "GG_BUILD_PLATFORMS=linux/amd64" in decision
    assert "GG_BUILD_PLATFORMS=linux/arm64" in decision
    assert decision.count("--print-common-revision") == 2
    assert 'GG_BUILD_VCS_REF="${published_revision}"' in decision
    assert 'security_refresh_epoch="$(date -u +%F)"' in decision
    assert 'compute_build_input_hash.sh "${security_refresh_epoch}"' in decision
    assert (
        check_job["outputs"]["security_refresh_epoch"]
        == "${{ steps.decision.outputs.security_refresh_epoch }}"
    )
    assert "container/scripts/compute_build_input_hash.sh" in decision
    assert "docker buildx imagetools inspect --format '{{json .Image}}'" in decision
    assert "container/scripts/check_published_build_input.py" in decision
    assert "needs.check-trigger.outputs.should_build == '1'" in workflow["jobs"]["prepare-build"]["if"]
    assert "upstream source resolution failed during preflight" in decision
    assert "build-input hash calculation failed during preflight" in decision
    assert "/actions/workflows/container-ghcr.yml/runs" not in decision
    assert "--gg-suite runtime workflow/tests" in runtime_contracts
    assert '"${IMAGE}@${DIGEST}"' in runtime_contracts
