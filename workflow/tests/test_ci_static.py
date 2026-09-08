import json
import os
import re
import subprocess
from pathlib import Path

import pytest
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
    "oras-project/setup-oras",
    "r-lib/actions/setup-r",
    "r-lib/actions/setup-r-dependencies",
}

SINGULARITY_CE_VERSION = "4.5.0"
SINGULARITY_CE_DEB_SHA256 = "85e6f7af5e7aad5b1bf28183ce333998bd37eb8f4769af352c47a5153f3373fb"
LOCAL_ACTIONS = {"./.github/actions/validate-sif"}


def validation_manifest():
    return json.loads(Path(__file__).with_name("validation_manifest.json").read_text())


def load_workflow(name: str) -> dict:
    path = GITHUB_WORKFLOWS_DIR / name
    workflow = yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    for job in workflow.get("jobs", {}).values():
        expanded = []
        for step in job.get("steps", []):
            expanded.append(step)
            if step.get("uses") in LOCAL_ACTIONS:
                action_path = GITHUB_WORKFLOWS_DIR.parent.parent / step["uses"] / "action.yml"
                action = yaml.load(action_path.read_text(), Loader=yaml.BaseLoader)
                assert action["runs"]["using"] == "composite"
                expanded.extend(action["runs"]["steps"])
        job["steps"] = expanded
    return workflow


def workflow_steps(workflow: dict):
    for job in workflow.get("jobs", {}).values():
        yield from job.get("steps", [])


def named_step(job: dict, name: str) -> dict:
    matches = [step for step in job.get("steps", []) if step.get("name") == name]
    assert len(matches) == 1, f"Expected exactly one CI step named {name!r}"
    return matches[0]


def step_run(job: dict, name: str) -> str:
    return str(named_step(job, name).get("run", ""))


@pytest.mark.parametrize("workflow_name,job_name,step_name", [
    ("container-ghcr.yml", "prepare-build", "Resolve image tags"),
    ("release-sif.yml", "prepare-release", "Resolve release, image tags, and source revisions"),
])
@pytest.mark.parametrize("resolution_fails", [False, True])
def test_publish_preparation_forwards_one_complete_source_snapshot(
    tmp_path, workflow_name, job_name, step_name, resolution_fails,
):
    variables = (
        "KFU52_AMALGKIT_REPO_SHA", "KFU52_CDSKIT_REPO_SHA", "KFU52_CSUBST_REPO_SHA",
        "KFU52_NWKIT_REPO_SHA", "BUSCO_REPO_SHA", "PAML_REPO_SHA", "KFL1OU_REPO_SHA",
        "KFFRACTBIAS_REPO_SHA", "KFTOOLS_REPO_SHA", "RKFTOOLS_REPO_SHA", "RADTE_REPO_SHA",
    )
    values = {variable: f"{index:040x}" for index, variable in enumerate(variables, start=1)}
    scripts = tmp_path / "container/scripts"
    scripts.mkdir(parents=True)
    resolver = scripts / "resolve_source_revisions.sh"
    resolver.write_text(
        "printf '%s\\n' " + " ".join(f"'{variable}={value}'" for variable, value in values.items())
        + f"\nexit {int(resolution_fails)}\n"
    )
    coverage = scripts / "check_env_coverage.sh"
    coverage.write_text("#!/bin/bash\nexit 0\n")
    coverage.chmod(0o755)
    (tmp_path / "VERSION").write_text("1.2.3\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    git = bin_dir / "git"
    git.write_text('#!/bin/bash\nif [[ "$1" == rev-list ]]; then printf "%040d\\n" 1; fi\nexit 0\n')
    git.chmod(0o755)
    output = tmp_path / "outputs"
    env = os.environ.copy()
    env.update({
        "PATH": f"{bin_dir}{os.pathsep}{env['PATH']}",
        "GITHUB_OUTPUT": str(output), "GITHUB_REPOSITORY_OWNER": "example",
        "GITHUB_REPOSITORY": "example/genegalleon", "GITHUB_SHA": "a" * 40,
        "GITHUB_EVENT_NAME": "release", "RELEASE_EVENT_TAG": "v1.2.3",
        "SECURITY_REFRESH_EPOCH": "2026-08-31",
    })
    step = step_run(load_workflow(workflow_name)["jobs"][job_name], step_name)
    completed = subprocess.run(["bash", "-c", step], cwd=tmp_path, env=env, capture_output=True, text=True)
    if resolution_fails:
        assert completed.returncode != 0
        assert not output.exists()
    else:
        assert completed.returncode == 0, completed.stderr
        actual = dict(line.split("=", 1) for line in output.read_text().splitlines())
        for variable, value in values.items():
            assert actual[variable.removeprefix("KFU52_").lower()] == value


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
        if value in LOCAL_ACTIONS:
            continue
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
    assert "io.genegalleon.build-target" in selection_step["run"]
    assert "io.genegalleon.security-refresh-epoch" in selection_step["run"]
    assert "--expected-value" in selection_step["run"]
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


def test_sif_runtime_validation_builds_current_dependency_corrected_runtime():
    sif_job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    validation = named_step(sif_job, "Validate exact runtime and reuse the shared SIF cache")

    assert "with" not in validation


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
    release_job = load_workflow("release-sif.yml")["jobs"]["build-platform"]
    runner_cleanup = step_run(release_job, "Reclaim amd64 runner disk space for SIF conversion")
    runtime_install = named_step(release_job, "Install Singularity runtime for amd64 SIF")
    conversion = named_step(release_job, "Build amd64 SIF from exact platform digest")
    disk_report = named_step(release_job, "Report amd64 SIF conversion disk usage")

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
    assert "always()" in disk_report["if"]
    assert "matrix.platform == 'linux/amd64'" in conversion["if"]
    assert '"docker://${IMAGE}@${DIGEST}"' in conversion["run"]
    assert "singularity-transport-tmp" in disk_report["run"]


def test_release_images_and_sif_are_validated_before_publication():
    jobs = load_workflow("release-sif.yml")["jobs"]
    release_resolution = step_run(
        jobs["prepare-release"],
        "Resolve release, image tags, and source revisions",
    )
    image_validation = named_step(jobs["build-platform"], "Validate release runtime contracts")
    conversion = named_step(jobs["build-platform"], "Build amd64 SIF from exact platform digest")
    identity = named_step(jobs["build-platform"], "Verify exact amd64 SIF identity")
    sif_validation = named_step(jobs["build-platform"], "Validate amd64 SIF runtime contracts")
    handoff = named_step(jobs["publish-manifest-and-sif"], "Verify qualified amd64 SIF handoff")

    assert '"${IMAGE}@${DIGEST}"' in image_validation["run"]
    assert '--user "$(id -u):$(id -g)"' in image_validation["run"]
    assert 'git show-ref --verify --quiet "refs/tags/${release_tag}"' in release_resolution
    assert 'git rev-list -n 1 "refs/tags/${release_tag}"' in release_resolution
    assert "workflow/tests/run_checks.py runtime" in image_validation["run"]
    assert conversion["env"]["DIGEST"] == "${{ steps.build.outputs.digest }}"
    assert "@${DIGEST}" in conversion["run"]
    assert "check_runtime_freshness.sh" in identity["run"]
    assert "--expected-hash" in identity["run"]
    assert "steps.build-input.outputs.runtime_value" in identity["env"]["EXPECTED_RUNTIME_INPUT"]
    assert "run_in_sif.sh" in sif_validation["run"]
    assert "workflow/tests/run_checks.py runtime" in sif_validation["run"]
    assert "sha256sum --check --strict" in handoff["run"]


def test_release_sif_is_published_early_as_durable_content_addressed_oci():
    jobs = load_workflow("release-sif.yml")["jobs"]
    build = jobs["build-platform"]
    publish = named_step(build, "Publish content-addressed amd64 SIF artifact")
    handoff = named_step(build, "Upload qualified SIF handoff")
    docker_build = next(
        step for step in build["steps"]
        if str(step.get("uses", "")).startswith("docker/build-push-action@")
    )

    assert "cache-to" not in docker_build["with"]
    assert "matrix.platform == 'linux/amd64'" in publish["if"]
    assert "oras push" in publish["run"]
    assert "sif-sha256-${sif_sha256}" in publish["run"]
    assert "oras resolve" in publish["run"]
    assert "@${manifest_digest}" in publish["run"]
    assert handoff["with"]["retention-days"] == "1"
    publish_job_text = json.dumps(jobs["publish-manifest-and-sif"], sort_keys=True)
    assert "singularity build" not in publish_job_text
    assert "Upload large SIF as workflow artifact" not in publish_job_text


def test_toolchain_dependent_r_integration_test_runs_in_sif_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_job = jobs["r-script-parse"]
    sif_job = jobs["sif-runtime-validation"]
    r_commands = "\n".join(str(step.get("run", "")) for step in r_job["steps"])
    sif_commands = "\n".join(str(step.get("run", "")) for step in sif_job["steps"])
    integration_test = "workflow/tests/test_orthogroup_copy_number_trait_pgls.R"

    assert not any(str(step.get("uses", "")).startswith("actions/setup-python@") for step in r_job["steps"])
    assert integration_test not in r_commands
    assert "run_in_sif.sh python workflow/tests/run_checks.py runtime" in sif_commands
    assert ["Rscript", integration_test] in validation_manifest()["r_commands"]


def test_runtime_python_suite_runs_in_authoritative_sif_job():
    sif_job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    commands = "\n".join(str(step.get("run", "")) for step in sif_job["steps"])

    assert "run_in_sif.sh python workflow/tests/run_checks.py runtime" in commands


def test_runtime_suite_contains_real_owned_upstream_contracts():
    contracts = (Path(__file__).parent / "test_owned_runtime_contracts.py").read_text(encoding="utf-8")

    assert "test_owned_runtime_contracts.py" in validation_manifest()["runtime_python_files"]
    assert "select_orthofinder_core_species.py" in contracts
    assert "SOURCE_REVISIONS" in contracts


def test_runtime_suite_guards_aster_and_ortools_loader_compatibility():
    assert "test_aster_ortools_runtime.py" in validation_manifest()["runtime_python_files"]


def test_treevis_package_validation_runs_only_in_the_container_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_commands = "\n".join(str(step.get("run", "")) for step in jobs["r-script-parse"]["steps"])
    sif_commands = "\n".join(str(step.get("run", "")) for step in jobs["sif-runtime-validation"]["steps"])

    assert "test_treevis_main.R" not in r_commands
    assert "workflow/tests/run_checks.py runtime" in sif_commands
    assert ["bash", "workflow/tests/check_treevis_package.sh"] in validation_manifest()["r_commands"]
    assert ["Rscript", "workflow/tests/test_treevis_main.R"] in validation_manifest()["r_commands"]


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
    assert "check_composite_actions.py" in actionlint_run
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
    assert decision.count("resolve_source_revisions.sh --format env --scope all") == 1
    assert "if ! resolved_sources=" in decision
    assert 'export "${variable}=${revision}"' in decision
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
    assert "workflow/tests/run_checks.py runtime" in runtime_contracts
    assert '"${IMAGE}@${DIGEST}"' in runtime_contracts


def test_daily_publisher_primes_the_same_exact_sif_cache_used_by_commit_checks():
    daily = load_workflow("container-ghcr.yml")["jobs"]
    prime = daily["prime-sif-cache"]
    commit = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    prepared = named_step(prime, "Prepare validated SIF for subsequent commit checks")
    assert prepared["uses"] == "./.github/actions/validate-sif"
    assert "immutable_tag" in prepared["with"]["image-ref"]
    assert prepared["with"]["runtime-input"] == "${{ needs.build-and-push.outputs.runtime_input }}"
    assert "steps.build-input.outputs.runtime_value" in daily["build-and-push"]["outputs"]["runtime_input"]
    assert named_step(prime, "Restore exact validation SIF")["with"] == named_step(
        commit, "Restore exact validation SIF"
    )["with"]
    assert "resolve_source_revisions.sh" not in step_run(prime, "Select published or local validation image")
    for job in (prime, commit):
        save = named_step(job, "Save validated SIF by exact runtime input")
        assert "github.event_name != 'pull_request'" in save["if"]
        assert "github.event.repository.default_branch" in save["if"]
        assert "restore-keys" not in named_step(job, "Restore exact validation SIF")["with"]


def test_parallel_python_lanes_install_the_same_prebuilt_offline_wheels():
    jobs = load_workflow("tests.yml")["jobs"]
    wheel_job = jobs["python-wheels"]
    artifact = named_step(wheel_job, "Share exact wheels with parallel test lanes")["with"]
    assert artifact["retention-days"] == "1"
    key = named_step(wheel_job, "Restore test wheels by resolved source and constraints")["with"]
    assert "steps.source.outputs.csubst_sha" in key["key"]
    assert "requirements.lock.txt" in key["key"]
    assert "restore-keys" not in key
    save = named_step(wheel_job, "Save trusted test wheels")
    assert "event_name != 'pull_request'" in save["if"]
    assert "github.event.repository.default_branch" in save["if"]
    for lane in ("python-fast", "python-heavy"):
        assert "python-wheels" in jobs[lane]["needs"]
        download = named_step(jobs[lane], "Download exact test wheels")["with"]
        assert download["name"] == artifact["name"]
        install = step_run(jobs[lane], "Install test dependencies")
        assert "--no-index" in install
        assert "--find-links" in install
        assert "install-requirements.txt" in install
        assert "git+" not in install


@pytest.mark.parametrize("invalid", [None, "mutable_tag", "sha_alias", "missing_hash", "missing_image", "unsafe_image"])
def test_prepared_sif_metadata_is_validated_without_resolving_new_upstream_tips(tmp_path, invalid):
    job = load_workflow("container-ghcr.yml")["jobs"]["prime-sif-cache"]
    run = step_run(job, "Resolve exact validation runtime input")
    output = tmp_path / "github-output"
    env = os.environ | {
        "PREPARED_IMAGE_REF": "ghcr.io/example/genegalleon:20260831-abcdef0-fedcba987654",
        "PREPARED_RUNTIME_INPUT": "a" * 64,
        "PREPARED_SECURITY_EPOCH": "2026-08-31",
        "GITHUB_OUTPUT": str(output),
    }
    if invalid == "mutable_tag":
        env["PREPARED_IMAGE_REF"] = "ghcr.io/example/genegalleon:latest"
    elif invalid == "sha_alias":
        env["PREPARED_IMAGE_REF"] = "ghcr.io/example/genegalleon:sha-abcdef0"
    elif invalid == "missing_hash":
        env["PREPARED_RUNTIME_INPUT"] = ""
    elif invalid == "missing_image":
        env["PREPARED_IMAGE_REF"] = ""
    elif invalid == "unsafe_image":
        env["PREPARED_IMAGE_REF"] += "\nsource=local"
    # This empty cwd has no resolver scripts: a valid prepared snapshot must
    # succeed without querying a potentially newer branch or rebuilding.
    result = subprocess.run(["bash", "-c", run], cwd=tmp_path, env=env,
                            capture_output=True, text=True, check=False)
    if invalid:
        assert result.returncode != 0
        assert not output.exists()
    else:
        assert result.returncode == 0, result.stderr
        assert "value=" + "a" * 64 in output.read_text()
        assert "published_tag=20260831-abcdef0-fedcba987654" in output.read_text()


def test_sif_checks_skip_re_resolution_only_after_verifying_exact_identity():
    job = load_workflow("tests.yml")["jobs"]["sif-runtime-validation"]
    steps = job["steps"]
    identity = named_step(job, "Verify exact validation SIF identity")
    validate = named_step(job, "Run SIF validation checks")
    assert steps.index(identity) < steps.index(validate)
    assert "--expected-hash" in identity["run"]
    assert validate["env"]["GG_RUNTIME_FRESHNESS"] == "off"
