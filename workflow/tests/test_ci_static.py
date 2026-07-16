from pathlib import Path

import yaml
from shell_static_helpers import GITHUB_WORKFLOWS_DIR


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


def test_workflows_use_current_checkout_and_setup_python_majors():
    uses = {
        step["uses"]
        for path in sorted(GITHUB_WORKFLOWS_DIR.glob("*.yml"))
        for step in workflow_steps(load_workflow(path.name))
        if "uses" in step
    }

    assert "actions/checkout@v6" in uses
    assert "actions/setup-python@v6" in uses
    assert not {value for value in uses if value.startswith("actions/checkout@") and value != "actions/checkout@v6"}
    assert not {
        value for value in uses if value.startswith("actions/setup-python@") and value != "actions/setup-python@v6"
    }


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

        assert 'source_revision_hash="$(' in metadata_run
        assert "| sha256sum | cut -c1-12" in metadata_run
        assert 'immutable_tag="${date_tag}-${short_sha}-${source_revision_hash}"' in metadata_run
        assert 'imagetools inspect "${IMAGE}:${IMMUTABLE_TAG}"' in publish_run
        assert "Refusing to overwrite immutable image tag" in publish_run


def test_sif_runtime_validation_builds_the_current_repository_image():
    jobs = load_workflow("tests.yml")["jobs"]
    sif_job = jobs["sif-runtime-validation"]
    build_run = step_run(sif_job, "Build current GeneGalleon image")
    conversion_run = step_run(sif_job, "Build validation SIF from current image")

    assert "IMAGE_SOURCE=local" in build_run
    assert "IMAGE=local/genegalleon" in build_run
    assert "TAG=ci" in build_run
    assert "FORCE_LOAD_IN_CI=1" in build_run
    assert "bash ./gg_container_build_entrypoint.sh" in build_run
    assert "SOURCE=docker-daemon" in conversion_run
    assert "latest" not in conversion_run.lower()


def test_toolchain_dependent_r_integration_test_runs_in_sif_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_job = jobs["r-script-parse"]
    sif_job = jobs["sif-runtime-validation"]
    r_commands = "\n".join(str(step.get("run", "")) for step in r_job["steps"])
    sif_commands = "\n".join(str(step.get("run", "")) for step in sif_job["steps"])
    integration_test = "workflow/tests/test_orthogroup_copy_number_trait_pgls.R"

    assert not any(step.get("uses") == "actions/setup-python@v6" for step in r_job["steps"])
    assert integration_test not in r_commands
    assert f"bash workflow/tests/run_in_sif.sh Rscript {integration_test}" in sif_commands


def test_treevis_package_validation_runs_only_in_the_container_job():
    jobs = load_workflow("tests.yml")["jobs"]
    r_commands = "\n".join(str(step.get("run", "")) for step in jobs["r-script-parse"]["steps"])
    sif_commands = "\n".join(
        str(step.get("run", "")) for step in jobs["sif-runtime-validation"]["steps"]
    )

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
    shellcheck_run = step_run(shell_job, "ShellCheck split workflow drivers")

    assert "github.com/rhysd/actionlint/cmd/actionlint@v1.7.12" in install_run
    assert "actionlint" in actionlint_run
    assert "shellcheck -S warning" in shellcheck_run
