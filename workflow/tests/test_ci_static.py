from shell_static_helpers import GITHUB_WORKFLOWS_DIR, read_text


def test_sif_runtime_validation_uses_sha_tag_and_warns_on_latest_fallback():
    workflow = read_text(GITHUB_WORKFLOWS_DIR / "tests.yml")

    assert 'sha_tag="sha-${GITHUB_SHA::7}"' in workflow
    assert 'TAG="${sha_tag}"' in workflow
    assert 'TAG="${GITHUB_SHA::7}"' not in workflow
    assert "::warning::GHCR image ${image}:${sha_tag} was unavailable" in workflow
    assert "falling back to ${image}:latest for dependency-runtime validation" in workflow
    assert 'echo "- warning: commit-specific image ${image}:${sha_tag} was unavailable; latest was used."' in workflow


def test_r_script_job_installs_python_dependencies_for_integration_test():
    workflow = read_text(GITHUB_WORKFLOWS_DIR / "tests.yml")
    r_job = workflow.split("\n  r-script-parse:\n", 1)[1].split("\n  sif-runtime-validation:\n", 1)[0]

    assert "uses: actions/setup-python@v5" in r_job
    assert 'python-version: "3.12"' in r_job
    assert "cache-dependency-path: workflow/tests/requirements.txt" in r_job
    assert "python -m pip install -r workflow/tests/requirements.txt" in r_job
    assert r_job.index("python -m pip install -r workflow/tests/requirements.txt") < r_job.index(
        "Rscript workflow/tests/test_orthogroup_copy_number_trait_pgls.R"
    )
