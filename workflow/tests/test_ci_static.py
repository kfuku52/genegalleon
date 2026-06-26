from shell_static_helpers import GITHUB_WORKFLOWS_DIR, read_text


def test_sif_runtime_validation_uses_sha_tag_and_warns_on_latest_fallback():
    workflow = read_text(GITHUB_WORKFLOWS_DIR / "tests.yml")

    assert 'sha_tag="sha-${GITHUB_SHA::7}"' in workflow
    assert 'TAG="${sha_tag}"' in workflow
    assert 'TAG="${GITHUB_SHA::7}"' not in workflow
    assert "::warning::GHCR image ${image}:${sha_tag} was unavailable" in workflow
    assert "falling back to ${image}:latest for dependency-runtime validation" in workflow
    assert 'echo "- warning: commit-specific image ${image}:${sha_tag} was unavailable; latest was used."' in workflow
