import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
KFUKU52_SHA_VARS = (
    "KFU52_AMALGKIT_REPO_SHA",
    "KFU52_CDSKIT_REPO_SHA",
    "KFU52_CSUBST_REPO_SHA",
    "KFU52_NWKIT_REPO_SHA",
    "KFL1OU_REPO_SHA",
    "KFTOOLS_REPO_SHA",
    "RKFTOOLS_REPO_SHA",
    "RADTE_REPO_SHA",
)


def test_kfuku52_sources_have_no_fixed_default_sha():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(
        encoding="utf-8"
    )

    assert not (REPO_ROOT / "container" / "source_pins.env").exists()
    assert "source_pins.env" not in dockerfile
    assert "GG_PIN_" not in dockerfile
    assert "GG_PIN_" not in buildx
    assert "GG_PIN_" not in apptainer
    for sha_var in KFUKU52_SHA_VARS:
        assert f'ARG {sha_var}=""' in dockerfile
        assert f"{sha_var}=${{{sha_var}:-}}" in buildx
        assert f"{sha_var}=${{{sha_var}:-}}" in apptainer


def test_all_container_build_paths_resolve_current_upstream_revisions():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(
        encoding="utf-8"
    )

    assert 'ARG KFU52_REPO_REF="master"' in dockerfile
    assert 'ARG KFU52_AMALGKIT_REPO_REF="master"' in dockerfile
    assert 'ARG KFU52_CSUBST_REPO_REF="master"' in dockerfile
    assert 'ARG KFL1OU_REPO_REF="main"' in dockerfile
    assert "KFL1OU_REPO_REF=${KFL1OU_REPO_REF:-main}" in buildx
    assert "KFL1OU_REPO_REF=${KFL1OU_REPO_REF:-main}" in apptainer
    assert buildx.count("resolve_source_sha ") == 8
    assert apptainer.count("resolve_source_sha ") == 8
    for source in ("amalgkit", "cdskit", "csubst", "nwkit", "kfl1ou", "kftools", "rkftools", "RADTE"):
        assert f" {source}" in buildx
        assert f" {source}" in apptainer
    assert "> /opt/pg/logs/source_revisions.tsv" in dockerfile


def test_container_build_paths_share_python_compatibility_constraints():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    apptainer_build = (
        REPO_ROOT / "container" / "apptainer_local_build.sh"
    ).read_text(encoding="utf-8")
    requirements_path = REPO_ROOT / "container" / "pip-compatibility.requirements.txt"
    requirements = requirements_path.read_text(encoding="utf-8")

    assert "-r /opt/pg/pip-compatibility.requirements.txt" in dockerfile
    assert "-r /opt/pg/pip-compatibility.requirements.txt" in apptainer_template
    assert (
        "@@STAGING_ROOT@@/pip-compatibility.requirements.txt "
        "/opt/pg/pip-compatibility.requirements.txt"
    ) in apptainer_template
    assert (
        'cp "${repo_root}/container/pip-compatibility.requirements.txt" '
        '"${staging_root}/"'
    ) in apptainer_build
    assert "python -m pip check" in dockerfile
    assert "python -m pip check" in apptainer_template
    for package in (
        "pillow",
        "lxml",
        "ete4",
        "genomepy",
        "mysql-connector-python",
        "ortools",
        "protobuf",
        "defusedxml",
    ):
        assert f"{package}==" in requirements


def test_native_apptainer_build_records_source_revisions():
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")

    assert "> /opt/pg/logs/source_revisions.tsv" in apptainer_template
    for source in (
        "amalgkit",
        "cdskit",
        "csubst",
        "nwkit",
        "kfl1ou",
        "kftools",
        "rkftools",
        "RADTE",
    ):
        assert source in apptainer_template


def test_treevis_package_is_part_of_every_repository_owned_image_context():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    apptainer_build = (
        REPO_ROOT / "container" / "apptainer_local_build.sh"
    ).read_text(encoding="utf-8")
    publish_workflow = (
        REPO_ROOT / ".github" / "workflows" / "container-ghcr.yml"
    ).read_text(encoding="utf-8")

    assert "COPY workflow/support/treevis /opt/pg/src/genegalleon.treevis" in dockerfile
    assert "R CMD INSTALL /opt/pg/src/genegalleon.treevis" in dockerfile
    assert "find workflow/support/treevis -type f -print" in buildx
    assert 'cp -R "${repo_root}/workflow/support/treevis" "${staging_root}/"' in apptainer_build
    assert "@@STAGING_ROOT@@/treevis /opt/pg/src/genegalleon.treevis" in apptainer_template
    assert '"workflow/support/treevis/"' in publish_workflow


def test_container_build_paths_pin_the_same_base_image_digest():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    base_pattern = re.compile(
        r"mambaorg/micromamba:noble@sha256:[0-9a-f]{64}"
    )

    docker_base = base_pattern.search(dockerfile)
    apptainer_base = base_pattern.search(apptainer_template)
    assert docker_base is not None
    assert apptainer_base is not None
    assert docker_base.group(0) == apptainer_base.group(0)


def test_container_build_paths_include_rar_extraction_runtime():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    required_commands = (
        REPO_ROOT / "container" / "spec" / "required_commands.tsv"
    ).read_text(encoding="utf-8")
    arm64_required_commands = (
        REPO_ROOT / "container" / "spec" / "required_commands.arm64.tsv"
    ).read_text(encoding="utf-8")

    assert "libarchive-tools" in dockerfile
    assert "libarchive-tools" in apptainer_template
    assert "base\tbsdtar" in required_commands
    assert "base\tbsdtar" in arm64_required_commands
