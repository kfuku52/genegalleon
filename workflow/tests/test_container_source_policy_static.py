import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PROGRAM_SHA_VARS = (
    "KFU52_AMALGKIT_REPO_SHA",
    "KFU52_CDSKIT_REPO_SHA",
    "KFU52_CSUBST_REPO_SHA",
    "KFU52_NWKIT_REPO_SHA",
    "BUSCO_REPO_SHA",
    "PAML_REPO_SHA",
    "KFL1OU_REPO_SHA",
    "KFFRACTBIAS_REPO_SHA",
    "KFTOOLS_REPO_SHA",
    "RKFTOOLS_REPO_SHA",
    "RADTE_REPO_SHA",
)


def test_program_source_defaults_are_moving_branches_not_commit_pins():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(encoding="utf-8")
    branches_path = REPO_ROOT / "container" / "source_branches.env"
    branches = branches_path.read_text(encoding="utf-8")

    assert not (REPO_ROOT / "container" / "source_pins.env").exists()
    assert 'source "${script_dir}/source_branches.env"' in buildx
    assert 'source "${script_dir}/source_branches.env"' in apptainer
    assert "GG_PIN_" not in branches
    assert "_REPO_SHA=" not in branches
    branch_assignments = re.findall(r"^GG_SOURCE_[A-Z0-9_]+_REPO_REF=(\S+)$", branches, re.MULTILINE)
    assert len(branch_assignments) == len(PROGRAM_SHA_VARS)
    assert set(branch_assignments) <= {"main", "master"}

    for sha_var in PROGRAM_SHA_VARS:
        assert f'ARG {sha_var}=""' in dockerfile
        assert f"{sha_var}=${{{sha_var}:-}}" in buildx
        assert f"{sha_var}=${{{sha_var}:-}}" in apptainer

    policy_files = [
        REPO_ROOT / "container" / "Dockerfile",
        REPO_ROOT / "container" / "buildx.sh",
        REPO_ROOT / "container" / "apptainer_local_build.sh",
        REPO_ROOT / "container" / "scripts" / "compute_build_input_hash.sh",
        REPO_ROOT / ".github" / "workflows" / "container-ghcr.yml",
        REPO_ROOT / ".github" / "workflows" / "release-sif.yml",
    ]
    for path in policy_files:
        text = path.read_text(encoding="utf-8")
        assert "source_pins" not in text
        assert "GG_PIN_" not in text
        assert not re.search(r'(?:ARG\s+)?[A-Z0-9_]+_REPO_SHA=(?:"?)[0-9a-f]{40}(?:"?)', text)


def test_all_container_build_paths_resolve_one_snapshot_per_build():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(encoding="utf-8")

    assert 'ARG KFU52_REPO_REF=""' in dockerfile
    assert 'ARG KFU52_AMALGKIT_REPO_REF=""' in dockerfile
    assert 'ARG KFU52_CSUBST_REPO_REF=""' in dockerfile
    assert 'ARG KFL1OU_REPO_REF=""' in dockerfile
    assert buildx.count("resolve_source_sha ") == 11
    assert apptainer.count("resolve_source_sha ") == 11
    for source in ("amalgkit", "cdskit", "csubst", "nwkit", "BUSCO", "paml", "kfl1ou", "kfFractBias", "kftools", "rkftools", "RADTE"):
        assert f" {source}" in buildx
        assert f" {source}" in apptainer
    assert "> /opt/pg/logs/source_revisions.tsv" in dockerfile


def test_container_build_paths_share_python_compatibility_constraints():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    apptainer_build = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(encoding="utf-8")
    requirements_path = REPO_ROOT / "container" / "pip-compatibility.requirements.txt"
    requirements = requirements_path.read_text(encoding="utf-8")

    assert "-r /opt/pg/pip-compatibility.requirements.txt" in dockerfile
    assert "-r /opt/pg/pip-compatibility.requirements.txt" in apptainer_template
    assert (
        "@@STAGING_ROOT@@/pip-compatibility.requirements.txt /opt/pg/pip-compatibility.requirements.txt"
    ) in apptainer_template
    assert ('cp "${repo_root}/container/pip-compatibility.requirements.txt" "${staging_root}/"') in apptainer_build
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
    assert "pypdf>=6.16.1" in requirements
    assert "setuptools<83" in requirements


def test_container_build_paths_refresh_system_security_packages_daily():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    build_hash = (
        REPO_ROOT / "container" / "scripts" / "compute_build_input_hash.sh"
    ).read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    apptainer_build = (
        REPO_ROOT / "container" / "apptainer_local_build.sh"
    ).read_text(encoding="utf-8")

    assert 'ARG SECURITY_REFRESH_EPOCH=""' in dockerfile
    assert "apt-get upgrade --with-new-pkgs -y" in dockerfile
    assert "apt-get --simulate dist-upgrade" in dockerfile
    assert "Upgradeable system packages remain after the security refresh" in dockerfile
    assert 'io.genegalleon.security-refresh-epoch="${SECURITY_REFRESH_EPOCH}"' in dockerfile
    assert "gg_version=${version};security_refresh_epoch=${1}" in build_hash
    assert 'SECURITY_REFRESH_EPOCH="$(date -u +%F)"' in buildx
    assert '--build-arg SECURITY_REFRESH_EPOCH="${SECURITY_REFRESH_EPOCH}"' in buildx
    assert 'SECURITY_REFRESH_EPOCH="$(date -u +%F)"' in apptainer_build
    assert "s|@@SECURITY_REFRESH_EPOCH@@|" in apptainer_build
    assert "io.genegalleon.security-refresh-epoch @@SECURITY_REFRESH_EPOCH@@" in apptainer_template
    assert "security_refresh_epoch.txt" in dockerfile
    assert "security_refresh_epoch.txt" in apptainer_template


def test_python_ci_follows_current_csubst_branch_without_a_commit_pin():
    requirements = (REPO_ROOT / "workflow/tests/requirements.txt").read_text(encoding="utf-8")
    lock = (REPO_ROOT / "workflow/tests/requirements.lock.txt").read_text(encoding="utf-8")

    assert "csubst @ git+https://github.com/kfuku52/csubst.git@master" in requirements
    assert not re.search(r"csubst[^\n]*@[0-9a-f]{40}", requirements)
    assert not re.search(r"csubst[^\n]*@[0-9a-f]{40}", lock)


def test_native_apptainer_build_records_source_revisions():
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")

    assert "> /opt/pg/logs/source_revisions.tsv" in apptainer_template
    for source in (
        "amalgkit",
        "cdskit",
        "csubst",
        "nwkit",
        "BUSCO",
        "paml",
        "kfl1ou",
        "kfFractBias",
        "kftools",
        "rkftools",
        "RADTE",
    ):
        assert source in apptainer_template


def test_treevis_package_is_part_of_every_repository_owned_image_context():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    apptainer_build = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(encoding="utf-8")
    publish_workflow = (REPO_ROOT / ".github" / "workflows" / "container-ghcr.yml").read_text(encoding="utf-8")

    assert "COPY workflow/support/treevis /opt/pg/src/genegalleon.treevis" in dockerfile
    assert "R CMD INSTALL /opt/pg/src/genegalleon.treevis" in dockerfile
    build_hash = (REPO_ROOT / "container" / "scripts" / "compute_build_input_hash.sh").read_text(encoding="utf-8")
    assert "compute_build_input_hash.sh" in buildx
    assert "find workflow/support/treevis -type f -print" in build_hash
    assert 'cp -R "${repo_root}/workflow/support/treevis" "${staging_root}/"' in apptainer_build
    assert "@@STAGING_ROOT@@/treevis /opt/pg/src/genegalleon.treevis" in apptainer_template
    assert "container/scripts/compute_build_input_hash.sh" in publish_workflow


def test_container_build_paths_pin_the_same_base_image_digest():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    base_pattern = re.compile(r"mambaorg/micromamba:noble@sha256:[0-9a-f]{64}")

    docker_base = base_pattern.search(dockerfile)
    apptainer_base = base_pattern.search(apptainer_template)
    assert docker_base is not None
    assert apptainer_base is not None
    assert docker_base.group(0) == apptainer_base.group(0)


def test_container_build_paths_include_rar_extraction_runtime():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    required_commands = (REPO_ROOT / "container" / "spec" / "required_commands.tsv").read_text(encoding="utf-8")
    arm64_required_commands = (REPO_ROOT / "container" / "spec" / "required_commands.arm64.tsv").read_text(
        encoding="utf-8"
    )

    assert "libarchive-tools" in dockerfile
    assert "libarchive-tools" in apptainer_template
    assert "base\tbsdtar" in required_commands
    assert "base\tbsdtar" in arm64_required_commands
