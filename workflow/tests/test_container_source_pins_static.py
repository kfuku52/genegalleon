import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PIN_FILE = REPO_ROOT / "container" / "source_pins.env"
SHA_PIN_NAMES = {
    "GG_PIN_AMALGKIT_REPO_SHA",
    "GG_PIN_CDSKIT_REPO_SHA",
    "GG_PIN_CSUBST_REPO_SHA",
    "GG_PIN_NWKIT_REPO_SHA",
    "GG_PIN_KFL1OU_REPO_SHA",
    "GG_PIN_KFTOOLS_REPO_SHA",
    "GG_PIN_RKFTOOLS_REPO_SHA",
    "GG_PIN_RADTE_REPO_SHA",
}


def read_pins():
    pins = {}
    for raw_line in PIN_FILE.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split("=", 1)
        pins[key] = value
    return pins


def test_all_github_source_revisions_are_pinned_to_commits():
    pins = read_pins()

    assert SHA_PIN_NAMES.issubset(pins)
    for name in SHA_PIN_NAMES:
        assert re.fullmatch(r"[0-9a-f]{40}", pins[name]), name
    assert re.fullmatch(r"\d+\.\d+\.\d+", pins["GG_PIN_CSUBST_VERSION"])


def test_all_container_build_paths_consume_authoritative_source_pins():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    buildx = (REPO_ROOT / "container" / "buildx.sh").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.sh").read_text(
        encoding="utf-8"
    )

    assert "COPY container/source_pins.env /opt/pg/source_pins.env" in dockerfile
    assert dockerfile.count("source /opt/pg/source_pins.env") == 2
    assert 'source "${script_dir}/source_pins.env"' in buildx
    assert 'source "${script_dir}/source_pins.env"' in apptainer
    assert "KFU52_CSUBST_REPO_SHA=\"${csubst_repo_sha}\"" in dockerfile
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
        "genomepy",
        "mysql-connector-python",
        "ortools",
        "protobuf",
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
