import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

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


def declared_context_inputs():
    return subprocess.check_output(
        [sys.executable, str(REPO_ROOT / "container/scripts/list_build_inputs.py")], text=True
    ).splitlines()


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
    for wrapper in (buildx, apptainer):
        assert "resolve_source_revisions.sh\" --format env --scope all" in wrapper
        assert "resolve_source_sha " not in wrapper
    installer = (REPO_ROOT / "container/scripts/install_source_artifacts.sh").read_text()
    assert "> /opt/pg/logs/source_revisions.tsv" in installer
    assert "install_source_artifacts.sh" in dockerfile


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
    assert "--stage-native-context" in apptainer_build
    assert "container/pip-compatibility.requirements.txt" in declared_context_inputs()
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
    assert '"security_refresh_epoch=${security_refresh_epoch}"' in build_hash
    assert 'hash_mode="runtime"' in build_hash
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

    assert "install_source_artifacts.sh" in apptainer_template
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


def test_shared_source_resolver_preserves_exact_overrides_and_owned_scope():
    resolver = REPO_ROOT / "container" / "scripts" / "resolve_source_revisions.sh"
    env = os.environ.copy()
    for index, variable in enumerate(PROGRAM_SHA_VARS, start=1):
        env[variable] = f"{index:040x}"

    all_sources = subprocess.run(
        ["bash", str(resolver), "--format", "env", "--scope", "all"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert all_sources.returncode == 0, all_sources.stderr
    assert all_sources.stdout.splitlines() == [f"{variable}={env[variable]}" for variable in PROGRAM_SHA_VARS]

    owned_sources = subprocess.run(
        ["bash", str(resolver), "--format", "tsv", "--scope", "owned"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert owned_sources.returncode == 0, owned_sources.stderr
    assert owned_sources.stdout.splitlines()[0] == "source\trevision"
    assert "BUSCO\t" not in owned_sources.stdout
    assert "paml\t" not in owned_sources.stdout
    assert "nwkit\t" in owned_sources.stdout
    assert "csubst\t" in owned_sources.stdout


@pytest.mark.parametrize("fail_source", ["", "csubst"])
def test_source_resolution_runs_concurrently_and_publishes_only_complete_snapshots(tmp_path, fail_source):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    barrier = tmp_path / "barrier"
    barrier.mkdir()
    git = bin_dir / "git"
    # Every lookup waits for all eleven to start. A serial implementation fails
    # at this barrier instead of relying on a timing-sensitive speed assertion.
    git.write_text(f"""#!{sys.executable}
import hashlib
import os
from pathlib import Path
import sys
import time

assert sys.argv[1:3] == ["ls-remote", "--exit-code"]
name = sys.argv[3].rsplit("/", 1)[-1].removesuffix(".git")
barrier = Path(os.environ["MOCK_GIT_BARRIER"])
(barrier / name).touch()
deadline = time.monotonic() + 10
while len(list(barrier.iterdir())) != 11:
    if time.monotonic() >= deadline:
        sys.exit("source lookups did not start concurrently")
    time.sleep(0.01)
if name == os.environ["MOCK_GIT_FAIL_SOURCE"]:
    sys.exit(1)
print(hashlib.sha1(name.encode()).hexdigest(), sys.argv[4], sep="\\t")
""")
    git.chmod(0o755)
    env = os.environ.copy()
    for variable in PROGRAM_SHA_VARS:
        env.pop(variable, None)
    env.update({
        "PATH": f"{bin_dir}{os.pathsep}{env['PATH']}",
        "MOCK_GIT_BARRIER": str(barrier),
        "MOCK_GIT_FAIL_SOURCE": fail_source,
    })
    completed = subprocess.run(
        ["bash", str(REPO_ROOT / "container/scripts/resolve_source_revisions.sh")],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=20,
    )
    assert len(list(barrier.iterdir())) == 11
    if fail_source:
        assert completed.returncode != 0
        assert completed.stdout == ""
        assert "One or more upstream source revisions could not be resolved" in completed.stderr
    else:
        assert completed.returncode == 0, completed.stderr
        records = dict(line.split("=", 1) for line in completed.stdout.splitlines())
        assert list(records) == list(PROGRAM_SHA_VARS)
        assert all(re.fullmatch("[0-9a-f]{40}", sha) for sha in records.values())


def test_source_stages_have_independent_revision_inputs_and_no_runtime_ancestry():
    dockerfile = (REPO_ROOT / "container/Dockerfile").read_text()
    declarations = list(re.finditer(r"^FROM (\S+) AS (\S+)$", dockerfile, re.MULTILINE))
    stages = {}
    for index, declaration in enumerate(declarations):
        end = declarations[index + 1].start() if index + 1 < len(declarations) else len(dockerfile)
        stages[declaration[2]] = (declaration[1], dockerfile[declaration.end():end])

    def ancestry(name, target="runtime"):
        while name in stages:
            yield name
            name = stages[name][0].replace("${BUILD_TARGET}", target)

    source_stages = {
        name: set(re.findall(r"^ARG ([A-Z0-9_]+_REPO_SHA)=", body, re.MULTILINE))
        for name, (_, body) in stages.items()
        if re.search(r"^ARG [A-Z0-9_]+_REPO_SHA=", body, re.MULTILINE)
    }
    assert len(source_stages) == len(PROGRAM_SHA_VARS)
    assert set.union(*source_stages.values()) == set(PROGRAM_SHA_VARS)
    for name, revisions in source_stages.items():
        assert len(revisions) == 1
        assert set(ancestry(name)) & source_stages.keys() == {name}
    for target in ("runtime", "development"):
        parents = set(ancestry(target, target))
        assert not parents.intersection(source_stages)
        assert not parents.intersection({"dependencies", "source-builder"})
        assert f"{target}-system" in parents
        if target == "runtime":
            assert "system" not in parents
        else:
            assert "system" in parents
        assert f'[[ "${{BUILD_TARGET}}" == "{target}" ]]' in stages[target][1]
    assert dockerfile.count("COPY --from=dependencies /opt/conda /opt/conda") == 1


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
    assert "list_build_inputs.py" in build_hash
    inputs = declared_context_inputs()
    assert "workflow/support/treevis/DESCRIPTION" in inputs
    assert "workflow/support/treevis/NAMESPACE" in inputs
    assert "--stage-native-context" in apptainer_build
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


def test_container_build_paths_share_the_fail_closed_aster_dependency_correction():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (
        REPO_ROOT / "container" / "apptainer_local_build.def.template"
    ).read_text(encoding="utf-8")
    script_name = "ensure_aster_loader_compatibility.sh"

    assert f"COPY container/scripts/{script_name}" in dockerfile
    assert f"bash /opt/pg/scripts/{script_name} /opt/conda" in dockerfile
    assert f"/opt/pg/scripts/{script_name} /opt/conda" in apptainer_template


def test_container_build_paths_include_archive_interoperability_commands():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer_template = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    required_commands = (REPO_ROOT / "container" / "spec" / "required_commands.tsv").read_text(encoding="utf-8")
    arm64_required_commands = (REPO_ROOT / "container" / "spec" / "required_commands.arm64.tsv").read_text(
        encoding="utf-8"
    )

    apt_runtime = (REPO_ROOT / "container/apt/runtime.txt").read_text().splitlines()
    assert "libarchive-tools" in apt_runtime
    assert "install_system_packages.sh runtime" in dockerfile
    assert "install_system_packages.sh development" in apptainer_template
    assert "base\tbsdtar" in required_commands
    assert "base\tbsdtar" in arm64_required_commands
    assert "apt-get install -y --no-install-recommends unzip" in dockerfile
    assert "unzip" in apt_runtime
    assert "base\tunzip" in required_commands
    assert "base\tunzip" in arm64_required_commands
