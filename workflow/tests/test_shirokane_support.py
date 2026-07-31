import hashlib
import os
import shlex
import subprocess
from pathlib import Path

import pytest
from shell_static_helpers import (
    WORKFLOW_DIR,
    entrypoint_scheduler_header,
    function_body,
    read_text,
)

SITE_RUNTIME_PATH = WORKFLOW_DIR / "support" / "gg_site_runtime.sh"
GG_UTIL_PATH = WORKFLOW_DIR / "support" / "gg_util.sh"
SUBMIT_HELPER_PATH = WORKFLOW_DIR / "gg_shirokane_submit.sh"
PREPARE_SIF_PATH = WORKFLOW_DIR / "gg_shirokane_prepare_sif.sh"


def run_bash(command: str, cwd: Path):
    return subprocess.run(
        ["bash", "-c", command],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def write_fake_apptainer(fake_bin: Path, arch: str = "x86_64") -> Path:
    fake_apptainer = fake_bin / "apptainer"
    fake_apptainer.write_text(
        "#!/usr/bin/env bash\n"
        'if [[ "${1:-}" == "exec" ]]; then\n'
        f"  printf '%s\\n' {shlex.quote(arch)}\n"
        "fi\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_apptainer.chmod(0o755)
    return fake_apptainer


def test_shirokane_site_profile_is_detected_from_sge_root(tmp_path):
    command = (
        "unset GG_SITE_PROFILE; "
        "export SGE_ROOT=/home/geadmin/N1GE; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        "gg_detect_site_profile"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "shirokane"


def test_shirokane_apptainer_module_loader_uses_site_module_directory(tmp_path):
    modulefiles_dir = tmp_path / "modulefiles"
    fake_bin = tmp_path / "bin"
    modulefiles_dir.mkdir()
    fake_bin.mkdir()
    fake_apptainer = fake_bin / "apptainer"
    fake_apptainer.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake_apptainer.chmod(0o755)

    command = (
        f"fake_bin={shlex.quote(str(fake_bin))}; "
        "module() { "
        'if [[ "$1" == "use" ]]; then return 0; fi; '
        'if [[ "$1" == "load" && "$2" == "apptainer" ]]; then '
        'export PATH="${fake_bin}:${PATH}"; return 0; fi; '
        "return 1; "
        "}; "
        f"export GG_SHIROKANE_MODULEFILES_DIR={shlex.quote(str(modulefiles_dir))}; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        "gg_shirokane_load_apptainer_module; "
        'printf "%s\\n" "$(command -v apptainer)"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(fake_apptainer)


def test_shirokane_submit_helper_defaults_to_ljob(tmp_path):
    entrypoint = WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh"
    completed = subprocess.run(
        [
            "bash",
            str(SUBMIT_HELPER_PATH),
            "--entrypoint",
            str(entrypoint),
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert "qsub" in completed.stdout
    assert "qsub -terse" in completed.stdout
    assert "-l ljob" in completed.stdout


@pytest.mark.parametrize(
    "tasks",
    [
        "0",
        "0-2",
        "2-1",
        "1-3:0",
        "1,,2",
    ],
)
def test_shirokane_submit_helper_rejects_invalid_task_ranges(tmp_path, tasks):
    entrypoint = WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh"
    completed = subprocess.run(
        [
            "bash",
            str(SUBMIT_HELPER_PATH),
            "--entrypoint",
            str(entrypoint),
            "--tasks",
            tasks,
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 2
    assert f"Invalid AGE task range: {tasks}" in completed.stderr


def test_shirokane_submit_helper_accepts_symlinked_repository_path(tmp_path):
    repo_link = tmp_path / "repo-link"
    repo_link.symlink_to(WORKFLOW_DIR.parent, target_is_directory=True)
    helper = repo_link / "workflow" / "gg_shirokane_submit.sh"

    completed = subprocess.run(
        [
            "bash",
            str(helper),
            "--entrypoint",
            "gg_gene_evolution_entrypoint.sh",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert "-l ljob" in completed.stdout


def test_shirokane_submit_helper_applies_resource_overrides(tmp_path):
    entrypoint = WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh"
    completed = subprocess.run(
        [
            "bash",
            str(SUBMIT_HELPER_PATH),
            "--entrypoint",
            str(entrypoint),
            "--tasks",
            "1-8",
            "--slots",
            "4",
            "--mem-per-slot",
            "32G",
            "--resource",
            "lmem",
            "--verify",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert "-verify" in completed.stdout
    assert "-pe def_slot 4" in completed.stdout
    assert "-l s_vmem=32G\\,lmem" in completed.stdout
    assert "-t 1-8" in completed.stdout


def test_shirokane_sif_preparation_job_requests_ljob():
    text = read_text(PREPARE_SIF_PATH)

    assert "#$ -S /bin/bash" in text
    assert "#$ -cwd" in text
    assert "#$ -l s_vmem=8G" in text
    assert "#$ -l ljob" in text
    assert "GG_SHIROKANE_PREBUILT_SIF is required" in text
    assert "GG_SHIROKANE_SIF_TAG is required" in text
    assert "GG_SHIROKANE_SIF_SHA256 is required" in text
    assert "gg_shirokane_load_apptainer_module" in text
    assert 'apptainer exec --cleanenv "${path}" uname -m' in text
    assert 'cp -- "${source_sif}" "${partial_path}"' in text
    assert "install.lock" in text
    assert 'mv -Tf -- "${stable_link_partial}" "${stable_link}"' in text


def test_shirokane_sif_preparation_installs_validated_prebuilt_image(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    write_fake_apptainer(fake_bin)
    source_sif = tmp_path / "prebuilt.sif"
    source_sif.write_bytes(b"validated-sif")
    digest = hashlib.sha256(source_sif.read_bytes()).hexdigest()
    sif_dir = tmp_path / "containers"
    stable_link = tmp_path / "genegalleon.sif"
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{fake_bin}:{env['PATH']}",
            "GG_SHIROKANE_PREBUILT_SIF": str(source_sif),
            "GG_SHIROKANE_SIF_TAG": "sha256-test",
            "GG_SHIROKANE_SIF_SHA256": digest,
            "GG_SHIROKANE_SIF_DIR": str(sif_dir),
            "GG_SHIROKANE_SIF_LINK": str(stable_link),
        }
    )

    completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    installed_sif = sif_dir / "genegalleon-sha256-test.sif"
    assert completed.returncode == 0, completed.stderr
    assert installed_sif.read_bytes() == source_sif.read_bytes()
    assert stable_link.is_symlink()
    assert os.readlink(stable_link) == "containers/genegalleon-sha256-test.sif"
    assert stable_link.resolve() == installed_sif.resolve()
    checksum_path = sif_dir / "genegalleon-sha256-test.sif.sha256"
    assert checksum_path.is_file()
    assert checksum_path.read_text(encoding="utf-8").split()[1] == "genegalleon-sha256-test.sif"

    env["GG_SHIROKANE_SIF_TAG"] = "sha256-next"
    second_completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert second_completed.returncode == 0, second_completed.stderr
    assert os.readlink(stable_link) == "containers/genegalleon-sha256-next.sif"
    assert stable_link.resolve() == (sif_dir / "genegalleon-sha256-next.sif").resolve()


def test_shirokane_sif_preparation_rejects_checksum_mismatch(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    write_fake_apptainer(fake_bin)
    source_sif = tmp_path / "prebuilt.sif"
    source_sif.write_bytes(b"invalid-sif")
    sif_dir = tmp_path / "containers"
    stable_link = tmp_path / "genegalleon.sif"
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{fake_bin}:{env['PATH']}",
            "GG_SHIROKANE_PREBUILT_SIF": str(source_sif),
            "GG_SHIROKANE_SIF_TAG": "sha256-test",
            "GG_SHIROKANE_SIF_SHA256": "0" * 64,
            "GG_SHIROKANE_SIF_DIR": str(sif_dir),
            "GG_SHIROKANE_SIF_LINK": str(stable_link),
        }
    )

    completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "SHA-256 mismatch" in completed.stderr
    assert not (sif_dir / "genegalleon-sha256-test.sif").exists()
    assert not stable_link.exists()


def test_shirokane_sif_preparation_requires_checksum(tmp_path):
    source_sif = tmp_path / "prebuilt.sif"
    source_sif.write_bytes(b"unchecked-sif")
    env = os.environ.copy()
    env.update(
        {
            "GG_SHIROKANE_PREBUILT_SIF": str(source_sif),
            "GG_SHIROKANE_SIF_TAG": "sha256-test",
            "GG_SHIROKANE_SIF_DIR": str(tmp_path / "containers"),
            "GG_SHIROKANE_SIF_LINK": str(tmp_path / "genegalleon.sif"),
        }
    )

    completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "GG_SHIROKANE_SIF_SHA256 is required" in completed.stderr


def test_shirokane_sif_preparation_rejects_non_amd64_image(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    write_fake_apptainer(fake_bin, arch="aarch64")
    source_sif = tmp_path / "prebuilt.sif"
    source_sif.write_bytes(b"arm64-sif")
    digest = hashlib.sha256(source_sif.read_bytes()).hexdigest()
    sif_dir = tmp_path / "containers"
    stable_link = tmp_path / "genegalleon.sif"
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{fake_bin}:{env['PATH']}",
            "GG_SHIROKANE_PREBUILT_SIF": str(source_sif),
            "GG_SHIROKANE_SIF_TAG": "sha256-test",
            "GG_SHIROKANE_SIF_SHA256": digest,
            "GG_SHIROKANE_SIF_DIR": str(sif_dir),
            "GG_SHIROKANE_SIF_LINK": str(stable_link),
        }
    )

    completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "SIF architecture mismatch" in completed.stderr
    assert not (sif_dir / "genegalleon-sha256-test.sif").exists()
    assert not stable_link.exists()


def test_shirokane_sif_preparation_rejects_concurrent_same_tag_install(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    write_fake_apptainer(fake_bin)
    source_sif = tmp_path / "prebuilt.sif"
    source_sif.write_bytes(b"validated-sif")
    digest = hashlib.sha256(source_sif.read_bytes()).hexdigest()
    sif_dir = tmp_path / "containers"
    sif_dir.mkdir()
    lock_dir = sif_dir / ".genegalleon-sha256-test.install.lock"
    lock_dir.mkdir()
    (lock_dir / "owner").write_text("job_id=another-job\n", encoding="utf-8")
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{fake_bin}:{env['PATH']}",
            "GG_SHIROKANE_PREBUILT_SIF": str(source_sif),
            "GG_SHIROKANE_SIF_TAG": "sha256-test",
            "GG_SHIROKANE_SIF_SHA256": digest,
            "GG_SHIROKANE_SIF_DIR": str(sif_dir),
            "GG_SHIROKANE_SIF_LINK": str(tmp_path / "genegalleon.sif"),
        }
    )

    completed = subprocess.run(
        ["bash", str(PREPARE_SIF_PATH)],
        cwd=WORKFLOW_DIR.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "Another SIF installation is using this tag" in completed.stderr
    assert "job_id=another-job" in completed.stderr


def test_uge_headers_default_to_shirokane_ljob():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints

    for entrypoint in entrypoints:
        header = entrypoint_scheduler_header(entrypoint)
        uge_header = header.split("## UGE", 1)[1].split("## PBS", 1)[0]
        assert "#$ -S /bin/bash" in uge_header
        assert "#$ -cwd" in uge_header
        assert "#$ -pe def_slot " in uge_header
        assert "#$ -l s_vmem=" in uge_header
        assert "#$ -l ljob" in uge_header
        assert "#$ -l mem_req=" not in uge_header
        assert "#$ -l epyc" not in uge_header
        assert "#$ -l d_rt=" not in uge_header
        assert "#$ -l s_rt=" not in uge_header


def test_sge_memory_value_conversion_rounds_up_to_whole_gib(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'printf "%s\\n" '
        '"$(gg_sge_memory_value_to_gb 8G)" '
        '"$(gg_sge_memory_value_to_gb 4096M)" '
        '"$(gg_sge_memory_value_to_gb 1.5G)" '
        '"$(gg_sge_memory_value_to_gb 1T)"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == ["8", "4", "2", "1024"]


def test_scheduler_normalization_reads_age_s_vmem(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_qstat = fake_bin / "qstat"
    fake_qstat.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' 'hard resource_list: ljob=TRUE,s_vmem=8G'\n",
        encoding="utf-8",
    )
    fake_qstat.chmod(0o755)

    command = (
        f"export PATH={shlex.quote(str(fake_bin))}:$PATH; "
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "unset GG_TASK_CPUS GG_JOB_ID GG_ARRAY_TASK_ID GG_MEM_PER_CPU_GB GG_MEM_TOTAL_GB; "
        "unset GG_MEM_TOOL_RESERVE_GB GG_MEM_TOOL_GB MEM_PER_SLOT MEM_PER_HOST SGE_TASK_ID; "
        "export SGE_ROOT=/home/geadmin/N1GE JOB_ID=42 NSLOTS=4; "
        "gg_normalize_scheduler_env >/dev/null; "
        'printf "cpus=%s\\nmem_per_cpu=%s\\ntotal=%s\\njob=%s\\n" '
        '"${GG_TASK_CPUS}" "${GG_MEM_PER_CPU_GB}" "${GG_MEM_TOTAL_GB}" "${GG_JOB_ID}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "cpus=4",
        "mem_per_cpu=8",
        "total=32",
        "job=42",
    ]


def test_scheduler_normalization_maps_age_undefined_task_to_one(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_qstat = fake_bin / "qstat"
    fake_qstat.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' 'hard_resource_list: ljob=TRUE,s_vmem=8G'\n",
        encoding="utf-8",
    )
    fake_qstat.chmod(0o755)

    command = (
        f"export PATH={shlex.quote(str(fake_bin))}:$PATH; "
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "unset GG_TASK_CPUS GG_JOB_ID GG_ARRAY_TASK_ID GG_MEM_PER_CPU_GB GG_MEM_TOTAL_GB; "
        "unset GG_MEM_TOOL_RESERVE_GB GG_MEM_TOOL_GB MEM_PER_SLOT MEM_PER_HOST; "
        "export SGE_ROOT=/home/geadmin/N1GE JOB_ID=42 NSLOTS=4 SGE_TASK_ID=undefined; "
        "gg_normalize_scheduler_env >/dev/null; "
        'printf "array=%s\\nsge=%s\\n" "${GG_ARRAY_TASK_ID}" "${SGE_TASK_ID}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == ["array=1", "sge=1"]


def test_shirokane_runtime_prefers_apptainer_over_singularity(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    for runtime in ("apptainer", "singularity"):
        path = fake_bin / runtime
        path.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)

    command = (
        f"export PATH={shlex.quote(str(fake_bin))}:$PATH; "
        "export SGE_ROOT=/home/geadmin/N1GE; "
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "gg_auto_enable_docker_runtime_if_available() { return 1; }; "
        "gg_detect_container_runtime_binary"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "apptainer"


def test_container_runtime_forwards_blas_thread_limits():
    util_text = read_text(GG_UTIL_PATH)
    body = function_body(util_text, "set_singularityenv")

    for prefix in ("SINGULARITYENV", "APPTAINERENV"):
        for name in ("OMP", "OPENBLAS", "MKL", "NUMEXPR"):
            assert f"export {prefix}_{name}_NUM_THREADS=${{GG_TASK_CPUS:-1}}" in body


def test_scheduler_normalization_precedes_duplicate_job_check():
    util_text = read_text(GG_UTIL_PATH)
    body = function_body(util_text, "gg_entrypoint_prepare_container_runtime")

    assert body.index("gg_normalize_scheduler_env") < body.index("exit_if_running_qstat")
    assert body.index("exit_if_running_qstat") < body.index("set_singularity_command")
