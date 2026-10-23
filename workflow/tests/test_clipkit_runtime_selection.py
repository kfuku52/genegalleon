"""Exercise the real native ClipKIT stage with a selected executable."""
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_gene_evolution_core.sh"


def stage_script():
    source = CORE.read_text()
    start = source.index('task="ClipKIT"')
    end = source.index('\ntask="cdskit stats for cleaned alignment"', start)
    return """set -euo pipefail
og_id=OG
dir_output_active="$PWD/result"
gg_workspace_dir="$PWD"
file_og_untrimmed_aln_analysis="$PWD/input.fa"
file_og_clipkit="$PWD/result/alignment.fa.gz"
file_og_clipkit_log="$PWD/result/alignment.log"
genetic_code=1
GG_TASK_CPUS=1
run_clipkit=${TEST_RUN_CLIPKIT:-1}
input_sequence_mode="$TEST_MODE"
clipkit_binary="${TEST_BINARY:-clipkit}"
clipkit_source_manifest="${TEST_MANIFEST:-}"
disable_if_no_input_file() { :; }
gg_artifact_prepare_stage() {
  local variable=$1
  shift 2
  printf '%s\\0' "$@" > "$PWD/provenance.args"
  printf -v "$variable" '%s' 1
}
gg_artifact_record() { :; }
gg_step_start() { :; }
gg_step_skip() { :; }
mv_out() { mv -- "$1" "$2"; }
cp_out() { cp -- "$1" "$2"; }
set_analysis_file() { :; }
""" + source[start:end]


def run_stage(tmp_path, mode="protein", binary=None, manifest=None, enabled=True):
    (tmp_path / "result").mkdir(exist_ok=True)
    (tmp_path / "input.fa").write_text(
        ">first full description\nATGACGACG\n>second\nATGACGACT\n"
    )
    env = dict(os.environ, TEST_MODE=mode, TEST_RUN_CLIPKIT=str(int(enabled)))
    if binary is not None:
        env["TEST_BINARY"] = str(binary)
    if manifest is not None:
        env["TEST_MANIFEST"] = str(manifest)
    return subprocess.run(["bash", "-c", stage_script()], cwd=tmp_path, env=env,
                          text=True, capture_output=True)


@pytest.mark.parametrize("mode", ["protein", "cds"])
def test_selected_clipkit_executable_and_manifest_are_bound(tmp_path, mode):
    native = shutil.which("clipkit")
    assert native, "The GeneGalleon runtime must provide ClipKIT"
    directory = tmp_path / "runtime with spaces"
    directory.mkdir()
    wrapper = directory / "clipkit tool"
    wrapper.write_text(
        "#!/usr/bin/env python\nimport json,os,sys\n"
        "from pathlib import Path\n"
        "Path('native-argv.json').write_text(json.dumps(sys.argv[1:]))\n"
        f"os.execv({native!r}, [{native!r}, *sys.argv[1:]])\n"
    )
    wrapper.chmod(0o700)
    manifest = directory / "source manifest.json"
    manifest.write_text('{"source": "complete-native-fixture"}\n')
    result = run_stage(tmp_path, mode, wrapper, manifest)
    assert result.returncode == 0, result.stdout + result.stderr
    args = json.loads((tmp_path / "native-argv.json").read_text())
    assert args[args.index("--mode") + 1] == "smart-gap"
    assert args[args.index("--sequence_type") + 1] == ("nt" if mode == "cds" else "aa")
    assert ("--codon" in args) == (mode == "cds")
    provenance = (tmp_path / "provenance.args").read_bytes().split(b"\0")
    assert f"clipkit_executable={wrapper}".encode() in provenance
    assert f"clipkit_source_manifest={manifest}".encode() in provenance
    assert (tmp_path / "result/alignment.fa.gz").is_file()
    assert (tmp_path / "result/alignment.log").is_file()


def test_default_clipkit_preserves_existing_provenance_inputs(tmp_path):
    result = run_stage(tmp_path)
    assert result.returncode == 0, result.stderr
    provenance = (tmp_path / "provenance.args").read_bytes()
    assert b"clipkit_executable=" not in provenance
    assert b"clipkit_source_manifest=" not in provenance


@pytest.mark.parametrize("missing", ["executable", "manifest"])
def test_missing_selected_runtime_stops_before_existing_outputs_change(tmp_path, missing):
    output = tmp_path / "result"
    output.mkdir()
    sentinel = output / "alignment.fa.gz"
    sentinel.write_bytes(b"preserved prior artifact")
    native = shutil.which("clipkit")
    assert native
    result = run_stage(tmp_path, binary=tmp_path / "absent" if missing == "executable" else native,
                       manifest=tmp_path / "absent-manifest")
    assert result.returncode != 0
    assert "Missing" in result.stderr
    assert sentinel.read_bytes() == b"preserved prior artifact"
    assert not (tmp_path / "provenance.args").exists()


def test_disabled_clipkit_does_not_require_unused_external_runtime(tmp_path):
    result = run_stage(tmp_path, binary=tmp_path / "absent", enabled=False)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "result/alignment.fa.gz").exists()
