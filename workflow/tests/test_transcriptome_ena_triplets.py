import gzip
import json
import os
import subprocess
from pathlib import Path

import pytest
from test_transcriptome_public_fastq_fallback import (
    CORE_PATH,
    _run_fallback,
    _run_manifest_validator,
)


@pytest.mark.parametrize("cached", [(), ("_1", "_2"), ("", "_1"), ("", "_1", "_2")])
def test_ena_triplet_download_and_restart_keep_all_reads(monkeypatch, tmp_path, cached):
    metadata = tmp_path / "metadata.tsv"
    output = tmp_path / "getfastq"
    run_dir = output / "SRR123"
    run_dir.mkdir(parents=True)
    # Singleton sorts before the mates lexically; provider order is immaterial.
    urls = {suffix: f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123{suffix}.fastq.gz"
            for suffix in ("", "_2", "_1")}
    metadata.write_text("run\tfastq_ftp\nSRR123\t" + ";".join(urls.values()) + "\n")
    payloads = {
        suffix: gzip.compress((f"@read{suffix or 'orphan'}\nACGT\n+\n!!!!\n" *
                               (2 if suffix else 1)).encode(), mtime=0)
        for suffix in urls
    }
    for suffix in cached:
        (run_dir / f"SRR123{suffix}.amalgkit.fastq.gz").write_bytes(payloads[suffix])
    attempts = _run_fallback(
        monkeypatch, metadata, output,
        {url: payloads[suffix] for suffix, url in urls.items() if suffix not in cached},
    )
    assert set(attempts) == {url for suffix, url in urls.items() if suffix not in cached}
    manifest = json.loads((output / "getfastq_completion.json").read_text())
    assert {entry["path"] for entry in manifest["runs"][0]["files"]} == {
        f"SRR123/SRR123{suffix}.amalgkit.fastq.gz" for suffix in urls
    }
    for suffix in urls:
        assert (run_dir / f"SRR123{suffix}.amalgkit.fastq.gz").read_bytes() == payloads[suffix]
    # Quant still consumes the two mates: do not inflate its read lengths with
    # orphan bases, even though the completion manifest preserves all three.
    assert "SRR123\t2\t16\t16\t16\n" in (run_dir / "getfastq_stats.tsv").read_text()
    result = _run_manifest_validator(metadata, output)
    assert result.returncode == 0, result.stdout + result.stderr
    _run_fallback(monkeypatch, metadata, output, {})
    metadata.write_text("run\nSRR123\n")
    _run_fallback(monkeypatch, metadata, output, {})


def test_ena_triplet_rejects_ambiguous_roles_before_download(monkeypatch, tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "run\tfastq_ftp\nSRR123\t" + ";".join(
            f"https://ftp.sra.ebi.ac.uk/reads/library{i}.fastq.gz" for i in range(3)
        ) + "\n"
    )
    with pytest.raises(SystemExit, match="Ambiguous public FASTQ triplet"):
        _run_fallback(monkeypatch, metadata, tmp_path / "getfastq", {})


def test_ena_triplet_still_rejects_unequal_mate_counts(monkeypatch, tmp_path):
    output = tmp_path / "getfastq"
    run_dir = output / "SRR123"
    run_dir.mkdir(parents=True)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("run\nSRR123\n")
    for suffix, count in (("", 1), ("_1", 2), ("_2", 1)):
        (run_dir / f"SRR123{suffix}.amalgkit.fastq.gz").write_bytes(
            gzip.compress(b"@read\nACGT\n+\n!!!!\n" * count)
        )
    with pytest.raises(SystemExit, match="paired FASTQ record counts differ"):
        _run_fallback(monkeypatch, metadata, output, {})
    assert not (output / "getfastq_completion.json").exists()


def test_short_read_classification_retains_triplet_singleton(tmp_path):
    run_dir = tmp_path / "reads"
    run_dir.mkdir()
    for suffix in ("", "_1", "_2"):
        (run_dir / f"SRR123{suffix}.amalgkit.fastq.gz").touch()
    classification = tmp_path / "classification.tsv"
    classification.write_text("run\tread_class\nSRR123\tshort_read\n")
    result = subprocess.run(
        ["bash", "-c", 'GG_CORE_SOURCE_ONLY=1 source "$1"; '
         'dir_amalgkit_getfastq_sp="$2"; load_classified_getfastq_files "$3"; '
         'printf "%s\\n" "${classified_short_single_fastq_files[@]}" '
         '"${classified_short_left_fastq_files[@]}" "${classified_short_right_fastq_files[@]}"',
         "test", str(CORE_PATH), str(run_dir), str(classification)],
        capture_output=True, text=True, check=True,
    )
    assert set(result.stdout.splitlines()) == {str(path) for path in run_dir.iterdir()}


@pytest.mark.parametrize("method,protocol", [("rnaspades", "same"), ("rnaspades", "mixed"),
                                              ("trinity", "same")])
@pytest.mark.parametrize("sample", [False, True])
@pytest.mark.parametrize("run_count", [1, 10])
def test_assembly_passes_companion_with_both_mates(tmp_path, method, protocol, sample, run_count):
    reads = tmp_path / "reads"
    reads.mkdir()
    for run_index in range(run_count):
        for suffix in ("", "_1", "_2"):
            (reads / f"SRR123{run_index}{suffix}.amalgkit.fastq.gz").write_bytes(
                gzip.compress(b"@read\nACGT\n+\n!!!!\n")
            )
    # Exercise the actual assembly selection/sampling/argument block. Only
    # external read statistics/sampling and the expensive assembler are stubs.
    fakebin = tmp_path / "bin"
    fakebin.mkdir()
    seqkit = fakebin / "seqkit"
    seqkit.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, shutil, sys\n"
        "args = sys.argv[1:]\n"
        "if args[0] == 'stats':\n"
        "    print('file\\tnum_seqs\\tsum_len')\n"
        "    for f in args[2:]: print(f + '\\t1\\t4')\n"
        "elif args[0] == 'sample':\n"
        "    shutil.copyfile(args[-1], args[args.index('--out-file') + 1])\n"
        "    with open('sampled.txt', 'a') as out: out.write(pathlib.Path(args[-1]).name + '\\n')\n"
        "else: raise SystemExit(args)\n"
    )
    seqkit.chmod(0o755)
    text = CORE_PATH.read_text()
    start = text.index('    mapfile -t files_right < <(find "${dir_amalgkit_getfastq_sp}"')
    stop = text.index('      if [[ -d "${dir_tmp}/rnaspades_output" ]]', start)
    block = text[start:stop] + '\n    fi\n'
    script = (
        'set -eo pipefail\nsource "${1%/core/*}/support/gg_util.sh"\n'
        'GG_CORE_SOURCE_ONLY=1 source "$1"\n'
        'dir_amalgkit_getfastq_sp="$2"\ndir_tmp="$PWD"\n'
        f'effective_assembly_method={method}\nprotocol_rna_seq={protocol}\n'
        f'max_assembly_input_fastq_size={6 if sample else 1000}\n'
        'assembly_cpus=1\nassembly_mem_gb=1\nbflyHeapSpaceMax=1\n'
        'Trinity() { printf "%s\\0" "$@" > args.bin; }\n' + block +
        ('printf "%s\\0" "${rnaspades_input_args[@]}" > args.bin\n' if method == "rnaspades" else '')
    )
    result = subprocess.run(
        ["bash", "-c", script, "test", str(CORE_PATH), str(reads)], cwd=tmp_path,
        env={**os.environ, "PATH": str(fakebin) + os.pathsep + os.environ["PATH"]},
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    args = (tmp_path / "args.bin").read_bytes().decode().strip("\0").split("\0")
    selected_runs = min(run_count, 9) if method == "rnaspades" and protocol == "mixed" else run_count
    if method == "rnaspades":
        flags = args[::2]
        for index in range(selected_runs):
            library = index + 1 if protocol == "mixed" else 1
            assert f"--pe{library}-1" in flags
            assert f"--pe{library}-2" in flags
            assert f"--pe{library}-s" in flags
        paths = args[1::2]
    else:
        paths = args[args.index("--left") + 1].split(",") + args[args.index("--right") + 1].split(",")
    names = {Path(path).name for path in paths}
    assert len(paths) == len(names) == selected_runs * 3
    assert names <= {path.name for path in reads.iterdir()}
    for name in names:
        if name.endswith("_1.amalgkit.fastq.gz"):
            assert name.replace("_1.amalgkit", "_2.amalgkit") in names
            assert name.replace("_1.amalgkit", ".amalgkit") in names
    assert (f"Total fastq length is {selected_runs * 12} bp" in result.stdout or
            f"Total paired fastq length is {selected_runs * 12} bp" in result.stdout)
    if sample:
        assert set((tmp_path / "sampled.txt").read_text().splitlines()) == names
