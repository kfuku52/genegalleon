from pathlib import Path
import csv
import os
import sqlite3
import subprocess
import sys
import textwrap


REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_PATH = REPO_ROOT / "workflow" / "core" / "gg_input_generation_core.sh"


def _write_text(path: Path, text: str, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(mode)


def _write_direct_species_fixture(root: Path) -> Path:
    input_dir = root / "raw_valid" / "Direct" / "species_wise_original"
    for species in ("Arabidopsis_thaliana", "Oryza_sativa"):
        species_dir = input_dir / species
        species_dir.mkdir(parents=True, exist_ok=True)
        _write_text(species_dir / f"{species}.cds.fa", ">gene1\nATGAAATTT\n")
        _write_text(
            species_dir / f"{species}.gff",
            "\n".join(
                [
                    "##gff-version 3",
                    "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                    "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                    "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=gene1.t1.CDS1;Parent=gene1.t1",
                    "",
                ]
            ),
        )
        _write_text(species_dir / f"{species}.genome.fa", ">chr1\nATGAAATTT\n")
    return input_dir


def _write_minimal_ete_taxonomy_db(workspace: Path) -> None:
    db_path = workspace / "downloads" / "ete_taxonomy" / "taxa.sqlite"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute("CREATE TABLE species (taxid INTEGER PRIMARY KEY, spname TEXT, common TEXT, rank TEXT, track TEXT)")
        cur.execute("CREATE TABLE synonym (taxid INTEGER, spname TEXT, common TEXT)")
        cur.executemany(
            "INSERT INTO species (taxid, spname, common, rank, track) VALUES (?, ?, ?, ?, ?)",
            [
                (3702, "Arabidopsis thaliana", "", "species", "1,131567,2759,33090,35493,3700,3701,3702"),
                (4530, "Oryza sativa", "", "species", "1,131567,2759,33090,35493,4527,4528,4530"),
            ],
        )
        cur.executemany(
            "INSERT INTO synonym (taxid, spname, common) VALUES (?, ?, ?)",
            [
                (3702, "Arabidopsis_thaliana", ""),
                (4530, "Oryza_sativa", ""),
            ],
        )
        conn.commit()


def _write_runtime_busco_dataset(workspace: Path, lineage: str = "eukaryota_odb12") -> None:
    lineage_dir = workspace / "downloads" / "busco_downloads" / "lineages" / lineage
    lineage_dir.mkdir(parents=True, exist_ok=True)
    _write_text(lineage_dir / "dataset.cfg", "name=test\n")
    _write_text(lineage_dir / ".download.ready", "ready\n")


def _install_fake_toolchain(root: Path) -> Path:
    bin_dir = root / "fake_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    _write_text(
        bin_dir / "conda",
        textwrap.dedent(
            """\
            #!/usr/bin/env bash
            set -euo pipefail
            if [[ $# -ge 2 && "$1" == "shell.bash" && "$2" == "hook" ]]; then
              cat <<'EOF'
            conda() {
              return 0
            }
            EOF
              exit 0
            fi
            exit 0
            """
        ),
        mode=0o755,
    )
    _write_text(
        bin_dir / "seqkit",
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import gzip
            import sys

            args = sys.argv[1:]
            if not args or args[0] != "seq":
                raise SystemExit(f"Unsupported seqkit invocation: {' '.join(sys.argv[1:])}")

            input_path = ""
            output_path = ""
            i = 1
            while i < len(args):
                arg = args[i]
                if arg == "--threads":
                    i += 2
                elif arg in ("--out-file", "-o"):
                    output_path = args[i + 1]
                    i += 2
                elif arg == "-":
                    input_path = arg
                    i += 1
                elif arg.startswith("-"):
                    raise SystemExit(f"Unsupported seqkit flag: {arg}")
                else:
                    input_path = arg
                    i += 1

            if input_path == "-":
                content = sys.stdin.read()
            elif input_path:
                opener = gzip.open if input_path.endswith(".gz") else open
                with opener(input_path, "rt", encoding="utf-8") as src:
                    content = src.read()
            else:
                content = sys.stdin.read()

            if output_path:
                opener = gzip.open if output_path.endswith(".gz") else open
                with opener(output_path, "wt", encoding="utf-8") as dst:
                    dst.write(content)
            else:
                sys.stdout.write(content)
            """
        ),
        mode=0o755,
    )
    _write_text(
        bin_dir / "busco",
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path

            def parse(flag, default=""):
                if flag in sys.argv:
                    idx = sys.argv.index(flag)
                    if idx + 1 < len(sys.argv):
                        return sys.argv[idx + 1]
                return default

            if "--download" in sys.argv:
                lineage = parse("--download")
                target = Path("busco_downloads") / "lineages" / lineage
                target.mkdir(parents=True, exist_ok=True)
                (target / "dataset.cfg").write_text("name=test\\n", encoding="utf-8")
                raise SystemExit(0)

            lineage_path = Path(parse("--lineage_dataset"))
            lineage = lineage_path.name if str(lineage_path) else "unknown_odb12"
            outdir = Path(parse("--out", "busco_tmp"))
            run_dir = outdir / f"run_{lineage}"
            run_dir.mkdir(parents=True, exist_ok=True)
            infile = Path(parse("--in"))
            header = "unknown_gene"
            if infile.exists():
                with infile.open("rt", encoding="utf-8") as handle:
                    for line in handle:
                        if line.startswith(">"):
                            header = line[1:].strip().split()[0]
                            break
            (run_dir / "full_table.tsv").write_text(
                "# BUSCO version is: fixture\\n"
                f"BUSCO1\\tComplete\\t{header}\\t1\\t100\\thttps://example.org/BUSCO1\\tfixture description\\n",
                encoding="utf-8",
            )
            (run_dir / "short_summary.txt").write_text(
                "C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1\\n",
                encoding="utf-8",
            )
            raise SystemExit(0)
            """
        ),
        mode=0o755,
    )
    _write_text(
        bin_dir / "Rscript",
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path

            def parse(flag, default=""):
                for arg in sys.argv[1:]:
                    if arg.startswith(flag + "="):
                        return arg.split("=", 1)[1]
                return default

            outdir = Path.cwd()
            busco_dir = Path(parse("--dir_species_cds_busco"))
            species = []
            if busco_dir.exists():
                for path in sorted(busco_dir.glob("*.busco.full.tsv")):
                    species.append(path.name.replace(".busco.full.tsv", ""))

            (outdir / "annotation_summary.tsv").write_text(
                "Species\\tbusco_cds_summary\\n"
                + "\\n".join(f"{sp}\\tC:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1" for sp in species)
                + ("\\n" if species else ""),
                encoding="utf-8",
            )
            (outdir / "busco_cds.pdf").write_text("fake pdf\\n", encoding="utf-8")
            (outdir / "busco_cds.svg").write_text("<svg/>\\n", encoding="utf-8")
            raise SystemExit(0)
            """
        ),
        mode=0o755,
    )
    return bin_dir


def _core_env(workspace: Path, input_dir: Path, fake_bin: Path, mode: str, task_id: int | None = None) -> dict[str, str]:
    env = {
        "HOME": os.environ["HOME"],
        "PATH": os.pathsep.join([str(fake_bin), str(Path(sys.executable).parent), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]),
        "TMPDIR": str(workspace / "tmp_runtime"),
        "GG_TASK_CPUS": "1",
        "gg_workspace_dir": str(workspace),
        "provider": "direct",
        "input_dir": str(input_dir),
        "input_generation_mode": mode,
        "run_format_inputs": "1",
        "run_validate_inputs": "1",
        "run_species_busco": "1",
        "run_multispecies_summary": "1",
        "run_generate_species_trait": "0",
        "busco_lineage": "eukaryota_odb12",
        "trait_profile": "none",
        "strict": "0",
        "overwrite": "1",
        "download_only": "0",
        "dry_run": "0",
        "download_timeout": "120",
        "download_manifest": "",
        "download_dir": "",
        "summary_output": "",
        "auth_bearer_token_env": "",
        "http_header": "",
        "species_cds_dir": "",
        "species_busco_full_dir": "",
        "species_busco_short_dir": "",
        "species_gff_dir": "",
        "species_genome_dir": "",
        "species_summary_output": "",
        "resolved_manifest_output": "",
        "species_trait_output": "",
        "task_plan_output": "",
        "trait_plan": "",
        "trait_database_sources": "",
        "trait_download_dir": "",
        "trait_download_timeout": "120",
        "trait_species_source": "download_manifest",
        "trait_databases": "auto",
    }
    if task_id is not None:
        env["GG_ARRAY_TASK_ID"] = str(task_id)
    return env


def _run_core(workspace: Path, input_dir: Path, fake_bin: Path, mode: str, task_id: int | None = None) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        ["bash", str(CORE_PATH)],
        cwd=REPO_ROOT,
        env=_core_env(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode=mode, task_id=task_id),
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    return completed


def _run_core_async(workspace: Path, input_dir: Path, fake_bin: Path, mode: str, task_id: int) -> subprocess.Popen[str]:
    return subprocess.Popen(
        ["bash", str(CORE_PATH)],
        cwd=REPO_ROOT,
        env=_core_env(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode=mode, task_id=task_id),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_expected_outputs(output_root: Path, expected_last_mode: str) -> None:
    summary_dir = output_root / "annotation_summary"
    summary = summary_dir / "annotation_summary.tsv"
    full_dir = output_root / "species_cds_busco_full"
    short_dir = output_root / "species_cds_busco_short"
    runs_path = output_root / "gg_input_generation_runs.tsv"

    assert (output_root / "species_cds" / "Arabidopsis_thaliana_cds.fa.gz").exists()
    assert (output_root / "species_cds" / "Oryza_sativa_cds.fa.gz").exists()
    assert (output_root / "species_gff" / "Arabidopsis_thaliana_gff.gff.gz").exists()
    assert (output_root / "species_gff" / "Oryza_sativa_gff.gff.gz").exists()
    assert (output_root / "species_genome" / "Arabidopsis_thaliana_genome.fa.gz").exists()
    assert (output_root / "species_genome" / "Oryza_sativa_genome.fa.gz").exists()
    assert (full_dir / "Arabidopsis_thaliana.busco.full.tsv").exists()
    assert (full_dir / "Oryza_sativa.busco.full.tsv").exists()
    assert (short_dir / "Arabidopsis_thaliana.busco.short.txt").exists()
    assert (short_dir / "Oryza_sativa.busco.short.txt").exists()

    arabidopsis_busco = (full_dir / "Arabidopsis_thaliana.busco.full.tsv").read_text(encoding="utf-8")
    oryza_busco = (full_dir / "Oryza_sativa.busco.full.tsv").read_text(encoding="utf-8")
    assert "Arabidopsis_thaliana_gene1" in arabidopsis_busco
    assert "Oryza_sativa_gene1" in oryza_busco

    assert (summary_dir / "busco_cds.pdf").exists()
    assert (summary_dir / "busco_cds.svg").exists()
    summary_text = summary.read_text(encoding="utf-8")
    assert "Arabidopsis_thaliana" in summary_text
    assert "Oryza_sativa" in summary_text

    runs = _read_tsv_rows(runs_path)
    assert runs[-1]["input_generation_mode"] == expected_last_mode
    assert runs[-1]["stage_species_busco_status"] == "ok"
    assert runs[-1]["stage_multispecies_summary_status"] == "ok"
    assert runs[-1]["num_species_busco_full"] == "2"
    assert runs[-1]["num_species_busco_short"] == "2"


def test_gg_input_generation_single_mode_end_to_end(tmp_path: Path):
    input_dir = _write_direct_species_fixture(tmp_path)
    workspace = tmp_path / "single_workspace"
    _write_minimal_ete_taxonomy_db(workspace)
    _write_runtime_busco_dataset(workspace)
    fake_bin = _install_fake_toolchain(tmp_path)

    _run_core(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode="single")

    _assert_expected_outputs(workspace / "output" / "input_generation", expected_last_mode="single")


def test_gg_input_generation_array_mode_end_to_end_with_parallel_workers(tmp_path: Path):
    input_dir = _write_direct_species_fixture(tmp_path)
    workspace = tmp_path / "array_workspace"
    _write_minimal_ete_taxonomy_db(workspace)
    fake_bin = _install_fake_toolchain(tmp_path)

    _run_core(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode="array_prepare")
    _write_runtime_busco_dataset(workspace)

    worker1 = _run_core_async(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode="array_worker", task_id=1)
    worker2 = _run_core_async(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode="array_worker", task_id=2)
    stdout1, stderr1 = worker1.communicate(timeout=180)
    stdout2, stderr2 = worker2.communicate(timeout=180)
    assert worker1.returncode == 0, stdout1 + "\n" + stderr1
    assert worker2.returncode == 0, stdout2 + "\n" + stderr2

    _run_core(workspace=workspace, input_dir=input_dir, fake_bin=fake_bin, mode="array_finalize")

    _assert_expected_outputs(workspace / "output" / "input_generation", expected_last_mode="array_finalize")
