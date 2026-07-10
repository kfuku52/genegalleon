import os
import re
import shlex
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Optional

import pandas
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
CORE_PATH = WORKFLOW_DIR / "core" / "gg_genome_evolution_core.sh"
ENTRYPOINT_PATH = WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh"
BASH_FOR_TESTS = shutil.which("bash") or "/bin/bash"
SYSTEM_BASH_MAJOR = int(
    subprocess.run(
        [BASH_FOR_TESTS, "-c", "printf '%s' \"${BASH_VERSINFO[0]}\""],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    or "0"
)


def _write_executable(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _prepare_stub_binaries(tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    capture_dir = tmp_path / "capture"
    capture_dir.mkdir(parents=True, exist_ok=True)

    _write_executable(
        bin_dir / "parallel",
        r"""#!/usr/bin/env bash
set -euo pipefail
func=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)
      shift 2
      ;;
    :::)
      shift
      break
      ;;
    *)
      if [[ -z "${func}" ]]; then
        func="$1"
      fi
      shift
      ;;
  esac
done
for item in "$@"; do
  bash -c "${func} \"\$1\"" _ "${item}"
done
""",
    )

    _write_executable(
        bin_dir / "conda",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "shell.bash" && "${2:-}" == "hook" ]]; then
  cat <<'EOF'
conda() {
  return 0
}
EOF
  exit 0
fi
exit 0
""",
    )

    _write_executable(
        bin_dir / "seqkit",
        r"""#!/usr/bin/env python3
import re
import sys


STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_fasta(lines):
    records = []
    header = None
    seq_chunks = []
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_chunks)))
            header = line[1:]
            seq_chunks = []
            continue
        seq_chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    return records


def load_records(path):
    if path:
        with open(path, "r", encoding="utf-8") as handle:
            return parse_fasta(handle)
    return parse_fasta(sys.stdin)


def write_records(records, handle):
    for header, sequence in records:
        handle.write(f">{header}\n{sequence}\n")


def apply_translation(records, transl_table):
    code = dict(STANDARD_CODE)
    if transl_table == "6":
        code["TAA"] = "Q"
        code["TAG"] = "Q"
    translated = []
    for header, sequence in records:
        sequence = sequence.upper()
        residues = []
        for idx in range(0, len(sequence) - len(sequence) % 3, 3):
            codon = sequence[idx:idx + 3]
            residues.append(code.get(codon, "X"))
        translated.append((header, "".join(residues)))
    return translated


def main():
    if len(sys.argv) < 2:
        raise SystemExit("seqkit stub requires a subcommand")
    subcommand = sys.argv[1]
    args = sys.argv[2:]

    if subcommand == "seq":
        input_path = None
        out_path = None
        remove_gaps = False
        name_only = False
        idx = 0
        while idx < len(args):
            arg = args[idx]
            if arg == "--threads":
                idx += 2
            elif arg == "--out-file":
                out_path = args[idx + 1]
                idx += 2
            elif arg == "--remove-gaps":
                remove_gaps = True
                idx += 1
            elif arg == "--name":
                name_only = True
                idx += 1
            elif arg.startswith("-"):
                raise SystemExit(f"Unsupported seqkit seq option: {arg}")
            else:
                input_path = arg
                idx += 1
        records = load_records(input_path)
        if remove_gaps:
            records = [(header, sequence.replace("-", "")) for header, sequence in records]
        output_handle = open(out_path, "w", encoding="utf-8") if out_path else sys.stdout
        try:
            if name_only:
                for header, _sequence in records:
                    output_handle.write(f"{header.split()[0]}\n")
            else:
                write_records(records, output_handle)
        finally:
            if out_path:
                output_handle.close()
        return

    if subcommand == "replace":
        input_path = None
        pattern = None
        replacement = ""
        by_seq = False
        ignore_case = False
        idx = 0
        while idx < len(args):
            arg = args[idx]
            if arg == "--threads":
                idx += 2
            elif arg == "--pattern":
                pattern = args[idx + 1]
                idx += 2
            elif arg == "--replacement":
                replacement = args[idx + 1]
                idx += 2
            elif arg == "--by-seq":
                by_seq = True
                idx += 1
            elif arg == "--ignore-case":
                ignore_case = True
                idx += 1
            elif arg.startswith("-"):
                raise SystemExit(f"Unsupported seqkit replace option: {arg}")
            else:
                input_path = arg
                idx += 1
        flags = re.IGNORECASE if ignore_case else 0
        records = load_records(input_path)
        updated = []
        for header, sequence in records:
            if by_seq:
                sequence = re.sub(pattern, replacement, sequence, flags=flags)
            else:
                header = re.sub(pattern, replacement, header, flags=flags)
            updated.append((header, sequence))
        write_records(updated, sys.stdout)
        return

    if subcommand == "translate":
        input_path = None
        transl_table = "1"
        idx = 0
        while idx < len(args):
            arg = args[idx]
            if arg == "--threads":
                idx += 2
            elif arg == "--transl-table":
                transl_table = args[idx + 1]
                idx += 2
            elif arg == "--allow-unknown-codon":
                idx += 1
            elif arg.startswith("-"):
                raise SystemExit(f"Unsupported seqkit translate option: {arg}")
            else:
                input_path = arg
                idx += 1
        write_records(apply_translation(load_records(input_path), transl_table), sys.stdout)
        return

    raise SystemExit(f"Unsupported seqkit subcommand: {subcommand}")


if __name__ == "__main__":
    main()
""",
    )

    _write_executable(
        bin_dir / "cdskit",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" != "pad" ]]; then
  echo "Unsupported cdskit subcommand: ${1:-}" >&2
  exit 1
fi
cat
""",
    )

    _write_executable(
        bin_dir / "ulimit",
        """#!/usr/bin/env bash
exit 0
""",
    )

    _write_executable(
        bin_dir / "orthofinder",
        f"""#!/usr/bin/env bash
set -euo pipefail
if [[ "${{1:-}}" == "-v" ]]; then
  echo "OrthoFinder:v2.5.5"
  exit 0
fi
capture_dir={shlex.quote(str(capture_dir))}
printf '%s\n' "$*" >> "${{capture_dir}}/orthofinder_args.txt"
input_dir=""
output_dir=""
run_name=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f)
      input_dir="$2"
      shift 2
      ;;
    -o)
      output_dir="$2"
      shift 2
      ;;
    -n)
      run_name="$2"
      shift 2
      ;;
    --assign)
      input_dir="$2"
      shift 2
      ;;
    --core)
      output_dir="$(dirname "$2")"
      shift 2
      ;;
    -s)
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done
if [[ -z "${{input_dir}}" || -z "${{output_dir}}" || -z "${{run_name}}" ]]; then
  echo "orthofinder stub requires input_dir, output_dir, and run_name" >&2
  exit 1
fi
results_dir="${{output_dir}}/Results_${{run_name}}"
mkdir -p "${{results_dir}}/Phylogenetic_Hierarchical_Orthogroups"
input_capture="${{capture_dir}}/input_files_${{run_name}}.txt"
proteins_capture="${{capture_dir}}/proteins_${{run_name}}.fasta"
find "${{input_dir}}" -maxdepth 1 -type f ! -name '.*' | sort > "${{input_capture}}"
cat "${{input_capture}}" > "${{capture_dir}}/input_files.txt"
>"${{proteins_capture}}"
>"${{capture_dir}}/proteins.fasta"
while IFS= read -r fasta; do
  cat "${{fasta}}" >> "${{proteins_capture}}"
  cat "${{fasta}}" >> "${{capture_dir}}/proteins.fasta"
done < "${{input_capture}}"
python - "${{input_capture}}" "${{results_dir}}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" <<'PY'
import pathlib
import sys

input_files = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
outfile = pathlib.Path(sys.argv[2])
species = []
for path in input_files:
    name = pathlib.Path(path).name
    if "_" in name:
        species.append(name.split("_", 2)[0] + "_" + name.split("_", 2)[1])
    else:
        species.append(pathlib.Path(name).stem)
with outfile.open("w", encoding="utf-8", newline="") as handle:
    cols = ["HOG", "OG", "Gene Tree Parent Clade", *species]
    handle.write("\\t".join(cols) + "\\n")
    genes = [f"{{sp}}_gene1" for sp in species]
    handle.write("\\t".join(["N0.HOG0000001", "OG0000001", "n0", *genes]) + "\\n")
PY
""",
    )

    _write_executable(
        bin_dir / "nwkit",
        f"""#!/usr/bin/env python3
import csv
from functools import cmp_to_key
from pathlib import Path
import sys

capture_dir = {str(capture_dir)!r}


def parse_specs(value):
    return [part for part in value.split(",") if part]


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def row_passes(row, spec):
    column, op, value = spec.split(":", 2)
    left = row.get(column, "")
    left_num = as_float(left)
    right_num = as_float(value)
    if op in {{"ge", "gt", "le", "lt"}}:
        if left_num is None or right_num is None:
            return False
        if op == "ge":
            return left_num >= right_num
        if op == "gt":
            return left_num > right_num
        if op == "le":
            return left_num <= right_num
        if op == "lt":
            return left_num < right_num
    if op == "eq":
        return left == value
    if op == "ne":
        return left != value
    return False


def compare_rows(rank_specs):
    parsed = [spec.split(":", 1) for spec in rank_specs]

    def compare(left, right):
        for column, direction in parsed:
            left_value = left.get(column, "")
            right_value = right.get(column, "")
            left_num = as_float(left_value)
            right_num = as_float(right_value)
            if left_num is not None and right_num is not None:
                if left_num == right_num:
                    continue
                result = -1 if left_num < right_num else 1
            else:
                if left_value == right_value:
                    continue
                result = -1 if left_value < right_value else 1
            return result if direction == "asc" else -result
        return -1 if left.get("leaf_name", "") < right.get("leaf_name", "") else 1

    return compare


def sample(args):
    Path(capture_dir, "nwkit_args.txt").write_text(" ".join(sys.argv[1:]) + "\\n", encoding="utf-8")
    trait = ""
    output_table = ""
    outfile = ""
    n = 0
    filters = []
    ranks = []
    idx = 0
    while idx < len(args):
        arg = args[idx]
        if arg == "--trait":
            trait = args[idx + 1]
            idx += 2
        elif arg == "--output-table":
            output_table = args[idx + 1]
            idx += 2
        elif arg == "--outfile":
            outfile = args[idx + 1]
            idx += 2
        elif arg == "--n":
            n = int(args[idx + 1])
            idx += 2
        elif arg == "--filter":
            filters.append(args[idx + 1])
            idx += 2
        elif arg == "--rank":
            ranks.append(args[idx + 1])
            idx += 2
        elif arg in {{"--infile", "--method", "--allow-fewer"}}:
            idx += 2
        else:
            idx += 1
    with open(trait, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\\t"))
    selected = [row for row in rows if all(row_passes(row, spec) for spec in filters)]
    if ranks:
        selected = sorted(selected, key=cmp_to_key(compare_rows(ranks)))
    selected = selected[:n]
    with open(output_table, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\\t", lineterminator="\\n")
        writer.writeheader()
        writer.writerows(selected)
    leaves = ",".join(row["leaf_name"] + ":0.1" for row in selected)
    Path(outfile).write_text("(" + leaves + ");\\n", encoding="utf-8")


def main():
    if len(sys.argv) < 2:
        raise SystemExit("nwkit stub requires a subcommand")
    if sys.argv[1] == "sample":
        sample(sys.argv[2:])
        return
    raise SystemExit("Unsupported nwkit subcommand: " + sys.argv[1])


if __name__ == "__main__":
    main()
""",
    )

    _write_executable(
        bin_dir / "omamer",
        f"""#!/usr/bin/env bash
set -euo pipefail
capture_dir={shlex.quote(str(capture_dir))}
subcommand="${{1:-}}"
if [[ "${{subcommand}}" != "search" ]]; then
  echo "Unsupported omamer subcommand: ${{subcommand}}" >&2
  exit 1
fi
shift
db=""
query=""
outfile=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)
      db="$2"
      shift 2
      ;;
    --query)
      query="$2"
      shift 2
      ;;
    --out)
      outfile="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done
if [[ -z "${{db}}" || -z "${{query}}" || -z "${{outfile}}" ]]; then
  echo "omamer search stub requires --db, --query, and --out" >&2
  exit 1
fi
if [[ ! -s "${{db}}" ]]; then
  echo "OMAmer DB not found: ${{db}}" >&2
  exit 1
fi
mkdir -p "$(dirname "${{outfile}}")"
printf '%s\\n' "${{db}}" >> "${{capture_dir}}/omamer_db_paths.txt"
printf '%s\\n' "${{query}}" >> "${{capture_dir}}/omamer_queries.txt"
printf '# query\\thog\\tscore\\n%s\\tHOG:0000001\\t100\\n' "$(basename "${{query}}")" > "${{outfile}}"
""",
    )

    _write_executable(
        bin_dir / "omark",
        """#!/usr/bin/env python3
import pathlib
import sys


def main():
    args = sys.argv[1:]
    omamer_file = ""
    db = ""
    outdir = ""
    idx = 0
    while idx < len(args):
        arg = args[idx]
        if arg == "-f":
            omamer_file = args[idx + 1]
            idx += 2
        elif arg == "-d":
            db = args[idx + 1]
            idx += 2
        elif arg == "-o":
            outdir = args[idx + 1]
            idx += 2
        else:
            idx += 1
    if not omamer_file or not db or not outdir:
        raise SystemExit("omark stub requires -f, -d, and -o")
    if not pathlib.Path(db).exists():
        raise SystemExit(f"OMArk DB not found: {db}")

    outdir_path = pathlib.Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    base = pathlib.Path(omamer_file).stem
    (outdir_path / f"{base}.sum").write_text(
        "#The selected clade was Viridiplantae\\n"
        "#Number of conserved HOGs is: 100\\n"
        "#Results on conserved HOGs is:\\n"
        "#S:Single:S, D:Duplicated[U:Unexpected,E:Expected],M:Missing\\n"
        "S:95,D:3[U:1,E:2],M:2\\n"
        "S:95.00%,D:3.00%[U:1.00%,E:2.00%],M:2.00%\\n"
        "#On the whole proteome, there are 1000 proteins\\n"
        "#Of which:\\n"
        "#A:Consistent (taxonomically)[P:Partial hits,F:Fragmented], I: Inconsistent (taxonomically)[P:Partial hits,F:Fragmented], C: Likely Contamination[P:Partial hits,F:Fragmented], U: Unknown \\n"
        "A:900[P:10,F:5],I:50[P:3,F:2],C:20[P:1,F:1],U:30\\n"
        "A:90.00%[P:1.00%,F:0.50%],I:5.00%[P:0.30%,F:0.20%],C:2.00%[P:0.10%,F:0.10%],U:3.00%\\n"
        "#From HOG placement, the detected species are:\\n"
        "#Clade\\tNCBI taxid\\tNumber of associated proteins\\tPercentage of proteome's total\\n"
        "Arabidopsis thaliana\\t3702\\t800\\t80.0%\\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
""",
    )

    return bin_dir


def _load_entrypoint_defaults() -> dict[str, str]:
    text = ENTRYPOINT_PATH.read_text(encoding="utf-8")
    in_block = False
    defaults: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("### Start: Modify this block"):
            in_block = True
            continue
        if line.startswith("### End: Modify this block"):
            break
        if not in_block or not line or line.startswith("#"):
            continue
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*?)(?:\s+#.*)?$", line)
        if not match:
            continue
        key = match.group(1)
        value = match.group(2).strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        defaults[key] = value
    return defaults


def _run_core(
    tmp_path: Path,
    extra_env: Optional[dict[str, str]] = None,
    *,
    missing_flags: tuple[str, ...] = (),
) -> subprocess.CompletedProcess[str]:
    workspace = tmp_path / "workspace"
    workspace.mkdir(parents=True, exist_ok=True)
    bin_dir = _prepare_stub_binaries(tmp_path)
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("CONDA") and not key.startswith("_CE_")
    }
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["gg_workspace_dir"] = str(workspace)
    env["GG_TASK_CPUS"] = "1"
    env["GG_MEM_PER_CPU_GB"] = "1"
    env["GG_MEM_TOTAL_GB"] = "1"
    env["GG_JOB_ID"] = "1"
    env["GG_ARRAY_TASK_ID"] = "1"
    env.update(_load_entrypoint_defaults())
    env.update(
        {
            "input_sequence_mode": "protein",
            "run_species_busco": "0",
            "run_build_species_busco_summary": "0",
            "run_species_omark": "0",
            "run_build_species_omark_summary": "1",
            "run_extract_species_tree_fasta": "0",
            "run_individual_mafft": "0",
            "run_individual_trimal": "0",
            "run_concat_alignment": "0",
            "run_concat_iqtree_protein": "0",
            "run_concat_iqtree_dna": "0",
            "run_individual_iqtree_pep": "0",
            "run_astral_pep": "0",
            "run_individual_iqtree_dna": "0",
            "run_astral_dna": "0",
            "run_plot_species_trees": "0",
            "run_constrained_tree": "0",
            "run_plot_constrained_tree": "0",
            "run_mcmctree1": "0",
            "run_mcmctree2": "0",
            "run_convert_tree_format": "0",
            "run_plot_mcmctreer": "0",
            "run_orthofinder": "1",
            "run_og_selection": "0",
            "run_orthogroup_method_comparison": "0",
            "run_single_copy_ortholog_decay_plot": "0",
            "run_busco_dupaware_extract_fasta": "0",
            "run_busco_dupaware_mafft": "0",
            "run_busco_dupaware_trimal": "0",
            "run_busco_dupaware_iqtree_dna": "0",
            "run_busco_dupaware_iqtree_pep": "0",
            "run_busco_dupaware_notung_root_dna": "0",
            "run_busco_dupaware_notung_root_pep": "0",
            "run_busco_dupaware_root_dna": "0",
            "run_busco_dupaware_root_pep": "0",
            "run_busco_dupaware_grampa_dna": "0",
            "run_busco_dupaware_grampa_pep": "0",
            "run_orthogroup_grampa": "0",
            "run_cafe": "0",
            "run_go_enrichment": "0",
            "change_direction_go": "increase",
            "target_branch_go": "",
        }
    )
    for flag in missing_flags:
        env.pop(flag, None)
    if extra_env:
        env.update(extra_env)

    command = (
        "set -euo pipefail; "
        "mapfile() { "
        "  if [[ ${1:-} == '-t' ]]; then shift; fi; "
        "  local target_array=${1:?}; "
        "  local line=''; "
        "  eval \"${target_array}=()\"; "
        "  while IFS= read -r line; do "
        "    eval \"${target_array}+=(\\\"\\${line}\\\")\"; "
        "  done; "
        "}; "
        "enable -n ulimit; "
        f"source {shlex.quote(str(CORE_PATH))}"
    )
    return subprocess.run(
        [BASH_FOR_TESTS, "-c", command],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_protein_mode_prefers_species_protein_inputs(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_cds_dir = workspace / "input" / "species_cds"
    species_code_dir = workspace / "input" / "species_genetic_code"
    species_protein_dir.mkdir(parents=True)
    species_cds_dir.mkdir(parents=True)
    species_code_dir.mkdir(parents=True)

    (species_protein_dir / "Tetrahymena_thermophila_pep.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_cds_dir / "Tetrahymena_thermophila_cds.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_code_dir / "species_genetic_code.tsv").write_text(
        "species\tgenetic_code\nTetrahymena_thermophila\t6\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "species_genetic_code.tsv is ignored because species_protein inputs are provided" in completed.stdout
    assert "Tetrahymena_thermophila.fa.gz" in completed.stdout
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_thermophila_gene1" in proteins
    assert "MPEP" in proteins
    assert "MQ" not in proteins
    assert "M*" not in proteins


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_protein_mode_translates_species_cds_with_species_specific_codes(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    species_code_dir = workspace / "input" / "species_genetic_code"
    species_cds_dir.mkdir(parents=True)
    species_code_dir.mkdir(parents=True)

    (species_cds_dir / "Tetrahymena_thermophila_cds.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_cds_dir / "Arabidopsis_thaliana_cds.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_code_dir / "species_genetic_code.tsv").write_text(
        "species\tgenetic_code\nTetrahymena_thermophila\t6\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Tetrahymena_thermophila.fa.gz" in completed.stdout
    assert "Arabidopsis_thaliana.fa.gz" in completed.stdout
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_thermophila_gene1" in proteins
    assert ">Arabidopsis_thaliana_gene1" in proteins
    assert "MQ" in proteins
    assert "M*" in proteins


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_protein_mode_preserves_qualified_species_labels_in_genetic_code_lookup(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    species_code_dir = workspace / "input" / "species_genetic_code"
    species_cds_dir.mkdir(parents=True)
    species_code_dir.mkdir(parents=True)

    (species_cds_dir / "Tetrahymena_cf_thermophila_cds.fa").write_text(
        ">Tetrahymena_cf_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_code_dir / "species_genetic_code.tsv").write_text(
        "species\tgenetic_code\nTetrahymena_cf_thermophila\t6\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Tetrahymena_cf_thermophila.fa.gz" in completed.stdout
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_cf_thermophila_gene1" in proteins
    assert "MQ" in proteins


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_core_defaults_shared_protein_flags_for_legacy_launchers(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    species_cds_dir.mkdir(parents=True)

    (species_cds_dir / "Tetrahymena_thermophila_cds.fa").write_text(
        ">Tetrahymena_thermophila_gene1\nATGTAA\n",
        encoding="utf-8",
    )

    completed = _run_core(
        tmp_path,
        missing_flags=("run_cds_translation", "run_species_omark", "run_build_species_omark_summary"),
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Tetrahymena_thermophila.fa.gz" in completed.stdout
    proteins = (tmp_path / "capture" / "proteins.fasta").read_text(encoding="utf-8")
    assert ">Tetrahymena_thermophila_gene1" in proteins
    assert "M*" in proteins
    assert "unbound variable" not in completed.stderr


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_accepts_legacy_species_tree_rooting_outgroup_lists(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_protein_dir.mkdir(parents=True)

    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_protein_dir / "Oryza_sativa_pep.fa").write_text(
        ">Oryza_sativa_gene1\nMPEP\n",
        encoding="utf-8",
    )

    completed = _run_core(
        tmp_path,
        {"species_tree_rooting": "Arabidopsis_thaliana,Oryza_sativa"},
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (
        "species_tree_rooting=Arabidopsis_thaliana,Oryza_sativa uses legacy outgroup-label syntax; "
        "interpreting it as outgroup,Arabidopsis_thaliana,Oryza_sativa."
    ) in completed.stdout
    assert "Resolved species_tree_rooting method: outgroup" in completed.stdout
    assert "Resolved species_tree_rooting value: Arabidopsis_thaliana,Oryza_sativa" in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_refuses_orthofinder_when_requested_species_tree_is_missing(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_protein_dir.mkdir(parents=True)
    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEP\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path, {"run_astral_pep": "1"})

    assert completed.returncode != 0
    assert "Refusing to run OrthoFinder without a species tree." in completed.stdout
    assert "Species-tree generation was requested, but no summary tree is available." in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_refuses_orthofinder_when_species_tree_lacks_input_species(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_tree_summary_dir = workspace / "output" / "species_tree" / "species_tree_summary"
    species_protein_dir.mkdir(parents=True)
    species_tree_summary_dir.mkdir(parents=True)
    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_protein_dir / "Oryza_sativa_pep.fa").write_text(
        ">Oryza_sativa_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_tree_summary_dir / "undated_species_tree.nwk").write_text(
        "(Arabidopsis_thaliana:0.1);\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path)

    assert completed.returncode != 0
    assert "Species tree is missing 1 species: Oryza_sativa" in completed.stdout
    assert "Refusing to run OrthoFinder without species tree constraints" in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_keeps_species_tree_outputs_when_generation_disabled_and_signature_changes(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_tree_dir = workspace / "output" / "species_tree"
    species_tree_summary_dir = species_tree_dir / "species_tree_summary"
    species_protein_dir.mkdir(parents=True)
    species_tree_summary_dir.mkdir(parents=True)
    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEP\n",
        encoding="utf-8",
    )
    existing_tree = species_tree_summary_dir / "undated_species_tree.nwk"
    existing_tree.write_text("(Arabidopsis_thaliana:0.1);\n", encoding="utf-8")
    stamp = species_tree_dir / ".shared_protein_input_signature"
    stamp.write_text("old-signature\n", encoding="utf-8")

    completed = _run_core(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert existing_tree.read_text(encoding="utf-8") == "(Arabidopsis_thaliana:0.1);\n"
    assert stamp.read_text(encoding="utf-8") == "old-signature\n"
    assert "Keeping existing species_tree outputs for reuse" in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_two_round_orthofinder_can_run_without_species_tree_when_generation_disabled(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_protein_dir.mkdir(parents=True)
    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEP\n",
        encoding="utf-8",
    )
    (species_protein_dir / "Oryza_sativa_pep.fa").write_text(
        ">Oryza_sativa_gene1\nMPEP\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path, {"max_orthofinder_core_species": "1"})

    assert completed.returncode == 0, completed.stdout + completed.stderr
    args = (tmp_path / "capture" / "orthofinder_args.txt").read_text(encoding="utf-8")
    assert "species_tree_core.nwk" not in args
    assert "No species tree summary was found. Species-tree generation flags are disabled" in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_two_round_orthofinder_uses_nwkit_sample_default_filters(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_tree_summary_dir = workspace / "output" / "species_tree" / "species_tree_summary"
    busco_short_dir = workspace / "output" / "species_protein_busco_short"
    species_protein_dir.mkdir(parents=True)
    species_tree_summary_dir.mkdir(parents=True)
    busco_short_dir.mkdir(parents=True)

    def write_protein(species: str, count: int) -> None:
        records = [f">{species}_gene{i}\nMPEP\n" for i in range(1, count + 1)]
        (species_protein_dir / f"{species}_pep.fa").write_text("".join(records), encoding="utf-8")

    def write_busco(species: str, complete_pct: float) -> None:
        (busco_short_dir / f"{species}.busco.short.txt").write_text(
            f"C:{complete_pct}%[S:{complete_pct}%,D:0.0%],F:0.0%,M:0.0%,n:100\n",
            encoding="utf-8",
        )

    write_protein("Arabidopsis_thaliana", 2)
    write_protein("Oryza_sativa", 1)
    write_protein("Nepenthes_gracilis", 1)
    write_protein("Dionaea_muscipula", 3)
    write_busco("Arabidopsis_thaliana", 95.0)
    write_busco("Oryza_sativa", 90.0)
    write_busco("Nepenthes_gracilis", 70.0)
    write_busco("Dionaea_muscipula", 99.0)
    (species_tree_summary_dir / "undated_species_tree.nwk").write_text(
        "((Arabidopsis_thaliana:0.1,Oryza_sativa:0.1):0.1,"
        "(Nepenthes_gracilis:0.1,Dionaea_muscipula:0.1):0.1);\n",
        encoding="utf-8",
    )

    completed = _run_core(tmp_path, {"max_orthofinder_core_species": "2"})

    assert completed.returncode == 0, completed.stdout + completed.stderr
    nwkit_args = (tmp_path / "capture" / "nwkit_args.txt").read_text(encoding="utf-8")
    assert "sample" in nwkit_args
    assert "--filter busco_complete_pct:ge:80" in nwkit_args
    assert "--filter num_seq:le:100000" in nwkit_args
    assert "--method max-pd" in nwkit_args
    core_inputs = (tmp_path / "capture" / "input_files_core.txt").read_text(encoding="utf-8")
    assert "Arabidopsis_thaliana.fa" in core_inputs
    assert "Oryza_sativa.fa" in core_inputs
    assert "Nepenthes_gracilis.fa" not in core_inputs
    assert "Dionaea_muscipula.fa" not in core_inputs
    all_inputs = (tmp_path / "capture" / "input_files_all.txt").read_text(encoding="utf-8")
    assert "Nepenthes_gracilis.fa" in all_inputs
    assert "Dionaea_muscipula.fa" in all_inputs
    args = (tmp_path / "capture" / "orthofinder_args.txt").read_text(encoding="utf-8")
    assert "species_tree_core.nwk" in args


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_reuses_existing_busco_outputs_without_retranslation(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    busco_full_dir = workspace / "output" / "species_protein_busco_full"
    busco_short_dir = workspace / "output" / "species_protein_busco_short"
    species_cds_dir.mkdir(parents=True)
    busco_full_dir.mkdir(parents=True)
    busco_short_dir.mkdir(parents=True)

    (species_cds_dir / "Arabidopsis_thaliana_cds.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (species_cds_dir / "Oryza_sativa_cds.fa").write_text(
        ">Oryza_sativa_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (busco_full_dir / "Arabidopsis_thaliana.busco.full.tsv").write_text(
        "# BUSCO version is: test\nbusco1\tComplete\tArabidopsis_thaliana_gene1\t1\t1\t100\n",
        encoding="utf-8",
    )
    (busco_full_dir / "Oryza_sativa.busco.full.tsv").write_text(
        "# BUSCO version is: test\nbusco1\tComplete\tOryza_sativa_gene1\t1\t1\t100\n",
        encoding="utf-8",
    )
    (busco_short_dir / "Arabidopsis_thaliana.busco.short.txt").write_text(
        "C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1\n",
        encoding="utf-8",
    )
    (busco_short_dir / "Oryza_sativa.busco.short.txt").write_text(
        "C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1\n",
        encoding="utf-8",
    )

    completed = _run_core(
        tmp_path,
        {
            "run_species_busco": "1",
            "run_build_species_busco_summary": "0",
            "run_build_species_omark_summary": "0",
            "run_orthofinder": "0",
        },
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Number of protein files for BUSCO: 2" in completed.stdout
    assert "Skipped BUSCO: Arabidopsis_thaliana_cds.fa" in completed.stdout
    assert "Skipped BUSCO: Oryza_sativa_cds.fa" in completed.stdout
    assert "Resolved BUSCO lineage for species set" not in completed.stdout
    assert "Translation started:" not in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_reuses_existing_omark_outputs_without_retranslation(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_cds_dir = workspace / "input" / "species_cds"
    omamer_dir = workspace / "output" / "genome_evolution" / "omamer_search"
    omark_species_dir = workspace / "output" / "genome_evolution" / "omark" / "Arabidopsis_thaliana"
    species_cds_dir.mkdir(parents=True)
    omamer_dir.mkdir(parents=True)
    omark_species_dir.mkdir(parents=True)

    (species_cds_dir / "Arabidopsis_thaliana_cds.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nATGTAA\n",
        encoding="utf-8",
    )
    (omamer_dir / "Arabidopsis_thaliana.omamer").write_text(
        "# query\thog\tscore\nArabidopsis_thaliana.query.fa\tHOG:0000001\t100\n",
        encoding="utf-8",
    )
    (omark_species_dir / "Arabidopsis_thaliana.sum").write_text(
        "#The selected clade was Viridiplantae\n",
        encoding="utf-8",
    )

    completed = _run_core(
        tmp_path,
        {
            "run_species_busco": "0",
            "run_species_omark": "1",
            "run_build_species_omark_summary": "0",
            "run_orthofinder": "0",
        },
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Number of protein files for OMArk: 1" in completed.stdout
    assert "Skipped OMArk: Arabidopsis_thaliana_cds.fa" in completed.stdout
    assert "Translation started:" not in completed.stdout


@pytest.mark.skipif(SYSTEM_BASH_MAJOR < 4, reason="gg_genome_evolution_core.sh requires bash 4+ features such as local -n")
def test_genome_evolution_omark_auto_downloads_database_and_summarizes_results(tmp_path: Path):
    workspace = tmp_path / "workspace"
    species_protein_dir = workspace / "input" / "species_protein"
    species_protein_dir.mkdir(parents=True)
    (species_protein_dir / "Arabidopsis_thaliana_pep.fa").write_text(
        ">Arabidopsis_thaliana_gene1\nMPEPTIDE\n",
        encoding="utf-8",
    )

    fake_db_source = tmp_path / "LUCA.h5"
    fake_db_source.write_text("fake omamer database\n", encoding="utf-8")

    completed = _run_core(
        tmp_path,
        extra_env={
            "run_species_omark": "1",
            "run_build_species_omark_summary": "1",
            "run_orthofinder": "0",
            "GG_OMARK_DB_URL": fake_db_source.resolve().as_uri(),
        },
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    runtime_db = workspace / "downloads" / "omark" / "LUCA.h5"
    assert runtime_db.exists()
    assert runtime_db.read_text(encoding="utf-8") == "fake omamer database\n"

    omark_summary = workspace / "output" / "genome_evolution" / "omark_summary_table" / "omark_summary.tsv"
    assert omark_summary.exists()
    summary_df = pandas.read_csv(omark_summary, sep="\t")
    assert summary_df.shape[0] == 1
    assert summary_df.loc[0, "species"] == "Arabidopsis_thaliana"
    assert summary_df.loc[0, "top_detected_taxid"] == 3702

    omamer_db_paths = (tmp_path / "capture" / "omamer_db_paths.txt").read_text(encoding="utf-8").splitlines()
    assert omamer_db_paths == [str(runtime_db)]
