import csv
import re
import shlex
import sqlite3
import subprocess
import sys
from pathlib import Path

from shell_static_helpers import read_text

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_SCRIPT = REPO_ROOT / "workflow" / "core" / "gg_transcriptome_generation_core.sh"
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"


def _taxonomy_db(tmp_path: Path) -> Path:
    db_path = tmp_path / "taxa.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.execute(
            "CREATE TABLE species (taxid INTEGER PRIMARY KEY, spname TEXT COLLATE NOCASE, rank TEXT, track TEXT)"
        )
        conn.execute(
            "CREATE TABLE synonym (taxid INTEGER, spname TEXT COLLATE NOCASE, PRIMARY KEY (spname, taxid))"
        )
        conn.executemany(
            "INSERT INTO species (taxid, spname, rank, track) VALUES (?, ?, 'species', '')",
            [
                (39324, "Dracocephalum officinale"),
                (100, "Species one"),
                (200, "Species two"),
                (300, "Accepted species"),
                (400, "Ambiguity one"),
                (401, "Ambiguity two"),
            ],
        )
        conn.executemany(
            "INSERT INTO synonym (taxid, spname) VALUES (?, ?)",
            [
                (39324, "Hyssopus officinalis"),
                (300, "Canonical species"),
                (400, "Ambiguous species"),
                (401, "Ambiguous species"),
            ],
        )
    return db_path


def _function_body(text: str, function_name: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(function_name)}\(\)\s*\{{", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"Function not found: {function_name}")
    start = match.start()
    next_match = re.search(r"^\s*[A-Za-z_][A-Za-z0-9_]*\(\)\s*\{", text[match.end():], re.MULTILINE)
    if next_match is None:
        return text[start:]
    return text[start:match.end() + next_match.start()]


def _embedded_python(function_name: str) -> str:
    text = read_text(CORE_SCRIPT)
    body = _function_body(text, function_name)
    match = re.search(r"<<'PY'\n(.*?)\nPY", body, re.DOTALL)
    if match is None:
        raise AssertionError(f"Embedded Python not found in function: {function_name}")
    return match.group(1)


def _function_script(function_name: str) -> str:
    text = read_text(CORE_SCRIPT)
    return _function_body(text, function_name)


def test_stage_quant_reference_fasta_aliases_uses_exact_metadata_scientific_name(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tscientific_name_original",
                "SRR1\tAsimitellaria furusei var. furusei\tSpecies one/original",
                "SRR2\tAsimitellaria furusei var. furusei\tSpecies one/original",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir), "Asimitellaria_furusei", str(SUPPORT_DIR)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    exact_scientific_name_path = output_dir / "Asimitellaria_furusei_var._furusei_for_kallisto_index.fasta"
    canonical_alias = output_dir / "Asimitellaria_furusei_for_kallisto_index.fasta"
    original_name_alias = output_dir / "Species_one_original_for_kallisto_index.fasta"

    assert exact_scientific_name_path.is_symlink()
    assert exact_scientific_name_path.resolve() == reference_path.resolve()
    assert not canonical_alias.exists()
    assert not original_name_alias.exists()


def test_stage_quant_reference_fasta_aliases_still_stages_canonical_prefix_without_metadata_columns(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tlib_layout",
                "SRR1\tpaired",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir), "Species_one", str(SUPPORT_DIR)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    canonical_alias = output_dir / "Species_one_for_kallisto_index.fasta"
    assert canonical_alias.is_symlink()
    assert canonical_alias.resolve() == reference_path.resolve()


def test_stage_quant_reference_fasta_aliases_stages_base_and_infraspecific_prefixes(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Abies_pinsapo_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name",
                "SRR1\tAbies pinsapo",
                "SRR2\tAbies pinsapo var. marocana",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir), "Abies_pinsapo", str(SUPPORT_DIR)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    base_alias = output_dir / "Abies_pinsapo_for_kallisto_index.fasta"
    var_alias = output_dir / "Abies_pinsapo_var._marocana_for_kallisto_index.fasta"
    assert base_alias.is_symlink()
    assert var_alias.is_symlink()
    assert base_alias.resolve() == reference_path.resolve()
    assert var_alias.resolve() == reference_path.resolve()


def test_stage_quant_reference_fasta_aliases_rejects_different_species_names(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name",
                "SRR1\tSpecies one",
                "SRR2\tSpecies two",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")
    taxonomy_db = _taxonomy_db(tmp_path)

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(output_dir),
            "Species_one",
            str(SUPPORT_DIR),
            str(taxonomy_db),
            "",
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "is not an alias of GeneGalleon species_key" in completed.stderr


def test_stage_quant_reference_fasta_aliases_accepts_cross_genus_shared_taxid_and_writes_audit(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Hyssopus_officinalis_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"
    audit_path = tmp_path / "reference_aliases.tsv"
    taxonomy_db = _taxonomy_db(tmp_path)
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\ttaxid\n"
        "SRR1\tDracocephalum officinale\t39324\t39324\n"
        "SRR2\tDracocephalum officinale\t39324\t39324\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(output_dir),
            "Hyssopus_officinalis",
            str(SUPPORT_DIR),
            str(taxonomy_db),
            str(audit_path),
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    alias_path = output_dir / "Dracocephalum_officinale_for_kallisto_index.fasta"
    assert alias_path.is_symlink()
    assert alias_path.resolve() == reference_path.resolve()
    with audit_path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["canonical_species_key"] == "Hyssopus_officinalis"
    assert rows[0]["declared_taxid"] == "39324"
    assert rows[0]["resolved_taxid"] == "39324"
    assert rows[0]["resolution_method"] == "shared_species_taxid"


def test_stage_quant_reference_fasta_aliases_rejects_declared_taxid_mismatch(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa"
    taxonomy_db = _taxonomy_db(tmp_path)
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\nSRR1\tDracocephalum officinale\t999\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(tmp_path / "fasta"),
            "Hyssopus_officinalis",
            str(SUPPORT_DIR),
            str(taxonomy_db),
            "",
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "does not match resolved taxid 39324" in completed.stderr


def test_stage_quant_reference_fasta_aliases_rejects_cross_name_when_taxonomy_db_is_missing(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa"
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\nSRR1\tDracocephalum officinale\t39324\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(tmp_path / "fasta"),
            "Hyssopus_officinalis",
            str(SUPPORT_DIR),
            str(tmp_path / "missing.sqlite"),
            "",
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "Taxonomy DB is required" in completed.stderr


def test_stage_quant_reference_fasta_aliases_rejects_ambiguous_name(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa"
    taxonomy_db = _taxonomy_db(tmp_path)
    metadata_path.write_text(
        "run\tscientific_name\nSRR1\tAmbiguous species\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(tmp_path / "fasta"),
            "Canonical_species",
            str(SUPPORT_DIR),
            str(taxonomy_db),
            "",
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "resolved ambiguously" in completed.stderr


def test_stage_amalgkit_merge_metadata_for_species_collapses_infraspecific_names(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    staged_path = tmp_path / "metadata.merge.tsv"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\texclusion",
                "SRR1\tAbies pinsapo\tno",
                "SRR2\tAbies pinsapo var. marocana\tno",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    command = "\n".join(
        [
            _function_script("stage_amalgkit_merge_metadata_for_species"),
            "stage_amalgkit_merge_metadata_for_species "
            f"{shlex.quote(str(metadata_path))} "
            f"{shlex.quote(str(staged_path))} "
            "Abies_pinsapo "
            f"{shlex.quote(str(SUPPORT_DIR))}",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    with staged_path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["scientific_name"] for row in rows] == ["Abies pinsapo", "Abies pinsapo"]


def test_stage_amalgkit_merge_metadata_for_species_collapses_cross_genus_shared_taxid(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    staged_path = tmp_path / "metadata.merge.tsv"
    taxonomy_db = _taxonomy_db(tmp_path)
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\texclusion\n"
        "SRR1\tDracocephalum officinale\t39324\tno\n",
        encoding="utf-8",
    )

    command = "\n".join(
        [
            _function_script("stage_amalgkit_merge_metadata_for_species"),
            "stage_amalgkit_merge_metadata_for_species "
            f"{shlex.quote(str(metadata_path))} "
            f"{shlex.quote(str(staged_path))} "
            "Hyssopus_officinalis "
            f"{shlex.quote(str(SUPPORT_DIR))} "
            f"{shlex.quote(str(taxonomy_db))}",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    with staged_path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["scientific_name"] for row in rows] == ["Hyssopus officinalis"]


def test_quant_alias_staging_failure_is_not_masked_by_process_substitution():
    text = read_text(CORE_SCRIPT)
    assert "if ! staged_quant_reference_aliases=$(stage_quant_reference_fasta_aliases" in text
    assert 'done <<< "${staged_quant_reference_aliases}"' in text
    assert "done < <(\n      stage_quant_reference_fasta_aliases" not in text


def test_resolve_amalgkit_merge_output_prefix_uses_exact_metadata_scientific_name(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    merge_dir = tmp_path / "merge"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tscientific_name_original",
                "SRR1\tAsimitellaria furusei var. furusei\tSpecies one/original",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    alias_dir = merge_dir / "Asimitellaria_furusei_var._furusei"
    alias_dir.mkdir(parents=True)
    (alias_dir / "Asimitellaria_furusei_var._furusei_eff_length.tsv").write_text("target_id\tSRR1\nx\t1\n", encoding="utf-8")

    command = "\n".join(
        [
            _function_script("resolve_amalgkit_merge_output_prefix"),
            "resolve_amalgkit_merge_output_prefix "
            f"{shlex.quote(str(metadata_path))} "
            f"{shlex.quote(str(merge_dir))} "
            "Species_one",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert completed.stdout.strip() == "Asimitellaria_furusei_var._furusei"


def test_resolve_amalgkit_merge_output_prefix_does_not_fall_back_to_canonical_prefix(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    merge_dir = tmp_path / "merge"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name",
                "SRR1\tSpecies one strain X",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    canonical_dir = merge_dir / "Species_one"
    canonical_dir.mkdir(parents=True)
    (canonical_dir / "Species_one_eff_length.tsv").write_text("target_id\tSRR1\nx\t1\n", encoding="utf-8")

    command = "\n".join(
        [
            _function_script("resolve_amalgkit_merge_output_prefix"),
            "resolve_amalgkit_merge_output_prefix "
            f"{shlex.quote(str(metadata_path))} "
            f"{shlex.quote(str(merge_dir))} "
            "Species_one",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
