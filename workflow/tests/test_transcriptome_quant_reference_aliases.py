import csv
import gzip
import hashlib
import json
import os
import random
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from shell_static_helpers import read_text

REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_SCRIPT = REPO_ROOT / "workflow" / "core" / "gg_transcriptome_generation_core.sh"
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
REFERENCE_TOKEN_CASES = [
    (
        "Sarracenia rosea-x-Sarracenia-psittacina",
        "Sarracenia_rosea-x-Sarracenia-psittacina",
        "",
        "Sarracenia_rosea-x-Sarracenia-psittacina",
    ),
    ("Species ONE", "Species_one", "", "Species_ONE"),
    ("Species one (Cultivar A)", "Species_one", "", "Species_one_Cultivar_A"),
    (
        "Sarracenia rosea × Sarracenia psittacina",
        "Sarracenia_rosea_x_Sarracenia_psittacina",
        "",
        "Sarracenia_rosea_Sarracenia_psittacina",
    ),
    ("Species one", "Species_one", "Reference_CUSTOM", "Reference_CUSTOM"),
]


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


@pytest.mark.parametrize(
    ("scientific_name", "species_key", "species_token", "expected_prefix"),
    REFERENCE_TOKEN_CASES,
)
def test_stage_quant_reference_fasta_aliases_preserves_amalgkit_search_token(
    tmp_path, scientific_name, species_key, species_token, expected_prefix
):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa.gz"
    output_dir = tmp_path / "fasta"
    audit_path = tmp_path / "reference_aliases.tsv"
    metadata_path.write_text(
        "run\tscientific_name\tspecies_token\n"
        f"SRR1\t{scientific_name}\t{species_token}\n"
        f"SRR2\t{scientific_name}\t{species_token}\n",
        encoding="utf-8",
    )
    original_metadata = metadata_path.read_bytes()
    reference_path.write_bytes(b">seq\nAAAA\n")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir),
         species_key, str(SUPPORT_DIR), "", str(audit_path)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode == 0, completed.stderr
    expected_name = expected_prefix + "_for_kallisto_index.fa.gz"
    assert completed.stdout.splitlines() == [expected_name]
    assert [path.name for path in output_dir.iterdir()] == [expected_name]
    assert (output_dir / expected_name).read_bytes() == reference_path.read_bytes()
    assert not (output_dir / expected_name).is_symlink()
    assert metadata_path.read_bytes() == original_metadata
    with audit_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["metadata_scientific_name"] == scientific_name
    assert rows[0]["metadata_prefix"] == expected_prefix
    assert rows[0]["alias_path"] == str(output_dir / expected_name)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ("SRR1\tSpecies one\tOne\nSRR2\tSpecies one\tTwo\n", "conflicting species_token"),
        ("SRR1\tSpecies one\tShared\nSRR2\tSpecies one var. X\tShared\n", "colliding"),
        ("SRR1\tSpecies one\tShared\nSRR2\tSpecies one var. X\tshared\n", "colliding"),
        ("SRR1\tSpecies one\t../reference\n", "safe path component"),
    ],
)
def test_stage_quant_reference_fasta_aliases_rejects_invalid_search_tokens(tmp_path, rows, message):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa"
    output_dir = tmp_path / "fasta"
    metadata_path.write_text("run\tscientific_name\tspecies_token\n" + rows, encoding="utf-8")
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir),
         "Species_one", str(SUPPORT_DIR)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode != 0
    assert message in completed.stderr
    assert not output_dir.exists()


def test_quant_reference_rejects_conflicting_taxids_after_missing_value(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "reference.fa"
    output_dir = tmp_path / "fasta"
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\n"
        "SRR1\tSpecies one\t\nSRR2\tSpecies one\t100\nSRR3\tSpecies one\t200\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir),
         "Species_one", str(SUPPORT_DIR)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode != 0
    assert "inconsistent species taxids" in completed.stderr
    assert not output_dir.exists()


@pytest.mark.parametrize("function_name", [
    "stage_quant_reference_fasta_aliases", "stage_amalgkit_merge_metadata_for_species",
])
def test_repeated_species_alias_revalidates_changed_taxid(tmp_path, function_name):
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "run\tscientific_name\ttaxid_species\n"
        "SRR1\tDracocephalum officinale\t39324\n"
        "SRR2\tDracocephalum officinale\t999\n",
        encoding="utf-8",
    )
    taxonomy_db = _taxonomy_db(tmp_path)
    output_path = tmp_path / "output"
    arguments = [str(metadata_path)]
    if function_name == "stage_quant_reference_fasta_aliases":
        reference_path = tmp_path / "reference.fa"
        reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")
        arguments.append(str(reference_path))
    arguments.extend([str(output_path), "Hyssopus_officinalis", str(SUPPORT_DIR), str(taxonomy_db)])
    completed = subprocess.run(
        [sys.executable, "-", *arguments], input=_embedded_python(function_name),
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode != 0
    assert "does not match resolved taxid 39324" in completed.stderr
    assert not output_path.exists()


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

    exact_scientific_name_path = output_dir / "Asimitellaria_furusei_var._furusei_for_kallisto_index.fa.gz"
    canonical_alias = output_dir / "Asimitellaria_furusei_for_kallisto_index.fa.gz"
    original_name_alias = output_dir / "Species_one_original_for_kallisto_index.fa.gz"

    assert exact_scientific_name_path.is_file()
    assert not exact_scientific_name_path.is_symlink()
    assert exact_scientific_name_path.read_bytes() == reference_path.read_bytes()
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
    canonical_alias = output_dir / "Species_one_for_kallisto_index.fa.gz"
    assert canonical_alias.is_file()
    assert not canonical_alias.is_symlink()
    assert canonical_alias.read_bytes() == reference_path.read_bytes()


@pytest.mark.parametrize(
    ("scientific_name", "species_key", "species_token", "expected_prefix"),
    REFERENCE_TOKEN_CASES,
)
def test_stage_quant_reference_fasta_aliases_materializes_discoverable_regular_file(
    tmp_path, scientific_name, species_key, species_token, expected_prefix
):
    quant = pytest.importorskip("amalgkit.quant")
    runtime_utils = pytest.importorskip("amalgkit.runtime_utils")
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"
    source_content = b">seq\nAAAA\n"

    metadata_path.write_text(
        "run\tscientific_name\tspecies_token\n"
        f"SRR1\t{scientific_name}\t{species_token}\n",
        encoding="utf-8",
    )
    reference_path.write_bytes(source_content)

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(output_dir),
            species_key,
            str(SUPPORT_DIR),
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    alias_path = output_dir / (expected_prefix + "_for_kallisto_index.fa.gz")
    assert alias_path.is_file()
    assert not alias_path.is_symlink()
    assert alias_path.read_bytes() == source_content
    assert reference_path.read_bytes() == source_content
    lookup_token = runtime_utils.resolve_species_token(scientific_name, species_token)
    assert lookup_token == expected_prefix
    assert quant.find_species_fasta_files(str(output_dir), lookup_token) == [str(alias_path)]


def test_stage_quant_reference_fasta_aliases_copies_across_filesystems(tmp_path):
    shared_memory_root = Path("/dev/shm")
    if not shared_memory_root.is_dir():
        pytest.skip("cross-filesystem fixture root is unavailable")
    if shared_memory_root.stat().st_dev == tmp_path.stat().st_dev:
        pytest.skip("cross-filesystem fixture root shares the test filesystem")

    metadata_path = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "fasta"
    source_content = b">seq\nAAAA\n"
    metadata_path.write_text(
        "run\tlib_layout\nSRR1\tpaired\n",
        encoding="utf-8",
    )

    with tempfile.TemporaryDirectory(
        prefix="genegalleon-quant-reference-",
        dir=shared_memory_root,
    ) as source_dir:
        reference_path = Path(source_dir) / "Species_one_longestCDS.fa.gz"
        reference_path.write_bytes(source_content)
        completed = subprocess.run(
            [
                sys.executable,
                "-",
                str(metadata_path),
                str(reference_path),
                str(output_dir),
                "Species_one",
                str(SUPPORT_DIR),
            ],
            input=_embedded_python("stage_quant_reference_fasta_aliases"),
            text=True,
            capture_output=True,
            check=False,
        )

        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
        alias_path = output_dir / "Species_one_for_kallisto_index.fa.gz"
        assert alias_path.is_file()
        assert not alias_path.is_symlink()
        assert alias_path.read_bytes() == source_content
        assert reference_path.read_bytes() == source_content
        assert os.stat(alias_path).st_dev != os.stat(reference_path).st_dev


def test_stage_quant_reference_fasta_aliases_rejects_colliding_regular_file(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"
    alias_path = output_dir / "Species_one_for_kallisto_index.fa.gz"
    metadata_path.write_text(
        "run\tlib_layout\nSRR1\tpaired\n",
        encoding="utf-8",
    )
    source_content = b">source\nAAAA\n"
    collision_content = b">collision\nCCCC\n"
    reference_path.write_bytes(source_content)
    output_dir.mkdir()
    alias_path.write_bytes(collision_content)

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(output_dir),
            "Species_one",
            str(SUPPORT_DIR),
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "already exists with different content" in completed.stderr
    assert reference_path.read_bytes() == source_content
    assert alias_path.read_bytes() == collision_content


def test_stage_quant_reference_fasta_aliases_rejects_symlink_collision(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"
    alias_path = output_dir / "Species_one_for_kallisto_index.fa.gz"
    metadata_path.write_text(
        "run\tlib_layout\nSRR1\tpaired\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")
    output_dir.mkdir()
    alias_path.symlink_to(reference_path)

    completed = subprocess.run(
        [
            sys.executable,
            "-",
            str(metadata_path),
            str(reference_path),
            str(output_dir),
            "Species_one",
            str(SUPPORT_DIR),
        ],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "must not be a symlink" in completed.stderr
    assert alias_path.is_symlink()


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
    base_alias = output_dir / "Abies_pinsapo_for_kallisto_index.fa.gz"
    var_alias = output_dir / "Abies_pinsapo_var._marocana_for_kallisto_index.fa.gz"
    assert base_alias.is_file()
    assert var_alias.is_file()
    assert not base_alias.is_symlink()
    assert not var_alias.is_symlink()
    assert base_alias.read_bytes() == reference_path.read_bytes()
    assert var_alias.read_bytes() == reference_path.read_bytes()


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
    alias_path = output_dir / "Dracocephalum_officinale_for_kallisto_index.fa.gz"
    assert alias_path.is_file()
    assert not alias_path.is_symlink()
    assert alias_path.read_bytes() == reference_path.read_bytes()
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


def test_stage_amalgkit_merge_metadata_collapses_explicit_tokens_with_species_aliases(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    staged_path = tmp_path / "metadata.merge.tsv"
    metadata_path.write_text(
        "run\tscientific_name\tspecies_token\n"
        "SRR1\tAbies pinsapo\tBase_CUSTOM\n"
        "SRR2\tAbies pinsapo var. marocana\tVariety_CUSTOM\n",
        encoding="utf-8",
    )
    original_metadata = metadata_path.read_bytes()
    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(staged_path), "Abies_pinsapo", str(SUPPORT_DIR), ""],
        input=_embedded_python("stage_amalgkit_merge_metadata_for_species"),
        text=True, capture_output=True, check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert metadata_path.read_bytes() == original_metadata
    with staged_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["scientific_name"] for row in rows] == ["Abies pinsapo", "Abies pinsapo"]
    assert [row["species_token"] for row in rows] == ["Abies_pinsapo", "Abies_pinsapo"]
    assert [row["run"] for row in rows] == ["SRR1", "SRR2"]


@pytest.mark.parametrize("species_token", ["", "Hybrid_CUSTOM"])
def test_amalgkit_quant_merge_with_mixed_case_hybrid_reference(tmp_path, species_token):
    metadata_utils = pytest.importorskip("amalgkit.metadata_utils")
    runtime_utils = pytest.importorskip("amalgkit.runtime_utils")
    pandas = pytest.importorskip("pandas")
    scientific_name = "Sarracenia rosea-x-Sarracenia-psittacina"
    species_key = scientific_name.replace(" ", "_")
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / (species_key + "_longestCDS.fa.gz")
    fasta_dir = tmp_path / "fasta"
    audit_path = tmp_path / "reference_aliases.tsv"
    rng = random.Random(47)
    sequence = "".join(rng.choices("ACGT", k=1000))
    with gzip.open(reference_path, "wt") as handle:
        handle.write(">tx1\n" + sequence + "\n")
    metadata_rows = []
    manifest_runs = []
    for run in ("R1", "R2"):
        fastq_path = tmp_path / "getfastq" / run / (run + ".amalgkit.fastq.gz")
        fastq_path.parent.mkdir(parents=True)
        with gzip.open(fastq_path, "wt") as handle:
            for index in range(2000):
                start = rng.randrange(len(sequence) - 100)
                handle.write(f"@read{index}\n{sequence[start:start + 100]}\n+\n" + "I" * 100 + "\n")
        stats = fastq_path.parent / "getfastq_stats.tsv"
        stats.write_text(
            f"run\tnum_written\tbp_fastp_in\tnum_fastp_in\tnum_fastp_out\n{run}\t2000\t200000\t2000\t1900\n",
            encoding="utf-8",
        )
        contracts = [
            {"path": path.relative_to(tmp_path / "getfastq").as_posix(), "size": path.stat().st_size,
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in (fastq_path, stats)
        ]
        manifest_runs.append({"run": run, "status": "complete", "files": contracts[:1], "stats": contracts[1]})
        metadata_rows.append(dict(
            run=run, scientific_name=scientific_name, species_token=species_token,
            lib_layout="single", total_spots=2000, total_bases=200000, spot_length=100,
            exclusion="no", is_sampled="yes", instrument="Illumina", taxid="1",
        ))
    metadata_utils.Metadata.from_DataFrame(pandas.DataFrame(metadata_rows)).df.to_csv(
        metadata_path, sep="\t", index=False,
    )
    original_metadata = metadata_path.read_bytes()
    manifest = tmp_path / "getfastq" / "getfastq_completion.json"
    manifest.write_text(json.dumps({"schema_version": 4, "read_source": "amalgkit", "status": "complete",
                                    "run_count": len(manifest_runs), "runs": manifest_runs}))

    def run_command(command, **kwargs):
        completed = subprocess.run(
            command, cwd=tmp_path, text=True, capture_output=True, check=False, timeout=90, **kwargs,
        )
        assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
        return completed

    run_command(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(fasta_dir), species_key,
         str(SUPPORT_DIR), "", str(audit_path)],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
    )
    run_command([
        "amalgkit", "quant", "--out_dir", str(tmp_path), "--metadata", str(metadata_path),
        "--threads", "1", "--clean_fastq", "no", "--fasta_dir", str(fasta_dir),
        "--build_index", "yes", "--quant_backend", "kallisto",
    ])
    lookup_token = runtime_utils.resolve_species_token(scientific_name, species_token)
    assert (tmp_path / "index" / (lookup_token + ".idx")).stat().st_size > 0
    run_command([
        sys.executable, str(SUPPORT_DIR / "validate_transcriptome_quant_outputs.py"),
        "--metadata", str(metadata_path), "--quant-root", str(tmp_path / "quant"),
    ])
    for run in ("R1", "R2"):
        run_dir = tmp_path / "quant" / run
        with (run_dir / (run + "_abundance.tsv")).open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["target_id"] for row in rows] == ["tx1"]
        assert float(rows[0]["est_counts"]) > 0
        assert json.loads((run_dir / (run + "_run_info.json")).read_text())["p_pseudoaligned"] > 0

    # Simulate one unfinished run after an interrupted species job.
    first = tmp_path / "quant" / "R1" / "R1_abundance.tsv"
    first_identity = first.read_bytes(), first.stat().st_mtime_ns
    index = tmp_path / "index" / (lookup_token + ".idx")
    index_mtime = index.stat().st_mtime_ns
    shutil.rmtree(tmp_path / "quant" / "R2")
    run_command([
        "bash", "-c", 'set -euo pipefail; GG_CORE_SOURCE_ONLY=1 source "$1"; source "$2"; '
        "quantification_only=1 prepare_transcriptome_quant_work_dirs",
        "_", str(CORE_SCRIPT), str(SUPPORT_DIR / "gg_util.sh"),
    ])
    proof = tmp_path / "verified_fastq.json"
    proof.touch(mode=0o600)
    run_command([
        "bash", "-c", 'set -euo pipefail; GG_CORE_SOURCE_ONLY=1 source "$1"; '
        'validate_amalgkit_getfastq_completion_manifest "$2" "$3" "$4"',
        "_", str(CORE_SCRIPT), str(manifest), str(metadata_path), str(proof),
    ])
    assert len(json.loads(proof.read_text())["files"]) == 2
    resume = run_command([
        "amalgkit", "quant", "--out_dir", str(tmp_path), "--metadata", str(metadata_path),
        "--threads", "1", "--clean_fastq", "no", "--fasta_dir", str(fasta_dir),
        "--build_index", "yes", "--quant_backend", "kallisto", "--redo", "no", "--verified_fastq_proof", str(proof),
    ])
    assert "The output will not be overwritten" in resume.stdout
    assert first_identity == (first.read_bytes(), first.stat().st_mtime_ns)
    assert index.stat().st_mtime_ns == index_mtime
    assert (tmp_path / "quant" / "R2" / "R2_abundance.tsv").is_file()

    retained = tmp_path / "retained"
    (tmp_path / "getfastq").rename(retained)
    for run in ("R1", "R2"):
        (retained / run / (run + ".amalgkit.fastq.gz")).unlink()
    run_command([
        "bash", "-c", 'set -euo pipefail; GG_CORE_SOURCE_ONLY=1 source "$1"; link_transcriptome_merge_getfastq "$2"',
        "_", str(CORE_SCRIPT), str(retained),
    ])
    merge_metadata = tmp_path / "metadata.merge.tsv"
    run_command(
        [sys.executable, "-", str(metadata_path), str(merge_metadata), species_key, str(SUPPORT_DIR), ""],
        input=_embedded_python("stage_amalgkit_merge_metadata_for_species"),
    )
    run_command([
        "amalgkit", "merge", "--out_dir", str(tmp_path), "--metadata", str(merge_metadata), "--threads", "1",
    ])
    resolved = run_command(
        [sys.executable, "-", str(merge_metadata), str(tmp_path / "merge"), species_key, str(SUPPORT_DIR)],
        input=_embedded_python("resolve_amalgkit_merge_output_prefix"),
    )
    merge_prefix = resolved.stdout.strip()
    assert merge_prefix == "Sarracenia_rosea-x-sarracenia-psittacina"
    for suffix in ("eff_length", "est_counts", "tpm"):
        with (tmp_path / "merge" / merge_prefix / f"{merge_prefix}_{suffix}.tsv").open(
            encoding="utf-8", newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert set(reader.fieldnames) == {"target_id", "R1", "R2"}
            rows = list(reader)
        assert [row["target_id"] for row in rows] == ["tx1"]
        assert all(float(rows[0][run]) > 0 for run in ("R1", "R2"))
    with (tmp_path / "merge" / "metadata.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["run"] for row in rows} == {"R1", "R2"}
    assert all(float(row["num_fastp_in"]) == 2000 and float(row["num_fastp_out"]) == 1900 for row in rows)
    assert metadata_path.read_bytes() == original_metadata

    # The native tool reads exactly the same files for the retained ZIP.
    merge_before = {
        path.relative_to(tmp_path / "merge"): path.read_bytes()
        for path in (tmp_path / "merge").rglob("*.tsv")
    }
    archive_tool = str(SUPPORT_DIR / "archive_transcriptome_quant.py")
    logical_quant = tmp_path / "output" / "transcriptome_quantification" / "amalgkit_quant" / species_key
    run_command([sys.executable, archive_tool, "pack", "--raw", str(logical_quant),
                 "--source", str(tmp_path / "quant")])
    shutil.rmtree(tmp_path / "quant")
    shutil.rmtree(tmp_path / "merge")
    run_command([sys.executable, archive_tool, "stage-merge", "--raw", str(logical_quant),
                 "--metadata", str(metadata_path), "--destination", str(tmp_path / "quant")])
    run_command([
        "amalgkit", "merge", "--out_dir", str(tmp_path), "--metadata", str(merge_metadata), "--threads", "1",
    ])
    assert {
        path.relative_to(tmp_path / "merge"): path.read_bytes()
        for path in (tmp_path / "merge").rglob("*.tsv")
    } == merge_before



def test_quant_alias_staging_failure_is_not_masked_by_process_substitution():
    text = read_text(CORE_SCRIPT)
    assert "if ! staged_quant_reference_aliases=$(stage_quant_reference_fasta_aliases" in text
    assert 'done <<< "${staged_quant_reference_aliases}"' in text
    assert "done < <(\n      stage_quant_reference_fasta_aliases" not in text


@pytest.mark.parametrize(
    ("scientific_name", "species_token", "expected_prefix"),
    [
        ("Asimitellaria furusei var. furusei", "", "Asimitellaria_furusei_var._furusei"),
        ("Sarracenia rosea-x-Sarracenia-psittacina", "", "Sarracenia_rosea-x-Sarracenia-psittacina"),
        ("Species one (Cultivar A)", "", "Species_one_Cultivar_A"),
        ("Species one", "Merged_CUSTOM", "Merged_CUSTOM"),
    ],
)
def test_resolve_amalgkit_merge_output_prefix_uses_exact_metadata_scientific_name(
    tmp_path, scientific_name, species_token, expected_prefix
):
    metadata_path = tmp_path / "metadata.tsv"
    merge_dir = tmp_path / "merge"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tspecies_token\tscientific_name_original",
                f"SRR1\t{scientific_name}\t{species_token}\tSpecies one/original",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    alias_dir = merge_dir / expected_prefix
    alias_dir.mkdir(parents=True)
    (alias_dir / (expected_prefix + "_eff_length.tsv")).write_text("target_id\tSRR1\nx\t1\n", encoding="utf-8")

    command = "\n".join(
        [
            _function_script("resolve_amalgkit_merge_output_prefix"),
            "resolve_amalgkit_merge_output_prefix "
            f"{shlex.quote(str(metadata_path))} "
            f"{shlex.quote(str(merge_dir))} "
            f"Species_one {shlex.quote(str(SUPPORT_DIR))}",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert completed.stdout.strip() == expected_prefix


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
            f"Species_one {shlex.quote(str(SUPPORT_DIR))}",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0


@pytest.mark.parametrize("substitute", ["", "_extra"])
def test_amalgkit_option_probe_is_exact_and_reuses_help(tmp_path, substitute):
    counts = tmp_path / "help_calls"
    command = r'''
source "$1"
counts=$2
suffix=$3
amalgkit() {
    printf '%s\n' "$*" >> "$counts"
    printf 'usage: amalgkit getfastq [--read_source_format%s {auto,fastq,sra}]\n' "$suffix"
}
if require_amalgkit_supported_options getfastq --read_source_format; then
    require_amalgkit_supported_options getfastq --read_source_format
else
    exit 1
fi
'''
    result = subprocess.run(
        ["bash", "-euc", command, "_", str(SUPPORT_DIR / "gg_amalgkit_runtime.sh"),
         str(counts), substitute],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == (1 if substitute else 0), result.stderr
    assert counts.read_text().splitlines() == ["getfastq --help"]


@pytest.mark.parametrize(("mode", "choices", "expected"), [
    ("auto", "auto|yes|no", 0),
    ("auto", "yes|no", 1),
    ("no", "yes|no", 0),
])
def test_amalgkit_option_probe_requires_auto_fallback_support(mode, choices, expected):
    command = r'''
source "$1"
amalgkit_allow_read_format_fallback=$2
choices=$3
required=()
gg_amalgkit_getfastq_options required
amalgkit() {
    printf '%s\n' --download_lock_dir --ncbi_metadata_max_concurrency --ncbi_download_max_concurrency \
        --aws_download_max_concurrency --gcp_download_max_concurrency --ena --rrna_filter_jobs \
        --rrna_filter_chunk_spots --rrna_filter_memory_limit --read_source_format
    printf '%s %s\n' --allow_read_format_fallback "$choices"
}
require_amalgkit_supported_options getfastq "${required[@]}"
'''
    result = subprocess.run(
        ["bash", "-euc", command, "_", str(SUPPORT_DIR / "gg_amalgkit_runtime.sh"), mode, choices],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == expected, result.stderr
