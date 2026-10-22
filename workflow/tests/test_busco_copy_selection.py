import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support/batch_extract_busco_fasta.py"


def run_extraction(tmp_path, records, rows, *, policy="best_score", mode="protein", extra=(), reverse=False, mutate=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    full = tmp_path / "full"
    full.mkdir(exist_ok=True)
    fasta_paths = []
    for species, sequences in records.items():
        fasta = tmp_path / f"{species}.fa.gz"
        with gzip.open(fasta, "wt") as handle:
            for gene, sequence in reversed(sequences) if reverse else sequences:
                handle.write(f">{gene} description\n{sequence}\n")
        fasta_paths.append(fasta)
        table = rows[species]
        (full / f"{species}.busco.full.tsv").write_text(
            "# Busco id\tStatus\tSequence\tScore\tLength\n" +
            "".join("\t".join(row) + "\n" for row in (reversed(table) if reverse else table)))
    summary = tmp_path / "summary.tsv"
    marker_ids = sorted({row[0] for table in rows.values() for row in table})
    with summary.open("w") as handle:
        handle.write("busco_id\torthodb_url\tdescription\t" + "\t".join(records) + "\n")
        for marker in marker_ids:
            cells = []
            for species in records:
                genes = [row[2].split(":")[0] for row in rows[species] if row[0] == marker and row[2] != "-"]
                cells.append(",".join(genes) if genes else "-")
            handle.write(marker + "\t-\t-\t" + "\t".join(cells) + "\n")
    fasta_list = tmp_path / "fastas.txt"
    fasta_list.write_text("".join(f"{path}\n" for path in (reversed(fasta_paths) if reverse else fasta_paths)))
    species_map = tmp_path / "species.tsv"
    species_map.write_text("species\tfasta\n" + "".join(f"{species}\t{path}\n"
                           for species, path in zip(records, fasta_paths, strict=True)))
    report = tmp_path / "selection.tsv"
    command = [sys.executable, str(SCRIPT), "--summary", str(summary), "--fasta-list", str(fasta_list),
               "--output-dir", str(tmp_path / "output"), "--suffix", ".fa",
               "--report", str(tmp_path / "counts.tsv"), "--mode", "single-copy",
               "--selection-report", str(report)]
    if policy is not None:
        command += ["--duplicate-policy", policy]
    if policy == "best_score":
        command += ["--busco-full-dir", str(full), "--species-fasta-map", str(species_map),
                    "--sequence-mode", mode, "--require-complete"]
    if mutate:
        mutate(tmp_path)
    result = subprocess.run(command + list(extra), capture_output=True, text=True)
    audit = []
    if report.exists():
        with report.open() as handle:
            audit = list(csv.DictReader(handle, delimiter="\t"))
    return result, audit


@pytest.mark.parametrize("policy", [None, "exclude"])
def test_default_excludes_duplicates_and_retains_existing_single_copy_behavior(tmp_path, policy):
    result, audit = run_extraction(tmp_path, {"A": [("single", "MAAA"), ("d1", "MCCC"), ("d2", "MDDD")]},
        {"A": [("b1", "Complete", "single", "20", "4"),
               ("b2", "Duplicated", "d1", "100", "4"), ("b2", "Duplicated", "d2", "90", "4")]},
        policy=policy)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output/b1.fa").read_text() == ">single description\nMAAA\n"
    assert not (tmp_path / "output/b2.fa").exists()
    assert audit == []


def test_quality_and_original_score_policies_cannot_be_enabled_together(tmp_path):
    result, audit = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"),
               ("b1", "Duplicated", "g2", "90", "4")]},
        extra=("--copy-policy", "best-score"))
    assert result.returncode == 2 and "conflicts" in result.stderr
    assert not (tmp_path / "output").exists() and audit == []


def test_original_score_policy_rejects_report_input_overlap(tmp_path):
    records = {"A": [("g1", "MAAA"), ("g2", "MCCC")]}
    rows = {"A": [("b1", "Duplicated", "g1", "100", "4"),
                  ("b1", "Duplicated", "g2", "90", "4")]}
    table = tmp_path / "full/A.busco.full.tsv"
    result, _ = run_extraction(tmp_path, records, rows, policy="exclude",
        extra=("--copy-policy", "best-score", "--busco-full-dir", str(tmp_path / "full"),
               "--selection-report", str(table)))
    assert result.returncode == 2 and "overlap" in result.stderr
    assert not (tmp_path / "output").exists()
    assert "g1\t100\t4" in table.read_text()


def test_best_score_selects_per_species_and_preserves_input_and_all_copy_summary(tmp_path):
    records = {"A": [("a1", "MAAA"), ("a2", "MCCC")], "B": [("b1", "MDDD"), ("b2", "MEEE")]}
    rows = {"A": [("b1", "Duplicated", "a1", "110", "4"), ("b1", "Duplicated", "a2", "100", "4")],
            "B": [("b1", "Duplicated", "b1", "90", "4"), ("b1", "Duplicated", "b2", "120", "4")]}
    result, audit = run_extraction(tmp_path, records, rows)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output/b1.fa").read_text() == ">a1 description\nMAAA\n>b2 description\nMEEE\n"
    chosen = [row for row in audit if row["selected"] == "1"]
    assert [(row["species"], row["gene_id"], row["reason"], row["score_margin"]) for row in chosen] == [
        ("A", "a1", "highest_score", "10.0"), ("B", "b2", "highest_score", "30.0")]
    assert all(row["status"] == "Duplicated" for row in audit)
    assert "a1,a2\tb1,b2" in (tmp_path / "summary.tsv").read_text()
    assert "a1\t110" in (tmp_path / "full/A.busco.full.tsv").read_text()


def test_quality_ties_and_audit_are_independent_of_input_order(tmp_path):
    records = {"A": [("stop", "MA*AA"), ("ambiguous", "MAXAA"), ("z", "MAAAA*"), ("a", "MAAAA")]}
    rows = {"A": [("b1", "Duplicated", gene, score, "5")
                   for gene, score in [("stop", "120"), ("ambiguous", "100"), ("z", "100"), ("a", "100")]]}
    first, audit = run_extraction(tmp_path / "first", records, rows)
    second, reordered = run_extraction(tmp_path / "second", records, rows, reverse=True)
    assert first.returncode == second.returncode == 0, first.stderr + second.stderr
    assert audit == reordered
    chosen = next(row for row in audit if row["selected"] == "1")
    assert chosen["gene_id"] == "a" and chosen["reason"] == "gene_id_tiebreak"
    assert {row["gene_id"]: row["reason"] for row in audit} == {
        "a": "gene_id_tiebreak", "ambiguous": "higher_ambiguous_fraction",
        "stop": "internal_stop", "z": "gene_id_tiebreak"}
    assert next(row for row in audit if row["gene_id"] == "z")["internal_stops"] == "0"


@pytest.mark.parametrize("genes,expected,reason", [
    ([("stop", "MA*AA", "120"), ("clean", "MAAAA", "100")], "clean", "no_internal_stop"),
    ([("x", "MAXAA", "100"), ("clean", "MAAAA", "100")], "clean", "lowest_ambiguous_fraction"),
    ([("stop1", "MA*AA", "120"), ("stop2", "MA**A", "100")], "stop1", "highest_score"),
])
def test_quality_preference_and_all_stop_candidates_remain_auditable(tmp_path, genes, expected, reason):
    result, audit = run_extraction(tmp_path, {"A": [(gene, sequence) for gene, sequence, _ in genes]},
        {"A": [("b1", "Duplicated", gene, score, "5") for gene, _, score in genes]})
    assert result.returncode == 0, result.stderr
    chosen = next(row for row in audit if row["selected"] == "1")
    assert (chosen["gene_id"], chosen["reason"]) == (expected, reason)


@pytest.mark.parametrize("mode", ["cds", "protein"])
def test_cds_and_protein_selection_use_same_source_gene_and_keep_original_sequence(tmp_path, mode):
    sequence = "ATGAAATAA" if mode == "cds" else "MK*"
    result, audit = run_extraction(tmp_path, {"A": [("g1", sequence), ("g2", sequence)]},
        {"A": [("b1", "Duplicated", "g1:1-9" if mode == "cds" else "g1", "100", "3"),
               ("b1", "Duplicated", "g2:1-9" if mode == "cds" else "g2", "90", "3")]}, mode=mode)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output/b1.fa").read_text() == f">g1 description\n{sequence}\n"
    assert audit[0]["gene_id"] == "g1" and audit[0]["internal_stops"] == "0"


def test_cds_quality_honors_species_specific_genetic_code(tmp_path):
    codes = tmp_path / "codes.tsv"
    codes.write_text("species\tgenetic_code\nA\t6\n")
    result, audit = run_extraction(tmp_path, {"A": [("g1", "ATGTAATAA"), ("g2", "ATGAAAAAA")]},
        {"A": [("b1", "Duplicated", "g1", "100", "3"), ("b1", "Duplicated", "g2", "90", "3")]},
        mode="cds", extra=("--species-genetic-codes", str(codes)))
    assert result.returncode == 0, result.stderr
    assert audit[0]["selected"] == "1" and audit[0]["internal_stops"] == "0"


def test_ambiguous_codon_that_encodes_stop_is_detected(tmp_path):
    result, audit = run_extraction(tmp_path, {"A": [("stop", "ATGTARAAA"), ("clean", "ATGAAAAAA")]},
        {"A": [("b1", "Duplicated", "stop", "100", "3"), ("b1", "Duplicated", "clean", "90", "3")]},
        mode="cds")
    assert result.returncode == 0, result.stderr
    assert next(row for row in audit if row["gene_id"] == "stop")["internal_stops"] == "1"
    assert next(row for row in audit if row["selected"] == "1")["gene_id"] == "clean"


def test_repeated_transcript_hits_select_one_source_gene_and_fragment_is_not_rescued(tmp_path):
    result, audit = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC"), ("fragment", "MDDD")]},
        {"A": [("b1", "Duplicated", "g1:1-9", "90", "3"),
               ("b1", "Duplicated", "g1:4-12", "100", "3"),
               ("b1", "Duplicated", "g2", "95", "4"),
               ("b1", "Fragmented", "fragment", "110", "4")]})
    assert result.returncode == 0, result.stderr
    assert len(audit) == 3
    chosen = next(row for row in audit if row["selected"] == "1")
    assert chosen["gene_id"] == "g1" and chosen["busco_sequence_ids"] == "g1:1-9,g1:4-12"
    assert next(row for row in audit if row["gene_id"] == "fragment")["reason"] == "not_complete"
    assert (tmp_path / "output/b1.fa").read_text().count(">") == 1


@pytest.mark.parametrize("score", ["-", "nan", "inf", ""])
def test_invalid_duplicate_score_fails_without_length_fallback(tmp_path, score):
    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", score, "4"), ("b1", "Duplicated", "g2", "90", "4")]})
    assert result.returncode != 0 and "Invalid BUSCO score" in result.stderr


@pytest.mark.parametrize("records", [[("g1", "MAAA")], [("g1", "MAAA"), ("g1", "MCCC"), ("g2", "MDDD")]])
def test_missing_or_repeated_candidate_id_fails(tmp_path, records):
    result, _ = run_extraction(tmp_path, {"A": records},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")]})
    assert result.returncode != 0
    assert "missing" in result.stderr or "Repeated candidate FASTA ID" in result.stderr


def test_strict_mode_conflicts_with_best_score(tmp_path):
    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA")]},
        {"A": [("b1", "Complete", "g1", "100", "4")]}, extra=("--strict",))
    assert result.returncode != 0 and "without --strict" in result.stderr


def test_fragment_only_duplicate_groups_are_not_rescued_and_single_id_behavior_is_preserved(tmp_path):
    result, audit = run_extraction(tmp_path, {"A": [("f1", "MAAA"), ("f2", "MCCC"), ("single", "MDDD")]},
        {"A": [("b1", "Fragmented", "f1", "90", "4"), ("b1", "Fragmented", "f2", "80", "4"),
               ("b2", "Fragmented", "single", "100", "4")]})
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "output/b1.fa").exists()
    assert (tmp_path / "output/b2.fa").is_file()
    assert all(row["selected"] == "0" and row["reason"] == "not_complete" for row in audit)


def test_summary_full_table_mismatch_fails(tmp_path):
    stale = tmp_path / "stale.tsv"
    stale.write_text("busco_id\torthodb_url\tdescription\tA\nb1\t-\t-\twrong,g2\n")
    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")]},
        extra=("--summary", str(stale)))
    assert result.returncode != 0 and "summary/full-table ID mismatch" in result.stderr


def test_species_map_mismatch_fails(tmp_path):
    mapping = tmp_path / "wrong.tsv"
    mapping.write_text(f"species\tfasta\nB\t{tmp_path / 'A.fa.gz'}\n")
    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")]},
        extra=("--species-fasta-map", str(mapping)))
    assert result.returncode != 0 and "unknown species FASTA mapping" in result.stderr


@pytest.mark.parametrize("omission", ["missing_cell", "missing_marker"])
def test_summary_cannot_silently_omit_full_table_candidates(tmp_path, omission):
    def omit(root):
        summary = root / "summary.tsv"
        text = summary.read_text()
        summary.write_text(text.replace("g1,g2", "-") if omission == "missing_cell" else text.splitlines()[0] + "\n")

    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")]}, mutate=omit)
    assert result.returncode != 0 and "summary/full-table ID mismatch" in result.stderr


@pytest.mark.parametrize("collision", ["summary", "fasta", "full", "hardlink", "counts", "marker", "input_directory"])
def test_output_overlap_is_rejected_before_any_input_or_existing_result_is_changed(tmp_path, collision):
    protected = {}

    def prepare(root):
        if collision == "hardlink":
            os.link(root / "full/A.busco.full.tsv", root / "alias.tsv")
        if collision == "input_directory":
            raw = root / "A.fa"
            raw.write_bytes(gzip.decompress((root / "A.fa.gz").read_bytes()))
            (root / "fastas.txt").write_text(f"{raw}\n")
            (root / "species.tsv").write_text(f"species\tfasta\nA\t{raw}\n")
        existing = root / "output/preserved.fa"
        existing.parent.mkdir()
        existing.write_text(">previous\nMAAA\n")
        protected.update({path: path.read_bytes() for path in root.rglob("*") if path.is_file()})

    target = {"summary": tmp_path / "summary.tsv", "fasta": tmp_path / "A.fa.gz",
              "full": tmp_path / "full/A.busco.full.tsv", "hardlink": tmp_path / "alias.tsv",
              "counts": tmp_path / "counts.tsv", "marker": tmp_path / "output/b1.fa"}
    extra = ("--output-dir", str(tmp_path)) if collision == "input_directory" else (
        "--selection-report", str(target[collision]))
    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")]},
        extra=extra, mutate=prepare)
    assert result.returncode != 0 and "Output paths overlap" in result.stderr
    assert all(path.is_file() and path.read_bytes() == value for path, value in protected.items())


@pytest.mark.parametrize("table,text", [
    ("species.tsv", "species\tfasta\nA\n"),
    ("species.tsv", "species\tfasta\nA\tanywhere\tunexpected\n"),
    ("codes.tsv", "species\tgenetic_code\nA\n"),
    ("codes.tsv", "species\tgenetic_code\nA\t\n"),
])
def test_malformed_mapping_tables_fail_with_context_instead_of_a_traceback(tmp_path, table, text):
    def malformed(root):
        (root / table).write_text(text)

    extra = ("--species-genetic-codes", str(tmp_path / "codes.tsv")) if table == "codes.tsv" else ()
    result, _ = run_extraction(tmp_path, {"A": [("g1", "ATGAAA"), ("g2", "ATGCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "2"), ("b1", "Duplicated", "g2", "90", "2")]},
        mode="cds", extra=extra, mutate=malformed)
    assert result.returncode != 0 and "Traceback" not in result.stderr
    assert "row" in result.stderr.lower()


@pytest.mark.parametrize("text", ["> \nMAAA\n", ">g1\n>g2\nMCCC\n"])
def test_malformed_candidate_fasta_fails_with_context(tmp_path, text):
    def malformed(root):
        with gzip.open(root / "A.fa.gz", "wt") as handle:
            handle.write(text)

    result, _ = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Complete", "g1", "100", "4"), ("b2", "Complete", "g2", "90", "4")]}, mutate=malformed)
    assert result.returncode != 0 and "Traceback" not in result.stderr
    assert "FASTA" in result.stderr


def test_explicit_species_headers_do_not_infer_species_from_arbitrary_or_reused_gene_ids(tmp_path):
    result, audit = run_extraction(tmp_path,
        {"A": [("g1", "MAAA"), ("g2", "MCCC")], "B": [("g1", "MDDD"), ("g2", "MEEE")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4")],
         "B": [("b1", "Duplicated", "g1", "90", "4"), ("b1", "Duplicated", "g2", "100", "4")]},
        extra=("--species-headers",))
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output/b1.fa").read_text() == ">A\nMAAA\n>B\nMEEE\n"
    assert {(row["species"], row["gene_id"]) for row in audit if row["selected"] == "1"} == {("A", "g1"), ("B", "g2")}


def test_missing_only_markers_can_be_absent_from_shared_summary(tmp_path):
    def remove_empty_marker(root):
        summary = root / "summary.tsv"
        summary.write_text("\n".join(line for line in summary.read_text().splitlines() if not line.startswith("missing\t")) + "\n")

    result, audit = run_extraction(tmp_path, {"A": [("g1", "MAAA"), ("g2", "MCCC")]},
        {"A": [("b1", "Duplicated", "g1", "100", "4"), ("b1", "Duplicated", "g2", "90", "4"),
               ("missing", "Missing", "-", "-", "-")]}, mutate=remove_empty_marker)
    assert result.returncode == 0, result.stderr
    assert len(audit) == 2 and audit[0]["selected"] == "1"


def test_literal_numeric_and_na_like_gene_ids_are_preserved(tmp_path):
    result, audit = run_extraction(tmp_path, {"A": [("001", "MAAA"), ("NA", "MCCC"), ("nan", "MDDD")]},
        {"A": [("b1", "Duplicated", gene, str(score), "4") for gene, score in [("001", 100), ("NA", 90), ("nan", 80)]]})
    assert result.returncode == 0, result.stderr
    assert {row["gene_id"] for row in audit} == {"001", "NA", "nan"}
    assert (tmp_path / "output/b1.fa").read_text().startswith(">001 ")


def test_output_pool_reopen_keeps_nofollow_protection_after_path_caching(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    from batch_extract_busco_fasta import OutputPool

    pool = OutputPool(tmp_path / "output", ".fa", 1)
    pool.output_dir.mkdir()
    protected = tmp_path / "protected.fa"
    protected.write_text("preserve\n")
    try:
        pool.write("b1", [">g1\n", "MAAA\n"])
        pool.write("b2", [">g2\n", "MCCC\n"])
        (pool.output_dir / "b1.fa").unlink()
        (pool.output_dir / "b1.fa").symlink_to(protected)
        with pytest.raises(OSError):
            pool.write("b1", [">g3\n", "MDDD\n"])
    finally:
        pool.close()
    assert protected.read_text() == "preserve\n"
