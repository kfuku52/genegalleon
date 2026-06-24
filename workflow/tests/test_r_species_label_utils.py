from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_PATH = REPO_ROOT / "workflow" / "support" / "species_label_utils.r"


def test_r_species_label_utils_preserves_taxonomic_qualifiers(tmp_path):
    test_script = tmp_path / "test_species_label_utils.R"
    test_script.write_text(
        f"""
source({str(UTILS_PATH)!r})

labels <- c(
  "Asimitellaria_furusei_var._furusei_longestCDS.fx2tab_cds.tsv",
  "Asimitellaria_furusei_var._subramosa_fx2tab_cds.tsv",
  "Asimitellaria_stylosa_var._stylosa_annotation.tsv",
  "Arisaema_sp._aooni_tpm.tsv",
  "/workspace/output/transcriptome_assembly/assembled_transcripts_with_isoforms/Arisaema_iyoanum_isoform.fa.gz",
  "/workspace/output/transcriptome_assembly/longest_cds/Asimitellaria_furusei_var._furusei_longestCDS.fa.gz",
  "Dictyostelium_cf._discoideum.longestCDS.fa.gz",
  "Asimitellaria_furusei_var._furusei.cds.fa.gz",
  "homo_sapiens.GRCh38.cds.all.fa.gz",
  "Arisaema_sp._aooni.dna.primary_assembly.fa.gz"
)
expected_keys <- c(
  "Asimitellaria_furusei_var._furusei",
  "Asimitellaria_furusei_var._subramosa",
  "Asimitellaria_stylosa_var._stylosa",
  "Arisaema_sp._aooni",
  "Arisaema_iyoanum",
  "Asimitellaria_furusei_var._furusei",
  "Dictyostelium_cf._discoideum",
  "Asimitellaria_furusei_var._furusei",
  "homo_sapiens",
  "Arisaema_sp._aooni"
)
actual_keys <- gg_species_label_from_filename(labels)
if (!identical(actual_keys, expected_keys)) {{
  stop(sprintf("Unexpected keys: %s", paste(actual_keys, collapse = ",")))
}}

expected_display <- c(
  "Asimitellaria furusei var. furusei",
  "Asimitellaria furusei var. subramosa",
  "Asimitellaria stylosa var. stylosa",
  "Arisaema sp. aooni"
)
actual_display <- gg_species_display_from_key(expected_keys[1:4])
if (!identical(actual_display, expected_display)) {{
  stop(sprintf("Unexpected display labels: %s", paste(actual_display, collapse = ",")))
}}

duplicate_message <- tryCatch(
  {{
    gg_stop_on_duplicate_species_keys(
      c("Asimitellaria_furusei_var._furusei", "Asimitellaria_furusei_var._furusei"),
      c("furusei_a.tsv", "furusei_b.tsv"),
      "unit-test files"
    )
    ""
  }},
  error = function(err) conditionMessage(err)
)
if (!grepl("unit-test files", duplicate_message, fixed = TRUE) ||
    !grepl("furusei_a.tsv", duplicate_message, fixed = TRUE) ||
    !grepl("furusei_b.tsv", duplicate_message, fixed = TRUE)) {{
  stop(sprintf("Duplicate guard did not report colliding files: %s", duplicate_message))
}}
""",
        encoding="utf-8",
    )

    completed = subprocess.run(
        ["Rscript", "--vanilla", str(test_script)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
