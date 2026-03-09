from importlib.util import module_from_spec, spec_from_file_location
import io
from pathlib import Path
import csv
import gzip
import json
import shutil
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def load_module():
    spec = spec_from_file_location("format_species_inputs_module", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


class FakeTextPipe:
    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)
        return len(text)

    def close(self):
        return None

    def getvalue(self):
        return "".join(self.parts)


def test_format_species_inputs_with_small_fixture_all_providers(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    stats_json = tmp_path / "stats.json"

    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--stats-output",
        str(stats_json),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    ensembl_cds = out_cds / "Ostreococcus_lucimarinus_ASM9206v1.cds.all.fa.gz"
    phycocosm_cds = out_cds / "Microglena_spYARC_MicrYARC1_GeneCatalog_CDS_20220803.fa.gz"
    phytozome_cds = out_cds / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa.gz"
    assert ensembl_cds.exists()
    assert phycocosm_cds.exists()
    assert phytozome_cds.exists()

    with gzip.open(ensembl_cds, "rt", encoding="utf-8") as handle:
        ensembl_text = handle.read()
    assert ensembl_text.count(">Ostreococcus_lucimarinus_OSTLU_25062") == 1
    assert ">Ostreococcus_lucimarinus_OSTLU_99999" in ensembl_text
    assert "ATGAAN" in ensembl_text

    with gzip.open(phycocosm_cds, "rt", encoding="utf-8") as handle:
        phycocosm_text = handle.read()
    assert ">Microglena_spYARC_mRNA.MigICE15955" in phycocosm_text
    assert ">Microglena_spYARC_mRNA.MigICE_00468" in phycocosm_text

    with gzip.open(phytozome_cds, "rt", encoding="utf-8") as handle:
        phytozome_text = handle.read()
    assert ">Hydrocotyle_leucocephala_HyleuH1.06G006800" in phytozome_text
    assert "ATGATGAN" in phytozome_text

    ensembl_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff.gz"
    phytozome_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff.gz"
    phytozome_exons = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene_exons.gff.gz"
    assert ensembl_gff.exists()
    assert phytozome_gff.exists()
    assert not phytozome_exons.exists()

    with gzip.open(ensembl_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "evm.model." not in gff_text
    assert "Oropetium_20150105_" not in gff_text

    with gzip.open(phytozome_gff, "rt", encoding="utf-8") as handle:
        phytozome_gff_text = handle.read()
    assert "evm_27.model." not in phytozome_gff_text
    stats = json.loads(stats_json.read_text(encoding="utf-8"))
    assert stats["species_processed"] == 3
    assert stats["num_species_cds_files"] == 3
    assert stats["num_species_gff_files"] == 3
    assert stats["cds_sequences_before"] >= stats["cds_sequences_after"]
    assert stats["cds_first_sequence_name"] != ""


def test_format_species_inputs_strict_mode_fails_on_missing_pair(tmp_path):
    dataset_copy = tmp_path / "small_gfe_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = dataset_copy / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    missing_gff.unlink()

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(dataset_copy / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--strict",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "missing GFF" in completed.stderr


def test_format_species_inputs_fernbase_prefers_primary_annotation_files(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATG\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.transcript.highconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATGAAA\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    formatted_gff = out_gff / "Azolla_filiculoides_gene_models.highconfidence_v1.1.gff.gz"
    formatted_genome = out_genome / "Azolla_filiculoides_genome_v1.2.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Azolla_filiculoides_Azfi_g1" in cds_text
    assert "ATGAAN" in cds_text
    assert "lowconfidence" not in formatted_gff.name


def test_format_species_inputs_fernbase_accepts_markerless_genome_fasta(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Ceratopteris_richardii"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Crichardii_676_v2.0_cds.fa").write_text(">Crich_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Crichardii_676_v2.1.gene.gff3").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Crich_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Crichardii_676_v2.0.fa").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_genome = out_genome / "Ceratopteris_richardii_Crichardii_676_v2.0.fa.gz"
    assert formatted_genome.exists()


def test_format_species_inputs_fernbase_uses_plain_gene_tag_for_aggregation(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Pteridium_aquilinum"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "pt_aq.cds.fa").write_text(
        (
            ">pteridium_mrna-253 gene=pteridium_gene-252 seq_id=scaff1 type=cds\n"
            "ATGAA\n"
            ">pteridium_mrna-303 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAA\n"
            ">pteridium_mrna-308 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAATAG\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq.gff3").write_text(
        (
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=pteridium_gene-252\n"
            "chr1\tsrc\tgene\t20\t40\t.\t+\t.\tID=pteridium_gene-295\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq_final_genome.fasta").write_text(">chr1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Pteridium_aquilinum_pt_aq.cds.fa.gz"
    assert formatted_cds.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Pteridium_aquilinum_pteridium_gene-252",
        ">Pteridium_aquilinum_pteridium_gene-295",
    ]

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["species_prefix"] == "Pteridium_aquilinum"
    assert rows[0]["aggregated_cds_removed"] == "1"
    assert rows[0]["cds_sequences_before"] == "3"
    assert rows[0]["cds_sequences_after"] == "2"


def test_write_fasta_records_gzip_prefers_seqkit(monkeypatch, tmp_path):
    mod = load_module()
    output_path = tmp_path / "species.fa.gz"
    calls = {}

    class FakeSeqkitProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.command = command
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            output_arg = self.command[self.command.index("-o") + 1]
            with gzip.open(output_arg, "wt", encoding="utf-8") as handle:
                handle.write(self._stdin_pipe.getvalue())
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "3")
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(mod.subprocess, "Popen", FakeSeqkitProcess)

    mod.write_fasta_records_gzip(output_path, [("seq1", "ATG"), ("seq2", "ATGA")])

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">seq1" in text
    assert ">seq2" in text
    assert calls["command"][:4] == ["/usr/bin/seqkit", "seq", "--threads", "3"]
    assert calls["command"][-1] == "-"


def test_write_gff_gzip_prefers_pigz(monkeypatch, tmp_path):
    mod = load_module()
    input_path = tmp_path / "input.gff3"
    input_path.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=evm.model.g1\n", encoding="utf-8")
    output_path = tmp_path / "output.gff.gz"
    calls = {}

    class FakePigzProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.stdout = kwargs["stdout"]
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            with gzip.GzipFile(fileobj=self.stdout, mode="wb") as handle:
                handle.write(self._stdin_pipe.getvalue().encode("utf-8"))
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "4")
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(mod.subprocess, "Popen", FakePigzProcess)

    line_count = mod.write_gff_gzip(input_path, output_path)

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert line_count == 1
    assert "evm.model." not in text
    assert calls["command"] == ["/usr/bin/pigz", "-p", "4", "-c"]


def test_resolve_provider_download_limits_keeps_fernbase_default_cap_at_two(monkeypatch):
    mod = load_module()
    monkeypatch.delenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FERNBASE", raising=False)
    limits = mod.resolve_provider_download_limits(8)
    assert limits["fernbase"] == 2


def test_format_species_inputs_fernbase_prefers_namespaced_transcript_gene_id_when_gene_tag_is_short(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Salvinia_molesta"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "sg2.cds").write_text(
        (
            ">sg2.g13251.t1 gene=g13251\n"
            "ATGAAATAG\n"
            ">sg2.g9.t1 gene=g9\n"
            "ATGAAATAA\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "sg2.gff3").write_text(
        (
            "Chr_1\tgmst\tgene\t1\t9\t.\t+\t.\tID=sg2.g13251\n"
            "Chr_1\tgmst\tgene\t20\t28\t.\t+\t.\tID=sg2.g9\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "sg2_genome.fasta").write_text(">Chr_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Salvinia_molesta_sg2.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Salvinia_molesta_sg2.g13251",
        ">Salvinia_molesta_sg2.g9",
    ]


def test_format_species_inputs_fernbase_strips_amt_suffix_when_gff_uses_base_gene_id(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        (
            ">Azfi_s0034.g025227.AMT2\n"
            "ATGAAATAG\n"
            ">Azfi_s0093.g043301.AMT2\n"
            "ATGAAATAA\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        (
            "SCAF_1\tAUGUSTUS\tgene\t1\t9\t.\t+\t.\tID=Azfi_s0034.g025227\n"
            "SCAF_1\tAUGUSTUS\tgene\t20\t28\t.\t+\t.\tID=Azfi_s0093.g043301\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">SCAF_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Azolla_filiculoides_Azfi_s0034.g025227",
        ">Azolla_filiculoides_Azfi_s0093.g043301",
    ]


def test_species_summary_is_incremental_and_persistent_across_runs(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"

    first = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert first.returncode == 0, first.stderr + "\n" + first.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_first = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_first) == 1
    assert rows_first[0]["provider"] == "ensemblplants"
    assert rows_first[0]["species_prefix"] == "Ostreococcus_lucimarinus"
    assert int(rows_first[0]["cds_sequences_before"]) >= int(rows_first[0]["cds_sequences_after"])
    assert rows_first[0]["cds_first_sequence_name"] != ""

    second = run_script(
        "--provider",
        "phycocosm",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "PhycoCosm" / "species_wise_original"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert second.returncode == 0, second.stderr + "\n" + second.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_second = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_second) == 2
    by_key = {(row["provider"], row["species_prefix"]): row for row in rows_second}
    assert ("ensemblplants", "Ostreococcus_lucimarinus") in by_key
    assert ("phycocosm", "Microglena_spYARC") in by_key
    assert by_key[("phycocosm", "Microglena_spYARC")]["cds_first_sequence_name"] != ""
