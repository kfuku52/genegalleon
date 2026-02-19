from pathlib import Path
import gzip
import subprocess
import sys

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "build_uniprot_metadata_from_dat.py"


def test_build_uniprot_metadata_from_dat_extracts_extended_fields(tmp_path):
    dat_gz = tmp_path / "uniprot_sprot.dat.gz"
    with gzip.open(dat_gz, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                "ID   TEST1_HUMAN              Reviewed;         100 AA.\n"
                "AC   P12345; Q99999;\n"
                "DE   RecName: Full=Protein kinase ABC;\n"
                "DE   EC=2.7.11.1;\n"
                "GN   Name=AKT1; Synonyms=PKB;\n"
                "OS   Homo sapiens (Human).\n"
                "OX   NCBI_TaxID=9606;\n"
                "CC   -!- SUBCELLULAR LOCATION: Cytoplasm.\n"
                "CC       Cell membrane {ECO:0000250}.\n"
                "KW   Kinase; Transferase.\n"
                "DR   GO; GO:0005524; F:ATP binding; IEA:InterPro.\n"
                "DR   GO; GO:0004672; P:protein kinase activity; IEA:UniProtKB-KW.\n"
                "DR   KEGG; hsa:1234; -.\n"
                "DR   KO; K00844; -.\n"
                "DR   InterPro; IPR000719; Protein kinase domain.\n"
                "DR   Pfam; PF00069; Pkinase.\n"
                "DR   Reactome; R-HSA-123456; Signal transduction.\n"
                "FT   SIGNAL          1..20\n"
                "FT   TRANSMEM        45..67\n"
                "FT   TRANSMEM        70..90\n"
                "//\n"
            )
        )

    outfile = tmp_path / "uniprot_sprot.meta.tsv.gz"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--uniprot_dat",
            str(dat_gz),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.shape[0] == 2
    assert sorted(out["accession"].tolist()) == ["P12345", "Q99999"]

    row = out.loc[out["accession"] == "P12345", :].iloc[0]
    assert row["signal_start"] == 1
    assert row["signal_end"] == 20
    assert row["transmem_count"] == 2
    assert row["transmem_aa"] == 44
    assert row["transmem_regions"] == "45-67;70-90"
    assert row["kegg_gene"] == "hsa:1234"
    assert row["kegg_orthology"] == "K00844"
    assert row["go_ids"] == "GO:0005524;GO:0004672"
    assert row["go_aspects"] == "F;P"
    assert row["go_terms"] == "ATP binding; protein kinase activity"
    assert row["go_evidence"] == "IEA"
    assert row["gene_name_primary"] == "AKT1"
    assert row["gene_name_synonyms"] == "PKB"
    assert row["ec_numbers"] == "2.7.11.1"
    assert row["subcellular_location"] == "Cytoplasm; Cell membrane"
    assert row["keywords"] == "Kinase;Transferase"
    assert row["interpro_ids"] == "IPR000719"
    assert row["pfam_ids"] == "PF00069"
    assert row["reactome_ids"] == "R-HSA-123456"
    assert row["organism"] == "Homo sapiens (Human)"
    assert row["taxid"] == 9606
