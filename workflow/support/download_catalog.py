"""Shared download catalog schema and lossless manifest/snapshot readers."""

import csv
import json
import re
from pathlib import Path

from format_species_manifest import read_xlsx_table

FERNBASE_CONFIDENCE_MODE_FIELD = "fernbase_confidence_mode"


MANIFEST_FIELDNAMES = (
    "provider",
    "id",
    "species_key",
    "cds_url",
    "gff_url",
    "gbff_url",
    "genome_url",
    "cds_archive_member",
    "gff_archive_member",
    "gbff_archive_member",
    "genome_archive_member",
    "cds_filename",
    "gff_filename",
    "gbff_filename",
    "genome_filename",
    "cds_url_template",
    "gff_url_template",
    "genome_url_template",
    "local_cds_path",
    "local_gff_path",
    "local_gbff_path",
    "local_genome_path",
    FERNBASE_CONFIDENCE_MODE_FIELD,
)


DIRECT_CATALOG_FIELDS = tuple(field for field in MANIFEST_FIELDNAMES if field not in ("provider", "id"))


ID_EXAMPLES_BY_PROVIDER = {
    "onekp": (("BMJR", "Adiantum raddianum (1KP BMJR)"),),
    "ensembl": (("homo_sapiens", "Homo sapiens"), ("mus_musculus", "Mus musculus")),
    "ensemblplants": (("Ostreococcus_lucimarinus", "Ostreococcus lucimarinus"),),
    "ensemblmetazoa": (("anopheles_gambiae", "Anopheles gambiae"),),
    "ensemblprotists": (("phytophthora_parasitica", "Phytophthora parasitica"),),
    "ncbi": (
        ("GCF_000001405.40", "Homo sapiens"),
        ("GCA_000001635.9", "Mus musculus"),
        ("GCF_049306965.1", "Danio rerio"),
        ("GCA_000001215.4", "Drosophila melanogaster"),
        ("GCF_000002985.6", "Caenorhabditis elegans"),
    ),
    "ddbj": (
        ("PRJDB15739", "Triphyophyllum peltatum"),
        ("BAAHMP000000000", "Triphyophyllum peltatum WGS master"),
    ),
    "coge": (
        ("24739", "Arabidopsis thaliana"),
        ("29177", "Oryza sativa"),
        ("42091", "Zea mays"),
        ("78085", "Drosophila melanogaster"),
        ("72417", "Caenorhabditis elegans"),
    ),
    "cngb": (
        ("CNA0012345", "Homo sapiens"),
        ("GCF_000001405.40", "Homo sapiens"),
        ("GCA_000001635.9", "Mus musculus"),
        ("GCF_049306965.1", "Danio rerio"),
        ("GCA_000001215.4", "Drosophila melanogaster"),
    ),
    "gwh": (
        ("GWHIGRM00000000.1", "Medicago sativa"),
        ("GWHCBHY00000000", "Allium sativum"),
    ),
    "citrusgenomedb": (
        ("https://www.citrusgenomedb.org/organism/5799", "Citrus australasica (CGD organism page)"),
        ("https://www.citrusgenomedb.org/Analysis/2530647", "Citrus australasica cv. AZM genome v1.0 (CGD analysis)"),
    ),
    "figshare": (
        (
            "https://figshare.com/articles/dataset/Construction_of_the_super_pan-genome_for_the_genus_Actinidia_reveals_structural_variations_linked_to_phenotypic_diversity/28759280",
            "Actinidia super pan-genome (figshare article)",
        ),
    ),
    "plantgarden": (
        ("https://plantgarden.jp/en/list/t64480/genome/t64480.G001", "Actinidia polygama (PlantGARDEN assembly page)"),
        ("https://plantgarden.jp/en/list/t385388/genome/t385388.G001", "Camellia oleifera (PlantGARDEN assembly page)"),
    ),
    "plantaedb": (
        (
            "https://plantaedb.com/taxa/phylum/angiosperms/order/asterales/family/asteraceae/subfamily/asteroideae/tribe/astereae/subtribe/conyzinae/genus/erigeron/species/erigeron-breviscapus",
            "Erigeron breviscapus (PlantaeDB page)",
        ),
        (
            "https://plantaedb.com/taxa/phylum/angiosperms/order/ranunculales/family/berberidaceae/genus/berberis/species/berberis-thunbergii",
            "Berberis thunbergii (PlantaeDB page)",
        ),
    ),
    "flybase": (("dmel_r6.61", "Drosophila melanogaster"),),
    "wormbase": (("celegans_prjna13758_ws290", "Caenorhabditis elegans"),),
    "vectorbase": (("anopheles_gambiae_pest", "Anopheles gambiae"),),
    "fernbase": (("Azolla_filiculoides", "Azolla filiculoides"), ("Salvinia_cucullata_v2", "Salvinia cucullata v2")),
    "veupathdb": (("EnuttalliP19", "Entamoeba nuttalli"),),
    "dictybase": (("Dictyostelium_discoideum", "Dictyostelium discoideum"),),
    "insectbase": (("IBG_00001", "Abrostola tripartita"),),
    "oryza_minuta": (("gramene_tetraploids", "Oryza minuta Gramene tetraploids"),),
    "direct": (("direct_example_species", "Direct URL manifest row"),),
    "local": (("/absolute/path/to/local/species_dir", "Local species directory"),),
}


def species_label_from_species_key(species_key: str) -> str:
    text = str(species_key or "").strip()
    if text == "":
        return ""
    parts = text.replace("_", " ").split()
    normalized = []
    for idx, part in enumerate(parts):
        if re.fullmatch(r"[a-z]+", part):
            if idx == 0:
                normalized.append(part.capitalize())
            else:
                normalized.append(part.lower())
        else:
            normalized.append(part)
    return " ".join(normalized)


def dedupe_preserve_order(items):
    out = []
    seen = set()
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def parse_snapshot_entries(raw_entries):
    out = []
    entries = raw_entries
    if isinstance(entries, dict):
        entries = [entries]
    if not isinstance(entries, list):
        return out
    for entry in entries:
        source_id = ""
        species_label = ""
        if isinstance(entry, str):
            source_id = entry
        elif isinstance(entry, dict):
            source_id = str(entry.get("id") or entry.get("source_id") or "").strip()
            species_label = str(entry.get("species") or entry.get("species_label") or entry.get("label") or "").strip()
        elif isinstance(entry, (tuple, list)) and len(entry) > 0:
            source_id = str(entry[0] or "").strip()
            if len(entry) > 1:
                species_label = str(entry[1] or "").strip()
        source_id = str(source_id or "").strip()
        if source_id == "":
            continue
        out.append((source_id, species_label))
    return dedupe_preserve_order(out)


def parse_direct_catalog_entries(raw_entries):
    out = []
    entries = raw_entries
    if isinstance(entries, dict):
        entries = [entries]
    if not isinstance(entries, list):
        return out
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        source_id = str(entry.get("id") or entry.get("source_id") or "").strip()
        species_label = str(entry.get("species") or entry.get("species_label") or entry.get("label") or "").strip()
        if source_id == "":
            continue
        key = (source_id, str(entry.get("species_key") or "").strip(), species_label)
        if key in seen:
            continue
        seen.add(key)
        record = {"id": source_id, "species": species_label}
        for field in DIRECT_CATALOG_FIELDS:
            record[field] = str(entry.get(field) or "").strip()
        out.append(record)
    return out


def load_snapshot_payload(path: Path):
    with open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    providers_blob = payload.get("providers", payload)
    if not isinstance(providers_blob, dict):
        raise ValueError("id options snapshot must contain object key 'providers'")
    return providers_blob


def load_manifest_rows(path):
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        try:
            from openpyxl import load_workbook
        except ImportError:
            raise RuntimeError("openpyxl is required to read .xlsx direct catalog manifests") from None
        workbook = load_workbook(path, read_only=True, data_only=False)
        try:
            return read_xlsx_table(workbook.active)[1]
        finally:
            workbook.close()
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def build_direct_catalog_from_manifest_rows(rows):
    out = []
    seen = set()
    for row in rows:
        provider = str(row.get("provider", "") or "").strip().lower()
        if provider != "direct":
            continue
        source_id = str(row.get("id", "") or "").strip()
        if source_id == "":
            continue
        species_key = str(row.get("species_key", "") or "").strip()
        species_label = species_label_from_species_key(species_key)
        key = (source_id, species_key, species_label)
        if key in seen:
            continue
        seen.add(key)
        record = {"id": source_id, "species": species_label}
        for field in DIRECT_CATALOG_FIELDS:
            record[field] = str(row.get(field, "") or "").strip()
        out.append(record)
    return sorted(out, key=lambda record: (str(record.get("id", "")).lower(), str(record.get("species_key", "")).lower()))
