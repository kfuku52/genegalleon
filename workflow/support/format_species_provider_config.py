#!/usr/bin/env python3
"""Provider registry metadata for format_species_inputs.py."""

from pathlib import Path

ORYZA_MINUTA_PROVIDER = "oryza_minuta"

DEFAULT_INPUT_RELATIVE_DIRS = {
    "ensembl": Path("Ensembl") / "original_files",
    "ensemblplants": Path("20230216_EnsemblPlants") / "original_files",
    "ensemblmetazoa": Path("EnsemblMetazoa") / "original_files",
    "ensemblprotists": Path("EnsemblProtists") / "original_files",
    "phycocosm": Path("PhycoCosm") / "species_wise_original",
    "phytozome": Path("Phytozome") / "species_wise_original",
    "ncbi": Path("NCBI_Genome") / "species_wise_original",
    "ddbj": Path("DDBJ") / "species_wise_original",
    "refseq": Path("NCBI_RefSeq") / "species_wise_original",
    "genbank": Path("NCBI_GenBank") / "species_wise_original",
    "coge": Path("CoGe") / "species_wise_original",
    "cngb": Path("CNGB") / "species_wise_original",
    "gwh": Path("GWH") / "species_wise_original",
    "citrusgenomedb": Path("CitrusGenomeDB") / "species_wise_original",
    "figshare": Path("Figshare") / "species_wise_original",
    "plantgarden": Path("PlantGARDEN") / "species_wise_original",
    "plantaedb": Path("PlantaeDB") / "species_wise_original",
    "flybase": Path("FlyBase") / "species_wise_original",
    "wormbase": Path("WormBase") / "species_wise_original",
    "vectorbase": Path("VectorBase") / "species_wise_original",
    "fernbase": Path("FernBase") / "species_wise_original",
    "veupathdb": Path("VEuPathDB") / "species_wise_original",
    "dictybase": Path("dictyBase") / "species_wise_original",
    "insectbase": Path("InsectBase") / "species_wise_original",
    ORYZA_MINUTA_PROVIDER: Path("OryzaMinuta") / "species_wise_original",
    "direct": Path("Direct") / "species_wise_original",
    "local": Path("Local") / "species_wise_original",
}

PROVIDERS = (
    "ensembl",
    "ensemblplants",
    "ensemblmetazoa",
    "ensemblprotists",
    "phycocosm",
    "phytozome",
    "ncbi",
    "ddbj",
    "refseq",
    "genbank",
    "coge",
    "cngb",
    "gwh",
    "citrusgenomedb",
    "figshare",
    "plantgarden",
    "plantaedb",
    "flybase",
    "wormbase",
    "vectorbase",
    "fernbase",
    "veupathdb",
    "dictybase",
    "insectbase",
    ORYZA_MINUTA_PROVIDER,
    "direct",
    "local",
)

DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS = (
    "ensembl",
    "ensemblplants",
    "ensemblmetazoa",
    "ensemblprotists",
    "ncbi",
    "ddbj",
    "refseq",
    "genbank",
    "coge",
    "cngb",
    "gwh",
    "citrusgenomedb",
    "figshare",
    "plantgarden",
    "plantaedb",
    "flybase",
    "wormbase",
    "vectorbase",
    "fernbase",
    "veupathdb",
    "dictybase",
    "insectbase",
    ORYZA_MINUTA_PROVIDER,
    "direct",
    "local",
)

DOWNLOAD_MANIFEST_BUILDER_PROVIDERS = tuple(
    provider
    for provider in DOWNLOAD_MANIFEST_SUPPORTED_PROVIDERS
    if provider not in ("refseq", "genbank")
)
DOWNLOAD_MANIFEST_LEGACY_NCBI_PROVIDER_ALIASES = ("refseq", "genbank")

ID_OPTIONS_FETCH_PROVIDERS = (
    "ensembl",
    "ensemblplants",
    "ensemblmetazoa",
    "ensemblprotists",
    "flybase",
    "wormbase",
    "vectorbase",
    "fernbase",
    "veupathdb",
    "insectbase",
    "local",
)

ENSEMBL_LIKE_PROVIDERS = (
    "ensembl",
    "ensemblplants",
    "ensemblmetazoa",
    "ensemblprotists",
)
