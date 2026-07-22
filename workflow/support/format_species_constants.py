"""Shared constants for species input formatting and downloads."""

import re
import threading
from collections import deque

from format_species_provider_config import ORYZA_MINUTA_PROVIDER

FERNBASE_CONFIDENCE_MODE_FIELD = "fernbase_confidence_mode"

FERNBASE_CONFIDENCE_MODE_HIGH_ONLY = "high-confidence only"

FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED = "high-low combined"

FERNBASE_CONFIDENCE_MODE_CHOICES = (
    FERNBASE_CONFIDENCE_MODE_HIGH_ONLY,
    FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED,
)

FERNBASE_COMBINED_FILENAME_MARKER = "highlowcombined"

ORYZA_MINUTA_DEFAULT_SPECIES_KEY = "Oryza_minuta"

ORYZA_MINUTA_CANONICAL_SOURCE_ID = "gramene_tetraploids"

ORYZA_MINUTA_SOURCE_ID_ALIASES = frozenset(
    (
        "",
        "oryza_minuta",
        "gramene",
        "tetraploids",
        "gramene_tetraploids",
        "oryza_minuta_gramene_tetraploids",
    )
)

ORYZA_MINUTA_SUBGENOME_BUNDLES = (
    {
        "label": "BB",
        "token": "oryza_minutabb",
        "cds_filename": "Oryza_minutabb.oryza_minutabb.cds.all.fa.gz",
        "gff_filename": "Oryza_minutabb.oryza_minutabb.gff3.gz",
        "genome_filename": "Oryza_minutabb.oryza_minutabb.dna.toplevel.fa.gz",
    },
    {
        "label": "CC",
        "token": "oryza_minutacc",
        "cds_filename": "Oryza_minutacc.oryza_minutacc.cds.all.fa.gz",
        "gff_filename": "Oryza_minutacc.oryza_minutacc.gff3.gz",
        "genome_filename": "Oryza_minutacc.oryza_minutacc.dna.toplevel.fa.gz",
    },
)

GENE_GROUPING_MODES = ("strict", "rescue_overlap")

RESCUE_SHARED_JUNCTION_MIN_SHORTER_OVERLAP = 0.70

RESCUE_SHARED_JUNCTION_MIN_LONGER_OVERLAP = 0.40

RESCUE_SAME_TERMINAL_MIN_SHORTER_OVERLAP = 0.90

RESCUE_SAME_TERMINAL_MIN_LONGER_OVERLAP = 0.60

FASTA_EXTENSIONS = (
    ".cds",
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".cds.gz",
    ".fa.bz2",
    ".fas.bz2",
    ".fasta.bz2",
    ".fna.bz2",
    ".fa.gz",
    ".fas.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fa.tar.gz",
    ".fas.tar.gz",
    ".fasta.tar.gz",
    ".fna.tar.gz",
    ".fa.tar.bz2",
    ".fas.tar.bz2",
    ".fasta.tar.bz2",
    ".fna.tar.bz2",
)

FASTA_ARCHIVE_EXTENSIONS = (
    ".fa.tar.gz",
    ".fas.tar.gz",
    ".fasta.tar.gz",
    ".fna.tar.gz",
    ".fa.tar.bz2",
    ".fas.tar.bz2",
    ".fasta.tar.bz2",
    ".fna.tar.bz2",
)

GFF_EXTENSIONS = (
    ".gff",
    ".gff3",
    ".gtf",
    ".gff.gz",
    ".gff3.gz",
    ".gtf.gz",
)

GENBANK_EXTENSIONS = (
    ".gb",
    ".gbk",
    ".gbff",
    ".genbank",
    ".gb.gz",
    ".gbk.gz",
    ".gbff.gz",
    ".genbank.gz",
)

SPECIES_KEY_FILENAME_SUFFIXES = (
    (
        ".derived.cds.fa.gz",
        "_derived.cds.fa.gz",
        ".derived.genome.fa.gz",
        "_derived.genome.fa.gz",
        ".derived.gff.gz",
        "_derived.gff.gz",
        ".cds.all.fa.gz",
        "_cds.all.fa.gz",
        ".cds.fa.gz",
        "_cds.fa.gz",
        ".cds.fna.gz",
        "_cds.fna.gz",
        ".genome.fa.gz",
        "_genome.fa.gz",
        ".genomic.fna.gz",
        "_genomic.fna.gz",
        ".dna.primary_assembly.fa.gz",
        ".dna.toplevel.fa.gz",
        ".pep.fa.gz",
        "_pep.fa.gz",
        ".protein.fa.gz",
        "_protein.fa.gz",
    )
    + FASTA_EXTENSIONS
    + GFF_EXTENSIONS
    + GENBANK_EXTENSIONS
)

ENSEMBL_CDS_PATTERN = re.compile(r"(?:^|[._-])cds(?:[._-]|$)", re.IGNORECASE)

INVALID_ID_CHARS = re.compile(r"[%/\+:;&\^\$#@!~=\'\"`*\(\)\{\}\[\]\|\?\s]+")

ENSEMBL_GENE_ID_PATTERN = re.compile(r"^ENS[A-Z0-9]*G[0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)

COGE_ID_HINT_PATTERN = re.compile(r"^coge[:_].+", re.IGNORECASE)

CNGB_ID_HINT_PATTERN = re.compile(r"^cngb[:_].+", re.IGNORECASE)

DDBJ_ID_HINT_PATTERN = re.compile(r"^ddbj[:_].+", re.IGNORECASE)

GWH_ID_HINT_PATTERN = re.compile(r"^gwh[:_].+", re.IGNORECASE)

CITRUSGENOMEDB_ID_HINT_PATTERN = re.compile(r"^citrusgenomedb[:_].+", re.IGNORECASE)

FIGSHARE_ID_HINT_PATTERN = re.compile(r"^figshare[:_].+", re.IGNORECASE)

PLANTGARDEN_ID_HINT_PATTERN = re.compile(r"^plantgarden[:_].+", re.IGNORECASE)

ENSEMBL_ID_HINT_PATTERN = re.compile(r"^ensembl[:_].+", re.IGNORECASE)

ENSEMBLMETAZOA_ID_HINT_PATTERN = re.compile(r"^ensemblmetazoa[:_].+", re.IGNORECASE)

ENSEMBLPROTISTS_ID_HINT_PATTERN = re.compile(r"^ensemblprotists[:_].+", re.IGNORECASE)

FERNBASE_ID_HINT_PATTERN = re.compile(r"^fernbase[:_].+", re.IGNORECASE)

VEUPATHDB_ID_HINT_PATTERN = re.compile(r"^veupathdb[:_].+", re.IGNORECASE)

DICTYBASE_ID_HINT_PATTERN = re.compile(r"^dictybase[:_].+", re.IGNORECASE)

INSECTBASE_ID_HINT_PATTERN = re.compile(r"^insectbase[:_].+", re.IGNORECASE)

PLANTAEDB_ID_HINT_PATTERN = re.compile(r"^plantaedb[:_].+", re.IGNORECASE)

ORYZA_MINUTA_ID_HINT_PATTERN = re.compile(r"^oryza_minuta[:_].+", re.IGNORECASE)

INSECTBASE_IBG_ID_PATTERN = re.compile(r"^IBG_[0-9]+$", re.IGNORECASE)

COGE_GID_PATTERN = re.compile(r"^[0-9]+$")

CNGB_ASSEMBLY_ACCESSION_PATTERN = re.compile(r"^CNA[0-9]+$", re.IGNORECASE)

DDBJ_BIOPROJECT_PATTERN = re.compile(r"(PRJDB[0-9]+)", re.IGNORECASE)

DDBJ_WGS_MASTER_ACCESSION_PATTERN = re.compile(r"^([A-Z]{4,6})0{8,9}(?:\.[0-9]+)?$", re.IGNORECASE)

DDBJ_WGS_ACCESSION_PATTERN = re.compile(r"^([A-Z]{4,6})[0-9]{8,9}(?:\.[0-9]+)?$", re.IGNORECASE)

DDBJ_WGS_PROJECT_PREFIX_PATTERN = re.compile(r"^[A-Z]{4,6}$", re.IGNORECASE)

GWH_ASSEMBLY_ACCESSION_PATTERN = re.compile(r"^GWH[A-Z0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)

GWH_ASSEMBLY_ACCESSION_SEARCH_PATTERN = re.compile(r"(GWH[A-Z0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)

CITRUSGENOMEDB_BINOMIAL_PATTERN = re.compile(r"\b([A-Z][a-z]+)\s+((?:x\s+)?[a-z][a-z-]+)\b")

FIGSHARE_ARTICLE_ID_PATTERN = re.compile(r"/articles(?:/[^/?#]+)*/([0-9]+)(?:$|[/?#])", re.IGNORECASE)

PLANTGARDEN_TAXID_PATTERN = re.compile(r"^t[0-9]+$", re.IGNORECASE)

PLANTGARDEN_GENOME_ID_PATTERN = re.compile(r"^(t[0-9]+)[.]G[0-9]+$", re.IGNORECASE)

DEFAULT_DOWNLOAD_LOCK_STALE_SECONDS = 900

DEFAULT_DOWNLOAD_LOCK_HEARTBEAT_SECONDS = 60

DEFAULT_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS = 3600

DEFAULT_DOWNLOAD_LOCK_POLL_SECONDS = 5.0

DEFAULT_DOWNLOAD_ATTEMPTS = 4

DEFAULT_DOWNLOAD_RETRY_BASE_SECONDS = 5.0

DOWNLOAD_DIAGNOSTIC_KEYS = (
    "cache_preexisting_partial_tmp",
    "cache_preexisting_corrupt",
    "cache_preexisting_locks",
    "cache_final_partial_tmp",
    "cache_final_corrupt",
    "cache_final_locks",
    "download_jobs",
    "failed_downloads",
    "transient_retries",
    "range_resumes",
    "corrupt_download_retries",
    "corrupt_cache_recoveries",
    "stale_locks_recovered",
    "lock_waits",
)

TRANSIENT_HTTP_STATUS_CODES = frozenset((408, 425, 429, 500, 502, 503, 504))

NCBI_DATASETS_INCLUDE_BY_LABEL = {
    "CDS": "CDS_FASTA",
    "GFF": "GENOME_GFF",
    "GBFF": "GENOME_GBFF",
    "GENOME": "GENOME_FASTA",
}

DOWNLOAD_LABELS = ("CDS", "GFF", "GENOME")

TAXONOMIC_PROXIMITY_QUALIFIERS = frozenset(("cf", "aff", "nr"))

TAXONOMIC_GENUS_ONLY_PLACEHOLDERS = frozenset(("sp", "spp"))

TAXONOMIC_HYBRID_CONNECTORS = frozenset(("x", "\u00d7", "hybrid"))

TAXONOMIC_INFRASPECIFIC_RANK_ALIASES = {
    "subsp": "subsp",
    "ssp": "subsp",
    "subspecies": "subsp",
    "var": "var",
    "variety": "var",
    "forma": "forma",
    "form": "forma",
    "f": "forma",
    "strain": "strain",
    "substrain": "substrain",
    "serovar": "serovar",
    "serotype": "serotype",
    "serogroup": "serogroup",
    "pathovar": "pathovar",
    "pv": "pathovar",
    "biovar": "biovar",
    "biotype": "biotype",
    "chemovar": "chemovar",
    "morphovar": "morphovar",
    "cultivar": "cultivar",
    "cv": "cultivar",
    "isolate": "isolate",
    "group": "group",
    "subgroup": "subgroup",
    "complex": "complex",
    "clade": "clade",
    "lineage": "lineage",
    "section": "section",
    "series": "series",
    "ecotype": "ecotype",
    "breed": "breed",
}

PROVIDER_DEFAULT_MAX_CONCURRENT_DOWNLOADS = {
    # NCBI E-utilities explicitly documents request-rate limits.
    # We keep concurrent file downloads conservative and separately enforce
    # E-utilities request throttling below.
    "ncbi": 2,
    "ddbj": 1,
    "refseq": 2,
    "genbank": 2,
    # No explicit numeric concurrency guidance was found for these sources;
    # keep defaults conservative to reduce ban/abuse risk.
    "ensembl": 2,
    "ensemblplants": 2,
    "ensemblmetazoa": 2,
    "ensemblprotists": 2,
    "phycocosm": 1,
    "phytozome": 1,
    "coge": 1,
    "cngb": 1,
    "gwh": 1,
    "citrusgenomedb": 2,
    "figshare": 2,
    "plantgarden": 2,
    "plantaedb": 2,
    "flybase": 2,
    "wormbase": 2,
    "vectorbase": 2,
    "fernbase": 2,
    "veupathdb": 1,
    "dictybase": 1,
    "insectbase": 2,
    ORYZA_MINUTA_PROVIDER: 2,
    "direct": 2,
    "local": 1,
}

DEFAULT_GLOBAL_DOWNLOAD_WORKERS = 1

DEFAULT_NCBI_EUTILS_RPS_NO_API_KEY = 3

DEFAULT_NCBI_EUTILS_RPS_WITH_API_KEY = 10

_ncbi_eutils_rate_lock = threading.Lock()

_ncbi_eutils_request_times = deque()

SPECIES_SUMMARY_COLUMNS = (
    "updated_utc",
    "run_started_utc",
    "provider",
    "species_key",
    "species_prefix",
    "taxid",
    "nuclear_genetic_code_id",
    "nuclear_genetic_code_name",
    "mitochondrial_genetic_code_id",
    "mitochondrial_genetic_code_name",
    "plastid_genetic_code_id",
    "plastid_genetic_code_name",
    "cds_input_path",
    "gff_input_path",
    "genome_input_path",
    "cds_output_path",
    "gff_output_path",
    "genome_output_path",
    "cds_status",
    "gff_status",
    "gff_repair_mode",
    "gff_repair_status",
    "gff_repair_audit_path",
    "gff_repaired_gene_ids",
    "gff_repaired_references",
    "gff_repair_ambiguous",
    "gff_repair_collisions",
    "genome_status",
    "cds_sequences_before",
    "cds_sequences_after",
    "cds_first_sequence_name",
    "aggregated_cds_removed",
    "gene_grouping_mode",
    "overwrite",
    "dry_run",
)

PLASTID_GENETIC_CODE_LINEAGE_DEFAULTS = {
    # NCBI taxonomy/taxdump exposes nuclear and mitochondrial codes directly, but
    # not plastid codes. Use the standard plastid code for well-established
    # plastid-bearing clades as a conservative best-effort fallback.
    "33090": "11",  # Viridiplantae
    "2763": "11",  # Rhodophyta
    "2830": "11",  # Haptophyta
    "3027": "11",  # Cryptophyceae
    "2696291": "11",  # Ochrophyta
}
