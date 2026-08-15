"""Command-line argument definitions for species input formatting."""

import argparse
import os
from pathlib import Path

from format_species_annotation.gff_repair import GFF_REPAIR_MODES
from format_species_constants import GENE_GROUPING_MODES
from format_species_provider_config import PROVIDERS


def default_output_path(*parts):
    root = Path(
        os.environ.get(
            "GG_INPUT_GENERATION_OUTPUT_ROOT",
            "workspace/output/input_generation",
        )
    )
    return str(root.joinpath(*parts))


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Format CDS/GFF inputs for genegalleon from provider-specific raw files. "
            "Accepts CDS-only inputs, GenBank/GBFF/EMBL flatfiles, or GFF plus genome FASTA inputs. "
            "When CDS is missing but both GFF and genome FASTA are present, CDS is derived automatically. "
            "Current targets are based on the established species_wise_original input layout."
        )
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + PROVIDERS,
        required=True,
        help="Input provider type. Use 'all' with --input-dir pointing to a provider-root directory.",
    )
    parser.add_argument(
        "--input-dir",
        default="",
        help=(
            "Provider input directory. For ensembl/ensemblplants/ensemblmetazoa/ensemblprotists: original_files/. "
            "For phycocosm/phytozome/ncbi/coge/cngb/gwh/citrusgenomedb/figshare/plantgarden/plantaedb/flybase/wormbase/vectorbase/fernbase/veupathdb/dictybase/insectbase/oryza_minuta/direct/local: species_wise_original/. "
            "Legacy aliases refseq/genbank are treated as ncbi. "
            "For --provider all, this must be the shared root containing all provider subdirectories."
        ),
    )
    parser.add_argument(
        "--species-cds-dir",
        default=default_output_path("species_cds"),
        help="Output directory for formatted species CDS FASTA files.",
    )
    parser.add_argument(
        "--species-gff-dir",
        default=default_output_path("species_gff"),
        help="Output directory for formatted species GFF files.",
    )
    parser.add_argument(
        "--species-genome-dir",
        default=default_output_path("species_genome"),
        help="Output directory for formatted species genome FASTA files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs.",
    )
    parser.add_argument(
        "--download-manifest",
        default="",
        help=(
            "Optional XLSX/TSV/CSV manifest for input download. "
            "Required columns: provider,id (in this order for XLSX templates). "
            "provider,id are required on every row; species_key is optional. "
            "CDS/GFF/genome can be set directly with cds_url/gff_url/genome_url or resolved from id "
            "(ncbi supports GCF/GCA/NCBI-URL auto-resolution; "
            "other supported providers support id-based template/index inference). "
            "Supported providers for --download-manifest: "
            "ensembl, ensemblplants, ensemblmetazoa, ensemblprotists, ncbi, coge, cngb, gwh, citrusgenomedb, figshare, plantgarden, plantaedb, flybase, wormbase, vectorbase, fernbase, veupathdb, dictybase, insectbase, oryza_minuta, direct, local "
            "(legacy aliases refseq/genbank are treated as ncbi). "
            "provider=citrusgenomedb accepts Citrus Genome Database organism or analysis URLs and resolves public genome/GFF/CDS downloads when available. "
            "provider=figshare accepts public article URLs or numeric article ids and can resolve file downloads from the figshare article API; use *_filename columns to disambiguate multi-file articles and *_archive_member columns for files inside .zip, .tar.*, or .rar archives. "
            "provider=plantgarden accepts PlantGARDEN species/genome/download URLs and resolves public genome/GFF/CDS-or-transcript downloads when available. "
            "provider=oryza_minuta downloads the public Gramene BB+CC tetraploid bundles and merges them into one species bundle. "
            "provider=direct requires explicit urls or an index-style id URL. "
            "cds_url may be provided by itself, or may be empty when both gff_url and genome_url are available; "
            "genegalleon will derive CDS during formatting. "
            "provider=local reads local files/directories (for example local phytozome files). "
            "Use --input-dir for direct local phycocosm/phytozome formatting. "
            "Optional columns: species_key,cds_archive_member,gff_archive_member,gbff_archive_member,genome_archive_member,cds_filename,gff_filename,genome_filename,"
            "cds_url_template,gff_url_template,genome_url_template,"
            "local_cds_path,local_gff_path,local_genome_path,"
            "fernbase_confidence_mode."
        ),
    )
    parser.add_argument(
        "--resolved-manifest-output",
        default=default_output_path("download_plan.resolved.tsv"),
        help=(
            "Output TSV path for resolved download-manifest rows. "
            "Rows selected for download are written with provider/species_key/URL/filename fields filled. "
            "Set empty string to disable."
        ),
    )
    parser.add_argument(
        "--download-dir",
        default=default_output_path("tmp", "input_download_cache"),
        help=("Directory for raw downloaded provider files. This can be reused as --input-dir for formatting."),
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Download raw files from manifest and exit without formatting.",
    )
    parser.add_argument(
        "--http-header",
        action="append",
        default=[],
        help="Additional HTTP header for download requests. Repeatable. Format: 'Key: Value'.",
    )
    parser.add_argument(
        "--auth-bearer-token-env",
        default="",
        help="Environment variable name containing bearer token for download Authorization header.",
    )
    parser.add_argument(
        "--auth-cookie-env",
        default="",
        help="Environment variable name containing raw Cookie header for download requests.",
    )
    parser.add_argument(
        "--jgi-login-env",
        default="",
        help="Environment variable name containing JGI login/email for Genome Portal sign-on.",
    )
    parser.add_argument(
        "--jgi-password-env",
        default="",
        help="Environment variable name containing JGI password for Genome Portal sign-on.",
    )
    parser.add_argument(
        "--download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) per download request.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=0,
        help=(
            "Maximum parallel download workers. "
            "If <=0, uses GG_TASK_CPUS when available, otherwise legacy NSLOTS, otherwise 1."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan actions without writing downloaded or formatted files.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Exit with error if any species lacks a usable source bundle or if any unexpected "
            "CDS-to-GFF mapping mismatch remains."
        ),
    )
    parser.add_argument(
        "--gene-grouping-mode",
        choices=GENE_GROUPING_MODES,
        default="rescue_overlap",
        help=(
            "How to group transcript/CDS records into gene-level representatives. "
            "'strict' uses annotation hierarchy and identifier normalization only. "
            "When both CDS and GFF are available, GFF transcript/gene aliases are used to group "
            "the provided CDS records. 'rescue_overlap' also merges highly similar overlapping "
            "GFF transcript models when annotation gene IDs appear inconsistent."
        ),
    )
    parser.add_argument(
        "--gff-repair-mode",
        choices=GFF_REPAIR_MODES,
        default="safe",
        help=(
            "How to repair formatted GFF gene identifiers against final CDS IDs. "
            "'safe' applies only unique, collision-free repairs; 'strict' also rejects "
            "ambiguous repair candidates; 'off' preserves the source GFF identifiers."
        ),
    )
    parser.add_argument(
        "--stats-output",
        default="",
        help=(
            "Optional JSON output path for run stats. "
            "Includes CDS counts before/after processing and first CDS sequence name."
        ),
    )
    parser.add_argument(
        "--species-summary-output",
        default=default_output_path("gg_input_generation_species.tsv"),
        help=(
            "Per-species TSV summary output path. "
            "Successful formatting results are persisted incrementally and retained across runs "
            "while referenced species output files exist."
        ),
    )
    return parser
