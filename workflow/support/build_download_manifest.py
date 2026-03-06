#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path
import re
import sys
from typing import Dict, Iterable, List, Tuple

try:
    from openpyxl import Workbook
    from openpyxl.utils import get_column_letter
    from openpyxl.workbook.defined_name import DefinedName
    from openpyxl.styles import Font
    from openpyxl.comments import Comment
    from openpyxl.worksheet.datavalidation import DataValidation
except Exception:  # pragma: no cover - exercised in runtime environments without openpyxl
    Workbook = None
    get_column_letter = None
    DefinedName = None
    Font = None
    Comment = None
    DataValidation = None


FASTA_EXTENSIONS = (
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fas.gz",
    ".fasta.gz",
    ".fna.gz",
)

GFF_EXTENSIONS = (
    ".gff",
    ".gff3",
    ".gtf",
    ".gff.gz",
    ".gff3.gz",
    ".gtf.gz",
)

PROVIDERS = (
    "ensembl",
    "ensemblplants",
    "ncbi",
    "coge",
    "cngb",
    "flybase",
    "wormbase",
    "vectorbase",
    "fernbase",
    "local",
)
LEGACY_NCBI_PROVIDER_ALIASES = ("refseq", "genbank")

DEFAULT_INPUT_RELATIVE_DIRS = {
    "ensembl": Path("Ensembl") / "original_files",
    "ensemblplants": Path("20230216_EnsemblPlants") / "original_files",
    "ncbi": Path("NCBI_Genome") / "species_wise_original",
    "coge": Path("CoGe") / "species_wise_original",
    "cngb": Path("CNGB") / "species_wise_original",
    "flybase": Path("FlyBase") / "species_wise_original",
    "wormbase": Path("WormBase") / "species_wise_original",
    "vectorbase": Path("VectorBase") / "species_wise_original",
    "fernbase": Path("FernBase") / "species_wise_original",
    "local": Path("Local") / "species_wise_original",
}
NCBI_MERGED_INPUT_RELATIVE_DIRS = (
    Path("NCBI_Genome") / "species_wise_original",
    Path("NCBI_RefSeq") / "species_wise_original",
    Path("NCBI_GenBank") / "species_wise_original",
)

MANIFEST_FIELDNAMES = (
    "provider",
    "id",
    "species_key",
    "cds_url",
    "gff_url",
    "genome_url",
    "cds_filename",
    "gff_filename",
    "genome_filename",
    "cds_url_template",
    "gff_url_template",
    "genome_url_template",
    "local_cds_path",
    "local_gff_path",
    "local_genome_path",
)

LARGE_ID_PROVIDERS = ("ncbi",)
SNAPSHOT_FULL_ID_PROVIDERS = ("ensembl", "ensemblplants", "flybase", "wormbase", "vectorbase", "fernbase", "local")
EXAMPLE_ONLY_PROVIDERS = LARGE_ID_PROVIDERS + ("coge", "cngb")

ID_EXAMPLES_BY_PROVIDER = {
    "ensembl": (("homo_sapiens", "Homo sapiens"), ("mus_musculus", "Mus musculus")),
    "ensemblplants": (("Ostreococcus_lucimarinus", "Ostreococcus lucimarinus"),),
    "ncbi": (
        ("GCF_000001405.40", "Homo sapiens"),
        ("GCA_000001635.9", "Mus musculus"),
        ("GCF_049306965.1", "Danio rerio"),
        ("GCA_000001215.4", "Drosophila melanogaster"),
        ("GCF_000002985.6", "Caenorhabditis elegans"),
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
    "flybase": (("dmel_r6.61", "Drosophila melanogaster"),),
    "wormbase": (("celegans_prjna13758_ws290", "Caenorhabditis elegans"),),
    "vectorbase": (("anopheles_gambiae_pest", "Anopheles gambiae"),),
    "fernbase": (("Azolla_filiculoides", "Azolla filiculoides"), ("Salvinia_cucullata_v2", "Salvinia cucullata v2")),
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


def format_id_choice(source_id: str, species_label: str) -> str:
    sid = str(source_id or "").strip()
    label = str(species_label or "").strip()
    if sid == "":
        return ""
    if label == "":
        return sid
    return "{} ({})".format(sid, label)


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


def load_id_options_snapshot(path: Path):
    with open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    providers_blob = payload.get("providers", payload)
    if not isinstance(providers_blob, dict):
        raise ValueError("id options snapshot must contain object key 'providers'")
    options = {}
    for provider in PROVIDERS:
        options[provider] = parse_snapshot_entries(providers_blob.get(provider, []))
    return options


def infer_coge_genome_id_from_files(species_key: str, files: List[Path], warnings: List[str]) -> str:
    matches = set()
    key_match = re.search(r"(?:^|[._-])gid([0-9]+)(?:[._-]|$)", species_key, flags=re.IGNORECASE)
    if key_match is not None:
        matches.add(key_match.group(1))
    for path in files:
        name = path.name
        for match in re.finditer(r"(?:^|[._-])gid([0-9]+)(?:[._-]|$)", name, flags=re.IGNORECASE):
            matches.add(match.group(1))
        for match in re.finditer(r"dsgid[=_-]?([0-9]+)", name, flags=re.IGNORECASE):
            matches.add(match.group(1))
    if len(matches) == 0:
        return ""
    ordered = sorted(matches)
    if len(ordered) > 1:
        warnings.append(
            "[coge] {}: multiple genome_id candidates detected {}. Using '{}'".format(
                species_key, ",".join(ordered), ordered[0]
            )
        )
    return ordered[0]


def build_provider_id_options(rows: List[Dict[str, str]], snapshot_options=None):
    if snapshot_options is None:
        snapshot_options = {}
    discovered: Dict[str, Dict[str, str]] = {provider: {} for provider in PROVIDERS}
    for row in rows:
        provider = str(row.get("provider", "") or "").strip().lower()
        source_id = str(row.get("id", "") or "").strip()
        species_label = species_label_from_species_key(str(row.get("species_key", "") or ""))
        if provider not in discovered or source_id == "":
            continue
        if source_id not in discovered[provider]:
            discovered[provider][source_id] = species_label

    options = {}
    for provider in PROVIDERS:
        values = []
        if provider in EXAMPLE_ONLY_PROVIDERS:
            values = [
                format_id_choice(source_id, species_label)
                for source_id, species_label in ID_EXAMPLES_BY_PROVIDER.get(provider, ())
            ]
        else:
            snapshot_entries = list(snapshot_options.get(provider, []))
            if provider in SNAPSHOT_FULL_ID_PROVIDERS and len(snapshot_entries) > 0:
                if provider == "local":
                    values = [str(source_id).strip() for source_id, _ in snapshot_entries if str(source_id).strip()]
                else:
                    values = [
                        format_id_choice(source_id, species_label if species_label != "" else species_label_from_species_key(source_id))
                        for source_id, species_label in snapshot_entries
                    ]
            if len(values) == 0:
                if provider == "local":
                    values = [str(source_id).strip() for source_id in sorted(discovered[provider].keys()) if str(source_id).strip()]
                    if len(values) == 0:
                        values = [str(source_id).strip() for source_id, _species_label in ID_EXAMPLES_BY_PROVIDER.get(provider, ())]
                else:
                    values = [
                        format_id_choice(source_id, discovered[provider].get(source_id, ""))
                        for source_id in sorted(discovered[provider].keys())
                    ]
                    if len(values) == 0:
                        values = [
                            format_id_choice(source_id, species_label)
                            for source_id, species_label in ID_EXAMPLES_BY_PROVIDER.get(provider, ())
                        ]
        if len(values) == 0:
            values = [""]
        options[provider] = dedupe_preserve_order(values)
    return options


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Build a download manifest from an existing local provider-root directory. "
            "Generated URLs are file:// URLs for reuse in format_species_inputs.py --download-manifest."
        )
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + PROVIDERS + LEGACY_NCBI_PROVIDER_ALIASES,
        default="all",
        help="Provider to scan (default: all).",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Root directory containing provider subdirectories used by format_species_inputs.py.",
    )
    parser.add_argument(
        "--output",
        default="workspace/input/input_generation/download_plan.xlsx",
        help="Output path (.xlsx recommended; .tsv/.csv also supported).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error when a species is missing CDS or GFF.",
    )
    parser.add_argument(
        "--id-options-snapshot",
        default="",
        help=(
            "Optional JSON snapshot for provider-specific id dropdown values. "
            "When set, non-large providers (ensembl/ensemblplants/flybase/wormbase/vectorbase/fernbase/local) "
            "prefer snapshot IDs over locally discovered IDs."
        ),
    )
    return parser


def is_fasta_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in FASTA_EXTENSIONS)


def is_gff_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in GFF_EXTENSIONS)


def is_probable_genome_filename(provider, name):
    if not is_fasta_filename(name):
        return False
    lower = name.lower()
    non_genome_markers = (
        "cds",
        "mrna",
        "trna",
        "rrna",
        "ncrna",
        "rna_from_genomic",
        "transcript",
        "cdna",
        "pep",
        "protein",
    )
    if any(marker in lower for marker in non_genome_markers):
        return False
    if provider == "ncbi":
        return "genomic" in lower
    if provider in ("ensembl", "ensemblplants"):
        return any(marker in lower for marker in ("dna", "genome", "toplevel", "primary_assembly", "chromosome"))
    return any(marker in lower for marker in ("genome", "assembly", "genomic", "dna", "chromosome", "scaffold"))


FERNBASE_GFF_EXCLUDE_PATTERN = re.compile(
    r"(?:^|[._-])(te|teanno|repeat|repeats|transpos(?:on)?)(?:[._-]|$)",
    re.IGNORECASE,
)


def is_probable_cds_filename(provider, name):
    lower = name.lower()
    if provider == "fernbase" and (lower.endswith(".cds") or lower.endswith(".cds.gz")):
        return True
    if not is_fasta_filename(name):
        return False
    return any(marker in lower for marker in ("cds", "transcript", "mrna", "cdna"))


def provider_candidate_sort_key(provider, label, name):
    lower = str(name or "").lower()
    label_upper = str(label or "").upper()
    if provider != "fernbase":
        return (lower,)
    if label_upper == "CDS":
        return (
            0 if "highconfidence" in lower else 1,
            1 if "lowconfidence" in lower else 0,
            0 if "cds" in lower else 1,
            1 if any(marker in lower for marker in ("transcript", "mrna", "cdna")) else 0,
            lower,
        )
    if label_upper == "GFF":
        return (
            1 if FERNBASE_GFF_EXCLUDE_PATTERN.search(lower) else 0,
            0 if "highconfidence" in lower else 1,
            1 if "lowconfidence" in lower else 0,
            lower,
        )
    if label_upper == "GENOME":
        return (
            1 if "chloroplast" in lower else 0,
            0 if any(marker in lower for marker in ("genome", "assembly", "chr", "chromosome")) else 1,
            lower,
        )
    return (lower,)


def pick_single_file(matches, provider, species_key, label, warnings):
    if len(matches) == 0:
        return None
    ordered = sorted(matches, key=lambda path: provider_candidate_sort_key(provider, label, path.name))
    if len(ordered) > 1:
        warnings.append(
            "[{}] {}: multiple {} files found. Using '{}'".format(
                provider, species_key, label, ordered[0].name
            )
        )
    return ordered[0]


def discover_ensembl_like(input_dir, provider):
    warnings = []
    errors = []
    rows = []

    species_to_cds = {}
    species_to_gff = {}
    species_to_genome = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        species_key = path.name.split(".", 1)[0]
        if species_key == "":
            continue
        if ".cds." in path.name.lower() and is_fasta_filename(path.name):
            species_to_cds.setdefault(species_key, []).append(path)
        elif is_gff_filename(path.name):
            species_to_gff.setdefault(species_key, []).append(path)
        elif is_probable_genome_filename(provider, path.name):
            species_to_genome.setdefault(species_key, []).append(path)

    for species_key in sorted(set(species_to_cds.keys()) | set(species_to_gff.keys()) | set(species_to_genome.keys())):
        cds_path = pick_single_file(species_to_cds.get(species_key, []), provider, species_key, "CDS", warnings)
        gff_path = pick_single_file(species_to_gff.get(species_key, []), provider, species_key, "GFF", warnings)
        genome_path = pick_single_file(
            species_to_genome.get(species_key, []), provider, species_key, "genome", warnings
        )
        if cds_path is None or gff_path is None:
            errors.append(
                "[{}] {}: missing {}".format(
                    provider,
                    species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        rows.append(
            {
                "provider": provider,
                "id": species_key,
                "species_key": species_key,
                "cds_url": cds_path.resolve().as_uri(),
                "gff_url": gff_path.resolve().as_uri(),
                "genome_url": genome_path.resolve().as_uri() if genome_path is not None else "",
                "cds_filename": cds_path.name,
                "gff_filename": gff_path.name,
                "genome_filename": genome_path.name if genome_path is not None else "",
            }
        )

    return rows, warnings, errors


def discover_species_dir_based(provider, input_dir):
    warnings = []
    errors = []
    rows = []
    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = species_dir.name
        files = [path for path in species_dir.iterdir() if path.is_file()]
        if provider == "phycocosm":
            cds_matches = [path for path in files if "fasta" in path.name.lower() and is_fasta_filename(path.name)]
            gff_matches = [path for path in files if "gff" in path.name.lower() and is_gff_filename(path.name)]
            genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]
        elif provider == "phytozome":
            cds_matches = [
                path
                for path in files
                if ("cds_" in path.name.lower() or ".cds." in path.name.lower()) and is_fasta_filename(path.name)
            ]
            gff_matches = [path for path in files if "gene.gff3" in path.name.lower() and is_gff_filename(path.name)]
            genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]
        elif provider == "ncbi":
            cds_matches = [
                path
                for path in files
                if "cds_from_genomic" in path.name.lower() and is_fasta_filename(path.name)
            ]
            gff_matches = [
                path
                for path in files
                if "genomic.gff" in path.name.lower() and is_gff_filename(path.name)
            ]
            genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]
        else:
            cds_matches = [
                path
                for path in files
                if is_probable_cds_filename(provider, path.name)
            ]
            if len(cds_matches) == 0:
                cds_matches = [
                    path
                    for path in files
                    if is_fasta_filename(path.name) and not is_probable_genome_filename(provider, path.name)
                ]
            gff_matches = [path for path in files if is_gff_filename(path.name)]
            genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]

        cds_path = pick_single_file(cds_matches, provider, species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, provider, species_key, "GFF", warnings)
        genome_path = pick_single_file(genome_matches, provider, species_key, "genome", warnings)
        if cds_path is None or gff_path is None:
            errors.append(
                "[{}] {}: missing {}".format(
                    provider, species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        source_id = species_key
        if provider == "coge":
            source_id = infer_coge_genome_id_from_files(species_key, files, warnings)
            if source_id == "":
                errors.append(
                    "[coge] {}: genome_id (gid) was not detected from file names. "
                    "Include 'gid<NUM>' in at least one input filename.".format(species_key)
                )
                continue
        rows.append(
            {
                "provider": provider,
                "id": source_id,
                "species_key": species_key,
                "cds_url": cds_path.resolve().as_uri(),
                "gff_url": gff_path.resolve().as_uri(),
                "genome_url": genome_path.resolve().as_uri() if genome_path is not None else "",
                "cds_filename": cds_path.name,
                "gff_filename": gff_path.name,
                "genome_filename": genome_path.name if genome_path is not None else "",
            }
        )
    return rows, warnings, errors


def discover(provider, input_root, allow_missing=False):
    if provider in LEGACY_NCBI_PROVIDER_ALIASES:
        provider = "ncbi"
    if provider == "ncbi":
        warnings = []
        errors = []
        rows = []
        seen_ids = set()
        existing_dirs = []
        for rel_dir in NCBI_MERGED_INPUT_RELATIVE_DIRS:
            input_dir = input_root / rel_dir
            if not input_dir.exists() or not input_dir.is_dir():
                continue
            existing_dirs.append(input_dir)
            sub_rows, sub_warnings, sub_errors = discover_species_dir_based("ncbi", input_dir)
            warnings.extend(sub_warnings)
            errors.extend(sub_errors)
            for row in sub_rows:
                source_id = str(row.get("id", "") or "").strip()
                if source_id != "" and source_id in seen_ids:
                    warnings.append(
                        "[ncbi] duplicate accession '{}' across NCBI_Genome/NCBI_RefSeq/NCBI_GenBank. Keeping first occurrence.".format(
                            source_id
                        )
                    )
                    continue
                if source_id != "":
                    seen_ids.add(source_id)
                rows.append(row)
        if len(existing_dirs) == 0:
            message = "[ncbi] input directories not found: {}".format(
                ", ".join(str((input_root / rel_dir)) for rel_dir in NCBI_MERGED_INPUT_RELATIVE_DIRS)
            )
            if allow_missing:
                return [], [message], []
            return [], [], [message]
        return rows, warnings, errors

    input_dir = input_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]
    if not input_dir.exists() or not input_dir.is_dir():
        message = "[{}] input directory not found: {}".format(provider, input_dir)
        if allow_missing:
            return [], [message], []
        return [], [], [message]
    if provider in ("ensembl", "ensemblplants"):
        return discover_ensembl_like(input_dir, provider)
    return discover_species_dir_based(provider, input_dir)


def write_manifest_tsv(rows: Iterable[Dict[str, str]], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(MANIFEST_FIELDNAMES),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: str(row.get(key, "") or "") for key in MANIFEST_FIELDNAMES})


def write_manifest_xlsx(rows: List[Dict[str, str]], output_path: Path, id_options_snapshot=None):
    if (
        Workbook is None
        or get_column_letter is None
        or DefinedName is None
        or Font is None
        or DataValidation is None
    ):
        raise RuntimeError(
            "openpyxl is required to write .xlsx manifests. Install openpyxl or use --output *.tsv."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "download_plan"
    sheet.append(list(MANIFEST_FIELDNAMES))
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    if Comment is not None:
        sheet["B1"].comment = Comment(
            "IDs not listed in the drop-down are still supported if they are valid in the provider database.",
            "genegalleon",
        )
    sheet.freeze_panes = "A2"

    for row in rows:
        sheet.append([str(row.get(key, "") or "") for key in MANIFEST_FIELDNAMES])

    list_sheet = workbook.create_sheet("_lists")
    for idx, provider in enumerate(PROVIDERS, start=1):
        list_sheet.cell(row=idx, column=1, value=provider)
    id_options_by_provider = build_provider_id_options(rows, snapshot_options=id_options_snapshot)
    for col_idx, provider in enumerate(PROVIDERS, start=2):
        values = list(id_options_by_provider.get(provider, [""]))
        for row_idx, value in enumerate(values, start=1):
            list_sheet.cell(row=row_idx, column=col_idx, value=value)
        col_letter = get_column_letter(col_idx)
        workbook.defined_names.add(
            DefinedName(
                name="id_opts_{}".format(provider),
                attr_text="'_lists'!${}$1:${}${}".format(col_letter, col_letter, len(values)),
            )
        )
    list_sheet.sheet_state = "hidden"

    provider_validation = DataValidation(
        type="list",
        formula1="=_lists!$A$1:$A${}".format(len(PROVIDERS)),
        allow_blank=False,
    )
    sheet.add_data_validation(provider_validation)
    provider_validation.add("A2:A5000")

    id_validation = DataValidation(
        type="list",
        formula1='=INDIRECT("id_opts_"&$A2)',
        allow_blank=False,
        errorStyle="information",
        errorTitle="ID candidates",
        error=(
            "This drop-down changes by provider. "
            "For provider=ncbi, five model-organism IDs are shown as examples. "
            "You can still type any other ID value."
        ),
    )
    id_validation.showErrorMessage = False
    sheet.add_data_validation(id_validation)
    id_validation.add("B2:B5000")

    workbook.save(output_path)


def write_manifest(rows: List[Dict[str, str]], output_path: Path, id_options_snapshot=None):
    suffix = output_path.suffix.lower()
    if suffix == ".xlsx":
        write_manifest_xlsx(rows, output_path, id_options_snapshot=id_options_snapshot)
        return
    write_manifest_tsv(rows, output_path)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    input_root = Path(args.input_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    id_options_snapshot = {}
    if str(args.id_options_snapshot or "").strip() != "":
        snapshot_path = Path(args.id_options_snapshot).expanduser().resolve()
        if not snapshot_path.exists():
            sys.stderr.write("Warning: id options snapshot was not found and will be ignored: {}\n".format(snapshot_path))
        else:
            try:
                id_options_snapshot = load_id_options_snapshot(snapshot_path)
            except Exception as exc:
                sys.stderr.write("Error: failed to read id options snapshot '{}': {}\n".format(snapshot_path, exc))
                return 1

    requested_provider = str(args.provider or "").strip().lower()
    if requested_provider in LEGACY_NCBI_PROVIDER_ALIASES:
        requested_provider = "ncbi"
    providers = PROVIDERS if requested_provider == "all" else (requested_provider,)
    all_rows = []
    all_warnings = []
    all_errors = []
    allow_missing = requested_provider == "all"
    for provider in providers:
        rows, warnings, errors = discover(provider, input_root, allow_missing=allow_missing)
        all_rows.extend(rows)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    provider_order = {provider: idx for idx, provider in enumerate(PROVIDERS)}
    all_rows = sorted(
        all_rows,
        key=lambda x: (
            provider_order.get(str(x.get("provider", "")), len(PROVIDERS)),
            str(x.get("species_key", "")),
        ),
    )
    write_manifest(all_rows, output_path, id_options_snapshot=id_options_snapshot)
    print("Manifest written: {} (rows={})".format(output_path, len(all_rows)))

    if args.strict and len(all_errors) > 0:
        return 1
    if len(all_rows) == 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
