#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
import sys


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

PROVIDERS = ("ensemblplants", "phycocosm", "phytozome")

DEFAULT_INPUT_RELATIVE_DIRS = {
    "ensemblplants": Path("20230216_EnsemblPlants") / "original_files",
    "phycocosm": Path("PhycoCosm") / "species_wise_original",
    "phytozome": Path("Phytozome") / "species_wise_original",
}


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Build a download manifest TSV from an existing local gfe_dataset-style directory. "
            "Generated URLs are file:// URLs for reuse in format_species_inputs.py --download-manifest."
        )
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + PROVIDERS,
        default="all",
        help="Provider to scan (default: all).",
    )
    parser.add_argument(
        "--dataset-root",
        required=True,
        help="Root of gfe_dataset.",
    )
    parser.add_argument(
        "--output",
        default="workspace/input/query_manifest/download_manifest.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error when a species is missing CDS or GFF.",
    )
    return parser


def is_fasta_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in FASTA_EXTENSIONS)


def is_gff_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in GFF_EXTENSIONS)


def pick_single_file(matches, provider, species_key, label, warnings):
    if len(matches) == 0:
        return None
    ordered = sorted(matches)
    if len(ordered) > 1:
        warnings.append(
            "[{}] {}: multiple {} files found. Using '{}'".format(
                provider, species_key, label, ordered[0].name
            )
        )
    return ordered[0]


def discover_ensemblplants(input_dir):
    warnings = []
    errors = []
    rows = []

    species_to_cds = {}
    species_to_gff = {}
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

    for species_key in sorted(set(species_to_cds.keys()) | set(species_to_gff.keys())):
        cds_path = pick_single_file(species_to_cds.get(species_key, []), "ensemblplants", species_key, "CDS", warnings)
        gff_path = pick_single_file(species_to_gff.get(species_key, []), "ensemblplants", species_key, "GFF", warnings)
        if cds_path is None or gff_path is None:
            errors.append(
                "[ensemblplants] {}: missing {}".format(
                    species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        rows.append(
            {
                "provider": "ensemblplants",
                "species_key": species_key,
                "cds_url": cds_path.resolve().as_uri(),
                "gff_url": gff_path.resolve().as_uri(),
                "cds_filename": cds_path.name,
                "gff_filename": gff_path.name,
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
        else:
            cds_matches = [
                path
                for path in files
                if ("cds_" in path.name.lower() or ".cds." in path.name.lower()) and is_fasta_filename(path.name)
            ]
            gff_matches = [path for path in files if "gene.gff3" in path.name.lower() and is_gff_filename(path.name)]

        cds_path = pick_single_file(cds_matches, provider, species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, provider, species_key, "GFF", warnings)
        if cds_path is None or gff_path is None:
            errors.append(
                "[{}] {}: missing {}".format(
                    provider, species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        rows.append(
            {
                "provider": provider,
                "species_key": species_key,
                "cds_url": cds_path.resolve().as_uri(),
                "gff_url": gff_path.resolve().as_uri(),
                "cds_filename": cds_path.name,
                "gff_filename": gff_path.name,
            }
        )
    return rows, warnings, errors


def discover(provider, dataset_root):
    input_dir = dataset_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]
    if not input_dir.exists() or not input_dir.is_dir():
        return [], [], ["[{}] input directory not found: {}".format(provider, input_dir)]
    if provider == "ensemblplants":
        return discover_ensemblplants(input_dir)
    return discover_species_dir_based(provider, input_dir)


def write_manifest(rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["provider", "species_key", "cds_url", "gff_url", "cds_filename", "gff_filename"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    dataset_root = Path(args.dataset_root).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    providers = PROVIDERS if args.provider == "all" else (args.provider,)
    all_rows = []
    all_warnings = []
    all_errors = []
    for provider in providers:
        rows, warnings, errors = discover(provider, dataset_root)
        all_rows.extend(rows)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    all_rows = sorted(all_rows, key=lambda x: (x["provider"], x["species_key"]))
    write_manifest(all_rows, output_path)
    print("Manifest written: {} (rows={})".format(output_path, len(all_rows)))

    if args.strict and len(all_errors) > 0:
        return 1
    if len(all_rows) == 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
