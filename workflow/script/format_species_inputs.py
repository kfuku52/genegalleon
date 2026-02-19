#!/usr/bin/env python

import argparse
from collections import defaultdict
import csv
import gzip
import os
from pathlib import Path
import re
import sys
from urllib.parse import urlparse
from urllib.request import Request, urlopen


COMMON_REPLACEMENTS = (
    ("evm_27.model.", ""),
    ("evm.model.", ""),
    ("Oropetium_20150105_", ""),
)

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

ENSEMBL_CDS_PATTERN = re.compile(r"(?:^|[._-])cds(?:[._-]|$)", re.IGNORECASE)
INVALID_ID_CHARS = re.compile(r"[%/\+:;&\^\$#@!~=\'\"`*\(\)\{\}\[\]\|\?\s]+")

DEFAULT_INPUT_RELATIVE_DIRS = {
    "ensemblplants": Path("20230216_EnsemblPlants") / "original_files",
    "phycocosm": Path("PhycoCosm") / "species_wise_original",
    "phytozome": Path("Phytozome") / "species_wise_original",
}

PROVIDERS = ("ensemblplants", "phycocosm", "phytozome")


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Format CDS/GFF inputs for genegalleon from provider-specific raw files. "
            "Current targets are based on legacy data_formatting*.sh scripts."
        )
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + PROVIDERS,
        required=True,
        help="Input provider type. Use 'all' with --dataset-root to process all providers.",
    )
    parser.add_argument(
        "--dataset-root",
        default="",
        help=(
            "Root of gfe_dataset. If --input-dir is omitted, provider-specific default "
            "subdirectories are used from this root."
        ),
    )
    parser.add_argument(
        "--input-dir",
        default="",
        help=(
            "Provider input directory. For ensemblplants: original_files/. "
            "For phycocosm/phytozome: species_wise_original/. "
            "Not allowed when --provider all is used."
        ),
    )
    parser.add_argument(
        "--species-cds-dir",
        default="workspace/input/species_cds",
        help="Output directory for formatted species CDS FASTA files.",
    )
    parser.add_argument(
        "--species-gff-dir",
        default="workspace/input/species_gff",
        help="Output directory for formatted species GFF files.",
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
            "Optional TSV/CSV manifest for input download. "
            "Required columns: provider,species_key,cds_url,gff_url. "
            "Optional columns: cds_filename,gff_filename."
        ),
    )
    parser.add_argument(
        "--download-dir",
        default="workspace/tmp/input_download_cache",
        help=(
            "Directory for raw downloaded provider files. "
            "This can be used as --dataset-root for formatting."
        ),
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
        "--download-timeout",
        type=float,
        default=120.0,
        help="Timeout (seconds) per download request.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan actions without writing downloaded or formatted files.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error if any species is missing either CDS or GFF.",
    )
    return parser


def open_text(path, mode):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def detect_manifest_delimiter(path):
    with open(path, "rt", encoding="utf-8") as handle:
        first_line = handle.readline()
    if "\t" in first_line:
        return "\t"
    return ","


def provider_raw_dir(provider, download_root, species_key):
    if provider == "ensemblplants":
        return download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]
    if provider in ("phycocosm", "phytozome"):
        return download_root / DEFAULT_INPUT_RELATIVE_DIRS[provider] / species_key
    raise ValueError("Unknown provider: {}".format(provider))


def default_download_filename(provider, species_key, label, url):
    parsed = urlparse(url)
    base = Path(parsed.path).name
    if base == "":
        ext = ".fa.gz" if label == "cds" else ".gff3.gz"
        base = "{}.{}{}".format(species_key, label, ext)
    if provider == "ensemblplants":
        if not base.startswith(species_key + "."):
            base = "{}.{}".format(species_key, base)
    return base


def parse_http_headers(http_header_values, auth_bearer_token_env):
    headers = {}
    for raw in http_header_values:
        text = raw.strip()
        if text == "":
            continue
        if ":" not in text:
            raise ValueError("Invalid --http-header value (missing ':'): {}".format(raw))
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()
        if key == "":
            raise ValueError("Invalid --http-header value (empty key): {}".format(raw))
        headers[key] = value
    if auth_bearer_token_env != "":
        token = os.environ.get(auth_bearer_token_env, "").strip()
        if token == "":
            raise ValueError(
                "Environment variable for bearer token is empty or undefined: {}".format(
                    auth_bearer_token_env
                )
            )
        headers["Authorization"] = "Bearer {}".format(token)
    return headers


def download_url_to_file(url, destination, headers, timeout, dry_run):
    if dry_run:
        return
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    request = Request(url, headers=headers)
    with urlopen(request, timeout=timeout) as response, open(tmp, "wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(destination)


def read_download_manifest(path):
    delimiter = detect_manifest_delimiter(path)
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = [row for row in reader]
    return rows


def download_from_manifest(
    manifest_path,
    download_root,
    provider_filter,
    overwrite,
    headers,
    timeout,
    dry_run,
):
    rows = read_download_manifest(manifest_path)
    warnings = []
    errors = []
    processed = 0
    downloaded = 0
    planned = 0

    required_cols = {"provider", "species_key", "cds_url", "gff_url"}
    if len(rows) == 0:
        errors.append("Download manifest is empty: {}".format(manifest_path))
        return {
            "warnings": warnings,
            "errors": errors,
            "processed": processed,
            "downloaded": downloaded,
            "planned": planned,
        }
    missing_cols = required_cols - set(rows[0].keys())
    if missing_cols:
        errors.append(
            "Download manifest missing required columns {}: {}".format(sorted(missing_cols), manifest_path)
        )
        return {
            "warnings": warnings,
            "errors": errors,
            "processed": processed,
            "downloaded": downloaded,
            "planned": planned,
        }

    for i, row in enumerate(rows, start=2):
        provider = (row.get("provider") or "").strip().lower()
        species_key = (row.get("species_key") or "").strip()
        cds_url = (row.get("cds_url") or "").strip()
        gff_url = (row.get("gff_url") or "").strip()
        cds_filename = (row.get("cds_filename") or "").strip()
        gff_filename = (row.get("gff_filename") or "").strip()

        if provider_filter != "all" and provider != provider_filter:
            continue
        processed += 1

        if provider not in PROVIDERS:
            errors.append("Manifest line {}: unsupported provider '{}'".format(i, provider))
            continue
        if species_key == "":
            errors.append("Manifest line {}: species_key is empty".format(i))
            continue
        if cds_url == "" or gff_url == "":
            errors.append("Manifest line {}: cds_url/gff_url is empty for {}".format(i, species_key))
            continue

        raw_dir = provider_raw_dir(provider, download_root, species_key)
        raw_dir.mkdir(parents=True, exist_ok=True)

        if cds_filename == "":
            cds_filename = default_download_filename(provider, species_key, "cds", cds_url)
        if gff_filename == "":
            gff_filename = default_download_filename(provider, species_key, "gff", gff_url)

        targets = [
            ("CDS", cds_url, raw_dir / cds_filename),
            ("GFF", gff_url, raw_dir / gff_filename),
        ]
        for label, url, target in targets:
            if target.exists() and target.stat().st_size > 0 and not overwrite:
                warnings.append(
                    "[download:{}] {} {} already exists. Skipping: {}".format(
                        provider, species_key, label, target
                    )
                )
                continue
            if dry_run:
                planned += 1
                continue
            try:
                download_url_to_file(url, target, headers=headers, timeout=timeout, dry_run=dry_run)
                downloaded += 1
            except Exception as exc:
                errors.append(
                    "[download:{}] failed {} {} from {} -> {} ({})".format(
                        provider, species_key, label, url, target, exc
                    )
                )

    if processed == 0:
        warnings.append("No manifest rows matched --provider {}.".format(provider_filter))

    return {
        "warnings": warnings,
        "errors": errors,
        "processed": processed,
        "downloaded": downloaded,
        "planned": planned,
    }


def is_fasta_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in FASTA_EXTENSIONS)


def is_gff_filename(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in GFF_EXTENSIONS)


def first_token(text):
    parts = text.split()
    if len(parts) == 0:
        return ""
    return parts[0]


def species_prefix_from_value(value):
    tokens = value.split("_")
    if len(tokens) < 2:
        return ""
    return "{}_{}".format(tokens[0], tokens[1])


def normalize_output_basename(source_name, species_prefix):
    if source_name.startswith(species_prefix + "_"):
        return source_name
    dotted_prefix = species_prefix + "."
    if source_name.startswith(dotted_prefix):
        return species_prefix + "_" + source_name[len(dotted_prefix):]
    return species_prefix + "_" + source_name


def apply_common_replacements(text):
    out = text
    for old, new in COMMON_REPLACEMENTS:
        out = out.replace(old, new)
    return out


def sanitize_identifier(identifier):
    out = identifier.replace("−", "-")
    out = apply_common_replacements(out)
    out = INVALID_ID_CHARS.sub("_", out)
    out = re.sub(r"_+", "_", out)
    out = out.strip("_")
    if out == "":
        out = "unnamed"
    return out


def pad_to_codon_length(seq):
    remainder = len(seq) % 3
    if remainder == 0:
        return seq
    return seq + ("N" * (3 - remainder))


def iter_fasta_records(path):
    header = None
    chunks = []
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
                continue
            chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        yield header, "".join(chunks)


def write_fasta_record(handle, record_id, sequence, width=80):
    handle.write(">{}\n".format(record_id))
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i + width] + "\n")


def extract_ensembl_id(header):
    match = re.search(r"(?:^|\s)gene:([^\s]+)", header)
    if match:
        return match.group(1)
    return first_token(header)


def extract_phycocosm_id(header):
    token = first_token(header)
    if token.startswith("jgi|"):
        parts = [part for part in token.split("|") if part != ""]
        if len(parts) >= 2:
            return parts[-1]
    if "|" in token:
        return token.split("|")[-1]
    return token


def extract_phytozome_id(header):
    return first_token(header)


def extract_provider_id(provider, header):
    if provider == "ensemblplants":
        return extract_ensembl_id(header)
    if provider == "phycocosm":
        return extract_phycocosm_id(header)
    return extract_phytozome_id(header)


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


def discover_ensemblplants_tasks(input_dir):
    warnings = []
    errors = []
    tasks = []

    cds_by_species = defaultdict(list)
    gff_by_species = defaultdict(list)

    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        name = path.name
        species_key = name.split(".", 1)[0]
        if species_key == "":
            continue
        if is_fasta_filename(name) and ENSEMBL_CDS_PATTERN.search(name):
            cds_by_species[species_key].append(path)
            continue
        if is_gff_filename(name):
            gff_by_species[species_key].append(path)

    for species_key in sorted(set(cds_by_species.keys()) | set(gff_by_species.keys())):
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append(
                "[ensemblplants] skipped '{}': unable to parse species prefix.".format(species_key)
            )
            continue

        cds_path = pick_single_file(
            cds_by_species.get(species_key, []),
            "ensemblplants",
            species_key,
            "CDS",
            warnings,
        )
        gff_path = pick_single_file(
            gff_by_species.get(species_key, []),
            "ensemblplants",
            species_key,
            "GFF",
            warnings,
        )
        if cds_path is None or gff_path is None:
            errors.append(
                "[ensemblplants] {}: missing {}".format(
                    species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        tasks.append(
            {
                "provider": "ensemblplants",
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
            }
        )

    return tasks, warnings, errors


def discover_phycocosm_tasks(input_dir):
    warnings = []
    errors = []
    tasks = []

    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = species_dir.name
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append(
                "[phycocosm] skipped '{}': unable to parse species prefix.".format(species_key)
            )
            continue

        files = [path for path in species_dir.iterdir() if path.is_file()]
        cds_matches = [path for path in files if "fasta" in path.name.lower() and is_fasta_filename(path.name)]
        gff_matches = [path for path in files if "gff" in path.name.lower() and is_gff_filename(path.name)]

        cds_path = pick_single_file(cds_matches, "phycocosm", species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, "phycocosm", species_key, "GFF", warnings)

        if cds_path is None or gff_path is None:
            errors.append(
                "[phycocosm] {}: missing {}".format(
                    species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        tasks.append(
            {
                "provider": "phycocosm",
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
            }
        )

    return tasks, warnings, errors


def discover_phytozome_tasks(input_dir):
    warnings = []
    errors = []
    tasks = []

    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = species_dir.name
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append(
                "[phytozome] skipped '{}': unable to parse species prefix.".format(species_key)
            )
            continue

        files = [path for path in species_dir.iterdir() if path.is_file()]
        cds_matches = [
            path
            for path in files
            if ("cds_" in path.name.lower() or ".cds." in path.name.lower()) and is_fasta_filename(path.name)
        ]
        gff_matches = [path for path in files if "gene.gff3" in path.name.lower() and is_gff_filename(path.name)]

        cds_path = pick_single_file(cds_matches, "phytozome", species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, "phytozome", species_key, "GFF", warnings)

        if cds_path is None or gff_path is None:
            errors.append(
                "[phytozome] {}: missing {}".format(
                    species_key, "CDS" if cds_path is None else "GFF"
                )
            )
            continue
        tasks.append(
            {
                "provider": "phytozome",
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
            }
        )

    return tasks, warnings, errors


def discover_tasks(provider, input_dir):
    if provider == "ensemblplants":
        return discover_ensemblplants_tasks(input_dir)
    if provider == "phycocosm":
        return discover_phycocosm_tasks(input_dir)
    if provider == "phytozome":
        return discover_phytozome_tasks(input_dir)
    raise ValueError("Unknown provider: {}".format(provider))


def format_cds(task, output_dir, overwrite, dry_run):
    output_name = normalize_output_basename(task["cds_path"].name, task["species_prefix"])
    output_path = output_dir / output_name
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        return {"status": "skip", "output_path": output_path, "written": 0, "duplicates": 0}
    if dry_run:
        return {"status": "dry-run", "output_path": output_path, "written": 0, "duplicates": 0}

    seen_ids = set()
    written = 0
    duplicates = 0

    with open_text(output_path, "wt") as fout:
        for header, sequence in iter_fasta_records(task["cds_path"]):
            extracted = extract_provider_id(task["provider"], header)
            sanitized = sanitize_identifier(extracted)
            prefixed = "{}_{}".format(task["species_prefix"], sanitized)
            final_id = sanitize_identifier(prefixed)
            if final_id in seen_ids:
                duplicates += 1
                continue
            seen_ids.add(final_id)
            seq = re.sub(r"\s+", "", sequence).upper()
            seq = pad_to_codon_length(seq)
            write_fasta_record(fout, final_id, seq)
            written += 1
    return {"status": "write", "output_path": output_path, "written": written, "duplicates": duplicates}


def format_gff(task, output_dir, overwrite, dry_run):
    output_name = normalize_output_basename(task["gff_path"].name, task["species_prefix"])
    output_path = output_dir / output_name
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        return {"status": "skip", "output_path": output_path, "lines": 0}
    if dry_run:
        return {"status": "dry-run", "output_path": output_path, "lines": 0}

    line_count = 0
    with open_text(task["gff_path"], "rt") as fin, open_text(output_path, "wt") as fout:
        for line in fin:
            fout.write(apply_common_replacements(line))
            line_count += 1
    return {"status": "write", "output_path": output_path, "lines": line_count}


def resolve_provider_inputs(args):
    if args.provider == "all":
        if args.input_dir != "":
            raise ValueError("--input-dir cannot be used with --provider all.")
        if args.dataset_root == "":
            raise ValueError("--dataset-root is required when --provider all is used.")
        dataset_root = Path(args.dataset_root).expanduser().resolve()
        return [
            (provider, (dataset_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]).resolve())
            for provider in PROVIDERS
        ]

    provider = args.provider
    if args.input_dir != "":
        return [(provider, Path(args.input_dir).expanduser().resolve())]
    if args.dataset_root != "":
        dataset_root = Path(args.dataset_root).expanduser().resolve()
        return [(provider, (dataset_root / DEFAULT_INPUT_RELATIVE_DIRS[provider]).resolve())]
    raise ValueError("Specify either --input-dir or --dataset-root.")


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    if args.download_only and args.download_manifest == "":
        parser.error("--download-only requires --download-manifest.")

    try:
        http_headers = parse_http_headers(args.http_header, args.auth_bearer_token_env)
    except ValueError as exc:
        parser.error(str(exc))

    if args.download_manifest != "":
        download_root = Path(args.download_dir).expanduser().resolve()
        download_root.mkdir(parents=True, exist_ok=True)
        report = download_from_manifest(
            manifest_path=Path(args.download_manifest).expanduser().resolve(),
            download_root=download_root,
            provider_filter=args.provider,
            overwrite=args.overwrite,
            headers=http_headers,
            timeout=float(args.download_timeout),
            dry_run=args.dry_run,
        )
        for warning in report["warnings"]:
            sys.stderr.write("Warning: {}\n".format(warning))
        for error in report["errors"]:
            sys.stderr.write("Error: {}\n".format(error))
        print(
            "Download stage: rows processed={}, files downloaded={}, planned downloads={}, download root={}, dry_run={}".format(
                report["processed"],
                report["downloaded"],
                report["planned"],
                download_root,
                int(args.dry_run),
            )
        )
        if args.download_only:
            return 0 if len(report["errors"]) == 0 else 1
        if len(report["errors"]) > 0 and args.strict:
            return 1

    if args.input_dir == "" and args.dataset_root == "" and args.download_manifest != "":
        args.dataset_root = str(Path(args.download_dir).expanduser().resolve())

    try:
        provider_inputs = resolve_provider_inputs(args)
    except ValueError as exc:
        parser.error(str(exc))

    output_cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    output_gff_dir = Path(args.species_gff_dir).expanduser().resolve()
    output_cds_dir.mkdir(parents=True, exist_ok=True)
    output_gff_dir.mkdir(parents=True, exist_ok=True)

    all_tasks = []
    all_warnings = []
    all_errors = []

    for provider, input_dir in provider_inputs:
        if not input_dir.exists() or not input_dir.is_dir():
            all_errors.append("[{}] input directory not found: {}".format(provider, input_dir))
            continue
        tasks, warnings, errors = discover_tasks(provider, input_dir)
        all_tasks.extend(tasks)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    if args.strict and len(all_errors) > 0:
        return 1
    if len(all_tasks) == 0:
        sys.stderr.write("No species pairs were discovered.\n")
        return 1

    processed = 0
    total_duplicates = 0
    for task in all_tasks:
        cds_result = format_cds(task, output_cds_dir, args.overwrite, args.dry_run)
        gff_result = format_gff(task, output_gff_dir, args.overwrite, args.dry_run)
        processed += 1
        total_duplicates += cds_result["duplicates"]
        print(
            "[{}] {}: CDS={} ({}, {}, duplicates={}), GFF={} ({}, lines={})".format(
                task["provider"],
                task["species_prefix"],
                task["cds_path"].name,
                cds_result["status"],
                cds_result["output_path"].name,
                cds_result["duplicates"],
                task["gff_path"].name,
                gff_result["status"],
                gff_result["lines"],
            )
        )

    print(
        "Finished. species processed: {}, duplicate CDS IDs skipped: {}, output CDS dir: {}, output GFF dir: {}, dry_run={}".format(
            processed,
            total_duplicates,
            output_cds_dir,
            output_gff_dir,
            int(args.dry_run),
        )
    )
    if len(all_errors) > 0:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
