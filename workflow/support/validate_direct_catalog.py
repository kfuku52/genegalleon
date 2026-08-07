#!/usr/bin/env python3
"""Validate curated direct-download bundles before publishing a download plan."""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from urllib.request import Request, urlopen

SUSPICIOUS_GFF_PATTERN = re.compile(
    r"(?:^|[._/-])(?:cds0|edta|reftrans(?:alignments?)?|repeat|teanno|transpos(?:on)?)(?:[._/-]|$)",
    re.IGNORECASE,
)
TRANSCRIPTOME_CDS_PATTERN = re.compile(
    r"(?:onekp|soapdenovo[-_.]?trans|transcriptome|trans[-_.]?assembly)",
    re.IGNORECASE,
)


class DownloadLimitReached(Exception):
    pass


class LimitedReader:
    def __init__(self, handle, limit):
        self.handle = handle
        self.remaining = int(limit)

    def read(self, size=-1):
        if self.remaining <= 0:
            raise DownloadLimitReached()
        requested = self.remaining if size is None or size < 0 else min(int(size), self.remaining)
        data = self.handle.read(requested)
        self.remaining -= len(data)
        return data

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)

    def readable(self):
        return True


def build_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--check-remote-gff", action="store_true")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--max-download-bytes", type=int, default=128 * 1024 * 1024)
    return parser


def nonempty(row, key):
    return str(row.get(key) or "").strip()


def validate_row_structure(row, line_number):
    errors = []
    species_key = nonempty(row, "species_key")
    label = species_key or "line {}".format(line_number)
    cds_url = nonempty(row, "cds_url") or nonempty(row, "local_cds_path")
    gff_url = nonempty(row, "gff_url") or nonempty(row, "local_gff_path")
    gbff_url = nonempty(row, "gbff_url") or nonempty(row, "local_gbff_path")
    genome_url = nonempty(row, "genome_url") or nonempty(row, "local_genome_path")
    if species_key == "":
        errors.append("line {}: missing species_key".format(line_number))
    if cds_url == "" and gbff_url == "" and not (gff_url != "" and genome_url != ""):
        errors.append("{}: no usable CDS, GBFF, or GFF+genome bundle".format(label))

    cds_descriptor = " ".join((nonempty(row, "cds_filename"), nonempty(row, "cds_url")))
    if cds_url != "" and gff_url == "" and gbff_url == "" and TRANSCRIPTOME_CDS_PATTERN.search(cds_descriptor):
        errors.append("{}: transcriptome assembly is not an accepted CDS-only bundle".format(label))

    gff_descriptor = " ".join((nonempty(row, "gff_filename"), nonempty(row, "gff_url")))
    if gff_url != "" and SUSPICIOUS_GFF_PATTERN.search(gff_descriptor):
        errors.append("{}: suspicious non-gene-model GFF source".format(label))
    return errors


def remote_gff_contains_cds(row, timeout, max_download_bytes):
    url = nonempty(row, "gff_url")
    if url == "":
        return "not_applicable", ""
    if nonempty(row, "gff_archive_member") != "":
        return "indeterminate", "archive member requires full archive validation"

    filename = nonempty(row, "gff_filename") or url.split("?", 1)[0]
    request = Request(url, headers={"User-Agent": "GeneGalleon-catalog-validator/1"})
    try:
        with urlopen(request, timeout=timeout) as response:
            limited = LimitedReader(response, max_download_bytes)
            stream = gzip.GzipFile(fileobj=limited, mode="rb") if filename.lower().endswith(".gz") else limited
            pending = b""
            while True:
                chunk = stream.read(1024 * 1024)
                if not chunk:
                    break
                pending += chunk
                lines = pending.split(b"\n")
                pending = lines.pop()
                for raw_line in lines:
                    if raw_line.startswith(b"#"):
                        continue
                    parts = raw_line.rstrip(b"\r").split(b"\t", 3)
                    if len(parts) >= 3 and parts[2].strip().lower() == b"cds":
                        return "cds", ""
            if pending:
                parts = pending.rstrip(b"\r").split(b"\t", 3)
                if len(parts) >= 3 and parts[2].strip().lower() == b"cds":
                    return "cds", ""
            return "no_cds", "complete GFF contained no CDS feature"
    except DownloadLimitReached:
        return "indeterminate", "download limit reached before a CDS feature"
    except Exception as exc:
        return "indeterminate", "{}: {}".format(type(exc).__name__, exc)


def main():
    args = build_arg_parser().parse_args()
    manifest = Path(args.manifest).expanduser().resolve()
    if not manifest.is_file():
        sys.stderr.write("Error: manifest not found: {}\n".format(manifest))
        return 2

    errors = []
    warnings = []
    rows = []
    seen_keys = set()
    with open(manifest, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, 2):
            rows.append(row)
            key = (nonempty(row, "provider"), nonempty(row, "species_key"))
            if key in seen_keys:
                errors.append("line {}: duplicate provider/species_key {}:{}".format(line_number, *key))
            seen_keys.add(key)
            errors.extend(validate_row_structure(row, line_number))

    remote_checked = 0
    if args.check_remote_gff:
        for row in rows:
            status, detail = remote_gff_contains_cds(row, args.timeout, args.max_download_bytes)
            if status == "not_applicable":
                continue
            remote_checked += 1
            species_key = nonempty(row, "species_key")
            if status == "no_cds":
                errors.append("{}: {}".format(species_key, detail))
            elif status == "indeterminate":
                warnings.append("{}: remote GFF check indeterminate ({})".format(species_key, detail))

    for warning in warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in errors:
        sys.stderr.write("Error: {}\n".format(error))
    print(
        "Direct catalog validation: rows={}, remote_gff_checked={}, warnings={}, errors={}".format(
            len(rows), remote_checked, len(warnings), len(errors)
        )
    )
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
