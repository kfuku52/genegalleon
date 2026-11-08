"""Resolve per-species genetic codes using the shared species-label contract."""

import argparse
import csv
import sys
from pathlib import Path

from species_labeling import extract_species_label

FASTA_SUFFIXES = (".fa", ".fa.gz", ".fas", ".fas.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz")


def write_table(cds_dir, default_code, outfile, input_table=None):
    species = sorted(
        extract_species_label(path.name, strip_extension=True) or path.stem
        for path in Path(cds_dir).iterdir()
        if path.is_file() and not path.name.startswith(".") and path.name.endswith(FASTA_SUFFIXES)
    )
    if not species:
        raise ValueError(f"No species CDS files were found in: {cds_dir}")
    overrides = {}
    if input_table is not None and Path(input_table).exists():
        with Path(input_table).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(
                (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
                delimiter="\t",
            )
            if reader.fieldnames is None or not {"species", "genetic_code"}.issubset(reader.fieldnames):
                raise ValueError(f"{input_table} must be a tab-delimited file with at least species and genetic_code columns.")
            for row in reader:
                name = (row.get("species") or "").strip().replace(" ", "_")
                code = (row.get("genetic_code") or "").strip()
                if not name or not code:
                    raise ValueError(f"Empty species or genetic_code entry was detected in {input_table}.")
                try:
                    code_int = int(code)
                except ValueError:
                    raise ValueError(f"Invalid genetic_code for {name}: {code}") from None
                if code_int <= 0:
                    raise ValueError(f"genetic_code must be a positive integer for {name}: {code}")
                if name in overrides:
                    raise ValueError(f"Duplicate species entry in {input_table}: {name}")
                overrides[name] = str(code_int)
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with outfile.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "genetic_code", "source"])
        for name in species:
            writer.writerow([name, overrides.get(name, str(default_code).strip()),
                             "species_genetic_code.tsv" if name in overrides else "default"])
    unknown = sorted(set(overrides).difference(species))
    if unknown:
        sys.stderr.write(
            "Warning: species_genetic_code.tsv entries were ignored because no matching species_cds file was found: "
            + ", ".join(unknown) + "\n"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cds_dir", type=Path)
    parser.add_argument("default_code")
    parser.add_argument("outfile", type=Path)
    parser.add_argument("input_table", nargs="?", default="")
    args = parser.parse_args()
    try:
        write_table(args.cds_dir, args.default_code, args.outfile, Path(args.input_table) if args.input_table else None)
    except ValueError as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
