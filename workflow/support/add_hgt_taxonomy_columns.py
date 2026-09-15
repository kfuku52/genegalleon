#!/usr/bin/env python3

"""Add taxonomy columns to already-written HGT candidate TSVs.

This is an additive migration for outputs made by an older GeneGalleon
version.  It preserves the existing columns and row order, and replaces the
three input files only after all enriched temporary files have been written.
"""

import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from score_hgt_candidates import (
    BRANCH_TAXONOMY_COLUMNS,
    GENE_TAXONOMY_COLUMNS,
    ORTHOGROUP_TAXONOMY_COLUMNS,
    REPRESENTATIVE_BRANCH_COLUMNS,
    TAXONOMIC_RANKS,
    TaxonomyResolver,
    join_unique_nonempty,
    representative_branch_annotation,
    resolve_taxonomy_annotation,
    taxonomy_lineage_column,
    taxonomy_rank_column,
    taxonomy_rank_list_column,
)

GENE_EXTRA_COLUMNS = list(GENE_TAXONOMY_COLUMNS)
BRANCH_EXTRA_COLUMNS = [*BRANCH_TAXONOMY_COLUMNS, *REPRESENTATIVE_BRANCH_COLUMNS]
ORTHOGROUP_EXTRA_COLUMNS = list(ORTHOGROUP_TAXONOMY_COLUMNS)


def split_gene_ids(value: object) -> List[str]:
    if value is None:
        return []
    return [token.strip() for token in str(value).split(";") if token.strip()]


def nonempty(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "<na>"}:
        return ""
    return text


def base_fieldnames(fieldnames: Sequence[str], extra_columns: Sequence[str]) -> List[str]:
    extra = set(extra_columns)
    return [field for field in fieldnames if field not in extra]


def read_fieldnames(path: Path) -> List[str]:
    if not path.is_file():
        raise FileNotFoundError(f"HGT table was not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        fieldnames = next(csv.reader(handle, delimiter="\t"), [])
    fieldnames = [str(field).strip() for field in fieldnames]
    if not fieldnames or any(field == "" for field in fieldnames):
        raise ValueError(f"HGT table has an invalid header: {path}")
    return fieldnames


def make_temp_path(path: Path, created: List[Path]) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(descriptor)
    temp_path = Path(name)
    created.append(temp_path)
    return temp_path


def write_row(writer: csv.DictWriter, row: Dict[str, object], fieldnames: Sequence[str], values: Dict[str, str]) -> None:
    output = {field: row.get(field, "") if row.get(field, "") is not None else "" for field in fieldnames}
    output.update(values)
    writer.writerow(output)


def aggregate_taxonomy_fields(records: Sequence[Dict[str, object]]) -> Dict[str, str]:
    """Collapse per-gene rank annotations into unique semicolon-separated lists."""
    values: Dict[str, str] = {}
    for rank in TAXONOMIC_RANKS:
        values[taxonomy_rank_list_column("recipient", rank)] = join_unique_nonempty(
            record.get(taxonomy_rank_column("recipient", rank), "") for record in records
        )
        values[taxonomy_rank_list_column("donor", rank)] = join_unique_nonempty(
            record.get(taxonomy_rank_column("donor", rank), "") for record in records
        )
    values[taxonomy_lineage_column("recipient", plural=True)] = join_unique_nonempty(
        record.get(taxonomy_lineage_column("recipient"), "") for record in records
    )
    values[taxonomy_lineage_column("donor", plural=True)] = join_unique_nonempty(
        record.get(taxonomy_lineage_column("donor"), "") for record in records
    )
    return values


def preserve_existing_values(row: Dict[str, object], values: Dict[str, str]) -> Dict[str, str]:
    """Keep already populated migration fields while filling newly added fields."""
    output = dict(values)
    for field in values:
        existing = nonempty(row.get(field))
        if existing:
            output[field] = existing
    return output


def enrich_hgt_tables(
    branch_path: Path,
    gene_path: Path,
    orthogroup_path: Path,
    taxonomy_dbfile: str,
) -> Dict[str, int]:
    """Append complete taxonomy fields to the three HGT tables and return row counts."""
    paths = {
        "branch": branch_path,
        "gene": gene_path,
        "orthogroup": orthogroup_path,
    }
    headers = {key: read_fieldnames(path) for key, path in paths.items()}
    required = {
        "branch": {"orthogroup", "candidate_genes"},
        "gene": {"orthogroup", "gene_id", "gene_taxon", "besthit_organism", "besthit_taxid"},
        "orthogroup": {"orthogroup"},
    }
    for key, required_columns in required.items():
        missing = sorted(required_columns - set(headers[key]))
        if missing:
            raise ValueError(f"{paths[key]} is missing required columns: {', '.join(missing)}")

    resolver = TaxonomyResolver(taxonomy_dbfile)
    if not resolver.enabled:
        print(
            "Warning: taxonomy DB could not be opened; existing taxonomy values will be preserved and missing values remain blank.",
            file=sys.stderr,
        )

    created: List[Path] = []
    gene_records_by_key: Dict[Tuple[str, str], Dict[str, object]] = {}
    row_counts: Dict[str, int] = {}
    try:
        gene_temp = make_temp_path(gene_path, created)
        gene_fields = base_fieldnames(headers["gene"], GENE_EXTRA_COLUMNS)
        with gene_path.open("r", encoding="utf-8", newline="") as source, gene_temp.open(
            "w", encoding="utf-8", newline=""
        ) as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(
                target,
                fieldnames=gene_fields + GENE_EXTRA_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            count = 0
            for row in reader:
                orthogroup = nonempty(row.get("orthogroup"))
                gene_id = nonempty(row.get("gene_id"))
                recipient_annotation = resolve_taxonomy_annotation(row.get("gene_taxon", ""), "", resolver)
                donor_annotation = resolve_taxonomy_annotation(
                    row.get("besthit_organism", ""), row.get("besthit_taxid", ""), resolver
                )
                taxonomy_values: Dict[str, str] = {}
                for rank in TAXONOMIC_RANKS:
                    recipient_field = taxonomy_rank_column("recipient", rank)
                    donor_field = taxonomy_rank_column("donor", rank)
                    taxonomy_values[recipient_field] = nonempty(row.get(recipient_field)) or recipient_annotation.get(
                        rank, ""
                    )
                    taxonomy_values[donor_field] = nonempty(row.get(donor_field)) or donor_annotation.get(rank, "")
                recipient_lineage_field = taxonomy_lineage_column("recipient")
                donor_lineage_field = taxonomy_lineage_column("donor")
                taxonomy_values[recipient_lineage_field] = nonempty(row.get(recipient_lineage_field)) or recipient_annotation.get(
                    "taxonomy", ""
                )
                taxonomy_values[donor_lineage_field] = nonempty(row.get(donor_lineage_field)) or donor_annotation.get(
                    "taxonomy", ""
                )
                gene_record = dict(row)
                gene_record.update(taxonomy_values)
                gene_records_by_key.setdefault((orthogroup, gene_id), gene_record)
                write_row(
                    writer,
                    row,
                    gene_fields,
                    taxonomy_values,
                )
                count += 1
            row_counts["gene"] = count

        branch_temp = make_temp_path(branch_path, created)
        branch_fields = base_fieldnames(headers["branch"], BRANCH_EXTRA_COLUMNS)
        with branch_path.open("r", encoding="utf-8", newline="") as source, branch_temp.open(
            "w", encoding="utf-8", newline=""
        ) as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(
                target,
                fieldnames=branch_fields + BRANCH_EXTRA_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            count = 0
            for row in reader:
                orthogroup = nonempty(row.get("orthogroup"))
                branch_gene_records = [
                    gene_records_by_key.get((orthogroup, gene_id), {})
                    for gene_id in split_gene_ids(row.get("candidate_genes"))
                ]
                representative = representative_branch_annotation(branch_gene_records)
                branch_values = aggregate_taxonomy_fields(branch_gene_records)
                branch_values.update(representative)
                branch_values = preserve_existing_values(row, branch_values)
                write_row(
                    writer,
                    row,
                    branch_fields,
                    branch_values,
                )
                count += 1
            row_counts["branch"] = count

        orthogroup_temp = make_temp_path(orthogroup_path, created)
        orthogroup_fields = base_fieldnames(headers["orthogroup"], ORTHOGROUP_EXTRA_COLUMNS)
        with orthogroup_path.open("r", encoding="utf-8", newline="") as source, orthogroup_temp.open(
            "w", encoding="utf-8", newline=""
        ) as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(
                target,
                fieldnames=orthogroup_fields + ORTHOGROUP_EXTRA_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            count = 0
            for row in reader:
                orthogroup = nonempty(row.get("orthogroup"))
                records = [
                    record
                    for (record_orthogroup, _gene_id), record in gene_records_by_key.items()
                    if record_orthogroup == orthogroup
                ]
                taxonomy_values = aggregate_taxonomy_fields(records)
                write_row(
                    writer,
                    row,
                    orthogroup_fields,
                    preserve_existing_values(row, taxonomy_values),
                )
                count += 1
            row_counts["orthogroup"] = count

        replacements = [
            (gene_temp, gene_path),
            (branch_temp, branch_path),
            (orthogroup_temp, orthogroup_path),
        ]
        for temporary, destination in replacements:
            os.replace(temporary, destination)
            created.remove(temporary)
    finally:
        for temporary in created:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
    return row_counts


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add taxonomy columns to existing GeneGalleon HGT tables.")
    parser.add_argument("--branch_tsv", required=True, type=Path)
    parser.add_argument("--gene_tsv", required=True, type=Path)
    parser.add_argument("--orthogroup_tsv", required=True, type=Path)
    parser.add_argument("--taxonomy_dbfile", default=None, type=str)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    taxonomy_dbfile = os.environ.get("GG_TAXONOMY_DBFILE", "") if args.taxonomy_dbfile is None else args.taxonomy_dbfile
    counts = enrich_hgt_tables(
        branch_path=args.branch_tsv,
        gene_path=args.gene_tsv,
        orthogroup_path=args.orthogroup_tsv,
        taxonomy_dbfile=taxonomy_dbfile,
    )
    print(
        "Added HGT taxonomy columns: "
        + ", ".join(f"{name}={counts.get(name, 0)} rows" for name in ("branch", "gene", "orthogroup"))
    )


if __name__ == "__main__":
    main()
