"""Select BUSCO representatives without changing the shared all-copy summary."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from busco_table_utils import BUSCO_ID_PATTERN, SOURCE_COORDINATE_PATTERN, normalize_species_label

SELECTION_COLUMNS = (
    "busco_id", "species", "gene_id", "busco_sequence_ids", "status", "score", "busco_length",
    "sequence_mode", "genetic_code", "sequence_length", "internal_stops", "ambiguous_fraction", "selected", "reason",
    "selected_gene_id", "runner_up_gene_id", "score_margin", "selection_rule",
)
SELECTION_RULE = "prefer_no_internal_stop;score_desc;ambiguous_fraction_asc;gene_id_asc"
AMBIGUITY_DELETE_TABLE = str.maketrans("", "", "ACDEFGHIKLMNPQRSTVWYOU*")


def source_gene_id(identifier):
    # Match collect_common_BUSCO_genes.py's transcriptome coordinate removal.
    return SOURCE_COORDINATE_PATTERN.sub("", identifier)


def finite_score(value, context):
    try:
        score = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid BUSCO score for {context}: {value!r}") from exc
    if not math.isfinite(score):
        raise ValueError(f"Invalid BUSCO score for {context}: {value!r}")
    return score


def sequence_quality(record, mode, genetic_code):
    sequence = "".join("".join(record[1:]).split()).upper()
    if not sequence:
        raise ValueError(f"Empty candidate sequence: {record[0].strip()}")
    if mode == "cds":
        from Bio.Data import CodonTable
        from Bio.Seq import Seq

        # Padding affects quality assessment only; extraction preserves the CDS.
        sequence = sequence.replace("X", "N")
        try:
            peptide = str(Seq(sequence + "N" * (-len(sequence) % 3)).translate(table=genetic_code))
        except CodonTable.TranslationError as exc:
            raise ValueError(f"Invalid candidate CDS: {record[0].strip()}: {exc}") from exc
    else:
        peptide = sequence
    stops = peptide[:-1].count("*")
    peptide = peptide.removesuffix("*")
    ambiguous = len(peptide.translate(AMBIGUITY_DELETE_TABLE))
    return len(sequence), stops, ambiguous / max(1, len(peptide))


class BestScoreSelection:
    def __init__(self, args, fasta_paths):
        self.args = args
        self.groups = defaultdict(dict)
        self.singles = defaultdict(dict)
        self.expected = defaultdict(int)
        self.report_rows = []
        with args.summary.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            species = header[3:]
            if not species or any(not label.strip() for label in species) or len(species) != len(set(species)):
                raise ValueError("BUSCO summary requires unique species columns")
            if args.species_headers and any(">" in label or any(char.isspace() for char in label) for label in species):
                raise ValueError("Species FASTA header labels must not contain whitespace or '>'")
            summary = {}
            for fields in reader:
                if not fields or not fields[0]:
                    continue
                if len(fields) != len(header) or fields[0] in summary:
                    raise ValueError(f"Malformed or repeated BUSCO summary row: {fields[0]}")
                if fields[0] in {".", ".."} or BUSCO_ID_PATTERN.fullmatch(fields[0]) is None:
                    raise ValueError(f"Invalid BUSCO identifier in summary: {fields[0]!r}")
                summary[fields[0]] = fields[3:]

        with args.species_fasta_map.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["species", "fasta"]:
                raise ValueError("Species FASTA map requires species and fasta columns")
            rows = list(reader)
        self.fasta_species = {}
        seen_species = set()
        for line_number, row in enumerate(rows, start=2):
            if None in row or not row["species"] or not row["fasta"]:
                raise ValueError(f"Malformed species FASTA map row {line_number}: {args.species_fasta_map}")
            label, path = row["species"], Path(row["fasta"])
            if label not in species or path in self.fasta_species or label in seen_species:
                raise ValueError(f"Ambiguous or unknown species FASTA mapping: {label}: {path}")
            self.fasta_species[path] = label
            seen_species.add(label)
        if set(self.fasta_species) != set(fasta_paths) or set(self.fasta_species.values()) != set(species):
            raise ValueError("Species FASTA map does not match FASTA inputs and BUSCO summary species")

        self.codes = {label: args.genetic_code for label in species}
        if args.sequence_mode == "cds":
            from Bio.Data import CodonTable

            if args.species_genetic_codes:
                with args.species_genetic_codes.open(encoding="utf-8", newline="") as handle:
                    seen = set()
                    reader = csv.DictReader(handle, delimiter="\t")
                    if not {"species", "genetic_code"}.issubset(reader.fieldnames or []):
                        raise ValueError("Species genetic codes require species and genetic_code columns")
                    for line_number, row in enumerate(reader, start=2):
                        if None in row or not row["species"] or not row["genetic_code"]:
                            raise ValueError(f"Malformed species genetic code row {line_number}: {args.species_genetic_codes}")
                        label = row["species"]
                        if label in seen:
                            raise ValueError(f"Repeated species genetic code: {label}")
                        seen.add(label)
                        if label in self.codes:
                            try:
                                self.codes[label] = int(row["genetic_code"])
                            except ValueError as exc:
                                raise ValueError(f"Invalid species genetic code row {line_number}: {args.species_genetic_codes}") from exc
            if any(code not in CodonTable.unambiguous_dna_by_id for code in self.codes.values()):
                raise ValueError("Unsupported NCBI genetic code for BUSCO copy selection")

        tables = {}
        for path in sorted(args.busco_full_dir.glob("*.tsv")):
            if not path.is_file() or path.name.startswith("."):
                continue
            label = normalize_species_label(path.name)
            if label in tables:
                raise ValueError(f"Ambiguous BUSCO full tables for species: {label}")
            tables[label] = path
        self.table_paths = list(tables.values())
        for species_index, label in enumerate(species):
            if label not in tables:
                raise ValueError(f"Missing BUSCO full table for species: {label}")
            full = defaultdict(lambda: defaultdict(list))
            with tables[label].open(encoding="utf-8", newline="") as handle:
                for fields in csv.reader(handle, delimiter="\t"):
                    if fields and fields[0] and not fields[0].startswith("#"):
                        if len(fields) < 2 or fields[1] not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
                            raise ValueError(f"Malformed BUSCO full-table row: {tables[label]}: {fields[0]}")
                        if fields[1] == "Missing":
                            if len(fields) > 2 and fields[2] not in {"", "-"}:
                                raise ValueError(f"Malformed Missing BUSCO row: {tables[label]}: {fields[0]}")
                            continue
                        if len(fields) < 3 or not fields[2] or fields[2] == "-":
                            raise ValueError(f"BUSCO full-table row lacks a sequence ID: {tables[label]}: {fields[0]}")
                        full[fields[0]][source_gene_id(fields[2])].append(fields)
            if set(full) - set(summary):
                raise ValueError(f"BUSCO summary/full-table ID mismatch: {label}: markers absent from summary")
            for marker, cells in summary.items():
                token = cells[species_index]
                matches_by_gene = full.get(marker, {})
                if token == "-":
                    if matches_by_gene:
                        raise ValueError(f"BUSCO summary/full-table ID mismatch: {label}: {marker}")
                    continue
                identifiers = set(token.split(","))
                if not identifiers or "" in identifiers or identifiers != set(matches_by_gene):
                    raise ValueError(f"BUSCO summary/full-table ID mismatch: {label}: {marker}")
                duplicated = "," in token or any(row[1] == "Duplicated" for matches in matches_by_gene.values() for row in matches)
                if not duplicated:
                    self.singles[label][token] = self.singles[label].get(token, []) + [marker]
                    self.expected[marker] += 1
                    continue
                candidates = {}
                for gene in sorted(identifiers):
                    matches = matches_by_gene[gene]
                    eligible = [row for row in matches if row[1] in {"Complete", "Duplicated"}]
                    scores = [finite_score(row[3] if len(row) > 3 else "", f"{label}: {marker}: {gene}")
                              for row in eligible]
                    candidates[gene] = {
                        "busco_id": marker, "species": label, "gene_id": gene,
                        "busco_sequence_ids": ",".join(sorted({row[2] for row in matches})),
                        "status": ",".join(sorted({row[1] for row in matches})),
                        "score": max(scores) if scores else "",
                        "busco_length": ",".join(sorted({row[4] for row in matches if len(row) > 4})),
                        "sequence_mode": args.sequence_mode,
                        "genetic_code": self.codes[label] if args.sequence_mode == "cds" else "",
                    }
                self.groups[label][marker] = candidates
                if any(row["score"] != "" for row in candidates.values()):
                    self.expected[marker] += 1

    def extract(self, fasta_paths, fasta_records, pool):
        found = defaultdict(int)
        for path in fasta_paths:
            species = self.fasta_species[path]
            groups = self.groups[species]
            genes = defaultdict(list)
            for marker, candidates in groups.items():
                for gene in candidates:
                    genes[gene].append(marker)
            best_records = {}
            observed = set()
            for identifier, record in fasta_records(path):
                if identifier not in genes and identifier not in self.singles[species]:
                    continue
                if identifier in observed:
                    raise ValueError(f"Repeated candidate FASTA ID: {species}: {identifier}")
                observed.add(identifier)
                for marker in self.singles[species].get(identifier, ()):
                    pool.write(marker, [f">{species}\n", *record[1:]] if self.args.species_headers else record)
                    found[marker] += 1
                if identifier not in genes:
                    continue
                length, stops, ambiguous = sequence_quality(record, self.args.sequence_mode, self.codes[species])
                for marker in genes[identifier]:
                    candidate = groups[marker][identifier]
                    candidate.update(sequence_length=length, internal_stops=stops, ambiguous_fraction=ambiguous)
                    if candidate["score"] == "":
                        continue
                    rank = (stops > 0, -candidate["score"], ambiguous, identifier)
                    if marker not in best_records or rank < best_records[marker][0]:
                        best_records[marker] = rank, record
            missing = (set(genes) | set(self.singles[species])) - observed
            if missing:
                raise ValueError(f"BUSCO candidate IDs missing from {path}: {', '.join(sorted(missing)[:20])}")
            for marker, candidates in groups.items():
                ranked = sorted((row for row in candidates.values() if row["score"] != ""),
                                key=lambda row: (row["internal_stops"] > 0, -row["score"],
                                                 row["ambiguous_fraction"], row["gene_id"]))
                winner = ranked[0] if ranked else None
                runner = ranked[1] if len(ranked) > 1 else None
                for candidate in candidates.values():
                    selected = candidate is winner
                    if candidate["score"] == "":
                        reason = "not_complete"
                    elif winner["internal_stops"] == 0 and candidate["internal_stops"] > 0:
                        reason = "internal_stop"
                    elif selected:
                        if runner is None:
                            reason = "only_complete_candidate"
                        elif candidate["internal_stops"] == 0 and runner["internal_stops"] > 0:
                            reason = "no_internal_stop"
                        elif candidate["score"] != runner["score"]:
                            reason = "highest_score"
                        elif candidate["ambiguous_fraction"] != runner["ambiguous_fraction"]:
                            reason = "lowest_ambiguous_fraction"
                        else:
                            reason = "gene_id_tiebreak"
                    elif candidate["score"] < winner["score"]:
                        reason = "lower_score"
                    elif candidate["ambiguous_fraction"] > winner["ambiguous_fraction"]:
                        reason = "higher_ambiguous_fraction"
                    else:
                        reason = "gene_id_tiebreak"
                    candidate.update(selected=int(selected), reason=reason,
                                     selected_gene_id=winner["gene_id"] if winner else "",
                                     runner_up_gene_id=runner["gene_id"] if runner else "",
                                     score_margin=winner["score"] - runner["score"] if runner else "",
                                     selection_rule=SELECTION_RULE)
                    self.report_rows.append(candidate)
                if winner:
                    record = best_records[marker][1]
                    pool.write(marker, [f">{species}\n", *record[1:]] if self.args.species_headers else record)
                    found[marker] += 1
        return found


def write_selection_report(handle, rows):
    writer = csv.DictWriter(handle, fieldnames=SELECTION_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda row: (row["busco_id"], row["species"], row["gene_id"])))
