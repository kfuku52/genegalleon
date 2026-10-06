"""Evidence-based normalisation of existing annotated CDS models.

Original FASTA/GFF records are never edited. Unsupported translations are
withheld; genomic reconstructions are used only with same-model evidence.
"""
import contextlib
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio.Seq import Seq

try:
    from performance_metrics import measure
except ImportError:
    from .performance_metrics import measure

try:
    from fasta_sequence_store import open_text
except ImportError:
    from .fasta_sequence_store import open_text


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def attributes(text):
    return {key: unquote(value) for item in text.split(";") if "=" in item
            for key, value in [item.split("=", 1)]}


def tokens(value):
    return {item.replace("GeneID:", "GeneID", 1) if item.startswith("GeneID:") else item
            for item in value.split(",") if item}


def phase_conventions(blocks):
    """Distinguish GFF3 skip-count phases from complementary frame fields."""
    if not blocks or any(row["phase"] not in {"0", "1", "2"} for row in blocks):
        return set()
    first, cumulative = int(blocks[0]["phase"]), 0
    conventions = {"gff3", "complementary"}
    for row in blocks:
        phase = int(row["phase"])
        if phase != (first - cumulative) % 3:
            conventions.discard("gff3")
        if phase != (first + cumulative) % 3:
            conventions.discard("complementary")
        cumulative += row["end"] - row["start"]
    return conventions


def protein(sequence, code, *, partial=False):
    if partial:
        sequence = sequence[:len(sequence) // 3 * 3]
    if not sequence or len(sequence) % 3:
        return None
    translated = str(Seq(sequence).translate(table=code)).removesuffix("*")
    return translated if translated and "*" not in translated else None


def agrees_with_padding(supplied, genomic):
    return supplied == genomic or any(supplied == genomic + "N" * n for n in (1, 2))


def ambiguity_to_n(sequence):
    # Legacy formatters collapse ambiguity codes to N. This is an identity
    # comparison only; neither the sequence nor its translation is masked.
    return "".join(base if base in "ACGTN" else "N" for base in sequence)


class CdsModelNormaliser:
    def __init__(self, source, directory, side, *, genome_records=None, reference_mapping=None,
                 gff_lines=None, attribute_parser=None, coding_only=False):
        self.source, self.directory, self.side = source, Path(directory), side
        self.genome_records, self.reference_mapping = genome_records, reference_mapping
        self.corrected = {}
        self.rows, self.models, self.nodes = [], defaultdict(list), defaultdict(list)
        self.exons, self.features = defaultdict(list), []
        self.lookup, self.genome, self.scratch = None, None, None
        self._genome_context = None
        self.reference_lengths = None
        wanted, retained = set(), set()
        with (contextlib.nullcontext(gff_lines()) if gff_lines is not None else open_text(Path(source["gff"]))) as handle:
            for line in handle:
                if line.strip() == "##FASTA":
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n\r").split("\t")
                if len(fields) != 9:
                    raise ValueError("Invalid GFF row in anchor admission")
                if coding_only and fields[2] not in {"CDS", "exon"}:
                    continue
                attr = (attribute_parser or attributes)(fields[8])
                row = {"seqid": fields[0], "source": fields[1], "feature": fields[2], "start": int(fields[3]) - 1,
                       "end": int(fields[4]), "strand": fields[6], "phase": fields[7], "attributes": attr,
                       "source_line": line.rstrip("\r\n")}
                self.features.append(row)
                if attr.get("ID"):
                    self.nodes[attr["ID"]].append(row)
                wanted.update(filter(None, attr.get("Parent", "").split(",")))
                if fields[2] in {"CDS", "exon"}:
                    parents = attr.get("Parent", attr.get("transcript_id", attr.get("ID", ""))).split(",")
                    for parent in filter(None, parents):
                        (self.models if fields[2] == "CDS" else self.exons)[parent].append(row)
        # Formatting uses CDS models and their declared ancestors, not millions
        # of independent alignment features. Resolve the full ancestor closure
        # with bounded passes, including unusual feature types and duplicate IDs.
        # The default retains all features for arbitrary anchor mappings.
        while coding_only and wanted - retained:
            pending = wanted - retained
            discovered = set()
            with (contextlib.nullcontext(gff_lines()) if gff_lines is not None else open_text(Path(source["gff"]))) as handle:
                for line in handle:
                    if line.strip() == "##FASTA":
                        break
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n\r").split("\t")
                    if len(fields) != 9:
                        raise ValueError("Invalid GFF row in anchor admission")
                    attr = (attribute_parser or attributes)(fields[8])
                    identifier = attr.get("ID")
                    if identifier not in pending:
                        continue
                    discovered.add(identifier)
                    wanted.update(filter(None, attr.get("Parent", "").split(",")))
                    if fields[2] in {"CDS", "exon"}:
                        continue  # Already retained during the initial pass.
                    row = {"seqid": fields[0], "source": fields[1], "feature": fields[2], "start": int(fields[3]) - 1,
                           "end": int(fields[4]), "strand": fields[6], "phase": fields[7], "attributes": attr,
                           "source_line": line.rstrip("\r\n")}
                    self.features.append(row)
                    self.nodes[identifier].append(row)
            # Missing parents retain the default resolver's behavior; never spin.
            retained.update(pending)
            if not discovered:
                break
        self.votes = Counter()
        for blocks in self.models.values():
            blocks.sort(key=lambda row: row["start"], reverse=blocks[0]["strand"] == "-")
            conventions = phase_conventions(blocks)
            if len(conventions) == 1 and all(row["source"] != "genegalleon_rescue" for row in blocks):
                self.votes.update(conventions)
        self.global_convention = next(iter(self.votes)) if len(self.votes) == 1 else None

    def ancestry(self, identifiers):
        found, todo = set(), list(identifiers)
        while todo:
            item = todo.pop()
            if item in found:
                continue
            found.add(item)
            for row in self.nodes.get(item, []):
                todo.extend(filter(None, row["attributes"].get("Parent", "").split(",")))
        return found

    def configure(self, mapping):
        aliases = defaultdict(set)
        for row in self.features:
            if row["feature"] != mapping.feature:
                continue
            attr = row["attributes"]
            roots = (set(filter(None, attr.get("Parent", "").split(",")))
                     if mapping.attribute == "Parent" else
                     {attr["ID"]} if attr.get("ID") else set(filter(None, attr.get("Parent", "").split(","))))
            for alias in tokens(attr.get(mapping.attribute, "")):
                aliases[alias].update(roots)
        by_root = defaultdict(set)
        for key, blocks in self.models.items():
            members = self.ancestry([key])
            members.update(row["attributes"]["ID"] for row in blocks if row["attributes"].get("ID"))
            for member in members:
                by_root[member].add(key)
        self.lookup = {alias: set().union(*(by_root[root] for root in roots)) for alias, roots in aliases.items()}

    def fetch(self, rows):
        if self.genome is None:
            with measure("genome_reconstruction_index"):
                self._open_genome()
        if len({(row["seqid"], row["strand"]) for row in rows}) != 1 or rows[0]["strand"] not in {"+", "-"}:
            raise ValueError("Inconsistent annotation strand/contig in anchor genome reconstruction")
        rows = sorted(rows, key=lambda row: row["start"], reverse=rows[0]["strand"] == "-")
        ascending = sorted(rows, key=lambda row: row["start"])
        if any(a["end"] > b["start"] for a, b in zip(ascending, ascending[1:], strict=False)):
            raise ValueError("Overlapping annotation blocks in anchor genome reconstruction")
        parts = []
        for row in rows:
            contig = (self.reference_mapping() if self.reference_mapping is not None else {}).get(row["seqid"], row["seqid"])
            length = self.reference_lengths.get(contig)
            if (length is None or row["start"] < 0 or row["end"] <= row["start"]
                    or row["end"] > length):
                raise ValueError("Annotation coordinates outside anchor genome: " + contig)
            sequence = self.genome.fetch(contig, row["start"], row["end"]).upper()
            parts.append(str(Seq(sequence).reverse_complement()) if row["strand"] == "-" else sequence)
        return "".join(parts)

    def _open_genome(self):
        import pysam
        if self.genome_records is None:
            from gene_model_catalog import indexed_genome
            self._genome_context = indexed_genome(Path(self.source["genome"]))
            self.genome = self._genome_context.__enter__()
            self.reference_lengths = dict(zip(self.genome.references, self.genome.lengths, strict=True))
            return
        self.scratch = tempfile.TemporaryDirectory(prefix=".anchor-genome-", dir=self.directory)
        path = Path(self.scratch.name) / "genome.fa"
        with path.open("w") as handle:
            for identifier, sequence in self.genome_records():
                handle.write(f">{identifier}\n{sequence}\n")
        result = subprocess.run([sys.executable, "-c", "import pysam,sys; pysam.faidx(sys.argv[1])", str(path)],
                                capture_output=True, text=True)
        if result.returncode or result.stderr.strip():
            raise ValueError("FASTA index warning or failure: " + result.stderr.strip())
        self.genome = pysam.FastaFile(str(path))
        self.reference_lengths = dict(zip(self.genome.references, self.genome.lengths, strict=True))

    def close(self):
        if self._genome_context is not None:
            self._genome_context.__exit__(None, None, None)
            self._genome_context = None
        elif self.genome is not None:
            self.genome.close()
        if self.scratch is not None:
            self.scratch.cleanup()

    def exceptional(self, key):
        rows = [*self.models[key], *(row for owner in self.ancestry([key]) for row in self.nodes.get(owner, []))]
        for row in rows:
            attr = {name.lower(): value for name, value in row["attributes"].items()}
            if (row["feature"].lower() == "pseudogene" or
                    any(attr.get(name, "").lower() not in {"", "false", "0", "no"} for name in ("pseudo", "pseudogene"))):
                return "annotated_pseudogene"
            if attr.get("transl_except") or attr.get("exception"):
                return "annotated_translation_exception"
            note = attr.get("note", "").lower()
            if "inserted" in note or "deleted" in note or "missing" in note:
                return "annotated_sequence_exception"
        return None

    def normalise(self, identifier, sequence, keys):
        keys = sorted(keys)
        code = self.source["genetic_code"]
        original = protein(sequence, code)
        record = {"original_id": identifier, "status": "excluded", "reason": "unsupported_translation",
                  "input_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "genetic_code": code,
                  "transcripts": keys, "evidence": []}
        self.rows.append(record)
        # Annotation exceptions are never translated by deleting or masking stops.
        exceptions = [self.exceptional(key) for key in keys]
        exception = next((value for value in exceptions if value), None)
        if exception:
            record["reason"] = exception
            return None
        if original:
            record.update(status="unchanged", reason="valid_original_translation",
                          protein_sha256=hashlib.sha256(original.encode()).hexdigest())
            return original
        compatible = []
        for key in keys:
            blocks = self.models[key]
            ascending = sorted(blocks, key=lambda row: row["start"])
            if (len({(row["seqid"], row["strand"]) for row in blocks}) != 1
                    or blocks[0]["strand"] not in {"+", "-"}
                    or any(a["end"] > b["start"] for a, b in zip(ascending, ascending[1:], strict=False))):
                record["evidence"].append({"transcript": key, "reason": "unsupported_cds_structure"})
            else:
                compatible.append(key)
        keys = compatible
        genomic_by_key = {key: self.fetch(self.models[key]) for key in keys}
        bound = [key for key in keys if agrees_with_padding(sequence, genomic_by_key[key])]
        # A supplied CDS matching an annotated CDS is bound to that model.
        # A shorter ORF in another isoform must not relabel its disrupted bases
        # as UTR and turn a genuine internal stop into a technical correction.
        record["bound_transcripts"] = bound
        if bound:
            keys = bound
        alternatives = {}
        for key in keys:
            blocks = self.models[key]
            genomic = genomic_by_key[key]
            first = int(blocks[0]["phase"]) if blocks[0]["phase"] in {"0", "1", "2"} else None
            conventions = phase_conventions(blocks)
            convention = next(iter(conventions)) if len(conventions) == 1 else (
                self.global_convention if self.global_convention in conventions else None)
            evidence = {"transcript": key, "first_phase": first, "phase_conventions": sorted(conventions)}
            record["evidence"].append(evidence)
            if agrees_with_padding(sequence, genomic):
                if first is None or not conventions:
                    evidence["reason"] = "inconsistent_or_missing_phase"
                    continue
                if first and convention is None:
                    evidence["reason"] = "ambiguous_phase_convention"
                    continue
                offset = (3 - first) % 3 if convention == "complementary" else first
                if any("," in row["attributes"].get("Parent", "") for row in blocks) and first:
                    evidence["reason"] = "shared_partial_cds_blocks"
                    continue
                translated = protein(genomic[offset:], code, partial=bool(first))
                if translated and (first or sequence != genomic):
                    evidence.update(reason="annotated_partial_frame" if first else "formatter_terminal_padding",
                                    phase_convention=convention, offset=offset,
                                    trailing_bases_omitted=(len(genomic) - offset) % 3 if first else 0)
                    repaired = genomic[offset:len(genomic) - evidence["trailing_bases_omitted"]]
                    alternatives.setdefault((translated, repaired), []).append(evidence)
                else:
                    evidence["reason"] = "genomic_internal_stop_or_incomplete_codon"
                continue
            exons = self.exons.get(key, [])
            if not exons or first != 0 or not conventions:
                evidence["reason"] = "input_genomic_cds_mismatch"
                continue
            transcript = self.fetch(exons)
            supplied = sequence
            # At most two formatter-added Ns can be ignored, and only when
            # their removal proves identity to this annotated transcript.
            matches = [supplied[:-n] if n else supplied for n in (0, 1, 2)
                       if (not n or supplied.endswith("N" * n))
                       and ambiguity_to_n(supplied[:-n] if n else supplied) in ambiguity_to_n(transcript)]
            translated = protein(genomic, code)
            if (translated and matches and genomic in transcript and
                    all(ambiguity_to_n(genomic) in ambiguity_to_n(value) for value in matches)
                    and all(len(value) > len(genomic) for value in matches)):
                evidence.update(reason="annotated_utr_in_cds", offset=0,
                                genomic_cds_sha256=hashlib.sha256(genomic.encode()).hexdigest())
                alternatives.setdefault((translated, genomic), []).append(evidence)
            else:
                evidence["reason"] = "input_genomic_cds_mismatch"
        if len(alternatives) == 1:
            (translated, repaired), evidence = next(iter(alternatives.items()))
            self.corrected[identifier] = {"sequence": repaired, "evidence": evidence}
            record["normalised_cds_sha256"] = hashlib.sha256(repaired.encode()).hexdigest()
            record.update(status="normalised", reason=evidence[0]["reason"],
                          protein_sha256=hashlib.sha256(translated.encode()).hexdigest(), selected_evidence=evidence)
            return translated
        record["reason"] = "ambiguous_genomic_translations" if alternatives else (
            record["evidence"][0]["reason"] if record["evidence"] else "no_bound_cds_model")
        return None

    def audit(self):
        summary = {"policy": "annotated_existing_anchor_v1", "original_records_modified": False,
                   "counts": dict(Counter(row["status"] for row in self.rows)),
                   "reasons": dict(Counter(row["reason"] for row in self.rows)),
                   "phase_convention_votes": dict(self.votes), "file_phase_convention": self.global_convention,
                   "records": [row for row in self.rows if row["status"] != "unchanged"]}
        (self.directory / f"{self.side}.anchor_admission.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        columns = ("original_id", "status", "reason", "genetic_code", "input_sha256", "protein_sha256")
        write_tsv(self.directory / f"{self.side}.anchor_admission.tsv", columns,
                  ([row.get(column, "") for column in columns] for row in self.rows))
        return {key: value for key, value in summary.items() if key != "records"}
