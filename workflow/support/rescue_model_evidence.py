#!/usr/bin/env python3
"""Audit frozen rescue models with optional target RNA, repeat and DNA evidence.

No model is rejected or rewritten by this CLI. Older rescue plans are read as
immutable evidence, never resumed using a different implementation identity.
"""
import argparse
import bisect
import csv
import gzip
import hashlib
import importlib.metadata
import json
import re
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    from fasta_sequence_store import open_text
    from gff_attribute_syntax import validate_attributes
    from input_generation_array_state import atomic_json, digest
    from rescue_gene_models import require_same_key, stage
    from rescue_model_quality import model_quality
except ImportError:
    from .fasta_sequence_store import open_text
    from .gff_attribute_syntax import validate_attributes
    from .input_generation_array_state import atomic_json, digest
    from .rescue_gene_models import require_same_key, stage
    from .rescue_model_quality import model_quality


def json_array(path, chunk_size=1024 * 1024):
    """Stream one JSON array without holding a whole-genome prediction set."""
    decoder = json.JSONDecoder()
    with open_text(Path(path)) as handle:
        buffer = ""
        eof = False
        def refill():
            nonlocal buffer, eof
            chunk = handle.read(chunk_size)
            buffer += chunk
            eof = not chunk
        def whitespace():
            nonlocal buffer
            while not buffer.strip() and not eof:
                refill()
            buffer = buffer.lstrip()
        whitespace()
        if not buffer.startswith("["):
            raise ValueError("Models must be a JSON array")
        buffer = buffer[1:]
        first = True
        while True:
            whitespace()
            if buffer.startswith("]"):
                buffer = buffer[1:]
                if buffer.strip() or handle.read().strip():
                    raise ValueError("Trailing data after models")
                return
            if not first:
                if not buffer.startswith(","):
                    raise ValueError("Missing model separator")
                buffer = buffer[1:]
                whitespace()
                if buffer.startswith("]"):
                    raise ValueError("Trailing model separator")
            while True:
                try:
                    value, end = decoder.raw_decode(buffer)
                    break
                except json.JSONDecodeError:
                    if eof:
                        raise ValueError("Truncated or invalid models JSON") from None
                    refill()
            if not isinstance(value, dict):
                raise ValueError("Model must be an object")
            yield value
            buffer = buffer[end:]
            first = False


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("Invalid half-open evidence interval")
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return result


def overlap_bases(blocks, intervals):
    merged = merge_intervals(intervals)
    total = 0
    for start, end, *_ in blocks:
        total += sum(max(0, min(end, right) - max(start, left)) for left, right in merged)
    return total


class IntervalIndex:
    def __init__(self, rows):
        self.data = {}
        for key, records in rows.items():
            records = sorted(records, key=lambda row: row[:2])
            maximum = 0
            prefix = []
            for row in records:
                maximum = max(maximum, row[1])
                prefix.append(maximum)
            self.data[key] = (records, [row[0] for row in records], prefix)

    def query(self, key, start, end):
        rows, starts, prefix = self.data.get(key, ([], [], []))
        index = bisect.bisect_left(starts, end) - 1
        result = []
        while index >= 0 and prefix[index] > start:
            if rows[index][1] > start:
                result.append(rows[index])
            index -= 1
        return result


def gff_rows(path, lengths):
    with open_text(Path(path)) as handle:
        for line in handle:
            if line.strip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or f[6] not in {"+", "-", "."}:
                raise ValueError("Invalid evidence GFF/GTF row")
            start, end = int(f[3]) - 1, int(f[4])
            if f[0] not in lengths or not 0 <= start < end <= lengths[f[0]]:
                raise ValueError("Evidence contig/coordinates disagree with genome")
            validate_attributes(f[8])
            yield f, start, end


def attributes(text):
    if re.match(r"^[^\s=;]+=", text.strip()):
        return {key: unquote(value) for token in text.rstrip(";").split(";") if "=" in token
                for key, value in [token.strip().split("=", 1)]}
    return dict(re.findall(r'(\w+)\s+"([^"\n]*)"', text))


def transcript_parents(text):
    if re.match(r"^[^\s=;]+=", text.strip()):
        raw = dict(token.strip().split("=", 1) for token in text.rstrip(";").split(";") if "=" in token)
        # Split encoded GFF lists before decoding literal commas in identifiers.
        return [unquote(value) for value in raw.get("Parent", "").split(",") if value]
    identifier = attributes(text).get("transcript_id")
    return [identifier] if identifier else []


def load_tracks(spec, lengths):
    junctions = defaultdict(set)
    transcripts = defaultdict(list)
    repeats = defaultdict(list)
    for track in spec.get("rna_junctions", []):
        for fields, start, end in gff_rows(track["path"], lengths):
            attr = attributes(fields[8])
            # Protein hints in BRAKER's mixed file are not RNA evidence.
            if fields[2] == "intron" and attr.get("src") == "E":
                junctions[fields[0], fields[6]].add((start, end))
    for number, track in enumerate(spec.get("rna_transcripts", [])):
        paths = defaultdict(list)
        owners = {}
        for fields, start, end in gff_rows(track["path"], lengths):
            if fields[2] != "exon":
                continue
            parents = transcript_parents(fields[8])
            if not parents:
                raise ValueError("RNA exon lacks transcript identity")
            for parent in parents:
                if parent in owners and owners[parent] != (fields[0], fields[6]):
                    raise ValueError("RNA transcript has conflicting contig/strand")
                owners[parent] = (fields[0], fields[6])
                paths[fields[0], fields[6], parent].append((start, end))
        for (seqid, strand, identifier), blocks in paths.items():
            blocks = sorted(set(blocks))
            if any(a[1] > b[0] for a, b in zip(blocks, blocks[1:], strict=False)):
                raise ValueError("Overlapping RNA exons in one transcript")
            transcripts[seqid, strand].append((blocks[0][0], blocks[-1][1], blocks,
                                               f"{number}:{identifier}", track["independence_group"]))
    for track in spec.get("repeats", []):
        with open_text(Path(track["path"])) as handle:
            for line in handle:
                fields = line.split()
                if not fields or not fields[0].isdigit():
                    continue
                if len(fields) < 14:
                    raise ValueError("Invalid RepeatMasker output row")
                seqid, start, end = fields[4], int(fields[5]) - 1, int(fields[6])
                if seqid not in lengths or not 0 <= start < end <= lengths[seqid]:
                    raise ValueError("Repeat coordinates disagree with genome")
                repeats[seqid].append((start, end, fields[10]))
    return junctions, IntervalIndex(transcripts), IntervalIndex(repeats)


def rna_evidence(model, spec, junctions, transcripts):
    blocks = sorted(model["cds"])
    introns = {(left[1], right[0]) for left, right in zip(blocks, blocks[1:], strict=False)}
    key = model["seqid"], model["strand"]
    supported = introns & junctions.get(key, set())
    unknown = introns & junctions.get((key[0], "."), set())
    compatible, unstranded = [], []
    for strand in (model["strand"], "."):
        for _start, _end, exons, identifier, group in transcripts.query((key[0], strand), blocks[0][0], blocks[-1][1]):
            chain_introns = {(a[1], b[0]) for a, b in zip(exons, exons[1:], strict=False)
                             if a[1] < blocks[-1][1] and b[0] > blocks[0][0]}
            if (chain_introns == introns and all(any(a <= start < end <= b for a, b in exons)
                                                for start, end, *_ in blocks)):
                (compatible if strand != "." else unstranded).append((identifier, group))
    available = bool(spec.get("rna_junctions") or spec.get("rna_transcripts"))
    return {"status": "available" if available else "not_provided", "junction_count": len(introns),
            "junctions_supported": len(supported) if spec.get("rna_junctions") else None,
            "unstranded_junctions": len(unknown) if spec.get("rna_junctions") else None,
            "exon_chain_transcripts": [item[0] for item in compatible] if spec.get("rna_transcripts") else None,
            "exon_chain_independence_groups": sorted({item[1] for item in compatible}) if spec.get("rna_transcripts") else None,
            "unstranded_exon_chain_transcripts": [item[0] for item in unstranded] if spec.get("rna_transcripts") else None,
            "translation_start_confirmed": False}


def repeat_evidence(model, spec, repeats):
    if not spec.get("repeats"):
        return {"status": "not_provided", "masked_fraction": None, "te_annotated_fraction": None, "classes": None}
    rows = repeats.query(model["seqid"], min(b[0] for b in model["cds"]), max(b[1] for b in model["cds"]))
    total = sum(end - start for start, end, *_ in model["cds"])
    te = {"LINE", "SINE", "LTR", "DNA", "RC", "Retroposon", "PLE"}
    return {"status": "available", "masked_fraction": overlap_bases(model["cds"], [row[:2] for row in rows]) / total,
            "te_annotated_fraction": overlap_bases(model["cds"], [row[:2] for row in rows if row[2].split("/")[0] in te]) / total,
            "classes": sorted({row[2] for row in rows if overlap_bases(model["cds"], [row[:2]])})}


def dna_evidence(model, spec, bam, genome):
    if bam is None:
        return {"status": "not_provided", "hq_depth_min": None, "hq_depth_median": None, "start_base_support": None}
    dna = spec["dna"]
    def usable(read):
        return (not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail)
                and read.mapping_quality >= dna["min_mapq"])
    depths = []
    # A split initiator may cross an exon boundary. Follow spliced CDS order.
    positions = []
    for start, end, *_ in sorted(model["cds"], reverse=model["strand"] == "-"):
        positions.extend(list(range(start, min(end, start + 3))) if model["strand"] == "+" else
                         list(range(end - 1, max(start - 1, end - 4), -1)))
        if len(positions) >= 3:
            positions = positions[:3]
            break
    calls = {}
    for start, end, *_ in model["cds"]:
        counts = bam.count_coverage(model["seqid"], start, end, quality_threshold=dna["min_baseq"], read_callback=usable)
        depths.extend(map(sum, zip(*counts, strict=True)))
        for position in positions:
            if start <= position < end:
                base = genome.fetch(model["seqid"], position, position + 1).upper()
                calls[position] = {"genomic_position_1based": position + 1, "genomic_base": base,
                                   "matching": counts["ACGT".index(base)][position - start],
                                   "total": sum(array[position - start] for array in counts)} if base in "ACGT" else None
    first_counts = [calls[position] for position in positions]
    return {"status": "available", "hq_depth_min": min(depths), "hq_depth_median": statistics.median(depths),
            "start_base_support": first_counts, "confirms_translation_initiation": False,
            "reference_identity": "manifest_asserted_with_exact_bam_contig_lengths"}


def check_model(model, genome, genetic_code=1):
    if model.get("strand") not in {"+", "-"} or not model.get("cds"):
        raise ValueError("Accepted model lacks CDS/strand")
    pieces = []
    cumulative = 0
    blocks = sorted(model["cds"], reverse=model["strand"] == "-")
    previous = None
    for start, end, phase in blocks:
        if (any(type(value) is not int for value in (start, end, phase))
                or not 0 <= start < end <= genome.get_reference_length(model["seqid"]) or phase != (-cumulative) % 3):
            raise ValueError("Accepted CDS coordinates/phase disagree with genome")
        if previous is not None and ((model["strand"] == "+" and start < previous[1])
                                     or (model["strand"] == "-" and end > previous[0])):
            raise ValueError("Accepted CDS blocks overlap")
        piece = genome.fetch(model["seqid"], start, end).upper()
        pieces.append(str(Seq(piece).reverse_complement()) if model["strand"] == "-" else piece)
        cumulative += end - start
        previous = (start, end)
    if "".join(pieces) != model["sequence"] or len(model["sequence"]) % 3:
        raise ValueError("Accepted model sequence disagrees with genomic CDS")
    sequence = model["sequence"]
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    if (set(sequence) - set("ACGT") or sequence[:3] not in table.start_codons
            or sequence[-3:] not in table.stop_codons or "*" in str(Seq(sequence).translate(table=genetic_code))[:-1]):
        raise ValueError("Accepted model lacks an intact code-compatible ORF")


def manifest_spec(path, name, genome_hash):
    if path is None:
        return {}
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != 1 or name not in manifest.get("species", {}):
        raise ValueError("Evidence manifest needs schema_version=1 and target species")
    spec = manifest["species"][name]
    if set(spec) - {"reference_genome_sha256", "rna_junctions", "rna_transcripts", "repeats", "dna"}:
        raise ValueError("Unknown evidence manifest fields")
    if spec.get("reference_genome_sha256") != genome_hash:
        raise ValueError("Evidence manifest reference differs from frozen genome")
    for kind in ("rna_junctions", "rna_transcripts", "repeats"):
        formats = {"rna_junctions": "braker_hints", "rna_transcripts": "exon_gff_gtf", "repeats": "repeatmasker_out"}
        if not isinstance(spec.get(kind, []), list):
            raise ValueError("Evidence tracks must be a list")
        seen_paths = set()
        for track in spec.get(kind, []):
            if not isinstance(track, dict) or set(track) - {"path", "format", "independence_group", "label", "tissue"}:
                raise ValueError("Unknown evidence track fields")
            if track.get("format") != formats[kind] or not Path(track["path"]).is_absolute():
                raise ValueError("Evidence track needs supported format and absolute path")
            if str(Path(track["path"]).resolve()) in seen_paths:
                raise ValueError("Duplicate evidence track path")
            seen_paths.add(str(Path(track["path"]).resolve()))
            if kind.startswith("rna") and not track.get("independence_group"):
                raise ValueError("RNA tracks require an explicit independence_group")
    if spec.get("dna"):
        dna = spec["dna"]
        if not isinstance(dna, dict) or set(dna) - {"path", "index", "format", "min_mapq", "min_baseq"}:
            raise ValueError("Unknown DNA evidence fields")
        if dna.get("format") != "bam" or not all(Path(dna[key]).is_absolute() for key in ("path", "index")):
            raise ValueError("DNA evidence needs BAM and explicit absolute index")
        if any(type(dna.get(key)) is not int or dna[key] < 0 for key in ("min_mapq", "min_baseq")):
            raise ValueError("DNA quality thresholds must be nonnegative integers")
    return spec


def audit(args):
    import pysam
    root = args.rescue_output.resolve()
    destination = args.output.resolve()
    if destination == root or root in destination.parents or destination in root.parents:
        raise ValueError("Audit output must be separate from the frozen rescue directory")
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    source = plan["request"]["sources"][args.species]
    worker = root / "rescued" / args.species
    models_path = worker / "models.json"
    receipt_path = worker / "receipt.json"
    frozen = json.loads(receipt_path.read_text())
    if frozen.get("key", {}).get("plan") != digest(plan_path):
        raise ValueError("Models receipt belongs to another plan")
    if frozen.get("files", {}).get("models.json") != digest(models_path):
        raise ValueError("Frozen models changed")
    genome_hash = plan["request"]["files"][source["genome"]]
    if digest(source["genome"]) != genome_hash:
        raise ValueError("Frozen genome changed")
    spec = manifest_spec(args.evidence_manifest, args.species, genome_hash)
    files = [plan_path, models_path, receipt_path, Path(source["genome"]), Path(__file__),
             *{Path(function.__code__.co_filename) for function in (model_quality, open_text, stage, atomic_json, validate_attributes)}]
    if args.evidence_manifest:
        files.append(args.evidence_manifest)
    for kind in ("rna_junctions", "rna_transcripts", "repeats"):
        files.extend(Path(track["path"]) for track in spec.get(kind, []))
    if spec.get("dna"):
        files.extend(Path(spec["dna"][key]) for key in ("path", "index"))
    if any(destination == p.resolve() or destination in p.resolve().parents for p in files):
        raise ValueError("Audit output must not contain a source file")
    def current_key():
        return {"schema": 1, "species": args.species, "inputs": {str(p.resolve()): digest(p) for p in files},
                "biopython": importlib.metadata.version("biopython"), "pysam": pysam.__version__}
    key = current_key()
    def build(tmp):
        # Only the genome needs indexing; preserve inputs and avoid loading every
        # source annotation or rerunning any protein alignment.
        genome_path = tmp / "reference.fa"
        if source["genome"].endswith(".gz"):
            with gzip.open(source["genome"], "rb") as incoming, genome_path.open("wb") as out:
                shutil.copyfileobj(incoming, out)
        else:
            shutil.copyfile(source["genome"], genome_path)
        names = set()
        with genome_path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    fields = line[1:].split()
                    if not fields or fields[0] in names:
                        raise ValueError("Empty or duplicate genome reference")
                    names.add(fields[0])
        if not names:
            raise ValueError("Genome has no references")
        pysam.faidx(str(genome_path))
        bam = None
        with pysam.FastaFile(str(genome_path)) as genome:
            lengths = dict(zip(genome.references, genome.lengths, strict=True))
            if len(lengths) != len(genome.references):
                raise ValueError("Duplicate genome references")
            junctions, transcripts, repeats = load_tracks(spec, lengths)
            if spec.get("dna"):
                bam = pysam.AlignmentFile(spec["dna"]["path"], "rb", index_filename=spec["dna"]["index"])
                if dict(zip(bam.references, bam.lengths, strict=True)) != lengths or not bam.has_index():
                    bam.close()
                    raise ValueError("DNA BAM reference/index disagrees with genome")
            records = []
            seen = set()
            try:
                for model in json_array(models_path):
                    if model.get("status") != "accepted":
                        continue
                    if model["model_id"] in seen:
                        raise ValueError("Duplicate accepted model identity")
                    seen.add(model["model_id"])
                    check_model(model, genome, source["genetic_code"])
                    quality = model_quality(model, source["genetic_code"])
                    rna = rna_evidence(model, spec, junctions, transcripts)
                    repeat = repeat_evidence(model, spec, repeats)
                    dna = dna_evidence(model, spec, bam, genome)
                    flags = set(quality["flags"])
                    if rna["exon_chain_transcripts"] == []:
                        flags.add("no_stranded_rna_exon_chain_detected")
                    if rna["junctions_supported"] is not None and rna["junctions_supported"] < rna["junction_count"]:
                        flags.add("incomplete_rna_junction_support")
                    if repeat["te_annotated_fraction"] is not None and repeat["te_annotated_fraction"] >= 0.5:
                        flags.add("te_overlap_ge_50pct")
                    if dna["hq_depth_min"] == 0:
                        flags.add("hq_dna_uncovered_cds_bases")
                    records.append({"model_id": model["model_id"], "seqid": model["seqid"], "strand": model["strand"],
                                    "cds": model["cds"], "cds_sha256": hashlib.sha256(model["sequence"].encode()).hexdigest(),
                                    "rescue_status": "accepted", "quality": quality, "rna": rna, "repeat": repeat, "dna": dna,
                                    "review_flags": sorted(flags), "decision": "unchanged", "functional_status": "not_established"})
            finally:
                if bam:
                    bam.close()
        atomic_json(tmp / "evidence.json", records)
        summary = {"species": args.species, "accepted_models": len(records), "sequence_changes": 0,
                   "flag_counts": dict(sorted(Counter(flag for row in records for flag in row["review_flags"]).items())),
                   "evidence_available": {kind: bool(spec.get(kind)) for kind in ("rna_junctions", "rna_transcripts", "repeats", "dna")},
                   "policy": "Advisory only; missing RNA, TE overlap and alternative starts are not automatic rejection criteria."}
        atomic_json(tmp / "summary.json", summary)
        with (tmp / "evidence.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["model_id", "seqid", "strand", "start_codon", "donor_species_count", "n_aligned", "c_aligned",
                             "rna_junctions_supported", "rna_junction_count", "rna_chain_groups", "te_fraction", "hq_dna_min", "hq_dna_median", "flags"])
            for row in records:
                q, r, rep, d = (row[k] for k in ("quality", "rna", "repeat", "dna"))
                writer.writerow([row["model_id"], row["seqid"], row["strand"], q["start_codon"], len(q["donor_species"]),
                                 q["terminal_alignment"]["n_aligned"], q["terminal_alignment"]["c_aligned"], r["junctions_supported"],
                                 r["junction_count"], None if r["exon_chain_independence_groups"] is None else len(r["exon_chain_independence_groups"]),
                                 rep["te_annotated_fraction"], d["hq_depth_min"], d["hq_depth_median"], ",".join(row["review_flags"])])
        genome_path.unlink()
        Path(str(genome_path) + ".fai").unlink()
    stage(destination.parent, Path(destination.name), key, build, lambda: require_same_key(key, current_key()))
    return json.loads((destination / "summary.json").read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescue-output", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--evidence-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args), sort_keys=True))


if __name__ == "__main__":
    main()
