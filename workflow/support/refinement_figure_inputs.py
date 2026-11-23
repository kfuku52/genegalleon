"""Read publication-bound inputs for the compact refinement figure."""

import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

from content_digest_cache import ContentDigestCache, cache
from fasta_sequence_store import open_text


def cached_locus_count(path, namespace, collect):
    path = Path(path)
    saved = cache().get(namespace, path)
    if saved is not None:
        return saved[1]
    signature = ContentDigestCache.signature(path)
    count = collect(path)
    if ContentDigestCache.signature(path) != signature:
        raise ValueError(f"Annotation changed while counting loci: {path}")
    cache().put(namespace, path, hashlib.sha256(str(count).encode()).hexdigest(), count,
                expected_signature=signature)
    return count


def _coding_loci(path):
    genes, coding, cds_parents = set(), set(), set()
    parents = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if len(fields) != 9:
                continue
            kind = fields[2]
            if kind not in ("gene", "CDS") and "RNA" not in kind and "transcript" not in kind:
                continue
            attributes = dict(part.split("=", 1) for part in fields[8].strip().split(";") if "=" in part)
            ident = attributes.get("ID")
            parent = attributes.get("Parent")
            if kind == "gene" and ident:
                genes.add(ident)
                if any(attributes.get(key) == "protein_coding" for key in ("gene_biotype", "gene_type", "biotype")):
                    coding.add(ident)
            elif kind == "CDS" and parent:
                cds_parents.update(parent.split(","))
            elif ident and parent:
                parents.setdefault(ident, set()).update(parent.split(","))
    pending = list(cds_parents)
    seen = set()
    while pending:
        parent = pending.pop()
        if parent in seen:
            continue
        seen.add(parent)
        if parent in genes:
            coding.add(parent)
        else:
            pending.extend(parents.get(parent, ()))
    return len(coding)


def _transcriptome_loci(path):
    with open_text(path, newline="") as handle:
        return len({row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")})


def _fasta_loci(path):
    with open_text(path) as handle:
        return sum(line.startswith(">") for line in handle)


def locus_totals(name, request, original_sources, source, entry, sequence_dir):
    assembly = request.get("transcriptome_sources", {}).get(name)
    if assembly:
        return {
            "input": cached_locus_count(assembly["transcript_gene_map"], "refinement_transcriptome_loci_v1",
                                       _transcriptome_loci),
            "refined": cached_locus_count(sequence_dir / entry["cds"], "refinement_selected_loci_v1", _fasta_loci),
        }
    if not source or not source.get("gff") or not entry.get("gff"):
        return {}
    input_gff = original_sources.get(name, source)["gff"]
    stages = {"input": input_gff, "refined": entry["gff"]}
    if name in original_sources:
        stages["rescued"] = source["gff"]
    return {stage: cached_locus_count(path, "refinement_coding_loci_v1", _coding_loci)
            for stage, path in stages.items()}


def published_member(effective, files, name):
    data = (effective / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != files[name]:
        raise ValueError(f"Refinement reporting input changed: {name}")
    return data


def refinement_annotations(effective, publication, analysed):
    files = publication["files"]
    if "utr_summary.tsv" not in files:
        return {}
    utr = csv.DictReader(io.StringIO(published_member(effective, files, "utr_summary.tsv").decode()), delimiter="\t")
    events = {"merges": Counter(), "splits": Counter()}
    for key, filename in (("merges", "split_model_merges.json"), ("splits", "fused_model_splits.json")):
        if filename in files:
            events[key] = Counter(row["species"] for row in json.loads(published_member(effective, files, filename)))
    repairs = Counter()
    if "model_repairs.json" in files:
        repairs.update((row["species"], bool(row["representative"]))
                       for row in json.loads(published_member(effective, files, "model_repairs.json")))
    counts = {}
    for row in utr:
        species = row["species"]
        if species not in analysed:
            continue
        total = int(row["exported_transcripts"])
        if total <= 0:
            continue
        raw = published_member(effective, files, f"utr_audit/{species}.json.gz")
        sides, origins = Counter(), Counter()
        for record in json.loads(gzip.decompress(raw)):
            ends = [key for key in ("utr5", "utr3") if record.get(key, {}).get("length")]
            if not ends:
                continue
            sources = {record[key]["source"] for key in ends}
            if not sources <= {"oriongeno", "existing", "selected"}:
                raise ValueError(f"Unknown UTR source for {species}: {sources}")
            sides["both" if len(ends) == 2 else ("five_only" if ends[0] == "utr5" else "three_only")] += 1
            origins["both" if len(sources) > 1 and "oriongeno" in sources
                    else "oriongeno_only" if sources == {"oriongeno"} else "original_gff_only"] += 1
        if (sides["five_only"] + sides["both"] != int(row["with_5prime_utr"])
                or sides["three_only"] + sides["both"] != int(row["with_3prime_utr"])
                or sides["both"] != int(row["with_both_utrs"])):
            raise ValueError(f"UTR counts differ from the published summary for {species}")
        sides["none"] = origins["none"] = total - sum(sides.values())
        counts[species] = {
            "repairs": (repairs[species, True], repairs[species, False]),
            "utr": {"final_genes": total, "sides": dict.fromkeys(("five_only", "three_only", "both", "none"), 0) | sides,
                    "origins": dict.fromkeys(("original_gff_only", "oriongeno_only", "both", "none"), 0) | origins},
            "architecture": {key: counter[species] for key, counter in events.items()},
        }
    return counts
