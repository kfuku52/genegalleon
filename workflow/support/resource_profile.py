#!/usr/bin/env python3
"""Build reusable, content-bound input-size evidence before scheduler submission.

Inputs are explicit: a profiler must not silently select scientific inputs or
count every historical workspace file as unfinished work. Only the Python
standard library is needed; gzip FASTA and FASTQ are streamed.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

WORKFLOWS = {
    "gg_genome_evolution", "gg_gene_evolution", "gg_genome_annotation",
    "gg_transcriptome_generation", "gg_input_generation",
    "gg_fractionation_bias", "gg_gene_summary", "gg_progress_summary",
}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def scan(path: Path, kind: str) -> dict:
    before = path.stat()
    checksum = hashlib.sha256()
    sequences = residues = longest = current = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        if kind == "fastq":
            while True:
                header = handle.readline()
                if not header:
                    break
                seq, plus, quality = (handle.readline() for _ in range(3))
                if (not header.startswith(b"@") or not plus.startswith(b"+")
                        or not seq.strip() or len(seq.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n"))):
                    raise ValueError(f"invalid four-line FASTQ: {path}")
                for line in (header, seq, plus, quality):
                    checksum.update(line)
                length = len(seq.rstrip(b"\r\n"))
                sequences += 1
                residues += length
                longest = max(longest, length)
        else:
            for line in handle:
                checksum.update(line)
                if line.startswith(b">"):
                    longest = max(longest, current)
                    current = 0
                    sequences += 1
                elif line.strip():
                    if not sequences:
                        raise ValueError(f"FASTA sequence before header: {path}")
                    length = len(b"".join(line.split()))
                    residues += length
                    current += length
            longest = max(longest, current)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"input changed during profiling: {path}")
    if not sequences:
        raise ValueError(f"empty sequence input: {path}")
    return {"path": str(path.resolve()), "kind": kind, "bytes": after.st_size,
            "mtime_ns": after.st_mtime_ns, "sha256": checksum.hexdigest(),
            "sequences": sequences, "residues": residues, "max_length": longest}


def profile(workflow, fasta=(), fastq=(), *, species=None, settings=None,
            stages=(), task_count=1):
    if workflow not in WORKFLOWS or type(task_count) is not int or not 1 <= task_count <= 75000:
        raise ValueError("unsupported workflow or task count")
    paths = [(Path(p).resolve(strict=True), kind)
             for kind, values in (("fasta", fasta), ("fastq", fastq)) for p in values]
    if not paths or len({p for p, _ in paths}) != len(paths):
        raise ValueError("provide nonempty, distinct input files")
    files = [scan(p, kind) for p, kind in sorted(paths)]
    if species is None:
        species = (max(f["sequences"] for f in files) if workflow == "gg_gene_evolution"
                   else len(fasta))
    if type(species) is not int or species < 0:
        raise ValueError("species must be a nonnegative integer")
    features = {"species": species, "files": len(files),
                "sequences": sum(f["sequences"] for f in files),
                "residues": sum(f["residues"] for f in files),
                "max_sequences": max(f["sequences"] for f in files),
                "max_length": max(f["max_length"] for f in files),
                "bytes": sum(f["bytes"] for f in files)}
    # Order/path/mtime are not scientific identity; contents and settings are.
    identity = {"workflow": workflow, "inputs": sorted((f["kind"], f["sha256"]) for f in files),
                "settings": settings or {}, "stages": sorted(set(stages)),
                "species": species, "task_count": task_count}
    return {"schema_version": 1, "workflow": workflow, "input_sha256": digest(identity),
            "features": features, "task_count": task_count,
            "stages": identity["stages"], "files": files,
            "settings_sha256": digest(settings or {})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow", required=True, choices=sorted(WORKFLOWS))
    parser.add_argument("--fasta", action="append", default=[])
    parser.add_argument("--fastq", action="append", default=[])
    parser.add_argument("--species", type=int)
    parser.add_argument("--task-count", type=int, default=1,
                        help="array members; supplied FASTA/FASTQ describe ONE representative member")
    parser.add_argument("--stage", action="append", default=[])
    parser.add_argument("--settings", type=Path, help="JSON scientific settings for comparable history")
    args = parser.parse_args()
    value = profile(args.workflow, args.fasta, args.fastq, species=args.species,
                    settings=json.loads(args.settings.read_text()) if args.settings else None,
                    stages=args.stage, task_count=args.task_count)
    print(json.dumps(value, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
