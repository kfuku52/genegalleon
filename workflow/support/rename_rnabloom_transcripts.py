#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename RNA-Bloom2 transcript FASTA headers into stable cluster-aware GeneGalleon IDs."
    )
    parser.add_argument("--input-fasta", required=True, help="Input RNA-Bloom2 transcript FASTA(.gz).")
    parser.add_argument(
        "--clusters",
        default="",
        help="Optional Corset clusters.txt file. When omitted, each transcript becomes its own singleton cluster.",
    )
    parser.add_argument("--species-prefix", required=True, help="Species prefix used in renamed FASTA headers.")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA path.")
    return parser.parse_args()


def open_maybe_gzip(path, mode):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def parse_fasta(path):
    current_id = None
    seq_chunks = []
    with open_maybe_gzip(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(seq_chunks)
                current_id = line[1:].split(None, 1)[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if current_id is not None:
        yield current_id, "".join(seq_chunks)


def load_clusters(path, *, strict=False):
    cluster_map = {}
    if not path:
        return cluster_map
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                if strict:
                    raise ValueError(f"Invalid Corset cluster assignment: {line}")
                continue
            transcript_id = fields[0].strip()
            cluster_id = fields[1].strip()
            if not transcript_id or not cluster_id:
                if strict:
                    raise ValueError(f"Invalid Corset cluster assignment: {line}")
                continue
            if strict and transcript_id in cluster_map and cluster_map[transcript_id] != cluster_id:
                raise ValueError(f"Conflicting Corset cluster assignments for {transcript_id}")
            cluster_map[transcript_id] = cluster_id
    return cluster_map


def wrap_sequence(sequence, width=80):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def cluster_transcripts(transcripts, cluster_map, species):
    cluster_indices = {}
    cluster_transcript_counts = {}
    for transcript_id, sequence in transcripts:
        cluster_key = cluster_map.get(transcript_id, f"singleton::{transcript_id}")
        cluster_index = cluster_indices.setdefault(cluster_key, len(cluster_indices) + 1)
        transcript_index = cluster_transcript_counts.get(cluster_key, 0) + 1
        cluster_transcript_counts[cluster_key] = transcript_index
        gene_id = f"{species}_c{cluster_index:06d}"
        yield transcript_id, gene_id, f"{gene_id}-i{transcript_index:06d}", sequence


def main():
    args = parse_args()
    input_path = Path(args.input_fasta)
    output_path = Path(args.output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster_map = load_clusters(args.clusters)
    renamed = cluster_transcripts(parse_fasta(input_path), cluster_map, args.species_prefix)

    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        for _, _, renamed_id, sequence in renamed:
            handle.write(f">{renamed_id}\n")
            for chunk in wrap_sequence(sequence):
                handle.write(f"{chunk}\n")


if __name__ == "__main__":
    main()
