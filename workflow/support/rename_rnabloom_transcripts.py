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


def load_clusters(path):
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
                continue
            transcript_id = fields[0].strip()
            cluster_id = fields[1].strip()
            if not transcript_id or not cluster_id:
                continue
            cluster_map[transcript_id] = cluster_id
    return cluster_map


def wrap_sequence(sequence, width=80):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def main():
    args = parse_args()
    input_path = Path(args.input_fasta)
    output_path = Path(args.output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cluster_map = load_clusters(args.clusters)
    cluster_indices = {}
    cluster_transcript_counts = {}
    next_cluster_index = 1

    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        for transcript_id, sequence in parse_fasta(input_path):
            cluster_key = cluster_map.get(transcript_id, f"singleton::{transcript_id}")
            cluster_index = cluster_indices.get(cluster_key)
            if cluster_index is None:
                cluster_index = next_cluster_index
                cluster_indices[cluster_key] = cluster_index
                next_cluster_index += 1
            transcript_index = cluster_transcript_counts.get(cluster_key, 0) + 1
            cluster_transcript_counts[cluster_key] = transcript_index

            renamed_id = f"{args.species_prefix}_c{cluster_index:06d}-i{transcript_index:06d}"
            handle.write(f">{renamed_id}\n")
            for chunk in wrap_sequence(sequence):
                handle.write(f"{chunk}\n")


if __name__ == "__main__":
    main()
