#!/usr/bin/env python3
"""Validated full-CDS inputs and portable CSUBST 3Di bundle metadata."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from Bio import Phylo, SeqIO

SCHEMA = "genegalleon-csubst-input-v2"
STRUCTURAL_DIR = "csubst.3di"
FIT_SUFFIXES = ("fasta", "nwk", "treefile", "state", "rate", "iqtree", "log", "ckp.gz")


def read_alignment(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    ids = [record.id for record in records]
    if not records or len(ids) != len(set(ids)):
        raise ValueError(f"Empty alignment or duplicate sequence IDs: {path}")
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1 or not next(iter(lengths)) or next(iter(lengths)) % 3:
        raise ValueError(f"Full CDS alignment must have equal, nonzero lengths in codon frame: {path}")
    return records


def prepare_full_alignment(source, tree, destination):
    records = read_alignment(source)
    tips = [tip.name for tip in Phylo.read(tree, "newick").get_terminals()]
    if any(not tip for tip in tips) or len(tips) != len(set(tips)):
        raise ValueError(f"Rooted tree requires unique named tips: {tree}")
    by_id = {record.id: record for record in records}
    missing = set(tips) - by_id.keys()
    if missing:
        raise ValueError(f"Full CDS alignment lacks tree tips: {sorted(missing)}")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([by_id[tip] for tip in tips], destination, "fasta")


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def structural_directory(bundle, verify=True):
    bundle = Path(bundle)
    path = bundle / "csubst.input.json"
    metadata = json.loads(path.read_text()) if path.is_file() else {}
    if metadata.get("schema") != SCHEMA or metadata.get("structural_directory") != STRUCTURAL_DIR:
        raise ValueError(f"3di20 requires a full-CDS bundle. Regenerate the iqtree_anc stage from the untrimmed CDS alignment: {bundle}")
    directory = bundle / STRUCTURAL_DIR
    if verify:
        expected = {f"csubst.{suffix}" for suffix in FIT_SUFFIXES} | {
            "csubst_3di_state_cache.npz", "inspect/csubst_alignment_3di.fa"}
        hashes = metadata.get("structural_sha256", {})
        if set(hashes) != expected:
            raise ValueError(f"Incomplete 3Di bundle metadata: {path}")
        for name in sorted(expected):
            file = directory / name
            if not file.is_file() or file.stat().st_size == 0 or _digest(file) != hashes[name]:
                raise ValueError(f"Missing or changed 3Di bundle input: {file}")
    return directory


def finalize(bundle, genetic_code):
    bundle = Path(bundle)
    directory = bundle / STRUCTURAL_DIR
    names = [f"csubst.{suffix}" for suffix in FIT_SUFFIXES] + [
        "csubst_3di_state_cache.npz", "inspect/csubst_alignment_3di.fa"]
    metadata = {"schema": SCHEMA, "genetic_code": int(genetic_code),
                "structural_directory": STRUCTURAL_DIR, "site_coordinates": "full_cds_codon_1based",
                "structural_sha256": {name: _digest(directory / name) for name in names}}
    (bundle / "csubst.input.json").write_text(json.dumps(metadata, indent=2) + "\n")
    structural_directory(bundle)


def write_structural_tip_alignment(codon_alignment, output):
    directory = Path(codon_alignment).parent
    structural_directory(directory.parent)
    records = read_alignment(codon_alignment)
    with open(directory / "inspect/csubst_alignment_3di.fa") as handle:
        structural = list(SeqIO.parse(handle, "fasta"))
    by_id = {record.id: record for record in structural}
    if len(by_id) != len(structural):
        raise ValueError("Duplicate IDs in CSUBST 3Di alignment")
    selected = []
    for record in records:
        state = by_id.get(record.id)
        if state is None or len(state.seq) * 3 != len(record.seq):
            raise ValueError(f"CSUBST 3Di alignment does not match full CDS: {record.id}")
        selected.append(state)
    SeqIO.write(selected, output, "fasta")
    return str(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare")
    for name in ("source", "tree", "destination"):
        prepare.add_argument("--" + name, required=True)
    finish = sub.add_parser("finalize")
    finish.add_argument("--bundle", required=True)
    finish.add_argument("--genetic-code", type=int, required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("--bundle", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare_full_alignment(args.source, args.tree, args.destination)
    elif args.command == "finalize":
        finalize(args.bundle, args.genetic_code)
    else:
        structural_directory(args.bundle)


if __name__ == "__main__":
    main()
