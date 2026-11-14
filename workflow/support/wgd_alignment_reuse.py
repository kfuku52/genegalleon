"""Reuse exact codon alignments from an explicitly supplied WGD archive."""

import csv
import json
import os
import platform
import re
import shutil
from pathlib import Path

try:
    from pairwise_synteny import digest
    from pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS
except ImportError:
    from .pairwise_synteny import digest
    from .pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS

ALIGNMENT_NAME = re.compile(r"cross_ks_aligned_pairs\.code([1-9][0-9]*)\.tsv\Z")
CODON_ALIGNMENT = re.compile(r"(?:[ACGTRYSWKMBDHVN]{3}|---)+\Z")
MAFFT_DEFAULT = re.compile(rb"^[ \t]*prefix=(/[a-zA-Z0-9_./-]+)[ \t]*$", re.MULTILINE)


def alignment_identity(ds):
    """Fingerprint MAFFT and its bundled binaries, not merely the wrapper."""
    executable = shutil.which("mafft")
    if not executable:
        return None, ()
    mafft = Path(executable).resolve()
    override = os.environ.get("MAFFT_BINARIES")
    defaults = set(MAFFT_DEFAULT.findall(mafft.read_bytes())) if not override else ()
    if not override and len(defaults) != 1:
        return None, ()
    root = Path(override or os.fsdecode(next(iter(defaults)))).resolve()
    if not (root / "tbfast").is_file():
        return None, ()
    artifacts = list(sorted(p for p in root.iterdir() if p.is_file()))
    helpers = {p.name: digest(p) for p in artifacts}
    related = {Path(path).name: value for path, value in ds["source_hashes"].items()
               if Path(path).name in {"backalign.py", "codonutil.py", "translate.py", "util.py", "pairwise_synteny_ds.py"}}
    identity = {"version": 1, "args": list(MAFFT_ARGS), "mafft": ds["mafft_version"],
                "wrapper": ds["mafft_executable_sha256"], "helpers": helpers, "backtranslation": related,
                "platform": [platform.system(), platform.machine(), *platform.libc_ver()]}
    return identity, [(str(mafft), identity["wrapper"]), *((str(p.resolve()), helpers[p.name]) for p in artifacts)]


def reuse_sources(summary, identity):
    archived_hash = digest(summary)
    previous = json.loads(summary.read_text(encoding="utf-8"))
    if identity is None or previous.get("plan", {}).get("alignment_identity") != identity:
        raise ValueError("Prior WGD alignments do not have the same verified MAFFT and backtranslation identity")
    sources, hashes = {}, {str(summary): archived_hash}
    for name, expected in previous.get("output_hashes", {}).items():
        match = ALIGNMENT_NAME.fullmatch(name)
        if match is None:
            continue
        source = (summary.parent / name).resolve()
        if digest(source) != expected:
            raise ValueError(f"Prior WGD alignment hash does not match its summary: {name}")
        sources[match[1]] = str(source)
        hashes[str(source)] = expected
    if not sources:
        raise ValueError("Prior WGD summary has no cross-block codon alignments")
    return sources, hashes


def matching_rows(path, selected, sequences):
    requests = {pair: info for pair, info in selected}
    rows, seen = {}, set()
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if tuple(reader.fieldnames or ()) != PAIR_COLUMNS:
            raise ValueError(f"Invalid prior WGD alignment header: {path}")
        for row in reader:
            if None in row or any(val is None for val in row.values()):
                raise ValueError(f"Invalid prior WGD alignment table: {path}")
            pair = row["target_gene"], row["query_gene"]
            if row["pair_id"] != "|".join(pair) or pair in seen:
                raise ValueError(f"Duplicate or invalid pair in prior WGD alignments: {row['pair_id']}")
            seen.add(pair)
            if pair not in requests:
                continue
            info = requests[pair]
            first, second = row["sequence_1"], row["sequence_2"]
            if (not first or len(first) != len(second) or len(first) % 3
                    or any(not CODON_ALIGNMENT.fullmatch(seq) for seq in (first, second))):
                raise ValueError(f"Malformed prior WGD codon alignment: {row['pair_id']}")
            if all(seq.replace("-", "") == sequences[species][gene][0]
                   for seq, species, gene in zip((first, second), info[:2], pair, strict=True)):
                rows[pair] = row
    return rows
