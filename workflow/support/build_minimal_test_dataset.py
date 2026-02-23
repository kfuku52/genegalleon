#!/usr/bin/env python3

"""Build a minimal workspace test dataset.

This script extracts a small, analysis-ready subset from an existing workspace:
- species: configurable (default: Arabidopsis, Cephalotus, Dionaea, Nepenthes, Oryza)
- genes: all AHA/YABBY homologs (DIAMOND-based) + first N complete BUSCO genes
- files: species_cds, species_expression, species_genome, species_gff

Genome extraction keeps only windows around target genes (gene +/- flank bp),
and GFF coordinates are shifted to match extracted window coordinates.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


DEFAULT_TARGET_SPECIES = [
    "Arabidopsis_thaliana",
    "Cephalotus_follicularis",
    "Dionaea_muscipula",
    "Nepenthes_gracilis",
    "Oryza_sativa",
]

FAMILY_QUERIES = ["AHA", "YABBY"]

DEFAULT_THRESHOLDS = {
    "evalue": 1e-40,
    "qcovhsp": 60.0,
    "pident": 35.0,
}

# Keep a tiny synthetic neighborhood signal for YABBY synteny panel testing.
SYNTENY_DUMMY_CONFIG = {
    "Cephalotus_follicularis": {
        "anchors": [
            {
                "anchor": "Cephalotus_follicularis_Cfol-v3-25951",
                "upstream": [
                    "Cephalotus_follicularis_Cfol-v3-90001",
                    "Cephalotus_follicularis_Cfol-v3-90002",
                    "Cephalotus_follicularis_Cfol-v3-90003",
                    "Cephalotus_follicularis_Cfol-v3-90004",
                    "Cephalotus_follicularis_Cfol-v3-90005",
                ],
                "downstream": [
                    "Cephalotus_follicularis_Cfol-v3-90006",
                    "Cephalotus_follicularis_Cfol-v3-90007",
                    "Cephalotus_follicularis_Cfol-v3-90008",
                    "Cephalotus_follicularis_Cfol-v3-90009",
                    "Cephalotus_follicularis_Cfol-v3-90010",
                ],
            },
            {
                "anchor": "Cephalotus_follicularis_Cfol-v3-04239",
                "upstream": [
                    "Cephalotus_follicularis_Cfol-v3-90011",
                    "Cephalotus_follicularis_Cfol-v3-90012",
                    "Cephalotus_follicularis_Cfol-v3-90013",
                    "Cephalotus_follicularis_Cfol-v3-90014",
                    "Cephalotus_follicularis_Cfol-v3-90015",
                ],
                "downstream": [
                    "Cephalotus_follicularis_Cfol-v3-90016",
                    "Cephalotus_follicularis_Cfol-v3-90017",
                    "Cephalotus_follicularis_Cfol-v3-90018",
                    "Cephalotus_follicularis_Cfol-v3-90019",
                    "Cephalotus_follicularis_Cfol-v3-90020",
                ],
            },
            {
                "anchor": "Cephalotus_follicularis_Cfol-v3-07401",
                "upstream": [
                    "Cephalotus_follicularis_Cfol-v3-90021",
                    "Cephalotus_follicularis_Cfol-v3-90022",
                    "Cephalotus_follicularis_Cfol-v3-90023",
                    "Cephalotus_follicularis_Cfol-v3-90024",
                    "Cephalotus_follicularis_Cfol-v3-90025",
                ],
                "downstream": [
                    "Cephalotus_follicularis_Cfol-v3-90026",
                    "Cephalotus_follicularis_Cfol-v3-90027",
                    "Cephalotus_follicularis_Cfol-v3-90028",
                    "Cephalotus_follicularis_Cfol-v3-90029",
                    "Cephalotus_follicularis_Cfol-v3-90030",
                ],
            },
            {
                "anchor": "Cephalotus_follicularis_Cfol-v3-04540",
                "upstream": [
                    "Cephalotus_follicularis_Cfol-v3-90031",
                    "Cephalotus_follicularis_Cfol-v3-90032",
                    "Cephalotus_follicularis_Cfol-v3-90033",
                    "Cephalotus_follicularis_Cfol-v3-90034",
                    "Cephalotus_follicularis_Cfol-v3-90035",
                ],
                "downstream": [
                    "Cephalotus_follicularis_Cfol-v3-90036",
                    "Cephalotus_follicularis_Cfol-v3-90037",
                    "Cephalotus_follicularis_Cfol-v3-90038",
                    "Cephalotus_follicularis_Cfol-v3-90039",
                    "Cephalotus_follicularis_Cfol-v3-90040",
                ],
            },
        ],
        # Shared groups are used to enforce reproducible cross-tip links.
        # Other dummy CDS are generated deterministically from dummy IDs.
        "shared_groups": {
            # Diagonal upstream links (offsets: -1, -2, -3, -4 across anchors).
            "diag_up": [
                "Cephalotus_follicularis_Cfol-v3-90001",
                "Cephalotus_follicularis_Cfol-v3-90012",
                "Cephalotus_follicularis_Cfol-v3-90023",
                "Cephalotus_follicularis_Cfol-v3-90034",
            ],
            # Diagonal downstream links (offsets: +4, +3, +2, +1 across anchors).
            "diag_down": [
                "Cephalotus_follicularis_Cfol-v3-90009",
                "Cephalotus_follicularis_Cfol-v3-90018",
                "Cephalotus_follicularis_Cfol-v3-90027",
                "Cephalotus_follicularis_Cfol-v3-90036",
            ],
        },
        "shared_group_sequences": {
            "diag_up": (
                "atgaaaactgctcttcaaaatgttcatggttatcctcttcgtgattctgttcaagctttc"
            ),
            "diag_down": (
                "atggacttccaagttgctaatcctggatttcatcgtaaggtttctgctcaagataactgg"
            ),
        },
        "dummy_cds_length_nt": 60,
        "upstream_base_bp": 120,
        "upstream_step_bp": 120,
        "downstream_base_bp": 120,
        "downstream_step_bp": 120,
    },
    "Nepenthes_gracilis": {
        "anchors": [
            {
                "anchor": "Nepenthes_gracilis_Nepgr000233",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91001",
                    "Nepenthes_gracilis_Nepgr91002",
                    "Nepenthes_gracilis_Nepgr91003",
                    "Nepenthes_gracilis_Nepgr91004",
                    "Nepenthes_gracilis_Nepgr91005",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91051",
                    "Nepenthes_gracilis_Nepgr91052",
                    "Nepenthes_gracilis_Nepgr91053",
                    "Nepenthes_gracilis_Nepgr91054",
                    "Nepenthes_gracilis_Nepgr91055",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr002443",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91006",
                    "Nepenthes_gracilis_Nepgr91007",
                    "Nepenthes_gracilis_Nepgr91008",
                    "Nepenthes_gracilis_Nepgr91009",
                    "Nepenthes_gracilis_Nepgr91010",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91056",
                    "Nepenthes_gracilis_Nepgr91057",
                    "Nepenthes_gracilis_Nepgr91058",
                    "Nepenthes_gracilis_Nepgr91059",
                    "Nepenthes_gracilis_Nepgr91060",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr006372",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91011",
                    "Nepenthes_gracilis_Nepgr91012",
                    "Nepenthes_gracilis_Nepgr91013",
                    "Nepenthes_gracilis_Nepgr91014",
                    "Nepenthes_gracilis_Nepgr91015",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91061",
                    "Nepenthes_gracilis_Nepgr91062",
                    "Nepenthes_gracilis_Nepgr91063",
                    "Nepenthes_gracilis_Nepgr91064",
                    "Nepenthes_gracilis_Nepgr91065",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr009191",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91016",
                    "Nepenthes_gracilis_Nepgr91017",
                    "Nepenthes_gracilis_Nepgr91018",
                    "Nepenthes_gracilis_Nepgr91019",
                    "Nepenthes_gracilis_Nepgr91020",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91066",
                    "Nepenthes_gracilis_Nepgr91067",
                    "Nepenthes_gracilis_Nepgr91068",
                    "Nepenthes_gracilis_Nepgr91069",
                    "Nepenthes_gracilis_Nepgr91070",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr011352",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91021",
                    "Nepenthes_gracilis_Nepgr91022",
                    "Nepenthes_gracilis_Nepgr91023",
                    "Nepenthes_gracilis_Nepgr91024",
                    "Nepenthes_gracilis_Nepgr91025",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91071",
                    "Nepenthes_gracilis_Nepgr91072",
                    "Nepenthes_gracilis_Nepgr91073",
                    "Nepenthes_gracilis_Nepgr91074",
                    "Nepenthes_gracilis_Nepgr91075",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr014538",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91026",
                    "Nepenthes_gracilis_Nepgr91027",
                    "Nepenthes_gracilis_Nepgr91028",
                    "Nepenthes_gracilis_Nepgr91029",
                    "Nepenthes_gracilis_Nepgr91030",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91076",
                    "Nepenthes_gracilis_Nepgr91077",
                    "Nepenthes_gracilis_Nepgr91078",
                    "Nepenthes_gracilis_Nepgr91079",
                    "Nepenthes_gracilis_Nepgr91080",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr022393",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91031",
                    "Nepenthes_gracilis_Nepgr91032",
                    "Nepenthes_gracilis_Nepgr91033",
                    "Nepenthes_gracilis_Nepgr91034",
                    "Nepenthes_gracilis_Nepgr91035",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91081",
                    "Nepenthes_gracilis_Nepgr91082",
                    "Nepenthes_gracilis_Nepgr91083",
                    "Nepenthes_gracilis_Nepgr91084",
                    "Nepenthes_gracilis_Nepgr91085",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr024231",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91036",
                    "Nepenthes_gracilis_Nepgr91037",
                    "Nepenthes_gracilis_Nepgr91038",
                    "Nepenthes_gracilis_Nepgr91039",
                    "Nepenthes_gracilis_Nepgr91040",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91086",
                    "Nepenthes_gracilis_Nepgr91087",
                    "Nepenthes_gracilis_Nepgr91088",
                    "Nepenthes_gracilis_Nepgr91089",
                    "Nepenthes_gracilis_Nepgr91090",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr028028",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91041",
                    "Nepenthes_gracilis_Nepgr91042",
                    "Nepenthes_gracilis_Nepgr91043",
                    "Nepenthes_gracilis_Nepgr91044",
                    "Nepenthes_gracilis_Nepgr91045",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91091",
                    "Nepenthes_gracilis_Nepgr91092",
                    "Nepenthes_gracilis_Nepgr91093",
                    "Nepenthes_gracilis_Nepgr91094",
                    "Nepenthes_gracilis_Nepgr91095",
                ],
            },
            {
                "anchor": "Nepenthes_gracilis_Nepgr029156",
                "upstream": [
                    "Nepenthes_gracilis_Nepgr91046",
                    "Nepenthes_gracilis_Nepgr91047",
                    "Nepenthes_gracilis_Nepgr91048",
                    "Nepenthes_gracilis_Nepgr91049",
                    "Nepenthes_gracilis_Nepgr91050",
                ],
                "downstream": [
                    "Nepenthes_gracilis_Nepgr91096",
                    "Nepenthes_gracilis_Nepgr91097",
                    "Nepenthes_gracilis_Nepgr91098",
                    "Nepenthes_gracilis_Nepgr91099",
                    "Nepenthes_gracilis_Nepgr91100",
                ],
            },
        ],
        "shared_groups": {
            "diag_up": [
                "Nepenthes_gracilis_Nepgr91001",
                "Nepenthes_gracilis_Nepgr91007",
                "Nepenthes_gracilis_Nepgr91013",
                "Nepenthes_gracilis_Nepgr91019",
                "Nepenthes_gracilis_Nepgr91025",
                "Nepenthes_gracilis_Nepgr91026",
                "Nepenthes_gracilis_Nepgr91032",
                "Nepenthes_gracilis_Nepgr91038",
                "Nepenthes_gracilis_Nepgr91044",
                "Nepenthes_gracilis_Nepgr91050",
            ],
            "diag_down": [
                "Nepenthes_gracilis_Nepgr91055",
                "Nepenthes_gracilis_Nepgr91059",
                "Nepenthes_gracilis_Nepgr91063",
                "Nepenthes_gracilis_Nepgr91067",
                "Nepenthes_gracilis_Nepgr91071",
                "Nepenthes_gracilis_Nepgr91080",
                "Nepenthes_gracilis_Nepgr91084",
                "Nepenthes_gracilis_Nepgr91088",
                "Nepenthes_gracilis_Nepgr91092",
                "Nepenthes_gracilis_Nepgr91096",
            ],
        },
        "shared_group_sequences": {
            "diag_up": (
                "atggctgaaatctctgctgttcaagattcctggactgactgcttatgttgaagccgatgc"
            ),
            "diag_down": (
                "atgcagttcaagccattggttgctgaagttgctgaagcctacgttcaagctgttggtcctgaa"
            ),
        },
        "dummy_cds_length_nt": 60,
        "upstream_base_bp": 120,
        "upstream_step_bp": 120,
        "downstream_base_bp": 120,
        "downstream_step_bp": 120,
    },
}


@dataclass
class Window:
    seqid: str
    start: int
    end: int
    cores: Set[str]


@dataclass
class EffectiveWindow:
    seqid: str
    start: int
    end: int
    window_id: str
    cores: Set[str]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def find_first_file(directory: Path, prefix: str) -> Path:
    matches = sorted(directory.glob(f"{prefix}*"))
    if not matches:
        raise FileNotFoundError(f"No file matching prefix '{prefix}' in {directory}")
    return matches[0]


def fasta_iter(path: Path) -> Iterator[Tuple[str, str, str]]:
    with open_text(path) as handle:
        header = None
        seq_buf: List[str] = []
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_buf)
                    yield header, header.split()[0], seq
                header = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line.strip())
        if header is not None:
            seq = "".join(seq_buf)
            yield header, header.split()[0], seq


def write_fasta_record(out_handle, header: str, seq: str, width: int = 60) -> None:
    out_handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        out_handle.write(seq[i : i + width] + "\n")


def run_cmd(cmd: Sequence[str], quiet: bool = False) -> None:
    kwargs = {}
    if quiet:
        kwargs = {"stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL}
    subprocess.run(cmd, check=True, **kwargs)


def load_query_ids(source_pg: Path) -> Dict[str, List[str]]:
    qdir = source_pg / "query_gene"
    out: Dict[str, List[str]] = {}
    for fam in FAMILY_QUERIES:
        path = qdir / fam
        ids = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
        if not ids:
            raise RuntimeError(f"No query IDs found in {path}")
        out[fam] = ids
    return out


def load_busco_first_n(
    source_pg: Path, species: Sequence[str], n: int
) -> Dict[str, List[str]]:
    busco_dir = source_pg / "species_cds_busco_full"
    out: Dict[str, List[str]] = {}
    for sp in species:
        path = busco_dir / f"{sp}.busco.full.tsv"
        picked: List[str] = []
        with open_text(path) as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                status = parts[1]
                seq = parts[2]
                if status != "Complete":
                    continue
                gid = seq.split(":", 1)[0]
                picked.append(gid)
                if len(picked) >= n:
                    break
        if len(picked) < n:
            raise RuntimeError(f"{sp}: only {len(picked)} BUSCO IDs found in {path}")
        out[sp] = picked
    return out


def extract_query_cds(
    arabidopsis_cds: Path, query_ids: Dict[str, List[str]], tmpdir: Path
) -> Dict[str, Path]:
    all_records: Dict[str, Tuple[str, str]] = {}
    for header, seqid, seq in fasta_iter(arabidopsis_cds):
        all_records[seqid] = (header, seq)

    fam_to_cds: Dict[str, Path] = {}
    for fam, ids in query_ids.items():
        out_cds = tmpdir / f"{fam}.query.cds.fa"
        with open(out_cds, "w") as out:
            for gid in ids:
                rec = all_records.get(gid)
                if rec is None:
                    # Source query lists can reference a different isoform than
                    # the current primary-transcript CDS set. Fallback to the
                    # longest available isoform for the same gene.
                    core = gid
                    if core.startswith("Arabidopsis_thaliana_"):
                        core = core[len("Arabidopsis_thaliana_") :]
                    gene_base = core.split(".", 1)[0]
                    pref = f"Arabidopsis_thaliana_{gene_base}"
                    candidates = [
                        sid
                        for sid in all_records
                        if sid == pref or sid.startswith(pref + ".")
                    ]
                    if not candidates:
                        raise RuntimeError(
                            f"Query ID '{gid}' not found in {arabidopsis_cds}"
                        )
                    best_sid = max(candidates, key=lambda sid: len(all_records[sid][1]))
                    rec = all_records[best_sid]
                header, seq = rec
                write_fasta_record(out, header, seq)
        fam_to_cds[fam] = out_cds
    return fam_to_cds


def prepare_species_diamond_db(cds_file: Path, tmpdir: Path, sp: str) -> Path:
    pep = tmpdir / f"{sp}.pep.fa"
    db = tmpdir / f"{sp}.db"
    run_cmd(["seqkit", "translate", str(cds_file), "-o", str(pep)], quiet=True)
    run_cmd(["diamond", "makedb", "--in", str(pep), "-d", str(db)], quiet=True)
    return db


def run_family_diamond_hits(
    species_to_cds: Dict[str, Path],
    fam_to_query_cds: Dict[str, Path],
    thresholds: Dict[str, float],
    tmpdir: Path,
) -> Dict[str, Dict[str, Set[str]]]:
    # translate query CDS once per family
    fam_to_query_pep: Dict[str, Path] = {}
    for fam, query_cds in fam_to_query_cds.items():
        query_pep = tmpdir / f"{fam}.query.pep.fa"
        run_cmd(["seqkit", "translate", str(query_cds), "-o", str(query_pep)], quiet=True)
        fam_to_query_pep[fam] = query_pep

    # prepare per-species db once
    species_to_db: Dict[str, Path] = {}
    for sp, cds_file in species_to_cds.items():
        species_to_db[sp] = prepare_species_diamond_db(cds_file, tmpdir, sp)

    out: Dict[str, Dict[str, Set[str]]] = {
        sp: {fam: set() for fam in FAMILY_QUERIES} for sp in species_to_cds
    }
    for sp, db in species_to_db.items():
        for fam, query_pep in fam_to_query_pep.items():
            blast_out = tmpdir / f"{sp}.{fam}.diamond.tsv"
            run_cmd(
                [
                    "diamond",
                    "blastp",
                    "--query",
                    str(query_pep),
                    "--db",
                    str(db),
                    "--out",
                    str(blast_out),
                    "--outfmt",
                    "6",
                    "qseqid",
                    "sseqid",
                    "pident",
                    "length",
                    "evalue",
                    "bitscore",
                    "qcovhsp",
                    "scovhsp",
                    "--max-target-seqs",
                    "2000",
                    "--threads",
                    "4",
                    "--sensitive",
                ],
                quiet=True,
            )
            with open(blast_out) as handle:
                for raw in handle:
                    qseqid, sseqid, pident, _length, evalue, _bits, qcov, _scov = raw.rstrip(
                        "\n"
                    ).split("\t")
                    # qseqid parsed for completeness (and to match outfmt shape)
                    _ = qseqid
                    if float(evalue) > thresholds["evalue"]:
                        continue
                    if float(qcov) < thresholds["qcovhsp"]:
                        continue
                    if float(pident) < thresholds["pident"]:
                        continue
                    out[sp][fam].add(sseqid)
    return out


def strip_species_prefix(sp: str, gene_id: str) -> str:
    pref = sp + "_"
    if gene_id.startswith(pref):
        return gene_id[len(pref) :]
    return gene_id


def token_variants(core: str) -> Set[str]:
    out = {core}
    out.add(core.replace("-", "_"))
    out.add(core.replace("_", "-"))
    return {x for x in out if x}


def build_token_pattern(core_ids: Iterable[str]) -> Tuple[re.Pattern, Dict[str, str]]:
    token_to_core: Dict[str, str] = {}
    for core in core_ids:
        for tok in token_variants(core):
            token_to_core[tok] = core
    # Longest tokens first to reduce accidental short-token matches.
    tokens = sorted(token_to_core.keys(), key=len, reverse=True)
    if not tokens:
        raise RuntimeError("No tokens to build regex pattern")
    pattern = re.compile("|".join(re.escape(t) for t in tokens))
    return pattern, token_to_core


def collect_synteny_dummy_ids(conf: Dict) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()
    for anchor_spec in conf.get("anchors", []):
        for direction in ("upstream", "downstream"):
            for dummy_id in anchor_spec.get(direction, []):
                if dummy_id in seen:
                    continue
                seen.add(dummy_id)
                out.append(dummy_id)
    return out


def make_deterministic_dummy_cds(dummy_id: str, length_nt: int = 60) -> str:
    length_nt = max(18, int(length_nt))
    if length_nt % 3 != 0:
        length_nt -= length_nt % 3
    codons = [
        "gct",
        "gcc",
        "gca",
        "gcg",
        "tct",
        "tcc",
        "tca",
        "tcg",
        "cct",
        "ccc",
        "cca",
        "ccg",
        "act",
        "acc",
        "aca",
        "acg",
        "gtt",
        "gtc",
        "gta",
        "gtg",
        "gat",
        "gac",
        "aat",
        "aac",
        "caa",
        "cag",
        "tat",
        "tac",
        "cat",
        "cac",
        "aaa",
        "aag",
        "cgt",
        "cgc",
        "cga",
        "cgg",
        "agt",
        "agc",
        "att",
        "atc",
        "ata",
        "ctt",
        "ctc",
        "cta",
        "ctg",
        "ttt",
        "ttc",
        "tgg",
        "tta",
        "ttg",
        "ggt",
        "ggc",
        "gga",
        "ggg",
    ]
    digest = hashlib.sha256(dummy_id.encode("utf-8")).digest()
    n_codon = length_nt // 3
    seq_codons = ["atg"]
    for i in range(1, n_codon):
        byte_val = digest[(i - 1) % len(digest)]
        idx = (byte_val + (i * 17)) % len(codons)
        seq_codons.append(codons[idx])
    return "".join(seq_codons)


def build_synteny_dummy_cds_map(conf: Dict) -> Dict[str, str]:
    dummy_ids = collect_synteny_dummy_ids(conf)
    length_nt = int(conf.get("dummy_cds_length_nt", 60))
    cds_map: Dict[str, str] = {}

    shared_groups = conf.get("shared_groups", {})
    shared_group_sequences = conf.get("shared_group_sequences", {})
    for group_name, members in shared_groups.items():
        seq = shared_group_sequences.get(group_name, "")
        seq = seq.strip().lower()
        if (not seq) or (len(seq) % 3 != 0):
            raise RuntimeError(
                f"Invalid shared_group_sequences entry for {group_name}: length must be non-zero and multiple of 3."
            )
        for dummy_id in members:
            cds_map[dummy_id] = seq

    used_sequences = set(cds_map.values())
    for dummy_id in dummy_ids:
        if dummy_id in cds_map:
            continue
        seq = make_deterministic_dummy_cds(dummy_id=dummy_id, length_nt=length_nt)
        salt = 0
        while seq in used_sequences:
            salt += 1
            seq = make_deterministic_dummy_cds(
                dummy_id=f"{dummy_id}#{salt}",
                length_nt=length_nt,
            )
        cds_map[dummy_id] = seq
        used_sequences.add(seq)

    missing_members = sorted(
        set(x for members in shared_groups.values() for x in members) - set(dummy_ids)
    )
    if missing_members:
        raise RuntimeError(
            "shared_groups include unknown dummy IDs: " + ",".join(missing_members)
        )
    return cds_map


def collect_gene_intervals(
    gff_file: Path,
    selected_gene_ids: Set[str],
    sp: str,
) -> Dict[str, List[Tuple[str, int, int]]]:
    core_ids = {strip_species_prefix(sp, gid) for gid in selected_gene_ids}
    pattern, token_to_core = build_token_pattern(core_ids)
    intervals: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)

    # First pass on gene/transcript-like rows.
    with open_text(gff_file) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in {"gene", "mRNA", "transcript"}:
                continue
            attr = parts[8]
            m = pattern.search(attr)
            if m is None:
                continue
            core = token_to_core[m.group(0)]
            start = int(parts[3])
            end = int(parts[4])
            if start > end:
                start, end = end, start
            intervals[core].append((parts[0], start, end))

    missing = sorted(core_ids - set(intervals))
    if missing:
        # Fallback pass on all feature rows for missing IDs.
        missing_pattern, missing_token_to_core = build_token_pattern(missing)
        with open_text(gff_file) as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attr = parts[8]
                m = missing_pattern.search(attr)
                if m is None:
                    continue
                core = missing_token_to_core[m.group(0)]
                start = int(parts[3])
                end = int(parts[4])
                if start > end:
                    start, end = end, start
                intervals[core].append((parts[0], start, end))

    return intervals


def build_windows(
    intervals: Dict[str, List[Tuple[str, int, int]]],
    flank: int,
) -> Dict[str, List[Window]]:
    raw_by_seqid: Dict[str, List[Window]] = defaultdict(list)
    for core, ivals in intervals.items():
        if not ivals:
            continue
        seq_count: Dict[str, int] = defaultdict(int)
        for seqid, _s, _e in ivals:
            seq_count[seqid] += 1
        best_seqid = max(seq_count.items(), key=lambda kv: kv[1])[0]
        starts = [s for seqid, s, _e in ivals if seqid == best_seqid]
        ends = [e for seqid, _s, e in ivals if seqid == best_seqid]
        gstart = min(starts)
        gend = max(ends)
        raw_by_seqid[best_seqid].append(
            Window(
                seqid=best_seqid,
                start=max(1, gstart - flank),
                end=gend + flank,
                cores={core},
            )
        )

    merged: Dict[str, List[Window]] = defaultdict(list)
    for seqid, wins in raw_by_seqid.items():
        wins_sorted = sorted(wins, key=lambda w: (w.start, w.end))
        cur = wins_sorted[0]
        for w in wins_sorted[1:]:
            if w.start <= cur.end + 1:
                cur.end = max(cur.end, w.end)
                cur.cores |= w.cores
            else:
                merged[seqid].append(cur)
                cur = Window(seqid=seqid, start=w.start, end=w.end, cores=set(w.cores))
        merged[seqid].append(cur)
    return merged


def extract_genome_windows(
    genome_file: Path,
    windows: Dict[str, List[Window]],
    out_fasta: Path,
) -> Dict[str, List[EffectiveWindow]]:
    effective: Dict[str, List[EffectiveWindow]] = defaultdict(list)
    with open(out_fasta, "w") as out:
        for _header, seqid, seq in fasta_iter(genome_file):
            if seqid not in windows:
                continue
            seqlen = len(seq)
            for w in windows[seqid]:
                start = max(1, w.start)
                end = min(seqlen, w.end)
                if start > end:
                    continue
                win_id = f"{seqid}:{start}-{end}"
                sub = seq[start - 1 : end]
                write_fasta_record(out, win_id, sub)
                effective[seqid].append(
                    EffectiveWindow(
                        seqid=seqid,
                        start=start,
                        end=end,
                        window_id=win_id,
                        cores=set(w.cores),
                    )
                )
    return effective


def write_shifted_gff(
    gff_file: Path,
    effective_windows: Dict[str, List[EffectiveWindow]],
    out_gff: Path,
) -> None:
    # Windows are non-overlapping per seqid. Keep sorted for deterministic mapping.
    sorted_wins: Dict[str, List[EffectiveWindow]] = {
        seqid: sorted(wins, key=lambda w: (w.start, w.end))
        for seqid, wins in effective_windows.items()
    }

    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        with open_text(gff_file) as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seqid = parts[0]
                if seqid not in sorted_wins:
                    continue
                start = int(parts[3])
                end = int(parts[4])
                if start > end:
                    start, end = end, start

                assigned: Optional[EffectiveWindow] = None
                for w in sorted_wins[seqid]:
                    if end < w.start:
                        break
                    if start > w.end:
                        continue
                    assigned = w
                    break
                if assigned is None:
                    continue

                ov_start = max(start, assigned.start)
                ov_end = min(end, assigned.end)
                if ov_start > ov_end:
                    continue

                parts[0] = assigned.window_id
                parts[3] = str(ov_start - assigned.start + 1)
                parts[4] = str(ov_end - assigned.start + 1)
                out.write("\t".join(parts) + "\n")


def write_trimmed_cds(
    cds_file: Path, selected_ids: Set[str], out_cds: Path
) -> Tuple[Set[str], int]:
    found: Set[str] = set()
    total = 0
    with open(out_cds, "w") as out:
        for header, seqid, seq in fasta_iter(cds_file):
            if seqid in selected_ids:
                write_fasta_record(out, header, seq)
                found.add(seqid)
                total += 1
    return found, total


def write_trimmed_expression(
    expression_file: Path, selected_ids: Set[str], out_expr: Path
) -> int:
    kept = 0
    with open_text(expression_file) as handle, open(out_expr, "w") as out:
        header = handle.readline()
        if header:
            out.write(header)
        for raw in handle:
            if not raw.strip():
                continue
            gid = raw.split("\t", 1)[0]
            if gid in selected_ids:
                out.write(raw)
                kept += 1
    return kept


def collect_core_gene_loci(gff_file: Path, core_ids: Set[str]) -> Dict[str, Tuple[str, int, int]]:
    loci: Dict[str, Tuple[str, int, int]] = {}
    if not core_ids:
        return loci
    pattern, token_to_core = build_token_pattern(core_ids)
    with open_text(gff_file) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            m = pattern.search(attr)
            if m is None:
                continue
            core = token_to_core[m.group(0)]
            start = int(parts[3])
            end = int(parts[4])
            if start > end:
                start, end = end, start
            seqid = parts[0]
            prev = loci.get(core)
            if prev is None:
                loci[core] = (seqid, start, end)
                continue
            prev_seqid, prev_start, prev_end = prev
            if prev_seqid != seqid:
                # Keep the longer locus span if multiple seqids are found.
                prev_span = prev_end - prev_start
                new_span = end - start
                if new_span > prev_span:
                    loci[core] = (seqid, start, end)
                continue
            loci[core] = (seqid, min(prev_start, start), max(prev_end, end))
    return loci


def append_species_synteny_dummies(
    species: str,
    out_cds: Path,
    out_gff: Path,
    found_ids: Set[str],
) -> Set[str]:
    conf = SYNTENY_DUMMY_CONFIG.get(species)
    if conf is None:
        return set()
    if (not out_cds.exists()) or (not out_gff.exists()):
        return set()
    anchor_gene_ids = [x.get("anchor", "") for x in conf.get("anchors", [])]
    anchor_gene_ids = [x for x in anchor_gene_ids if x]
    if not anchor_gene_ids:
        return set()
    core_ids = {
        strip_species_prefix(species, gid).replace("-", "_")
        for gid in anchor_gene_ids
    }
    dummy_cds = build_synteny_dummy_cds_map(conf)
    upstream_base_bp = int(conf.get("upstream_base_bp", 120))
    upstream_step_bp = int(conf.get("upstream_step_bp", 120))
    downstream_base_bp = int(conf.get("downstream_base_bp", 120))
    downstream_step_bp = int(conf.get("downstream_step_bp", 120))
    loci = collect_core_gene_loci(out_gff, core_ids)
    missing_anchors = []
    dummy_ids: Set[str] = set()
    with open(out_cds, "a") as cds_out, open(out_gff, "a") as gff_out:
        for anchor_spec in conf.get("anchors", []):
            anchor = anchor_spec.get("anchor", "")
            if not anchor:
                continue
            if anchor not in found_ids:
                missing_anchors.append(anchor)
                continue
            core = strip_species_prefix(species, anchor).replace("-", "_")
            locus = loci.get(core)
            if locus is None:
                missing_anchors.append(anchor)
                continue
            seqid, gstart, gend = locus
            for direction in ("upstream", "downstream"):
                dummy_list = list(anchor_spec.get(direction, []))
                for rank, dummy_id in enumerate(dummy_list, start=1):
                    seq = dummy_cds.get(dummy_id, "")
                    if not seq:
                        raise RuntimeError(
                            f"Dummy CDS sequence not defined: {species} {anchor} {dummy_id}"
                        )
                    if direction == "upstream":
                        shift = upstream_base_bp + ((rank - 1) * upstream_step_bp)
                        dstart = max(1, gstart - shift)
                    else:
                        shift = downstream_base_bp + ((rank - 1) * downstream_step_bp)
                        dstart = gend + shift
                    dend = dstart + len(seq) - 1
                    dcore = strip_species_prefix(species, dummy_id).replace("-", "_")
                    write_fasta_record(cds_out, dummy_id, seq)
                    gff_out.write(
                        "\t".join(
                            [
                                seqid,
                                "CoGe",
                                "CDS",
                                str(dstart),
                                str(dend),
                                ".",
                                "+",
                                ".",
                                (
                                    f"Parent={dcore}.mRNA1;ID={dcore}.CDS1;"
                                    f"Name={dcore}.CDS1;Alias={dcore};CDS={dcore}.CDS1;"
                                    f"coge_fid=gg_dummy_{dcore.split('_')[-1]}"
                                ),
                            ]
                        )
                        + "\n"
                    )
                    dummy_ids.add(dummy_id)
    if missing_anchors:
        print(
            f"Warning: skipped synteny dummies for missing anchors in {species}: "
            + ",".join(sorted(missing_anchors))
        )
    return dummy_ids


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pg", required=True, help="Source workspace directory")
    parser.add_argument("--out-pg", required=True, help="Output minimal workspace directory")
    parser.add_argument(
        "--species",
        default=",".join(DEFAULT_TARGET_SPECIES),
        help=(
            "Comma-separated target species list "
            f"(default: {','.join(DEFAULT_TARGET_SPECIES)})"
        ),
    )
    parser.add_argument(
        "--query-reference-species",
        default="Arabidopsis_thaliana",
        help=(
            "Species prefix used to extract AHA/YABBY query CDS from source "
            "(default: Arabidopsis_thaliana)"
        ),
    )
    parser.add_argument(
        "--busco-n",
        type=int,
        default=5,
        help="Number of complete BUSCO genes per species (default: 5)",
    )
    parser.add_argument(
        "--flank-bp",
        type=int,
        default=5000,
        help="Flanking bp added upstream/downstream of each target gene (default: 5000)",
    )
    parser.add_argument(
        "--skip-genome-gff",
        action="store_true",
        help="Skip species_genome/species_gff extraction (CDS/manifest only).",
    )
    parser.add_argument(
        "--diamond-evalue",
        type=float,
        default=DEFAULT_THRESHOLDS["evalue"],
        help="DIAMOND evalue threshold (default: 1e-40)",
    )
    parser.add_argument(
        "--diamond-qcovhsp",
        type=float,
        default=DEFAULT_THRESHOLDS["qcovhsp"],
        help="DIAMOND qcovhsp threshold (default: 60)",
    )
    parser.add_argument(
        "--diamond-pident",
        type=float,
        default=DEFAULT_THRESHOLDS["pident"],
        help="DIAMOND pident threshold (default: 35)",
    )
    args = parser.parse_args()

    source_pg = Path(args.source_pg).resolve()
    out_pg = Path(args.out_pg).resolve()
    target_species = [s.strip() for s in args.species.split(",") if s.strip()]
    if not target_species:
        raise RuntimeError("No target species specified via --species")

    if out_pg.exists():
        shutil.rmtree(out_pg)
    out_pg.mkdir(parents=True, exist_ok=True)

    out_dirs = {
        "species_cds": out_pg / "species_cds",
        "species_expression": out_pg / "species_expression",
        "species_genome": out_pg / "species_genome",
        "species_gff": out_pg / "species_gff",
        "query_gene": out_pg / "query_gene",
        "manifest": out_pg / "manifest",
    }
    for d in out_dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    src_dirs = {
        "species_cds": source_pg / "species_cds",
        "species_expression": source_pg / "species_expression",
        "species_genome": source_pg / "species_genome",
        "species_gff": source_pg / "species_gff",
    }

    species_to_files: Dict[str, Dict[str, Path]] = {}
    for sp in target_species:
        species_to_files[sp] = {
            "cds": find_first_file(src_dirs["species_cds"], sp),
            "genome": find_first_file(src_dirs["species_genome"], sp),
            "gff": find_first_file(src_dirs["species_gff"], sp),
        }
        expr = sorted(src_dirs["species_expression"].glob(f"{sp}*.tsv"))
        species_to_files[sp]["expression"] = expr[0] if expr else None

    query_ids = load_query_ids(source_pg)
    for fam, ids in query_ids.items():
        (out_dirs["query_gene"] / fam).write_text("\n".join(ids) + "\n")

    busco_ids = load_busco_first_n(source_pg, target_species, args.busco_n)

    thresholds = {
        "evalue": float(args.diamond_evalue),
        "qcovhsp": float(args.diamond_qcovhsp),
        "pident": float(args.diamond_pident),
    }

    with tempfile.TemporaryDirectory(prefix="gg_miniset_") as td:
        tmpdir = Path(td)
        query_ref_cds = find_first_file(src_dirs["species_cds"], args.query_reference_species)
        fam_to_query_cds = extract_query_cds(
            query_ref_cds, query_ids, tmpdir
        )
        family_hits = run_family_diamond_hits(
            {sp: files["cds"] for sp, files in species_to_files.items()},
            fam_to_query_cds,
            thresholds,
            tmpdir,
        )

    selected_ids: Dict[str, Set[str]] = {}
    sources: Dict[str, Dict[str, Set[str]]] = {}
    for sp in target_species:
        fam_union = set()
        for fam in FAMILY_QUERIES:
            fam_union |= family_hits[sp][fam]
        selected = fam_union | set(busco_ids[sp])
        selected_ids[sp] = selected
        sources[sp] = {
            "AHA": set(family_hits[sp]["AHA"]),
            "YABBY": set(family_hits[sp]["YABBY"]),
            "BUSCO": set(busco_ids[sp]),
        }

    # Write per-species subsets.
    manifest_rows: List[Tuple[str, str, str]] = []
    summary_rows: List[Tuple[str, int, int, int, int, int, int, int]] = []
    for sp in target_species:
        files = species_to_files[sp]
        selected = set(selected_ids[sp])

        out_cds = out_dirs["species_cds"] / files["cds"].name
        found_ids, cds_count = write_trimmed_cds(files["cds"], selected, out_cds)
        out_gff = None

        expr_count = -1
        expr_file = files["expression"]
        if expr_file is not None:
            out_expr = out_dirs["species_expression"] / expr_file.name
            expr_count = write_trimmed_expression(expr_file, found_ids, out_expr)

        if not args.skip_genome_gff:
            intervals = collect_gene_intervals(files["gff"], found_ids, sp)
            windows = build_windows(intervals, flank=args.flank_bp)

            genome_base = files["genome"].name
            if genome_base.endswith(".gz"):
                genome_base = genome_base[: -len(".gz")]
            out_genome = out_dirs["species_genome"] / genome_base
            effective_windows = extract_genome_windows(files["genome"], windows, out_genome)

            gff_base = files["gff"].name
            if gff_base.endswith(".gz"):
                gff_base = gff_base[: -len(".gz")]
            out_gff = out_dirs["species_gff"] / gff_base
            write_shifted_gff(files["gff"], effective_windows, out_gff)

        if out_gff is not None:
            dummy_ids = append_species_synteny_dummies(
                species=sp,
                out_cds=out_cds,
                out_gff=out_gff,
                found_ids=found_ids,
            )
            if dummy_ids:
                selected |= dummy_ids
                found_ids |= dummy_ids
                cds_count += len(dummy_ids)
                sources[sp]["DUMMY"] = set(dummy_ids)
        missing_ids = sorted(selected - found_ids)

        aha_n = len(sources[sp]["AHA"])
        yabby_n = len(sources[sp]["YABBY"])
        busco_n = len(sources[sp]["BUSCO"])
        merged_n = len(selected)
        summary_rows.append(
            (
                sp,
                aha_n,
                yabby_n,
                busco_n,
                merged_n,
                cds_count,
                len(missing_ids),
                expr_count,
            )
        )

        for gid in sorted(selected):
            src = []
            for label in ("AHA", "YABBY", "BUSCO", "DUMMY"):
                if gid in sources[sp].get(label, set()):
                    src.append(label)
            manifest_rows.append((sp, gid, ";".join(src)))

    # Manifest outputs
    manifest_path = out_dirs["manifest"] / "selected_genes.tsv"
    with open(manifest_path, "w") as out:
        out.write("species\tgene_id\tsource\n")
        for sp, gid, src in sorted(manifest_rows):
            out.write(f"{sp}\t{gid}\t{src}\n")

    summary_path = out_dirs["manifest"] / "summary.tsv"
    with open(summary_path, "w") as out:
        out.write(
            "species\tAHA_hits\tYABBY_hits\tBUSCO_ids\tselected_union\tcds_written\tcds_missing\texpression_rows_written\n"
        )
        for row in summary_rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"Done. Output: {out_pg}")
    print(f"Summary: {summary_path}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
