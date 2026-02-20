#!/usr/bin/env python3

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import pandas

try:
    import fcntl

    HAS_FCNTL = True
except ImportError:
    HAS_FCNTL = False


FASTA_EXTENSIONS = (
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fas",
    ".fas.gz",
    ".fna",
    ".fna.gz",
)

GFF_EXTENSIONS = (
    ".gff",
    ".gff.gz",
    ".gff3",
    ".gff3.gz",
    ".gtf",
    ".gtf.gz",
)


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--focal_cds_fasta", required=True, type=str)
    parser.add_argument("--dir_sp_cds", required=True, type=str)
    parser.add_argument("--dir_sp_gff", required=True, type=str)
    parser.add_argument("--cache_dir", required=True, type=str)
    parser.add_argument("--gff2genestat_script", required=True, type=str)
    parser.add_argument("--window", default=5, type=int)
    parser.add_argument("--evalue", default=0.01, type=float)
    parser.add_argument("--genetic_code", default=1, type=int)
    parser.add_argument("--threads", default=1, type=int)
    parser.add_argument("--outfile", required=True, type=str)
    return parser


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fasta_ids(path):
    seq_ids = []
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seq_id = line[1:].strip().split()[0]
            if seq_id:
                seq_ids.append(seq_id)
    return seq_ids


def guess_species_name(gene_id):
    parts = gene_id.split("_")
    if len(parts) < 2:
        return ""
    return parts[0] + "_" + parts[1]


def find_species_file(directory, species_name, extensions):
    if not os.path.isdir(directory):
        return ""
    candidates = []
    for name in os.listdir(directory):
        if not name.startswith(species_name):
            continue
        lower_name = name.lower()
        if any(lower_name.endswith(ext) for ext in extensions):
            candidates.append(name)
    if not candidates:
        return ""
    candidates.sort()
    return os.path.join(directory, candidates[0])


def run_cmd(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        cmd_txt = " ".join(cmd)
        raise RuntimeError(
            "Command failed (exit code {}): {}\nSTDOUT:\n{}\nSTDERR:\n{}".format(
                proc.returncode, cmd_txt, proc.stdout, proc.stderr
            )
        )
    return proc


def load_gene_info(path):
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return pandas.DataFrame(columns=["gene_id", "chromosome", "start", "end", "strand"])
    df = pandas.read_csv(path, sep="\t", header=0)
    required_cols = {"gene_id", "chromosome", "start", "end"}
    if not required_cols.issubset(set(df.columns)):
        return pandas.DataFrame(columns=["gene_id", "chromosome", "start", "end", "strand"])
    out = df.loc[:, [c for c in ["gene_id", "chromosome", "start", "end", "strand"] if c in df.columns]].copy()
    if "strand" not in out.columns:
        out.loc[:, "strand"] = "+"
    out.loc[:, "gene_id"] = out["gene_id"].astype(str)
    out.loc[:, "chromosome"] = out["chromosome"].fillna("").astype(str)
    out.loc[:, "start"] = pandas.to_numeric(out["start"], errors="coerce")
    out.loc[:, "end"] = pandas.to_numeric(out["end"], errors="coerce")
    out = out.dropna(subset=["start", "end"])
    out.loc[:, "start"] = out[["start", "end"]].min(axis=1).astype(int)
    out.loc[:, "end"] = out[["start", "end"]].max(axis=1).astype(int)
    out = out.loc[out["chromosome"] != "", :]
    out = out.drop_duplicates(subset=["gene_id"], keep="first")
    out = out.sort_values(["chromosome", "start", "end", "gene_id"], kind="mergesort").reset_index(drop=True)
    return out


def ensure_species_gene_cache(
    species_name,
    species_cds_path,
    dir_sp_gff,
    cache_dir,
    gff2genestat_script,
    threads,
):
    ensure_parent_dir(os.path.join(cache_dir, "dummy"))
    out_path = os.path.join(cache_dir, species_name + ".gff_info.tsv")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    lock_path = out_path + ".lock"
    with open(lock_path, "a") as lock_handle:
        if HAS_FCNTL:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                return out_path
            tmp_path = out_path + ".tmp." + str(os.getpid())
            cmd = [
                sys.executable,
                gff2genestat_script,
                "--mode",
                "gene_delim",
                "--dir_gff",
                dir_sp_gff,
                "--seqfile",
                species_cds_path,
                "--feature",
                "CDS",
                "--multiple_hits",
                "longest",
                "--ncpu",
                str(max(1, int(threads))),
                "--outfile",
                tmp_path,
            ]
            run_cmd(cmd)
            if (not os.path.exists(tmp_path)) or os.path.getsize(tmp_path) == 0:
                raise RuntimeError("gff2genestat did not produce an output: {}".format(tmp_path))
            os.replace(tmp_path, out_path)
            return out_path
        finally:
            if HAS_FCNTL:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def parse_fasta_subset(path, wanted_ids):
    seqs = {}
    if len(wanted_ids) == 0:
        return seqs
    current_id = ""
    current_seq = []
    keep_current = False
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if keep_current and current_id and current_seq:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                keep_current = current_id in wanted_ids
                current_seq = []
                continue
            if keep_current:
                current_seq.append(line)
    if keep_current and current_id and current_seq:
        seqs[current_id] = "".join(current_seq)
    return seqs


def write_fasta(records, path):
    with open(path, "w") as handle:
        for seq_id in sorted(records.keys()):
            seq = records[seq_id]
            if not seq:
                continue
            handle.write(">{}\n".format(seq_id))
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


class UnionFind:
    def __init__(self, nodes):
        self.parent = {n: n for n in nodes}
        self.rank = {n: 0 for n in nodes}

    def find(self, x):
        parent = self.parent[x]
        if parent != x:
            self.parent[x] = self.find(parent)
        return self.parent[x]

    def union(self, x, y):
        rx = self.find(x)
        ry = self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1


def cluster_neighbors_by_similarity(cds_fasta, evalue_cutoff, genetic_code, threads, tmpdir):
    genes = read_fasta_ids(cds_fasta)
    genes = sorted(set(genes))
    if len(genes) == 0:
        return {}, {}
    uf = UnionFind(genes)
    if len(genes) == 1:
        return {genes[0]: "SG000001"}, {"SG000001": 1}
    pep_fasta = os.path.join(tmpdir, "neighbors.pep.fasta")
    db_prefix = os.path.join(tmpdir, "neighbors.pep")
    blast_out = os.path.join(tmpdir, "neighbors.blast.tsv")
    run_cmd(
        [
            "seqkit",
            "translate",
            "--allow-unknown-codon",
            "--transl-table",
            str(genetic_code),
            "--threads",
            str(max(1, int(threads))),
            cds_fasta,
            "--out-file",
            pep_fasta,
        ]
    )
    if (not os.path.exists(pep_fasta)) or os.path.getsize(pep_fasta) == 0:
        gene_to_group = {}
        group_size = {}
        for i, gene in enumerate(genes, start=1):
            gid = "SG{:06d}".format(i)
            gene_to_group[gene] = gid
            group_size[gid] = 1
        return gene_to_group, group_size
    run_cmd(["diamond", "makedb", "--in", pep_fasta, "--db", db_prefix])
    run_cmd(
        [
            "diamond",
            "blastp",
            "--query",
            pep_fasta,
            "--db",
            db_prefix,
            "--out",
            blast_out,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "evalue",
            "--evalue",
            str(evalue_cutoff),
            "--max-target-seqs",
            "2000",
            "--threads",
            str(max(1, int(threads))),
        ]
    )
    if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
        with open(blast_out, "r") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                qseqid = fields[0]
                sseqid = fields[1]
                if (qseqid == sseqid) or (qseqid not in uf.parent) or (sseqid not in uf.parent):
                    continue
                try:
                    evalue = float(fields[2])
                except ValueError:
                    continue
                if evalue <= evalue_cutoff:
                    uf.union(qseqid, sseqid)
    root_to_genes = defaultdict(list)
    for gene in genes:
        root_to_genes[uf.find(gene)].append(gene)
    components = []
    for members in root_to_genes.values():
        members_sorted = sorted(set(members))
        components.append(members_sorted)
    components.sort(key=lambda m: (-len(m), m[0]))
    gene_to_group = {}
    group_size = {}
    for i, members in enumerate(components, start=1):
        gid = "SG{:06d}".format(i)
        group_size[gid] = len(members)
        for gene in members:
            gene_to_group[gene] = gid
    return gene_to_group, group_size


def write_empty_output(path):
    ensure_parent_dir(path)
    cols = [
        "node_name",
        "species",
        "direction",
        "offset",
        "neighbor_gene",
        "group_id",
        "group_size",
    ]
    pandas.DataFrame(columns=cols).to_csv(path, sep="\t", index=False)


def main():
    args = build_arg_parser().parse_args()
    args.window = max(1, int(args.window))
    args.threads = max(1, int(args.threads))
    focal_gene_ids = sorted(set(read_fasta_ids(args.focal_cds_fasta)))
    if len(focal_gene_ids) == 0:
        write_empty_output(args.outfile)
        return
    focal_by_species = defaultdict(set)
    for gene_id in focal_gene_ids:
        species_name = guess_species_name(gene_id)
        if species_name:
            focal_by_species[species_name].add(gene_id)
    if len(focal_by_species) == 0:
        write_empty_output(args.outfile)
        return
    species_cds_paths = {}
    species_neighbor_ids = defaultdict(set)
    neighbor_rows = []
    for species_name in sorted(focal_by_species.keys()):
        sp_cds = find_species_file(args.dir_sp_cds, species_name, FASTA_EXTENSIONS)
        if not sp_cds:
            print("species CDS file not found: {}".format(species_name), file=sys.stderr)
            continue
        species_cds_paths[species_name] = sp_cds
        try:
            cache_path = ensure_species_gene_cache(
                species_name=species_name,
                species_cds_path=sp_cds,
                dir_sp_gff=args.dir_sp_gff,
                cache_dir=args.cache_dir,
                gff2genestat_script=args.gff2genestat_script,
                threads=args.threads,
            )
        except Exception as exc:
            print("failed to prepare gene cache for {}: {}".format(species_name, exc), file=sys.stderr)
            continue
        gene_info = load_gene_info(cache_path)
        if gene_info.shape[0] == 0:
            print("gene info is empty for species: {}".format(species_name), file=sys.stderr)
            continue
        chrom_genes = {}
        gene_locus = {}
        for chromosome, df_chr in gene_info.groupby("chromosome", sort=False):
            df_chr = df_chr.sort_values(["start", "end", "gene_id"], kind="mergesort").reset_index(drop=True)
            genes = df_chr["gene_id"].astype(str).tolist()
            strands = df_chr["strand"].fillna("+").astype(str).tolist()
            chrom_genes[chromosome] = {"genes": genes, "strands": strands}
            for i, gene_id in enumerate(genes):
                if gene_id not in gene_locus:
                    gene_locus[gene_id] = (chromosome, i, strands[i])
        for focal_gene in sorted(focal_by_species[species_name]):
            if focal_gene not in gene_locus:
                continue
            chromosome, focal_idx, strand = gene_locus[focal_gene]
            genes = chrom_genes[chromosome]["genes"]
            if strand == "-":
                upstream_sign = +1
                downstream_sign = -1
            else:
                upstream_sign = -1
                downstream_sign = +1
            for step in range(1, args.window + 1):
                up_idx = focal_idx + (upstream_sign * step)
                if 0 <= up_idx < len(genes):
                    neighbor_gene = genes[up_idx]
                    neighbor_rows.append(
                        {
                            "node_name": focal_gene,
                            "species": species_name,
                            "direction": "upstream",
                            "offset": -step,
                            "neighbor_gene": neighbor_gene,
                        }
                    )
                    species_neighbor_ids[species_name].add(neighbor_gene)
                down_idx = focal_idx + (downstream_sign * step)
                if 0 <= down_idx < len(genes):
                    neighbor_gene = genes[down_idx]
                    neighbor_rows.append(
                        {
                            "node_name": focal_gene,
                            "species": species_name,
                            "direction": "downstream",
                            "offset": step,
                            "neighbor_gene": neighbor_gene,
                        }
                    )
                    species_neighbor_ids[species_name].add(neighbor_gene)
    if len(neighbor_rows) == 0:
        write_empty_output(args.outfile)
        return
    neighbor_seqs = {}
    for species_name in sorted(species_neighbor_ids.keys()):
        if species_name not in species_cds_paths:
            continue
        ids = species_neighbor_ids[species_name]
        if len(ids) == 0:
            continue
        seqs = parse_fasta_subset(species_cds_paths[species_name], ids)
        neighbor_seqs.update(seqs)
    if len(neighbor_seqs) == 0:
        write_empty_output(args.outfile)
        return
    tmpdir = tempfile.mkdtemp(prefix="gg_synteny_neighbors_")
    try:
        cds_fasta = os.path.join(tmpdir, "neighbors.cds.fasta")
        write_fasta(neighbor_seqs, cds_fasta)
        gene_to_group, group_size = cluster_neighbors_by_similarity(
            cds_fasta=cds_fasta,
            evalue_cutoff=args.evalue,
            genetic_code=args.genetic_code,
            threads=args.threads,
            tmpdir=tmpdir,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    out_rows = []
    for row in neighbor_rows:
        neighbor_gene = row["neighbor_gene"]
        if neighbor_gene not in gene_to_group:
            continue
        row_out = dict(row)
        gid = gene_to_group[neighbor_gene]
        row_out["group_id"] = gid
        row_out["group_size"] = int(group_size.get(gid, 1))
        out_rows.append(row_out)
    ensure_parent_dir(args.outfile)
    if len(out_rows) == 0:
        write_empty_output(args.outfile)
        return
    df_out = pandas.DataFrame(out_rows)
    df_out = df_out.drop_duplicates(subset=["node_name", "direction", "offset", "neighbor_gene"])
    df_out = df_out.sort_values(["node_name", "direction", "offset", "neighbor_gene"], kind="mergesort")
    df_out.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
