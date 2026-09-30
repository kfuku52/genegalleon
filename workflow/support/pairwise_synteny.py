#!/usr/bin/env python3
"""Prepare and render pairwise JCVI synteny without species-tree dependencies."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq
from kffractbias.io import annotation_to_genes, natural_key, select_isoforms, write_bed

try:
    from fasta_sequence_store import fasta_records
    from species_labeling import extract_species_label
except ImportError:  # package imports in tests
    from .fasta_sequence_store import fasta_records
    from .species_labeling import extract_species_label

REQUIRED = ("analysis_id", "target_species", "query_species")
OPTIONAL = tuple(f"{side}_{field}" for side in ("target", "query") for field in
                 ("fasta", "gff", "feature", "attribute", "seqids"))
FASTA_SUFFIXES = (".fa", ".fas", ".fasta", ".fna", ".faa")
COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def safe_token(value, label):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"{label} must be a safe identifier: {value!r}")
    return value


def source_file(workspace, explicit, directory, species, suffixes):
    if explicit:
        path = Path(explicit)
        path = path if path.is_absolute() else workspace / path
    else:
        root = workspace / "input" / directory
        matches = [p for p in root.iterdir() if p.is_file()
                   and extract_species_label(p.name) == species
                   and p.name.removesuffix(".gz").endswith(suffixes)] if root.is_dir() else []
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one {directory} source for {species}; found {len(matches)}. Specify an explicit path in the pair table.")
        path = matches[0]
    if not path.is_file() or any(c in str(path) for c in "\n\r\0"):
        raise ValueError(f"Invalid source file: {path}")
    return str(path.resolve())


def tool_identity():
    import kffractbias.io

    versions = {name: importlib.metadata.version(name) for name in ("jcvi", "kfFractBias", "biopython")}
    versions["annotation_source_sha256"] = digest(kffractbias.io.__file__)
    versions["diamond"] = subprocess.check_output(["diamond", "version"], text=True).strip()
    return versions


def build_plan(args):
    if not 0 < args.cscore <= 1 or not 0 < args.minimum_mapping_fraction <= 1:
        raise ValueError("cscore and minimum-mapping-fraction must be in (0, 1]")
    if args.min_anchors < 2 or args.distance < 1:
        raise ValueError("min-anchors must be at least 2 and distance must be positive")
    formats = args.formats.split(",")
    if not formats or len(set(formats)) != len(formats) or set(formats) - {"pdf", "svg", "png"}:
        raise ValueError("formats must be a unique comma-separated subset of pdf,svg,png")
    workspace = args.workspace.resolve()
    pairs_file = args.pairs if args.pairs.is_absolute() else workspace / args.pairs
    with pairs_file.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or set(REQUIRED) - set(fields) or set(fields) - set(REQUIRED + OPTIONAL):
            raise ValueError("Pair table requires analysis_id,target_species,query_species and only documented optional columns")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Pair table row width does not match its header")
            if any(value.strip() for value in row.values()):
                rows.append({k: v.strip() for k, v in row.items()})
    if not rows:
        raise ValueError("Pair table has no analyses")
    identifiers = [row["analysis_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("analysis_id values must be unique")
    codes = {}
    code_file = workspace / "input/species_genetic_code/species_genetic_code.tsv"
    if code_file.is_file():
        with code_file.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["species"] in codes:
                    raise ValueError(f"Duplicate genetic-code species: {row['species']}")
                codes[row["species"]] = int(row["genetic_code"])
    if any(code not in CodonTable.generic_by_id for code in (args.genetic_code, *codes.values())):
        raise ValueError("Unknown NCBI genetic code")
    pairs = []
    for row in rows:
        pair = {key: safe_token(row[key], key) for key in REQUIRED}
        if pair["target_species"] == pair["query_species"]:
            raise ValueError("Pairwise synteny requires different species; use the existing self-synteny workflow for self comparisons")
        for side in ("target", "query"):
            species = pair[f"{side}_species"]
            mode = args.sequence_mode
            if mode == "auto":
                if row.get(f"{side}_fasta"):
                    raise ValueError("Explicit FASTA paths require sequence-mode protein or cds")
                protein_dir = workspace / "input/species_protein"
                mode = "protein" if protein_dir.is_dir() and any(
                    p.is_file() and p.name.removesuffix(".gz").endswith(FASTA_SUFFIXES)
                    and extract_species_label(p.name) == species for p in protein_dir.iterdir()
                ) else "cds"
            feature, attribute = row.get(f"{side}_feature", ""), row.get(f"{side}_attribute", "")
            if bool(feature) != bool(attribute):
                raise ValueError("GFF feature and attribute must be specified together")
            pair[side] = {
                "species": species, "mode": mode,
                "fasta": source_file(workspace, row.get(f"{side}_fasta"), f"species_{mode}", species, FASTA_SUFFIXES),
                "gff": source_file(workspace, row.get(f"{side}_gff"), "species_gff", species, (".gff", ".gff3", ".gtf")),
                "feature": feature, "attribute": attribute,
                "genetic_code": codes.get(species, args.genetic_code) if mode == "cds" else None,
            }
            pair[f"{side}_seqids"] = row.get(f"{side}_seqids", "")
        pairs.append(pair)
    inputs = {str(Path(__file__).resolve()): digest(__file__)}
    for pair in pairs:
        for side in ("target", "query"):
            for kind in ("fasta", "gff"):
                path = pair[side][kind]
                inputs[path] = digest(path)
    return {
        "workspace": str(workspace), "pairs": sorted(pairs, key=lambda p: p["analysis_id"]),
        "parameters": {"cscore": args.cscore, "min_anchors": args.min_anchors, "distance": args.distance,
                       "minimum_mapping_fraction": args.minimum_mapping_fraction, "quota": None,
                       "isoform_policy": "longest"},
        "formats": formats, "tools": tool_identity(), "input_hashes": inputs, "schema_version": 1,
    }


def verify_inputs(plan):
    for path, expected in plan["input_hashes"].items():
        if digest(path) != expected:
            raise ValueError(f"Input changed during pairwise synteny; previous results were preserved: {path}")


def phase_root(plan, phase):
    return Path(plan["workspace"]) / "output/genome_evolution/synteny" / phase


def contract_args(plan, phase):
    workspace = Path(plan["workspace"])
    provenance = workspace / "output/artifact_provenance/genome_evolution" / f"pairwise_synteny.{phase}.json"
    result = ["--manifest", str(provenance), "--step", f"pairwise_synteny_{phase}", "--family-id", "all_pairs",
              "--logical-root", str(workspace / "output/.gg_global_artifacts"), "--workspace-root", str(workspace),
              "--output", f"results={phase_root(plan, phase)}", "--input", f"implementation={Path(__file__).resolve()}"]
    if phase == "analysis":
        pairs = [{k: v for k, v in pair.items() if not k.endswith("_seqids")} for pair in plan["pairs"]]
        for pair in pairs:
            for side in ("target", "query"):
                for kind in ("fasta", "gff"):
                    result.extend(("--input", f"{pair['analysis_id']}.{side}.{kind}={pair[side][kind]}"))
        parameters = {"pairs": pairs, **plan["parameters"], "tools": plan["tools"]}
    else:
        result.extend(("--input", f"analysis={phase_root(plan, 'analysis')}"))
        parameters = {"formats": plan["formats"], "coordinate_system": "gene_rank",
                      "pairs": [{k: pair[k] for k in ("analysis_id", "target_seqids", "query_seqids")} for pair in plan["pairs"]],
                      "jcvi": plan["tools"]["jcvi"]}
    result.extend(("--parameter", "configuration=" + json.dumps(parameters, sort_keys=True)))
    return result


def prepare_genome(source, directory, side, minimum_mapping_fraction):
    sequences = {}
    aliases = {}
    species = source["species"]
    for identifier, _header, sequence in fasta_records(Path(source["fasta"])):
        safe_token(identifier, "FASTA identifier")
        if identifier in sequences or not sequence:
            raise ValueError(f"Duplicate or empty FASTA record: {identifier}")
        alphabet = "ACGTURYSWKMBDHVN" if source["mode"] == "cds" else "ACDEFGHIKLMNPQRSTVWYBXZJUO*"
        if set(sequence.upper()) - set(alphabet):
            raise ValueError(f"Invalid {source['mode']} sequence: {identifier}")
        sequences[identifier] = sequence.upper()
        for alias in {identifier, identifier.removeprefix(species + "_")}:
            if alias in aliases and aliases[alias] != identifier:
                raise ValueError(f"Ambiguous species-prefix alias: {alias}")
            aliases[alias] = identifier
    if not sequences:
        raise ValueError(f"Empty FASTA: {source['fasta']}")
    mapping = annotation_to_genes(source["gff"], set(aliases), feature=source["feature"] or None,
                                  attribute=source["attribute"] or None)
    genes = tuple(replace(gene, gene_id=aliases[gene.gene_id]) for gene in mapping.genes)
    if len({gene.gene_id for gene in genes}) != len(genes):
        raise ValueError("Multiple GFF identifiers map to the same FASTA record")
    mapping = replace(mapping, genes=genes, fasta_gene_count=len(sequences),
                      locus_by_id={aliases[key]: value for key, value in mapping.locus_by_id.items()})
    if len(genes) / len(sequences) < minimum_mapping_fraction:
        raise ValueError(f"Only {len(genes)}/{len(sequences)} FASTA identifiers mapped to {source['gff']}")
    mapping = select_isoforms(mapping, {key: len(value) for key, value in sequences.items()}, "longest")
    selected = {gene.gene_id for gene in mapping.genes}
    matched = {gene.gene_id for gene in genes}
    names = {identifier: species + "_" + identifier.removeprefix(species + "_") for identifier in selected}
    normalized = tuple(replace(gene, gene_id=names[gene.gene_id]) for gene in mapping.genes)
    write_bed(normalized, directory / f"{side}.bed")
    with (directory / f"{side}.pep").open("w", encoding="utf-8") as handle:
        for gene in mapping.genes:
            sequence = sequences[gene.gene_id]
            if source["mode"] == "cds":
                if len(sequence) % 3:
                    raise ValueError(f"CDS length is not divisible by three: {gene.gene_id}")
                sequence = str(Seq(sequence).translate(table=source["genetic_code"]))
            sequence = sequence.removesuffix("*")
            if not sequence or "*" in sequence:
                raise ValueError(f"Empty translation or internal stop: {gene.gene_id}")
            handle.write(f">{names[gene.gene_id]}\n{sequence}\n")
    write_tsv(directory / f"{side}.id_map.tsv", ("original_id", "locus_id", "jcvi_id", "status"),
              ((identifier, mapping.locus_by_id.get(identifier, ""), names.get(identifier, ""),
                "selected" if identifier in selected else "isoform_excluded" if identifier in matched else "unmapped")
               for identifier in sorted(sequences, key=natural_key)))
    metadata = {**mapping.metadata(), "source": source, "fasta_sha256": digest(source["fasta"]),
                "gff_sha256": digest(source["gff"]), "unmapped_count": len(sequences) - len(genes)}
    return normalized, metadata


def run_command(command, cwd, commands, label):
    logs = cwd / "logs"
    logs.mkdir(exist_ok=True)
    log = logs / f"{len(commands) + 1:02d}.{label}.log"
    print(f"Pairwise synteny: {label}; log: {log}", file=sys.stderr, flush=True)
    env = {**os.environ, "MPLCONFIGDIR": str(cwd / ".mplconfig"), "MPLBACKEND": "Agg"}
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT, check=False)
    commands.append({"argv": command, "returncode": result.returncode, "log": str(log.relative_to(cwd))})
    write_json(cwd / "commands.json", commands)
    if result.returncode:
        raise RuntimeError(f"{label} failed ({result.returncode}):\n{log.read_text(encoding='utf-8')[-8000:]}")


def summarize_anchors(directory, genomes):
    lookup = [{g.gene_id: g for g in genome} for genome in genomes]
    blocks = []
    block = []
    for line in (directory / "target.query.lifted.anchors").read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            if block:
                blocks.append(block)
            block = []
        elif line.strip():
            fields = line.split()
            if len(fields) < 3 or fields[0] not in lookup[0] or fields[1] not in lookup[1]:
                raise ValueError(f"Invalid JCVI anchor: {line}")
            block.append(tuple(fields[:3]))
    if block:
        blocks.append(block)
    if not blocks:
        raise ValueError("No syntenic blocks were found; no completed plots will be published")
    rows, block_rows = [], []
    for number, block in enumerate(blocks, 1):
        target_genes = [lookup[0][left] for left, _, _ in block]
        query_genes = [lookup[1][right] for _, right, _ in block]
        if len({g.seqid for g in target_genes}) != 1 or len({g.seqid for g in query_genes}) != 1:
            raise ValueError("A JCVI block spans multiple chromosomes")
        orientation = "+" if (target_genes[-1].start - target_genes[0].start) * (query_genes[-1].start - query_genes[0].start) >= 0 else "-"
        block_rows.append((number, target_genes[0].seqid, min(g.start for g in target_genes), max(g.end for g in target_genes),
                           query_genes[0].seqid, min(g.start for g in query_genes), max(g.end for g in query_genes), orientation, len(block)))
        rows.extend((number, left, right, score) for left, right, score in block)
    write_tsv(directory / "anchors.tsv", ("block_id", "target_gene", "query_gene", "score"), rows)
    write_tsv(directory / "blocks.tsv", ("block_id", "target_seqid", "target_start", "target_end", "query_seqid", "query_start", "query_end", "orientation", "anchor_count"), block_rows)
    return {"block_count": len(blocks), "anchor_count": len(rows),
            "coordinate_system": "BED: 0-based half-open", "syntenic_genes":
            [len({row[index] for row in rows}) for index in (1, 2)],
            "syntenic_gene_fraction": [len({row[index + 1] for row in rows}) / len(genome) for index, genome in enumerate(genomes)]}


def analyze(plan, output, cpus):
    output.mkdir(parents=True)
    for pair in plan["pairs"]:
        directory = output / pair["analysis_id"]
        directory.mkdir()
        target, target_info = prepare_genome(pair["target"], directory, "target", plan["parameters"]["minimum_mapping_fraction"])
        query, query_info = prepare_genome(pair["query"], directory, "query", plan["parameters"]["minimum_mapping_fraction"])
        if {g.gene_id for g in target} & {g.gene_id for g in query}:
            raise ValueError("Target and query normalized identifiers overlap")
        parameters = plan["parameters"]
        commands = []
        run_command([sys.executable, "-m", "jcvi.compara.catalog", "ortholog", "target", "query", "--dbtype=prot",
                     "--align_soft=diamond_blastp", "--no_strip_names", "--no_dotplot", f"--cpus={cpus}",
                     f"--cscore={parameters['cscore']}", f"--min_size={parameters['min_anchors']}", f"--dist={parameters['distance']}"],
                    directory, commands, "mcscan")
        counts = summarize_anchors(directory, (target, query))
        run_command([sys.executable, "-m", "jcvi.compara.synteny", "screen", "--minspan=0", "--minsize=0", "--simple",
                     "target.query.lifted.anchors", "target.query.screened.anchors"], directory, commands, "blocks")
        write_json(directory / "summary.json", {"schema_version": 1, "analysis_id": pair["analysis_id"], "target": target_info,
                   "query": query_info, "parameters": parameters, "tools": plan["tools"], **counts})


def selected_seqids(value, genes):
    available = sorted({gene.seqid for gene in genes}, key=natural_key)
    selected = value.split(",") if value else available
    if len(selected) != len(set(selected)) or set(selected) - set(available):
        raise ValueError(f"Display seqids must be unique annotated chromosomes; requested {selected}, available {available}")
    if any(any(char in seqid for char in ",\n\r") or seqid.endswith("-") for seqid in selected):
        raise ValueError("JCVI display seqids cannot contain commas/newlines or end with '-' (an orientation control)")
    return selected


def render(plan, output):
    from kffractbias.io import read_bed

    output.mkdir(parents=True)
    for pair in plan["pairs"]:
        source = phase_root(plan, "analysis") / pair["analysis_id"]
        directory = output / pair["analysis_id"]
        directory.mkdir()
        for name in ("target.bed", "query.bed", "target.query.lifted.anchors", "target.query.screened.simple"):
            shutil.copyfile(source / name, directory / name)
        genomes = [read_bed(directory / f"{side}.bed") for side in ("target", "query")]
        selected = [selected_seqids(pair[f"{side}_seqids"], genes) for side, genes in zip(("target", "query"), genomes, strict=True)]
        target_by_id = {gene.gene_id: gene.seqid for gene in genomes[0]}
        colors = {seqid: COLORS[index % len(COLORS)] for index, seqid in enumerate(sorted(set(target_by_id.values()), key=natural_key))}
        with (directory / "colored.simple").open("w", encoding="utf-8") as handle:
            for line in (directory / "target.query.screened.simple").read_text(encoding="utf-8").splitlines():
                if line.strip():
                    handle.write(colors[target_by_id[line.split()[0]]] + "*" + line + "\n")
        (directory / "seqids").write_text("\n".join(",".join(seqids) for seqids in selected) + "\n", encoding="utf-8")
        (directory / "layout").write_text(
            "# y, xstart, xend, rotation, color, label, va, bed, label_va\n"
            f"0.7,0.12,0.92,0,,{pair['target_species'].replace('_', ' ')} (gene rank),top,target.bed,top\n"
            f"0.3,0.12,0.92,0,,{pair['query_species'].replace('_', ' ')} (gene rank),bottom,query.bed,bottom\n"
            "# edges\ne,0,1,colored.simple\n", encoding="utf-8")
        write_tsv(directory / "display.tsv", ("species", "seqid", "selected_for_karyotype", "gene_count"),
                  ((pair[f"{side}_species"], seqid, int(seqid in seqids), count)
                   for side, genes, seqids in zip(("target", "query"), genomes, selected, strict=True)
                   for seqid, count in sorted(Counter(g.seqid for g in genes).items(), key=lambda item: natural_key(item[0]))))
        commands = []
        for fmt in plan["formats"]:
            common = ["--notex", f"--format={fmt}", "--seed=1"]
            run_command([sys.executable, "-m", "jcvi.graphics.dotplot", "target.query.lifted.anchors", "--nosort", "--nochpf", "--colororientation",
                         f"--nmax={json.loads((source / 'summary.json').read_text(encoding='utf-8'))['anchor_count']}",
                         f"--genomenames={pair['target_species'].replace('_', ' ')}_{pair['query_species'].replace('_', ' ')}",
                         "--title=Pairwise synteny (gene rank)", "--style=white", f"--outfile=dotplot.{fmt}", *common], directory, commands, f"dotplot-{fmt}")
            run_command([sys.executable, "-m", "jcvi.graphics.karyotype", "seqids", "layout", "--keep-chrlabels",
                         "--figsize=12x7", f"--outfile=karyotype.{fmt}", *common], directory, commands, f"karyotype-{fmt}")
            for name in ("dotplot", "karyotype"):
                if not (directory / f"{name}.{fmt}").is_file() or not (directory / f"{name}.{fmt}").stat().st_size:
                    raise ValueError(f"Missing or empty {name}.{fmt}")
        shutil.rmtree(directory / ".mplconfig", ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    planning = subparsers.add_parser("plan")
    planning.add_argument("--workspace", required=True, type=Path)
    planning.add_argument("--pairs", required=True, type=Path)
    planning.add_argument("--sequence-mode", choices=("auto", "protein", "cds"), default="auto")
    planning.add_argument("--genetic-code", type=int, default=1)
    planning.add_argument("--cscore", type=float, default=0.7)
    planning.add_argument("--min-anchors", type=int, default=4)
    planning.add_argument("--distance", type=int, default=20)
    planning.add_argument("--minimum-mapping-fraction", type=float, default=1)
    planning.add_argument("--formats", default="pdf,svg,png")
    planning.add_argument("--outfile", required=True, type=Path)
    contract = subparsers.add_parser("contract")
    contract.add_argument("--plan", required=True, type=Path)
    contract.add_argument("--phase", choices=("analysis", "plots"), required=True)
    for phase in ("analysis", "plots"):
        execution = subparsers.add_parser(phase)
        execution.add_argument("--plan", required=True, type=Path)
        execution.add_argument("--output", required=True, type=Path)
        execution.add_argument("--cpus", type=int, default=1)
    verification = subparsers.add_parser("verify")
    verification.add_argument("--plan", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "plan":
            write_json(args.outfile, build_plan(args))
        else:
            plan = json.loads(args.plan.read_text(encoding="utf-8"))
            if args.command == "verify":
                verify_inputs(plan)
            elif args.command == "contract":
                sys.stdout.buffer.write(b"\0".join(value.encode() for value in contract_args(plan, args.phase)) + b"\0")
            elif args.cpus < 1:
                raise ValueError("cpus must be positive")
            elif args.command == "analysis":
                analyze(plan, args.output, args.cpus)
            else:
                render(plan, args.output)
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"Pairwise synteny failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
