#!/usr/bin/env python3
"""Prepare and render pairwise JCVI synteny without species-tree dependencies."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq
from kffractbias.io import annotation_to_genes, natural_key, select_isoforms, write_bed

try:
    from fasta_sequence_store import fasta_records
    from pairwise_synteny_dotplot import chromosome_lengths, prepare_dotplot
    from pairwise_synteny_karyotype import chromosome_colors
    from pairwise_synteny_layout import FIGSIZE, order_by_ribbon_length
    from representative_selection import (
        identifier_aliases,
        load_effective_inputs,
        load_representative_map,
        matching_identifier,
    )
    from species_labeling import extract_species_label
except ImportError:  # package imports in tests
    from .fasta_sequence_store import fasta_records
    from .pairwise_synteny_dotplot import chromosome_lengths, prepare_dotplot
    from .pairwise_synteny_karyotype import chromosome_colors
    from .pairwise_synteny_layout import FIGSIZE, order_by_ribbon_length
    from .representative_selection import (
        identifier_aliases,
        load_effective_inputs,
        load_representative_map,
        matching_identifier,
    )
    from .species_labeling import extract_species_label

REQUIRED = ("analysis_id", "target_species", "query_species")
OPTIONAL = tuple(f"{side}_{field}" for side in ("target", "query") for field in
                 ("fasta", "gff", "feature", "attribute", "seqids", "cds", "genome", "sizes"))
FASTA_SUFFIXES = (".fa", ".fas", ".fasta", ".fna", ".faa")


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
    dotplot_color = getattr(args, "dotplot_color", "orientation")
    ds_color_max = getattr(args, "ds_color_max", 2.0)
    minimum_length = getattr(args, "dotplot_min_length", 1_000_000)
    dotplot_sort = getattr(args, "dotplot_sort", "homoeolog")
    karyotype_color = getattr(args, "karyotype_color", "chromosome")
    if karyotype_color not in {"chromosome", "homoeolog"}:
        raise ValueError("karyotype-color must be chromosome or homoeolog")
    karyotype_scale = getattr(args, "karyotype_scale", "shared")
    if karyotype_scale not in {"shared", "independent"}:
        raise ValueError("karyotype-scale must be shared or independent")
    karyotype_track_order = getattr(args, "karyotype_track_order", "target-query")
    if karyotype_track_order not in {"target-query", "query-target"}:
        raise ValueError("karyotype-track-order must be target-query or query-target")
    if minimum_length < 0 or dotplot_sort not in {"karyotype", "homoeolog", "none"}:
        raise ValueError("dotplot-min-length must be nonnegative; dotplot-sort must be karyotype, homoeolog or none")
    if dotplot_color not in {"orientation", "ds"} or not math.isfinite(ds_color_max) or ds_color_max <= 0:
        raise ValueError("dotplot-color must be orientation or ds; ds-color-max must be finite and positive")
    if not 0 < args.cscore <= 1 or not 0 < args.minimum_mapping_fraction <= 1:
        raise ValueError("cscore and minimum-mapping-fraction must be in (0, 1]")
    if args.min_anchors < 2 or args.distance < 1:
        raise ValueError("min-anchors must be at least 2 and distance must be positive")
    formats = args.formats.split(",")
    if not formats or len(set(formats)) != len(formats) or set(formats) - {"pdf", "svg", "png"}:
        raise ValueError("formats must be a unique comma-separated subset of pdf,svg,png")
    workspace = args.workspace.resolve()
    effective_manifest = getattr(args, "representative_inputs", "")
    effective_inputs = {}
    if effective_manifest:
        effective_manifest = Path(effective_manifest)
        effective_manifest = (effective_manifest if effective_manifest.is_absolute()
                              else workspace / effective_manifest).resolve()
        effective_inputs = load_effective_inputs(effective_manifest)
        if any(int(row["genetic_code"]) not in CodonTable.generic_by_id for row in effective_inputs.values()):
            raise ValueError("Effective inputs manifest has an unknown NCBI genetic code")
    representative_path = getattr(args, "representative_map", "")
    if representative_path:
        representative_path = Path(representative_path)
        representative_path = (representative_path if representative_path.is_absolute()
                               else workspace / representative_path).resolve()
        load_representative_map(representative_path)
    pairs_file = args.pairs if args.pairs.is_absolute() else workspace / args.pairs
    pairs_digest = digest(pairs_file)
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
    code_digest = None
    code_file = workspace / "input/species_genetic_code/species_genetic_code.tsv"
    if code_file.is_file():
        code_digest = digest(code_file)
        with code_file.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if len(reader.fieldnames or ()) != 2 or set(reader.fieldnames or ()) != {"species", "genetic_code"}:
                raise ValueError("Genetic-code table requires species and genetic_code columns")
            for row in reader:
                if None in row or any(value is None or not value.strip() for value in row.values()):
                    raise ValueError("Invalid genetic-code table row")
                row = {key: value.strip() for key, value in row.items()}
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
            effective = effective_inputs.get(species)
            if effective_inputs and effective is None:
                raise ValueError(f"Effective inputs manifest lacks species: {species}")
            if effective and any(row.get(f"{side}_{role}") for role in ("fasta", "gff", "cds", "genome")):
                raise ValueError("Effective inputs manifest cannot be combined with explicit pair source paths")
            mode = args.sequence_mode
            if mode == "auto":
                if row.get(f"{side}_fasta"):
                    raise ValueError("Explicit FASTA paths require sequence-mode protein or cds")
                protein_dir = workspace / "input/species_protein"
                mode = "protein" if effective or (protein_dir.is_dir() and any(
                    p.is_file() and p.name.removesuffix(".gz").endswith(FASTA_SUFFIXES)
                    and extract_species_label(p.name) == species for p in protein_dir.iterdir()
                )) else "cds"
            feature, attribute = row.get(f"{side}_feature", ""), row.get(f"{side}_attribute", "")
            if bool(feature) != bool(attribute):
                raise ValueError("GFF feature and attribute must be specified together")
            pair[side] = {
                "species": species, "mode": mode,
                "fasta": (effective.get('analysis_cds', effective['cds']) if mode == 'cds' else effective[mode]) if effective else source_file(workspace, row.get(f"{side}_fasta"), f"species_{mode}", species, FASTA_SUFFIXES),
                "gff": (effective.get('analysis_gff', effective['gff']) if mode == 'cds' else effective['gff']) if effective else source_file(workspace, row.get(f"{side}_gff"), "species_gff", species, (".gff", ".gff3", ".gtf")),
                "feature": feature, "attribute": attribute,
                "genetic_code": (int(effective["genetic_code"]) if effective else codes.get(species, args.genetic_code)) if mode == "cds" else None,
            }
            if representative_path or effective:
                if effective and representative_path and str(representative_path) != effective["representative_map"]:
                    raise ValueError("Explicit representative map differs from effective inputs manifest")
                pair[side]["representative_map"] = effective["representative_map"] if effective else str(representative_path)
            if effective and mode == 'cds':
                # Raw biological CDS can be partial or deliberately excluded
                # from translation. Use the bundle's admitted coding view.
                pair[side]['admitted_protein'] = effective['protein']
            pair[f"{side}_seqids"] = row.get(f"{side}_seqids", "")
            genome, sizes = row.get(f"{side}_genome"), row.get(f"{side}_sizes")
            genome = effective["genome"] if effective else genome
            if genome and sizes:
                raise ValueError("Specify either genome or sizes for each species, not both")
            kind, length_path = "gff", pair[side]["gff"]
            if sizes or genome:
                kind = "sizes" if sizes else "genome"
                length_path = source_file(workspace, sizes or genome, "species_genome", species, FASTA_SUFFIXES)
            elif minimum_length:
                root = workspace / "input/species_genome"
                matches = [p for p in root.iterdir() if p.is_file() and extract_species_label(p.name) == species
                           and p.name.removesuffix(".gz").endswith(FASTA_SUFFIXES)] if root.is_dir() else []
                if matches:
                    kind = "genome"
                    length_path = source_file(workspace, "", "species_genome", species, FASTA_SUFFIXES)
            pair.setdefault("dotplot_lengths", {})[side] = {"kind": kind, "path": length_path}
        if dotplot_color == "ds":
            pair["ds"] = {}
            for side in ("target", "query"):
                source = pair[side]
                effective = effective_inputs.get(source["species"])
                path = effective.get('analysis_cds', effective['cds']) if effective else source_file(workspace, row.get(f"{side}_cds") or (source["fasta"] if source["mode"] == "cds" else ""),
                                   "species_cds", source["species"], FASTA_SUFFIXES)
                pair["ds"][side] = {"fasta": path, "genetic_code": int(effective["genetic_code"]) if effective else codes.get(source["species"], args.genetic_code)}
            if pair["ds"]["target"]["genetic_code"] != pair["ds"]["query"]["genetic_code"]:
                raise ValueError("dS estimation requires the same genetic code for both species; use orientation coloring for mixed-code pairs")
        pairs.append(pair)
    layout_source = Path(__file__).with_name("pairwise_synteny_layout.py").resolve()
    inputs = {str(Path(__file__).resolve()): digest(__file__), str(layout_source): digest(layout_source),
              str(pairs_file.resolve()): pairs_digest}
    dotplot_source = Path(__file__).with_name("pairwise_synteny_dotplot.py").resolve()
    inputs[str(dotplot_source)] = digest(dotplot_source)
    for helper in ("pairwise_synteny_karyotype.py", "pairwise_synteny_style.py", "pairwise_synteny_layout.py"):
        path = Path(__file__).with_name(helper).resolve()
        inputs[str(path)] = digest(path)
    reader_source = Path(__file__).with_name("fasta_sequence_store.py").resolve()
    inputs[str(reader_source)] = digest(reader_source)
    if representative_path or effective_inputs:
        if representative_path:
            inputs[str(representative_path)] = digest(representative_path)
        if effective_manifest:
            inputs[str(effective_manifest)] = digest(effective_manifest)
            for sources in effective_inputs.values():
                for role in ("cds", "protein", "gff", "genome", "representative_map"):
                    inputs[sources[role]] = sources[role + "_sha256"]
                for role in ('analysis_cds', 'analysis_gff'):
                    if role in sources:
                        inputs[sources[role]] = sources[role + '_sha256']
        for helper in ("representative_selection.py", "gff2genestat.py", "gff_feature_structure.py"):
            path = Path(__file__).with_name(helper).resolve()
            inputs[str(path)] = digest(path)
    if code_digest is not None:
        inputs[str(code_file.resolve())] = code_digest
    for pair in pairs:
        for side in ("target", "query"):
            for kind in ("fasta", "gff"):
                path = pair[side][kind]
                inputs[path] = digest(path)
            path = pair["dotplot_lengths"][side]["path"]
            inputs[path] = digest(path)
    ds_tools = {}
    if dotplot_color == "ds":
        try:
            from pairwise_synteny_ds import ds_tool_identity
        except ImportError:
            from .pairwise_synteny_ds import ds_tool_identity
        ds_tools = ds_tool_identity()
        inputs.update(ds_tools["source_hashes"])
        for pair in pairs:
            for side in ("target", "query"):
                inputs[pair["ds"][side]["fasta"]] = digest(pair["ds"][side]["fasta"])
    plan = {
        "workspace": str(workspace), "pairs": sorted(pairs, key=lambda p: p["analysis_id"]),
        "parameters": {"cscore": args.cscore, "min_anchors": args.min_anchors, "distance": args.distance,
                       "minimum_mapping_fraction": args.minimum_mapping_fraction, "quota": None,
                       "isoform_policy": "representative_map" if representative_path or effective_inputs else "longest"},
        "formats": formats, "karyotype_sort": args.karyotype_sort, "karyotype_color": karyotype_color,
        "karyotype_scale": karyotype_scale,
        "karyotype_track_order": karyotype_track_order,
        "dotplot_color": dotplot_color, "ds_color_max": ds_color_max, "ds_tools": ds_tools,
        "dotplot_min_length": minimum_length, "dotplot_sort": dotplot_sort,
        "tools": tool_identity(), "input_hashes": inputs, "schema_version": 1,
    }
    if effective_manifest:
        plan["representative_inputs"] = str(effective_manifest)
    return plan


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
              "--output", f"results={phase_root(plan, phase)}", "--input", f"implementation={Path(__file__).resolve()}",
              "--input", f"sequence_reader={Path(__file__).with_name('fasta_sequence_store.py').resolve()}"]
    if phase == "analysis":
        if plan.get("representative_inputs"):
            result.extend(("--input", "representative_inputs=" + plan["representative_inputs"]))
        pairs = [{k: v for k, v in pair.items() if not k.endswith("_seqids") and k not in {"ds", "dotplot_lengths"}} for pair in plan["pairs"]]
        for pair in pairs:
            for side in ("target", "query"):
                for kind in ("fasta", "gff"):
                    result.extend(("--input", f"{pair['analysis_id']}.{side}.{kind}={pair[side][kind]}"))
                if pair[side].get('admitted_protein'):
                    result.extend(('--input', f"{pair['analysis_id']}.{side}.admitted_protein="
                                   + pair[side]['admitted_protein']))
                if pair[side].get("representative_map"):
                    result.extend(("--input", f"{pair['analysis_id']}.{side}.representative_map="
                                   + pair[side]["representative_map"]))
        if plan["parameters"].get("isoform_policy") == "representative_map":
            for helper in ("representative_selection.py", "gff2genestat.py", "gff_feature_structure.py"):
                result.extend(("--input", f"selection_implementation.{helper}="
                               + str(Path(__file__).with_name(helper).resolve())))
        parameters = {"pairs": pairs, **plan["parameters"], "tools": plan["tools"]}
    elif phase == "ds":
        result.extend(("--input", f"analysis={phase_root(plan, 'analysis')}"))
        for path in plan["ds_tools"]["source_hashes"]:
            result.extend(("--input", f"ds_implementation.{Path(path).name}={path}"))
        for pair in plan["pairs"]:
            for side in ("target", "query"):
                result.extend(("--input", f"{pair['analysis_id']}.{side}.cds={pair['ds'][side]['fasta']}"))
        parameters = {"tools": plan["ds_tools"], "pairs": [{"analysis_id": p["analysis_id"], "ds": p["ds"]} for p in plan["pairs"]]}
    else:
        result.extend(("--input", f"analysis={phase_root(plan, 'analysis')}"))
        result.extend(("--input", f"layout_implementation={Path(__file__).with_name('pairwise_synteny_layout.py').resolve()}"))
        result.extend(("--input", f"dotplot_implementation={Path(__file__).with_name('pairwise_synteny_dotplot.py').resolve()}"))
        for helper in ("pairwise_synteny_karyotype.py", "pairwise_synteny_style.py"):
            result.extend(("--input", f"{helper}={Path(__file__).with_name(helper).resolve()}"))
        for pair in plan["pairs"]:
            for side, source in pair.get("dotplot_lengths", {}).items():
                result.extend(("--input", f"{pair['analysis_id']}.{side}.lengths={source['path']}"))
        parameters = {"formats": plan["formats"], "coordinate_system": "gene_rank",
                      "karyotype_sort": plan.get("karyotype_sort", "both_length"), "karyotype_figsize": list(FIGSIZE),
                      "karyotype_color": plan.get("karyotype_color", "chromosome"),
                      "karyotype_scale": plan.get("karyotype_scale", "shared"),
                      "karyotype_track_order": plan.get("karyotype_track_order", "target-query"),
                      "pairs": [{k: pair[k] for k in ("analysis_id", "target_seqids", "query_seqids")} for pair in plan["pairs"]],
                      "jcvi": plan["tools"]["jcvi"]}
        parameters.update(dotplot_color=plan.get("dotplot_color", "orientation"), ds_color_max=plan.get("ds_color_max", 2.0))
        parameters.update(dotplot_min_length=plan.get("dotplot_min_length", 1_000_000),
                          dotplot_sort=plan.get("dotplot_sort", "homoeolog"),
                          length_sources=[p.get("dotplot_lengths", {}) for p in plan["pairs"]])
        if plan.get("dotplot_color") == "ds":
            result.extend(("--input", f"ds={phase_root(plan, 'ds')}"))
            result.extend(("--input", f"ds_renderer={Path(__file__).with_name('pairwise_synteny_ds.py').resolve()}"))
    result.extend(("--parameter", "configuration=" + json.dumps(parameters, sort_keys=True)))
    return result


def select_declared_isoforms(mapping, representative_map, source, *, coding_traits=None):
    """Bind both candidate FASTAs and locus FASTAs to the same declared CDS."""
    try:
        from gff2genestat import process_single_gff
    except ImportError:  # package imports in tests
        from .gff2genestat import process_single_gff
    species = source["species"]
    loci = {}
    for gene in mapping.genes:
        loci.setdefault(mapping.locus_by_id.get(gene.gene_id, gene.gene_id), []).append(gene)
    selected, choices, requested = [], {}, []
    for locus, genes in loci.items():
        row = representative_map.choice(locus, species)
        if len({gene.seqid for gene in genes}) > 1 or len({gene.strand for gene in genes} - {"."}) > 1:
            raise ValueError(f"Isoforms of locus {locus} have incompatible coordinates")
        transcript = matching_identifier(row["source_transcript_id"], {gene.gene_id for gene in genes}, species)
        if transcript is None:
            # Effective FASTAs are keyed by locus rather than source transcript.
            locus_ids = [gene.gene_id for gene in genes
                         if identifier_aliases(gene.gene_id, species)
                         & identifier_aliases(row["gene_id"], species)]
            if len(locus_ids) != 1:
                raise ValueError(f"Representative transcript absent from FASTA for {locus}: "
                                 + row["source_transcript_id"])
            transcript = locus_ids[0]
        gene = next(gene for gene in genes if gene.gene_id == transcript)
        selected.append(gene)
        canonical = species + "_" + row["gene_id"].removeprefix(species + "_")
        requested.append(canonical)
        choices[gene.gene_id] = canonical
    path = Path(source["gff"])
    columns = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    out_columns = ["gene_id", "feature_size", "num_intron", "intron_positions", "chromosome",
                   "start", "end", "strand", "feature_blocks", "feature_type", 'cds_first_phase',
                   'phase_status', 'splice_mode']
    traits = process_single_gff(path.name, str(path.parent), requested, "CDS", "longest", columns,
                                out_columns, "report", "report", representative_map=representative_map)
    by_id = {row.gene_id: row for row in traits.itertuples(index=False)}
    normalized = []
    for gene in selected:
        row = by_id[choices[gene.gene_id]]
        if coding_traits is not None:
            coding_traits[gene.gene_id] = row
        if not math.isfinite(float(row.start)) or not math.isfinite(float(row.end)):
            raise ValueError(f"Representative CDS cannot be placed in a BED interval: {gene.gene_id}")
        if gene.seqid != row.chromosome or (gene.strand != "." and gene.strand != row.strand):
            raise ValueError(f"Representative CDS disagrees with locus coordinates: {gene.gene_id}")
        normalized.append(replace(gene, start=int(row.start) - 1, end=int(row.end), strand=row.strand))
    return replace(mapping, genes=tuple(normalized), isoform_policy="representative_map",
                   collapsed_isoform_count=len(mapping.genes) - len(normalized))


def prepare_genome(source, directory, side, minimum_mapping_fraction, *, protein_transform=None, annotation_mapper=None,
                   representative_map=None):
    sequences = {}
    aliases = {}
    species = source["species"]
    for identifier, _header, sequence in fasta_records(Path(source["fasta"])):
        safe_token(identifier, "FASTA identifier")
        if identifier in sequences or not sequence:
            raise ValueError(f"Duplicate or empty FASTA record: {identifier}")
        alphabet = "ACGTURYSWKMBDHVN" if source["mode"] == "cds" else "ACDEFGHIKLMNPQRSTVWYBXZJUO*"
        if not sequence.isascii() or set(sequence.upper()) - set(alphabet):
            raise ValueError(f"Invalid {source['mode']} sequence: {identifier}")
        sequences[identifier] = sequence.upper()
        for alias in {identifier, identifier.removeprefix(species + "_")}:
            if alias in aliases and aliases[alias] != identifier:
                raise ValueError(f"Ambiguous species-prefix alias: {alias}")
            aliases[alias] = identifier
    if not sequences:
        raise ValueError(f"Empty FASTA: {source['fasta']}")
    mapping = (annotation_mapper(source, set(aliases)) if annotation_mapper is not None else
               annotation_to_genes(source["gff"], set(aliases), feature=source["feature"] or None,
                                   attribute=source["attribute"] or None))
    genes = tuple(replace(gene, gene_id=aliases[gene.gene_id]) for gene in mapping.genes)
    if len({gene.gene_id for gene in genes}) != len(genes):
        raise ValueError("Multiple GFF identifiers map to the same FASTA record")
    mapping = replace(mapping, genes=genes, fasta_gene_count=len(sequences),
                      locus_by_id={aliases[key]: value for key, value in mapping.locus_by_id.items()})
    if len(genes) / len(sequences) < minimum_mapping_fraction:
        raise ValueError(f"Only {len(genes)}/{len(sequences)} FASTA identifiers mapped to {source['gff']}")
    map_path = representative_map or source.get("representative_map")
    representative_map = load_representative_map(map_path)
    coding_traits = {}
    mapping = (select_declared_isoforms(mapping, representative_map, source, coding_traits=coding_traits) if representative_map is not None
               else select_isoforms(mapping, {key: len(value) for key, value in sequences.items()}, "longest"))
    admitted = None
    if source.get('admitted_protein'):
        admitted = {}
        for identifier, _header, sequence in fasta_records(Path(source['admitted_protein'])):
            if identifier in admitted or identifier not in sequences:
                raise ValueError('Duplicate or absent admitted protein identifier: ' + identifier)
            sequence = sequence.upper().removesuffix('*')
            if not sequence or not sequence.isascii() or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
                raise ValueError('Invalid admitted protein: ' + identifier)
            admitted[identifier] = sequence
        if set(admitted) - {gene.gene_id for gene in mapping.genes}:
            raise ValueError('Admitted protein does not map to a selected CDS')
    proteins = {}
    for gene in mapping.genes:
        sequence = sequences[gene.gene_id]
        if protein_transform is not None:
            sequence = protein_transform(gene, sequence, mapping)
            if sequence is None:
                continue
        elif admitted is not None:
            if gene.gene_id not in admitted:
                continue
            if source['mode'] == 'cds':
                if representative_map is not None:
                    translated = translate_selected_cds(gene.gene_id, sequence, source['genetic_code'], coding_traits)
                    if translated is not None:
                        translated = translated.removesuffix('*')
                else:
                    if len(sequence) % 3:
                        raise ValueError('Admitted CDS requires an explicit phase-aware representative map: ' + gene.gene_id)
                    translated = str(Seq(sequence).translate(table=source['genetic_code'])).removesuffix('*')
                if translated != admitted[gene.gene_id]:
                    raise ValueError('Admitted protein disagrees with phase-aware CDS translation: ' + gene.gene_id)
            sequence = admitted[gene.gene_id]
        elif source['mode'] == 'cds' and representative_map is not None:
            sequence = translate_selected_cds(gene.gene_id, sequence, source['genetic_code'], coding_traits)
            if sequence is None:
                continue
        elif source["mode"] == "cds":
            if len(sequence) % 3:
                raise ValueError(f"CDS length is not divisible by three: {gene.gene_id}")
            sequence = str(Seq(sequence).translate(table=source["genetic_code"]))
        sequence = sequence.removesuffix("*")
        if not sequence or "*" in sequence:
            raise ValueError(f"Empty translation or internal stop: {gene.gene_id}")
        proteins[gene.gene_id] = sequence
    selected = set(proteins)
    excluded = {gene.gene_id for gene in mapping.genes} - selected
    mapping = replace(mapping, genes=tuple(g for g in mapping.genes if g.gene_id in selected))
    matched = {gene.gene_id for gene in genes}
    names = {identifier: species + "_" + identifier.removeprefix(species + "_") for identifier in selected}
    normalized = tuple(replace(gene, gene_id=names[gene.gene_id]) for gene in mapping.genes)
    write_bed(normalized, directory / f"{side}.bed")
    with (directory / f"{side}.pep").open("w", encoding="utf-8") as handle:
        for gene in mapping.genes:
            handle.write(f">{names[gene.gene_id]}\n{proteins[gene.gene_id]}\n")
    write_tsv(directory / f"{side}.id_map.tsv", ("original_id", "locus_id", "jcvi_id", "status"),
              ((identifier, mapping.locus_by_id.get(identifier, ""), names.get(identifier, ""),
                "selected" if identifier in selected else "translation_excluded" if identifier in excluded else
                "isoform_excluded" if identifier in matched else "unmapped")
               for identifier in sorted(sequences, key=natural_key)))
    metadata = {**mapping.metadata(), "source": source, "fasta_sha256": digest(source["fasta"]),
                "gff_sha256": digest(source["gff"]), "unmapped_count": len(sequences) - len(genes)}
    if representative_map is not None:
        metadata["representative_map_sha256"] = digest(representative_map.path)
        metadata["representative_selection_sha256"] = digest(Path(__file__).with_name("representative_selection.py"))
    if source.get('admitted_protein'):
        metadata['admitted_protein_sha256'] = digest(source['admitted_protein'])
    return normalized, metadata


def translate_selected_cds(gene_id, sequence, genetic_code, coding_traits):
    """Verify and translate the exact selected path without changing raw CDS."""
    trait = coding_traits[gene_id]
    if trait.splice_mode == 'pseudogene':
        return None
    if trait.phase_status != 'consistent' or not math.isfinite(float(trait.cds_first_phase)):
        raise ValueError('Representative CDS has missing or conflicting phase: ' + gene_id)
    if len(sequence) != int(trait.feature_size):
        raise ValueError('Representative CDS length disagrees with its annotated path: ' + gene_id)
    coding = sequence[int(trait.cds_first_phase):]
    return str(Seq(coding[:len(coding) // 3 * 3]).translate(table=genetic_code))


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


def order_karyotype(selected, genomes, anchors, mode, simple=None, scale_mode="shared"):
    """Order chromosomes without changing chromosome or gene orientations."""
    if mode not in {"none", "target", "query", "target_length", "query_length", "both_length"}:
        raise ValueError("karyotype-sort must be none, target, query, target_length, query_length or both_length")
    if scale_mode not in {"shared", "independent"}:
        raise ValueError("karyotype-scale must be shared or independent")
    if mode.endswith("_length"):
        if simple is None:
            raise ValueError("Ribbon-length sorting requires the JCVI simple blocks")
        moving = None if mode == "both_length" else (0 if mode == "target_length" else 1)
        ordered, metadata = order_by_ribbon_length(selected, genomes, simple, moving, scale_mode)
        metadata["mode"] = mode
        return ordered, metadata
    ordered = [list(seqids) for seqids in selected]
    metadata = {"mode": mode, "method": "input_order" if mode == "none" else "dominant_anchor_partner",
                "input_order": dict(zip(("target", "query"), selected, strict=True)),
                "display_order": dict(zip(("target", "query"), ordered, strict=True)),
                "orientation_changed": False, "scale_mode": scale_mode, "chromosomes": []}
    if mode == "none":
        return ordered, metadata
    moving = 0 if mode == "target" else 1
    fixed = 1 - moving
    metadata.update(fixed_side=("target", "query")[fixed], weight="unique_anchor_pairs")
    fixed_order = {seqid: index for index, seqid in enumerate(selected[fixed])}
    lookup = [{gene.gene_id: gene for gene in genome} for genome in genomes]
    fixed_ranks = {}
    for seqid in selected[fixed]:
        genes = sorted((gene for gene in genomes[fixed] if gene.seqid == seqid),
                       key=lambda gene: (gene.start, gene.end, natural_key(gene.gene_id)))
        fixed_ranks.update((gene.gene_id, rank) for rank, gene in enumerate(genes))
    counts = {seqid: Counter() for seqid in selected[moving]}
    rank_sums = {seqid: Counter() for seqid in selected[moving]}
    seen = set()
    with Path(anchors).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2 or any(fields[index] not in lookup[index] for index in (0, 1)):
                raise ValueError(f"Invalid JCVI anchor for chromosome ordering: {line.rstrip()}")
            pair = tuple(fields[:2])
            if pair in seen:
                continue
            seen.add(pair)
            moving_seqid = lookup[moving][pair[moving]].seqid
            fixed_seqid = lookup[fixed][pair[fixed]].seqid
            if moving_seqid not in counts or fixed_seqid not in fixed_order:
                continue
            counts[moving_seqid][fixed_seqid] += 1
            rank_sums[moving_seqid][fixed_seqid] += fixed_ranks[pair[fixed]]
    sort_keys = {}
    for index, seqid in enumerate(selected[moving]):
        support = counts[seqid]
        partner = min(support, key=lambda other: (-support[other], fixed_order[other])) if support else None
        count = support[partner] if partner is not None else 0
        center = Fraction(rank_sums[seqid][partner], count) if count else None
        sort_keys[seqid] = (fixed_order[partner], center, index) if count else (len(fixed_order), 0, index)
        metadata["chromosomes"].append({"seqid": seqid, "dominant_partner": partner,
                                        "dominant_anchor_count": count, "anchor_count": sum(support.values()),
                                        "partner_gene_rank_mean": float(center) if center is not None else None})
    ordered[moving] = sorted(ordered[moving], key=sort_keys.__getitem__)
    metadata["display_order"] = dict(zip(("target", "query"), ordered, strict=True))
    return ordered, metadata


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
        selected, ordering = order_karyotype(selected, genomes, directory / "target.query.lifted.anchors",
                                            plan.get("karyotype_sort", "both_length"), directory / "target.query.screened.simple",
                                            plan.get("karyotype_scale", "shared"))
        write_json(directory / "karyotype_order.json", ordering)
        dotplot = prepare_dotplot(directory, pair, genomes, selected, plan.get("dotplot_min_length", 1_000_000),
                                  plan.get("dotplot_sort", "homoeolog"))
        write_json(directory / "dotplot_order.json", dotplot)
        target_by_id = {gene.gene_id: gene.seqid for gene in genomes[0]}
        color_map = chromosome_colors(selected, genomes, directory / "target.query.lifted.anchors",
                                      plan.get("karyotype_color", "chromosome"))
        write_json(directory / "karyotype_colors.json", color_map)
        colors = color_map["chromosomes"]["target"]
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
            if plan.get("dotplot_color") == "ds":
                try:
                    from pairwise_synteny_ds import render_ds_dotplot
                except ImportError:
                    from .pairwise_synteny_ds import render_ds_dotplot
                render_ds_dotplot(directory, phase_root(plan, "ds") / pair["analysis_id"], pair, fmt, plan["ds_color_max"], filtered=True)
            else:
                try:
                    from pairwise_synteny_dotplot import render_orientation_pdf
                except ImportError:
                    from .pairwise_synteny_dotplot import render_orientation_pdf
                render_orientation_pdf(directory, pair, fmt=fmt)
            run_command([sys.executable, str(Path(__file__).with_name("pairwise_synteny_karyotype.py").resolve()),
                         "--directory", str(directory.resolve()), "--target-species", pair["target_species"],
                         "--query-species", pair["query_species"], "--format", fmt,
                         "--scale", plan.get("karyotype_scale", "shared"),
                         "--track-order", plan.get("karyotype_track_order", "target-query"),
                         "--analysis", str(source.resolve())], directory, commands, f"karyotype-{fmt}")
            for name in ("dotplot", "karyotype"):
                if not (directory / f"{name}.{fmt}").is_file() or not (directory / f"{name}.{fmt}").stat().st_size:
                    raise ValueError(f"Missing or empty {name}.{fmt}")
        shutil.rmtree(directory / ".mplconfig", ignore_errors=True)


def validate_dotplot_lengths(plan):
    """Reject missing lengths or an empty display before expensive dS alignment."""
    from kffractbias.io import read_bed

    try:
        from pairwise_synteny_ds import anchor_pairs
    except ImportError:
        from .pairwise_synteny_ds import anchor_pairs

    for pair in plan["pairs"]:
        minimum = plan.get("dotplot_min_length", 1_000_000)
        all_ids, selected_ids = [], []
        for side in ("target", "query"):
            source = pair.get("dotplot_lengths", {}).get(side, {"kind": "gff", "path": pair[side]["gff"]})
            genes = read_bed(phase_root(plan, "analysis") / pair["analysis_id"] / f"{side}.bed")
            lengths, _metadata = chromosome_lengths(source, genes, minimum)
            all_ids.append({gene.gene_id for gene in genes})
            selected_ids.append({gene.gene_id for gene in genes if not minimum or lengths[gene.seqid] >= minimum})
        pairs = anchor_pairs(phase_root(plan, "analysis") / pair["analysis_id"] / "target.query.lifted.anchors")
        if any(genes[i] not in all_ids[i] for genes in pairs for i in (0, 1)):
            raise ValueError("Invalid dotplot anchor")
        if not any(all(genes[i] in selected_ids[i] for i in (0, 1)) for genes in pairs):
            raise ValueError(f"No anchors connect chromosomes meeting the dotplot minimum length ({minimum} bp)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    planning = subparsers.add_parser("plan")
    planning.add_argument("--workspace", required=True, type=Path)
    planning.add_argument("--pairs", required=True, type=Path)
    planning.add_argument("--sequence-mode", choices=("auto", "protein", "cds"), default="auto")
    planning.add_argument("--genetic-code", type=int, default=1)
    planning.add_argument("--representative-map", default="", type=Path)
    planning.add_argument("--representative-inputs", default="", type=Path)
    planning.add_argument("--cscore", type=float, default=0.7)
    planning.add_argument("--min-anchors", type=int, default=4)
    planning.add_argument("--distance", type=int, default=20)
    planning.add_argument("--minimum-mapping-fraction", type=float, default=1)
    planning.add_argument("--formats", default="pdf,svg,png")
    planning.add_argument("--dotplot-color", choices=("orientation", "ds"), default="orientation")
    planning.add_argument("--dotplot-min-length", type=int, default=1_000_000, help="Minimum assembly chromosome length in bp; 0 disables filtering")
    planning.add_argument("--dotplot-sort", choices=("homoeolog", "karyotype", "none"), default="homoeolog",
                          help="Adjacent supported 2x2 chromosome groups, ribbon order, or original BED order")
    planning.add_argument("--ds-color-max", type=float, default=2.0, help="Upper display color limit only; anchors are never filtered")
    planning.add_argument("--karyotype-sort", choices=("none", "target", "query", "target_length", "query_length", "both_length"), default="both_length",
                          help="Minimize width-weighted ribbon length (default both_length); target/query use dominant-anchor partners")
    planning.add_argument("--karyotype-color", choices=("chromosome", "homoeolog"), default="chromosome",
                          help="Soft chromosome colors (default); homoeolog shares supported 2x2 group colors")
    planning.add_argument("--karyotype-track-order", choices=("target-query", "query-target"), default="target-query",
                          help="Top-to-bottom ribbon display order, independent of analysis target/query")
    planning.add_argument("--karyotype-scale", choices=("shared", "independent"), default="shared",
                          help="Shared gene width and one scale bar (default); independent normalizes each track separately")
    planning.add_argument("--outfile", required=True, type=Path)
    contract = subparsers.add_parser("contract")
    contract.add_argument("--plan", required=True, type=Path)
    contract.add_argument("--phase", choices=("analysis", "ds", "plots"), required=True)
    for phase in ("analysis", "ds", "plots"):
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
            elif args.command == "ds":
                try:
                    from pairwise_synteny_ds import estimate_ds
                except ImportError:
                    from .pairwise_synteny_ds import estimate_ds
                validate_dotplot_lengths(plan)
                estimate_ds(plan, args.output, args.cpus)
            else:
                render(plan, args.output)
    except (ValueError, OSError, RuntimeError, csv.Error, subprocess.SubprocessError) as exc:
        print(f"Pairwise synteny failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
