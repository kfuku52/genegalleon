#!/usr/bin/env python3
"""Native WGD candidate inference and conservative gene-node origin evidence."""

import argparse
import csv
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, nullcontext
from functools import partial
from pathlib import Path

try:
    from pairwise_synteny import digest, prepare_genome, safe_token, source_file, write_json
    from pairwise_synteny_ds import DS_STATUSES, PAIR_COLUMNS, align_pair, ds_tool_identity, load_cds
    from parallel_io import ordered_bounded_results
    from wgd_alignment_reuse import alignment_identity, matching_cached_rows, reuse_bundle
    from wgd_cross_ks import CrossKs, cross_ks_cds_sources, plot_cross_ks
    from wgd_cross_overlap import OVERLAP_FIELDS, duplication_synteny_overlap
    from wgd_cross_phase import (
        ORIENTATION_PHASE_FIELDS,
        PHASE_FIELDS,
        orientation_phase_overlap,
        phase_rows,
        require_paired_orientations,
        trace_has_unique_phase,
        unique_phase_windows,
    )
    from wgd_cross_regions import NESTED_FIELDS, REGION_DETAIL, REGION_SUMMARY, nested_copy_overlap, summarize_regions
    from wgd_cross_replication import REPLICATION_FIELDS, replicated_copy_rows
    from wgd_cross_support import (
        DUPLICATION_SUPPORT_FIELDS,
        SUPPORT_FIELDS,
        DuplicationTreeSupport,
    )
    from wgd_cross_synteny import (
        COVERAGE_OUTPUT,
        NONREUSE_OUTPUT,
        ORIENTATION_OUTPUT,
        PAIR_OUTPUT,
        WINDOW_OUTPUT,
        CrossSynteny,
        archive_window_floor,
        cross_synteny_inputs,
        orientation_rows,
        read_pairs,
    )
    from wgd_evidence import KsBranchPlacer, combine_node, number, read_table, summarize_events, valid_position
    from wgd_mul_diagnostics import run_diagnostics
    from wgd_rooted_copy_ks import RootedCopyKs
    from wgd_taxon_holdout import HOLDOUT_FIELDS
    from wgd_terminal_ks import TerminalKs
    from wgd_tree_scan import COVERAGE_FIELDS as DUPLICATION_COVERAGE_FIELDS
    from wgd_tree_scan import COPY_SPLIT_FIELDS, PANEL_FIELDS, PANEL_SUMMARY_FIELDS, ROOT_FIELDS, FamilyPanels, scan_duplications
except ImportError:
    from .pairwise_synteny import digest, prepare_genome, safe_token, source_file, write_json
    from .pairwise_synteny_ds import DS_STATUSES, PAIR_COLUMNS, align_pair, ds_tool_identity, load_cds
    from .parallel_io import ordered_bounded_results
    from .wgd_alignment_reuse import alignment_identity, matching_cached_rows, reuse_bundle
    from .wgd_cross_ks import CrossKs, cross_ks_cds_sources, plot_cross_ks
    from .wgd_cross_overlap import OVERLAP_FIELDS, duplication_synteny_overlap
    from .wgd_cross_phase import (
        ORIENTATION_PHASE_FIELDS,
        PHASE_FIELDS,
        orientation_phase_overlap,
        phase_rows,
        require_paired_orientations,
        trace_has_unique_phase,
        unique_phase_windows,
    )
    from .wgd_cross_regions import NESTED_FIELDS, REGION_DETAIL, REGION_SUMMARY, nested_copy_overlap, summarize_regions
    from .wgd_cross_replication import REPLICATION_FIELDS, replicated_copy_rows
    from .wgd_cross_support import (
        DUPLICATION_SUPPORT_FIELDS,
        SUPPORT_FIELDS,
        DuplicationTreeSupport,
    )
    from .wgd_cross_synteny import (
        COVERAGE_OUTPUT,
        NONREUSE_OUTPUT,
        ORIENTATION_OUTPUT,
        PAIR_OUTPUT,
        WINDOW_OUTPUT,
        CrossSynteny,
        archive_window_floor,
        cross_synteny_inputs,
        orientation_rows,
        read_pairs,
    )
    from .wgd_evidence import KsBranchPlacer, combine_node, number, read_table, summarize_events, valid_position
    from .wgd_mul_diagnostics import run_diagnostics
    from .wgd_rooted_copy_ks import RootedCopyKs
    from .wgd_taxon_holdout import HOLDOUT_FIELDS
    from .wgd_terminal_ks import TerminalKs
    from .wgd_tree_scan import COVERAGE_FIELDS as DUPLICATION_COVERAGE_FIELDS
    from .wgd_tree_scan import COPY_SPLIT_FIELDS, PANEL_FIELDS, PANEL_SUMMARY_FIELDS, ROOT_FIELDS, FamilyPanels, scan_duplications


CROSS_DS_METHODS = ("yn00", "gy94", "gy94-regularized")
CROSS_GY94_STATUSES = frozenset(("identical", "kappa_boundary", "rate_boundary", "not_converged",
                                "frequency_not_converged"))
CROSS_DS_REPORT_LABELS = {"gy94": "GY94_CF3x4_pairwise", "gy94-regularized": "GY94_CF3x4_log_ratio_penalty"}


def write_table(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def species_tree(path):
    from nwkit.rooting_state import require_rooted
    from nwkit.util import read_tree
    from nwkit.wgd_count import _count_tree

    tree = read_tree(str(path), "auto", True, rooted="auto")
    require_rooted(tree, "Native WGD evidence requires a rooted species tree.")
    _count_tree(tree)
    return tree


def owned_identity():
    import kffractbias.io
    import kffractbias.selfevidence
    import kffractbias.selfscan
    import nwkit
    import nwkit.ksrate
    import nwkit.ksrate_cli
    import nwkit.ksrate_model
    import nwkit.wgd_count
    import nwkit.wgd_count_cli
    import nwkit.wgd_count_fit
    import nwkit.wgd_count_model
    import nwkit.wgd_tree
    import nwkit.wgd_tree_cli
    import nwkit.wgd_tree_model

    modules = (kffractbias.io, kffractbias.selfevidence, kffractbias.selfscan, nwkit.ksrate,
               nwkit.ksrate_cli, nwkit.ksrate_model, nwkit.wgd_count, nwkit.wgd_count_cli,
               nwkit.wgd_count_fit, nwkit.wgd_count_model, nwkit.wgd_tree, nwkit.wgd_tree_cli, nwkit.wgd_tree_model)
    files = [Path(module.__file__).resolve() for module in modules]
    # Shared parsers/tree utilities also affect branch IDs and scientific fits.
    for root in (Path(nwkit.__file__).parent, Path(kffractbias.io.__file__).parent):
        files.extend(root.rglob("*.py"))
    files += [Path(__file__).resolve(), Path(__file__).with_name("wgd_evidence.py").resolve(),
              Path(__file__).with_name("wgd_alignment_reuse.py").resolve(),
              Path(__file__).with_name("wgd_mul_diagnostics.py").resolve(),
              Path(__file__).with_name("wgd_tree_scan.py").resolve(),
              Path(__file__).with_name("wgd_taxon_holdout.py").resolve(),
              Path(__file__).with_name("wgd_cross_ks.py").resolve(),
              Path(__file__).with_name("wgd_cross_overlap.py").resolve(),
              Path(__file__).with_name("wgd_cross_phase.py").resolve(),
              Path(__file__).with_name("wgd_cross_regions.py").resolve(),
              Path(__file__).with_name("wgd_cross_replication.py").resolve(),
              Path(__file__).with_name("wgd_cross_support.py").resolve(),
              Path(__file__).with_name("wgd_cross_synteny.py").resolve(),
              Path(__file__).with_name("wgd_ks_outgroups.py").resolve(),
              Path(__file__).with_name("wgd_quartet_ks.py").resolve(),
              Path(__file__).with_name("wgd_rooted_copy_ks.py").resolve(),
              Path(__file__).with_name("wgd_terminal_ks.py").resolve(),
              Path(__file__).with_name("pairwise_synteny.py").resolve()]
    ds = ds_tool_identity()
    return {"nwkit": nwkit.__version__, "ds": ds,
            "source_hashes": {str(path): digest(path) for path in files},
            "last_version": subprocess.check_output(["lastal", "--version"], text=True).strip()}


def make_plan(args):
    workspace = args.workspace.resolve()

    def resolve(value, default):
        path = Path(value) if value else workspace / default
        return path.resolve() if path.is_absolute() else (workspace / path).resolve()

    tree_path = resolve(args.species_tree, "output/species_tree/species_tree_summary/dated_species_tree.nwk")
    if not args.species_tree and not tree_path.is_file():
        tree_path = workspace / "output/species_tree/species_tree_summary/undated_species_tree.nwk"
    tree = species_tree(tree_path)
    names = sorted(leaf.name for leaf in tree.leaves())
    counts = resolve(args.counts, "output/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv")
    members = resolve(args.members, "output/orthofinder/Orthogroups/Orthogroups.tsv")
    inputs = {str(path): digest(path) for path in (tree_path, counts, members)}
    tree_archive = resolve(args.gene_trees, "") if getattr(args, "gene_trees", "") else None
    if tree_archive:
        inputs[str(tree_archive)] = digest(tree_archive)
    cross_manifest = getattr(args, "cross_synteny", "")
    cross_pairs = None
    if cross_manifest:
        if tree_archive is None:
            raise ValueError("Cross synteny requires --gene-trees for exact copy-level tree agreement")
        cross_path = resolve(cross_manifest, "")
        inputs[str(cross_path)] = digest(cross_path)

        def hash_input(path):
            if str(path) not in inputs:
                inputs[str(path)] = digest(path)

        cross_pairs = read_pairs(cross_path, lambda path: resolve(path, ""), set(names), hash_input)
    support_archive = ""
    if getattr(args, "cross_tree_support", ""):
        if tree_archive is None:
            raise ValueError("Tree split support requires --gene-trees")
        support_archive = str(resolve(args.cross_tree_support, ""))
        inputs[support_archive] = digest(support_archive)
    genomes_path = resolve(args.genomes, "input/wgd_genomes.tsv")
    if args.genomes or genomes_path.is_file():
        inputs[str(genomes_path)] = digest(genomes_path)
        genomes = read_table(genomes_path)
        if genomes and ("species" not in genomes[0] or set(genomes[0]) - {
                "species", "mode", "fasta", "gff", "cds", "feature", "attribute"}):
            raise ValueError("Invalid genome manifest columns")
    else:
        genomes = []
    codes = {}
    code_file = workspace / "input/species_genetic_code/species_genetic_code.tsv"
    if code_file.is_file():
        inputs[str(code_file)] = digest(code_file)
        for row in read_table(code_file):
            if row["species"] in codes:
                raise ValueError("Duplicate species genetic code")
            codes[row["species"]] = int(row["genetic_code"])
    cross_ks_manifest = getattr(args, "cross_ks_cds", "")
    cross_ks_sources = {}
    if cross_ks_manifest:
        if cross_pairs is None:
            raise ValueError("Cross-block Ks requires --cross-synteny")
        path = resolve(cross_ks_manifest, "")
        inputs[str(path)] = digest(path)
        cross_ks_sources = cross_ks_cds_sources(path, workspace, set(names), codes, args.genetic_code)
        for source in cross_ks_sources.values():
            inputs[source["cds"]] = digest(source["cds"])
    seen = set()
    for genome in genomes:
        name = safe_token(genome["species"], "species")
        if name not in names or name in seen:
            raise ValueError("Genome manifest requires unique species-tree tips")
        seen.add(name)
        mode = genome.get("mode") or args.sequence_mode
        if mode not in {"cds", "protein"}:
            raise ValueError("Genome sequence mode must be cds or protein")
        genome.update(mode=mode, genetic_code=codes.get(name, args.genetic_code),
                      feature=genome.get("feature", ""), attribute=genome.get("attribute", ""))
        for field, directory, suffixes in (("fasta", "species_cds" if mode == "cds" else "species_protein",
                                            (".fa", ".fas", ".fasta", ".fna", ".faa")),
                                           ("gff", "species_gff", (".gff", ".gff3"))):
            genome[field] = source_file(workspace, genome.get(field, ""), directory, name, suffixes)
        cds = genome.get("cds", "")
        if mode == "cds" and cds and source_file(workspace, cds, "species_cds", name, (".fa", ".fas", ".fasta", ".fna")) != genome["fasta"]:
            raise ValueError("CDS mode requires cds and fasta to refer to the same source")
        genome["cds"] = genome["fasta"] if mode == "cds" else (
            source_file(workspace, cds, "species_cds", name, (".fa", ".fas", ".fasta", ".fna")) if cds else "")
        for field in ("fasta", "gff", "cds"):
            if genome[field]:
                inputs[genome[field]] = digest(genome[field])
    parameters = {key: getattr(args, key) for key in (
        "count_bootstrap", "ks_bootstrap", "seed", "max_pairs", "max_ks_families", "max_count_families",
        "diagonal_bound", "cscore", "self_hit_percent", "minimum_mapping_fraction", "alpha",
        "min_coverage", "min_blocks", "max_states", "max_iterations", "multiplicity")}
    parameters["min_duplication_score"] = getattr(args, "min_duplication_score", 0.01)
    parameters["duplication_taxon_holdout"] = getattr(args, "duplication_taxon_holdout", 0)
    if parameters["duplication_taxon_holdout"] not in (0, 1):
        raise ValueError("Duplication taxon holdout must be 0 or 1")
    if parameters["duplication_taxon_holdout"] and tree_archive is None:
        raise ValueError("Duplication taxon holdout requires --gene-trees")
    for audit in ("duplication_coverage_audit", "duplication_root_audit", "duplication_copy_split_audit"):
        if getattr(args, audit, 0):
            if tree_archive is None:
                raise ValueError(f"{audit.replace('_', ' ').capitalize()} requires --gene-trees")
            parameters[audit] = 1
    panel_count = getattr(args, "duplication_family_panels", 0)
    if panel_count and (not 2 <= panel_count <= 64 or tree_archive is None):
        raise ValueError("Family panels require --gene-trees and a panel count of 2 to 64, or 0 to disable")
    if panel_count:
        parameters["duplication_family_panels"] = panel_count
    parameters["min_cross_windows"] = getattr(args, "min_cross_windows", 3)
    parameters["cross_orientation_audit"] = getattr(args, "cross_orientation_audit", 0)
    if parameters["cross_orientation_audit"] not in (0, 1):
        raise ValueError("Cross-synteny orientation audit must be 0 or 1")
    if parameters["cross_orientation_audit"] and cross_pairs is None:
        raise ValueError("Cross-synteny orientation audit requires --cross-synteny")
    if getattr(args, "cross_require_unique_phase", 0):
        if cross_pairs is None:
            raise ValueError("Unique cross-synteny phase requires --cross-synteny")
        parameters["cross_require_unique_phase"] = 1
    if getattr(args, "cross_require_both_orientations", 0):
        if not (parameters.get("cross_require_unique_phase") and parameters["cross_orientation_audit"]):
            raise ValueError("Shared cross-synteny paths require --cross-require-unique-phase 1 and --cross-orientation-audit 1")
        parameters["cross_require_both_orientations"] = 1
    parameters["count_bootstrap_workers"] = getattr(args, "count_bootstrap_workers", 1)
    parameters["count_scan_workers"] = getattr(args, "count_scan_workers", 1)
    for key in ("count_event_fractions", "count_family_gamma_shape"):
        value = getattr(args, key, "")
        if value:
            if tree_archive is not None:
                raise ValueError("Count event fractions and family gamma shape require gene-count mode")
            parameters[key] = value
    count_inputs = {}
    for kind in ("detection", "rate_groups"):
        value = getattr(args, f"count_{kind}_tsv", "")
        if value:
            if tree_archive is not None:
                raise ValueError("Count detection and rate-group tables require gene-count mode")
            path = resolve(value, "")
            inputs[str(path)] = digest(path)
            count_inputs[f"count_{kind}_tsv"] = str(path)
    checkpoint = getattr(args, "count_bootstrap_checkpoint", "")
    if checkpoint:
        if not args.count_bootstrap or tree_archive is not None:
            raise ValueError("A bootstrap checkpoint requires gene-count mode and positive --count-bootstrap")
        checkpoint = str(resolve(checkpoint, ""))
    parameters["max_ks_quartets"] = getattr(args, "max_ks_quartets", 0)
    if getattr(args, "ks_rooted_copies", 0):
        if parameters["max_ks_quartets"] < 1 or not cross_ks_sources:
            raise ValueError("Rooted-copy dS requires --cross-ks-cds and positive --max-ks-quartets")
        parameters["ks_rooted_copies"] = 1
    outgroup_limit = getattr(args, "max_ks_outgroup_families", 0)
    if outgroup_limit < 0 or (outgroup_limit and not cross_ks_sources):
        raise ValueError("Common-outgroup dS requires --cross-ks-cds and a nonnegative family cap")
    if outgroup_limit:
        parameters["max_ks_outgroup_families"] = outgroup_limit
    if getattr(args, "ks_sister_panels", 0):
        if not outgroup_limit:
            raise ValueError("Sister dS panels require a positive --max-ks-outgroup-families")
        parameters["ks_sister_panels"] = 1
    terminal_limit = getattr(args, "max_terminal_ks_families", 0)
    if terminal_limit < 0 or (terminal_limit and not cross_ks_sources):
        raise ValueError("Terminal dS requires --cross-ks-cds and a nonnegative family cap")
    if terminal_limit:
        parameters["max_terminal_ks_families"] = terminal_limit
    if getattr(args, "ks_component_holdout", 0):
        if not cross_ks_sources:
            raise ValueError("Copy dS component holdout requires --cross-ks-cds")
        parameters["ks_component_holdout"] = 1
    ds_method = getattr(args, "cross_ks_method", "yn00")
    if ds_method != "yn00":
        if ds_method not in CROSS_DS_METHODS or not cross_ks_sources:
            raise ValueError("Alternative dS estimation requires --cross-ks-cds and a supported --cross-ks-method")
        parameters["cross_ks_method"] = ds_method
        import cdskit.cli
        import cdskit.dnds_mle
        import scipy

        modules = (cdskit.cli, cdskit.dnds_mle, scipy)
        inputs.update({str(Path(module.__file__).resolve()): digest(Path(module.__file__)) for module in modules})
    if parameters["max_ks_quartets"] and not cross_ks_sources:
        raise ValueError("Four-copy dS diagnostics require --cross-ks-cds")
    if (args.count_bootstrap < 0 or args.ks_bootstrap < 0 or args.seed < 0
            or parameters["max_ks_quartets"] < 0
            or any(parameters[key] < 1 for key in ("max_pairs", "max_ks_families", "max_count_families",
                                                  "diagonal_bound", "min_blocks", "max_iterations",
                                                  "count_bootstrap_workers", "count_scan_workers"))
            or not 0 < args.cscore <= 1 or not 0 < args.minimum_mapping_fraction <= 1
            or not 0 < args.self_hit_percent <= 100 or not 0 < args.alpha < 1
            or not 0 < parameters["min_duplication_score"] <= 1
            or parameters["min_cross_windows"] < 1
            or not 0 < args.min_coverage <= 1 or args.max_states < 16 or args.multiplicity < 2):
        raise ValueError("Invalid WGD analysis parameters")
    tools = owned_identity()
    inputs.update(tools["source_hashes"])
    inputs.update(tools["ds"]["source_hashes"])
    reuse = getattr(args, "cross_ks_reuse_summary", ())
    if isinstance(reuse, str):
        reuse = [reuse]
    reuse = [summary for summary in reuse if summary]
    if reuse and not cross_ks_sources:
        raise ValueError("Cross-block alignment reuse requires --cross-ks-cds")
    alignment, files = alignment_identity(tools["ds"]) if cross_ks_sources else (None, ())
    inputs.update(files)
    reused_alignments = {}
    if reuse:
        reused_alignments, hashes = reuse_bundle((resolve(summary, "") for summary in reuse), alignment)
        managed_output = workspace / "output/genome_evolution/wgd_ssd"
        protect_output(managed_output, hashes)
        inputs.update(hashes)
    if checkpoint:
        protect_checkpoint(checkpoint, inputs | {str(args.outfile.resolve()): ""}, workspace)
    protect_output(args.outfile, inputs, workspace)
    return {"schema_version": 1, "workspace": str(workspace), "species_tree": str(tree_path),
            "counts": str(counts), "members": str(members), "gene_trees": str(tree_archive) if tree_archive else "",
            "cross_pairs": cross_pairs, "cross_ks_sources": cross_ks_sources, "cross_tree_support": support_archive,
            "alignment_identity": alignment, "reused_cross_ks_alignments": reused_alignments,
            "count_bootstrap_checkpoint": checkpoint, **count_inputs,
            "genomes": genomes, "species": names,
            "parameters": parameters, "tools": tools, "input_hashes": inputs,
            "absent_inputs": ([str(genomes_path)] if not genomes_path.is_file() and not args.genomes else [])
            + ([str(code_file)] if not code_file.is_file() else [])}


def verify(plan):
    identity = plan.get("alignment_identity")
    if identity is not None and alignment_identity(plan["tools"]["ds"])[0] != identity:
        raise ValueError("MAFFT or backtranslation identity changed during WGD analysis")
    for path, expected in plan["input_hashes"].items():
        if digest(path) != expected:
            raise ValueError(f"WGD input changed during analysis: {path}")
    if any(Path(path).exists() for path in plan.get("absent_inputs", [])):
        raise ValueError("A previously absent WGD input appeared during analysis")


def protect_output(path, inputs, workspace=None):
    target = Path(path).resolve()
    sources = [Path(source).resolve() for source in inputs]
    if any(source == target or target in source.parents or (source.is_dir() and source in target.parents)
           for source in sources):
        raise ValueError("WGD output must not replace or contain an input")
    if workspace and (target == Path(workspace) / "input" or Path(workspace) / "input" in target.parents):
        raise ValueError("WGD outputs must not be written into curated workspace/input")


def protect_checkpoint(path, inputs, workspace):
    protected = [workspace / "output/genome_evolution/wgd_ssd", workspace / "output/tmp/wgd_ssd",
                 workspace / "output/artifact_provenance", workspace / "output/.gg_global_artifacts"]
    checkpoint = Path(path).resolve()
    if not checkpoint.parent.is_dir():
        raise ValueError("WGD bootstrap checkpoint needs an existing parent directory")
    if any(parent == checkpoint or parent in checkpoint.parents for parent in protected):
        raise ValueError("WGD bootstrap checkpoint must be outside managed WGD outputs and task scratch")
    protect_output(checkpoint, inputs, workspace)


def validated_members(plan):
    rows, key, _ = _validated_family_tables(plan)
    return rows, key


def _validated_family_tables(plan):
    rows = read_table(plan["members"])
    count_rows = dict(validated_counts(plan))
    if not rows:
        raise ValueError("Empty family-membership or gene-count table")
    key = "family_id" if "family_id" in rows[0] else "Orthogroup"
    if set(rows[0]) != {key, *plan["species"]}:
        raise ValueError("Family membership species must exactly match the full species tree")
    prefixes = [(name, name + "_") for name in plan["species"]]
    seen_families, seen_genes = set(), set()
    for row in rows:
        family = safe_token(row[key], "family ID")
        if family in seen_families or family not in count_rows:
            raise ValueError("Duplicate or unmatched family-membership row")
        seen_families.add(family)
        for name, prefix in prefixes:
            raw_genes = row[name].strip()
            genes = [gene.strip() for gene in raw_genes.split(",")] if raw_genes else []
            if "" in genes:
                raise ValueError("Empty gene identifier in family membership")
            for gene in genes:
                safe_token(gene, "Gene identifier")
                if not gene.startswith(prefix) or gene == prefix:
                    raise ValueError("Family-membership gene identifier must have its species prefix")
            if (len(genes) > 1 and len(genes) != len(set(genes))) or not seen_genes.isdisjoint(genes):
                raise ValueError("Gene identifiers must be unique across family membership")
            seen_genes.update(genes)
            count = count_rows[family][name]
            if count is not None and count != len(genes):
                raise ValueError("Family membership and gene count disagree")
    if seen_families != set(count_rows):
        raise ValueError("Family membership and gene counts have different families")
    return rows, key, count_rows


def validate_annotation_counts(plan, members):
    from kffractbias.io import annotation_to_genes

    for source in plan["genomes"]:
        name = source["species"]
        identifiers = {gene.strip() for row in members for gene in row[name].split(",") if gene.strip()}
        if not identifiers:
            continue
        aliases = {}
        for identifier in identifiers:
            for alias in (identifier, identifier.removeprefix(name + "_")):
                if alias in aliases and aliases[alias] != identifier:
                    raise ValueError("Ambiguous species-prefix alias in counted genes")
                aliases[alias] = identifier
        mapping = annotation_to_genes(source["gff"], set(aliases), feature=source["feature"] or None,
                                      attribute=source["attribute"] or None)
        seen_loci = {}
        for identifier, locus in mapping.locus_by_id.items():
            gene = aliases[identifier]
            if locus in seen_loci and seen_loci[locus] != gene:
                raise ValueError(f"Family gene counts include multiple isoforms of annotated locus {name}:{locus}; "
                                 "provide representative-locus membership and counts")
            seen_loci[locus] = gene


def contract(plan):
    root = Path(plan["workspace"])
    args = ["--manifest", str(root / "output/artifact_provenance/genome_evolution/wgd_ssd.json"),
            "--step", "wgd_ssd", "--family-id", "all_families", "--logical-root",
            str(root / "output/.gg_global_artifacts"), "--workspace-root", str(root),
            "--output", f"results={root / 'output/genome_evolution/wgd_ssd'}"]
    for index, path in enumerate(sorted(plan["input_hashes"])):
        args += ["--input", f"source{index}={path}"]
    return args + ["--parameter", "configuration=" + json.dumps(plan, sort_keys=True)]


def run_tool(argv, directory, label):
    logs = directory / "logs"
    logs.mkdir(exist_ok=True)
    with (logs / f"{label}.log").open("w", encoding="utf-8") as handle:
        process = subprocess.run([str(value) for value in argv], cwd=directory, stdout=handle,
                                 stderr=subprocess.STDOUT, env={**os.environ, "MPLBACKEND": "Agg"})
    if process.returncode:
        log = logs / (label + ".log")
        raise RuntimeError(f"{label} failed ({process.returncode}); see {log}\n{log.read_text(errors='replace')[-3000:]}")


def validated_counts(plan):
    names = plan["species"]
    rows = read_table(plan["counts"])
    if not rows:
        raise ValueError("Empty gene-count table")
    key = "family_id" if "family_id" in rows[0] else "Orthogroup"
    if key not in rows[0] or set(rows[0]) - {key, "Total", *names} or set(names) - set(rows[0]):
        raise ValueError("Gene-count species must exactly match the full species tree")
    result, seen = [], set()
    for row in rows:
        family = safe_token(row[key], "family ID")
        if family in seen:
            raise ValueError("Duplicate gene-count family")
        seen.add(family)
        values = {name: number(row[name]) for name in names}
        if any(value is not None and value != int(value) for value in values.values()):
            raise ValueError("Gene counts must be integers")
        if any(row[name] not in {"NA", "nan", "", "."} and values[name] is None for name in names):
            raise ValueError("Invalid gene count")
        if "Total" in row and row["Total"] not in {"NA", "nan", "", "."}:
            total = number(row["Total"])
            observed = sum(value for value in values.values() if value is not None)
            if (total is None or total != int(total) or total < observed
                    or (all(value is not None for value in values.values()) and total != observed)):
                raise ValueError("Gene-count Total disagrees with species counts")
        result.append((family, values))
    return result


def normalized_counts(plan, directory, tree=None, count_rows=None, minimum_families=2):
    names = plan["species"]
    tree = tree if tree is not None else species_tree(plan["species_tree"])
    clades = [{leaf.name for leaf in child.leaves()} for child in tree.children]
    selected, audit = [], []
    for family, values in validated_counts(plan) if count_rows is None else count_rows.items():
        retained = all(any(values[name] is not None and values[name] > 0 for name in group) for group in clades)
        audit.append({"family_id": family, "selection": "root_clades_present" if retained else "excluded_root_clade_absence"})
        if retained:
            selected.append({"family_id": family, **{name: "NA" if values[name] is None else int(values[name]) for name in names}})
    limit = plan["parameters"]["max_count_families"]
    if len(selected) > limit:
        selected = sorted(random.Random(plan["parameters"]["seed"]).sample(selected, limit), key=lambda row: row["family_id"])
    chosen = {row["family_id"] for row in selected}
    for row in audit:
        if row["selection"] == "root_clades_present" and row["family_id"] not in chosen:
            row["selection"] = "not_sampled"
    if len(selected) < minimum_families:
        raise ValueError("At least two root-clade-present families are required")
    write_table(directory / "counts.tsv", selected, ["family_id", *names])
    write_table(directory / "count_family_selection.tsv", audit, ["family_id", "selection"])


def raw_anchors(path, species):
    result, block = [], 0
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            block += 1
        elif line.strip():
            fields = line.split()
            if len(fields) < 2 or block == 0:
                raise ValueError("Malformed unquota anchors")
            a, b = sorted(fields[:2])
            result.append({"species": species, "block_id": f"block{block:06d}", "gene_a": a, "gene_b": b})
    return result


def prepare_sources(plan, output, cpus):
    from kffractbias.io import annotation_locus_order, read_bed

    positions, sequences, summaries, anchors = [], {}, {}, []
    parameters = plan["parameters"]
    for source in plan["genomes"]:
        name = source["species"]
        directory = output / "self_synteny" / name
        directory.mkdir(parents=True)
        genes, metadata = prepare_genome(source, directory, "self", parameters["minimum_mapping_fraction"])
        genes = read_bed(directory / "self.bed")
        loci = {row["jcvi_id"]: row["locus_id"] for row in read_table(directory / "self.id_map.tsv") if row["status"] == "selected"}
        annotation = annotation_locus_order(source["gff"], feature=metadata["feature"], attribute=metadata["attribute"])
        ranks = defaultdict(int)
        locus_ranks = {}
        for locus in annotation:
            ranks[locus.seqid] += 1
            locus_ranks[locus.gene_id] = (locus, ranks[locus.seqid])
        for gene in genes:
            locus, rank = locus_ranks[loci[gene.gene_id]]
            if locus.seqid != gene.seqid:
                raise ValueError("Selected gene and annotated locus chromosomes differ")
            positions.append({"gene_id": gene.gene_id, "species": name, "locus_id": loci[gene.gene_id],
                              "seqid": gene.seqid, "rank": rank, "start": locus.start, "end": locus.end})
        write_json(directory / "mapping.json", metadata)
        run_tool([sys.executable, "-m", "kffractbias.selfscan", "align", "--aligner=last", f"--cpus={cpus}",
                  f"--cscore={parameters['cscore']}", f"--self-hit-percent={parameters['self_hit_percent']}",
                  "--sequence-type=prot"], directory, "align")
        filtered = list(directory.glob("self.self.last*.inverse.filtered"))
        if len(filtered) != 1:
            raise RuntimeError("Expected one filtered self alignment")
        run_tool([sys.executable, "-m", "kffractbias.selfscan", "scan", filtered[0], filtered[0].with_suffix(""),
                  directory / "self.bed", directory / "unquota.anchors", f"--diagonal-bound={parameters['diagonal_bound']}",
                  "--screening=none", "--allow-empty"], directory, "scan")
        anchors.extend(raw_anchors(directory / "self.self.lifted.anchors", name))
        summaries[name] = json.loads((directory / "self.self.raw.summary.json").read_text())
        summaries[name]["num_annotated_loci"] = len(annotation)
        sequences[name] = load_cds(source["cds"], directory, "self", name, source["genetic_code"]) if source["cds"] else {}
    write_table(output / "gene_positions.tsv", positions,
                ["gene_id", "species", "locus_id", "seqid", "rank", "start", "end"])
    if len({row["gene_id"] for row in positions}) != len(positions):
        raise ValueError("Duplicate canonical gene identifiers across prepared genomes")
    return positions, sequences, summaries, anchors


def estimate_ks(plan, output, sequences, anchors, cpus, family_tables=None):
    rng = random.Random(plan["parameters"]["seed"])
    codes = {source["species"]: source["genetic_code"] for source in plan["genomes"]}
    requests, selection = {}, []
    sequence_species = sorted(sequences)
    anchor_pairs = defaultdict(set)
    for row in anchors:
        anchor_pairs[row["species"]].add(tuple(sorted((row["gene_a"], row["gene_b"]))))
    for name in sequence_species:
        all_pairs = sorted(anchor_pairs[name])
        pairs = [pair for pair in all_pairs if all(gene in sequences[name] for gene in pair)]
        limit = plan["parameters"]["max_pairs"]
        if len(pairs) > limit:
            pairs = sorted(rng.sample(pairs, limit))
        for pair in pairs:
            requests[pair] = (name, name, "", "self_anchor")
        chosen_pairs = set(pairs)
        for pair in all_pairs:
            if pair not in chosen_pairs:
                status = "cds_unavailable" if any(gene not in sequences[name] for gene in pair) else "not_sampled"
                selection.append({"pair_id": "|".join(pair), "species_a": name, "species_b": name,
                                  "family_id": "", "pair_role": "self_anchor", "status": status, "dS": "NA"})
    members, family_key, counts = family_tables or _validated_family_tables(plan)
    single_copy, audit = [], []
    for row in members:
        family = row[family_key]
        family_counts = counts[family]
        # Validation ensures commas separate nonempty identifiers.
        if any((count is not None and count > 1) or "," in row[name] for name, count in family_counts.items()):
            audit.append({"family_id": family, "selection": "known_multicopy_family"})
            continue
        available = {}
        for name in sequence_species:
            gene = row[name].strip()
            if family_counts[name] == 1 and gene in sequences[name]:
                available[name] = gene
        if len(available) >= 2:
            single_copy.append((family, available))
        audit.append({"family_id": family, "selection": "eligible" if len(available) >= 2 else "insufficient_single_copy_cds_species"})
    del members, counts
    limit = plan["parameters"]["max_ks_families"]
    if len(single_copy) > limit:
        single_copy = sorted(rng.sample(single_copy, limit))
    chosen = {family for family, _ in single_copy}
    for row in audit:
        if row["selection"] == "eligible":
            row["selection"] = "selected" if row["family_id"] in chosen else "not_sampled"
    for family, available in single_copy:
        for a, b in itertools.combinations(sorted(available), 2):
            if codes[a] == codes[b]:
                if (available[a], available[b]) in requests:
                    raise ValueError("Duplicate Ks pair request")
                requests[(available[a], available[b])] = (a, b, family, "single_copy_family_contrast")
            else:
                selection.append({"pair_id": "|".join((available[a], available[b])), "species_a": a,
                                  "species_b": b, "family_id": family, "pair_role": "single_copy_family_contrast",
                                  "status": "genetic_code_mismatch", "dS": "NA"})
    reports = estimate_pair_reports(output, sequences, codes, requests, cpus)
    contrasts = []
    for pair, info in requests.items():
        row = reports["|".join(pair)]
        selection.append({"pair_id": row["pair_id"], "species_a": info[0], "species_b": info[1],
                          "family_id": info[2], "pair_role": info[3], "status": row["status"], "dS": row["dS"]})
        if info[3] == "single_copy_family_contrast" and row["status"] == "ok":
            contrasts.append({"species_a": info[0], "species_b": info[1], "family_id": info[2], "ks": row["dS"]})
    write_table(output / "pair_selection.tsv", selection, ["pair_id", "species_a", "species_b", "family_id", "pair_role", "status", "dS"])
    write_table(output / "ks_family_selection.tsv", audit, ["family_id", "selection"])
    write_table(output / "family_ks_contrasts.tsv", contrasts, ["species_a", "species_b", "family_id", "ks"])
    return reports, contrasts


def estimate_pair_reports(output, sequences, codes, requests, cpus, prefix="", *, reuse=None, method="yn00"):
    # MAFFT's many short-lived files need not traverse a shared workspace mount.
    if prefix and not os.environ.get("MAFFT_TMPDIR") and requests:
        with tempfile.TemporaryDirectory(prefix="gg-wgd-mafft-", dir="/tmp") as scratch:
            return _estimate_pair_reports(output, sequences, codes, requests, cpus, prefix,
                                          {**os.environ, "MAFFT_TMPDIR": scratch}, reuse, method)
    return _estimate_pair_reports(output, sequences, codes, requests, cpus, prefix, None, reuse, method)


def _estimate_pair_reports(output, sequences, codes, requests, cpus, prefix, align_env, reuse, method):
    reports = {}
    by_code = defaultdict(list)
    for pair, info in requests.items():
        by_code[codes[info[0]]].append((pair, info))
    options = {"env": align_env} if align_env is not None else {}
    for code in sorted(set(codes.values())):
        selected = by_code[code]
        if not selected:
            continue
        pairs_file, ds_file = output / f"{prefix}aligned_pairs.code{code}.tsv", output / f"{prefix}pair_ds.code{code}.tsv"
        cached = matching_cached_rows(reuse[str(code)], selected, sequences) if reuse and str(code) in reuse else {}
        if cached:
            print(f"Cross-block dS code {code}: reusing {len(cached):,}/{len(selected):,} codon alignments", flush=True)

        def align(item, code=code, cached=cached):
            pair, info = item
            if pair in cached:
                return cached[pair]
            return align_pair(pair, [sequences[info[0]], sequences[info[1]]], code, **options)

        with ThreadPoolExecutor(max_workers=cpus) as pool:
            with closing(ordered_bounded_results(pool, align, ((item,) for item in selected), max_pending=cpus * 8)) as aligned:
                write_table(pairs_file, aligned, PAIR_COLUMNS)
        method_args = ["--method", method] if method != "yn00" else []
        run_tool(["cdskit", "dnds", "--pairs_file", pairs_file, *method_args,
                  "--codon_table", str(code), "--threads", str(cpus), "--out_file", ds_file], output, f"{prefix}dnds.code{code}")
        expected = {"|".join(pair) for pair, _ in selected}
        observed = set()
        statuses = DS_STATUSES if method == "yn00" else DS_STATUSES | CROSS_GY94_STATUSES
        for row in read_table(ds_file):
            pair_id = row.get("pair_id")
            if pair_id not in expected or pair_id in observed:
                raise ValueError("dS report contains an unrequested or duplicate Ks pair")
            if method in CROSS_DS_REPORT_LABELS and row.get("method") != CROSS_DS_REPORT_LABELS[method]:
                raise ValueError("CDSKIT did not report the requested cross-block dS model")
            if (not {"pair_id", "status", "dS"}.issubset(row) or row.get("status") not in statuses
                    or (row["status"] == "ok" and number(row.get("dS")) is None)):
                raise ValueError("dS report contains an invalid Ks status or value")
            observed.add(pair_id)
            reports[pair_id] = row
        if observed != expected:
            raise ValueError("dS report does not cover all requested Ks pairs")
    return reports


FOCAL_DS_STYLES = (
    ("Copy-pair dS", "#262E3D", "o"),
    ("Younger split", "#237A75", "|"),
    ("Older split", "#9B5937", "|"),
)


def draw_focal_ks_intervals(ax, y, row, matched=()):
    estimates = [("ok", (row["median_ds"], row["ci_lower"], row["ci_upper"]))]
    estimates.extend((row[f"{name}_boundary_status"], tuple(row[f"{name}_{field}"] for field in ("ds", "lower", "upper")))
                     for name in ("younger", "older"))
    has_estimate = False
    for shift, (_, color, marker), (status, estimate) in zip((0, .23, -.23), FOCAL_DS_STYLES, estimates, strict=True):
        center, lo, hi = map(number, estimate)
        diagnostic = status == "nonmonotone_or_invalid_interval"
        if status != "ok" and not diagnostic:
            continue
        if lo is not None and hi is not None and lo <= hi:
            ax.plot((lo, hi), (y + shift,) * 2, color=color, lw=1.4, marker="|" if not shift else None, ms=5,
                    ls="--" if diagnostic else "-")
            has_estimate = True
        if center is not None:
            ax.plot(center, y + shift, marker="x" if diagnostic else marker, color=color, ms=3 if not shift else 4,
                    ls="none", zorder=4)
            has_estimate = True
    for y_adjust, (_, color, _), estimate in zip((.40, -.40), FOCAL_DS_STYLES[1:], matched, strict=False):
        center, lo, hi = (number(estimate.get(key)) for key in ("matched_median_ds", "matched_lower", "matched_upper"))
        lane = y + y_adjust
        if lo is not None and hi is not None and lo <= hi:
            ax.plot((lo, hi), (lane, lane), color=color, lw=.9, ls=":")
        if center is not None:
            ax.plot(center, lane, marker="D", mfc="white", mec=color, ms=3, ls="none", zorder=4)
    if not has_estimate:
        ax.text(.02, y, "No dS estimate", transform=ax.get_yaxis_transform(), va="center", fontsize=8, color="#646D78")


def focal_ks_label(row):
    components = f", c={int(row['num_window_components']):,}" if "num_window_components" in row else ""
    return f"{row['focal']} (n={int(row['num_estimated_pairs']):,}{components})"


def plot_candidates(rows, output, tree=None, alpha=0.05, duplication_threshold=None, focal_rows=(), region_rows=(),
                    matched_trios=(), phase_diagnostics=(), common_outgroups=(), terminal_rows=(),
                    orientation_rows=(), copy_holdout=(), ds_method="yn00"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import ConnectionPatch, Patch

    tree_mode = duplication_threshold is not None
    cross_mode = any(row.get("synteny_method") == "cross_block_correspondence" for row in rows)
    block_definitions = {row.get("synteny_block_definition") for row in rows}
    both_orientations = "uniquely_concordant_four_copy_paths_in_both_orientations" in block_definitions
    strict_phase = both_orientations or "same_original_blocks_with_uniquely_concordant_four_copy_paths" in block_definitions
    paired_regions = {row["species_event_id"]: row for row in region_rows if row["mode"] == "paired"}
    baseline_cross_support = {row["species_event_id"] for row in region_rows
                              if row["mode"] == "all" and int(row["supports_before_holdout"])}
    phase_windows = unique_phase_windows(phase_diagnostics)
    shared_windows = defaultdict(int)
    for row in orientation_rows:
        if row["mode"] == "unique_phase" and row["status"] == "assessed_both":
            shared_windows[row["species_event_id"]] += min(row["shared_windows_a"], row["shared_windows_b"])
    age_styles = {"supports_branch": ("o", "#15504B", "dS supports placement"),
                  "partial_support": ("o", "#9F620A", "dS supports one daughter"),
                  "age_conflict": ("x", "#AD3232", "dS age conflict"),
                  "discordant_focals": ("x", "#AD3232", "dS age conflict")}
    metric = "duplication_score" if tree_mode else "likelihood_ratio_statistic"
    cutoff = (f"duplication score >= {duplication_threshold:g} (exploratory)" if tree_mode else
              f"search-wide P <= {alpha:g}")
    if strict_phase:
        cutoff += "; unique paths in both orientations" if both_orientations else "; unique paired paths"
    styles = {"duplication-and-synteny-candidate" if tree_mode else "WGD-supported": (
                  "#176B67", "", "Duplications + synteny" if tree_mode else "Counts + synteny"),
              "duplication-only-candidate" if tree_mode else "count-only-candidate": (
                  "#A34D22", "//", "Duplication-only candidate" if tree_mode else "Count-only candidate"),
              "unresolved": ("#777777", "", "Unresolved")}
    one_copy_style = "#176B67", "..", "Combined candidate; requires one-copy paths"
    phase_limited_style = "#826326", "//", "Assessed; unique-path support below floor"
    if both_orientations:
        phase_limited_style = "#826326", "//", "Assessed; shared unique-path support below floor"

    def phase_limited(row):
        return (strict_phase and row["event_support"] == "unresolved"
                and row.get("tree_support") == "tree_score_supported"
                and row.get("synteny_status") == "assessed_below_threshold"
                and row["species_event_id"] in baseline_cross_support)

    def row_style(row):
        state = row["event_support"]
        if phase_limited(row):
            return phase_limited_style
        region = paired_regions.get(row.get("species_event_id"))
        if (state == "duplication-and-synteny-candidate" and region is not None
                and not int(region["supports_before_holdout"])):
            return one_copy_style
        return styles[state]

    active_styles = [*styles.values()]
    if any(row_style(row) == one_copy_style for row in rows):
        active_styles.append(one_copy_style)
    if any(phase_limited(row) for row in rows):
        active_styles.append(phase_limited_style)
    shown = sorted(rows, key=lambda row: float(row[metric]), reverse=True)[:40]
    labels = []
    for row in reversed(shown):
        taxa = next(csv.reader([row["descendant_taxa"]]))
        label = f"B{row['branch_id']}: " + ", ".join(taxa[:2])
        if len(taxa) > 2:
            label += f" (+{len(taxa) - 2})"
        labels.append(label if len(label) <= 55 else label[:52] + "...")
    fig, ax = plt.subplots(figsize=(9, max(3, len(shown) * 0.25)))
    for i, row in enumerate(reversed(shown)):
        color, hatch, _ = row_style(row)
        ax.barh(i, float(row[metric]), color=color, hatch=hatch, edgecolor=color)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("Supporting / eligible gene families (not a P-value)" if tree_mode else
                  "Conditional count-model likelihood-ratio statistic")
    if tree_mode:
        ax.axvline(duplication_threshold, color="#555555", ls=":", lw=1)
    prefix = f"Top {len(shown)} of {len(rows)} " if len(shown) < len(rows) else ""
    ax.set_title(prefix + f"WGD candidates ({cutoff})", fontsize=11)
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(handles=[Patch(facecolor=c, hatch=h, edgecolor=c, label=label)
                       for c, h, label in active_styles], fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output / "wgd_candidates.pdf")
    plt.close(fig)
    if tree is None:
        return

    from nwkit.clade_index import CladeIndex

    clades = CladeIndex(tree)
    nodes = list(tree.traverse("preorder"))
    index = {clades.clade_id_for_node(node): node for node in nodes}
    events = {index[row["species_event_id"]]: row for row in rows}
    if len(events) != len(rows) or tree in events:
        raise ValueError("WGD candidates must map uniquely to non-root species branches")
    x = {tree: 0.0}
    for node in nodes[1:]:
        x[node] = x[node.up] + node.dist
    leaves = list(tree.leaves())
    y = {node: len(leaves) - i - 1 for i, node in enumerate(leaves)}
    for node in reversed(nodes):
        if not node.is_leaf:
            children = [y[child] for child in node.children]
            y[node] = (min(children) + max(children)) / 2
    by_node = defaultdict(list)
    trios = {(row["focal"], row["species_event_id"]): row for row in matched_trios}
    closed = {(row["focal"], row["species_event_id"]): row for row in copy_holdout
              if row["status"] == "ok" and row["drop_component_min"] < row["drop_component_max"]}
    gap_rows = defaultdict(list)
    for row in common_outgroups:
        if row["median_signed_gap"] is not None:
            gap_rows[row["focal"], row["species_event_id"]].append(row)
    matched_seen = False
    for row in focal_rows:
        by_node[index[row["species_event_id"]]].append(row)
    terminal = {index[row["species_event_id"]]: row for row in terminal_rows if row["median_copy_ds"] is not None}

    def leaf_label(leaf):
        row = terminal.get(leaf)
        return (leaf.name if row is None else
                f"{leaf.name} (n={row['six_pair_families']:,}, b={row['original_blocks']:,})")

    if by_node or terminal:
        # Insert room around each annotated branch without changing root distances.
        offset = 0
        for node in sorted(nodes, key=y.get):
            reserve = len(by_node.get(node, ())) * .58
            y[node] += offset + reserve
            offset += 2 * reserve
        names = [leaf_label(leaf) for leaf in leaves] + [focal_ks_label(row) for row in focal_rows]
        widths = (8.25, max(2.85, 0.087 * max(map(len, names))), 4.1)
        if gap_rows:
            gap_name_width = max(len(row["common_outgroup"]) for group in gap_rows.values() for row in group)
            widths += (3.15, max(1.4, .08 * (gap_name_width + 12)))
        height = max(4.8, (max(y.values()) + 1) * .26 + 2)
    else:
        widths = (8.25, 2.75)
        height = max(4, len(leaves) * .29 + 2)
    fig, axes = plt.subplots(ncols=len(widths), figsize=(sum(widths), height),
                             gridspec_kw={"width_ratios": widths, "wspace": 0.03})
    ax, labels_ax = axes[:2]
    labels_ax.sharey(ax)
    if by_node or terminal:
        ds_ax = axes[2]
        ds_ax.sharey(ax)
        axis_title = "dS (focal species scale)" if ds_method == "yn00" else f"dS ({ds_method}, focal species scale)"
        ds_ax.set(xlabel=axis_title)
        ds_ax.set_axisbelow(True)
        ds_ax.grid(axis="x", color="#E5E8EA", linewidth=.7)
        ds_ax.tick_params(axis="x", labelsize=8)
        ds_ax.tick_params(axis="y", left=False, labelleft=False)
        for side in ("top", "left", "right"):
            ds_ax.spines[side].set_visible(False)
    if gap_rows and by_node:
        gap_ax, gap_labels = axes[3:]
        gap_ax.sharey(ax)
        gap_labels.sharey(ax)
        gap_labels.set(xlim=(0, 1))
        gap_labels.set_axis_off()
        gap_ax.set(xlabel="Common outgroup: older − younger dS")
        gap_ax.axvline(0, color="#A3ABB2", lw=.85, ls=":")
        gap_ax.tick_params(axis="x", labelsize=8)
        for side in ("top", "left", "right"):
            gap_ax.spines[side].set_visible(False)
    for node in reversed(nodes):
        if not node.is_leaf:
            children = [y[child] for child in node.children]
            ax.vlines(x[node], min(y[node], *children), max(y[node], *children), color="#9AA5AA", lw=1)
        else:
            labels_ax.text(0, y[node], leaf_label(node), va="center", fontsize=8.5)
            if node in terminal:
                record = terminal[node]
                if record["copy_lower"] is not None:
                    ds_ax.hlines(y[node], record["copy_lower"], record["copy_upper"], color="#405A83", lw=1.25)
                ds_ax.plot(record["median_copy_ds"], y[node], "^", color="#405A83", ms=3.5, clip_on=False)
        if node is tree:
            continue
        row = events.get(node)
        state = row["event_support"] if row else "unresolved"
        color, hatch, _ = row_style(row) if row else styles[state]
        unobserved = (cross_mode and row and state == "unresolved" and
                      row["synteny_status"] == "insufficient_observability")
        if unobserved:
            color = "#6B537F"
        ax.plot((x[node.up], x[node]), (y[node], y[node]), color=color,
                ls=":" if unobserved else "-." if hatch == ".." else "--" if hatch else "-",
                lw=1.7 if unobserved else 1 if state == "unresolved" else 2.6)
        age = age_styles.get(row.get("cross_ks_status")) if row else None
        midpoint = (x[node.up] + x[node]) / 2
        if age:
            marker, edge, _ = age
            ax.plot(midpoint, y[node], marker=marker, mec=edge, mfc="white", ms=5.5, mew=1.2, zorder=4)
        if row and (state != "unresolved" or unobserved or phase_limited(row)):
            value = f"S={float(row[metric]):.2%}" if tree_mode else f"P={float(row['p_value']):.3g}"
            region = paired_regions.get(row["species_event_id"])
            windows = (f"  W={region['matched_windows']}" if region is not None
                       and (state == "duplication-and-synteny-candidate" or phase_limited(row)) else "")
            if windows and row["species_event_id"] in phase_windows:
                windows += f"  U={phase_windows[row['species_event_id']]}"
            if windows and both_orientations and row["species_event_id"] in shared_windows:
                windows += f"  R={shared_windows[row['species_event_id']]}"
            ax.annotate(f"B{row['branch_id']}  {value}{windows}",
                        (midpoint, y[node]), xytext=(0, 5),
                        textcoords="offset points", ha="center", fontsize=7.5, color=color)
        if node in by_node:
            focals = sorted(by_node[node], key=lambda r: r["focal"])
            # Continue the branch's own row to its dS tracks; this is not a tree edge.
            leader = ConnectionPatch((midpoint, y[node]), (.12, y[node]), "data",
                                     labels_ax.get_yaxis_transform(), axesA=ax, axesB=labels_ax,
                                     color="#AFB8BE", lw=.65, linestyle=":", zorder=1, clip_on=False)
            ax.add_artist(leader)
            lanes = [y[node] + (len(focals) / 2 - .5 - i) for i in range(len(focals))]
            labels_ax.plot((.12, .12), (lanes[-1], lanes[0]), color="#AFB8BE", lw=.75)
            for lane, focal in zip(lanes, focals, strict=True):
                labels_ax.plot((.12, .16), (lane, lane), color="#AFB8BE", lw=.75)
                labels_ax.text(.20, lane, focal_ks_label(focal), va="center", fontsize=8, color="#47545A")
                matched = ()
                if any(focal[f"{name}_boundary_status"] == "nonmonotone_or_invalid_interval" for name in ("younger", "older")):
                    matched = tuple(trios.get((focal["focal"], clades.clade_id_for_node(boundary_node)), {})
                                    for boundary_node in (node, node.up))
                    matched_seen |= any(number(estimate.get("matched_median_ds")) is not None for estimate in matched)
                draw_focal_ks_intervals(ds_ax, lane, focal, matched)
                if sensitivity := closed.get((focal["focal"], focal["species_event_id"])):
                    ds_ax.hlines(lane - .11, sensitivity["drop_component_min"], sensitivity["drop_component_max"],
                                 color="#78536F", lw=1.3)
                diagnostics = gap_rows.get((focal["focal"], focal["species_event_id"]), ())
                for j, diagnostic in enumerate(diagnostics):
                    ordinate = lane + ((len(diagnostics) - 1) / 2 - j) * .42
                    point = float(diagnostic["median_signed_gap"])
                    lo, hi = diagnostic["gap_lower"], diagnostic["gap_upper"]
                    color = {"ordered": "#15504B", "reversed": "#AD3232"}.get(diagnostic["status"], "#586471")
                    if lo is not None and hi is not None:
                        gap_ax.hlines(ordinate, lo, hi, color=color, lw=1.15)
                    gap_ax.plot(point, ordinate, ".", ms=5, color=color)
                    name, count = diagnostic["common_outgroup"], diagnostic["num_measured_all_outgroups"]
                    gap_labels.text(0, ordinate, f"{name}  n={count}", va="center", fontsize=6.5, color=color)
    span = max(x.values()) or 1
    labels_ax.set(xlim=(0, 1))
    labels_ax.set_axis_off()
    ax.set(xlim=(-span * .02, span * 1.02), ylim=(-1, max(y.values()) + 1), yticks=[],
           xlabel="Distance from root (input branch units)")
    for side in ("top", "left", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(f"WGD candidate branches | {cutoff}", loc="left")
    handles = [Line2D([], [], color=c, ls="-." if h == ".." else "--" if h else "-", label=label)
               for c, h, label in active_styles]
    if cross_mode:
        handles.append(Line2D([], [], color="#6B537F", ls=":", label="Synteny unobserved; score below cutoff"))
    if phase_windows:
        handles.append(Line2D([], [], color="none", label="W: paired windows; U: uniquely concordant windows"))
    if both_orientations and shared_windows:
        handles.append(Line2D([], [], color="none", label="R: shared unique windows (smaller reference count per pair)"))
    shown_age = {row.get("cross_ks_status") for row in rows}
    age_handles = [Line2D([], [], marker=marker, mec=edge, mfc="white", ls="none", label=label)
                   for marker, edge, label in dict.fromkeys(style for state, style in age_styles.items() if state in shown_age)]
    handles.extend(age_handles)
    legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -.11),
                       ncols=2 if cross_mode else 3, fontsize=8, frameon=False)
    if by_node or terminal:
        ds_ax.set_xlim(left=0, right=max(ds_ax.get_xlim()[1], .1) if ds_ax.has_data() else 1)
        ds_handles = ([Line2D([], [], color=color, marker=marker, lw=1.6, label=label)
                       for label, color, marker in FOCAL_DS_STYLES] if by_node else [])
        if any(row[f"{name}_boundary_status"] == "nonmonotone_or_invalid_interval"
               and number(row[f"{name}_ds"]) is not None for row in focal_rows for name in ("younger", "older")):
            ds_handles.append(Line2D([], [], color="#646D78", marker="x", ls="--", label="Unusable split"))
        if matched_seen:
            ds_handles.append(Line2D([], [], color="#646D78", marker="D", mfc="white", ls=":", label="Matched-family split"))
        if terminal:
            ds_handles.append(Line2D([], [], color="#405A83", marker="^", label="Terminal self anchors (b: blocks)"))
        if closed:
            ds_handles.append(Line2D([], [], color="#78536F", label="Range omitting one component"))
        ds_ax.legend(handles=ds_handles, loc="lower left", bbox_to_anchor=(0, 1.0),
                     ncols=2 if matched_seen or terminal or closed else 3, fontsize=8, frameon=False)
    top = min(.92, 1 - .95 / height) if closed else min(.92, 1 - .80 / height) if terminal else (
        min(.92, 1 - .64 / height) if matched_seen else .92)
    bottom = 0.19 if age_handles else 0.15
    fig.subplots_adjust(left=0.07, right=0.97, bottom=bottom, top=top)
    if both_orientations:
        fig.canvas.draw()
        low = legend.get_window_extent(fig.canvas.get_renderer()).y0 / fig.bbox.height
        if low < .02:
            # The legend is anchored to an axis-relative position below the
            # panel. Increasing the bottom margin also shortens that panel.
            fig.subplots_adjust(bottom=bottom + (.02 - low) / 1.11)
    for extension in ("pdf", "png"):
        fig.savefig(output / f"wgd_species_tree.{extension}", dpi=175)
    plt.close(fig)


def observable_synteny(placer, summaries, anchors):
    species_with_ks = {row["species"] for row in anchors
                       if row["ks_status"] == "ok" and number(row["ks"]) is not None}
    opportunity = defaultdict(set)
    for name, summary in summaries.items():
        # A successful zero-block scan is an observation, unlike missing
        # sequences or blocks that cannot be placed on the species tree.
        if summary.get("num_blocks") == 0:
            branches = placer.lineage_branches(name)
        elif name in species_with_ks:
            branches = placer.observable_branches(name)
        else:
            continue
        for branch in branches:
            opportunity[branch].add(name)
    return opportunity


def run_analysis(plan, output, cpus):
    if cpus < 1:
        raise ValueError("cpus must be positive")
    p = plan["parameters"]
    tree_mode = bool(plan.get("gene_trees"))
    count_input_args = []
    for key, option in (("count_detection_tsv", "--detection-tsv"), ("count_rate_groups_tsv", "--rate-groups-tsv")):
        if plan.get(key):
            if tree_mode or plan[key] not in plan["input_hashes"]:
                raise ValueError("Count observation and rate-group tables require hashed inputs in gene-count mode")
            count_input_args.extend((option, plan[key]))
    for key, option in (("count_event_fractions", "--event-fractions"), ("count_family_gamma_shape", "--family-gamma-shape")):
        if p.get(key):
            if tree_mode:
                raise ValueError("Count event fractions and family gamma shape require gene-count mode")
            count_input_args.extend((option, p[key]))
    bootstrap_workers = p.get("count_bootstrap_workers", 1)
    scan_workers = p.get("count_scan_workers", 1)
    if not tree_mode and scan_workers > cpus:
        raise ValueError("WGD scan workers exceed the allocated CPUs")
    if not tree_mode and p["count_bootstrap"] and bootstrap_workers > cpus:
        raise ValueError("WGD bootstrap workers exceed the allocated CPUs")
    verify(plan)
    protect_output(output, plan["input_hashes"], Path(plan["workspace"]))
    checkpoint = plan.get("count_bootstrap_checkpoint", "")
    if checkpoint:
        if tree_mode or not p["count_bootstrap"]:
            raise ValueError("A bootstrap checkpoint requires gene-count mode and positive --count-bootstrap")
        protect_checkpoint(checkpoint, plan["input_hashes"] | {str(output.resolve()): ""}, Path(plan["workspace"]))
        if output.resolve() in Path(checkpoint).resolve().parents:
            raise ValueError("WGD bootstrap checkpoint must be outside the output bundle")
    family_tables = _validated_family_tables(plan)
    members, family_key, count_rows = family_tables
    validate_annotation_counts(plan, members)
    tree = species_tree(plan["species_tree"])
    output.mkdir(parents=True, exist_ok=False)
    normalized_counts(plan, output, tree, count_rows, minimum_families=0 if tree_mode else 2)
    cross_evidence = None
    cross_ks = terminal = rooted_copy = None
    region_summary = ()
    phase_diagnostics = ()
    traced_ks_hits = ()
    orientation = ()
    plot_names = ["wgd_candidates.pdf", "wgd_species_tree.pdf", "wgd_species_tree.png"]
    if tree_mode:
        age_sources = plan.get("cross_ks_sources", {})
        cross = reverse = None
        pairs = plan.get("cross_pairs")
        if pairs is not None:
            audit_orientation = p.get("cross_orientation_audit")
            terminal_limit = p.get("max_terminal_ks_families", 0)
            inputs = cross_synteny_inputs(pairs) if audit_orientation or terminal_limit else None
            cross = CrossSynteny(tree, pairs, p["min_cross_windows"], age_sources,
                                 trace_branches=None, audit_coverage=True, audit_phase=True, inputs=inputs,
                                 reverse_phase=bool(p.get("cross_require_unique_phase")) if audit_orientation else None)
            if audit_orientation:
                reverse, cross.reverse = cross.reverse, None
            if terminal_limit:
                terminal = TerminalKs(tree, members, family_key, age_sources, inputs, terminal_limit, p["seed"])
            del inputs
        by_family = {row[family_key]: row for row in members} if cross is not None else {}

        def family_members(family):
            row = by_family.get(family)
            return (None if row is None else
                    {gene.strip(): name for name in plan["species"] if row[name].strip() for gene in row[name].split(",")})

        cross_ks = CrossKs(tree, age_sources, p, family_tables) if age_sources else None
        if p.get("ks_rooted_copies"):
            rooted_copy = RootedCopyKs(tree, cross_ks.codes)

        def observe_family(family, leaves):
            cross.observe_family(family, leaves)
            if cross_ks is not None:
                cross_ks.observe_family(family, leaves)

        support_file = plan.get("cross_tree_support")
        paired_hits = set() if cross is not None and (support_file or rooted_copy is not None) else None

        def observe_cross(branch, family, node):
            cross.observe(branch, family, node, paired_hits)
            if reverse is not None:
                reverse.observe(branch, family, node)

        def observe_support(family, gene_tree, duplications, membership):
            if support:
                support.observe(family, gene_tree, duplications, paired_hits or ())
            if paired_hits:
                if rooted_copy is not None:
                    rooted_copy.observe(family, gene_tree, membership, paired_hits)
                paired_hits.clear()

        holdout_rows = [] if p.get("duplication_taxon_holdout") else None
        duplication_coverage = [] if p.get("duplication_coverage_audit") else None
        root_audit = [] if p.get("duplication_root_audit") else None
        copy_splits = [] if p.get("duplication_copy_split_audit") else None
        panel_count = p.get("duplication_family_panels", 0)
        panels = FamilyPanels(panel_count, p["seed"]) if panel_count else None
        panel_output = []
        with DuplicationTreeSupport(support_file) if support_file else nullcontext() as support:
            raw_candidates, num_trees = scan_duplications(
                tree, plan["gene_trees"], p["min_duplication_score"],
                observe=observe_cross if cross is not None else None,
                family_members=family_members if cross is not None else None,
                observe_family=observe_family if cross is not None else None,
                observe_duplications=observe_support if support or rooted_copy is not None else None, holdout_output=holdout_rows,
                coverage_output=duplication_coverage, root_output=root_audit,
                observe_topology=terminal.observe if terminal is not None else None,
                panels=panels, panel_output=panel_output, copy_split_output=copy_splits)
        if holdout_rows is not None:
            write_table(output / "duplication_taxon_holdout.tsv", holdout_rows, HOLDOUT_FIELDS)
        if duplication_coverage is not None:
            write_table(output / "duplication_coverage.tsv", duplication_coverage, DUPLICATION_COVERAGE_FIELDS)
        if root_audit is not None:
            write_table(output / "duplication_root_orientations.tsv", root_audit, ROOT_FIELDS)
        if copy_splits is not None:
            write_table(output / "duplication_copy_splits.tsv", copy_splits, COPY_SPLIT_FIELDS)
        if panels is not None:
            detail, summary = panel_output
            write_table(output / "duplication_family_panels.tsv", detail, PANEL_FIELDS)
            write_table(output / "duplication_family_panel_summary.tsv", summary, PANEL_SUMMARY_FIELDS)
        if support:
            write_table(output / "duplication_tree_support.tsv", support.rows(raw_candidates), DUPLICATION_SUPPORT_FIELDS)
        del by_family
        if cross is not None:
            if reverse is None:
                coverage_rows = list(cross.coverage_rows())
            else:
                coverage_rows, nonreuse_rows, reverse_floor = cross.nonreused_coverage_rows()
                write_table(output / "cross_synteny_anchor_nonreuse.tsv", nonreuse_rows, NONREUSE_OUTPUT)
            cross_evidence, pair_rows, window_rows = cross.summarize(
                archive_floor=archive_window_floor(coverage_rows),
                require_unique_phase=p.get("cross_require_unique_phase", 0))
            write_table(output / "cross_synteny_pairs.tsv", pair_rows, PAIR_OUTPUT)
            phase_diagnostics = list(phase_rows(raw_candidates, pair_rows, cross.phases))
            write_table(output / "cross_synteny_copy_phase.tsv", phase_diagnostics, PHASE_FIELDS)
            both_hits = None
            reverse_windows = None
            if reverse is not None:
                reverse_evidence, reverse_pairs, _ = reverse.summarize(include_windows=False, archive_floor=reverse_floor)
                write_table(output / "cross_synteny_orientation.tsv",
                            orientation_rows(raw_candidates, pair_rows, reverse_pairs, p["min_cross_windows"]),
                            ORIENTATION_OUTPUT)
                orientation = list(orientation_phase_overlap(raw_candidates, pair_rows, reverse_pairs, cross, reverse,
                                                             unique=bool(p.get("cross_require_unique_phase"))))
                write_table(output / "cross_synteny_orientation_quartets.tsv", orientation, ORIENTATION_PHASE_FIELDS)
                if p.get("cross_require_both_orientations"):
                    reverse_windows = {}
                    both_hits = require_paired_orientations(cross_evidence, cross.daughters, orientation, cross,
                                                           reverse, p["min_cross_windows"], reverse_windows=reverse_windows,
                                                           reverse_evidence=reverse_evidence)
                del reverse
            write_table(output / "cross_synteny_windows.tsv", window_rows, WINDOW_OUTPUT)
            write_table(output / "duplication_synteny_family_overlap.tsv",
                        duplication_synteny_overlap(raw_candidates, window_rows), OVERLAP_FIELDS)
            write_table(output / "cross_synteny_archive_coverage.tsv", coverage_rows, COVERAGE_OUTPUT)
            unique_hits = ({hit for hit in cross.traced_hits if trace_has_unique_phase(cross.phases, hit)}
                           if p.get("cross_require_unique_phase") else None)
            write_table(output / "cross_synteny_copy_replication.tsv",
                        replicated_copy_rows(raw_candidates, pair_rows, cross.traced_hits, unique_hits), REPLICATION_FIELDS)
            region_summary, region_detail = summarize_regions(raw_candidates, cross.traced_hits, cross.daughters,
                                                               p["min_cross_windows"], unique_phase_hits=unique_hits,
                                                               reverse_windows=reverse_windows)
            if reverse_windows is not None:
                for row in region_summary:
                    if (row["mode"] == "both_orientations" and
                            bool(row["supports_before_holdout"]) != cross_evidence[row["species_event_id"]]["supported"]):
                        raise ValueError("Both-reference region audit changed the pre-holdout support")
            write_table(output / "cross_synteny_region_summary.tsv", region_summary, REGION_SUMMARY)
            write_table(output / "cross_synteny_regions.tsv", region_detail, REGION_DETAIL)
            write_table(output / "cross_synteny_nested_copies.tsv",
                        nested_copy_overlap(tree, raw_candidates, cross.traced_hits), NESTED_FIELDS)
            if support:
                write_table(output / "cross_synteny_tree_support.tsv", support.paired_rows(raw_candidates), SUPPORT_FIELDS)
            if cross_ks:
                age_hits = unique_hits if both_hits is None else both_hits
                age_pairs = (cross.age_pairs if age_hits is None else
                             CrossKs.age_pairs_from_trace(age_hits, cross.age_species))
                cross_ks.prepare(age_pairs, raw_candidates, cross_evidence)
                if p.get("max_ks_quartets"):
                    traced_ks_hits = cross.traced_hits if age_hits is None else age_hits
            del cross
        if not raw_candidates:
            raise ValueError("Duplication scoring needs at least one non-root internal species branch")
        write_table(output / "duplication_candidates.tsv", raw_candidates, list(raw_candidates[0]))
    else:
        worker_args = ["--bootstrap-workers", str(bootstrap_workers)] if p["count_bootstrap"] and bootstrap_workers > 1 else []
        if scan_workers > 1:
            worker_args.extend(("--scan-workers", str(scan_workers)))
        if checkpoint:
            worker_args.extend(("--bootstrap-checkpoint", checkpoint))
        run_tool([sys.executable, "-m", "nwkit", "wgd-count", "--infile", plan["species_tree"],
                  "--input-rooted", "auto", "--counts", output / "counts.tsv", "--ascertainment", "root-clades",
                  "--bootstrap", p["count_bootstrap"], "--seed", p["seed"], "--alpha", p["alpha"],
                  "--max-states", p["max_states"], "--max-iterations", p["max_iterations"],
                  "--multiplicity", p["multiplicity"], "--outfile", output / "count_candidates.tsv",
                  "--model-out", output / "count_model.json", *worker_args, *count_input_args], output, "wgd-count")
        raw_candidates = read_table(output / "count_candidates.tsv")
        num_trees = None
    positions, sequences, summaries, anchors = prepare_sources(plan, output, cpus)
    reports, contrasts = estimate_ks(plan, output, sequences, anchors, cpus, family_tables)
    del family_tables, members, count_rows
    bounds = {}
    if contrasts:
        run_tool([sys.executable, "-m", "nwkit", "ksrate", "--infile", plan["species_tree"], "--input-rooted", "auto",
                  "--ks-tsv", output / "family_ks_contrasts.tsv", "--bootstrap", p["ks_bootstrap"], "--seed", p["seed"],
                  "--ci-method", "pair-median-bonferroni",
                  "--outfile", output / "ks_boundaries.tsv", "--trios-out", output / "ks_trios.tsv",
                  "--model-out", output / "ks_model.json"], output, "ksrate")
        bounds = {(row["focal"], row["species_event_id"]): row for row in read_table(output / "ks_boundaries.tsv")}
    ks_placer = KsBranchPlacer(tree, bounds)
    for row in anchors:
        report = reports.get(row["gene_a"] + "|" + row["gene_b"])
        row["ks"] = report["dS"] if report else "NA"
        row["ks_status"] = report["status"] if report else (
            "cds_unavailable" if any(gene not in sequences.get(row["species"], {}) for gene in (row["gene_a"], row["gene_b"]))
            else "not_sampled")
        branch, status = ks_placer.place(row["species"], number(row["ks"])) if row["ks_status"] == "ok" else (None, "missing_ks")
        row.update(species_event_id=branch or "", placement_status=status)
    write_table(output / "anchor_evidence.tsv", anchors,
                ["species", "block_id", "gene_a", "gene_b", "ks", "ks_status", "species_event_id", "placement_status"])
    daughter_taxa = None
    if tree_mode:
        from nwkit.clade_index import CladeIndex

        clades = CladeIndex(tree)
        daughter_taxa = {clades.clade_id_for_node(node): tuple(set(child.leaf_names()) for child in node.children)
                         for node in tree.descendants() if not node.is_leaf}
    candidates = summarize_events(raw_candidates, anchors, summaries,
                                  p["min_coverage"], p["min_blocks"], {row["gene_id"]: row for row in positions},
                                  synteny_opportunity=observable_synteny(ks_placer, summaries, anchors),
                                  daughter_taxa=daughter_taxa, cross_evidence=cross_evidence)
    if cross_evidence is not None:
        for row in candidates:
            row["synteny_method"] = "cross_block_correspondence"
            row["minimum_cross_windows"] = p["min_cross_windows"]
            if p.get("cross_require_both_orientations"):
                row["synteny_block_definition"] = "uniquely_concordant_four_copy_paths_in_both_orientations"
            elif p.get("cross_require_unique_phase"):
                row["synteny_block_definition"] = "same_original_blocks_with_uniquely_concordant_four_copy_paths"
    age_summary, focal_rows = {}, ()
    if cross_ks:
        ds_method = p.get("cross_ks_method", "yn00")
        estimate = partial(estimate_pair_reports, reuse=plan.get("reused_cross_ks_alignments"),
                           method=ds_method)
        focal_rows = cross_ks.run(output, cpus, candidates, estimate, write_table, traced_ks_hits, terminal=terminal,
                                  rooted_copy=rooted_copy, run_tool=run_tool,
                                  method_label=CROSS_DS_REPORT_LABELS.get(ds_method, plan["tools"]["ds"]["method"]))
        plot_names += plot_cross_ks(focal_rows, output)
        for row in candidates:
            status = row["cross_ks_status"]
            age_summary[status] = age_summary.get(status, 0) + 1
    write_table(output / "wgd_events.tsv", candidates, list(candidates[0]))
    plot_candidates(candidates, output, tree=tree, alpha=p["alpha"],
                    duplication_threshold=p["min_duplication_score"] if tree_mode else None,
                    focal_rows=focal_rows, region_rows=region_summary,
                    matched_trios=read_table(output / "cross_ks_matched_trios.tsv") if cross_ks else (),
                    phase_diagnostics=phase_diagnostics, common_outgroups=cross_ks.outgroup_rows if cross_ks else (),
                    terminal_rows=cross_ks.terminal_rows if cross_ks else (), orientation_rows=orientation,
                    copy_holdout=cross_ks.copy_holdout_rows if cross_ks else (), ds_method=p.get("cross_ks_method", "yn00"))
    verify(plan)
    write_json(output / "summary.json", {"schema_version": 1, "plan": plan,
        "output_hashes": {path.name: digest(path) for path in sorted(output.iterdir())
                          if path.is_file() and path.suffix in {".tsv", ".json"}},
        "plot_hashes": {name: digest(output / name) for name in plot_names},
        "num_positioned_genes": len(positions),
        "num_raw_anchor_rows": len(anchors), "num_ks_contrasts": len(contrasts), "species_synteny": summaries,
        "ks_boundary_ci_methods": sorted({row.get("ci_method", "not_reported") for row in bounds.values()}),
        "num_WGD_supported_candidates": sum(row["event_support"] == "WGD-supported" for row in candidates),
        "num_count_only_candidates": sum(row["event_support"] == "count-only-candidate" for row in candidates),
        "num_duplication_only_candidates": sum(row["event_support"] == "duplication-only-candidate" for row in candidates),
        "num_duplication_and_synteny_candidates": sum(
            row["event_support"] == "duplication-and-synteny-candidate" for row in candidates),
        "num_input_gene_trees": num_trees,
        **({"cross_ks_branch_status_counts": age_summary} if cross_ks else {}),
        "synteny_method": "cross_block_correspondence" if cross_evidence is not None else "self_synteny_ks",
        "candidate_model": "rooted_gene_tree_duplication_score" if tree_mode else "native_gene_count",
        "interpretation": ("Exploratory rooted gene-tree duplication scores. " if tree_mode else
                           "Experimental one-event-at-a-time gene-count scan. ") +
                          "Branch candidates are not gene-node origin labels or posterior probabilities. "
                          "Missing synteny is not SSD evidence."})


def same_species_cross_child_pairs(node, positions):
    """Select positioned pairs in the original left-leaf/right-leaf order."""
    right_by_species = defaultdict(list)
    for leaf in node.children[1].leaves():
        if leaf.name in positions:
            right_by_species[positions[leaf.name]["species"]].append(leaf.name)
    pairs = []
    for leaf in node.children[0].leaves():
        if leaf.name in positions:
            pairs.extend((leaf.name, right) for right in right_by_species.get(positions[leaf.name]["species"], ()))
    return pairs


def _load_family_evidence(evidence, gene_ids):
    """Validate complete tables while retaining only this family's evidence."""
    def rows(path, keys):
        with Path(path).open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            fields = next(reader, [])
            if not fields or len(set(fields)) != len(fields) or any(not field for field in fields):
                raise ValueError(f"Invalid TSV header: {path}")
            indices = None
            for row in reader:
                # DictReader skips physically empty rows, but a quoted empty
                # field remains a real row and must still pass the width check.
                if not row:
                    continue
                if len(row) != len(fields):
                    raise ValueError(f"Invalid TSV row width: {path}")
                if indices is None:
                    # Match DictReader's key errors only once a row is present,
                    # after its width has been validated.
                    for key in keys:
                        if key not in fields:
                            raise KeyError(key)
                    indices = tuple(fields.index(key) for key in keys)
                yield fields, row, indices

    positions, seen = {}, set()
    duplicate_gene = False
    for fields, values, indices in rows(evidence / "gene_positions.tsv", ("gene_id",)):
        gene = values[indices[0]]
        duplicate_gene |= gene in seen
        seen.add(gene)
        if gene in gene_ids:
            positions[gene] = dict(zip(fields, values, strict=True))
    event_rows = read_table(evidence / "wgd_events.tsv")
    events = {row["species_event_id"]: row for row in event_rows}
    if duplicate_gene or len(events) != len(event_rows):
        raise ValueError("Duplicate gene or species event identifier in WGD evidence")
    del seen
    anchors = defaultdict(list)
    for fields, values, indices in rows(evidence / "anchor_evidence.tsv", ("gene_a", "gene_b")):
        left, right = values[indices[0]], values[indices[1]]
        # Invalid or unmatched anchors also participate in conflict detection.
        if left in gene_ids and right in gene_ids:
            anchors[tuple(sorted((left, right)))].append(dict(zip(fields, values, strict=True)))
    return positions, events, anchors


def classify(args):
    from nwkit.clade_index import CladeIndex
    from nwkit.rooting_state import require_rooted
    from nwkit.util import _is_missing_support_value, assign_branch_ids, read_tree, write_tree

    mul_enabled = bool(getattr(args, "mul_diagnostics", 0))
    if mul_enabled and (not args.mul_h1.strip() or min(args.mul_max_candidates, args.mul_max_state_pairs, args.mul_max_maps) < 1):
        raise ValueError("MUL diagnostics require explicit H1 and positive resource limits")
    evidence = args.evidence.resolve()
    # Ensure the candidate IDs refer to this full tree, not a family-pruned tree.
    summary = json.loads((evidence / "summary.json").read_text())
    if digest(args.species_tree) != summary["plan"]["input_hashes"][str(Path(summary["plan"]["species_tree"]))]:
        raise ValueError("WGD evidence and classification species trees differ")
    verify(summary["plan"])
    required = {"gene_positions.tsv", "anchor_evidence.tsv", "wgd_events.tsv"}
    if args.native_tree_likelihood:
        if summary.get("candidate_model") == "rooted_gene_tree_duplication_score":
            raise ValueError("Native tree likelihood requires a fitted count model, not duplication-score candidates")
        required.add("count_model.json")
    if not required.issubset(summary.get("output_hashes", {})):
        raise ValueError("WGD evidence lacks output integrity hashes; regenerate it")
    for name, expected in summary["output_hashes"].items():
        if Path(name).name != name or digest(evidence / name) != expected:
            raise ValueError("WGD evidence output changed or has an invalid hash path")
    snapshot_paths = [args.gene_tree, args.species_tree, evidence / "summary.json"]
    if args.species_map:
        snapshot_paths.append(args.species_map)
    snapshot = {str(path): digest(path) for path in snapshot_paths}
    protect_output(args.output, [args.gene_tree, args.species_tree, evidence],
                   Path(summary["plan"]["workspace"]) if "workspace" in summary["plan"] else None)
    species = species_tree(args.species_tree)
    if mul_enabled:
        from nwkit.mul_reconcile_model import validate_binary

        validate_binary(species, "MUL diagnostic species tree")
    tree = read_tree(str(args.gene_tree), "auto", True, rooted="auto")
    require_rooted(tree, "WGD classification requires a rooted gene tree.")
    args.output.mkdir(parents=True, exist_ok=False)
    for node in tree.traverse():
        for prop in ("duplication_origin", "conditional_wgd_probability", "mul_mapping_status", "mul_dl_duplication",
                     "mul_gene_node", "mul_best_hypotheses", "mul_optimal_mappings"):
            node.props.pop(prop, None)
    event_source = "nhx" if any("D" in node.props or "H" in node.props for node in tree.traverse()) else "lca"
    command = [sys.executable, "-m", "nwkit", "reconcile", "--infile", args.gene_tree,
               "--species-tree", args.species_tree, "--event-source", event_source, "--unmatched", "ignore",
               "--tree-id", args.family_id, "--species-parser", args.species_parser, "--outfile", args.output / "reconciliation.tsv"]
    if args.species_regex:
        command += ["--species-regex", args.species_regex]
    if args.species_map:
        command += ["--species-map-tsv", args.species_map]
    run_tool(command, args.output, "reconcile")
    native = {}
    if args.native_tree_likelihood:
        native_command = [sys.executable, "-m", "nwkit", "wgd-tree", "--infile", args.gene_tree,
                          "--input-rooted", "auto", "--species-tree", args.species_tree,
                          "--count-model", evidence / "count_model.json", "--tree-id", args.family_id,
                          "--species-parser", args.species_parser, "--outfile", args.output / "native_tree_origins.tsv",
                          "--model-out", args.output / "native_tree_likelihood.json"]
        if args.species_regex:
            native_command += ["--species-regex", args.species_regex]
        if args.species_map:
            native_command += ["--species-map-tsv", args.species_map]
        run_tool(native_command, args.output, "wgd-tree")
        native = {(row["species_event_id"], row["gene_clade_id"]): row
                  for row in read_table(args.output / "native_tree_origins.tsv")}
    positions, events, anchors = _load_family_evidence(evidence, {leaf.name for leaf in tree.leaves()})
    clades, branches = CladeIndex(tree), assign_branch_ids(tree)
    nodes = {clades.clade_id_for_node(node): node for node in tree.traverse()}
    result, pair_rows = [], []
    reconciliation = read_table(args.output / "reconciliation.tsv")
    leaf_mapping = {row["gene_name"]: row for row in reconciliation if row["event_type"] == "leaf"}
    species_clades = CladeIndex(species)
    terminal_branches = {species_clades.clade_id_for_node(leaf) for leaf in species.leaves()}
    for row in reconciliation:
        if row["event_type"] != "duplication":
            continue
        node = nodes[row["gene_clade_id"]]
        pairs = same_species_cross_child_pairs(node, positions)
        status, reason, features, supported = combine_node(row["species_event_id"], pairs, positions, anchors, events,
            args.proximal_distance, terminal_cherry=len(list(node.leaves())) == 2 and row["species_event_id"] in terminal_branches)
        if row["mapping_status"] != "mapped" or row["event_status"] != "resolved":
            status, reason = "unresolved", "unresolved_reconciliation"
        if any(leaf.name not in positions or not valid_position(positions[leaf.name]) for leaf in node.leaves()):
            status, reason = "unresolved", "incomplete_genomic_mapping"
        elif any(leaf_mapping.get(leaf.name, {}).get("mapping_status") != "mapped"
                 or leaf_mapping[leaf.name]["species_name"] != positions[leaf.name]["species"] for leaf in node.leaves()):
            status, reason = "unresolved", "genomic_species_mapping_conflict"
        else:
            loci = [(positions[leaf.name]["species"], positions[leaf.name]["locus_id"]) for leaf in node.leaves()]
            if len(loci) != len(set(loci)):
                status, reason = "unresolved", "same_locus_annotation_ambiguity"
        node.add_prop("duplication_origin", status)
        probability = native.get((row["species_event_id"], row["gene_clade_id"]), {}).get("conditional_wgd_probability", "NA")
        if probability != "NA":
            node.add_prop("conditional_wgd_probability", probability)
        support = node.props.get("support")
        if support is None or _is_missing_support_value(support):
            support = "NA"
        result.append({"family_id": args.family_id, "gene_branch_id": branches[node], "gene_clade_id": row["gene_clade_id"],
                       "species_event_id": row["species_event_id"], "classification": status, "reason": reason,
                       "event_source": event_source, "num_same_species_cross_child_pairs": len(pairs),
                       "num_anchor_rows": len(supported), "input_tree_support": support,
                       "native_conditional_wgd_probability": probability,
                       "score_meaning": "evidence_rule_not_posterior"})
        pair_rows.extend({"gene_clade_id": row["gene_clade_id"], "gene_a": a, "gene_b": b,
                          "position_feature": feature, "gene_rank_distance": distance if distance is not None else "NA"}
                         for a, b, feature, distance in features)
    fields = ["family_id", "gene_branch_id", "gene_clade_id", "species_event_id", "classification", "reason", "event_source",
              "num_same_species_cross_child_pairs", "num_anchor_rows", "input_tree_support", "native_conditional_wgd_probability", "score_meaning"]
    write_table(args.output / "duplication_origins.tsv", result, fields)
    write_table(args.output / "pair_evidence.tsv", pair_rows, ["gene_clade_id", "gene_a", "gene_b", "position_feature", "gene_rank_distance"])
    mul_summary = run_diagnostics(args, tree, species, reconciliation, result, run_tool, read_table, write_table) if mul_enabled else {"enabled": False}
    properties = {key for node in tree.traverse() for key in node.props} - {"name", "dist", "support"}
    write_tree(tree, argparse.Namespace(outfile=str(args.output / "classified_gene_tree.nhx")),
               "auto", quiet=True, props=sorted(properties))
    plot_origins(tree, args.output)
    verify(summary["plan"])
    verify({"input_hashes": snapshot})
    verify({"input_hashes": {str(evidence / name): expected for name, expected in summary["output_hashes"].items()}})
    write_json(args.output / "summary.json", {"schema_version": 1, "num_duplications": len(result),
        "class_counts": {status: sum(row["classification"] == status for row in result) for status in ("WGD-supported", "SSD-supported", "unresolved")},
        "native_tree_likelihood": bool(args.native_tree_likelihood),
        "mul_diagnostics": mul_summary,
        "interpretation": "Conditional on one input gene tree. Origin labels are inspectable evidence support, not probabilities. Tandem adjacency supports an SSD hypothesis but does not exclude WGD-derived copies relocated by rearrangement. Broad segmental duplication can mimic combined WGD evidence. Optional native probabilities condition on fixed topology, supplied count parameters and one event, not WGD occurrence or parameter/tree uncertainty."})


def plot_origins(tree, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {"WGD-supported": "#197c53", "SSD-supported": "#be4d31", "unresolved": "#777777"}
    leaves = list(tree.leaves())
    y = {leaf: index for index, leaf in enumerate(reversed(leaves))}
    x = {tree: 0}
    for node in tree.traverse("preorder"):
        for child in node.children:
            x[child] = x[node] + 1
    longest_name = max(len(leaf.name) for leaf in leaves)
    fig, ax = plt.subplots(figsize=(max(10, longest_name * 0.07 + 5), max(3, len(leaves) * 0.22)))
    for node in tree.traverse("postorder"):
        if node.children:
            y[node] = sum(y[child] for child in node.children) / len(node.children)
            ax.plot([x[node], x[node]], [min(y[child] for child in node.children), max(y[child] for child in node.children)], color="#999999", linewidth=0.8)
        if node is not tree:
            ax.plot([x[node.up], x[node]], [y[node], y[node]], color=colors.get(node.props.get("duplication_origin"), "#333333"), linewidth=1.4)
        if node.is_leaf:
            ax.annotate(node.name, (x[node], y[node]), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
        elif "duplication_origin" in node.props:
            ax.plot(x[node], y[node], "o", color=colors[node.props["duplication_origin"]], markersize=4)
    depth = max(1, max(x.values()))
    ax.set_xlim(-0.03 * depth, depth * (1 + max(0.25, min(1.5, longest_name / 60))))
    ax.set_ylim(-1, len(leaves))
    ax.axis("off")
    ax.legend(handles=[Line2D([], [], color=color, marker="o", label=status) for status, color in colors.items()], loc="upper left", bbox_to_anchor=(0, 1.08), frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "duplication_origins.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--workspace", required=True, type=Path)
    for field in ("species-tree", "counts", "members", "genomes", "gene-trees", "cross-synteny", "cross-ks-cds",
                  "cross-tree-support", "count-bootstrap-checkpoint", "count-detection-tsv",
                  "count-rate-groups-tsv", "count-event-fractions", "count-family-gamma-shape"):
        plan.add_argument("--" + field, default="")
    plan.add_argument("--cross-ks-reuse-summary", action="append", default=[], metavar="SUMMARY_JSON")
    plan.add_argument("--sequence-mode", choices=("cds", "protein"), default="cds")
    plan.add_argument("--duplication-taxon-holdout", type=int, choices=(0, 1), default=0)
    plan.add_argument("--duplication-coverage-audit", type=int, choices=(0, 1), default=0)
    plan.add_argument("--duplication-root-audit", type=int, choices=(0, 1), default=0)
    plan.add_argument("--duplication-copy-split-audit", type=int, choices=(0, 1), default=0)
    plan.add_argument("--duplication-family-panels", type=int, default=0)
    plan.add_argument("--cross-orientation-audit", type=int, choices=(0, 1), default=0)
    plan.add_argument("--cross-require-unique-phase", type=int, choices=(0, 1), default=0)
    plan.add_argument("--cross-require-both-orientations", type=int, choices=(0, 1), default=0)
    plan.add_argument("--ks-component-holdout", type=int, choices=(0, 1), default=0)
    plan.add_argument("--ks-sister-panels", type=int, choices=(0, 1), default=0)
    plan.add_argument("--ks-rooted-copies", type=int, choices=(0, 1), default=0)
    plan.add_argument("--cross-ks-method", choices=CROSS_DS_METHODS, default="yn00")
    for field, default in (("genetic-code", 1), ("count-bootstrap", 199), ("count-bootstrap-workers", 1), ("count-scan-workers", 1),
                           ("ks-bootstrap", 199), ("seed", 1),
                           ("max-pairs", 5000), ("max-ks-families", 500), ("max-ks-quartets", 0), ("max-ks-outgroup-families", 0),
                           ("max-terminal-ks-families", 0),
                           ("max-count-families", 10000),
                           ("diagonal-bound", 300), ("min-blocks", 3), ("min-cross-windows", 3), ("max-states", 256),
                           ("max-iterations", 200), ("multiplicity", 2)):
        plan.add_argument("--" + field, type=int, default=default)
    for field, default in (("cscore", 0.7), ("self-hit-percent", 98), ("minimum-mapping-fraction", 1),
                           ("alpha", 0.05), ("min-coverage", 0.2), ("min-duplication-score", 0.01)):
        plan.add_argument("--" + field, type=float, default=default)
    plan.add_argument("--outfile", required=True, type=Path)
    for name in ("contract", "verify", "run"):
        command = commands.add_parser(name)
        command.add_argument("--plan", required=True, type=Path)
        if name == "run":
            command.add_argument("--output", required=True, type=Path)
            command.add_argument("--cpus", type=int, default=1)
    node = commands.add_parser("classify")
    for field in ("gene-tree", "species-tree", "evidence", "output"):
        node.add_argument("--" + field, required=True, type=Path)
    node.add_argument("--family-id", required=True)
    node.add_argument("--species-parser", default="extended-taxonomic")
    node.add_argument("--species-regex", default="")
    node.add_argument("--species-map", default="")
    node.add_argument("--proximal-distance", type=int, default=10)
    node.add_argument("--native-tree-likelihood", type=int, choices=(0, 1), default=0)
    node.add_argument("--mul-diagnostics", type=int, choices=(0, 1), default=0)
    node.add_argument("--mul-h1", default="")
    node.add_argument("--mul-h2", default="")
    node.add_argument("--mul-max-candidates", type=int, default=10000)
    node.add_argument("--mul-max-state-pairs", type=int, default=10000000)
    node.add_argument("--mul-max-maps", type=int, default=100000)
    args = parser.parse_args()
    if args.command == "plan":
        write_json(args.outfile, make_plan(args))
    elif args.command == "classify":
        classify(args)
    else:
        payload = json.loads(args.plan.read_text())
        if args.command == "verify":
            verify(payload)
        elif args.command == "contract":
            sys.stdout.buffer.write(b"\0".join(value.encode() for value in contract(payload)) + b"\0")
        else:
            run_analysis(payload, args.output.resolve(), args.cpus)


if __name__ == "__main__":
    main()
