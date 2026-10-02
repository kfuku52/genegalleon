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
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    from pairwise_synteny import digest, prepare_genome, safe_token, source_file, write_json
    from pairwise_synteny_ds import DS_STATUSES, PAIR_COLUMNS, align_pair, ds_tool_identity, load_cds
    from wgd_evidence import branch_for_ks, combine_node, number, read_table, summarize_events, valid_position
except ImportError:
    from .pairwise_synteny import digest, prepare_genome, safe_token, source_file, write_json
    from .pairwise_synteny_ds import DS_STATUSES, PAIR_COLUMNS, align_pair, ds_tool_identity, load_cds
    from .wgd_evidence import branch_for_ks, combine_node, number, read_table, summarize_events, valid_position


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
    if (args.count_bootstrap < 0 or args.ks_bootstrap < 0 or args.seed < 0
            or any(parameters[key] < 1 for key in ("max_pairs", "max_ks_families", "max_count_families",
                                                  "diagonal_bound", "min_blocks", "max_iterations"))
            or not 0 < args.cscore <= 1 or not 0 < args.minimum_mapping_fraction <= 1
            or not 0 < args.self_hit_percent <= 100 or not 0 < args.alpha < 1
            or not 0 < args.min_coverage <= 1 or args.max_states < 16 or args.multiplicity < 2):
        raise ValueError("Invalid WGD analysis parameters")
    tools = owned_identity()
    inputs.update(tools["source_hashes"])
    inputs.update(tools["ds"]["source_hashes"])
    protect_output(args.outfile, inputs, workspace)
    return {"schema_version": 1, "workspace": str(workspace), "species_tree": str(tree_path),
            "counts": str(counts), "members": str(members), "genomes": genomes, "species": names,
            "parameters": parameters, "tools": tools, "input_hashes": inputs,
            "absent_inputs": ([str(genomes_path)] if not genomes_path.is_file() and not args.genomes else [])
            + ([str(code_file)] if not code_file.is_file() else [])}


def verify(plan):
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


def validated_members(plan):
    rows = read_table(plan["members"])
    counts = validated_counts(plan)
    if not rows:
        raise ValueError("Empty family-membership or gene-count table")
    key = "family_id" if "family_id" in rows[0] else "Orthogroup"
    if set(rows[0]) != {key, *plan["species"]}:
        raise ValueError("Family membership species must exactly match the full species tree")
    count_rows = dict(counts)
    seen_families, seen_genes = set(), set()
    for row in rows:
        family = safe_token(row[key], "family ID")
        if family in seen_families or family not in count_rows:
            raise ValueError("Duplicate or unmatched family-membership row")
        seen_families.add(family)
        for name in plan["species"]:
            genes = [gene.strip() for gene in row[name].split(",") if gene.strip()]
            if row[name].strip() and any(not gene.strip() for gene in row[name].split(",")):
                raise ValueError("Empty gene identifier in family membership")
            for gene in genes:
                safe_token(gene, "Gene identifier")
                if not gene.startswith(name + "_") or gene == name + "_":
                    raise ValueError("Family-membership gene identifier must have its species prefix")
            if len(genes) != len(set(genes)) or seen_genes.intersection(genes):
                raise ValueError("Gene identifiers must be unique across family membership")
            seen_genes.update(genes)
            count = count_rows[family][name]
            if count is not None and count != len(genes):
                raise ValueError("Family membership and gene count disagree")
    if seen_families != set(count_rows):
        raise ValueError("Family membership and gene counts have different families")
    return rows, key


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


def normalized_counts(plan, directory):
    names = plan["species"]
    tree = species_tree(plan["species_tree"])
    clades = [{leaf.name for leaf in child.leaves()} for child in tree.children]
    selected, audit = [], []
    for family, values in validated_counts(plan):
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
    if len(selected) < 2:
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


def estimate_ks(plan, output, sequences, anchors, cpus):
    rng = random.Random(plan["parameters"]["seed"])
    codes = {source["species"]: source["genetic_code"] for source in plan["genomes"]}
    requests, selection = {}, []
    for name in sorted(sequences):
        all_pairs = sorted({tuple(sorted((row["gene_a"], row["gene_b"]))) for row in anchors if row["species"] == name})
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
    members, family_key = validated_members(plan)
    counts = dict(validated_counts(plan))
    single_copy, audit = [], []
    for row in members:
        family = row[family_key]
        if any((counts[family][name] is not None and counts[family][name] > 1)
               or len([gene for gene in row[name].split(",") if gene.strip()]) > 1 for name in plan["species"]):
            audit.append({"family_id": family, "selection": "known_multicopy_family"})
            continue
        available = {}
        for name in sorted(sequences):
            genes = [gene.strip() for gene in row[name].split(",") if gene.strip()]
            if counts[family][name] == 1 and len(genes) == 1 and genes[0] in sequences[name]:
                available[name] = genes[0]
        if len(available) >= 2:
            single_copy.append((row[family_key], available))
        audit.append({"family_id": row[family_key], "selection": "eligible" if len(available) >= 2 else "insufficient_single_copy_cds_species"})
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
    reports, contrasts = {}, []
    for code in sorted(set(codes.values())):
        selected = [(pair, info) for pair, info in requests.items() if codes[info[0]] == code]
        if not selected:
            continue
        pairs_file, ds_file = output / f"aligned_pairs.code{code}.tsv", output / f"pair_ds.code{code}.tsv"
        with ThreadPoolExecutor(max_workers=cpus) as pool:
            aligned = pool.map(lambda item, code=code: align_pair(item[0], [sequences[item[1][0]], sequences[item[1][1]]], code), selected)
            write_table(pairs_file, aligned, PAIR_COLUMNS)
        run_tool(["cdskit", "dnds", "--pairs_file", pairs_file,
                  "--codon_table", str(code), "--threads", str(cpus), "--out_file", ds_file], output, f"dnds.code{code}")
        expected = {"|".join(pair) for pair, _ in selected}
        observed = set()
        for row in read_table(ds_file):
            pair_id = row.get("pair_id")
            if pair_id not in expected or pair_id in observed:
                raise ValueError("dS report contains an unrequested or duplicate Ks pair")
            if (not {"pair_id", "status", "dS"}.issubset(row) or row.get("status") not in DS_STATUSES
                    or (row["status"] == "ok" and number(row.get("dS")) is None)):
                raise ValueError("dS report contains an invalid Ks status or value")
            observed.add(pair_id)
            reports[pair_id] = row
        if observed != expected:
            raise ValueError("dS report does not cover all requested Ks pairs")
        for pair, info in selected:
            row = reports["|".join(pair)]
            selection.append({"pair_id": row["pair_id"], "species_a": info[0], "species_b": info[1],
                              "family_id": info[2], "pair_role": info[3], "status": row["status"], "dS": row["dS"]})
            if info[3] == "single_copy_family_contrast" and row["status"] == "ok" and number(row["dS"]) is not None:
                contrasts.append({"species_a": info[0], "species_b": info[1], "family_id": info[2], "ks": row["dS"]})
    write_table(output / "pair_selection.tsv", selection, ["pair_id", "species_a", "species_b", "family_id", "pair_role", "status", "dS"])
    write_table(output / "ks_family_selection.tsv", audit, ["family_id", "selection"])
    write_table(output / "family_ks_contrasts.tsv", contrasts, ["species_a", "species_b", "family_id", "ks"])
    return reports, contrasts


def plot_candidates(rows, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shown = sorted(rows, key=lambda row: float(row["likelihood_ratio_statistic"]), reverse=True)[:40]
    labels = []
    for row in reversed(shown):
        taxa = next(csv.reader([row["descendant_taxa"]]))
        label = f"B{row['branch_id']}: " + ", ".join(taxa[:2])
        if len(taxa) > 2:
            label += f" (+{len(taxa) - 2})"
        labels.append(label if len(label) <= 55 else label[:52] + "...")
    fig, ax = plt.subplots(figsize=(9, max(3, len(shown) * 0.25)))
    ax.barh(labels,
            [float(row["likelihood_ratio_statistic"]) for row in reversed(shown)],
            color=["#197c53" if row["event_support"] == "WGD-supported" else "#777777" for row in reversed(shown)])
    ax.set_xlabel("Conditional count-model likelihood-ratio statistic")
    prefix = f"Top {len(shown)} of {len(rows)} " if len(shown) < len(rows) else ""
    ax.set_title(prefix + "WGD candidates (green: combined evidence support)", fontsize=11)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(output / "wgd_candidates.pdf")
    plt.close(fig)


def run_analysis(plan, output, cpus):
    if cpus < 1:
        raise ValueError("cpus must be positive")
    verify(plan)
    protect_output(output, plan["input_hashes"], Path(plan["workspace"]))
    members, _ = validated_members(plan)
    validate_annotation_counts(plan, members)
    tree = species_tree(plan["species_tree"])
    output.mkdir(parents=True, exist_ok=False)
    normalized_counts(plan, output)
    p = plan["parameters"]
    run_tool([sys.executable, "-m", "nwkit", "wgd-count", "--infile", plan["species_tree"],
              "--input-rooted", "auto", "--counts", output / "counts.tsv", "--ascertainment", "root-clades",
              "--bootstrap", p["count_bootstrap"], "--seed", p["seed"], "--alpha", p["alpha"],
              "--max-states", p["max_states"], "--max-iterations", p["max_iterations"],
              "--multiplicity", p["multiplicity"], "--outfile", output / "count_candidates.tsv",
              "--model-out", output / "count_model.json"], output, "wgd-count")
    positions, sequences, summaries, anchors = prepare_sources(plan, output, cpus)
    reports, contrasts = estimate_ks(plan, output, sequences, anchors, cpus)
    bounds = {}
    if contrasts:
        run_tool([sys.executable, "-m", "nwkit", "ksrate", "--infile", plan["species_tree"], "--input-rooted", "auto",
                  "--ks-tsv", output / "family_ks_contrasts.tsv", "--bootstrap", p["ks_bootstrap"], "--seed", p["seed"],
                  "--ci-method", "pair-median-bonferroni",
                  "--outfile", output / "ks_boundaries.tsv", "--trios-out", output / "ks_trios.tsv",
                  "--model-out", output / "ks_model.json"], output, "ksrate")
        bounds = {(row["focal"], row["species_event_id"]): row for row in read_table(output / "ks_boundaries.tsv")}
    for row in anchors:
        report = reports.get(row["gene_a"] + "|" + row["gene_b"])
        row["ks"] = report["dS"] if report else "NA"
        row["ks_status"] = report["status"] if report else (
            "cds_unavailable" if any(gene not in sequences.get(row["species"], {}) for gene in (row["gene_a"], row["gene_b"]))
            else "not_sampled")
        branch, status = branch_for_ks(row["species"], number(row["ks"]), tree, bounds) if row["ks_status"] == "ok" else (None, "missing_ks")
        row.update(species_event_id=branch or "", placement_status=status)
    write_table(output / "anchor_evidence.tsv", anchors,
                ["species", "block_id", "gene_a", "gene_b", "ks", "ks_status", "species_event_id", "placement_status"])
    candidates = summarize_events(read_table(output / "count_candidates.tsv"), anchors, summaries,
                                  p["min_coverage"], p["min_blocks"], {row["gene_id"]: row for row in positions})
    write_table(output / "wgd_events.tsv", candidates, list(candidates[0]))
    plot_candidates(candidates, output)
    verify(plan)
    write_json(output / "summary.json", {"schema_version": 1, "plan": plan,
        "output_hashes": {path.name: digest(path) for path in sorted(output.iterdir())
                          if path.is_file() and path.suffix in {".tsv", ".json"}},
        "num_positioned_genes": len(positions),
        "num_raw_anchor_rows": len(anchors), "num_ks_contrasts": len(contrasts), "species_synteny": summaries,
        "ks_boundary_ci_methods": sorted({row.get("ci_method", "not_reported") for row in bounds.values()}),
        "num_WGD_supported_candidates": sum(row["event_support"] == "WGD-supported" for row in candidates),
        "interpretation": "Experimental one-event-at-a-time count scan plus branch-matched unquota synteny and Ks. Not a posterior, joint multi-event history, ploidy estimate, or validated allopolyploid model. Single-copy family contrasts may retain hidden paralogy. Missing synteny is not SSD evidence."})


def classify(args):
    from nwkit.clade_index import CladeIndex
    from nwkit.rooting_state import require_rooted
    from nwkit.util import _is_missing_support_value, assign_branch_ids, read_tree, write_tree

    evidence = args.evidence.resolve()
    # Ensure the candidate IDs refer to this full tree, not a family-pruned tree.
    summary = json.loads((evidence / "summary.json").read_text())
    if digest(args.species_tree) != summary["plan"]["input_hashes"][str(Path(summary["plan"]["species_tree"]))]:
        raise ValueError("WGD evidence and classification species trees differ")
    verify(summary["plan"])
    required = {"gene_positions.tsv", "anchor_evidence.tsv", "wgd_events.tsv"}
    if args.native_tree_likelihood:
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
    tree = read_tree(str(args.gene_tree), "auto", True, rooted="auto")
    require_rooted(tree, "WGD classification requires a rooted gene tree.")
    args.output.mkdir(parents=True, exist_ok=False)
    for node in tree.traverse():
        for prop in ("duplication_origin", "conditional_wgd_probability"):
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
    position_rows = read_table(evidence / "gene_positions.tsv")
    event_rows = read_table(evidence / "wgd_events.tsv")
    positions = {row["gene_id"]: row for row in position_rows}
    events = {row["species_event_id"]: row for row in event_rows}
    if len(positions) != len(position_rows) or len(events) != len(event_rows):
        raise ValueError("Duplicate gene or species event identifier in WGD evidence")
    anchors = defaultdict(list)
    for row in read_table(evidence / "anchor_evidence.tsv"):
        anchors[tuple(sorted((row["gene_a"], row["gene_b"])))].append(row)
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
        pairs = [(a.name, b.name) for a in node.children[0].leaves() for b in node.children[1].leaves()
                 if a.name in positions and b.name in positions and positions[a.name]["species"] == positions[b.name]["species"]]
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
    for field in ("species-tree", "counts", "members", "genomes"):
        plan.add_argument("--" + field, default="")
    plan.add_argument("--sequence-mode", choices=("cds", "protein"), default="cds")
    for field, default in (("genetic-code", 1), ("count-bootstrap", 199), ("ks-bootstrap", 199), ("seed", 1),
                           ("max-pairs", 5000), ("max-ks-families", 500), ("max-count-families", 10000),
                           ("diagonal-bound", 300), ("min-blocks", 3), ("max-states", 256),
                           ("max-iterations", 200), ("multiplicity", 2)):
        plan.add_argument("--" + field, type=int, default=default)
    for field, default in (("cscore", 0.7), ("self-hit-percent", 98), ("minimum-mapping-fraction", 1),
                           ("alpha", 0.05), ("min-coverage", 0.2)):
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
    node.add_argument("--species-parser", default="taxonomic")
    node.add_argument("--species-regex", default="")
    node.add_argument("--species-map", default="")
    node.add_argument("--proximal-distance", type=int, default=10)
    node.add_argument("--native-tree-likelihood", type=int, choices=(0, 1), default=0)
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
