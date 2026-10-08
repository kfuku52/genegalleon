#!/usr/bin/env python3
"""Recover unannotated CDS models from sparse multi-species and self synteny.

Independent CLI. Input files and the first BUSCO/reference selection are frozen;
no orthogroup, inferred species tree, or family-loss calls are required.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import importlib.util
import inspect
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote, unquote

from Bio import Phylo
from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    import busco_reference_quality as busco_reference_implementation
    from busco_reference_quality import COMPARABLE_QUALITY, busco_quality, patristic_distances
    from cds_model_normalisation import CdsModelNormaliser
    from fasta_sequence_store import exclusive_lock, fasta_records, open_text
    from gene_model_species_profiles import parameters_for, read_profiles
    from gff_attribute_syntax import validate_gff
    from input_generation_array_state import atomic_json, digest, digest_paths
    from pairwise_synteny import prepare_genome, safe_token, write_tsv
    from rescue_additional_candidates import nominate_genome_only_candidates, reassess_unanchored_models
    from rescue_anchor_admission import prepare_rescue_genome
    from rescue_coding_paths import (
        OwnershipIndex,
        alignment_support,
        coding_shape,
        original_ownership,
        resolve_coding_paths,
    )
    from rescue_model_quality import model_quality
    from rescue_model_store import export_legacy_models, iter_accepted_models, write_model_store
    from rescue_prediction_cache import (
        MINIPROT_MAX_SECONDARY,
        MINIPROT_OUTPUT_SCORE_RATIO,
        RAW_FIELDS,
        SEARCH_CONTRACT_FILE,
        frozen_prediction_cache_key,
        prediction_search_contract,
        verify_prediction_cache,
    )
    from rescue_prediction_cache import (
        stream_json_array as stream_json_array,
    )
    from rescue_raw_validation import RawValidationMemo
    from rescue_search_inputs import (
        LOCAL_RECORDS,
        export_inputs,
        genome_aliases,
        local_record,
        write_genome_records,
        write_record,
    )
    from rescue_terminal_completion import complete_terminals
    from species_labeling import extract_species_label
    from stage_output_hashes import hash_outputs, hash_paths
except ImportError:
    from . import busco_reference_quality as busco_reference_implementation
    from .busco_reference_quality import COMPARABLE_QUALITY, busco_quality, patristic_distances
    from .cds_model_normalisation import CdsModelNormaliser
    from .fasta_sequence_store import exclusive_lock, fasta_records, open_text
    from .gene_model_species_profiles import parameters_for, read_profiles
    from .gff_attribute_syntax import validate_gff
    from .input_generation_array_state import atomic_json, digest, digest_paths
    from .pairwise_synteny import prepare_genome, safe_token, write_tsv
    from .rescue_additional_candidates import nominate_genome_only_candidates, reassess_unanchored_models
    from .rescue_anchor_admission import prepare_rescue_genome
    from .rescue_coding_paths import (
        OwnershipIndex,
        alignment_support,
        coding_shape,
        original_ownership,
        resolve_coding_paths,
    )
    from .rescue_model_quality import model_quality
    from .rescue_model_store import export_legacy_models, iter_accepted_models, write_model_store
    from .rescue_prediction_cache import (
        MINIPROT_MAX_SECONDARY,
        MINIPROT_OUTPUT_SCORE_RATIO,
        RAW_FIELDS,
        SEARCH_CONTRACT_FILE,
        frozen_prediction_cache_key,
        prediction_search_contract,
        verify_prediction_cache,
    )
    from .rescue_prediction_cache import (
        stream_json_array as stream_json_array,
    )
    from .rescue_raw_validation import RawValidationMemo
    from .rescue_search_inputs import (
        LOCAL_RECORDS,
        export_inputs,
        genome_aliases,
        local_record,
        write_genome_records,
        write_record,
    )
    from .rescue_terminal_completion import complete_terminals
    from .species_labeling import extract_species_label
    from .stage_output_hashes import hash_outputs, hash_paths

SCHEMA = 1
PARAMETERS = ("common_references", "nearest_references", "minimum_busco", "cscore",
              "min_anchors", "distance", "diagonal_bound", "max_interval", "padding",
              "minimum_coverage", "minimum_identity", "max_intron", "genome_fallback",
              "max_genome_queries", "unanchored_min_species", "terminal_max_extension",
              "terminal_max_unaligned_c_overhang")


def run(command, directory, label, stdout=None):
    command = [str(x) for x in command]
    logs = directory / "logs"
    logs.mkdir(exist_ok=True)
    (logs / (label + ".command.json")).write_text(json.dumps(command) + "\n")
    with (logs / (label + ".log")).open("w") as log:
        with Path(stdout).open("w") if stdout else (logs / (label + ".stdout")).open("w") as out:
            result = subprocess.run(command, cwd=directory, stdout=out, stderr=log, check=False,
                                    env={**os.environ, "MPLBACKEND": "Agg"})
    if result.returncode:
        raise RuntimeError(f"{label} failed ({result.returncode}): {logs / (label + '.log')}")


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def package_source_hashes(name):
    """Bind support modules and compiled extensions, including editable sources."""
    spec = importlib.util.find_spec(name)
    hashes = {}
    for index, location in enumerate(spec.submodule_search_locations):
        root = Path(location)
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.suffix in {".py", ".so", ".pyd"}:
                hashes[f"{name}.files/{index}/{path.relative_to(root)}"] = digest(path)
    if not hashes:
        raise ValueError("No comparison package implementation found: " + name)
    return hashes


def identities():
    versions = {name: importlib.metadata.version(name) for name in
                ("nwkit", "kfFractBias", "jcvi", "biopython", "pysam", "numpy", "natsort", "more-itertools")}
    versions["python"] = sys.version
    for tool, arg in (("miniprot", "--version"), ("diamond", "version"), ("lastal", "--version"), ("lastdb", "--version")):
        versions[tool] = subprocess.check_output([tool, arg], text=True).strip()
        versions[tool + "_sha256"] = digest(shutil.which(tool))
    # Editable upstream sources can change without changing their version.
    modules = ("nwkit.sample", "kffractbias.io", "kffractbias.selfscan", "kffractbias.selfevidence", "jcvi.compara.catalog",
               "jcvi.compara.synteny", "jcvi.compara.blastfilter", "jcvi.apps.align")
    versions["source_hashes"] = {name: digest(importlib.util.find_spec(name).origin) for name in modules}
    for package in ("kffractbias", "jcvi"):
        versions["source_hashes"].update(package_source_hashes(package))
    versions["implementation"] = digest(__file__)
    versions["prediction_search_contract"] = prediction_search_contract()
    versions["attribute_syntax_implementation"] = digest(sys.modules[validate_gff.__module__].__file__)
    versions["ownership_attribute_parser_implementation"] = digest(Path(__file__).with_name("format_species_annotation") / "common.py")
    versions["mapping_implementation"] = digest(sys.modules[prepare_genome.__module__].__file__)
    versions["anchor_admission_implementation"] = digest(sys.modules[prepare_rescue_genome.__module__].__file__)
    versions["cds_normalisation_implementation"] = digest(sys.modules[CdsModelNormaliser.__module__].__file__)
    versions["reader_implementation"] = digest(sys.modules[fasta_records.__module__].__file__)
    versions["state_implementation"] = digest(sys.modules[atomic_json.__module__].__file__)
    versions["busco_quality_implementation"] = digest(busco_reference_implementation.__file__)
    versions["quality_implementation"] = digest(sys.modules[model_quality.__module__].__file__)
    versions["species_profiles_implementation"] = digest(sys.modules[read_profiles.__module__].__file__)
    versions["genome_index_implementation"] = digest(Path(__file__).with_name("gene_model_catalog.py"))
    versions["output_hashes_implementation"] = digest(sys.modules[hash_outputs.__module__].__file__)
    versions["search_inputs_implementation"] = digest(sys.modules[local_record.__module__].__file__)
    versions["model_store_implementation"] = digest(sys.modules[write_model_store.__module__].__file__)
    versions["raw_validation_implementation"] = digest(sys.modules[RawValidationMemo.__module__].__file__)
    for name, module in (("resolve_coding_paths", "rescue_coding_paths"),
                         ("complete_terminals", "rescue_terminal_completion"),
                         ("nominate_genome_only_candidates", "rescue_additional_candidates"),
                         ("verify_prediction_cache", "rescue_prediction_cache")):
        versions[name + "_implementation"] = digest(Path(__file__).with_name(module + ".py"))
    return versions


def guide_evidence(path, tree_path, expected_inputs=None, nearest_references=None):
    """Read and freeze the guide's exact inputs, outputs and panel diagnostics."""
    raw = path.read_bytes()
    payload = json.loads(raw)
    request = payload.get("request", {})
    if request.get("schema") != SCHEMA:
        raise ValueError("Invalid BUSCO guide receipt schema")
    files = payload.get("files", {})
    if (not isinstance(files, dict) or not {"guide_tree.nwk", "stability.json", "markers.json"} <= set(files)
            or any(Path(p).name != p for p in files) or (path.parent / "guide_tree.nwk").resolve() != tree_path.resolve()):
        raise ValueError("Invalid BUSCO guide receipt/tree path")
    sources = payload.get("request", {}).get("files", {})
    expected = {str(path.parent / p): h for p, h in files.items()}
    expected.update(sources)
    if digest_paths(expected) != expected or digest(path) != hashlib.sha256(raw).hexdigest():
        raise ValueError("Frozen BUSCO guide inputs/outputs changed")
    if expected_inputs is not None:
        guide_sources = payload["request"].get("sources", {})
        if set(guide_sources) != set(expected_inputs):
            raise ValueError("BUSCO guide species cohort differs from rescue; rebuild the guide for these species")
        for name, hashes in expected_inputs.items():
            source = guide_sources.get(name, {})
            if tuple(sources.get(source.get(key)) for key in ("cds", "short")) != hashes:
                raise ValueError("BUSCO guide does not use the rescue CDS/BUSCO inputs: " + name)
    stability = json.loads((path.parent / "stability.json").read_text())
    names = set(payload["request"].get("sources", {}))
    count = request.get("parameters", {}).get("nearest")
    if type(count) is not int or count < 0 or (nearest_references is not None and count != nearest_references):
        raise ValueError("BUSCO guide nearest-reference count differs from rescue; rebuild the guide with this count")
    count = min(count, len(names)-1)
    if not isinstance(stability, dict) or set(stability) != names:
        raise ValueError("Invalid BUSCO guide stability species")
    for name, evidence in stability.items():
        if not isinstance(evidence, dict):
            raise ValueError("Invalid BUSCO guide stability evidence: " + name)
        panels = [evidence.get(key) for key in ("nearest", "panel_0", "panel_1")]
        if (any(not isinstance(panel, list) or len(panel) != count
                or any(not isinstance(value, str) for value in panel) or len(panel) != len(set(panel))
                or not set(panel) <= names - {name} for panel in panels)
                or evidence.get("stable") is not (set(panels[0]) == set(panels[1]) == set(panels[2]))):
            raise ValueError("Invalid BUSCO guide stability evidence: " + name)
    expected[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
    return expected, stability


def build_plan(args):
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(root / ".plan.lock"):
        if (args.max_genome_queries < 0 or args.unanchored_min_species < 2
                or args.terminal_max_extension < 1):
            raise ValueError("Invalid bounded genome-search or terminal-completion parameters")
        species = sorted({extract_species_label(p.name) for p in args.cds_dir.iterdir()
                          if p.is_file() and p.name.removesuffix(".gz").endswith((".fa", ".fas", ".fasta", ".fna"))})
        if not species:
            raise ValueError("No CDS species found")
        tree_bytes = args.tree.read_bytes()
        tree = Phylo.read(io.StringIO(tree_bytes.decode("utf-8")), "newick")
        names = [tip.name for tip in tree.get_terminals()]
        if len(names) != len(set(names)) or not set(species) <= set(names):
            raise ValueError("Tree leaves must be unique and contain every CDS species")
        for clade in tree.find_clades():
            if clade.branch_length is not None and (not math.isfinite(clade.branch_length) or clade.branch_length < 0):
                raise ValueError("Tree branch lengths must be finite and nonnegative")
        topology_only = not any((c.branch_length or 0) > 0 for c in tree.find_clades() if c is not tree.root)
        if not topology_only and any(c.branch_length is None for c in tree.find_clades() if c is not tree.root):
            raise ValueError("Partial branch lengths: provide lengths on every non-root edge or a topology-only tree")
        if topology_only:
            for clade in tree.find_clades():
                clade.branch_length = 1.0 if clade is not tree.root else 0.0
        for name in set(names) - set(species):
            tree.prune(name)
        codes = {}
        code_hash = None
        if args.genetic_codes:
            code_hash = digest(args.genetic_codes)
            with args.genetic_codes.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if len(reader.fieldnames or ()) != 2 or set(reader.fieldnames or ()) != {"species", "genetic_code"}:
                    raise ValueError("Genetic-code table requires species and genetic_code columns")
                rows = list(reader)
            for row in rows:
                if None in row or any(value is None or not value.strip() for value in row.values()):
                    raise ValueError("Invalid genetic-code table row")
                row = {k: v.strip() for k, v in row.items()}
                if row["species"] in codes:
                    raise ValueError("Duplicate genetic-code species")
                codes[row["species"]] = int(row["genetic_code"])
                if codes[row["species"]] not in CodonTable.unambiguous_dna_by_id:
                    raise ValueError("Unknown genetic code: " + row["species"])
        sources = {}
        for name in species:
            safe_token(name, "species")
            def discover(directory, suffixes, name=name):
                matches = [p for p in directory.iterdir() if p.is_file() and extract_species_label(p.name) == name
                           and p.name.removesuffix(".gz").endswith(suffixes)]
                if len(matches) != 1:
                    raise ValueError(f"Expected one {name} source in {directory}, found {len(matches)}")
                return str(matches[0].resolve())
            cds = discover(args.cds_dir, (".fa", ".fas", ".fasta", ".fna"))
            gff = discover(args.gff_dir, (".gff", ".gff3"))
            validate_gff(gff)
            genome = discover(args.genome_dir, (".fa", ".fas", ".fasta", ".fna"))
            busco = discover(args.busco_dir, (".busco.short.txt",))
            code = codes.get(name, args.genetic_code)
            if code not in CodonTable.unambiguous_dna_by_id:
                raise ValueError(f"Unknown genetic code for {name}")
            sources[name] = {"species": name, "fasta": cds, "gff": gff, "genome": genome,
                             "busco": busco, "quality": busco_quality(busco), "mode": "cds",
                             "feature": args.feature, "attribute": args.attribute, "genetic_code": code}
        comparable = {tuple(s["quality"][k] for k in COMPARABLE_QUALITY) for s in sources.values()}
        if len(comparable) != 1:
            raise ValueError("Reference BUSCO scores must use the same lineage, version, mode, dataset date and marker count")
        eligible = {n for n, s in sources.items() if s["quality"]["complete_pct"] >= args.minimum_busco}
        if len(eligible) < args.common_references:
            raise ValueError(f"Need {args.common_references} eligible common references; found {len(eligible)}")
        if args.gemoma_jar and any(s["genetic_code"] != 1 for s in sources.values()):
            raise ValueError("Optional GeMoMa refinement currently requires genetic code 1")
        files = [str(args.tree.resolve())] + [s[k] for s in sources.values() for k in ("fasta", "gff", "genome", "busco")]
        if args.genetic_codes:
            files.append(str(args.genetic_codes.resolve()))
        params = {key: getattr(args, key) for key in PARAMETERS}
        profiles_path = getattr(args, "species_profiles", None)
        profiles = read_profiles(profiles_path, sources)
        if profiles_path:
            files.append(str(profiles_path.resolve()))
        request = {"schema": SCHEMA, "sources": sources, "files": digest_paths(files), "parameters": params,
                   "species_profiles": profiles,
                   "tools": identities(), "tree_metric": "unit_edges" if topology_only else "patristic_distance",
                   "gemoma_jar": str(args.gemoma_jar.resolve()) if args.gemoma_jar else None,
                   "gemoma_java": None}
        request["output_storage"] = {"format": getattr(args, "model_storage", "compact"),
                                     "retain_search_inputs": bool(getattr(args, "retain_search_inputs", False))}
        prediction_cache = getattr(args, "prediction_cache", None)
        if prediction_cache:
            request["prediction_cache"] = frozen_prediction_cache_key(prediction_cache, sources)
        guide_stability = None
        if getattr(args, "guide_tree_receipt", None):
            if topology_only:
                raise ValueError("A BUSCO k-mer guide must have informative branch lengths")
            guide_files, guide_stability = guide_evidence(
                args.guide_tree_receipt.resolve(), args.tree,
                {n: (request["files"][sources[n]["fasta"]], request["files"][sources[n]["busco"]]) for n in species},
                nearest_references=args.nearest_references)
            if not set(species) <= set(guide_stability):
                raise ValueError("BUSCO guide diagnostics lack rescue species")
            request["files"].update(guide_files)
            request["guide_tree_receipt"] = str(args.guide_tree_receipt.resolve())
        if request["files"][str(args.tree.resolve())] != hashlib.sha256(tree_bytes).hexdigest():
            raise ValueError("Initial tree changed during planning")
        if args.genetic_codes and request["files"][str(args.genetic_codes.resolve())] != code_hash:
            raise ValueError("Genetic-code table changed during planning")
        if any(busco_quality(s["busco"]) != s["quality"] for s in sources.values()):
            raise ValueError("BUSCO quality changed during planning")
        if args.gemoma_jar:
            request["files"][request["gemoma_jar"]] = digest(args.gemoma_jar)
            executable = shutil.which(args.gemoma_java)
            if not executable:
                raise ValueError("GeMoMa Java executable not found: " + args.gemoma_java)
            request["gemoma_java"] = str(Path(executable).resolve())
            request["files"][request["gemoma_java"]] = digest(request["gemoma_java"])
            request["gemoma_java_version"] = subprocess.check_output(
                [request["gemoma_java"], "-version"], stderr=subprocess.STDOUT, text=True).strip()
        path = root / "plan.json"
        if path.exists():
            old = json.loads(path.read_text())
            if old["request"] != request:
                raise ValueError("Frozen plan differs from current inputs/settings/tools. Use a new output directory.")
            return old
        with tempfile.TemporaryDirectory(prefix=".selection-", dir=root) as scratch:
            tmp = Path(scratch)
            Phylo.write(tree, tmp / "selection.nwk", "newick", format_branch_length="%1.17g")
            write_tsv(tmp / "quality.tsv", ("leaf_name", "busco_complete_pct"),
                      [(n, sources[n]["quality"]["complete_pct"]) for n in species])
            run(["nwkit", "sample", "--infile", tmp / "selection.nwk", "--outfile", tmp / "references.nwk",
                 "--n", args.common_references, "--method", "max-pd", "--trait", tmp / "quality.tsv",
                 "--filter", f"busco_complete_pct:ge:{args.minimum_busco}", "--rank", "busco_complete_pct:desc",
                 "--report", tmp / "references.tsv"], tmp, "select")
            refs = [row["leaf_name"] for row in table(tmp / "references.tsv")]
            if len(refs) != args.common_references or len(set(refs)) != len(refs) or not set(refs) <= eligible:
                raise ValueError("NWKIT did not select the requested eligible reference set")
            distances = patristic_distances(tree)
            neighbors = {n: sorted((m for m in species if m != n),
                                   key=lambda m: (distances[n][m], -sources[m]["quality"]["complete_pct"], m))[:args.nearest_references]
                         for n in species}
            if guide_stability:
                # Retain alternative supported neighbourhoods within a bounded
                # extra budget; marker-panel agreement is not a bootstrap value.
                for name in species:
                    evidence = guide_stability[name]
                    if not evidence["stable"]:
                        alternatives = (set(evidence["panel_0"]) | set(evidence["panel_1"])) & set(species) - {name}
                        ordered = sorted(alternatives - set(neighbors[name]),
                                         key=lambda m: (distances[name][m], -sources[m]["quality"]["complete_pct"], m))
                        neighbors[name] += ordered[:args.nearest_references]
            donors = {n: sorted((set(refs) | set(neighbors[n])) - {n}) for n in species}
            pairs = sorted({tuple(sorted((n, m))) for n in species for m in donors[n]})
            jobs = [{"a": a, "b": b, "kind": "pair"} for a, b in pairs]
            jobs += [{"a": n, "b": n, "kind": "self"} for n in species]
            plan = {"request": request, "common_references": refs, "nearest_references": neighbors,
                    "donors": donors, "species": species,
                    "synteny_jobs": [{**job, "index": i, "id": f"comparison_{i:06d}"} for i, job in enumerate(jobs, 1)]}
            if guide_stability:
                plan["guide_tree_stability"] = {n: guide_stability[n] for n in species}
            if guide_stability and digest_paths(request["files"]) != request["files"]:
                raise ValueError("Guide evidence changed during reference selection")
            for file in ("selection.nwk", "quality.tsv", "references.nwk", "references.tsv"):
                shutil.copyfile(tmp / file, root / file)
            atomic_json(path, plan, immutable=True)
        return plan


def load(root, *, inspect_only=False):
    plan = json.loads((root / "plan.json").read_text())
    if (plan["request"]["schema"] != SCHEMA
            or (not inspect_only and plan["request"]["tools"] != identities())):
        raise ValueError("Rescue schema/tools changed; use a new output directory")
    plan_digest(root, plan)
    if plan["request"].get("guide_tree_receipt"):
        path = Path(plan["request"]["guide_tree_receipt"])
        if digest(path) != plan["request"]["files"][str(path)]:
            raise ValueError("Frozen BUSCO guide receipt changed")
    return plan


def plan_digest(root, plan):
    """Refuse to stamp a stale in-memory plan with the current file's hash."""
    path = root / "plan.json"
    raw = path.read_bytes()
    result = hashlib.sha256(raw).hexdigest()
    if json.loads(raw) != plan or digest(path) != result:
        raise ValueError("Frozen rescue plan changed during execution")
    return result


def verify_sources(plan, names, keys):
    paths = [plan["request"]["sources"][n][k] for n in names for k in keys]
    current = digest_paths(paths)
    if any(plan["request"]["files"][p] != value for p, value in current.items()):
        raise ValueError("Frozen rescue input changed")


def verified(directory, key, hash_workers=1):
    try:
        receipt = json.loads((directory / "receipt.json").read_text())
        if not (isinstance(receipt, dict) and receipt.get("key") == key
                and isinstance(receipt.get("files"), dict) and bool(receipt["files"])
                and all(isinstance(p, str) and isinstance(value, str)
                        for p, value in receipt["files"].items())):
            return False
        paths = [directory / p for p in receipt["files"]]
        # digest itself refuses nonregular files and fences content reads with
        # before/after stat identities. Every listed file is still read in full.
        return hash_paths(directory, paths, workers=hash_workers, hash_function=digest) == receipt["files"]
    except (OSError, ValueError):
        return False


def require_same_key(expected, current):
    if expected != current:
        raise ValueError("Stage dependencies changed during execution; retry after they are complete")


def recover_publication(dest, journal, token):
    """Recover a replacement interrupted by SIGKILL, under the job's lock."""
    if not journal.exists():
        return
    state = json.loads(journal.read_text())
    if (not isinstance(state, dict) or set(state) != {"destination", "backup", "temporary", "key"}
            or state["destination"] != dest.name):
        raise ValueError("Invalid publication journal: " + str(journal))
    paths = {}
    for key, prefix in (("backup", ".previous-"), ("temporary", ".working-")):
        name = state[key]
        if not isinstance(name, str) or Path(name).name != name or not name.startswith(prefix + token + "-"):
            raise ValueError("Unsafe publication journal path: " + str(journal))
        paths[key] = dest.parent / name
        if paths[key].is_symlink():
            raise ValueError("Publication directory must not be a symlink: " + str(paths[key]))
    backup, temporary = paths["backup"], paths["temporary"]
    if backup.exists():
        if dest.exists():
            if verified(dest, state["key"]):
                shutil.rmtree(backup)
            else:
                quarantine = Path(tempfile.mkdtemp(prefix=".interrupted-", dir=dest.parent))
                quarantine.rmdir()
                dest.rename(quarantine)
                backup.rename(dest)
        else:
            backup.rename(dest)
    if temporary.exists():
        failed = dest.parent / (dest.name + ".failed")
        if failed.exists():
            shutil.rmtree(failed)
        temporary.rename(failed)
    journal.unlink()


def stage(root, relative, key, builder, guard=None, hash_workers=1):
    """Publish only complete jobs; retain their diagnostics after failed attempts."""
    dest = root / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    lock = root / ".locks" / (str(relative).replace("/", "__") + ".lock")
    lock.parent.mkdir(exist_ok=True)
    journal = lock.with_suffix(".publish.json")
    token = hashlib.sha256(str(relative).encode()).hexdigest()[:12]
    with exclusive_lock(lock):
        recover_publication(dest, journal, token)
        if guard:
            guard()
        if verified(dest, key, hash_workers=hash_workers):
            print(f"Reused {relative}", flush=True)
            return dest
        tmp = Path(tempfile.mkdtemp(prefix=".working-" + token + "-", dir=dest.parent))
        try:
            builder(tmp)
            if guard:
                guard()
            files = hash_outputs(tmp, workers=hash_workers, hash_function=digest)
            if not files:
                raise ValueError("Empty stage outputs")
            atomic_json(tmp / "receipt.json", {"key": key, "files": files})
            # Retain the previous complete directory until replacement succeeds.
            backup = None
            if dest.exists():
                backup = Path(tempfile.mkdtemp(prefix=".previous-" + token + "-", dir=dest.parent))
                backup.rmdir()
                atomic_json(journal, {"destination": dest.name, "backup": backup.name, "temporary": tmp.name, "key": key})
                dest.rename(backup)
            try:
                tmp.rename(dest)
            except BaseException:
                if backup is not None:
                    backup.rename(dest)
                raise
            if backup is not None:
                shutil.rmtree(backup)
                journal.unlink()
        except BaseException:
            failed = dest.parent / (dest.name + ".failed")
            if tmp.exists():
                if failed.exists():
                    shutil.rmtree(failed)
                tmp.rename(failed)
            raise
    print(f"Completed {relative}", flush=True)
    return dest


def prepared(root, plan, name):
    source = plan["request"]["sources"][name]
    key = {"plan": plan_digest(root, plan), "species": name}
    def guard():
        require_same_key(key["plan"], plan_digest(root, plan))
        verify_sources(plan, [name], ["fasta", "gff", "genome"])
    guard()
    def build(tmp):
        genes, meta = prepare_rescue_genome(source, tmp, "genes", 1.0)
        atomic_json(tmp / "mapping.json", meta)
        atomic_json(tmp / "positions.json", [vars(g) for g in genes])
    return stage(root, Path("prepared") / name, key, build, guard)


def parse_anchors(path):
    blocks, current = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            if current:
                blocks.append(current)
            current = []
        elif line.strip():
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("Invalid anchor row")
            current.append(fields[:2])
    if current:
        blocks.append(current)
    return blocks


def align_self(cpus, cscore):
    """Retain high-identity paralogs; remove true same-ID self hits only."""
    from jcvi.apps.align import last
    from jcvi.compara.blastfilter import main as filter_synteny_hits

    last(["self.pep", "self.pep", f"--cpus={cpus}"], "prot")
    with Path("self.self.last").open() as handle, Path("self.self.last.noself").open("w") as out:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 12:
                raise ValueError("Malformed LAST protein alignment")
            if fields[0] != fields[1]:
                out.write(line)
    filter_synteny_hits(["self.self.last.noself", f"--cscore={cscore}", "--tandem_Nmax=0", "--no_strip_names"])


def build_comparison(tmp, job, dirs, params, cpus):
    """Only the BED/protein files, comparison options and toolchain affect this stage."""
    for name, side in ((job["a"], "target"), (job["b"], "query")) if job["kind"] == "pair" else ((job["a"], "self"),):
        for ext in ("bed", "pep"):
            shutil.copyfile(dirs[name] / ("genes." + ext), tmp / (side + "." + ext))
    if job["kind"] == "pair":
        run([sys.executable, "-m", "jcvi.compara.catalog", "ortholog", "target", "query",
             "--dbtype=prot", "--align_soft=diamond_blastp", "--no_strip_names", "--no_dotplot",
             "--tandem_Nmax=0", "--ignore_zero_anchor",
             f"--cpus={cpus}", f"--cscore={params['cscore']}", f"--min_size={params['min_anchors']}",
             f"--dist={params['distance']}"], tmp, "pair")
        anchors = tmp / "target.query.lifted.anchors"
        if not anchors.exists():
            raw = tmp / "target.query.anchors"
            if not raw.is_file() or raw.read_text().strip():
                raise ValueError("Pair scan did not produce a complete anchor result")
            anchors.write_text("")
    else:
        run([sys.executable, Path(__file__).resolve(), "self-align", "--cpus", cpus,
             "--cscore", params["cscore"]], tmp, "self_align")
        filtered = tmp / "self.self.last.noself.filtered"
        run([sys.executable, "-m", "kffractbias.selfscan", "scan", filtered, filtered.with_suffix(""),
             tmp / "self.bed", tmp / "unquota.anchors", f"--diagonal-bound={params['diagonal_bound']}",
             "--screening=none", "--allow-empty"], tmp, "self_scan")
        anchors = tmp / "self.self.lifted.anchors"
    atomic_json(tmp / "blocks.json", parse_anchors(anchors))


def comparison_cache_key(root, plan, job):
    names = sorted({job["a"], job["b"]})
    tools = plan["request"]["tools"]
    owners = ("kfFractBias", "jcvi", "biopython", "numpy", "natsort", "more-itertools", "python",
              "diamond", "diamond_sha256", "lastal", "lastal_sha256", "lastdb", "lastdb_sha256")
    modules = {k: v for k, v in tools["source_hashes"].items() if k.startswith(("kffractbias.", "jcvi."))}
    algorithm = "\n".join(inspect.getsource(f) for f in (build_comparison, align_self, parse_anchors, run))
    return {"schema": 1, "job": {k: job[k] for k in ("a", "b", "kind")},
            "inputs": {n: {ext: digest(root / "prepared" / n / ("genes." + ext)) for ext in ("bed", "pep")} for n in names},
            "parameters": {k: plan["request"]["parameters"][k] for k in ("cscore", "min_anchors", "distance", "diagonal_bound")},
            "tools": {k: tools[k] for k in owners}, "source_hashes": modules,
            "algorithm_sha256": hashlib.sha256(algorithm.encode()).hexdigest()}


def copy_verified_comparison(source, destination, key):
    """Copy hashed files, never link writable outputs to another plan's cache."""
    receipt_hash = digest(source / "receipt.json")
    if not verified(source, key):
        raise ValueError("Comparison cache incomplete or corrupted")
    receipt = json.loads((source / "receipt.json").read_text())
    for relative, expected in receipt["files"].items():
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Unsafe comparison cache path")
        target = destination / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / path, target)
        if digest(target) != expected:
            raise ValueError("Comparison cache changed during copying")
    if digest(source / "receipt.json") != receipt_hash or not verified(source, key):
        raise ValueError("Comparison cache changed during copying")


def synteny(root, plan, index, cpus, comparison_cache=None):
    jobs = plan["synteny_jobs"]
    if not 1 <= index <= len(jobs):
        raise ValueError("Synteny index outside frozen plan")
    job = jobs[index - 1]
    dirs = {n: prepared(root, plan, n) for n in {job["a"], job["b"]}}
    params = plan["request"]["parameters"]
    key = comparison_key(root, job)
    cache_key = comparison_cache_key(root, plan, job)
    cache_root = (comparison_cache or root.parent / "gene_model_rescue_comparison_cache").resolve()
    if cache_root == root or root in cache_root.parents:
        raise ValueError("Comparison cache must be outside the frozen rescue output")
    cache_id = hashlib.sha256(json.dumps(cache_key, sort_keys=True).encode()).hexdigest()
    def guard():
        plan_digest(root, plan)
        require_same_key(key, comparison_key(root, job))
        require_same_key(cache_key, comparison_cache_key(root, plan, job))
    def build(tmp):
        cached = stage(cache_root, Path("comparisons") / cache_id, cache_key,
                       lambda output: build_comparison(output, job, dirs, params, cpus), guard)
        copy_verified_comparison(cached, tmp, cache_key)
        atomic_json(tmp / "job.json", job)
        atomic_json(tmp / "cache.json", {"cache_key": cache_key, "cache_receipt_sha256": digest(cached / "receipt.json")})
    return stage(root, Path("synteny") / job["id"], key, build, guard)


def comparison_key(root, job):
    plan_hash = digest(root / "plan.json")
    names = sorted({job["a"], job["b"]})
    for name in names:
        if not verified(root / "prepared" / name, {"plan": plan_hash, "species": name}):
            raise ValueError("Prepared annotation incomplete or corrupted: " + name)
    return {"plan": plan_hash, "job": job,
            "prepared": {n: digest(root / "prepared" / n / "receipt.json") for n in names}}


def local_flanks(left, right):
    """Match parallel donor paths in rank order, retaining all equal-cost ties."""
    if len(left) == len(right):
        return list(zip(left, right, strict=True))
    swapped = len(left) > len(right)
    small, large = (right, left) if swapped else (left, right)
    costs = [[0] * (len(large) + 1)] + [[math.inf] * (len(large) + 1) for _ in small]
    for i in range(1, len(small) + 1):
        for j in range(1, len(large) + 1):
            costs[i][j] = min(costs[i][j - 1], costs[i - 1][j - 1] + abs(small[i - 1][0] - large[j - 1][0]))
    pending, visited, matches = [(len(small), len(large))], set(), set()
    while pending:
        i, j = pending.pop()
        if not i or not j or (i, j) in visited:
            continue
        visited.add((i, j))
        if costs[i][j] == costs[i][j - 1]:
            pending.append((i, j - 1))
        if costs[i][j] == costs[i - 1][j - 1] + abs(small[i - 1][0] - large[j - 1][0]):
            pair = (small[i - 1], large[j - 1])
            matches.add(tuple(reversed(pair)) if swapped else pair)
            pending.append((i - 1, j - 1))
    return sorted(matches)


def candidate_intervals(target, donor, blocks, target_positions, donor_positions, params):
    """Two flanking anchors per block, independent of hits in other WGD copies."""
    by_chr = defaultdict(list)
    for gene in donor_positions.values():
        by_chr[gene["seqid"]].append(gene)
    rank = {}
    for genes in by_chr.values():
        genes.sort(key=lambda g: (g["start"], g["end"], g["gene_id"]))
        rank.update({g["gene_id"]: i for i, g in enumerate(genes)})
    candidates = []
    for block_index, block in enumerate(blocks, 1):
        if any(a not in target_positions or b not in donor_positions for a, b in block):
            raise ValueError("Anchor ID missing from prepared annotation")
        chromosomes = defaultdict(lambda: defaultdict(set))
        for a, b in block:
            chromosomes[(target_positions[a]["seqid"], donor_positions[b]["seqid"])][a].add(b)
        flanks = []
        for _, groups in sorted(chromosomes.items()):
            ordered = sorted(groups, key=lambda a: (target_positions[a]["start"], target_positions[a]["end"], a))
            for a, c in zip(ordered, ordered[1:], strict=False):
                left = sorted((rank[b], b) for b in groups[a])
                right = sorted((rank[d], d) for d in groups[c])
                # Lifted blocks can merge parallel WGD copies. Follow adjacent
                # target anchors and parallel donor paths, rather than
                # interleaving copies by donor rank alone.
                flanks.extend(((a, b), (c, d)) for (_, b), (_, d) in local_flanks(left, right))
        for (a, b), (c, d) in flanks:
            left, right = target_positions[a], target_positions[c]
            first, last = donor_positions[b], donor_positions[d]
            if first["seqid"] != last["seqid"] or left["seqid"] != right["seqid"] or a == c:
                continue
            lo, hi = sorted((left, right), key=lambda g: g["start"])
            start, end = lo["end"], hi["start"]
            if not 0 < end - start <= params["max_interval"]:
                continue
            lower, upper = sorted((rank[b], rank[d]))
            for gene in by_chr[first["seqid"]][lower + 1:upper]:
                if target == donor and gene["seqid"] == lo["seqid"] and gene["start"] < end and gene["end"] > start:
                    continue
                value = {"target": target, "donor": donor, "query": gene["gene_id"], "seqid": lo["seqid"],
                         "start": max(0, start - params["padding"]), "end": end + params["padding"],
                         "expected_start": start, "expected_end": end,
                         "block": block_index, "left_anchor": lo["gene_id"], "right_anchor": hi["gene_id"],
                         "donor_left_anchor": b, "donor_right_anchor": d,
                         "orientation": "+" if rank[b] < rank[d] else "-"}
                value["expected_strand"] = gene["strand"] if value["orientation"] == "+" else {"+": "-", "-": "+"}.get(gene["strand"], ".")
                value["id"] = "region_" + hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:20]
                candidates.append(value)
    return candidates


def candidates(root, plan, name):
    params = parameters_for(plan["request"], name)
    position = {}
    def positions(n):
        if n not in position:
            position[n] = {g["gene_id"]: g for g in json.loads((root / "prepared" / n / "positions.json").read_text())}
        return position[n]
    result = []
    for job in plan["synteny_jobs"]:
        a, b = job["a"], job["b"]
        if name not in {a, b} or (a != b and (b if name == a else a) not in plan["donors"][name]):
            continue
        directory = root / "synteny" / job["id"]
        if not verified(directory, comparison_key(root, job)):
            raise ValueError("Comparison incomplete or corrupted: " + job["id"])
        blocks = json.loads((directory / "blocks.json").read_text())
        for target, donor, reverse in ((a, b, False), (b, a, True)):
            if target != name:
                continue
            oriented = [[[y, x] for x, y in block] for block in blocks] if reverse else blocks
            found = candidate_intervals(target, donor, oriented, positions(target), positions(donor), params)
            for row in found:
                row["comparison"] = job["id"]
            result.extend(found)
    return list({r["id"]: r for r in result}.values())


def rescue_key(root, plan, name):
    jobs = [j for j in plan["synteny_jobs"] if name in {j["a"], j["b"]}
            and (j["a"] == j["b"] or (j["b"] if name == j["a"] else j["a"]) in plan["donors"][name])]
    plan_hash = plan_digest(root, plan)
    for job in jobs:
        if not verified(root / "synteny" / job["id"], comparison_key(root, job)):
            raise ValueError("Comparison incomplete or corrupted: " + job["id"])
    donors = [name, *plan["donors"][name]]
    for donor in donors:
        if not verified(root / "prepared" / donor, {"plan": plan_hash, "species": donor}):
            raise ValueError("Prepared annotation incomplete or corrupted: " + donor)
    cache = plan["request"].get("prediction_cache")
    if cache:
        expected = {**cache, "species": {name: cache["species"][name]} if name in cache["species"] else {}}
        require_same_key(expected, frozen_prediction_cache_key(cache["root"], [name]))
    return {"plan": plan_hash, "species": name,
            "comparisons": {j["id"]: digest(root / "synteny" / j["id"] / "receipt.json") for j in jobs},
            "prepared": {n: digest(root / "prepared" / n / "receipt.json") for n in donors}}


def attributes(text):
    return {k: unquote(v) for field in text.split(";") if "=" in field for k, v in [field.split("=", 1)]}


def read_miniprot(path):
    return list(iter_miniprot(miniprot_lines(path)))


def iter_miniprot(lines):
    """Yield completed predictions without retaining another expanded array."""
    current = None
    paf = None
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith("##PAF\t"):
            if current is not None:
                yield current
                current = None
            f = line.split("\t")[1:]
            if len(f) < 12:
                raise ValueError("Invalid embedded miniprot PAF")
            if f[5] == "*":  # -u emits PAF records for unmapped proteins too.
                paf = None
                continue
            tags = {v[:2]: v[5:] for v in f[12:]}
            cigar = tags.get("cg", "")
            if not re.fullmatch(r"(?:[1-9]\d*[MIDFGNUV])+", cigar):
                raise ValueError("Mapped miniprot record lacks valid protein CIGAR")
            operations = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDFGNUV])", cigar)]
            span = int(f[3]) - int(f[2])
            consumed = sum(n if op in "MI" else 1 if op in "GUV" else 0 for n, op in operations)
            aligned = sum(n if op == "M" else 1 if op in "GUV" else 0 for n, op in operations)
            if not 0 <= int(f[2]) < int(f[3]) <= int(f[1]) or consumed != span:
                raise ValueError("Inconsistent miniprot query coordinates and CIGAR")
            paf = {"coverage": aligned / int(f[1]), "query_span_coverage": span / int(f[1]), "query_start": int(f[2]),
                   "query_end": int(f[3]), "query_length": int(f[1]),
                   "frameshift": any(op in "FG" for _, op in operations) or int(tags.get("fs", 0)) > 0,
                   "paf_query": f[0], "paf_seqid": f[5], "paf_strand": f[4], "paf": line}
        elif not line.startswith("#") and line.strip():
            f = line.split("\t")
            if len(f) != 9:
                raise ValueError("Invalid miniprot GFF row")
            attr = attributes(f[8])
            if f[2] == "mRNA":
                if paf is None:
                    raise ValueError("miniprot GFF requires embedded PAF records")
                current = {**paf, "query": attr["Target"].split()[0], "seqid": f[0], "strand": f[6],
                           "identity": float(attr.get("Identity", 0)), "cds": [], "id": attr["ID"]}
                if (current["query"], f[0], f[6]) != (paf["paf_query"], paf["paf_seqid"], paf["paf_strand"]):
                    raise ValueError("miniprot GFF and PAF do not describe the same model")
                paf = None
            elif f[2] == "CDS":
                if current is None or attr.get("Parent") != current["id"]:
                    raise ValueError("Unbound miniprot CDS")
                if (f[0], f[6]) != (current["seqid"], current["strand"]):
                    raise ValueError("miniprot CDS has inconsistent contig/strand")
                current["cds"].append([int(f[3]) - 1, int(f[4]), int(f[7])])
    if current is not None:
        yield current


def miniprot_lines(path):
    """Stream large predictor output without a second full text/list copy."""
    with Path(path).open() as handle:
        for line in handle:
            yield line.rstrip("\r\n")


def validate_model(model, genome, code, params):
    """Check assembled CDS in transcription order without repairing disruptions."""
    problems = []
    if model["strand"] not in {"+", "-"}:
        raise ValueError("Predicted CDS has unknown strand")
    exons = sorted([list(e) for e in model["cds"]], reverse=model["strand"] == "-")
    if not exons:
        return {**model, "problems": ["no_cds"], "sequence": ""}
    if model["frameshift"]:
        problems.append("frameshift")
    if any(not math.isfinite(model[key]) or not 0 <= model[key] <= 1 for key in ("coverage", "identity")):
        raise ValueError("Invalid alignment coverage or identity")
    if model["coverage"] < params["minimum_coverage"]:
        problems.append("low_coverage")
    if model["identity"] < params["minimum_identity"]:
        problems.append("low_identity")
    seq = ""
    cumulative = 0
    for i, (start, end, phase) in enumerate(exons):
        if start < 0 or end > genome.get_reference_length(model["seqid"]) or end <= start:
            raise ValueError("Predicted CDS outside genome")
        if phase != (3 - cumulative % 3) % 3:
            problems.append("invalid_phase")
        piece = genome.fetch(model["seqid"], start, end).upper()
        seq += str(Seq(piece).reverse_complement()) if model["strand"] == "-" else piece
        cumulative += end - start
        if i:
            prev = exons[i - 1]
            intron_start, intron_end = (prev[1], start) if model["strand"] == "+" else (end, prev[0])
            if intron_end - intron_start < 4:
                problems.append("invalid_splice")
            else:
                if intron_end - intron_start > params.get("max_intron", math.inf):
                    problems.append("intron_too_long")
                donor = genome.fetch(model["seqid"], intron_start, intron_start + 2).upper()
                acceptor = genome.fetch(model["seqid"], intron_end - 2, intron_end).upper()
                splice = (donor, acceptor) if model["strand"] == "+" else (
                    str(Seq(acceptor).reverse_complement()), str(Seq(donor).reverse_complement()))
                if splice not in {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}:
                    problems.append("noncanonical_splice")
    if set(seq) - set("ACGT"):
        problems.append("assembly_gap_or_ambiguity")
    span = genome.fetch(model["seqid"], min(e[0] for e in exons), max(e[1] for e in exons)).upper()
    gap_bases = len(span) - span.count("A") - span.count("C") - span.count("G") - span.count("T")
    if gap_bases:
        problems.append("assembly_gap_within_model_span")
    if len(seq) % 3:
        problems.append("incomplete_frame")
    else:
        codons = CodonTable.unambiguous_dna_by_id[code]
        if seq[:3] not in codons.start_codons:
            problems.append("missing_start")
        # If the predictor omits the terminal stop, extend only if the genome
        # supplies that exact next codon; never manufacture or mask residues.
        if seq[-3:] not in codons.stop_codons:
            last = exons[-1]
            start, end = (last[1], last[1] + 3) if model["strand"] == "+" else (last[0] - 3, last[0])
            stop = genome.fetch(model["seqid"], max(0, start), min(end, genome.get_reference_length(model["seqid"]))).upper()
            if model["strand"] == "-":
                stop = str(Seq(stop).reverse_complement())
            if len(stop) == 3 and stop in codons.stop_codons and start >= 0:
                seq += stop
                if model["strand"] == "+":
                    last[1] += 3
                else:
                    last[0] -= 3
            else:
                problems.append("missing_stop")
        translation = str(Seq(seq).translate(table=code))
        if "*" in translation[:-1] or ("missing_stop" in problems and "*" in translation):
            problems.append("internal_stop")
    return {**model, "cds": exons, "sequence": seq, "assembly_ambiguous_bases": gap_bases,
            "problems": sorted(set(problems))}


def overlaps(a, b):
    return a["seqid"] == b["seqid"] and a["start"] < b["end"] and b["start"] < a["end"]


def check_interval(model):
    region = model["evidence"]
    if region.get("genome_only"):
        model["problems"].append("unanchored_genome_search")
        return model
    if (model["seqid"] != region["seqid"] or not model["cds"]
            or min(s for s, _, _ in model["cds"]) < region["expected_start"]
            or max(e for _, e, _ in model["cds"]) > region["expected_end"]):
        model["problems"].append("outside_expected_synteny_interval")
    return model


def search_intervals(tmp, windows, proteins, genome, code, max_intron, cpus, interval_workers=None,
                     retain_inputs=False):
    """Bound pending work and retain input order, with total threads <= cpus.

    Only the submitting thread fetches from pysam's seekable FASTA handle.
    Alignments remain independent so other WGD intervals cannot compete.
    """
    workers = min(cpus, len(windows)) if interval_workers is None else interval_workers
    if cpus < 1 or (windows and not 1 <= workers <= cpus):
        raise ValueError("Interval workers must be positive and no greater than cpus")
    interval_dir = tmp / "intervals"
    interval_dir.mkdir()
    if not windows:
        return []
    threads = cpus // workers
    def predict(directory):
        run(["miniprot", "-u", "-T", code, "-t", threads, "-G", max_intron,
             "-N", MINIPROT_MAX_SECONDARY, f"--outs={MINIPROT_OUTPUT_SCORE_RATIO}",
             "--gff", directory / "region.fa", directory / "queries.fa"],
            directory, "miniprot", directory / "models.gff")
        result = read_miniprot(directory / "models.gff")
        if not retain_inputs:
            (directory / "region.fa").unlink()
            (directory / "queries.fa").unlink()
        return result
    pending, predictions = deque(), []
    with gzip.open(tmp / LOCAL_RECORDS, "wt", compresslevel=1) as metadata, ThreadPoolExecutor(max_workers=workers) as executor:
        for i, ((seqid, start, end), queries) in enumerate(windows.items(), 1):
            directory = interval_dir / str(i)
            directory.mkdir()
            sequence = genome.fetch(seqid, start, end)
            (directory / "region.fa").write_text(f">interval\n{sequence}\n")
            write_record(metadata, local_record(i, (seqid, start, end), queries, proteins, sequence))
            with (directory / "queries.fa").open("w") as out:
                for region in queries:
                    out.write(f">{region['id']}\n{proteins[region['donor']][region['query']]}\n")
            pending.append(executor.submit(predict, directory))
            if len(pending) >= workers * 2:
                predictions.extend(pending.popleft().result())
        for future in pending:
            predictions.extend(future.result())
    return predictions


def write_local_search_inputs(directory, windows, proteins, genome):
    """Materialize diagnostics only for local queries requiring a new search.

    Candidate evidence and genome-search coverage live in their verified JSON/
    TSV receipts. These FASTAs are optional local-search diagnostics, so a fully
    covered cache needs neither file nor a repeated fetch of every old interval.
    """
    if not windows:
        return
    with (directory / "regions.fa").open("w") as out, (directory / "queries.fa").open("w") as queries:
        for (seqid, start, end), regions in windows.items():
            sequence = genome.fetch(seqid, start, end)
            for region in regions:
                queries.write(f">{region['id']}\n{proteins[region['donor']][region['query']]}\n")
                out.write(f">{region['id']}\n{sequence}\n")


def write_unique_queries(regions, proteins, path):
    """Exact sequence equality only; map every original candidate to its first representative."""
    representatives, by_sequence = {}, {}
    with path.open("w") as out:
        for region in regions:
            sequence = proteins[region["donor"]][region["query"]]
            if region["id"] in representatives:
                raise ValueError("Duplicate candidate ID in genome search")
            representative = by_sequence.get(sequence)
            if representative is None:
                representative = region["id"]
                by_sequence[sequence] = representative
                out.write(f">{representative}\n{sequence}\n")
            representatives[region["id"]] = representative
    return representatives


def expand_miniprot_queries(source, destination, representatives):
    """Restore query order, names, PAF and GFF IDs before the existing model reader/QC."""
    lines = iter(expand_miniprot_lines(source, representatives))
    first = next(lines, None)
    with destination.open("w") as out:
        if first is not None:
            out.write(first)
        out.writelines(lines)


def expand_miniprot_lines(source, representatives):
    """Restore the legacy GFF stream without an expanded persistent GFF copy."""
    blocks, headers = defaultdict(list), []
    query, lines = None, []
    known = set(representatives.values())
    with source.open() as handle:
        for line in handle:
            if line.startswith("##PAF\t"):
                if query is not None:
                    blocks[query].append("".join(lines))
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 13 or fields[1] not in known:
                    raise ValueError("Genome search returned an invalid/unknown representative")
                query, lines = fields[1], [line]
            elif query is None:
                if not line.startswith("#") and line.strip():
                    raise ValueError("Genome search GFF lacks embedded PAF")
                headers.append(line)
            else:
                lines.append(line)
    if query is not None:
        blocks[query].append("".join(lines))
    if set(blocks) != known:
        raise ValueError("Genome search did not report every representative with -u")
    number = 0
    yield from headers
    for original, representative in representatives.items():
        for block in blocks[representative]:
            identifiers = {}
            for line in block.splitlines(keepends=True):
                if line.startswith("##PAF\t"):
                    fields = line.rstrip("\n").split("\t")
                    fields[1] = original
                    yield "\t".join(fields) + "\n"
                elif not line.startswith("#") and line.strip():
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        raise ValueError("Invalid representative miniprot GFF row")
                    values = fields[8].split(";")
                    if fields[2] == "mRNA":
                        old = attributes(fields[8])["ID"]
                        match = re.fullmatch(r"MP(\d+)", old)
                        if not match:
                            raise ValueError("Unexpected miniprot model ID")
                        number += 1
                        identifiers[old] = "MP" + str(number).zfill(len(match[1]))
                    rewritten = []
                    for value in values:
                        if value.startswith(("ID=", "Parent=")):
                            key, token = value.split("=", 1)
                            if token not in identifiers:
                                raise ValueError("Unbound representative miniprot feature")
                            value = key + "=" + identifiers[token]
                        elif value.startswith("Target="):
                            target = value[len("Target="):].split(" ", 1)
                            if target[0] != representative or len(target) != 2:
                                raise ValueError("Inconsistent representative miniprot Target")
                            value = "Target=" + original + " " + target[1]
                        rewritten.append(value)
                    fields[8] = ";".join(rewritten)
                    yield "\t".join(fields) + "\n"
                else:
                    yield line


def output_storage(plan):
    value = plan["request"].get("output_storage", {"format": "legacy", "retain_search_inputs": True})
    if (not isinstance(value, dict) or set(value) != {"format", "retain_search_inputs"}
            or not isinstance(value["format"], str) or value["format"] not in {"compact", "legacy"}
            or type(value["retain_search_inputs"]) is not bool):
        raise ValueError("Unknown or malformed rescue output-storage configuration")
    return value


def rescue(root, plan, name, cpus, interval_workers=None):
    if name not in plan["species"]:
        raise ValueError("Unknown rescue species")
    storage = output_storage(plan)
    donors = [name, *plan["donors"][name]]
    for donor in donors:
        prepared(root, plan, donor)
    key = rescue_key(root, plan, name)
    regions = candidates(root, plan, name)
    verify_sources(plan, [name], ["genome", "fasta", "gff"])
    source = plan["request"]["sources"][name]
    params = parameters_for(plan["request"], name)
    additional, nomination = nominate_genome_only_candidates(root, plan, name, regions, params)
    regions.extend(additional)
    def build(tmp):
        from gene_model_catalog import indexed_genome
        expanded_genome_sha256 = None
        atomic_json(tmp / SEARCH_CONTRACT_FILE, prediction_search_contract())
        atomic_json(tmp / "genome_search_nomination.json", nomination)
        if not regions:
            # No copy/index/predictor is needed, but invalid source annotations
            # must still fail. Stream contig lengths without materializing DNA.
            lengths = {}
            identifier = None
            with open_text(Path(source["genome"])) as handle:
                for line in handle:
                    if line.startswith(">"):
                        token = line[1:].split()
                        if not token or token[0] in lengths:
                            raise ValueError("FASTA index warning: missing or duplicate contig")
                        identifier = token[0]
                        lengths[identifier] = 0
                    elif line.strip():
                        if identifier is None:
                            raise ValueError("FASTA index warning: sequence before header")
                        lengths[identifier] += len(line.strip())
            with open_text(Path(source["gff"])) as handle:
                for line in handle:
                    if line.strip() == "##FASTA":
                        break
                    if line.startswith("#") or not line.strip():
                        continue
                    f = line.rstrip().split("\t")
                    if (len(f) == 9 and f[2] in {"gene", "mRNA", "transcript", "CDS"}
                            and (f[0] not in lengths or not 0 <= int(f[3]) - 1 < int(f[4]) <= lengths[f[0]])):
                        raise ValueError("Original annotation outside or absent from genome")
            atomic_json(tmp / "candidates.json", [])
            if storage["format"] == "compact":
                write_model_store(tmp, [], revisions=[])
            else:
                atomic_json(tmp / "models.json", [])
                atomic_json(tmp / "revision_candidates.json", [])
                atomic_json(tmp / "partial_models.json", [])
            write_tsv(tmp / "audit.tsv", ("candidate", "donor_gene", "status", "reasons", "coverage", "identity", "model_id"), [])
            write_tsv(tmp / "quality_flags.tsv", ("candidate", "model_id", "status", "start_codon",
                      "donor_n_terminus_aligned", "donor_c_terminus_aligned", "donor_species", "flags",
                      "donor_aligned_query_fraction", "donor_internal_unaligned_query_fraction"), [])
            return
        with open_text(Path(source["gff"])) as handle:
            existing = original_ownership(handle, attributes)
        placement_ownership = OwnershipIndex(existing)
        proteins = {donor: {i: s for i, _, s in fasta_records(root / "prepared" / donor / "genes.pep")} for donor in donors}
        frozen_cache = plan["request"].get("prediction_cache")
        with indexed_genome(source["genome"]) as genome:
            lengths = dict(zip(genome.references, genome.lengths, strict=True))
            positions = json.loads((root / "prepared" / name / "positions.json").read_text())
            for feature in [*existing, *positions]:
                if (feature["seqid"] not in lengths or feature["start"] < 0
                        or feature["end"] <= feature["start"] or feature["end"] > lengths[feature["seqid"]]):
                    raise ValueError("Original annotation outside or absent from genome: " + str(feature))
            for region in regions:
                if not region.get("genome_only"):
                    region["end"] = min(region["end"], genome.get_reference_length(region["seqid"]))
            cached = (verify_prediction_cache(frozen_cache["root"], root, plan, name, regions, params, frozen_cache)
                      if frozen_cache else None)
            cached_ids = cached.candidate_ids() if cached else set()
            cached_genome_ids = cached.genome_candidate_ids() if cached else set()
            atomic_json(tmp / "candidates.json", regions)
            windows = defaultdict(list)
            for region in regions:
                if region.get("genome_only") or region["id"] in cached_ids:
                    continue
                windows[(region["seqid"], region["start"], region["end"])].append(region)
            if storage["retain_search_inputs"]:
                write_local_search_inputs(tmp, windows, proteins, genome)
            predictions = search_intervals(tmp, windows, proteins, genome, source["genetic_code"],
                                           params["max_intron"], cpus, interval_workers,
                                           retain_inputs=storage["retain_search_inputs"])
            by_id = {r["id"]: r for r in regions}
            validated = []
            raw_validation = RawValidationMemo(validate_model, genome, source["genetic_code"], params)
            def checked_prediction(model):
                # Preserve predictor coordinates and metrics before any genomic
                # completion, so a later plan rechecks its own terminal bounds.
                raw_prediction = {key: model[key] for key in RAW_FIELDS if key in model}
                raw_prediction["cds"] = [list(block) for block in model["cds"]]
                raw_prediction["evidence"] = dict(model["evidence"])
                checked = raw_validation.validate(model)
                evidence = model["evidence"]
                checked = complete_terminals(checked, genome, source["genetic_code"], params,
                                             proteins[evidence["donor"]][evidence["query"]], validate_model)
                checked["raw_prediction"] = raw_prediction
                return check_interval(checked)
            if cached:
                for model in cached.iter_models():
                    validated.append(checked_prediction(model))
            for model in predictions:
                region = by_id[model["query"]]
                model["seqid"] = region["seqid"]
                model["cds"] = [[s + region["start"], e + region["start"], p] for s, e, p in model["cds"]]
                model["evidence"] = region
                model["search"] = "synteny_interval"
                validated.append(checked_prediction(model))
            def unresolved_regions():
                preview = [{**m, "problems": list(m["problems"])} for m in validated]
                reassess_unanchored_models(preview, params["unanchored_min_species"], ownership=placement_ownership)
                preview = consolidate(preview, existing, name, plan.get("nearest_references", {}).get(name, ()))
                resolved = {m["query"] for m in preview if m["status"] in {"accepted", "duplicate_support", "accepted_alternative_path"}}
                return [r for r in regions if r["id"] not in resolved]
            unresolved = unresolved_regions()
            # Old local-search coverage and old whole-genome coverage are
            # separate, because acceptance can change under the new policy.
            unresolved = [r for r in unresolved if r["id"] not in cached_genome_ids]
            searched_genome_ids = set(cached_genome_ids)
            if params["genome_fallback"] and unresolved:
                if storage["retain_search_inputs"]:
                    with (tmp / "unresolved.fa").open("w") as out:
                        for r in unresolved:
                            out.write(f">{r['id']}\n{proteins[r['donor']][r['query']]}\n")
                representatives = write_unique_queries(unresolved, proteins, tmp / "unresolved.unique.fa")
                write_genome_records(tmp, unresolved, representatives, proteins)
                write_tsv(tmp / "genome_prediction_query_mapping.tsv", ("candidate", "representative"), representatives.items())
                searched_genome_ids.update(representatives)
                (tmp / "genome.fa").symlink_to(os.fsdecode(genome.filename))
                run(["miniprot", "-T", source["genetic_code"], "-t", cpus, "-d", tmp / "genome.mpi", tmp / "genome.fa"], tmp, "miniprot_index")
                run(["miniprot", "-u", "-t", cpus, "-G", params["max_intron"],
                     "-N", MINIPROT_MAX_SECONDARY, f"--outs={MINIPROT_OUTPUT_SCORE_RATIO}", "--gff",
                     tmp / "genome.mpi", tmp / "unresolved.unique.fa"], tmp, "miniprot_genome", tmp / "genome.unique.gff")
                if storage["retain_search_inputs"]:
                    expand_miniprot_queries(tmp / "genome.unique.gff", tmp / "genome.gff", representatives)
                expanded_digest = hashlib.sha256()
                def hashed_expansion():
                    for line in expand_miniprot_lines(tmp / "genome.unique.gff", representatives):
                        expanded_digest.update(line.encode())
                        yield line
                for model in iter_miniprot(hashed_expansion()):
                    model["evidence"] = by_id[model["query"]]
                    model["search"] = "genome_fallback"
                    validated.append(checked_prediction(model))
                expanded_genome_sha256 = expanded_digest.hexdigest()
                unresolved = unresolved_regions()
            if plan["request"]["gemoma_jar"] and unresolved:
                # Selected transcript IDs and target coordinates constrain the
                # optional refinement. Predictions still pass the same QC.
                gemoma_regions = [r for r in unresolved if not r.get("genome_only")]
                if gemoma_regions:
                    if not (tmp / "genome.fa").exists():
                        (tmp / "genome.fa").symlink_to(os.fsdecode(genome.filename))
                    refine_gemoma(tmp, root, plan, source, gemoma_regions, genome, validated, cpus)
            if cached:
                cached.check()
            if searched_genome_ids:
                covered = [r for r in regions if r["id"] in searched_genome_ids]
                coverage_mapping = write_unique_queries(covered, proteins, tmp / "genome.covered.unique.fa")
                write_tsv(tmp / "genome_query_mapping.tsv", ("candidate", "representative"), coverage_mapping.items())
            atomic_json(tmp / "raw_validation_memo.json", raw_validation.diagnostics())
        # Recheck sources before publication, including the target genome.
        verify_sources(plan, [name], ["genome", "fasta", "gff"])
        placement = reassess_unanchored_models(validated, params["unanchored_min_species"], ownership=placement_ownership)
        models = consolidate(validated, existing, name, plan.get("nearest_references", {}).get(name, ()))
        for model in models:
            model["quality_evidence"] = model_quality(model, source["genetic_code"])
        atomic_json(tmp / "placement_audit.json", placement)
        revisions = {}
        for model in models:
            if len(model.get("revision_owner_ids", [])) != 1 or set(model["problems"]) - {"overlap_existing_annotation"}:
                continue
            shape = (*coding_shape(model), model["revision_owner_ids"][0])
            evidence = {**model["evidence"], "alignment": {"coverage": model["coverage"],
                        "identity": model["identity"], "problems": list(model["problems"])}}
            if shape not in revisions:
                revisions[shape] = {**model, "support": []}
            revisions[shape]["support"].append(evidence)
        if storage["format"] == "compact":
            write_model_store(tmp, models, revisions=revisions.values())
        else:
            atomic_json(tmp / "models.json", models)
            atomic_json(tmp / "revision_candidates.json", list(revisions.values()))
            atomic_json(tmp / "partial_models.json", [m for m in models if m.get("partial_evidence", {}).get("partial")])
        atomic_json(tmp / "search_inputs.json", {
            "schema": 1, "plan_sha256": plan_digest(root, plan),
            "source_files": {source[k]: plan["request"]["files"][source[k]] for k in ("genome", "fasta", "gff")},
            "prepared_receipts": {donor: digest(root / "prepared" / donor / "receipt.json") for donor in donors},
            "retained_search_inputs": storage["retain_search_inputs"],
            "local_records": LOCAL_RECORDS if (tmp / LOCAL_RECORDS).exists() else None,
            "expanded_genome_gff_sha256": expanded_genome_sha256,
            "prediction_search_contract": prediction_search_contract(),
            "prediction_cache": plan["request"].get("prediction_cache"),
        })
        atomic_json(tmp / "prediction_reuse.json", {"cached_local_queries": len(cached_ids),
                    "cached_genome_queries": len(cached_genome_ids), "searched_local_windows": len(windows),
                    "additional_genome_queries": len(additional), "all_queries": len(regions),
                    "search_contract": prediction_search_contract(),
                    "cached_search_contract": cached.search_contract if cached else None,
                    "cached_local_search_compatible": cached.local_search_compatible if cached else None,
                    "cached_genome_search_compatible": cached.genome_search_compatible if cached else None})
        write_tsv(tmp / "quality_flags.tsv", ("candidate", "model_id", "status", "start_codon",
                  "donor_n_terminus_aligned", "donor_c_terminus_aligned", "donor_species", "flags",
                  "donor_aligned_query_fraction", "donor_internal_unaligned_query_fraction"),
                  [(m["query"], m.get("model_id", ""), m["status"], m["quality_evidence"]["start_codon"],
                    m["quality_evidence"]["terminal_alignment"]["n_aligned"],
                    m["quality_evidence"]["terminal_alignment"]["c_aligned"],
                    ",".join(m["quality_evidence"]["donor_species"]), ",".join(m["quality_evidence"]["flags"]),
                    m["quality_evidence"]["query_alignment"]["aligned_query_fraction"],
                    m["quality_evidence"]["query_alignment"]["internal_unaligned_query_fraction"])
                   for m in models])
        detected = {m["query"] for m in validated}
        audit = [{"candidate": m["query"], "donor_gene": m["evidence"]["query"], "status": m["status"],
                  "reasons": ",".join(m["problems"]), "coverage": m["coverage"], "identity": m["identity"],
                  "model_id": m.get("model_id", "")} for m in models]
        audit += [{"candidate": r["id"], "donor_gene": r["query"], "status": "unresolved", "reasons": "no_alignment",
                   "coverage": "", "identity": "", "model_id": ""} for r in regions if r["id"] not in detected]
        write_tsv(tmp / "audit.tsv", ("candidate", "donor_gene", "status", "reasons", "coverage", "identity", "model_id"),
                  [list(row.values()) for row in audit])
        # Huge indexes are execution scratch, not reusable unverified outputs.
        (tmp / "genome.mpi").unlink(missing_ok=True)
        (tmp / "genome.fa").unlink(missing_ok=True)
        if not storage["retain_search_inputs"]:
            (tmp / "unresolved.unique.fa").unlink(missing_ok=True)
            (tmp / "genome.covered.unique.fa").unlink(missing_ok=True)
    return stage(root, Path("rescued") / name, key, build,
                 lambda: require_same_key(key, rescue_key(root, plan, name)), hash_workers=cpus)


def consolidate(validated, existing, species, nearest=()):
    ownership = OwnershipIndex(existing)
    unique = {}
    for model in validated:
        shape = coding_shape(model)
        model["start"] = min((e[0] for e in model["cds"]), default=0)
        model["end"] = max((e[1] for e in model["cds"]), default=0)
        owners = ownership.overlapping(model)
        if owners:
            model["problems"] = sorted(set(model["problems"]) | {"overlap_existing_annotation"})
            identifiers = sorted({o["gene_id"] for o in owners if o.get("gene_id")})
            model["revision_owner_ids"] = identifiers
            model["revision_disposition"] = "existing_model_refinement" if len(identifiers) == 1 else "ambiguous_split_or_merge"
        model["status"] = "unresolved" if model["problems"] else "accepted"
        if model["status"] == "accepted":
            model["model_id"] = species + "_ggrescue_" + hashlib.sha256(repr(shape).encode()).hexdigest()[:16]
            evidence = {**model.get("evidence", {}), "alignment": {
                "coverage": model.get("coverage", 0), "identity": model.get("identity", 0),
                "problems": list(model["problems"])}}
            if shape in unique:
                unique[shape]["support"].append(evidence)
                model["status"] = "duplicate_support"
            else:
                model["support"] = [evidence]
                unique[shape] = model
    accepted = sorted(unique.values(), key=lambda m: (m["seqid"], m["start"], m["end"]))
    resolve_coding_paths(accepted, nearest)
    for model in validated:
        shape = coding_shape(model)
        if model["status"] == "duplicate_support" and unique[shape]["status"] == "unresolved":
            model["status"] = "unresolved"
            model["problems"] = sorted(set(model["problems"]) | set(unique[shape]["problems"]))
    return validated


def refine_gemoma(tmp, root, plan, source, regions, genome, validated, cpus):
    jar = plan["request"]["gemoma_jar"]
    if digest(jar) != plan["request"]["files"][jar]:
        raise ValueError("GeMoMa jar changed")
    java = plan["request"]["gemoma_java"]
    if digest(java) != plan["request"]["files"][java]:
        raise ValueError("GeMoMa Java executable changed")
    params = parameters_for(plan["request"], source.get("species", ""))
    for donor in sorted({r["donor"] for r in regions}):
        ref = plan["request"]["sources"][donor]
        verify_sources(plan, [donor], ["genome", "gff"])
        mapping = {r["jcvi_id"]: r for r in table(root / "prepared" / donor / "genes.id_map.tsv") if r["status"] == "selected"}
        admission = {r["original_id"]: r for r in json.loads(
            (root / "prepared" / donor / "genes.anchor_admission.json").read_text())["records"]}
        transcripts = {}
        with open_text(Path(ref["gff"])) as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) == 9 and fields[2] in {"mRNA", "transcript"}:
                    attr = attributes(fields[8])
                    if "ID" in attr:
                        transcripts[attr["ID"]] = attr.get("Parent", "").split(",")
        subset = [r for r in regions if r["donor"] == donor]
        directory = tmp / ("gemoma_" + donor)
        directory.mkdir()
        proteins = None
        # Separate queries avoid merging duplicate transcript IDs across WGD
        # intervals in GeMoMa's selected-file lookup.
        for r in subset:
            row = mapping[r["query"]]
            evidence = admission.get(row["original_id"], {}).get("selected_evidence", [])
            if any(e.get("phase_convention") == "complementary" for e in evidence):
                (directory / (r["id"] + ".skipped.json")).write_text(json.dumps(
                    {"reason": "unsupported_reference_phase_convention", "evidence": evidence}) + "\n")
                continue
            aliases = {row["original_id"], row["original_id"].removeprefix(donor + "_")}
            matches = sorted(set(transcripts) & aliases)
            if not matches:
                matches = sorted(t for t, parents in transcripts.items() if row["locus_id"] in parents)
            if len(matches) != 1 or r["expected_strand"] not in {"+", "-"}:
                (directory / (r["id"] + ".skipped.json")).write_text(json.dumps({"reason": "ambiguous_reference_transcript_or_strand", "transcripts": matches}) + "\n")
                continue
            selected = directory / (r["id"] + ".tsv")
            selected.write_text(f"{matches[0]}\t{r['seqid']}\t{r['expected_strand']}\t{r['start'] + 1}\t{r['end']}\n")
            out = directory / r["id"]
            run([java, "-jar", jar, "CLI", "GeMoMaPipeline", "s=own", f"i={donor}",
                 f"a={ref['gff']}", f"g={ref['genome']}", f"t={tmp / 'genome.fa'}", f"selected={selected}",
                 "tblastn=true", "AnnotationFinalizer.r=NO", "GeMoMa.Score=ReAlign",
                 f"threads={cpus}", f"outdir={out}"], directory, r["id"])
            result = out / "final_annotation.gff"
            if not result.exists():
                raise RuntimeError("GeMoMa did not write final_annotation.gff")
            groups = {}
            for line in result.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                if len(f) != 9:
                    raise ValueError("Invalid GeMoMa GFF row")
                attr = attributes(f[8])
                if f[2] in {"mRNA", "transcript"}:
                    groups[attr["ID"]] = {"query": r["id"], "seqid": f[0], "strand": f[6], "cds": [],
                                         "frameshift": False, "coverage": 0.0, "identity": 0.0,
                                         "evidence": r, "search": "gemoma", "id": attr["ID"]}
                elif f[2] == "CDS":
                    for parent in attr["Parent"].split(","):
                        if parent not in groups or (f[0], f[6]) != (groups[parent]["seqid"], groups[parent]["strand"]):
                            raise ValueError("Unbound or inconsistent GeMoMa CDS")
                        groups[parent]["cds"].append([int(f[3]) - 1, int(f[4]), int(f[7])])
            # Re-align each refined CDS to its donor protein for comparable QC.
            for model in groups.values():
                checked = validate_model(model, genome, source["genetic_code"],
                                         {**params, "minimum_coverage": 0, "minimum_identity": 0})
                protein = str(Seq(checked["sequence"][:len(checked["sequence"]) // 3 * 3]).translate(table=source["genetic_code"])).rstrip("*")
                if proteins is None:
                    proteins = {i: s for i, _, s in fasta_records(root / "prepared" / donor / "genes.pep")}
                reference = proteins[r["query"]]
                over_budget = len(reference) * len(protein) > params.get("terminal_max_alignment_cells", 25_000_000)
                alignment = None
                if protein and not over_budget:
                    from Bio.Align import PairwiseAligner
                    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                                              open_gap_score=-5, extend_gap_score=-1)
                    alignment = aligner.align(reference, protein)[0]
                aligned = sum(b - a for a, b in alignment.aligned[0]) if alignment else 0
                matches = sum(reference[a + k] == protein[c + k] for (a, b), (c, d) in zip(*alignment.aligned, strict=True)
                              for k in range(b - a)) if alignment else 0
                checked["coverage"] = aligned / len(reference)
                checked["identity"] = matches / aligned if aligned else 0
                if aligned:
                    checked.update(query_start=int(alignment.aligned[0][0][0]),
                                   query_end=int(alignment.aligned[0][-1][1]), query_length=len(reference))
                    checked["query_span_coverage"] = (checked["query_end"] - checked["query_start"]) / len(reference)
                if checked["coverage"] < params["minimum_coverage"]:
                    checked["problems"].append("low_coverage")
                if checked["identity"] < params["minimum_identity"]:
                    checked["problems"].append("low_identity")
                if over_budget:
                    checked["problems"].append("donor_alignment_budget_exceeded")
                raw_prediction = {key: checked[key] for key in RAW_FIELDS if key in checked}
                raw_prediction["cds"] = [list(block) for block in model["cds"]]
                raw_prediction["evidence"] = dict(model["evidence"])
                checked = complete_terminals(checked, genome, source["genetic_code"], params, reference, validate_model)
                checked["raw_prediction"] = raw_prediction
                validated.append(check_interval(checked))


def rescue_selection_attributes(model, path=None):
    """Export locus corroboration without assigning it to every coding path."""
    selection = model.get("path_selection")
    if not selection:
        return ""
    locus = model["locus_support"]
    values = [("representative_status", selection["representative_status"]),
              ("representative_selection_policy", selection["selection_policy"]),
              ("locus_independent_donor_species_count", len(locus["independent_donor_species"])),
              ("rescue_orthology", locus["orthology"]), ("rescue_expected_copy", locus["expected_copy"])]
    if path is not None:
        targets = {row["target"] for item in [model, *model.get("alternative_coding_paths", [])]
                   for row in [item.get("evidence", {}), *item.get("support", [])] if row.get("target")}
        donors = {donor for donor, _, _ in alignment_support(path)} - targets
        values.append(("path_independent_donor_species_count", len(donors)))
    return ";" + ";".join(f"{key}={quote(str(value), safe='')}" for key, value in values)


def finalize(root, plan, names=None, destination=Path("augmented"), hash_workers=1):
    names = names or plan["species"]
    key = {"plan": plan_digest(root, plan), "rescue_receipts": {
        n: digest(root / "rescued" / n / "receipt.json") for n in names}}
    def guard():
        for name in names:
            if not verified(root / "rescued" / name, rescue_key(root, plan, name), hash_workers=hash_workers):
                raise ValueError("Rescue incomplete or corrupted: " + name)
        verify_sources(plan, names, ["fasta", "gff", "genome", "busco"])
        require_same_key(key, {"plan": plan_digest(root, plan), "rescue_receipts": {
            n: digest(root / "rescued" / n / "receipt.json") for n in names}})
    def build(tmp):
        for subdir in ("species_cds", "species_gff", "anchor_admission"):
            (tmp / subdir).mkdir()
        rows = []
        admission_summaries = {}
        alternative_counts = {}
        for name in names:
            source = plan["request"]["sources"][name]
            models = list(iter_accepted_models(root / "rescued" / name))
            alternative_counts[name] = sum(len(m.get("alternative_coding_paths", [])) for m in models)
            cds = tmp / "species_cds" / (name + ".rescue.cds.fa")
            gff = tmp / "species_gff" / (name + ".rescue.gff3")
            existing_ids = set()
            with cds.open("w") as out:
                for identifier, header, sequence in fasta_records(Path(source["fasta"])):
                    if identifier in existing_ids:
                        raise ValueError("Duplicate original CDS ID")
                    existing_ids.add(identifier)
                    out.write(f">{header}\n{sequence}\n")
                for model in models:
                    if model["model_id"] in existing_ids:
                        raise ValueError("Rescued ID collides with original")
                    out.write(f">{model['model_id']}\n{model['sequence']}\n")
            with open_text(Path(source["gff"])) as handle:
                original = handle.read()
            lines = original.splitlines(keepends=True)
            boundary = next((i for i, line in enumerate(lines) if line.strip() == "##FASTA"), len(lines))
            annotation, embedded = "".join(lines[:boundary]), "".join(lines[boundary:])
            gff_ids = {attributes(line.rstrip().split("\t")[8])["ID"] for line in lines[:boundary]
                       if not line.startswith("#") and len(line.rstrip().split("\t")) == 9
                       and "ID" in attributes(line.rstrip().split("\t")[8])}
            mapping = json.loads((root / "prepared" / name / "mapping.json").read_text())
            with gff.open("w") as out:
                out.write(annotation)
                if annotation and not annotation.endswith("\n"):
                    out.write("\n")
                for model in models:
                    identifier = model["model_id"]
                    common = [model["seqid"], "genegalleon_rescue", "", str(model["start"] + 1),
                              str(model["end"]), ".", model["strand"], ".", ""]
                    feature, attribute = mapping["feature"], mapping["attribute"]
                    gene_level = feature == "gene" or (feature in {"mRNA", "transcript"} and attribute == "Parent")
                    gene = identifier if gene_level else identifier + ".gene"
                    transcript = identifier + ".t1" if gene_level or (feature == "CDS" and attribute == "ID") else identifier
                    cds_ids = [identifier if feature == "CDS" and attribute == "ID" else f"{identifier}.cds{i}"
                               for i in range(1, len(model["cds"]) + 1)]
                    new_ids = {gene, transcript, *cds_ids}
                    if new_ids & gff_ids:
                        raise ValueError("Rescued GFF ID collides with original or another new model")
                    gff_ids.update(new_ids)
                    def fields(kind, attr, feature=feature, attribute=attribute, identifier=identifier):
                        if kind == feature and attribute not in {"ID", "Parent"}:
                            attr += f";{attribute}={identifier}"
                        return attr
                    alternative_paths = model.get("alternative_coding_paths", [])
                    gene_start = min([model["start"], *[min(b[0] for b in a["cds"]) for a in alternative_paths]])
                    gene_end = max([model["end"], *[max(b[1] for b in a["cds"]) for a in alternative_paths]])
                    locus_attributes = rescue_selection_attributes(model)
                    out.write("\t".join([*common[:2], "gene", str(gene_start + 1), str(gene_end), *common[5:8], fields("gene", "ID=" + gene + locus_attributes)]) + "\n")
                    transcript_feature = "transcript" if feature == "transcript" else "mRNA"
                    out.write("\t".join([*common[:2], transcript_feature, *common[3:8], fields(transcript_feature, f"ID={transcript};Parent={gene}" + rescue_selection_attributes(model, model))]) + "\n")
                    for cds_id, (start, end, phase) in zip(cds_ids, model["cds"], strict=True):
                        out.write("\t".join([model["seqid"], "genegalleon_rescue", "CDS", str(start + 1), str(end), ".",
                                            model["strand"], str(phase), fields("CDS", f"ID={cds_id};Parent={transcript}")]) + "\n")
                    for alternative in alternative_paths:
                        alt_id = identifier + ".alt_" + alternative["model_id"].rsplit("_", 1)[-1]
                        start, end = min(b[0] for b in alternative["cds"]), max(b[1] for b in alternative["cds"])
                        ids = {alt_id, *[f"{alt_id}.cds{i}" for i in range(1, len(alternative["cds"]) + 1)]}
                        if ids & gff_ids:
                            raise ValueError("Alternative coding path GFF ID collides")
                        gff_ids.update(ids)
                        out.write("\t".join([model["seqid"], "genegalleon_rescue", transcript_feature, str(start + 1), str(end),
                                            ".", model["strand"], ".", f"ID={alt_id};Parent={gene};support=homology_coding_path" + rescue_selection_attributes(model, alternative)]) + "\n")
                        for i, (s, e, phase) in enumerate(alternative["cds"], 1):
                            out.write("\t".join([model["seqid"], "genegalleon_rescue", "CDS", str(s + 1), str(e), ".",
                                                model["strand"], str(phase), f"ID={alt_id}.cds{i};Parent={alt_id}"]) + "\n")
                out.write(embedded)
            # Exercise the same full-ID annotation mapping used by downstream synteny.
            check = tmp / ("check_" + name)
            check.mkdir()
            mapping = json.loads((root / "prepared" / name / "mapping.json").read_text())
            _, checked_mapping = prepare_rescue_genome(
                {**source, "fasta": str(cds), "gff": str(gff), "feature": mapping["feature"],
                 "attribute": mapping["attribute"]}, check, "check", 1.0,
                required_ids=[m["model_id"] for m in models])
            admission_summaries[name] = checked_mapping["anchor_admission"]
            for suffix in ("json", "tsv"):
                shutil.copyfile(check / ("check.anchor_admission." + suffix),
                                tmp / "anchor_admission" / (name + "." + suffix))
            shutil.rmtree(check)
            rows.append((name, len(existing_ids), len(models), len(existing_ids) + len(models),
                         source["quality"]["complete_pct"], str(root / destination / "species_cds" / cds.name),
                         str(root / destination / "species_gff" / gff.name), source["genome"]))
        write_tsv(tmp / "inputs.tsv", ("species", "original_cds", "rescued_models", "augmented_cds", "busco_before_pct",
                                        "cds", "gff", "genome"), rows)
        atomic_json(tmp / "summary.json", {"counting_unit": "new gene model", "species": len(rows),
                                           "rescued_models": sum(r[2] for r in rows), "changed_species": [r[0] for r in rows if r[2]],
                                           "common_references": plan["common_references"],
                                           "alternative_coding_paths": alternative_counts,
                                           "anchor_admission": admission_summaries,
                                           "gene_loss_calls": False, "plan_sha256": plan_digest(root, plan)})
    return stage(root, destination, key, build, guard, hash_workers=hash_workers)


def qc_work_items(root, plan, indices, hash_workers=1):
    """Verify the exact effective inputs before dispatching or reusing BUSCO."""
    result = []
    for index in indices:
        name = plan["species"][index - 1]
        if not verified(root / "rescued" / name, rescue_key(root, plan, name), hash_workers=hash_workers):
            raise ValueError("Rescue incomplete or corrupted: " + name)
        effective_key = {"plan": plan_digest(root, plan), "rescue_receipts": {
            name: digest(root / "rescued" / name / "receipt.json")}}
        if not verified(root / "effective" / name, effective_key, hash_workers=hash_workers):
            raise ValueError("Effective inputs incomplete or corrupted: " + name)
        rows = table(root / "effective" / name / "inputs.tsv")
        if len(rows) != 1 or rows[0]["species"] != name:
            raise ValueError("Effective input table has the wrong species")
        done = verified(root / "workers" / name, {"plan": plan_digest(root, plan), "species": name}, hash_workers=hash_workers)
        result.append((index, name, rows[0]["cds"], rows[0]["rescued_models"], int(done)))
    return result


def export_destination(root, destination):
    """Keep review exports outside the frozen producer and its future stages."""
    destination = Path(destination).resolve()
    if destination.is_relative_to(Path(root).resolve()):
        raise ValueError("Review export destination must be outside the frozen rescue output")
    if destination.exists():
        raise FileExistsError("Review export destination already exists")
    return destination


def export_worker_inputs(root, plan, name, destination, *, combined=False, cpus=1):
    """Read a completed producer without rerunning it or changing its receipt."""
    try:
        from gene_model_catalog import indexed_genome
    except ImportError:
        from .gene_model_catalog import indexed_genome
    destination = export_destination(root, destination)
    directory = root / "rescued" / name
    key = rescue_key(root, plan, name)
    if not verified(directory, key, hash_workers=cpus):
        raise ValueError("Predictor input export requires a complete verified rescue")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=".working-search-inputs-", dir=destination.parent))
    tmp.rmdir()
    source = plan["request"]["sources"][name]
    donors = [name, *plan["donors"][name]]
    try:
        prepared_receipts = {}
        def verify_prepared():
            for donor in donors:
                folder = root / "prepared" / donor
                expected = {"plan": plan_digest(root, plan), "species": donor}
                if not verified(folder, expected, hash_workers=cpus):
                    raise ValueError("Predictor input export requires verified prepared donors")
                current = digest(folder / "receipt.json")
                if donor in prepared_receipts and prepared_receipts[donor] != current:
                    raise ValueError("Prepared donor changed during predictor input export")
                prepared_receipts[donor] = current
        verify_prepared()
        proteins = {donor: {i: s for i, _, s in fasta_records(root / "prepared" / donor / "genes.pep")} for donor in donors}
        verify_sources(plan, [name], ["genome", "fasta", "gff"])
        if (directory / LOCAL_RECORDS).exists() or (directory / "genome_search_inputs.jsonl.gz").exists():
            with indexed_genome(source["genome"]) as genome:
                result = export_inputs(directory, tmp, genome, proteins, combined=combined)
            if combined and (directory / "genome.unique.gff").is_file():
                expand_miniprot_queries(directory / "genome.unique.gff", tmp / "genome.gff", genome_aliases(directory))
                metadata = json.loads((directory / "search_inputs.json").read_text())
                if digest(tmp / "genome.gff") != metadata.get("expanded_genome_gff_sha256"):
                    raise ValueError("Regenerated genome GFF differs from frozen prediction stream")
        else:
            # Older publications retained the actual inputs. Copy only those
            # explicitly requested; never alter their original receipt members.
            tmp.mkdir()
            count = 0
            for interval in sorted((directory / "intervals").glob("*"), key=lambda p: int(p.name)):
                if not interval.is_dir():
                    continue
                paths = [interval / "region.fa", interval / "queries.fa"]
                if any(not path.is_file() for path in paths):
                    raise ValueError("Legacy predictor inputs lack reconstruction records")
                target = tmp / "intervals" / interval.name
                target.mkdir(parents=True)
                for path in paths:
                    shutil.copyfile(path, target / path.name)
                count += 1
            for filename in ["unresolved.unique.fa", *(["regions.fa", "queries.fa", "unresolved.fa", "genome.gff"] if combined else [])]:
                if (directory / filename).is_file():
                    shutil.copyfile(directory / filename, tmp / filename)
            result = {"local_windows": count, "combined_diagnostics": combined, "legacy_retained_inputs": True}
        verify_sources(plan, [name], ["genome", "fasta", "gff"])
        verify_prepared()
        if not verified(directory, key, hash_workers=cpus):
            raise ValueError("Rescue changed during predictor input export")
        atomic_json(tmp / "export_receipt.json", {
            "origin_receipt_sha256": digest(directory / "receipt.json"), "plan_sha256": plan_digest(root, plan),
            "prepared_receipts": prepared_receipts,
            "files": hash_outputs(tmp, workers=cpus, hash_function=digest), **result})
        tmp.rename(destination)
        return destination
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    self_parser = sub.add_parser("self-align", help="Internal LAST preparation preserving high-identity paralogs")
    self_parser.add_argument("--cpus", type=int, default=1)
    self_parser.add_argument("--cscore", type=float, default=0.7)
    plan = sub.add_parser("plan", help="Freeze sources, five common references, nearest donors and self jobs")
    for name in ("cds-dir", "gff-dir", "genome-dir", "busco-dir", "tree", "output"):
        plan.add_argument("--" + name, type=Path, required=True)
    plan.add_argument("--common-references", type=int, default=5)
    plan.add_argument("--guide-tree-receipt", type=Path, help="Optional frozen pre-rescue BUSCO k-mer guide receipt")
    plan.add_argument("--nearest-references", type=int, default=3)
    plan.add_argument("--minimum-busco", type=float, default=90)
    plan.add_argument("--genetic-code", type=int, default=1)
    plan.add_argument("--genetic-codes", type=Path)
    plan.add_argument("--species-profiles", type=Path, help="Explicit target-species prediction parameter overrides (TSV)")
    plan.add_argument("--feature", default="")
    plan.add_argument("--attribute", default="")
    plan.add_argument("--cscore", type=float, default=0.7)
    plan.add_argument("--min-anchors", type=int, default=4, help="Minimum anchors for pair comparisons; selfscan uses 4")
    plan.add_argument("--distance", type=int, default=20, help="Pair chaining distance; selfscan uses 20")
    plan.add_argument("--diagonal-bound", type=int, default=20)
    plan.add_argument("--max-interval", type=int, default=200000)
    plan.add_argument("--padding", type=int, default=0)
    plan.add_argument("--minimum-coverage", type=float, default=0.95)
    plan.add_argument("--minimum-identity", type=float, default=0.5)
    plan.add_argument("--max-intron", type=int, default=20000)
    plan.add_argument("--genome-fallback", type=int, choices=(0, 1), default=1)
    plan.add_argument("--prediction-cache", type=Path, help="Frozen, checksum-verified previous rescue output; raw predictions are revalidated")
    plan.add_argument("--model-storage", choices=("compact", "legacy"), default="compact",
                      help="Versioned shared prediction storage, or explicit legacy JSON arrays")
    plan.add_argument("--retain-search-inputs", action="store_true",
                      help="Retain diagnostic and completed predictor FASTAs (default: regenerate on demand)")
    plan.add_argument("--max-genome-queries", type=int, default=20000)
    plan.add_argument("--unanchored-min-species", type=int, default=2)
    plan.add_argument("--terminal-max-extension", type=int, default=300)
    plan.add_argument("--terminal-max-unaligned-c-overhang", type=int, choices=(0, 1, 2), default=2)
    plan.add_argument("--gemoma-jar", type=Path, help="Optional GeMoMa refinement; requires Java and tblastn")
    plan.add_argument("--gemoma-java", default="java", help="Java executable compatible with the supplied GeMoMa jar")
    for name in ("synteny", "rescue", "finalize", "run", "status", "qc", "worker-complete", "qc-inputs", "export-search-inputs", "export-models"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--output", type=Path, required=True)
        cmd.add_argument("--cpus", type=int, default=1, help="Predictor and full-checksum CPU budget")
        if name in {"rescue", "run"}:
            cmd.add_argument("--interval-workers", type=int, help="Concurrent independent intervals (default: cpus); total threads stay within cpus")
        if name in {"synteny", "run"}:
            cmd.add_argument("--comparison-cache", type=Path, help="Shared comparison cache (default: OUTPUT parent/gene_model_rescue_comparison_cache)")
        if name in {"synteny", "rescue", "worker-complete", "qc-inputs"}:
            cmd.add_argument("--task-index", type=int, help="One-based frozen job/species index; omit to run all")
        if name in {"export-search-inputs", "export-models"}:
            cmd.add_argument("--task-index", type=int, required=True, help="One-based frozen species index")
            cmd.add_argument("--destination", type=Path, required=True, help="New export directory; producer files are untouched")
        if name == "export-search-inputs":
            cmd.add_argument("--combined", action="store_true", help="Also recreate repeated regions.fa/queries.fa diagnostics")
        if name == "qc":
            cmd.add_argument("--busco-dir", type=Path, required=True)
    return p


def main():
    p = parser()
    args = p.parse_args()
    if args.command == "self-align":
        if args.cpus < 1 or not 0 < args.cscore <= 1:
            p.error("Invalid self alignment parameters")
        align_self(args.cpus, args.cscore)
        return
    if args.command == "plan":
        if (args.common_references < 1 or args.nearest_references < 0 or not 0 <= args.minimum_busco <= 100
                or not 0 < args.cscore <= 1 or args.min_anchors < 2 or args.distance < 1
                or args.max_interval < 1 or args.padding < 0 or args.max_intron < 1 or args.diagonal_bound < 1
                or not 0 < args.minimum_coverage <= 1 or not 0 < args.minimum_identity <= 1
                or bool(args.feature) != bool(args.attribute)):
            p.error("Invalid rescue parameters")
        plan = build_plan(args)
        print(json.dumps({"species": len(plan["species"]), "comparisons": len(plan["synteny_jobs"]),
                          "common_references": plan["common_references"]}))
        return
    root = args.output.resolve()
    plan = load(root, inspect_only=args.command in {"export-search-inputs", "export-models"})
    if getattr(args, "cpus", 1) < 1:
        p.error("cpus must be positive")
    if getattr(args, "interval_workers", None) is not None and not 1 <= args.interval_workers <= args.cpus:
        p.error("interval-workers must be positive and no greater than cpus")
    if getattr(args, "task_index", None) is not None:
        count = len(plan["synteny_jobs"]) if args.command == "synteny" else len(plan["species"])
        if not 1 <= args.task_index <= count:
            p.error("task index outside frozen plan")
    if args.command == "export-search-inputs":
        print(export_worker_inputs(root, plan, plan["species"][args.task_index - 1], args.destination,
                                   combined=args.combined, cpus=args.cpus))
        return
    if args.command == "export-models":
        destination = export_destination(root, args.destination)
        name = plan["species"][args.task_index - 1]
        directory = root / "rescued" / name
        if not verified(directory, rescue_key(root, plan, name), hash_workers=args.cpus):
            raise ValueError("Legacy export requires a complete verified rescue")
        export_legacy_models(directory, destination, directory=True)
        print(destination)
        return
    if args.command == "qc-inputs":
        indices = [args.task_index] if args.task_index is not None else range(1, len(plan["species"]) + 1)
        for row in qc_work_items(root, plan, indices, hash_workers=args.cpus):
            print(*row, sep="\t")
    if args.command == "status":
        pending_pairs = []
        for job in plan["synteny_jobs"]:
            try:
                complete = verified(root / "synteny" / job["id"], comparison_key(root, job), hash_workers=args.cpus)
            except (OSError, ValueError):
                complete = False
            if not complete:
                pending_pairs.append(job["index"])
        pending_species = []
        for i, n in enumerate(plan["species"], 1):
            try:
                complete = verified(root / "rescued" / n, rescue_key(root, plan, n), hash_workers=args.cpus)
            except (OSError, ValueError):
                complete = False
            if not complete or pending_pairs:
                pending_species.append(i)
        print(json.dumps({"synteny_tasks": pending_pairs, "rescue_tasks": pending_species}))
    if args.command == "worker-complete":
        if not args.task_index or not 1 <= args.task_index <= len(plan["species"]):
            p.error("Worker index outside frozen plan")
        name = plan["species"][args.task_index - 1]
        if not verified(root / "rescued" / name, rescue_key(root, plan, name), hash_workers=args.cpus):
            raise ValueError("Rescue worker has no verified models")
        effective_key = {"plan": plan_digest(root, plan), "rescue_receipts": {
            name: digest(root / "rescued" / name / "receipt.json")}}
        if not verified(root / "effective" / name, effective_key, hash_workers=args.cpus):
            raise ValueError("Rescue worker has no verified exported inputs")
        files = [root / "rescued" / name / "receipt.json", root / "effective" / name / "receipt.json",
                 root / "qc/species_cds_busco_full" / (name + ".busco.full.tsv"),
                 root / "qc/species_cds_busco_short" / (name + ".busco.short.txt")]
        summary = files[-1]
        for directory in (root / "rescued" / name, root / "effective" / name):
            files += [directory / p for p in json.loads((directory / "receipt.json").read_text())["files"]]
        worker_dir = root / "workers" / name
        def current_files():
            hashes = hash_paths(root, files, workers=args.cpus, hash_function=digest)
            return {os.path.relpath(root / p, worker_dir): value for p, value in hashes.items()}
        frozen_files = current_files()
        quality = busco_quality(summary)
        initial = plan["request"]["sources"][name]["quality"]
        if any(quality[k] != initial[k] for k in COMPARABLE_QUALITY):
            raise ValueError("Worker BUSCO is not comparable to the initial run")
        effective = table(root / "effective" / name / "inputs.tsv")[0]
        if int(effective["rescued_models"]) == 0 and quality != initial:
            raise ValueError("Unchanged species BUSCO differs from the initial run")
        require_same_key(frozen_files, current_files())
        if (not verified(root / "rescued" / name, rescue_key(root, plan, name), hash_workers=args.cpus)
                or not verified(root / "effective" / name, effective_key, hash_workers=args.cpus)):
            raise ValueError("Worker dependencies changed during execution")
        atomic_json(worker_dir / "receipt.json", {
            "key": {"plan": plan_digest(root, plan), "species": name},
            "files": frozen_files})
    if args.command == "qc":
        finalize(root, plan, hash_workers=args.cpus)
        files = {n: args.busco_dir / (n + ".busco.short.txt") for n in plan["species"]}
        def current_key():
            return {"plan": plan_digest(root, plan), "augmented": digest(root / "augmented" / "receipt.json"),
                    "inputs": digest(root / "augmented" / "inputs.tsv"), "busco": digest_paths(files.values())}
        key = current_key()
        rows = table(root / "augmented" / "inputs.tsv")
        def build(tmp):
            report = []
            for row in rows:
                before = plan["request"]["sources"][row["species"]]["quality"]
                after = busco_quality(files[row["species"]])
                if any(before[k] != after[k] for k in COMPARABLE_QUALITY):
                    raise ValueError("Post-rescue BUSCO must use the initial lineage/version/mode/dataset date/marker count")
                report.append((row["species"], row["rescued_models"], before["complete_pct"], after["complete_pct"],
                               after["complete_pct"] - before["complete_pct"]))
            write_tsv(tmp / "before_after.tsv", ("species", "rescued_models", "busco_before_pct", "busco_after_pct", "delta_pct"), report)
        stage(root, Path("qc_report"), key, build, lambda: require_same_key(key, current_key()))
    if args.command in {"synteny", "run"}:
        indices = [args.task_index] if getattr(args, "task_index", None) else range(1, len(plan["synteny_jobs"]) + 1)
        for index in indices:
            synteny(root, plan, index, args.cpus, args.comparison_cache)
    if args.command in {"rescue", "run"}:
        indices = [args.task_index] if getattr(args, "task_index", None) else range(1, len(plan["species"]) + 1)
        for index in indices:
            if not 1 <= index <= len(plan["species"]):
                p.error("rescue task index outside frozen plan")
            rescue(root, plan, plan["species"][index - 1], args.cpus, args.interval_workers)
            name = plan["species"][index - 1]
            finalize(root, plan, [name], Path("effective") / name, hash_workers=args.cpus)
    if args.command in {"finalize", "run"}:
        finalize(root, plan, hash_workers=args.cpus)


if __name__ == "__main__":
    main()
