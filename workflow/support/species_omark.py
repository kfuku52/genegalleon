"""Species-wise OMAmer and OMArk results shared by both workflow stages."""
import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import os
import shutil
import subprocess
import tempfile
from contextlib import nullcontext
from pathlib import Path

from content_digest_cache import ContentDigestCache, cache, cached_sha256_file, configure
from fasta_sequence_store import exclusive_lock, fasta_records, write_json_atomic
from input_generation_array_state import safe_component

EMPTY_QUERY = hashlib.sha256(b"").hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def tool_identity(name):
    executable = shutil.which(name)
    if not executable:
        raise ValueError(f"Required program not found: {name}")
    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        version = None
    return [version, cached_sha256_file(Path(executable))[0]]


def query_bytes(records):
    """Canonical OMAmer input: order, wrapping and FASTA descriptions are immaterial."""
    proteins = {}
    for identifier, sequence in records:
        if not identifier or identifier in proteins or any(c.isspace() for c in identifier):
            raise ValueError(f"Empty, ambiguous or duplicated protein ID: {identifier!r}")
        sequence = sequence.upper().removesuffix("*")
        if not sequence:
            raise ValueError(f"Empty protein sequence: {identifier}")
        proteins[identifier] = sequence
    return b"".join(f">{key}\n{proteins[key]}\n".encode() for key in sorted(proteins))


def result_directory(cache_dir, name, contract, key):
    if not safe_component(name):
        raise ValueError(f"Invalid species name: {name!r}")
    return cache_dir / name / fingerprint([contract, key])


class SpeciesOMArk:
    def __init__(self, directory, database, scratch, keep_query=False, threads=1, taxonomy=None):
        self.directory = directory.resolve()
        configure(self.directory / "digests.sqlite")
        self.db = database.resolve(strict=True)
        self.db_stat = ContentDigestCache.signature(self.db)
        self.tools = {name: tool_identity(name) for name in ("omamer", "omark")}
        # Use the workspace ETE database when supplied, otherwise use OMArk's ETE3 default.
        workspace_taxonomy = os.environ.get("GG_TAXONOMY_DBFILE")
        if taxonomy is None and workspace_taxonomy and Path(workspace_taxonomy).is_file():
            taxonomy = workspace_taxonomy
        self.taxonomy = Path(taxonomy).resolve(strict=True) if taxonomy else Path.home() / ".etetoolkit" / "taxa.sqlite"
        # Array jobs can arrive together on the first database read.
        with exclusive_lock(self.directory / ".database.lock"):
            self.db_sha256 = cached_sha256_file(self.db)[0]
            self.bind_taxonomy()
        self.searches = 0
        self.threads = threads
        self.scratch = scratch
        self.keep_query = keep_query

    def taxonomy_signature(self):
        try:
            return ContentDigestCache.signature(self.taxonomy)
        except FileNotFoundError:
            return None

    def bind_taxonomy(self):
        self.taxonomy_stat = self.taxonomy_signature()
        self.taxonomy_sha256 = cached_sha256_file(self.taxonomy)[0] if self.taxonomy_stat else None
        self.contract = fingerprint([1, self.db_sha256, self.taxonomy_sha256, self.tools])

    def check_database(self):
        if ContentDigestCache.signature(self.db) != self.db_stat or self.taxonomy_signature() != self.taxonomy_stat:
            raise ValueError("OMArk database changed during evaluation")

    def result(self, name, key):
        return result_directory(self.directory, name, self.contract, key)

    def available(self, name, key):
        result = self.result(name, key)
        try:
            return ((result / f"{name}.omamer").stat().st_size > 0
                    and (result / "omark" / f"{name}.sum").stat().st_size > 0)
        except FileNotFoundError:
            return False

    def evaluate(self, name, data):
        key = hashlib.sha256(data).hexdigest()
        result = self.result(name, key)
        if not data:
            return EMPTY_QUERY
        if self.available(name, key):
            return key
        with exclusive_lock(result.parent / ".lock"):
            if self.taxonomy_stat is None and self.taxonomy.is_file():
                self.bind_taxonomy()
            if self.available(name, key):
                return key
            self.scratch.mkdir(parents=True, exist_ok=True)
            query = self.scratch / f"{name}.query.fa"
            if self.scratch.is_symlink() or query.is_symlink():
                raise ValueError(f"Refusing symlinked OMAmer query scratch: {query}")
            query.write_bytes(data)
            with tempfile.TemporaryDirectory(prefix=".run-", dir=result.parent) as tmp:
                stage = Path(tmp) / "result"
                out = stage / "omark"
                out.mkdir(parents=True)
                hits = stage / f"{name}.omamer"
                subprocess.run(["omamer", "search", "--db", self.db, "--query", query, "--out", hits,
                                "--nthreads", str(self.threads)], check=True)
                if not hits.stat().st_size:
                    raise ValueError(f"Empty OMAmer search result for {name}")
                # On the first run OMArk can initialise its own taxonomy DB.
                # Serialize just that case and publish under the resulting DB identity.
                initialising = self.taxonomy_stat is None
                with exclusive_lock(self.directory / ".database.lock") if initialising else nullcontext():
                    if initialising:
                        self.bind_taxonomy()
                    command = ["omark", "-f", hits, "-d", self.db, "-o", out]
                    if self.taxonomy_stat:
                        command.extend(["-e", self.taxonomy])
                    subprocess.run(command, check=True)
                    if self.taxonomy_stat is None:
                        self.bind_taxonomy()
                if not (out / f"{name}.sum").stat().st_size:
                    raise ValueError(f"Empty OMArk summary for {name}")
                self.check_database()
                result = self.result(name, key)
                if not self.available(name, key):
                    if result.exists():
                        result.rename(Path(tmp) / "incomplete")
                    stage.rename(result)
                self.searches += 1
            if not self.keep_query:
                query.unlink()
        return key

    def protein(self, name, path):
        path = path.resolve(strict=True)
        namespace = "omark_query_v1"
        saved = cache().get(namespace, path)
        if saved and (saved[0] == EMPTY_QUERY or self.available(name, saved[0])):
            return saved[0]
        before = ContentDigestCache.signature(path)
        data = query_bytes((i, seq) for i, _, seq in fasta_records(path))
        if ContentDigestCache.signature(path) != before:
            raise ValueError(f"Protein input changed during reading: {path}")
        key = self.evaluate(name, data)
        cache().put(namespace, path, key, len(data), expected_signature=before)
        return key


def link_result(source, target):
    if target.is_symlink() and target.resolve() == source:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".link-", dir=target.parent) as tmp:
        link = Path(tmp) / "link"
        link.symlink_to(os.path.relpath(source, target.parent))
        os.replace(link, target)


def genome(args, runner):
    sources = dict(args.input)
    if not sources or len(sources) != len(args.input) or not all(map(safe_component, sources)):
        raise ValueError("OMArk requires unique species inputs")
    output = args.genome
    state_file = output / ".omark-inputs.json"
    with exclusive_lock(output / ".omark.lock"):
        old = json.loads(state_file.read_text()) if state_file.exists() else {}
        if (args.cds_signature and old.get("cds_signature") == args.cds_signature
                and old.get("contract") == runner.contract and set(old["queries"]) == set(sources)
                and all(runner.available(name, key) for name, key in old["queries"].items())):
            queries = old["queries"]
            print("Skipped OMArk: matching species queries")
        elif args.cds_signature and not args.translated:
            return 75  # Ask the existing shell translation stage for its effective proteins.
        else:
            queries = {}
            for name, source in sources.items():
                protein = args.translated / (name + ".fa.gz") if args.translated else Path(source)
                queries[name] = runner.protein(name, protein)
                if old.get("contract") == runner.contract and old.get("queries", {}).get(name) == queries[name]:
                    print(f"Skipped OMArk: {name}")
        runner.check_database()
        if EMPTY_QUERY in queries.values():
            raise ValueError("Genome evolution OMArk requires nonempty protein inputs")
        for name, key in queries.items():
            result = runner.result(name, key)
            link_result(result / (name + ".omamer"), output / "omamer_search" / (name + ".omamer"))
            link_result(result / "omark", output / "omark" / name)
        for name in old.get("queries", {}).keys() - queries.keys():
            if not safe_component(name):
                raise ValueError("Invalid species in previous OMArk selection")
            for path in (output / "omark" / name, output / "omamer_search" / (name + ".omamer")):
                if path.is_symlink():
                    path.unlink()
        state = {"cds_signature": args.cds_signature, "contract": runner.contract, "queries": queries,
                 "cache": os.path.relpath(runner.directory, output)}
        if old != state:
            write_json_atomic(state_file, state)
    print(f"OMArk evaluated {runner.searches} new protein sets for {len(sources)} species.")
    return 0


def read_rows(path, expected):
    # These small published tables identify the protein files and rescue
    # boundaries. Check the same bytes that will be parsed.
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Published input table changed: {path}")
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))


def appended_ids(row):
    """Rescue finalize appends exactly one CDS per newly admitted locus."""
    if not int(row["rescued_models"]):
        return set()
    from fasta_sequence_store import open_text
    found = set()
    count = 0
    original = int(row["original_cds"])
    with open_text(Path(row["cds"])) as handle:
        for line in handle:
            if line.startswith(">"):
                if count >= original:
                    found.add(line[1:].split()[0])
                count += 1
    if len(found) != int(row["rescued_models"]) or count != original + len(found):
        raise ValueError(f"Rescue CDS and published counts disagree: {row['species']}")
    return found


def rescued_ids(locus, added):
    aliases = {locus.get("source_gene_id")}
    for candidate in locus["candidates"]:
        aliases.add(candidate.get("source_transcript_id"))
        aliases.update(candidate.get("source_fasta_ids", ()))
    return aliases & added


def baseline_keys(runner, name, catalog, augmentation, assembly, input_only=False):
    from gene_model_selection import _baseline

    phases = ("input", "rescued") if augmentation and not input_only else ("input",)
    # Only the two consumed inputs can change baseline membership or sequence.
    augmentation_id = None
    if augmentation:
        augmentation_id = [augmentation["original_cds"], augmentation["rescued_models"]]
        if int(augmentation["rescued_models"]):
            augmentation_id.append(ContentDigestCache.signature(Path(augmentation["cds"])))
    prefix = "omark_baseline_v1_" + fingerprint([augmentation_id, assembly])
    saved = {p: cache().get(prefix + p, catalog) for p in phases}
    if all(s and (s[0] == EMPTY_QUERY or runner.available(name, s[0])) for s in saved.values()):
        return {p: s[0] for p, s in saved.items()}

    before = ContentDigestCache.signature(catalog)
    added = appended_ids(augmentation) if augmentation else set()
    if added and ContentDigestCache.signature(Path(augmentation["cds"])) != augmentation_id[-1]:
        raise ValueError(f"Rescue CDS changed while extracting added loci: {name}")
    seen = set()
    records = []
    with catalog.open() as handle:
        for line in handle:
            locus = json.loads(line)
            candidates = locus["candidates"]
            candidate = _baseline(locus, candidates)
            new = set()
            if added:
                new = rescued_ids(locus, added)
                if len(new) > 1 or new & seen:
                    raise ValueError(f"Ambiguous rescue membership at {locus['gene_id']}")
                seen.update(new)
            quality = candidate["quality"]
            if (quality["usable"] and candidate["protein"] and "*" not in candidate["protein"]
                    and (assembly or (candidate["blocks"] and not locus.get("ambiguous_coordinates")
                         and not quality.get("structure_problem") and not quality.get("sequence_mismatch")))):
                records.append((locus["gene_id"], candidate["protein"], bool(new)))
    missing = added - seen
    if missing:
        excluded = catalog.with_name("excluded_loci.jsonl")
        if excluded.is_file():
            with excluded.open() as handle:
                for line in handle:
                    missing -= rescued_ids(json.loads(line), missing)
                    if not missing:
                        break
        if missing:
            raise ValueError(f"Unmapped newly rescued CDS in {name}: {sorted(missing)[:5]}")
    if ContentDigestCache.signature(catalog) != before:
        raise ValueError(f"Refinement catalogue changed: {catalog}")
    keys = {}
    for phase in phases:
        if phase == "rescued" and not added:
            keys[phase] = keys["input"]
        else:
            data = query_bytes((gene, seq) for gene, seq, new in records if phase == "rescued" or not new)
            keys[phase] = runner.evaluate(name, data)
        cache().put(prefix + phase, catalog, keys[phase], len(records), expected_signature=before)
    return keys


def comparison_fields(summary):
    if summary is None:
        return dict(representative_proteins=0, consistent_full=0, consistent_full_pct="", fragment=0,
                    partial=0, selected_clade="", conserved_hogs="", missing_hogs="")
    full = int(summary["consistent_count"] - summary["consistent_partial_count"] - summary["consistent_fragment_count"])
    total = int(summary["protein_count"])
    return dict(representative_proteins=total, consistent_full=full, consistent_full_pct=100 * full / total,
                fragment=sum(int(summary[f"{g}_fragment_count"]) for g in ("consistent", "inconsistent", "contamination")),
                partial=sum(int(summary[f"{g}_partial_count"]) for g in ("consistent", "inconsistent", "contamination")),
                selected_clade=summary["selected_clade"], conserved_hogs=int(summary["conserved_hogs"]),
                missing_hogs=int(summary["missing_count"]))


def augmentation_rows(request):
    if request.get("rescue_output"):
        manifest = Path(request["rescue_output"]) / "augmented/inputs.tsv"
        # Rescue can also have supplied a guide tree only.
        if str(manifest) in request.get("files", {}):
            return {r["species"]: r for r in read_rows(manifest, request["files"][str(manifest)])}
    return None


def catalog_keys(runner, root, request, augmented, name, input_only=False):
    if not safe_component(name):
        raise ValueError(f"Invalid species name: {name!r}")
    assembly = name in request.get("transcriptome_sources", {})
    rescue = augmented.get(name) if augmented is not None and not assembly else None
    if augmented is not None and not assembly and not rescue:
        raise ValueError(f"Missing rescue augmentation for {name}")
    if rescue and Path(request["sources"][name]["fasta"]).resolve() != Path(rescue["cds"]).resolve():
        raise ValueError(f"Rescue publication differs from refinement input: {name}")
    return baseline_keys(runner, name, root / "catalog" / name / "loci.jsonl", rescue, assembly, input_only)


def baseline(args, runner):
    root = args.baseline.resolve()
    plan = json.loads((root / "plan.json").read_bytes())
    species = plan["request"].get("target_species", plan["species"])
    if not 1 <= args.task_index <= len(species):
        raise ValueError(f"Baseline task index {args.task_index} is outside 1..{len(species)}")
    name = species[args.task_index - 1]
    request = plan["request"]
    keys = catalog_keys(runner, root, request, augmentation_rows(request), name, input_only=True)
    runner.check_database()
    print(f"{name} input protein_set_sha256={keys['input']} new_searches={runner.searches}")


def refinement(args, runner):
    root = args.refinement.resolve()
    from refinement_figure_inputs import locus_totals
    from summarize_omark import parse_sum_file

    plan_data = (root / "plan.json").read_bytes()
    plan = json.loads(plan_data)
    request = plan["request"]
    effective = root / "effective"
    receipt = effective / "receipt.json"
    generation = ContentDigestCache.signature(receipt)
    receipt_data = receipt.read_bytes()
    publication = json.loads(receipt_data)
    if publication["key"]["plan"] != hashlib.sha256(plan_data).hexdigest():
        raise ValueError(f"Completed effective publication belongs to a different refinement plan: {root}")
    mixed = bool(request.get("transcriptome_sources"))
    table = effective / "sequence/inputs.tsv" if mixed else effective / "inputs.tsv"
    augmented = augmentation_rows(request)
    original_sources = {}
    if augmented:
        rescue_plan = Path(request["rescue_output"]) / "plan.json"
        rescue_data = rescue_plan.read_bytes()
        expected = request["files"][str(rescue_plan)]
        if hashlib.sha256(rescue_data).hexdigest() != expected:
            raise ValueError(f"Frozen rescue plan changed: {rescue_plan}")
        original_sources = json.loads(rescue_data)["request"]["sources"]
    genomic = {}
    if mixed and "inputs.tsv" in publication["files"]:
        genomic = {row["species"]: row for row in read_rows(effective / "inputs.tsv", publication["files"]["inputs.tsv"])}
    rows = []
    inputs = read_rows(table, publication["files"][str(table.relative_to(effective))])
    species = [row["species"] for row in inputs]
    if not species or len(set(species)) != len(inputs) or set(species) != set(request.get("target_species", plan["species"])):
        raise ValueError("Final OMArk input species differ from the refinement plan")
    publication_hash = hashlib.sha256(receipt_data).hexdigest()
    for entry in inputs:
        name = entry["species"]
        keys = catalog_keys(runner, root, request, augmented, name)
        keys["refined"] = runner.protein(name, table.parent / entry["protein"])
        counts = locus_totals(name, request, original_sources, request.get("sources", {}).get(name),
                              genomic.get(name, entry), table.parent)
        summaries = {}
        for phase, key in keys.items():
            result = runner.result(name, key) / "omark" / (name + ".sum") if key != EMPTY_QUERY else None
            if key not in summaries:
                summaries[key] = comparison_fields(parse_sum_file(result) if result else None)
            stats = summaries[key]
            baseline = summaries[keys["input"]]
            comparable = (stats["selected_clade"] not in ("", "-") and stats["selected_clade"] == baseline["selected_clade"]
                          and stats["conserved_hogs"] == baseline["conserved_hogs"])
            rows.append({"species": name, "stage": phase, **stats, "hogs_comparable_to_input": int(comparable),
                         "delta_missing_hogs_from_input": stats["missing_hogs"] - baseline["missing_hogs"] if comparable else "",
                         "annotated_gene_loci": counts.get(phase) or "",
                         "refinement_plan_sha256": publication["key"]["plan"],
                         "effective_receipt_sha256": publication_hash,
                         "protein_set_sha256": key, "database_sha256": runner.db_sha256,
                         "taxonomy_sha256": runner.taxonomy_sha256 or "",
                         **{name + "_version": version or sha for name, (version, sha) in runner.tools.items()},
                         "sum_file": str(result) if result else ""})
    runner.check_database()
    if ContentDigestCache.signature(receipt) != generation:
        raise ValueError("Effective publication changed during OMArk comparison")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.report.parent) as tmp:
        temp = Path(tmp) / "comparison.tsv"
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, rows[0], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        temp.replace(args.report)
    print(f"OMArk evaluated {runner.searches} new protein sets for {len(inputs)} species. Comparison: {args.report}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--genome", type=Path, help="Genome evolution output directory")
    target.add_argument("--refinement", type=Path, help="Completed refinement publication")
    target.add_argument("--baseline", type=Path, help="Refinement plan with a completed species catalogue")
    parser.add_argument("--task-index", type=int, help="One-based species index in the baseline plan")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, help="Shared ETE taxonomy database (default: GG_TAXONOMY_DBFILE or ETE3 default)")
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--keep-query", action="store_true")
    parser.add_argument("--threads", type=int, default=1, help="Allocated OMAmer search CPUs")
    parser.add_argument("--input", nargs=2, action="append", metavar=("SPECIES", "FASTA"))
    parser.add_argument("--cds-signature", help="Effective CDS input and translation policy identity")
    parser.add_argument("--translated", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if bool(args.baseline) != (args.task_index is not None):
        parser.error("--baseline requires --task-index, which is exclusive to baseline evaluation")
    if args.cds_signature and not args.genome:
        parser.error("--cds-signature requires --genome")
    if args.genome and not args.input:
        parser.error("--genome requires species inputs")
    if args.refinement and not args.report:
        parser.error("--refinement requires --report")
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.cds_signature:
        if not args.translated and not (args.genome / ".omark-inputs.json").exists():
            return 75
        configure(args.cache.resolve() / "digests.sqlite")
        # Use the same translation policy as prepare_species_protein_tmp. A
        # changed translator may still produce an existing canonical query.
        args.cds_signature = fingerprint([1, args.cds_signature, tool_identity("cdskit"), tool_identity("seqkit")])
    runner = SpeciesOMArk(args.cache, args.db, args.scratch, args.keep_query, args.threads, args.taxonomy)
    if args.baseline:
        return baseline(args, runner)
    return genome(args, runner) if args.genome else refinement(args, runner)


if __name__ == "__main__":
    raise SystemExit(main())
