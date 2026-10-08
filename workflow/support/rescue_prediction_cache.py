"""Verified predictor reuse, independently of old acceptance decisions.

The small frozen key is cheap to compare repeatedly. Verification reads source,
query and result contents before a consumer may reuse predictions. JSON arrays
are streamed so a gigabyte result does not require a gigabyte temporary string.
"""
import csv
import hashlib
import json
import math
import os
import stat
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

try:
    from fasta_sequence_store import fasta_records
    from input_generation_array_state import FreshDigestBatch
except ImportError:
    from .fasta_sequence_store import fasta_records
    from .input_generation_array_state import FreshDigestBatch


SEARCH_CONTRACT_FILE = "prediction_search_contract.json"
MINIPROT_OUTPUT_SCORE_RATIO = 0.5
MINIPROT_MAX_SECONDARY = 30


def prediction_search_contract():
    """Effective output bounds; other defaults are bound to the miniprot SHA."""
    scope = {"output_score_ratio": MINIPROT_OUTPUT_SCORE_RATIO, "max_secondary": MINIPROT_MAX_SECONDARY}
    return {"schema": 1, "local": dict(scope), "genome": dict(scope)}


# Repository-owned, archived source schemas, not upstream version pins. The
# .172 source was verified byte-for-byte against its frozen producer-plan SHA.
# .189 added verified inherited predictions; its cache-verifier SHA is bound
# below, and its parent plan/worker receipt must prove the inherited contract.
# Both explicitly searched local intervals at .5 and genome fallback at .99.
LEGACY_DIRECT_IMPLEMENTATION = "608902a388c9c7a9530fce9371e01e7cbed24fa8428bca3e3bb9c457cf699e98"
LEGACY_INHERITED_IMPLEMENTATION = "db0e6657b081331e0cdfef84e66b553517f5dfa3e87cdccf0fc7e3e75458a6d7"
LEGACY_INHERITED_VERIFIER = "339c210589a679560520405d83e0fdd8fa27e1d11fd08a66bea45382767be111"
# Those sources left -N at the executable default. The independently audited
# historical executable has N=30; unknown executables cannot inherit that fact.
LEGACY_MINIPROT_SHA256 = "f8822f41eceb53a6ca611bd6606faf035aeb913030d1dede7fd7889ab37d4f52"


def _validated_search_contract(value):
    if (not isinstance(value, dict) or set(value) != {"schema", "local", "genome"}
            or type(value["schema"]) is not int or value["schema"] != 1):
        raise ValueError("Malformed prediction search contract")
    for kind in ("local", "genome"):
        scope = value[kind]
        if (not isinstance(scope, dict) or set(scope) != {"output_score_ratio", "max_secondary"}
                or isinstance(scope["output_score_ratio"], bool)
                or not isinstance(scope["output_score_ratio"], (int, float))
                or not math.isfinite(scope["output_score_ratio"]) or not 0 < scope["output_score_ratio"] <= 1
                or isinstance(scope["max_secondary"], bool) or not isinstance(scope["max_secondary"], int)
                or scope["max_secondary"] < 1):
            raise ValueError("Malformed prediction search contract scope")
    return value


def _verified_search_contract(root, plan, name, receipt, expect, ancestors=None):
    """Prove legacy inherited flags using only frozen metadata, without logs.

    Unknown implementations abstain. A forged or changed declared ancestor is
    an integrity error, not a silently ignored cache miss. No ancestor models,
    commands or scratch inputs are consumed by this contract proof.
    """
    root = Path(root)
    ancestors = set() if ancestors is None else ancestors
    token = (str(root.resolve(strict=True)), name)
    if token in ancestors or len(ancestors) >= 64:
        raise ValueError("Cyclic or excessive prediction cache ancestry")
    ancestors = ancestors | {token}
    tools = plan["request"]["tools"]
    if SEARCH_CONTRACT_FILE in receipt["files"]:
        path = _safe_member(root / "rescued" / name, SEARCH_CONTRACT_FILE)
        expect(path, receipt["files"][SEARCH_CONTRACT_FILE])
        declared = _validated_search_contract(json.loads(path.read_text()))
        if tools.get("prediction_search_contract") != declared:
            raise ValueError("Prediction search contract is not bound to its producer plan")
        return declared
    if "prediction_search_contract" in tools:
        raise ValueError("Prediction producer omitted its frozen search contract")
    implementation = tools.get("implementation")
    if (implementation not in {LEGACY_DIRECT_IMPLEMENTATION, LEGACY_INHERITED_IMPLEMENTATION}
            or tools.get("miniprot_sha256") != LEGACY_MINIPROT_SHA256):
        return {"schema": 1, "local": None, "genome": None}
    result = {"schema": 1, "local": {"output_score_ratio": 0.5, "max_secondary": 30},
              "genome": {"output_score_ratio": 0.99, "max_secondary": 30}}
    parent = plan["request"].get("prediction_cache")
    if not parent:
        return result
    if implementation != LEGACY_INHERITED_IMPLEMENTATION or tools.get("verify_prediction_cache_implementation") != LEGACY_INHERITED_VERIFIER:
        raise ValueError("Legacy inherited prediction verifier is unproven")
    if (not isinstance(parent, dict) or parent.get("schema") != 1 or not isinstance(parent.get("species"), dict)
            or not isinstance(parent.get("root"), str)):
        raise ValueError("Malformed frozen prediction ancestor")
    if name not in parent["species"]:
        return result  # This worker had no inherited producer to reuse.
    parent_root = Path(parent["root"])
    if str(parent_root.resolve(strict=True)) != parent["root"]:
        raise ValueError("Frozen prediction ancestor root differs")
    expect(parent_root / "plan.json", parent.get("plan_sha256"))
    parent_plan = json.loads((parent_root / "plan.json").read_text())
    if name not in parent_plan["species"]:
        raise ValueError("Prediction ancestor lacks its frozen worker species")
    directory = parent_root / "rescued" / name
    frozen = parent["species"][name]
    if not isinstance(frozen, dict) or not isinstance(frozen.get("files"), dict):
        raise ValueError("Malformed frozen ancestor worker")
    if not {"models.json", "candidates.json"} <= set(frozen["files"]):
        raise ValueError("Frozen prediction ancestor lacks result digests")
    for member in frozen["files"]:
        _safe_member(directory, member, resolve=False)
    expect(directory / "receipt.json", frozen.get("receipt_sha256"))
    parent_receipt = _receipt(directory)
    if (parent_receipt["key"].get("plan") != parent["plan_sha256"] or parent_receipt["key"].get("species") != name
            or any(parent_receipt["files"].get(member) != digest for member, digest in frozen["files"].items())):
        raise ValueError("Prediction ancestor worker differs from its frozen generation")
    inherited = _verified_search_contract(parent_root, parent_plan, name, parent_receipt, expect, ancestors)
    for kind in ("local", "genome"):
        if inherited[kind] != result[kind]:
            result[kind] = None  # Mixed or unknown search scopes are not complete coverage.
    return result


def _safe_member(directory, member, *, resolve=True):
    if not isinstance(member, str):
        raise ValueError("Unsafe prediction cache member")
    relative = PurePosixPath(member)
    if (not member or relative.is_absolute()
            or ".." in relative.parts or "\\" in member or any(ord(c) < 32 for c in member)):
        raise ValueError("Unsafe prediction cache member")
    path = Path(directory) / member
    if resolve and not path.resolve(strict=True).is_relative_to(Path(directory).resolve(strict=True)):
        raise ValueError("Prediction cache member escapes its producer directory")
    return path


def _json_snapshot(path):
    """One content read supplies both JSON and SHA, fenced against replacement."""
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("Prediction metadata must be a regular file")
        raw = handle.read()
        after = os.fstat(handle.fileno())
        current = os.stat(path)
    def identity(info):
        return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns
    if identity(before) != identity(after) or identity(before) != identity(current):
        raise OSError("Prediction metadata changed while reading: " + str(path))
    return json.loads(raw), hashlib.sha256(raw).hexdigest(), identity(current)


def _receipt(directory, data=None):
    if data is None:
        data, _, _ = _json_snapshot(Path(directory) / "receipt.json")
    if not isinstance(data, dict) or not isinstance(data.get("key"), dict) or not isinstance(data.get("files"), dict):
        raise ValueError("Malformed prediction producer receipt")
    for member, value in data["files"].items():
        # Only consumed prediction/prepared members need filesystem resolution.
        # A producer may have millions of interval scratch paths; their bytes
        # are not used by prediction reuse and must not cause metadata storms.
        _safe_member(directory, member, resolve=False)
        if not isinstance(value, str) or len(value) != 64 or set(value) - set("0123456789abcdef"):
            raise ValueError("Invalid prediction receipt digest")
    return data


def frozen_prediction_cache_key(root, names=None):
    """Capture immutable plan/worker receipt digests without hashing large models."""
    root = Path(root).resolve(strict=True)
    plan_path = root / "plan.json"
    plan, plan_hash, plan_identity = _json_snapshot(plan_path)
    selected = list(plan["species"] if names is None else names)
    species = {}
    for name in selected:
        if name not in plan["species"] or "/" in name or "\\" in name or name.startswith("."):
            raise ValueError("Unsafe or unknown prediction cache species")
        directory = root / "rescued" / name
        if not (directory / "receipt.json").exists():
            continue
        receipt_data, receipt_hash, _ = _json_snapshot(directory / "receipt.json")
        receipt = _receipt(directory, receipt_data)
        if receipt["key"].get("species") != name or receipt["key"].get("plan") != plan_hash:
            raise ValueError("Prediction worker receipt is not bound to its original plan")
        required = ("models.json", "candidates.json")
        if any(member not in receipt["files"] for member in required):
            raise ValueError("Prediction producer receipt lacks required result files")
        selected_files = {member: receipt["files"][member] for member in required}
        if "genome_query_mapping.tsv" in receipt["files"]:
            selected_files["genome_query_mapping.tsv"] = receipt["files"]["genome_query_mapping.tsv"]
        if SEARCH_CONTRACT_FILE in receipt["files"]:
            selected_files[SEARCH_CONTRACT_FILE] = receipt["files"][SEARCH_CONTRACT_FILE]
        species[name] = {"receipt_sha256": receipt_hash, "files": selected_files}
    current = plan_path.stat()
    if plan_identity != (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns, current.st_ctime_ns):
        raise OSError("Prediction plan changed while capturing worker metadata")
    return {"schema": 1, "root": str(root), "plan_sha256": plan_hash, "species": species}


def stream_json_array(path, *, chunk_size=1024 * 1024, max_record_bytes=64 * 1024 * 1024, hasher=None):
    """Read a strict JSON array with bounded buffering and no extra dependency."""
    def invalid_constant(value):
        raise ValueError("Nonfinite value in prediction cache JSON: " + value)
    def unique_keys(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate key in prediction cache JSON: " + key)
            result[key] = value
        return result
    decoder = json.JSONDecoder(parse_constant=invalid_constant, object_pairs_hook=unique_keys)
    with Path(path).open(encoding="utf-8", newline="") as handle:
        buffer, position, eof = "", 0, False
        def fill():
            nonlocal buffer, position, eof
            buffer = buffer[position:]
            position = 0
            block = handle.read(chunk_size)
            if hasher is not None:
                hasher.update(block.encode("utf-8"))
            eof = not block
            buffer += block
        def token():
            nonlocal position
            while True:
                while position < len(buffer) and buffer[position].isspace():
                    position += 1
                if position < len(buffer) or eof:
                    return buffer[position:position + 1]
                fill()
        if token() != "[":
            raise ValueError("Prediction cache requires a JSON array")
        position += 1
        first = True
        while True:
            next_token = token()
            if next_token == "]":
                position += 1
                break
            if not first:
                if next_token != ",":
                    raise ValueError("Invalid prediction array separator")
                position += 1
                if token() == "]":
                    raise ValueError("Trailing prediction array comma")
            elif not next_token:
                raise ValueError("Truncated prediction JSON array")
            first = False
            while True:
                try:
                    value, end = decoder.raw_decode(buffer, position)
                    position = end
                    break
                except json.JSONDecodeError as error:
                    if eof:
                        raise ValueError("Malformed or truncated prediction record") from error
                    if len(buffer) - position > max_record_bytes:
                        raise ValueError("Prediction cache record exceeds bounded buffer") from error
                    fill()
                    token()
            yield value
        if token():
            raise ValueError("Trailing data after prediction JSON array")


RAW_FIELDS = ("coverage", "query_span_coverage", "query_start", "query_end", "query_length", "frameshift",
              "paf_query", "paf_seqid", "paf_strand", "paf", "query", "seqid", "strand", "identity", "cds", "id", "search")


def _raw_prediction(model, regions):
    if not isinstance(model, dict) or model.get("query") not in regions:
        return None
    completion = model.get("terminal_completion")
    if ("raw_prediction" not in model and isinstance(completion, dict)
            and (completion.get("status") in {"completed", "accepted"} or "selected_cds" in completion)):
        raise ValueError("Completed terminal model lacks pristine cached prediction")
    if "raw_prediction" in model:
        original = model["raw_prediction"]
        if (not isinstance(original, dict) or any(original.get(key) != model.get(key)
                                                for key in ("query", "evidence", "search", "seqid", "strand"))):
            raise ValueError("Cached raw prediction does not match its model provenance")
        model = original
    region = regions[model["query"]]
    old_evidence = model.get("evidence", {})
    if old_evidence != region:
        raise ValueError("Cached prediction candidate evidence changed")
    if model.get("search") not in {"synteny_interval", "genome_fallback"}:
        raise ValueError("Unsupported cached predictor search")
    if (not isinstance(model.get("seqid"), str) or not model["seqid"] or any(c.isspace() for c in model["seqid"])
            or model.get("strand") not in {"+", "-"} or not isinstance(model.get("frameshift"), bool)
            or not isinstance(model.get("cds"), list) or not model["cds"]):
        raise ValueError("Malformed cached genomic prediction")
    for exon in model["cds"]:
        if (not isinstance(exon, list) or len(exon) != 3 or any(isinstance(x, bool) or not isinstance(x, int) for x in exon)
                or not 0 <= exon[0] < exon[1] or exon[2] not in {0, 1, 2}):
            raise ValueError("Invalid cached prediction CDS")
    for key in ("coverage", "identity"):
        if (isinstance(model.get(key), bool) or not isinstance(model.get(key), (int, float))
                or not math.isfinite(model[key]) or not 0 <= model[key] <= 1):
            raise ValueError("Invalid cached prediction alignment metric")
    result = {key: model[key] for key in RAW_FIELDS if key in model}
    result["cds"] = [list(exon) for exon in model["cds"]]
    result["evidence"] = region
    result["prediction_cache"] = True
    # No old sequence, QC, support, selection, model ID or acceptance is reused.
    return result


@dataclass
class VerifiedPredictionCache:
    directory: Path
    regions: dict
    frozen: dict
    batch: FreshDigestBatch
    old_root: Path = None
    local_search_compatible: bool = True
    genome_search_compatible: bool = True
    search_contract: dict = field(default_factory=dict)
    _candidate_ids: set = None
    _old_genome_regions: dict = field(default_factory=dict)
    _genome_mapping: dict = field(default_factory=dict)
    _genome_bindings: dict = field(default_factory=dict)
    genome_reuse: dict = field(default_factory=dict)

    def check(self):
        self.batch.check()

    def iter_models(self):
        self.check()
        hasher = hashlib.sha256()
        for model in stream_json_array(self.directory / "models.json", hasher=hasher):
            if isinstance(model, dict) and model.get("search") == "genome_fallback":
                if not self.genome_search_compatible:
                    continue
                if model.get("query") not in self._old_genome_regions:
                    raise ValueError("Cached genome prediction lacks searched-query provenance")
                bindings = self._genome_bindings.get(model.get("query"), ())
                if not bindings:
                    continue  # Expanded old duplicates are supplied once below.
                original = _raw_prediction(model, self._old_genome_regions)
                for candidate in bindings:
                    row = {**original, "cds": [list(exon) for exon in original["cds"]],
                           "query": candidate, "evidence": self.regions[candidate]}
                    if "paf_query" in row:
                        row["paf_query"] = candidate
                    if "paf" in row:
                        paf = row["paf"].split("\t")
                        if len(paf) < 13 or paf[:2] != ["##PAF", model["query"]]:
                            raise ValueError("Cached genome prediction PAF does not match its query")
                        paf[1] = candidate
                        row["paf"] = "\t".join(paf)
                    if "id" in row:
                        row["id"] = "cached_" + candidate + "_" + hashlib.sha256(str(row["id"]).encode()).hexdigest()[:12]
                    row["prediction_cache_source_query"] = model["query"]
                    yield row
            else:
                if not self.local_search_compatible:
                    continue
                raw = _raw_prediction(model, self.regions)
                if raw is not None:
                    yield raw
        if hasher.hexdigest() != self.frozen["files"]["models.json"]:
            raise ValueError("Prediction model content changed after verification")
        self.check()

    def candidate_ids(self):
        """Queries with no old alignment are reused as verified empty results."""
        self.check()
        if self._candidate_ids is not None:
            return set(self._candidate_ids)
        ids = set()
        seen = set()
        self._read_genome_mapping()
        for row in stream_json_array(self.directory / "candidates.json"):
            if not isinstance(row, dict) or not isinstance(row.get("id"), str):
                raise ValueError("Malformed cached candidate")
            if row["id"] in seen:
                raise ValueError("Duplicate cached candidate")
            seen.add(row["id"])
            if row["id"] in self.regions:
                if row != self.regions[row["id"]]:
                    raise ValueError("Cached candidate definition changed")
                if self.local_search_compatible and not row.get("genome_only"):
                    ids.add(row["id"])
            if row["id"] in self._genome_mapping:
                if row["id"] in self._old_genome_regions:
                    raise ValueError("Duplicate cached genome candidate")
                self._old_genome_regions[row["id"]] = row
        if set(self._old_genome_regions) != set(self._genome_mapping):
            raise ValueError("Cached genome mapping references absent candidate evidence")
        self._candidate_ids = ids
        self.check()
        return set(ids)

    def genome_candidate_ids(self):
        """Old genome-search coverage, independently of its acceptance outcome."""
        self.check()
        if not self.genome_search_compatible:
            return set()
        if self.genome_reuse:
            return {candidate for bindings in self._genome_bindings.values() for candidate in bindings}
        return set(self._read_genome_mapping()) & set(self.regions)

    def genome_query_mapping(self):
        """Current searched representatives, suitable for transitive receipts."""
        self.check()
        result = {}
        for ids in self._genome_bindings.values():
            representative = min(ids)
            result.update({candidate: representative for candidate in ids})
        return result

    def _read_genome_mapping(self):
        if self._genome_mapping or "genome_query_mapping.tsv" not in self.frozen["files"]:
            return self._genome_mapping
        with (self.directory / "genome_query_mapping.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["candidate", "representative"]:
                raise ValueError("Malformed cached genome query mapping")
            for row in reader:
                if row["candidate"] in self._genome_mapping or not row["representative"]:
                    raise ValueError("Duplicate or malformed cached genome query")
                self._genome_mapping[row["candidate"]] = row["representative"]
        if any(self._genome_mapping.get(rep) != rep for rep in self._genome_mapping.values()):
            raise ValueError("Cached genome representative lacks its own searched query")
        self.check()
        return self._genome_mapping

    def bind_exact_genome_sequences(self):
        """Rebind only genome searches whose full prepared protein is identical."""
        self.candidate_ids()
        if not self.genome_search_compatible:
            self.genome_reuse = {"covered_current_queries": 0, "covered_unique_sequences": 0,
                                 "new_genome_only_queries": 0, "policy": "search_contract_mismatch"}
            return
        required = {}
        for row in [*self._old_genome_regions.values(), *self.regions.values()]:
            required.setdefault(row["donor"], set()).add(row["query"])
        sequences = {}
        for donor, genes in required.items():
            for identifier, _, sequence in fasta_records(self.old_root / "prepared" / donor / "genes.pep"):
                if identifier in genes:
                    sequences[donor, identifier] = sequence
            if any((donor, gene) not in sequences for gene in genes):
                raise ValueError("Cached genome query lacks verified prepared protein")
        sources = {}
        for identifier, row in self._old_genome_regions.items():
            sequence = sequences[row["donor"], row["query"]]
            representative = self._genome_mapping[identifier]
            source_row = self._old_genome_regions[representative]
            if sequences[source_row["donor"], source_row["query"]] != sequence:
                raise ValueError("Cached genome representative protein differs")
            sources[sequence] = representative
        for candidate, row in self.regions.items():
            source = sources.get(sequences[row["donor"], row["query"]])
            if source is not None:
                self._genome_bindings.setdefault(source, []).append(candidate)
        self.genome_reuse = {"covered_current_queries": sum(map(len, self._genome_bindings.values())),
                             "covered_unique_sequences": len(self._genome_bindings),
                             "new_genome_only_queries": sum(self.regions[c]["genome_only"] for ids in self._genome_bindings.values()
                                                            for c in ids if self.regions[c].get("genome_only")),
                             "policy": "exact_prepared_protein_identity_genome_search_only"}
        self.check()


def verify_prediction_cache(old_root, new_root, plan, name, regions, params, frozen=None):
    """Verify compatible prediction inputs once; changed generation must fail."""
    old_root, new_root = Path(old_root), Path(new_root)
    frozen = frozen if frozen is not None else frozen_prediction_cache_key(old_root, [name])
    if frozen.get("schema") != 1 or str(old_root.resolve(strict=True)) != frozen.get("root"):
        raise ValueError("Frozen prediction cache root differs")
    if name not in frozen["species"]:
        return None
    batch = FreshDigestBatch()
    def expect(path, expected):
        if batch.read([path])[str(path)] != expected:
            raise ValueError("Frozen prediction cache content changed: " + str(path))
    expect(old_root / "plan.json", frozen["plan_sha256"])
    old_plan = json.loads((old_root / "plan.json").read_text())
    directory = old_root / "rescued" / name
    worker = frozen["species"][name]
    expect(directory / "receipt.json", worker["receipt_sha256"])
    receipt = _receipt(directory)
    if (receipt["key"].get("plan") != frozen["plan_sha256"] or receipt["key"].get("species") != name
            or any(receipt["files"].get(p) != v for p, v in worker["files"].items())):
        raise ValueError("Cached worker is not bound to frozen plan/results")
    for member, expected in worker["files"].items():
        expect(_safe_member(directory, member), expected)
    for tool in ("miniprot", "miniprot_sha256"):
        if not old_plan["request"]["tools"].get(tool) or old_plan["request"]["tools"][tool] != plan["request"]["tools"].get(tool):
            raise ValueError("Prediction cache miniprot identity changed")
    try:
        from gene_model_species_profiles import parameters_for
    except ImportError:
        from .gene_model_species_profiles import parameters_for
    old_params = parameters_for(old_plan["request"], name)
    # All search/nominated-region parameters must match. Acceptance/QC changes
    # are permitted because the caller recomputes them against the genome.
    for key in ("max_intron", "genome_fallback", "max_interval", "padding", "cscore", "min_anchors", "distance", "diagonal_bound"):
        if old_params.get(key) != params.get(key):
            raise ValueError("Prediction cache search parameter changed: " + key)
    donors = {name, *(row["donor"] for row in regions), *receipt["key"].get("prepared", {})}
    for donor in sorted(donors):
        old_source = old_plan["request"]["sources"].get(donor)
        source = plan["request"]["sources"].get(donor)
        if not old_source or not source or old_source["genetic_code"] != source["genetic_code"]:
            raise ValueError("Prediction cache source or genetic code changed")
        for source_field in ("genome", "fasta", "gff"):
            old_hash = old_plan["request"]["files"][old_source[source_field]]
            if plan["request"]["files"][source[source_field]] != old_hash:
                raise ValueError("Prediction cache frozen source changed: " + donor + "/" + source_field)
            expect(Path(old_source[source_field]), old_hash)
            expect(Path(source[source_field]), old_hash)
        old_prepared = old_root / "prepared" / donor
        old_receipt = _receipt(old_prepared)
        if (old_receipt["key"] != {"plan": frozen["plan_sha256"], "species": donor}
                or receipt["key"].get("prepared", {}).get(donor) != batch.read([old_prepared / "receipt.json"])[str(old_prepared / "receipt.json")]):
            raise ValueError("Prediction cache prepared producer changed")
        for member, expected in old_receipt["files"].items():
            expect(_safe_member(old_prepared, member), expected)
        expected_protein = old_receipt["files"]["genes.pep"]
        expect(new_root / "prepared" / donor / "genes.pep", expected_protein)
    mapping = {row["id"]: row for row in regions}
    if len(mapping) != len(regions):
        raise ValueError("Duplicate current prediction candidate IDs")
    old_contract = _verified_search_contract(old_root, old_plan, name, receipt, expect)
    expected_contract = prediction_search_contract()
    result = VerifiedPredictionCache(directory, mapping, worker, batch, old_root=old_root,
                                     local_search_compatible=old_contract["local"] == expected_contract["local"],
                                     genome_search_compatible=old_contract["genome"] == expected_contract["genome"],
                                     search_contract=old_contract)
    result.candidate_ids()  # Validate every intersecting empty/aligned query.
    result.bind_exact_genome_sequences()
    result.check()
    return result
