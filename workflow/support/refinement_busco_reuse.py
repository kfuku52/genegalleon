"""Reuse completed BUSCO leaves only under a verified computation identity."""
from __future__ import annotations

import ast
import copy
import hashlib
import json
import shutil
from pathlib import Path

from input_generation_array_state import FreshDigestBatch, atomic_json, digest, safe_component

KERNELS = ("phases", "protein_genetic_code", "read_genetic_codes", "translate_cds",
           "_validate_admitted_input", "protein_input_identity", "verify_translation",
           "contracted_result", "read_result", "paired_result", "staged_result", "run_one")
DEPENDENCIES = ("busco_quality_metadata.py", "fasta_sequence_store.py", "refinement_busco_inputs.py")


def _kernel_identity(raw, *, current=False):
    """Inspect producer source without executing historical Python code."""
    nodes = {}
    for node in ast.parse(raw).body:
        if isinstance(node, ast.FunctionDef):
            name = "run_one" if current and node.name == "_run_one" else node.name
            if name in KERNELS and not (current and node.name == "run_one"):
                value = copy.deepcopy(node)
                value.name = name
                nodes[name] = ast.dump(value, include_attributes=False)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PROTEIN_TRANSLATION_POLICY"
                                                 for t in node.targets):
            nodes["PROTEIN_TRANSLATION_POLICY"] = ast.dump(node, include_attributes=False)
    if set(nodes) != {*KERNELS, "PROTEIN_TRANSLATION_POLICY"}:
        raise ValueError("Previous producer lacks the complete protein BUSCO computation kernel")
    return hashlib.sha256(json.dumps(nodes, sort_keys=True).encode()).hexdigest()


def prepare_reuse(previous, producer, contract, implementation):
    """Bind all parameters/tools plus an equivalent historical computation kernel.

    Only report implementation identity and physical lineage/download locators
    can differ. The lineage contents, translated protein bytes, all numerical
    settings and tool/wrapper identities remain mandatory.
    """
    if previous is None or producer is None:
        raise ValueError("Previous report and exact producer source are both required")
    if contract.get("mode") != "proteins" or contract.get("translation", {}).get("input_kind") != "admitted":
        raise ValueError("Previous-report reuse requires admitted protein inputs")
    previous, producer, implementation = Path(previous).resolve(), Path(producer).resolve(), Path(implementation).resolve()
    boundary = FreshDigestBatch()
    manifest = previous / "contract.json"
    paths = [manifest, producer, implementation]
    paths.extend(root.parent / name for root in (producer, implementation) for name in DEPENDENCIES)
    hashes = boundary.read(paths)
    raw = producer.read_bytes()
    current_raw = implementation.read_bytes()
    frozen = json.loads(manifest.read_text())["contract"]
    if (hashlib.sha256(raw).hexdigest() != hashes[str(producer)]
            or frozen.get("implementation") != hashes[str(producer)]
            or contract.get("implementation") != hashes[str(implementation)]
            or hashlib.sha256(current_raw).hexdigest() != hashes[str(implementation)]):
        raise ValueError("Previous BUSCO producer does not match its recorded implementation")
    previous_kernel = _kernel_identity(raw)
    current_kernel = _kernel_identity(current_raw, current=True)
    if previous_kernel != current_kernel:
        raise ValueError("Previous BUSCO computation kernel differs from the current implementation")
    for name in DEPENDENCIES:
        if hashes[str(producer.parent / name)] != hashes[str(implementation.parent / name)]:
            raise ValueError("Previous BUSCO dependency differs: " + name)
    locators = {"implementation", "lineage_path", "download_path"}
    if {k:v for k,v in frozen.items() if k not in locators} != {k:v for k,v in contract.items() if k not in locators}:
        raise ValueError("Previous BUSCO computation contract differs")
    boundary.check()
    return {"report": previous, "contract": frozen, "boundary": boundary,
            "proof": {"producer_source": str(producer), "producer_sha256": hashes[str(producer)],
                      "contract_sha256": hashes[str(manifest)], "computation_kernel_sha256": current_kernel,
                      "dependency_sha256": {name: hashes[str(implementation.parent / name)] for name in DEPENDENCIES}}}


def reuse_previous(pair, phase, report, contract, reuse, helpers):
    """Copy identical verified proteins/scores; bind a fresh current-input receipt."""
    import fcntl

    previous, report = reuse["report"], Path(report).resolve()
    if report == previous or report.is_relative_to(previous) or previous.is_relative_to(report):
        raise ValueError("Previous report must be separate from the current report")
    if not safe_component(pair["species"]) or phase not in helpers["phases"](pair):
        raise ValueError("Invalid BUSCO reuse species/phase")
    source = previous / "runs" / pair["species"] / phase
    receipt_path = source / "receipt.json"
    if not receipt_path.is_file():
        return False
    reuse["boundary"].check()
    source_hash = digest(Path(pair[phase]))
    key = helpers["protein_input_identity"](pair, phase, contract, source_hash)
    descriptor = key["protein_input"]
    code = contract["translation"]["genetic_codes"][pair["species"]]
    helpers["_validate_admitted_input"](descriptor, source_hash, code, species=pair["species"], phase=phase)
    boundary = FreshDigestBatch()
    members = {"summary.txt":"summary_sha256", "full_table.tsv":"full_table_sha256",
               "protein.faa":"protein_sha256", "translation_audit.json":"translation_audit_sha256"}
    hashes = boundary.read([receipt_path, *[source / name for name in members]])
    saved = json.loads(receipt_path.read_text())
    if saved.get("key", {}).get("contract") != reuse["contract"]:
        raise ValueError("Previous BUSCO leaf belongs to another computation contract")
    for name, field in members.items():
        if hashes[str(source / name)] != saved.get(field):
            raise ValueError("Previous BUSCO leaf changed: " + name)
    helpers["verify_translation"](source, saved)
    helpers["contracted_result"](source / "summary.txt", contract)
    if saved["protein_sha256"] != descriptor["protein_sha256"]:
        boundary.check()
        return False
    directory = report / "runs" / pair["species"] / phase
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (directory / "receipt.json").exists():
            return False
        for name in ("summary.txt", "full_table.tsv"):
            shutil.copyfile(source / name, directory / name)
        shutil.copyfile(descriptor["protein_path"], directory / "protein.faa")
        shutil.copyfile(descriptor["audit_path"], directory / "translation_audit.json")
        expected = {"summary.txt":saved["summary_sha256"], "full_table.tsv":saved["full_table_sha256"],
                    "protein.faa":descriptor["protein_sha256"], "translation_audit.json":descriptor["audit_sha256"]}
        for name, sha in expected.items():
            if digest(directory / name) != sha:
                raise OSError("BUSCO reuse input changed during copying: " + name)
        helpers["_validate_admitted_input"](descriptor, source_hash, code, species=pair["species"], phase=phase)
        if digest(Path(pair[phase])) != source_hash:
            raise OSError("Current BUSCO DNA input changed during reuse")
        boundary.check()
        reuse["boundary"].check()
        proof = dict(reuse["proof"], source_receipt=str(receipt_path), source_receipt_sha256=hashes[str(receipt_path)],
                     basis="Identical admitted protein bytes and verified full computation contract/kernel")
        atomic_json(directory / "receipt.json", {"key":key,
                    **{members[name]:sha for name,sha in expected.items()}, "reused_verified_previous":proof})
    return True
