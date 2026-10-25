"""Count qualifying nomination routes for accepted refinement coding paths.

These are predictor origins, not homology donor groups. A Miniprot structural
confirmation attached to native evidence does not constitute a second route.
"""

import copy
import json
from collections import Counter
from pathlib import Path

from input_generation_array_state import FreshDigestBatch

METHODS = ("miniprot_only", "oriongeno_only", "both")
LABELS = ("Miniprot only", "OrionGeno only", "Both methods")
COLOURS = ("#a6cbe5", "#4b91c4", "#154f79")


def classify_path_methods(models, native_identity=None):
    """Partition published accepted candidate IDs, preserving multiple isoforms."""
    paths = {}
    for model in models:
        if model.get("status") != "accepted":
            continue
        candidate = model["candidate"]
        identifier = candidate["candidate_id"]
        if identifier in paths or model["change_type"] not in {"model_revision", "isoform_addition"}:
            raise ValueError("Duplicate accepted coding path or unknown change type")
        methods, donors = set(), set()
        for alignment in model["alignments"]:
            if type(alignment.get("supports_path", True)) is not bool:
                raise ValueError("Ambiguous coding-path support status")
            if not alignment.get("supports_path", True):
                continue
            if alignment.get("annevo") or alignment.get("rescue_revision_model_id"):
                raise ValueError(
                    "Three-route figure requires directly recorded Miniprot/OrionGeno nominations; "
                    "imported rescue revisions need their original nomination provenance"
                )
            donors.add(alignment["donor_species"])
            native = alignment.get("oriongeno")
            if native:
                if (
                    not isinstance(native, dict)
                    or not native.get("native_transcript_id")
                    or not native.get("native_gene_id")
                    or (native_identity is not None and any(native.get(k) != v for k, v in native_identity.items()))
                ):
                    raise ValueError("OrionGeno path origin differs from its frozen native identity")
                methods.add("oriongeno")
            else:
                # The ordinary producer's unlabelled, non-imported alignments
                # originate in its explicit Miniprot interval-search loop.
                methods.add("miniprot")
        if not methods or donors != set(model["donors"]) or donors != set(candidate["support"]["donors"]):
            raise ValueError("Accepted coding path lacks consistent qualifying nomination evidence")
        group = "both" if len(methods) == 2 else next(iter(methods)) + "_only"
        paths[identifier] = dict(gene_id=model["gene_id"], change_type=model["change_type"], category=group)
    counts = Counter(row["category"] for row in paths.values())
    return dict(counts={key: counts[key] for key in METHODS}, paths=paths)


def collect_path_methods(data, root, *, verification=None):
    """Reclassify a summary from freshly verified, immutable prediction members.

    This neither predicts nor adopts a model. The caller can share the fresh digest
    boundary with its other report inputs and check it again before publication.
    """
    batch = verification if verification is not None else FreshDigestBatch()
    root = Path(root).resolve()

    def read(path, expected=None):
        raw = path.read_bytes()
        import hashlib

        actual = hashlib.sha256(raw).hexdigest()
        if expected is not None and actual != expected:
            raise ValueError("Refinement report input changed: " + str(path))
        batch.read([path], expected={str(path): actual})
        return json.loads(raw), actual

    plan, plan_sha = read(root / "plan.json")
    if data.get("refinement_plan_sha256") != plan_sha:
        raise ValueError("Figure counts must declare the same frozen refinement plan SHA")
    names = {row["species"] for row in data["species"] if row["analysed"]}
    if names != set(plan["request"]["sources"]) or any(
        Path(name).name != name or name in {".", ".."} or "\\" in name for name in names
    ):
        raise ValueError("Method summary species differ from the refinement plan")
    derived, evidence = copy.deepcopy(data), {}
    for row in derived["species"]:
        if not row["analysed"]:
            continue
        name = row["species"]
        directory = root / "predictions" / name
        receipt, receipt_sha = read(directory / "receipt.json")
        if receipt["key"]["plan"] != plan_sha:
            raise ValueError("Coding-path prediction receipt belongs to another plan")
        models, _ = read(directory / "predictions.json", receipt["files"]["predictions.json"])
        native = plan["request"].get("oriongeno", {}).get("predictions", {}).get(name, {}).get("identity")
        if not native and any(
            a.get("oriongeno")
            for m in models
            if m.get("status") == "accepted"
            for a in m["alignments"]
            if a.get("supports_path", True)
        ):
            raise ValueError("Native coding-path origin is absent from the frozen plan")
        summary = classify_path_methods(models, native)
        counts = Counter(path["change_type"] for path in summary["paths"].values())
        if (
            counts["model_revision"] != row["coding_paths"]["repairs"]
            or counts["isoform_addition"] != row["coding_paths"]["isoforms"]
        ):
            raise ValueError("Published coding paths differ from the supplied figure totals")
        row["coding_paths"].update(counts=[summary["counts"][key] for key in METHODS], colours=list(COLOURS))
        evidence[name] = dict(prediction_receipt_sha256=receipt_sha, **summary)
    derived["coding_path_legend"] = dict(labels=list(LABELS), colours=list(COLOURS))
    batch.check()
    return derived, dict(
        refinement_plan_sha256=plan_sha,
        species=evidence,
        unit="accepted refinement candidate IDs",
        predictors_executed=False,
    )
