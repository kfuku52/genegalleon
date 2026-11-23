"""Admit source-verified OrionGeno models at uncovered coding loci."""

import json
import shutil
from collections import defaultdict
from pathlib import Path

from annevo_gene_models import read_paths, verify_inputs

BASIS = "oriongeno_without_donor"
AUDIT = BASIS + ".json"


def contract():
    from rescue_gene_models import digest

    return {"schema": 1, "implementation_sha256": digest(Path(__file__))}


def admit(frozen, species, genome, code, params, original, models):
    """Fill uncovered loci after ordinary rescue and audit every native gene."""
    from rescue_gene_models import OwnershipIndex, consolidate, require_same_key, validate_model

    prediction = frozen["predictions"].get(species)
    if prediction is None:
        return {"native_prediction_available": False, "genes": {}}

    verify_inputs(frozen, [species])
    paths_by_gene = defaultdict(list)
    for path in read_paths(prediction["gff"], genome, engine="oriongeno"):
        paths_by_gene[path["native_gene_id"]].append(path)

    occupied, represented = [], set()
    for model in models:
        if model["status"] != "accepted":
            continue
        for path in (model, *model.get("alternative_coding_paths", [])):
            occupied.append({"seqid": path["seqid"], "strand": path["strand"], "cds": path["cds"],
                             "gene_id": model["model_id"],
                             "start": min(a for a, _, _ in path["cds"]),
                             "end": max(b for _, b, _ in path["cds"])})
            represented.update(row["oriongeno"]["native_gene_id"] for row in path["support"]
                               if row.get("oriongeno"))
    rescued = OwnershipIndex(occupied)
    qc_params = {**params, "minimum_coverage": 0, "minimum_identity": 0, "minimum_broad_coverage": 0}
    attempted, genes = [], {}
    for gene, paths in paths_by_gene.items():
        audit = genes[gene] = {}
        if any(original.overlapping(path) for path in paths):
            audit["status"] = "original_cds_overlap"
        elif gene in represented or any(rescued.overlapping(path) for path in paths):
            audit["status"] = "already_rescued"
        else:
            valid, problems = [], set()
            for path in paths:
                candidate = validate_model(
                    {**path, "coverage": 0, "identity": 0, "frameshift": False}, genome, code, qc_params)
                if candidate["problems"]:
                    problems.update(candidate["problems"])
                else:
                    candidate.update(start=min(a for a, _, _ in candidate["cds"]),
                                     end=max(b for _, b, _ in candidate["cds"]))
                    valid.append(candidate)
            if not valid:
                audit.update(status="invalid_cds", reasons=sorted(problems))
            elif any(original.overlapping(path) for path in valid):
                audit["status"] = "original_cds_overlap"
            elif any(rescued.overlapping(path) for path in valid):
                audit["status"] = "already_rescued"
            else:
                path = min(valid, key=lambda p: (-len(p["sequence"]), p["native_transcript_id"]))
                attempted.append(path)
                audit["native_transcript_id"] = path["native_transcript_id"]

    # A conflicting same-strand CDS from another native locus is ambiguous.
    contenders = OwnershipIndex([
        {key: path[key] for key in ("seqid", "strand", "cds", "start", "end")} | {"gene_id": path["native_gene_id"]}
        for path in attempted
    ])
    additions = []
    for path in attempted:
        gene, transcript = path["native_gene_id"], path["native_transcript_id"]
        if contenders.overlapping(path, owner=gene):
            genes[gene]["status"] = "conflicting_native_locus"
            continue
        path.update(coverage=None, identity=None, query=f"{BASIS}:{transcript}", search=BASIS,
                    evidence={"target": species, "query": "", "admission_basis": BASIS,
                              "oriongeno": {**prediction["identity"], "native_gene_id": gene,
                                           "native_transcript_id": transcript}},
                    placement_evidence={"classification": BASIS, "independent_donor_species": [],
                                        "orthology": "unassigned", "expected_copy": "unassigned"})
        additions.append(path)
    # Use the same genomic IDs and accepted-model format as ordinary rescue.
    consolidate(additions, [], species)
    for model in additions:
        require_same_key("accepted", model["status"])
        genes[model["native_gene_id"]].update(status="admitted", model_id=model["model_id"])
    verify_inputs(frozen, [species])
    models.extend(additions)
    return {"native_prediction_available": True, "genes": genes}


def extend_worker(root, plan, species, cpus, dependencies):
    """Append validated native paths to a frozen rescue without replaying donors."""
    import rescue_admission_replay as replay
    import rescue_diagnostics as diagnostics
    import rescue_gene_models as rescue
    import rescue_model_store as store
    from gene_model_catalog import indexed_genome
    from rescue_model_store_patch import patch_model_store

    previous, prior = replay.previous_context(root, plan)
    cached = plan["request"]["admission_cache"]
    producer = previous / "rescued" / species
    expected_receipt = cached["species"][species]["receipt_sha256"]
    receipt, receipt_sha = store._receipt(producer)
    rescue.require_same_key(expected_receipt, receipt_sha)
    rescue.require_same_key(receipt["key"], rescue.rescue_key(previous, prior, species, dependency_snapshot=dependencies))
    key = rescue.rescue_key(root, plan, species, dependency_snapshot=dependencies)
    source = plan["request"]["sources"][species]
    source_proof = rescue.RescueSourceBatch(plan)
    proof = rescue.FreshDigestBatch()
    expected = {str(producer / "receipt.json"): receipt_sha,
                str(previous / "plan.json"): cached["plan_sha256"],
                str(root / "plan.json"): key["plan"]}
    rescue.require_same_key(expected, proof.read(expected))

    def guard():
        proof.check()
        dependencies.check_all()
        source_proof.check()

    def build(tmp):
        rescue.verify_sources(plan, [species], ["genome", "fasta", "gff"], verification=source_proof)
        verify_inputs(plan["request"]["oriongeno"], [species], verification=source_proof)
        frozen = store.frozen_model_store_key(producer)
        accepted = list(store.iter_accepted_models(producer, verification=proof))
        initial = len(accepted)
        with rescue.open_text(Path(source["gff"])) as stream:
            original = rescue.OwnershipIndex(rescue.original_ownership(stream, rescue.attributes))
        with indexed_genome(source["genome"]) as genome:
            audit = admit(plan["request"]["oriongeno"], species, genome, source["genetic_code"],
                          rescue.parameters_for(plan["request"], species), original, accepted)
        additions = accepted[initial:]
        rescue.atomic_json(tmp / AUDIT, audit)
        if additions:
            for model in additions:
                model["quality_evidence"] = rescue.model_quality(model, source["genetic_code"])
            patch_model_store(producer, tmp, {}, appended=additions, frozen_key=frozen,
                              revisions=store.iter_revision_models(producer), publication_workers=min(cpus, 4))

        mutable = {"search_inputs.json", diagnostics.SUMMARY, "audit.tsv", "audit.tsv.gz",
                   "quality_flags.tsv", "quality_flags.tsv.gz"}

        def copy_member(src):
            relative = src.relative_to(producer).as_posix()
            target = tmp / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, target)
            rescue.require_same_key(receipt["files"][relative], rescue.digest(target))

        rescue.hash_paths(producer, (store._member(producer, name) for name in receipt["files"]
                                    if name not in mutable and not (additions and name.startswith("model_store/"))),
                          workers=min(cpus, 4), hash_function=copy_member)

        def read_sidecar(name):
            path = producer / name
            rescue.require_same_key(receipt["files"][name], proof.read([path])[str(path)])
            return json.loads(path.read_text())

        if diagnostics.SUMMARY in receipt["files"]:
            summary = read_sidecar(diagnostics.SUMMARY)
            updated = diagnostics.replace_counts(summary, diagnostics.new_counts(), diagnostics.summarize(additions))
            rescue.atomic_json(tmp / diagnostics.SUMMARY, updated)
        else:
            n = store._manifest(producer / "model_store")[0]["counts"]["models"]
            for name, row, allow_tail in (("audit.tsv", diagnostics.audit_row, True),
                                          ("quality_flags.tsv", diagnostics.quality_row, False)):
                path = diagnostics.legacy_table(producer, name)
                if path is not None:
                    rescue.require_same_key(receipt["files"][path.name], proof.read([path])[str(path)])
                    replay._rewrite_table(path, tmp / path.name, {}, additions, row, n, allow_tail=allow_tail)
        if "search_inputs.json" in receipt["files"]:
            search = read_sidecar("search_inputs.json")
            search.update(plan_sha256=key["plan"], prepared_receipts=key["prepared"],
                          prediction_cache=plan["request"]["prediction_cache"])
            rescue.atomic_json(tmp / "search_inputs.json", search)
        rescue.atomic_json(tmp / "oriongeno_admission_parent.json",
                           {"plan_sha256": cached["plan_sha256"], "receipt_sha256": receipt_sha,
                            "source_worker": str(producer), "appended_models": len(additions)})

    return rescue.stage(root, Path("rescued") / species, key, build, guard,
                        publication_guard=guard, hash_workers=min(cpus, 4))


def extend(previous, output, cpus):
    """Publish an augmented rescue using the same frozen inputs and donors."""
    import rescue_admission_replay as replay
    import rescue_gene_models as rescue

    root, previous = output.resolve(), previous.resolve(strict=True)
    if (root / "plan.json").is_file():
        plan = rescue.load(root)
        cache = plan["request"].get("admission_cache", {})
        rescue.require_same_key(str(previous), cache.get("root"))
        rescue.require_same_key("oriongeno_additions", cache.get("contract", {}).get("scope"))
    else:
        plan = replay.build_plan(previous, None, root, oriongeno_additions=True)
    replay.prepare(root, plan, cpus)
    for species in plan["species"]:
        rescue.rescue(root, plan, species, cpus)
    return rescue.finalize(root, plan, hash_workers=cpus)
