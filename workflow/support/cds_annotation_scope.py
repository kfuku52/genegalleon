"""Keep the original gene catalogue separate from CDS annotation eligibility."""
import gzip
import hashlib
import json
from pathlib import Path

SCHEMA = "annotation_cds_scope_v1"
QC_COLUMNS = ("cds_qc_status", "cds_qc_decision", "cds_qc_exclusion_reasons",
              "cds_qc_candidates_json", "cds_qc_selected_source", "functional_annotation_status")


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ids_digest(ids):
    return hashlib.sha256(canonical(sorted(ids)).encode()).hexdigest()


def fasta_ids(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        ids = [line[1:].split()[0] for line in handle if line.startswith(">")]
    if len(set(ids)) != len(ids) or any(not item for item in ids):
        raise ValueError("Original/resolved CDS IDs must be nonempty and unique")
    return ids


def load_scope(catalogue, resolved, report_path):
    """Validate existing resolution evidence; do not evaluate or relax CDS QC."""
    paths = [Path(value) for value in (catalogue, resolved, report_path)]
    guards = [(p.stat().st_dev, p.stat().st_ino, p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in paths]
    original_ids, accepted_ids = fasta_ids(catalogue), fasta_ids(resolved)
    report = json.loads(Path(report_path).read_text())
    contract = report["contract"]
    catalogue_sha = sha256(catalogue)
    if (contract["sources"]["cds"]["sha256"] != catalogue_sha
            or report["output_sha256"] != sha256(resolved)):
        raise ValueError("CDS catalogue/resolved sequence SHA differs from resolution evidence")
    decisions = report["decisions"]
    by_id = {item["gene_id"]: item for item in decisions}
    accepted = {key for key, item in by_id.items() if item["selected_source"] is not None}
    if (len(by_id) != len(decisions) or set(by_id) != set(original_ids)
            or accepted != set(accepted_ids) or report["input_count"] != len(original_ids)
            or report["accepted_count"] != len(accepted)
            or report["excluded_count"] != len(original_ids) - len(accepted)):
        raise ValueError("CDS resolution catalogue/eligibility coverage mismatch")
    rows = {}
    for gene, item in by_id.items():
        selected = item["selected_source"]
        candidates = item["candidates"]
        if selected is not None and (selected not in candidates or not candidates[selected]["accepted"]):
            raise ValueError("Selected CDS candidate was not accepted by resolution QC")
        if selected is None and (not candidates or any(v["accepted"] for v in candidates.values())):
            raise ValueError("Excluded CDS requires rejected candidate evidence")
        reasons = {key: value["reason"] for key, value in candidates.items()} if selected is None else {}
        rows[gene] = dict(zip(QC_COLUMNS, (
            "accepted" if selected is not None else "excluded", item["reason"], canonical(reasons),
            canonical(candidates), selected or "",
            "eligible" if selected is not None else "not_assessed_cds_qc_excluded"), strict=True))
    scope = {"schema": SCHEMA, "catalogue_sha256": catalogue_sha,
             "resolved_cds_sha256": report["output_sha256"], "resolution_sha256": sha256(report_path),
             "resolution_contract": contract, "catalogue_gene_ids_sha256": ids_digest(original_ids),
             "accepted_gene_ids_sha256": ids_digest(accepted_ids), "input_count": len(original_ids),
             "accepted_count": len(accepted), "excluded_count": len(original_ids) - len(accepted)}
    if guards != [(p.stat().st_dev, p.stat().st_ino, p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in paths]:
        raise ValueError("CDS scope input changed while reading")
    return original_ids, accepted, rows, scope


def validate_qc_row(row, product, go_columns):
    """Absent legacy QC is unknown, never implicitly accepted."""
    if not set(QC_COLUMNS) & row.keys():
        return "legacy_unspecified"
    if not set(QC_COLUMNS) <= row.keys():
        raise ValueError("Incomplete CDS QC annotation columns")
    status = row["cds_qc_status"]
    candidates = json.loads(row["cds_qc_candidates_json"])
    reasons = json.loads(row["cds_qc_exclusion_reasons"])
    selected = row["cds_qc_selected_source"]
    if status == "accepted":
        if (row["functional_annotation_status"] != "eligible" or reasons or selected not in candidates
                or candidates[selected].get("accepted") is not True):
            raise ValueError("Contradictory accepted CDS QC annotation")
    elif status == "excluded":
        if (row["functional_annotation_status"] != "not_assessed_cds_qc_excluded" or selected
                or not candidates or any(v.get("accepted") is not False for v in candidates.values())
                or reasons != {key: value["reason"] for key, value in candidates.items()}):
            raise ValueError("Excluded CDS requires explicit rejection reasons and unassessed status")
        predicted = {product, *go_columns.values(), *(key for key in row if key.startswith(("sprot_", "cdskit_localize_")))}
        if any(row.get(key) for key in predicted):
            raise ValueError("QC-excluded CDS cannot carry functional predictions")
    else:
        raise ValueError("Unknown CDS QC status")
    return status
