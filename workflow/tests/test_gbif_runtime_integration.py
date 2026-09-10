"""Local GBIF acquisition through the real trait consumers and NWKIT runtime."""

import json
import subprocess
import sys
from pathlib import Path

import pandas

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def run(*args, check=True):
    result = subprocess.run([str(arg) for arg in args], check=False, text=True, capture_output=True)
    if check:
        assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_local_gbif_alias_preserves_meaning_through_copy_number_and_rsc(tmp_path):
    species = [f"Genus_{letter}" for letter in "abcdef"]
    manifest = tmp_path / "download.tsv"
    manifest.write_text("provider\tid\tspecies_key\n" + "".join(f"local\t{i}\t{name}\n" for i, name in enumerate(species)))
    mapping = tmp_path / "taxa.tsv"
    mapping.write_text("species\ttaxon_key\tscientific_name\n" + "".join(f"{name}\t{i + 100}\t{name.replace('_', ' ')}\n" for i, name in enumerate(species)))
    records = tmp_path / "records.tsv"
    latitudes = [12, 19, 17, 32, 28, 45]
    records.write_text("gbifID\tspeciesKey\tdecimalLatitude\tdecimalLongitude\toccurrenceStatus\tissue\n" + "".join(
        f"{i + 1}\t{i + 100}\t{lat}\t140\tPRESENT\t\n" for i, lat in enumerate(latitudes)))
    download = tmp_path / "download.json"
    download.write_text(json.dumps({"status": "SUCCEEDED", "key": "fixture", "doi": "10.example/fixture", "totalRecords": 6,
                                    "request": {"predicate": {"type": "in", "key": "TAXON_KEY", "values": [str(i + 100) for i in range(6)]}}}))
    plan = tmp_path / "traits.plan.tsv"
    metric = "gbif_observed_northern_limit_lat"
    plan.write_text(f"database\tsource_column\toutput_trait\tvalue_type\ngbif\t{metric}\tobserved_north\tnumeric\n")
    traits = tmp_path / "traits.tsv"
    run(sys.executable, SUPPORT / "generate_species_trait.py", "--download-manifest", manifest,
        "--trait-plan", plan, "--database-sources", tmp_path / "unused.tsv", "--databases", "gbif",
        "--gbif-occurrence-file", records, "--gbif-taxon-map", mapping, "--gbif-download-metadata", download,
        "--gbif-require-complete", "yes", "--downloads-dir", tmp_path / "cache", "--output", traits)
    generated = pandas.read_csv(traits, sep="\t")
    assert generated["observed_north"].tolist() == latitudes
    assert all(name == "species" or name == "observed_north" or name.startswith("gbif_observed_") for name in generated)
    source_meta = json.loads(Path(str(traits) + ".metadata.json").read_text())
    assert source_meta["traits"]["observed_north"]["role"] == "observation"

    tree = tmp_path / "species.nwk"
    tree.write_text("(((Genus_a:1,Genus_b:1):1,Genus_c:2):1,((Genus_d:1,Genus_e:1):1,Genus_f:2):1);")
    counts = tmp_path / "counts.tsv"
    counts.write_text("Orthogroup\t" + "\t".join(species) + "\nOG1\t1\t3\t2\t4\t6\t5\n")
    output = tmp_path / "pgls"
    arguments = ["Rscript", SUPPORT / "orthogroup_copy_number_trait_pgls.r", f"--file_orthogroup_copy_number={counts}",
                 f"--file_sptree={tree}", f"--file_trait={traits}", f"--outdir={output}"]
    automatic = run(*arguments, check=False)
    assert automatic.returncode != 0
    assert not output.exists()
    run(*arguments, "--trait=observed_north")
    stats = pandas.read_csv(output / "orthogroup_copy_number_trait_pgls.tsv", sep="\t")
    assert stats["trait"].tolist() == ["observed_north"]
    assert stats["status"].tolist() == ["ok"]
    audit = json.loads((output / "species_trait_input.json").read_text())
    assert audit["traits"]["observed_north"]["source_column"] == metric
    assert "not evidence" in audit["interpretation"]
    assert audit["gbif"]["records_sha256"] == source_meta["gbif"]["records_sha256"]
    assert (output / "selected_species_traits.tsv.metadata.json").exists()
    assert json.loads((output / "selected_species_traits.tsv.schema.json").read_text())["traits"] == {"observed_north": "numeric"}
    assert "GBIF responses describe retained occurrence records" in (output / "orthogroup_copy_number_trait_pgls.summary.svg").read_text()

    expression = tmp_path / "expression.tsv"
    pandas.DataFrame({"gene": [name + "_g1" for name in species], "expression": [2, 5, 3, 8, 9, 7]}).to_csv(expression, sep="\t", index=False)
    prepared = tmp_path / "prepared.tsv"
    metadata = tmp_path / "preparation.tsv"
    run(sys.executable, SUPPORT / "reconciled_speciation_contrast.py", "prepare", "--expression", expression,
        "--species-traits", traits, "--predictors", "observed_north", "--expression-output", tmp_path / "expression.prepared.tsv",
        "--species-traits-output", prepared, "--analysis-plan-output", tmp_path / "analysis.tsv", "--metadata-output", metadata)
    derived = json.loads(Path(str(prepared) + ".metadata.json").read_text())
    assert derived["traits"]["observed_north"]["role"] == "observation"
    assert derived["gbif"]["records_sha256"] == source_meta["gbif"]["records_sha256"]
    preparation = dict(pandas.read_csv(metadata, sep="\t", keep_default_na=False).itertuples(index=False, name=None))
    assert json.loads(preparation["predictor_input_contract"])["selected"] == ["observed_north"]

    reconciliation = tmp_path / "reconciliation.tsv"
    pandas.DataFrame({"node_class": ["tip"] * 6, "gene_name": [name + "_g1" for name in species],
                      "species_name": species}).to_csv(reconciliation, sep="\t", index=False)
    outputs = {field: tmp_path / (field + (".jsonl" if field == "audit" else ".tsv")) for field in (
        "native", "comparison", "status", "audit", "expression-summary", "expression-audit",
        "response-tip-summary", "response-sampling-covariance", "predictor-tip-summary", "predictor-sampling-covariance")}
    output_args = [value for field, path in outputs.items() for value in ("--" + field + "-out", path)]
    run(sys.executable, SUPPORT / "species_tree_pgls.py", "--methods", "species-nwkit", "--tree-id", "OG1",
        "--species-tree", tree, "--reconciliation", reconciliation, "--expression", tmp_path / "expression.prepared.tsv",
        "--species-traits", prepared, "--analysis-plan", tmp_path / "analysis.tsv", "--metadata", metadata,
        "--expression-value-type", "identity", "--response-evolution-model", "brownian",
        "--predictor-evolution-model", "brownian", *output_args)
    method_audit = [json.loads(line) for line in outputs["audit"].read_text().splitlines()]
    input_audit = next(record for record in method_audit if record.get("stage") == "species_trait_input")
    assert input_audit["traits"]["observed_north"]["role"] == "observation"
    assert input_audit["gbif"]["records_sha256"] == source_meta["gbif"]["records_sha256"]
    method_status = pandas.read_csv(outputs["status"], sep="\t")
    assert method_status.loc[method_status["analysis_method"].eq("species_nwkit"), "status"].tolist() == ["ok"]
    # The default RSC-only path skips species PGLS, but must retain the audit
    # after the core removes its temporary prepared table and metadata.
    run(sys.executable, SUPPORT / "species_tree_pgls.py", "--methods", "rsc", "--tree-id", "OG1",
        "--species-tree", tree, "--reconciliation", reconciliation, "--expression", tmp_path / "expression.prepared.tsv",
        "--species-traits", prepared, "--analysis-plan", tmp_path / "analysis.tsv", "--metadata", metadata,
        "--response-evolution-model", "brownian", "--predictor-evolution-model", "brownian", *output_args)
    rsc_audit = [json.loads(line) for line in outputs["audit"].read_text().splitlines()]
    assert next(record for record in rsc_audit if record.get("stage") == "species_trait_input")["gbif"]["records_sha256"] == source_meta["gbif"]["records_sha256"]




def test_copy_number_selection_explicit_observation_masks_partial_and_carries_audit(tmp_path):
    import numpy as np
    sys.path.insert(0, str(SUPPORT))
    from species_trait_contract import write_trait_bundle
    rng = np.random.default_rng(11)
    leaves = [f"sp{i}" for i in range(18)]
    tree, counts, folds, traits = [tmp_path / name for name in ("tree.nwk", "counts.tsv", "folds.tsv", "traits.tsv")]
    tree.write_text("[&R] (" + ",".join(f"{leaf}:1" for leaf in leaves) + ");")
    copy = rng.poisson(2, (2, 18))
    pandas.DataFrame({"Orthogroup": ["OG1", "OG2"], **dict(zip(leaves, copy.T, strict=True))}).to_csv(counts, sep="\t", index=False)
    pandas.DataFrame({"leaf_name": leaves, "fold": np.repeat(["a", "b", "c"], 6)}).to_csv(folds, sep="\t", index=False)
    definition = {"observed_north": {"role": "observation", "source": "gbif", "source_column": "gbif_observed_northern_limit_lat"}}
    quality = [{"species": leaf, "status": "complete_search", "analysis_eligible": True} for leaf in leaves]
    quality[0].update(status="capped_partial", analysis_eligible=False, termination_reason="record_limit")
    write_trait_bundle(pandas.DataFrame({"species": leaves, "observed_north": 30 + copy[0] + rng.normal(size=18)}),
                       traits, definition, {"quality": quality, "observations": []}, [])
    out = tmp_path / "selection"
    args = [sys.executable, SUPPORT / "orthogroup_copy_number_trait_selection.py", "--copy-number", counts,
            "--tree", tree, "--traits", traits, "--folds", folds, "--outdir", out, "--strengths", "0.5", "--l1-ratios", "0.5"]
    automatic = run(*args, check=False)
    assert automatic.returncode != 0 and "automatically eligible" in automatic.stderr
    assert not out.exists()
    run(*args, "--trait", "observed_north")
    manifest = pandas.read_csv(out / "manifest.tsv", sep="\t")
    assert manifest["n_species"].tolist() == [17]
    audit = json.loads((out / "species_trait_input.json").read_text())
    assert audit["masked_species"]["observed_north"][0]["species"] == "sp0"
    selected = pandas.read_csv(out / "selected_species_traits.tsv", sep="\t")
    assert pandas.isna(selected.loc[selected.species.eq("sp0"), "observed_north"]).all()
    metadata = json.loads((out / "trait_0001.metadata.json").read_text())
    assert metadata["trait_input_audit"]["traits"] == definition
    assert "not evidence" in metadata["trait_input_audit"]["interpretation"]
