#!/usr/bin/env python3
"""Compare observation summaries from one verified saved GBIF record snapshot.

This measures sensitivity to record selection, not uncertainty in true ranges.
No network access or species-distribution model is used.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
from pathlib import Path

import pandas
from gbif_observations import (
    METRIC_DEFINITIONS,
    SCHEMA_VERSION,
    build_gbif_distribution_metrics,
    effective_config,
    file_sha256,
    gbif_row_to_point,
    occupied_cell,
    read_snapshot,
    summarize_records,
)


def choices(value: str) -> list[str]:
    return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))


def run(args: argparse.Namespace) -> int:
    metadata = json.loads(args.metadata.read_text())
    bundle = metadata.get("gbif", metadata)
    if bundle.get("identity", {}).get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported GBIF snapshot schema")
    snapshot = Path(bundle["records_path"])
    if file_sha256(snapshot) != bundle["records_sha256"]:
        raise ValueError("GBIF snapshot hash differs from its acquisition metadata")
    config = effective_config(bundle["identity"]["configuration"])
    grids = [config["gbif_grid_degrees"]] if args.grid_degrees == "source" else [float(value) for value in choices(args.grid_degrees)]
    windows, sizes, selections = choices(args.year_windows), choices(args.sample_sizes), choices(args.selection)
    establishments = choices(args.establishment_means)
    countries = list(dict.fromkeys(value.strip() for value in args.countries.split(";") if value.strip()))
    missing_dates = choices(args.missing_date)
    exclusions = [None, *choices(args.exclude_datasets)]
    if not all([grids, windows, sizes, selections, establishments, countries, missing_dates]) or args.replicates < 1:
        raise ValueError("Sensitivity dimensions and replicates must be nonempty/positive")
    if set(selections) - {"records", "unique_location", "one_per_cell"}:
        raise ValueError("selection must be records, unique_location or one_per_cell")
    if set(missing_dates) - {"source", "keep", "exclude"}:
        raise ValueError("missing-date must be source, keep or exclude")
    if any(size != "all" and (not size.isdigit() or int(size) < 1) for size in sizes):
        raise ValueError("sample-sizes must contain all or positive integers")
    acquisitions = {}
    for entry in read_snapshot(snapshot):
        species = entry["species"]
        if "acquisition" in entry:
            acquisitions[species] = {**entry["acquisition"], "rows": []}
        else:
            acquisitions[species]["rows"].append(entry["record"])
    dimensions = list(itertools.product(grids, windows, establishments, countries, missing_dates, exclusions, selections, sizes))
    results = []
    for species, acquisition in sorted(acquisitions.items()):
        for grid, window, establishment, country, missing_date, excluded_dataset, selection, size in dimensions:
            scenario_config = {**config, "gbif_grid_degrees": grid}
            if window != "source":
                if window == "all":
                    scenario_config.update(gbif_year_min=None, gbif_year_max=None)
                else:
                    start, end = window.split(":")
                    scenario_config.update(gbif_year_min=int(start), gbif_year_max=int(end))
            if establishment != "source":
                scenario_config["gbif_include_establishment_means"] = "" if establishment == "all" else establishment
            if country != "source":
                scenario_config["gbif_countries"] = "" if country == "all" else country
            if missing_date != "source":
                scenario_config["gbif_missing_date"] = missing_date
            scenario_config = effective_config(scenario_config)
            scenario_acquisition = dict(acquisition)
            scenario_acquisition["rows"] = [row for row in acquisition["rows"] if not excluded_dataset or row.get("datasetKey") != excluded_dataset]
            _metrics, quality, retained = summarize_records(scenario_acquisition, scenario_config)
            retained.sort(key=lambda row: (str(row.get("gbifID") or ""), json.dumps(row, sort_keys=True)))
            replicates = args.replicates if size != "all" or selection != "records" else 1
            scenario = json.dumps([grid, window, establishment, country, missing_date, excluded_dataset, selection, size], separators=(",", ":"))
            for replicate in range(replicates):
                seed = int.from_bytes(hashlib.sha256(f"{args.seed}|{species}|{scenario}|{replicate}".encode()).digest()[:8], "big")
                rng = random.Random(seed)
                selected = retained.copy()
                if selection != "records":
                    groups = {}
                    for row in retained:
                        point = gbif_row_to_point(row, scenario_config)
                        key = point[:2] if selection == "unique_location" else occupied_cell(point, grid)
                        groups.setdefault(key, []).append(row)
                    selected = [rng.choice(rows) for _key, rows in sorted(groups.items())]
                available_units = len(selected)
                insufficient = size != "all" and len(selected) < int(size)
                if insufficient:
                    selected = []
                elif size != "all":
                    selected = rng.sample(selected, int(size))
                points = [gbif_row_to_point(row, scenario_config) for row in selected]
                results.append({
                    "species": species, "grid_degrees": grid, "year_window": window,
                    "countries": country, "missing_date": missing_date,
                    "establishment_means": establishment, "excluded_dataset": excluded_dataset or "",
                    "selection": selection, "sample_size_requested": size, "replicate": replicate,
                    "seed": seed, "acquisition_status": acquisition["status"],
                    "status": "insufficient_records" if insufficient else ("no_retained_records" if not selected else "descriptive"),
                    "retained_before_subsampling": len(retained), "selected_records": len(selected),
                    "available_sampling_units": available_units,
                    "selected_ids_sha256": hashlib.sha256(json.dumps(sorted(str(row.get("gbifID") or "") for row in selected), separators=(",", ":")).encode()).hexdigest(),
                    "unknown_date_count_before_subsampling": quality["unknown_date_count"],
                    **build_gbif_distribution_metrics(points, grid),
                })
    if not results:
        raise ValueError("GBIF snapshot contains no species")
    from species_trait_contract import trait_output_transaction
    sidecar = Path(str(args.output) + ".metadata.json")
    outputs = [args.output, sidecar]
    with trait_output_transaction(outputs, [snapshot, args.metadata]) as staged:
        pandas.DataFrame.from_records(results).to_csv(staged[args.output], sep="\t", index=False, na_rep="NA")
        report = {"schema_version": 1, "source_metadata_sha256": file_sha256(args.metadata),
                  "records_path": str(snapshot), "records_sha256": bundle["records_sha256"],
                  "source_identity": bundle["identity"], "table_sha256": file_sha256(staged[args.output]),
                  "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
                  "metric_definitions": METRIC_DEFINITIONS,
                  "interpretation": "Descriptive sensitivity to selection within this acquisition universe, not bias correction or confidence intervals for true ranges",
                  "sampling": "Without replacement after GBIF-ID deduplication and filtering; each location/cell first supplies one uniformly chosen record. Unit inclusion probability is selected_records/available_sampling_units; within a selected group each record has probability 1/group_size. Species with fewer units than requested are not silently given a smaller sample"}
        Path(staged[sidecar]).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {len(results)} descriptive sensitivity rows to {args.output}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path, help="Generated species_trait.tsv.metadata.json or saved acquisition cache JSON")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--grid-degrees", default="source", help="source or comma-separated degree grid sizes")
    parser.add_argument("--year-windows", default="source", help="Comma-separated source, all, or MIN:MAX years; all remains limited to the original acquisition universe")
    parser.add_argument("--establishment-means", default="source", help="Comma-separated source, all, native, introduced, etc.; missing is never assumed native")
    parser.add_argument("--countries", default="source", help="Semicolon-separated country-code sets, source or all; e.g. 'source;JP;US,CA'")
    parser.add_argument("--missing-date", default="source", help="Comma-separated source, keep, exclude; applies with active year windows")
    parser.add_argument("--exclude-datasets", default="", help="Compare baseline with excluding each listed dataset key in turn")
    parser.add_argument("--selection", default="records", help="Comma-separated records, unique_location, one_per_cell")
    parser.add_argument("--sample-sizes", default="all", help="Comma-separated all or positive integers")
    parser.add_argument("--replicates", default=20, type=int)
    parser.add_argument("--seed", default=1, type=int)
    return run(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
