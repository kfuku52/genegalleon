#!/usr/bin/env python3
"""NWKIT OU adapter: replicate means, sampling errors and artifacts."""

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def aggregate_replicates(table, separator="_"):
    if table.empty or table.shape[1] < 2:
        raise ValueError("Expression input needs named genes and at least one trait.")
    names = table.iloc[:, 0]
    if names.isna().any() or names.astype(str).duplicated().any():
        raise ValueError("Expression gene names must be nonmissing and unique.")
    groups = {}
    for column in table.columns[1:]:
        trait = column.rsplit(separator, 1)[0] if separator and separator in column else column
        if not trait or "," in trait:
            raise ValueError("Trait names must be nonempty and cannot contain commas.")
        groups.setdefault(trait, []).append(column)
    aggregated = pd.DataFrame({"leaf_name": names.astype(str)})
    audit = []
    for index, (trait, columns) in enumerate(groups.items()):
        values = table[columns].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
        if np.isinf(values).any():
            raise ValueError("Replicate values cannot contain infinity.")
        count = np.isfinite(values).sum(axis=1)
        mean = np.divide(np.nansum(values, axis=1), count, out=np.full(len(values), np.nan), where=count > 0)
        ss = np.nansum((values - mean[:, None]) ** 2, axis=1)
        variance = np.divide(ss, count * (count - 1), out=np.zeros(len(values)), where=count > 1)
        error_column = f"__native_se_{index}"
        if error_column in groups or trait == "leaf_name":
            raise ValueError("Trait name collides with an adapter output column.")
        aggregated[trait] = mean
        aggregated[error_column] = np.sqrt(variance)
        audit.extend(
            {
                "leaf_name": str(name),
                "trait": trait,
                "replicate_count": int(n),
                "mean": float(value) if n else None,
                "sampling_variance_of_mean": float(error) if n else None,
                "sampling_variance_status": "replicate_estimate"
                if n > 1
                else ("unavailable_single_observation" if n else "missing"),
            }
            for name, n, value, error in zip(names, count, mean, variance, strict=True)
        )
    return aggregated, list(groups), [f"__native_se_{i}" for i in range(len(groups))], audit


def branch_summary(model):
    if model.get("selection") != "native" or model.get("completion_status") != "complete":
        raise ValueError("Native branch summary requires a completed NWKIT model.")
    support = (model.get("selection_support") or {}).get("branch_frequencies", [])
    frequencies = {row["branch_id"]: row["frequency"] for row in support}
    shifts = set(model["shift_branch_ids"])
    rows = [
        {
            "branch_id": row["branch_id"],
            "ou_native_regime": row["regime"],
            "ou_native_is_shift": int(row["branch_id"] in shifts),
            "ou_native_selection_frequency": frequencies.get(row["branch_id"]),
            "ou_native_research_only": True,
        }
        for row in model["branches"]
    ]
    return pd.DataFrame(rows)


def plot_native_model(model, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D

    branches = {row["branch_id"]: row for row in model["branches"]}
    children = {branch: [] for branch in branches}
    for branch, row in branches.items():
        if row["parent"] != -1:
            children[row["parent"]].append(branch)
    order, stack = [], [0]
    while stack:
        branch = stack.pop()
        order.append(branch)
        stack.extend(reversed(children[branch]))
    tips = [branch for branch in order if not children[branch]]
    y = {branch: len(tips) - i for i, branch in enumerate(tips)}
    for branch in reversed(order):
        if children[branch]:
            y[branch] = np.mean([y[child] for child in children[branch]])
    x = {0: 0.0}
    for branch in order[1:]:
        x[branch] = x[branches[branch]["parent"]] + branches[branch]["dist"]
    regimes = sorted({row["regime"] for row in branches.values()})
    palette = plt.get_cmap("tab10")
    colors = {regime: palette(i % 10) for i, regime in enumerate(regimes)}
    colors[branches[0]["regime"]] = "0.25"
    with PdfPages(path) as pdf:
        for trait in model["trait_names"]:
            fig, axes = plt.subplots(
                1, 2, figsize=(12, max(5, 0.18 * len(tips))), sharey=True, gridspec_kw={"width_ratios": [1.5, 1]}
            )
            for branch in order:
                row = branches[branch]
                if children[branch]:
                    axes[0].plot(
                        [x[branch]] * 2,
                        [min(y[c] for c in children[branch]), max(y[c] for c in children[branch])],
                        color="0.6",
                        linewidth=0.7,
                    )
                if branch:
                    axes[0].plot(
                        [x[row["parent"]], x[branch]], [y[branch]] * 2, color=colors[row["regime"]], linewidth=1.5
                    )
                    if branch in model["shift_branch_ids"]:
                        axes[0].scatter((x[row["parent"]] + x[branch]) / 2, y[branch], marker="|", color="black", s=70)
            rows = {row["branch_id"]: row for row in model["tip_predictions"] if row["trait"] == trait}
            for branch in tips:
                row = rows[branch]
                if row["observed"] is not None:
                    axes[1].errorbar(
                        row["observed"], y[branch], xerr=row["standard_error"], fmt="o", color="0.3", markersize=3
                    )
                axes[1].scatter(row["predicted"], y[branch], marker="x", color=colors[row["regime"]], s=25)
            axes[0].set_yticks([y[b] for b in tips], [branches[b]["name"] for b in tips], fontsize=7)
            axes[0].set_xlabel("Time from root")
            axes[1].set_xlabel(f"{trait}: observed ± sampling SE; × fitted mean")
            for ax in axes:
                ax.spines[["top", "right"]].set_visible(False)
            fig.suptitle(f"NWKIT OU shifts | {trait}")
            legend = [Line2D([0], [0], color=colors[regime], label=regime) for regime in regimes]
            legend.append(Line2D([0], [0], color="black", linestyle="", marker="|", label="shift location"))
            fig.legend(handles=legend, loc="lower center", ncol=min(5, len(legend)), fontsize=8)
            fig.tight_layout(rect=(0, 0.07, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--traits", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--replicate-separator", default="_")
    parser.add_argument("--resume-model")
    parser.add_argument("--regime-map", help="Optional fixed layout; skips shift selection.")
    parser.add_argument("--criterion", choices=["AIC", "AICc", "BIC", "pBIC", "bootstrap"], default="AICc")
    parser.add_argument(
        "--max-shifts", default="auto", help="Nonnegative integer or auto; NWKIT resolves the search cap."
    )
    parser.add_argument("--calibration-replicates", type=int, default=199)
    parser.add_argument("--calibration-level", type=float, default=0.05)
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--bootstrap-seed", type=int, default=2)
    parser.add_argument("--convergence", action="store_true")
    parser.add_argument("--estimate-measurement-error", choices=["yes", "no"], default="yes")
    parser.add_argument("--root-model", choices=["OUfixedRoot", "OUrandomRoot"], default="OUfixedRoot")
    parser.add_argument("--trait-covariance", choices=["diagonal", "full"], default="diagonal")
    parser.add_argument("--alpha-model", choices=["trait-specific", "shared"], default="trait-specific")
    parser.add_argument("--search-strategy", choices=["auto", "exhaustive", "lasso", "native-path"], default="auto")
    for flag in ("candidate-pool", "refit-budget", "screening-budget", "beam-width"):
        parser.add_argument("--" + flag, type=int, help="Optional override; omission uses the NWKIT default.")
    parser.add_argument("--alpha")
    parser.add_argument("--process-tip-variance")
    args = parser.parse_args(argv)
    from nwkit.cli import main as nwkit_main
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction, validate_output_targets

    prefix = Path(args.output_prefix).resolve()
    suffixes = [
        "model.json",
        "regime-map.tsv",
        "effects.tsv",
        "regimes.tsv",
        "tips.tsv",
        "replicates.tsv",
        "branch-summary.tsv",
        "pdf",
    ]
    destinations = {suffix: str(prefix) + "." + suffix for suffix in suffixes}
    validate_output_targets(list(destinations.values()))
    validate_outputs_do_not_replace_inputs(
        [
            (label, value)
            for label, value in [
                ("tree", args.tree),
                ("traits", args.traits),
                ("resume", args.resume_model),
                ("map", args.regime_map),
            ]
            if value
        ],
        list(destinations.items()),
    )
    table = pd.read_csv(args.traits, sep="\t", dtype={0: str})
    aggregated, traits, errors, audit = aggregate_replicates(table, args.replicate_separator)
    adapter = {
        "name": "genegalleon_native_ou",
        "raw_expression_sha256": hashlib.sha256(Path(args.traits).read_bytes()).hexdigest(),
        "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "replicate_separator": args.replicate_separator,
        "sampling_variance": "sample_variance/n; plug-in; zero known component for n=1, explicitly marked unavailable",
        "replicate_audit": audit,
    }
    if args.resume_model:
        previous = json.loads(Path(args.resume_model).read_text())
        if previous.get("genegalleon_adapter") != adapter:
            raise ValueError("Native resume rejected: replicate inputs or adapter changed.")
    with tempfile.TemporaryDirectory(prefix="gg-native-ou-") as directory:
        working = Path(directory)
        aggregated.to_csv(working / "traits.tsv", sep="\t", index=False, na_rep="NA")
        command = [
            "shift",
            "--selection",
            "native",
            "--infile",
            args.tree,
            "--format",
            "1",
            "--trait",
            str(working / "traits.tsv"),
            "--state-column",
            ",".join(traits),
            "--standard-error-column",
            ",".join(errors),
        ]
        for flag, suffix in [
            ("--model-out", "model.json"),
            ("--outfile", "regime-map.tsv"),
            ("--effects-out", "effects.tsv"),
            ("--regime-parameters-out", "regimes.tsv"),
            ("--tip-summary-out", "tips.tsv"),
        ]:
            command.extend([flag, str(working / suffix)])
        for name in (
            "max_shifts",
            "calibration_replicates",
            "calibration_level",
            "bootstrap",
            "seed",
            "bootstrap_seed",
            "root_model",
            "trait_covariance",
            "alpha_model",
            "search_strategy",
            "candidate_pool",
            "refit_budget",
            "screening_budget",
            "beam_width",
            "alpha",
            "process_tip_variance",
            "resume_model",
            "regime_map",
        ):
            value = getattr(args, name)
            if value is not None:
                command.extend(["--" + name.replace("_", "-"), str(value)])
        if args.criterion != "bootstrap":
            command.extend(["--criterion", args.criterion])
        if args.convergence:
            command.append("--convergence")
        if args.estimate_measurement_error == "yes":
            command.append("--estimate-measurement-error")
        nwkit_main(command)
        model = json.loads((working / "model.json").read_text())
        model["genegalleon_adapter"] = adapter
        (working / "model.json").write_text(json.dumps(model, indent=2, allow_nan=False) + "\n")
        pd.DataFrame(audit).to_csv(working / "replicates.tsv", sep="\t", index=False, na_rep="NA")
        branch_summary(model).to_csv(working / "branch-summary.tsv", sep="\t", index=False, na_rep="NA")
        plot_native_model(model, working / "pdf")
        with output_transaction(list(destinations.values())) as staged:
            for suffix, destination in destinations.items():
                shutil.copyfile(working / suffix, staged[destination])


if __name__ == "__main__":
    main()
