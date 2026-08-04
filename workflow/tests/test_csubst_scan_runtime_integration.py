import ast
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

csubst = pytest.importorskip("csubst")
from csubst import ete, substitution_scan, tree  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
DB_SCRIPT = REPO_ROOT / "workflow" / "support" / "generate_orthogroup_database.py"
CSUBST_CORE_SCRIPT = REPO_ROOT / "workflow" / "core" / "gg_gene_evolution_core.sh"
CSUBST_SITE_WRAPPER = REPO_ROOT / "workflow" / "support" / "csubst_site_wrapper.py"
BASELINE_SCAN_COLUMNS = {
    "site_rate_categorized",
    "q_rate_enrichment_empirical",
    "q_rate_enrichment_empirical_by_trait",
    "q_rate_enrichment_empirical_by_trait_match",
}


def shell_csubst_command_options(subcommand):
    lines = CSUBST_CORE_SCRIPT.read_text(encoding="utf-8").splitlines()
    command_start = f"csubst {subcommand} \\"
    for start, line in enumerate(lines):
        if line.strip() == command_start:
            break
    else:
        raise AssertionError(f"CSUBST command not found: {subcommand}")

    options = set()
    for line in lines[start + 1 :]:
        match = re.match(r"\s*(--[A-Za-z][A-Za-z0-9_]*)", line)
        if match is not None:
            options.add(match.group(1))
        if not line.rstrip().endswith("\\"):
            break
    if subcommand == "search":
        options.update({"--foreground", "--fg_format"})
    return options


def csubst_sites_wrapper_options():
    tree = ast.parse(
        CSUBST_SITE_WRAPPER.read_text(encoding="utf-8"),
        filename=str(CSUBST_SITE_WRAPPER),
    )
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "build_csubst_sites_command"
    )
    return {
        node.value
        for node in ast.walk(function)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and re.fullmatch(r"--[A-Za-z][A-Za-z0-9_]*", node.value)
    }


@pytest.mark.parametrize("subcommand", ["search", "scan", "sites"])
def test_current_csubst_cli_supports_genegalleon_options(subcommand):
    proc = subprocess.run(
        ["csubst", subcommand, "--help-advanced"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    supported_options = set(
        re.findall(r"--[A-Za-z][A-Za-z0-9_]*", proc.stdout)
    )
    if subcommand == "sites":
        used_options = csubst_sites_wrapper_options()
    else:
        used_options = shell_csubst_command_options(subcommand)
    assert used_options <= supported_options, (
        f"GeneGalleon passes unsupported csubst {subcommand} options: "
        f"{sorted(used_options - supported_options)}"
    )


def set_state(state, branch_id, site, state_id):
    state[int(branch_id), int(site), :] = 0.0
    state[int(branch_id), int(site), int(state_id)] = 1.0


def toy_scan_context():
    phylogeny = tree.add_numerical_node_labels(
        ete.PhyloNode("((A:1,B:1)X:1,(C:1,D:1)Y:1)R;", format=1)
    )
    labels = {
        node.name: int(ete.get_prop(node, "numerical_label"))
        for node in phylogeny.traverse()
    }
    num_node = max(labels.values()) + 1
    for node in phylogeny.traverse():
        ete.set_prop(node, "SNdist", 0.25)
        ete.set_prop(node, "Ndist", 0.10)

    fg_leaf_names = {"trait": [["A"], ["C"]]}
    for index, names in enumerate(fg_leaf_names["trait"], start=1):
        name_set = set(names)
        for node in phylogeny.traverse():
            node_leaf_names = set(ete.get_leaf_names(node))
            ete.add_features(
                node,
                **{
                    f"is_lineage_fg_trait_{index}": node_leaf_names.issubset(
                        name_set
                    )
                },
            )
    for node in phylogeny.traverse():
        ete.add_features(node, is_fg_trait=False)

    state_nsy = np.zeros((num_node, 1, 2), dtype=float)
    state_pep = np.zeros((num_node, 1, 2), dtype=float)
    for node_id in labels.values():
        set_state(state_nsy, node_id, 0, 0)
        set_state(state_pep, node_id, 0, 0)
    for name in ["A", "C"]:
        set_state(state_nsy, labels[name], 0, 1)
        set_state(state_pep, labels[name], 0, 1)
        node = next(node for node in phylogeny.traverse() if node.name == name)
        ete.add_features(node, is_fg_trait=True)

    on_tensor = np.zeros((num_node, 1, 1, 2, 2), dtype=float)
    on_tensor[labels["A"], 0, 0, 0, 1] = 0.9
    on_tensor[labels["C"], 0, 0, 0, 1] = 0.8
    context = {
        "tree": phylogeny,
        "fg_df": pd.DataFrame({"name": ["A", "C"], "trait": [1, 2]}),
        "fg_leaf_names": fg_leaf_names,
        "fg_ids": {
            "trait": np.array([labels["A"], labels["C"]], dtype=np.int64)
        },
        "fg_stem_only": True,
        "scan_sister_stem_only": True,
        "state_nsy": state_nsy,
        "state_pep": state_pep,
        "nonsyn_state_orders": np.array(["A", "K"], dtype=object),
        "amino_acid_orders": np.array(["A", "K"], dtype=object),
        "iqtree_rate_values": np.array([0.25], dtype=float),
        "float_tol": 1e-12,
        "nonsyn_recode": "no",
        "scan_unit_mode": "clade",
        "scan_match": "any2spe",
        "scan_min_event_pp": 0.5,
        "scan_min_support": "2",
        "scan_rate_length": "raw",
        "scan_rate_exposure": "state_aware",
        "scan_rate_event_mode": "posterior_sum",
        "scan_other_scope": "all",
        "scan_pvalue_calibration": "none",
        "scan_n_permutations": 0,
        "scan_permutation_seed": 1,
        "scan_permutation_sample_original": False,
        "scan_permutation_retry_sample_original": True,
        "min_clade_bin_count": 1,
    }
    return context, on_tensor


def test_current_csubst_scan_output_imports_into_gene_family_database(tmp_path):
    context, on_tensor = toy_scan_context()
    scan_df, units_df = substitution_scan.scan_substitutions(
        g=context, ON_tensor=on_tensor
    )

    assert scan_df.shape[0] == 1
    assert BASELINE_SCAN_COLUMNS.issubset(scan_df.columns)
    assert "fg_clade_branch_ids" in units_df.columns
    assert scan_df.iloc[0]["site_rate"] == pytest.approx(0.25)

    stat_tree = tmp_path / "stat_tree"
    stat_branch = tmp_path / "stat_branch"
    scan_dir = tmp_path / "csubst_scan"
    units_dir = tmp_path / "csubst_scan_units"
    for directory in [stat_tree, stat_branch, scan_dir, units_dir]:
        directory.mkdir()
    pd.DataFrame(
        {
            "num_branch": [3],
            "num_spe": [1],
            "num_dup": [0],
            "num_sp": [2],
            "tree_metric": [1.0],
        }
    ).to_csv(
        stat_tree / "OG0001_stat.tree.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "branch_id": [0],
            "node_name": ["n0"],
            "num_sp": [2],
            "so_event": ["S"],
            "branch_metric": [1.0],
        }
    ).to_csv(
        stat_branch / "OG0001_stat.branch.tsv", sep="\t", index=False
    )
    scan_df.to_csv(scan_dir / "OG0001_csubst_scan.tsv", sep="\t", index=False)
    units_df.to_csv(
        units_dir / "OG0001_csubst_scan_units.tsv", sep="\t", index=False
    )

    db_path = tmp_path / "gg_orthogroup.db"
    proc = subprocess.run(
        [
            sys.executable,
            str(DB_SCRIPT),
            "--overwrite",
            "1",
            "--dbpath",
            str(db_path),
            "--dir_stat_tree",
            str(stat_tree),
            "--dir_stat_branch",
            str(stat_branch),
            "--dir_csubst_aa_change",
            str(scan_dir),
            "--dir_csubst_aa_change_unit",
            str(units_dir),
            "--ncpu",
            "1",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr

    with sqlite3.connect(db_path) as conn:
        db_scan = pd.read_sql_query("SELECT * FROM aa_change", conn)
        db_units = pd.read_sql_query("SELECT * FROM aa_change_unit", conn)
    assert db_scan.shape[0] == scan_df.shape[0]
    assert db_units.shape[0] == units_df.shape[0]
    assert set(scan_df.columns).issubset(db_scan.columns)
    assert set(units_df.columns).issubset(db_units.columns)
