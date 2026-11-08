"""Adapt ape's numbered edges to NWKIT's native split-based property transfer."""

import argparse
import json
import sys

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.numbered_tree import ape_tree, numbered_clade_keys
from nwkit.output_transaction import output_transaction
from nwkit.transfer import transfer_properties


def transfer_labels(data):
    target, target_nodes = ape_tree(data["target"])
    source, _ = ape_tree(data["source"])
    result = transfer_properties(
        target, source, [("name", "gg_support_label")], target_class="intnode",
        taxon_mode="intersection" if data["allow_subset"] else "exact",
        align_roots=False, allow_target_reroot=False, exclude_root=True,
        match_basis="split", root_edge_policies={"gg_support_label": "equal-only"},
        collect_report=False,
    )
    if result["status_counts"].get("ambiguous", 0):
        raise ValueError("Ambiguous support labels: NWKIT could not resolve equal values for a unique split")
    return [node.props.get("gg_support_label", "") for node in target_nodes]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", help="JSON tree-transfer or clade-index request")
    parser.add_argument("response", help="Output JSON path")
    args = parser.parse_args()
    validate_outputs_do_not_replace_inputs([("request", args.request)], [("response", args.response)])
    with open(args.request, encoding="utf-8-sig") as handle:
        data = json.load(handle)
    labels = (numbered_clade_keys(data["tree"], data.get("tip_set"), compact=True)
              if data.get("operation") == "clade_keys" else transfer_labels(data))
    with output_transaction([args.response], create_parents=True) as staged:
        with open(staged[args.response], "w", encoding="utf-8") as handle:
            json.dump(labels, handle)


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
