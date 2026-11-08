"""Publish matching FASTA/Newick names using CDSKIT maps and NWKIT parsing."""

import argparse
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from cdskit.label_mapping import iter_fasta_headers, map_fasta_headers, shorten_labels, write_label_mapping
from nwkit.cli import parser
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.rename import rename_main
from nwkit.util import read_tree


def shorten(input_fasta, output_fasta, input_newick, output_newick, max_length):
    inputs = [Path(input_fasta), Path(input_newick)]
    outputs = [Path(output_fasta), Path(output_newick)]
    validate_outputs_do_not_replace_inputs(
        [("input", p) for p in inputs], [("output", p) for p in outputs])
    with inputs[0].open(encoding="utf-8", newline="") as fasta:
        labels = list(iter_fasta_headers(fasta))
        mapping = shorten_labels(labels, int(max_length))
        tree = read_tree(str(inputs[1]), 1, True, quiet=True)
        tips = list(tree.leaf_names())
        if len(tips) != len(set(tips)) or set(tips) != set(labels):
            raise ValueError("FASTA headers and tree tips must have matching unique labels")
        with output_transaction(outputs, create_parents=True) as staged:
            if not mapping:
                for source, output in zip(inputs, outputs, strict=True):
                    shutil.copyfile(source, staged[output])
                return
            fasta.seek(0)
            with Path(staged[outputs[0]]).open("w", encoding="utf-8", newline="") as output:
                output.writelines(map_fasta_headers(fasta, mapping))
            with TemporaryDirectory(prefix="gg-label-map-") as temporary:
                table = Path(temporary) / "mapping.tsv"
                write_label_mapping(mapping, table)
                args = parser.parse_args([
                    "rename", "--infile", str(inputs[1]), "--outfile", str(staged[outputs[1]]),
                    "--format", "1", "--outformat", "1", "--name-tsv", str(table),
                    "--target", "leaf", "--require-all-old-names", "yes"])
                rename_main(args)
            result = read_tree(str(staged[outputs[1]]), 1, True, quiet=True)
            expected = [mapping.get(label, label) for label in labels]
            if set(result.leaf_names()) != set(expected):
                raise ValueError("Renamed FASTA and tree disagree")


def main():
    arguments = argparse.ArgumentParser(description=__doc__)
    arguments.add_argument("input_fasta")
    arguments.add_argument("output_fasta")
    arguments.add_argument("input_newick")
    arguments.add_argument("output_newick")
    arguments.add_argument("max_length", type=int)
    args = arguments.parse_args()
    shorten(args.input_fasta, args.output_fasta, args.input_newick, args.output_newick, args.max_length)


if __name__ == "__main__":
    main()
