#!/usr/bin/env python3
# coding: utf-8

"""Convert ASTRAL extended support annotations to a single selected label.

Example input annotation:
  '[q1=0.5;q2=0.25;q3=0.25;f1=1.0;...;pp1=0.45;...;QC=2;EN=2.0]'

This script replaces each annotation with one selected key value (e.g., q1).
It is robust to field order changes across ASTRAL versions.
"""

import argparse
import sys

from nwkit.astral import extract_support_labels as extract_support_labels
from nwkit.astral import parse_annotation as parse_annotation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", required=True, type=str, help="Input ASTRAL tree file")
    parser.add_argument("--outfile", required=True, type=str, help="Output tree file")
    parser.add_argument("--label_key", required=True, type=str, help="Label key to extract (e.g., q1, pp1)")
    args = parser.parse_args()

    with open(args.infile, "r", encoding="utf-8") as handle:
        text = handle.read()

    converted, missing_count = extract_support_labels(text, args.label_key)

    with open(args.outfile, "w", encoding="utf-8") as handle:
        handle.write(converted)

    if missing_count > 0:
        print(
            f"Warning: {missing_count} annotation blocks did not contain key '{args.label_key}'.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
