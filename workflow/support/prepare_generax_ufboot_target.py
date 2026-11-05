#!/usr/bin/env python3
"""Write an unrooted, annotation-free GeneRax support-mapping target."""

import argparse


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, metavar="PATH")
    parser.add_argument("--output", required=True, metavar="PATH")
    return parser


from nwkit.tool_trees import prepare_target as prepare_target


def main():
    args = build_arg_parser().parse_args()
    prepare_target(args.input, args.output)


if __name__ == "__main__":
    main()
