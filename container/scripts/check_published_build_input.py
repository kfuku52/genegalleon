#!/usr/bin/env python3
"""Compare published multi-architecture image labels with expected build inputs."""

import argparse
import json
import re
import sys

BUILD_INPUT_LABEL = "io.genegalleon.build-input"
REVISION_LABEL = "org.opencontainers.image.revision"
SHA256_RE = re.compile(r"[0-9a-f]{64}")
GIT_REVISION_RE = re.compile(r"[0-9a-f]{40}")


def parse_expected(values: list[str]) -> dict[str, str]:
    expected: dict[str, str] = {}
    for value in values:
        platform, separator, digest = value.partition("=")
        if not separator or not platform or not SHA256_RE.fullmatch(digest):
            raise ValueError(f"Expected build input must use PLATFORM=<64-character lowercase SHA-256>: {value!r}")
        if platform in expected:
            raise ValueError(f"Duplicate expected platform: {platform}")
        expected[platform] = digest
    if not expected:
        raise ValueError("At least one expected build input is required.")
    return expected


def compare_build_inputs(
    image_payload: object,
    expected: dict[str, str],
    label: str = BUILD_INPUT_LABEL,
) -> list[str]:
    if not isinstance(image_payload, dict):
        return ["Published image metadata is not a platform mapping."]

    problems: list[str] = []
    for platform, expected_digest in sorted(expected.items()):
        image = image_payload.get(platform)
        if not isinstance(image, dict):
            problems.append(f"Published image is missing platform {platform}.")
            continue
        config = image.get("config")
        labels = config.get("Labels") if isinstance(config, dict) else None
        observed = labels.get(label) if isinstance(labels, dict) else None
        if not isinstance(observed, str) or not SHA256_RE.fullmatch(observed):
            problems.append(f"Published image {platform} lacks a valid {label} label.")
        elif observed != expected_digest:
            problems.append(
                f"Published image {platform} build input changed: published={observed} expected={expected_digest}."
            )
    return problems


def common_revision(image_payload: object, platforms: list[str]) -> tuple[str | None, list[str]]:
    if not isinstance(image_payload, dict):
        return None, ["Published image metadata is not a platform mapping."]

    problems: list[str] = []
    revisions: dict[str, str] = {}
    for platform in sorted(set(platforms)):
        image = image_payload.get(platform)
        if not isinstance(image, dict):
            problems.append(f"Published image is missing platform {platform}.")
            continue
        config = image.get("config")
        labels = config.get("Labels") if isinstance(config, dict) else None
        observed = labels.get(REVISION_LABEL) if isinstance(labels, dict) else None
        if not isinstance(observed, str) or not GIT_REVISION_RE.fullmatch(observed):
            problems.append(f"Published image {platform} lacks a valid {REVISION_LABEL} label.")
        else:
            revisions[platform] = observed

    if problems:
        return None, problems
    if len(set(revisions.values())) != 1:
        return None, ["Published image revisions differ across expected platforms."]
    return next(iter(revisions.values())), []


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read `docker buildx imagetools inspect --format '{{json .Image}}'` "
            "JSON from stdin and compare io.genegalleon.build-input labels."
        )
    )
    parser.add_argument(
        "--label",
        default=BUILD_INPUT_LABEL,
        help=f"SHA-256 image label to compare (default: {BUILD_INPUT_LABEL}).",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--expected",
        action="append",
        metavar="PLATFORM=SHA256",
        help="Expected build-input label for one platform; repeat for every platform.",
    )
    mode.add_argument(
        "--print-common-revision",
        action="append",
        metavar="PLATFORM",
        help="Print the shared OCI revision for the listed platforms.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        image_payload = json.load(sys.stdin)
    except json.JSONDecodeError as exc:
        print(f"Published image build-input comparison failed: {exc}", file=sys.stderr)
        return 1

    if args.print_common_revision:
        revision, problems = common_revision(image_payload, args.print_common_revision)
        if problems:
            for problem in problems:
                print(problem, file=sys.stderr)
            return 1
        print(revision)
        return 0

    try:
        expected = parse_expected(args.expected)
    except ValueError as exc:
        print(f"Published image build-input comparison failed: {exc}", file=sys.stderr)
        return 1

    problems = compare_build_inputs(image_payload, expected, label=args.label)
    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1

    print(
        f"Published image {args.label} values match: " + ", ".join(sorted(expected)),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
