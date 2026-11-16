"""Write the SPAdes YAML dataset (using YAML-compatible JSON syntax)."""

import argparse
import json
from pathlib import Path


def dataset(left, right, single):
    if len(left) != len(right):
        raise ValueError("rnaSPAdes input has unequal mate inventories")
    singles = {Path(path).resolve() for path in single}
    libraries = []
    for first, second in zip(left, right, strict=True):
        first, second = Path(first).resolve(), Path(second).resolve()
        library = {
            "type": "paired-end",
            "orientation": "fr",
            "left reads": [str(first)],
            "right reads": [str(second)],
        }
        companion = first.with_name(first.name.removesuffix("_1.amalgkit.fastq.gz") + ".amalgkit.fastq.gz")
        if companion in singles:
            library["single reads"] = [str(companion)]
            singles.remove(companion)
        libraries.append(library)
    libraries.extend({"type": "single", "single reads": [str(path)]} for path in sorted(singles))
    if not libraries:
        raise ValueError("rnaSPAdes input dataset is empty")
    return libraries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("left_count", type=int)
    parser.add_argument("right_count", type=int)
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()
    paths, left = args.paths, args.left_count
    right = left + args.right_count
    libraries = dataset(paths[:left], paths[left:right], paths[right:])
    args.output.write_text(json.dumps(libraries, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
