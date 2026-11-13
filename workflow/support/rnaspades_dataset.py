"""Write the SPAdes YAML dataset (using YAML-compatible JSON syntax)."""

import json
import sys
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


if __name__ == "__main__":
    output, left_count, right_count, *paths = sys.argv[1:]
    left_count = int(left_count)
    right_end = left_count + int(right_count)
    Path(output).write_text(json.dumps(dataset(paths[:left_count], paths[left_count:right_end], paths[right_end:]),
                                      indent=2) + "\n", encoding="utf-8")
