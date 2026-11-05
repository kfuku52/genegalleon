"""Calculate source identity afresh at each workflow command boundary."""

import argparse
import hashlib
import importlib.metadata
import importlib.util
from pathlib import Path


def package_identity(package, revision_file):
    spec = importlib.util.find_spec(package)
    if spec is None or spec.origin is None:
        identity = "unavailable"
    else:
        try:
            version = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            version = "unpackaged"
        root = Path(spec.origin).parent
        digest = hashlib.sha256()
        suffixes = {".py", ".R", ".r", ".tsv", ".nwk", ".json", ".txt", ".cpp"}
        for path in sorted(path for path in root.rglob("*") if path.is_file() and path.suffix in suffixes):
            digest.update(str(path.relative_to(root)).encode() + b"\0")
            digest.update(path.read_bytes())
        identity = f"{package} {version} source_sha256={digest.hexdigest()}"
    if revision_file.is_file() and revision_file.stat().st_size:
        for line in revision_file.read_text().splitlines():
            fields = line.split("\t")
            if fields[0] == package:
                identity += " " + (fields[1] if len(fields) > 1 else "")
                break
        else:
            identity += " "
    return identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package")
    parser.add_argument("revision_file", type=Path)
    args = parser.parse_args()
    print(package_identity(args.package, args.revision_file))


if __name__ == "__main__":
    main()
