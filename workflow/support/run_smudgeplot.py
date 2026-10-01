#!/usr/bin/env python3
"""Count read k-mers with FastK and run Smudgeplot in its isolated runtime."""

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path


def positive_int(value):
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def validate_reads(values):
    reads = [Path(value).resolve(strict=True) for value in values]
    identities = set()
    for path in reads:
        stat = path.stat()
        if not path.is_file() or stat.st_size == 0:
            raise ValueError(f"every read input must be a nonempty file: {path}")
        identity = (stat.st_dev, stat.st_ino)
        if identity in identities:
            raise ValueError(f"duplicate read inputs would count the same reads twice: {path}")
        identities.add(identity)
    return reads


def discover_reads(directory):
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise ValueError(f"read input directory is not a directory: {root}")
    visited = set()
    candidates = []

    def fail(error):
        raise error

    for current, directories, filenames in os.walk(root, followlinks=True, onerror=fail):
        stat = Path(current).stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity in visited:
            raise ValueError(f"repeated or cyclic read input directory: {current}")
        visited.add(identity)
        directories[:] = sorted(name for name in directories if not name.startswith("."))
        for name in sorted(filenames):
            if not name.startswith(".") and name.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
                candidates.append(Path(current) / name)
    return sorted(validate_reads(candidates))


def validate_native_paths(paths):
    # Native FastK/Smudgeplot commands interpolate paths into shell commands.
    # Fail before invoking them rather than rewriting or patching dependencies.
    for path in paths:
        if shlex.quote(str(path)) != str(path):
            raise ValueError(
                "FastK/Smudgeplot native tools require paths without whitespace or "
                f"shell metacharacters: {path}"
            )


def run_analysis(args):
    reads = validate_reads(args.reads)
    read_stats = [path.stat() for path in reads]
    output = Path(args.output_dir).resolve()
    database_dir = Path(args.database_dir).resolve() if args.database_dir else None
    validate_native_paths([*reads, output, *([database_dir] if database_dir else [])])
    for name in ("FastK", "Histex", "Logex", "Symmex", "Fastrm", "smudgeplot"):
        if shutil.which(name) is None:
            raise ValueError(f"required command is missing: {name}; rebuild the GeneGalleon runtime")
    if database_dir == output:
        raise ValueError("output and database directories must be different")
    for path, label in ((output, "output"), (database_dir, "database")):
        if path is not None and path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"{label} directory must be empty: {path}")
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, MPLBACKEND="Agg")
    commands = []
    # Claim the empty output atomically, so concurrent runs cannot overwrite
    # one another's log or pair tables after passing the directory check.
    with (output / "commands.log").open("x") as log:

        def execute(command, stdout=None):
            commands.append(command)
            log.write(json.dumps(command) + "\n")
            log.flush()
            return subprocess.run(
                command, cwd=output, env=environment, check=True,
                stdin=subprocess.DEVNULL, stdout=stdout or log, stderr=log, text=True,
            )

        if database_dir:
            database_dir.mkdir(parents=True, exist_ok=True)
            workspace = nullcontext(str(database_dir))
        else:
            workspace = tempfile.TemporaryDirectory(prefix="fastk-", dir=output)
        with workspace as temporary:
            database = str(Path(temporary) / "reads")
            execute([
                "FastK", "-v", "-t1", f"-k{args.kmer_length}",
                f"-M{args.memory_gb}", f"-T{args.threads}",
                f"-P{temporary}", f"-N{database}", *map(str, reads),
            ])
            with (output / "fastk.histo.tsv").open("w") as histogram:
                execute(["Histex", "-G", database], stdout=histogram)
            lower = args.lower_count
            if lower == "auto":
                result = execute(
                    ["smudgeplot", "cutoff", "fastk.histo.tsv", "L"],
                    stdout=subprocess.PIPE,
                )
                try:
                    lower = positive_int(result.stdout.strip())
                except argparse.ArgumentTypeError as error:
                    raise ValueError("Smudgeplot cutoff did not return a positive integer") from error
            else:
                lower = positive_int(lower)
            execute([
                "smudgeplot", "hetmers", "-L", str(lower), "-t", str(args.threads),
                "-tmp", temporary, "-o", "kmerpairs", "--verbose", "--json_report", database,
            ])
        pairs = output / "kmerpairs.smu"
        if not pairs.is_file() or pairs.stat().st_size == 0:
            raise ValueError("Smudgeplot found no k-mer pairs; inspect the read coverage and cutoff")
        try:
            execute([
                "smudgeplot", "all", "kmerpairs.smu", "-o", "smudgeplot",
                "-d", str(args.aggregation_distance), "--format", "pdf", "--json_report",
                "--title", args.title,
            ])
        except subprocess.CalledProcessError as error:
            if args.lower_count == "auto":
                raise ValueError(
                    f"inference failed at automatic cutoff {lower}; inspect commands.log and "
                    "fastk.histo.tsv, assess read coverage and consider an explicit --lower-count; "
                    f"{error}"
                ) from error
            raise
    for filename in ("smudgeplot_smudgeplot.pdf", "smudgeplot_smudgeplot_log10.pdf", "smudgeplot_smudgeplot_report.json"):
        if not (output / filename).is_file() or (output / filename).stat().st_size == 0:
            raise ValueError(f"Smudgeplot did not produce the required output: {filename}")
    for filename in ("smudgeplot_smudgeplot.pdf", "smudgeplot_smudgeplot_log10.pdf"):
        with (output / filename).open("rb") as handle:
            if handle.read(5) != b"%PDF-":
                raise ValueError(f"Smudgeplot produced an invalid PDF: {filename}")
    report = json.loads((output / "smudgeplot_smudgeplot_report.json").read_text())
    if not isinstance(report, dict):
        raise ValueError("Smudgeplot inference report must be a JSON object")
    for path, original in zip(reads, read_stats, strict=True):
        current = path.stat()
        if (original.st_dev, original.st_ino, original.st_size, original.st_mtime_ns, original.st_ctime_ns) != (
            current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns, current.st_ctime_ns
        ):
            raise ValueError(f"read input changed during analysis: {path}")
    metadata = {
        "schema_version": 1,
        "reads": [{"path": str(path), "size_bytes": stat.st_size,
                   "mtime_ns": stat.st_mtime_ns} for path, stat in zip(reads, read_stats, strict=True)],
        "kmer_length": args.kmer_length, "lower_count_requested": args.lower_count,
        "lower_count_used": lower, "aggregation_distance": args.aggregation_distance,
        "threads": args.threads, "memory_gb": args.memory_gb,
        "commands": commands,
        "interpretation": "K-mer pair structure is supporting evidence; low read coverage can limit ploidy inference.",
    }
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--reads", nargs="+")
    inputs.add_argument("--list-reads", metavar="DIRECTORY", help="validate and list canonical FASTQ paths, separated by NUL")
    parser.add_argument("--output-dir")
    parser.add_argument("--database-dir", help="retain FastK files here for cutoff sensitivity analyses")
    parser.add_argument("--kmer-length", type=positive_int, default=21)
    parser.add_argument("--lower-count", default="auto", help="auto or positive count; inspect low-coverage data explicitly")
    parser.add_argument("--aggregation-distance", type=positive_int, default=2)
    parser.add_argument("--threads", type=positive_int, default=4)
    parser.add_argument("--memory-gb", type=positive_int, default=12)
    parser.add_argument("--title", default="Smudgeplot")
    args = parser.parse_args()
    if not args.list_reads and not args.output_dir:
        parser.error("--output-dir is required for analysis")
    if args.lower_count != "auto":
        try:
            positive_int(args.lower_count)
        except (ValueError, argparse.ArgumentTypeError) as error:
            parser.error(str(error))
    try:
        if args.list_reads:
            for path in discover_reads(args.list_reads):
                sys.stdout.buffer.write(os.fsencode(path) + b"\0")
        else:
            run_analysis(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Smudgeplot failed: {error}\n")


if __name__ == "__main__":
    main()
