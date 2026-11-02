#!/usr/bin/env python3
"""Plan and retain GenomeScope fits with a traceable Smudgeplot initial lambda."""

import argparse
import hashlib
import json
import math
import shlex
import subprocess
import zipfile
from pathlib import Path


def integer_coverage(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError("Smudgeplot haploid_coverage must be a finite positive number")
    rounded = max(1, math.floor(value + 0.5))
    if rounded > 2147483647:
        raise ValueError("initial coverage exceeds GenomeScope's integer range")
    return rounded


def smudgeplot_coverage(archive, species, kmer_length, reads):
    prefix = f"{species}.smudgeplot/"
    with zipfile.ZipFile(archive) as handle:
        names = handle.namelist()
        for name in ("run_metadata.json", "smudgeplot_smudgeplot_report.json"):
            if names.count(prefix + name) != 1:
                raise ValueError(f"Smudgeplot archive must contain exactly one {name}")
        metadata = json.loads(handle.read(prefix + "run_metadata.json"))
        report = json.loads(handle.read(prefix + "smudgeplot_smudgeplot_report.json"))
    if not isinstance(metadata, dict) or not isinstance(report, dict):
        raise ValueError("Smudgeplot metadata and inference report must be JSON objects")
    if metadata.get("kmer_length") != kmer_length:
        raise ValueError("Smudgeplot and GenomeScope must use the same k-mer length")
    expected = []
    for read in reads:
        path = Path(read).resolve(strict=True)
        stat = path.stat()
        expected.append(dict(path=str(path), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))
    if not expected or metadata.get("reads") != expected:
        raise ValueError("Smudgeplot and GenomeScope must use the same unchanged canonical reads")
    h = report.get("haploid_coverage")
    integer_coverage(h)
    with Path(archive).open("rb") as archive_handle:
        digest = hashlib.file_digest(archive_handle, "sha256").hexdigest()
    return dict(haploid_coverage=h, archive=str(Path(archive).resolve()),
                archive_sha256=digest,
                kmer_length=kmer_length, reads=expected,
                lower_count_requested=metadata.get("lower_count_requested"),
                lower_count_used=metadata.get("lower_count_used"))


def plan(initial, source, archive, species, kmer_length, reads):
    if initial != -1 and not 1 <= initial <= 2147483647:
        raise ValueError("initial coverage must be -1 or a positive GenomeScope integer")
    if source not in ("smudgeplot", "auto"):
        raise ValueError("coverage source must be smudgeplot or auto")
    result = dict(schema_version=1, requested_initial_coverage=initial, coverage_source=source)
    if initial != -1:
        result.update(selected_source="manual", initial_coverage=initial, fits=[dict(label="primary", initial_coverage=initial)])
    elif source == "auto":
        result.update(selected_source="auto", initial_coverage=-1, fits=[dict(label="primary", initial_coverage=-1)])
    else:
        if not archive:
            raise ValueError("Smudgeplot-guided initialization requires matching reads and its archive; select auto or a manual lambda for histogram-only input")
        evidence = smudgeplot_coverage(archive, species, kmer_length, reads)
        h = evidence["haploid_coverage"]
        initial = integer_coverage(h)
        fits = [dict(label="primary", initial_coverage=initial)]
        for label, value in (("half_H", h / 2), ("twice_H", h * 2), ("auto", -1)):
            value = -1 if value == -1 else integer_coverage(value)
            if value not in [fit["initial_coverage"] for fit in fits]:
                fits.append(dict(label=label, initial_coverage=value))
        result.update(selected_source="smudgeplot", initial_coverage=initial, smudgeplot=evidence, fits=fits)
    return result


def validate_report(path):
    report = json.loads(path.read_text())
    if not isinstance(report, dict):
        raise ValueError("GenomeScope report must be a JSON object")
    summary = report.get("model_summary")
    length = report.get("genome_haploid_length")
    if not isinstance(summary, dict) or not summary or not isinstance(length, dict):
        raise ValueError("GenomeScope report lacks a fitted model summary or haploid length")
    values = (report.get("kcov"), length.get("avg"))
    if not all(isinstance(value, (int, float)) and not isinstance(value, bool)
               and math.isfinite(value) and value > 0 for value in values):
        raise ValueError("GenomeScope fitted coverage and haploid length must be finite and positive")


def run_fits(selection, output, histogram, program, kmer_length, ploidy, maximum):
    """The primary H fit stays at the archive root; alternatives never replace it."""
    if selection["selected_source"] == "smudgeplot":
        for record in selection["smudgeplot"]["reads"]:
            stat = Path(record["path"]).stat()
            if stat.st_size != record["size_bytes"] or stat.st_mtime_ns != record["mtime_ns"]:
                raise ValueError("Smudgeplot reads changed during GenomeScope counting; existing archive was preserved")
    directories = [output if fit["label"] == "primary" else output / "sensitivity" / fit["label"]
                   for fit in selection["fits"]]
    required = ("summary.txt", "transformed_linear_plot.png", "model.txt", "report.json")
    # Check every fit before executing any program. Old model files must never
    # make a zero-exit invocation that produced no output appear successful.
    for directory in directories:
        if any((directory / name).exists() for name in (*required, "command.txt", "fit.log")):
            raise ValueError(f"GenomeScope fit output already exists: {directory}; use a fresh output directory")
    results = []
    for fit, directory in zip(selection["fits"], directories, strict=True):
        directory.mkdir(parents=True, exist_ok=True)
        command = [program, "-i", str(histogram.resolve()), "-o", str(directory.resolve()),
                   "-k", str(kmer_length), "-p", str(ploidy), "-l", str(fit["initial_coverage"]),
                   "-m", str(maximum), "--json_report"]
        (directory / "command.txt").write_text(shlex.join(command) + "\n")
        result = dict(**fit, status="failed", directory=str(directory.relative_to(output)))
        try:
            with (directory / "fit.log").open("w") as log:
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
            result["returncode"] = completed.returncode
            if completed.returncode != 0:
                raise ValueError(f"GenomeScope exited with status {completed.returncode}; inspect fit.log")
            missing = [name for name in required if not (directory / name).is_file() or not (directory / name).stat().st_size]
            if missing:
                raise ValueError("GenomeScope did not produce fitted outputs: " + ", ".join(missing))
            validate_report(directory / "report.json")
            result["status"] = "success"
        except (OSError, ValueError) as error:
            result["error"] = str(error)
        results.append(result)
    selection.update(kmer_length=kmer_length, ploidy=ploidy, max_kmercov=maximum, results=results)
    (output / "coverage_selection.json").write_text(json.dumps(selection, indent=2, allow_nan=False) + "\n")
    if results[0]["status"] != "success":
        raise ValueError("GenomeScope primary fit failed; inspect fit.log and coverage_selection.json; existing archive was preserved")
    for result in results[1:]:
        if result["status"] != "success":
            print(f"GenomeScope sensitivity fit {result['label']} failed; diagnostic retained")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    configure = subparsers.add_parser("plan")
    configure.add_argument("--initial", type=int, required=True)
    configure.add_argument("--source", choices=("smudgeplot", "auto"), required=True)
    configure.add_argument("--archive")
    configure.add_argument("--species", required=True)
    configure.add_argument("--kmer-length", type=int, required=True)
    configure.add_argument("--reads", nargs="*", default=[])
    configure.add_argument("--output", type=Path, required=True)
    fit = subparsers.add_parser("fit")
    fit.add_argument("--selection", type=Path, required=True)
    fit.add_argument("--output", type=Path, required=True)
    fit.add_argument("--histogram", type=Path, required=True)
    fit.add_argument("--program", required=True)
    fit.add_argument("--kmer-length", type=int, required=True)
    fit.add_argument("--ploidy", type=int, required=True)
    fit.add_argument("--max-kmercov", type=int, required=True)
    args = parser.parse_args()
    if args.command == "plan":
        selection = plan(args.initial, args.source, args.archive, args.species, args.kmer_length, args.reads)
        args.output.write_text(json.dumps(selection, indent=2, allow_nan=False) + "\n")
    else:
        run_fits(json.loads(args.selection.read_text()), args.output, args.histogram,
                 args.program, args.kmer_length, args.ploidy, args.max_kmercov)


if __name__ == "__main__":
    main()
