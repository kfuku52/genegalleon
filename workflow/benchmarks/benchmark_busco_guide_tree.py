#!/usr/bin/env python3
"""Synthetic BUSCO guide scaling; run in a GeneGalleon runtime.
Excludes BUSCO execution. Synthetic topology cannot validate real plant relationships.
"""
import argparse
import csv
import gzip
import json
import math
import platform
import random
import resource
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.support import busco_guide_tree as guide
from workflow.support.input_generation_array_state import digest, digest_paths

AA = "ACDEFGHIKLMNPQRSTVWY"


def fixture(root, species, markers):
    for folder in ("cds", "full/single_copy", "short"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    rng = random.Random(92317)
    bases = ["M" + "".join(rng.choices(AA, k=299)) for _ in range(markers)]
    codons = dict(zip(AA, ("GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT",
                          "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"), strict=True))
    for i in range(species):
        name = f"Synthetic_species{i:04d}"
        records, dna, rows = {}, [], []
        for g, base in enumerate(bases):
            shared, private = random.Random(g*1000+i//5), random.Random(g*100000+i)
            seq = list(base)
            for position in shared.sample(range(1, 300), 24):
                seq[position] = shared.choice(AA)
            for position in private.sample(range(1, 300), 4):
                seq[position] = private.choice(AA)
            seq = "".join(seq)
            marker, identifier = f"BUSCO{g:05d}", f"{name}_g{g}"
            records[marker] = dict(sequence=seq, busco_sequence_id=identifier, protein_id=identifier, protein_header=identifier)
            rows.append(f"{marker}\tComplete\t{identifier}\t1\t300\n")
            dna.append(">"+identifier+"\n"+"".join(codons[c] for c in seq)+"TAA\n")
        cds, full, short = (root/"cds"/(name+".fa"), root/"full"/(name+".busco.full.tsv"), root/"short"/(name+".busco.short.txt"))
        cds.write_text("".join(dna))
        full.write_text("".join(rows))
        short.write_text("# BUSCO version is: synthetic\n"
                         "# The lineage dataset is: synthetic_odb12 (Creation date: 2026-01-01)\n"
                         "# BUSCO was run in mode: transcriptome\n"
                         f"C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:{markers}\n")
        with gzip.open(root/"full/single_copy"/(name+".json.gz"), "wt") as handle:
            json.dump(dict(schema=guide.SCHEMA, species=name, quality=guide.busco_quality(short), records=records,
                           marker_ids=sorted(records), source_hashes=digest_paths([cds, full, short])), handle)


def read_sketch(path):
    raw, position, values = path.read_bytes(), 32, []
    k, size, markers = struct.unpack_from("<QQQ", raw, 8)
    for _ in range(markers):
        count, = struct.unpack_from("<Q", raw, position)
        position += 8
        values.append(set(struct.unpack_from("<"+"Q"*count, raw, position)))
        position += count*8
    return k, size, values


def baseline(root, repeats, cpus):
    """Identical packed hashes and bottom-k union estimator, Python versus C++."""
    root.mkdir()
    rng, paths = random.Random(9), []
    for i in range(50):
        text, binary = root/f"{i}.txt", root/f"{i}.bin"
        text.write_text("40\n"+"\n".join("".join(rng.choices(AA, k=300)) for _ in range(40))+"\n")
        subprocess.run(["gg-kmer-distance", "sketch", str(text), str(binary), "5", "256"], check=True)
        paths.append(binary)
    listing = root/"sketches.txt"
    listing.write_text("\n".join(map(str, paths))+"\n")
    k, size, _ = read_sketch(paths[0])
    data = [read_sketch(path)[2] for path in paths]
    def distance(a, b):
        sample = set(sorted(a | b)[:size])
        j = len(a & b & sample)/len(sample)
        return min(1., -math.log(2*j/(1+j))/k) if j else 1.
    def native(threads=1):
        subprocess.run(["gg-kmer-distance", "compare", str(listing), str(root/"pairs.tsv"), str(threads)], check=True)
    def python():
        return {(i,j): sum(distance(a,b) for a,b in zip(data[i], data[j], strict=True))/40
                for i in range(len(data)) for j in range(i+1, len(data))}
    native()
    expected = python()
    with (root/"pairs.tsv").open() as handle:
        actual = {(int(row["i"]), int(row["j"])): float(row["distance"]) for row in csv.DictReader(handle, delimiter="\t")}
    error = max(abs(expected[key]-actual[key]) for key in expected)
    if error > 1e-12:
        raise AssertionError(f"Unequal distances: {error}")
    times = {}
    for label, run in (("python_1cpu", python), ("native_1cpu", native), ("native_parallel", lambda: native(cpus))):
        samples = []
        for _ in range(repeats):
            start = time.perf_counter()
            run()
            samples.append(time.perf_counter()-start)
        times[label] = samples
    return dict(species=50, markers=40, parallel_cpus=cpus, maximum_distance_error=error, seconds=times)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--species", type=int, nargs="+", default=[100, 250, 500])
    parser.add_argument("--markers", type=int, default=200)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cpus", type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(scope="synthetic proteins; excludes BUSCO; does not validate biological topology",
                  platform=platform.platform(), python=sys.version, cpus=args.cpus,
                  tools={t: dict(path=shutil.which(t), sha256=digest(shutil.which(t))) for t in ("rapidnj", "gg-kmer-distance")},
                  baseline=baseline(args.output/"equivalence", args.repeats, args.cpus), scaling={})
    for count in args.species:
        root = args.output/str(count)
        fixture(root, count, args.markers)
        runs = []
        for repeat in range(args.repeats):
            cache = root/f"cache{repeat}"
            for mode in ("cold", "warm"):
                options = guide.parser().parse_args(["build", "--cds-dir", str(root/"cds"),
                          "--full-dir", str(root/"full"), "--short-dir", str(root/"short"),
                          "--output", str(root/f"guide{repeat}-{mode}"), "--cache", str(cache),
                          "--markers", str(args.markers), "--cpus", str(args.cpus)])
                start = time.perf_counter()
                receipt = guide.build(options)
                runs.append(dict(repeat=repeat, mode=mode, wall_seconds=time.perf_counter()-start, **receipt["performance"]))
                print(json.dumps(runs[-1]), flush=True)
        result["scaling"][str(count)] = runs
        result["process_peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["child_peak_rss_kib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        (args.output/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print("Results:", args.output/"results.json")


if __name__ == "__main__":
    main()
