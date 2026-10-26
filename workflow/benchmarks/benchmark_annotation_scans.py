#!/usr/bin/env python3
"""Compare full source-identity/organelle scans with identical generated annotations."""

import argparse
import gzip
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def worker(a):
    sys.path.insert(0, str(a.support))
    from format_species_annotation.organelle import gff_organelle_seqids
    from format_species_annotation.source_identity import source_annotation_path

    start = time.perf_counter()
    cpu = time.process_time()
    result = gff_organelle_seqids(a.input)
    view = source_annotation_path(a.input)
    elapsed = time.perf_counter() - start
    cpu = time.process_time() - cpu
    payload = dict(seqids=sorted(result), view_sha256=sha(Path(view)))
    start = time.perf_counter()
    assert gff_organelle_seqids(a.input) == result
    cached = time.perf_counter() - start
    print(
        json.dumps(
            dict(
                wall_s=elapsed,
                cpu_s=cpu,
                cached_s=cached,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                outcome=payload,
            )
        )
    )


def fixtures(root, rows):
    cases = []
    for case, n, gtf in [("small", 32, False), ("ncbi", rows, False), ("gtf", max(32, rows // 20), True)]:
        p = root / (case + ".gff")
        data = ["##gff-version 3\n"]
        for i in range(n):
            feature = ["gene", "mRNA", "exon", "CDS"][i % 4]
            gene = i // 4
            start = gene * 12 + 1
            end = start + 8
            attrs = (
                f'gene_id "GENE{gene}"; transcript_id "RNA{gene}"; db_xref "GeneID:{gene},RefSeq:NM_{gene}"; description "ordinary; description, escaped "text"";'
                if gtf
                else f"ID={feature}{i};Parent=RNA{gene};Dbxref=GeneID:{gene},Genbank:NM_{gene};Name=model_{gene};gbkey={feature};gene=GENE{gene};product=ordinary%20protein%20{gene};inference=similar%20to%20RNA%20sequence;Note=ordinary%20note%2Cno%20author%20identity"
            )
            data.append(f"chr1\tbenchmark\t{feature}\t{start}\t{end}\t.\t+\t0\t{attrs}\n")
        data.append("chrM\tbenchmark\tregion\t1\t99\t.\t+\t.\tID=mt;genome=mitochondrion\n")
        data.append("chrP\tbenchmark\tCDS\t1\t9\t.\t+\t0\tID=cp;genome_type=chloroplast\n")
        p.write_text("".join(data))
        for codec in ["plain", "gzip"]:
            target = p if codec == "plain" else p.with_suffix(".gff.gz")
            if codec == "gzip":
                target.write_bytes(gzip.compress(p.read_bytes(), mtime=0))
            cases.append((case + "-" + codec, target, n + 2))
    return cases


def main(a):
    sources = {"baseline": a.baseline_support}
    if not a.baseline_only:
        sources["candidate"] = a.support
    trials = []
    summary = {}
    inputs = []
    with tempfile.TemporaryDirectory(prefix="gg-annotation-scans-") as tmp:
        cases = [] if a.existing_only else fixtures(Path(tmp), a.rows)
        if a.existing_gff:
            cases.append(("existing-gff", a.existing_gff, None))
        for name, path, rows in cases:
            inputs.append(dict(case=name, rows=rows, bytes=path.stat().st_size, sha256=sha(path)))
            expected = None
            for rep in range(-1, a.repeats):
                for label in list(sources) if rep % 2 == 0 else list(sources)[::-1]:
                    x = json.loads(
                        subprocess.check_output(
                            [
                                sys.executable,
                                __file__,
                                "--worker",
                                "--support",
                                str(sources[label]),
                                "--input",
                                str(path),
                            ],
                            text=True,
                        )
                    )
                    if name != "existing-gff":
                        assert x["outcome"]["seqids"] == ["chrM", "chrP"]
                    if expected is None:
                        expected = x["outcome"]
                    else:
                        assert x["outcome"] == expected, (name, label)
                    x.update(case=name, variant=label, warmup=rep < 0)
                    trials.append(x)
                    print(json.dumps(x), flush=True)
            summary[name] = {
                label: {
                    key: statistics.median(
                        x[key] for x in trials if x["case"] == name and x["variant"] == label and not x["warmup"]
                    )
                    for key in ["wall_s", "cpu_s", "cached_s", "peak_rss_bytes"]
                }
                for label in sources
            }
    implementations = {
        label: {
            str(Path("format_species_annotation") / name): sha(root / "format_species_annotation" / name)
            for name in ["common.py", "source_identity.py", "organelle.py"]
        }
        for label, root in sources.items()
    }
    report = dict(
        implementations=implementations,
        environment=dict(python=sys.version, platform=platform.platform()),
        configuration=dict(rows=a.rows, repeats=a.repeats, warmups=1),
        scope="Full native source identity and organelle scans only, cold Python caches in fresh alternating processes; warm filesystem cache. Complete source-view bytes and organelle sets compared; no whole scientific workflow throughput claim.",
        inputs=inputs,
        summary=summary,
        trials=trials,
    )
    a.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--support", type=Path, default=Path(__file__).resolve().parents[1] / "support")
    p.add_argument("--baseline-support", type=Path)
    p.add_argument("--baseline-only", action="store_true")
    p.add_argument("--rows", type=int, default=100000)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--output", type=Path)
    p.add_argument("--existing-gff", type=Path)
    p.add_argument("--existing-only", action="store_true")
    p.add_argument("--worker", action="store_true")
    p.add_argument("--input", type=Path)
    a = p.parse_args()
    if a.worker:
        worker(a)
    else:
        if a.baseline_support is None or a.output is None or min(a.rows, a.repeats) < 1:
            p.error("baseline support, output, and positive row/repeat counts are required")
        if a.existing_only and a.existing_gff is None:
            p.error("existing-only requires existing-gff")
        main(a)
