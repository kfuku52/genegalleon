"""Reuse exact codon alignments from an explicitly supplied WGD archive."""

import csv
import hashlib
import json
import os
import platform
import re
import shutil
import sys
from importlib.util import find_spec
from pathlib import Path

try:
    from pairwise_synteny import digest
    from pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS
except ImportError:
    from .pairwise_synteny import digest
    from .pairwise_synteny_ds import MAFFT_ARGS, PAIR_COLUMNS

ALIGNMENT_NAME = re.compile(r"cross_ks_aligned_pairs\.code([1-9][0-9]*)\.tsv\Z")
CODON_ALIGNMENT = re.compile(r"(?:[ACGTRYSWKMBDHVN]{3}|---)+\Z")
MAFFT_DEFAULT = re.compile(rb"^[ \t]*prefix=(/[a-zA-Z0-9_./-]+)[ \t]*$", re.MULTILINE)


def likelihood_identity(method):
    import cdskit
    import numpy
    import scipy
    import scipy.linalg
    import scipy.optimize

    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        return None, {}

    executable = shutil.which("cdskit")
    if executable is None:
        raise ValueError("Cross-block likelihood dS needs the CDSKIT executable")
    root = Path(cdskit.__file__).resolve().parent
    paths = [Path(executable).resolve(), Path(sys.executable).resolve(), *sorted(root.rglob("*.py"))]
    # Exact optimization/spectral kernels can change without a version change
    # in locally rebuilt numerical packages.
    for module in ("scipy.optimize._lbfgsb", "scipy.linalg._flapack", "scipy.linalg._matfuncs_expm"):
        spec = find_spec(module)
        if spec is None or spec.origin is None:
            return None, {}
        paths.append(Path(spec.origin).resolve())
    kernels = paths[-3:]
    libraries = threadpool_info()
    if not any(record.get("user_api") == "blas" and record.get("filepath") for record in libraries):
        return None, {}
    paths.extend(Path(record["filepath"]).resolve() for record in libraries if record.get("filepath"))
    files = {str(path): digest(path) for path in dict.fromkeys(paths)}
    identity = {"version": 1, "method": method, "python": sys.version, "numpy": numpy.__version__,
                "scipy": scipy.__version__, "cdskit": cdskit.__version__,
                "code": {str(path.relative_to(root)): files[str(path)] for path in paths if path.is_relative_to(root)},
                "platform": [platform.system(), platform.machine(), *platform.libc_ver()],
                "executable": files[str(paths[0])], "interpreter": files[str(paths[1])],
                "kernels": {path.name: files[str(path)] for path in kernels},
                "libraries": sorted(((record.get("internal_api"), record.get("architecture", ""), record.get("num_threads"),
                                      files[str(Path(record["filepath"]).resolve())])
                                     for record in libraries if record.get("filepath")), key=str),
                "environment": {key: os.environ.get(key) for key in ("OPENBLAS_CORETYPE", "NPY_DISABLE_CPU_FEATURES",
                                "MKL_CBWR", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")}}
    return json.loads(json.dumps(identity)), files


def alignment_identity(ds):
    """Fingerprint MAFFT and its bundled binaries, not merely the wrapper."""
    executable = shutil.which("mafft")
    if not executable:
        return None, ()
    mafft = Path(executable).resolve()
    override = os.environ.get("MAFFT_BINARIES")
    defaults = set(MAFFT_DEFAULT.findall(mafft.read_bytes())) if not override else ()
    if not override and len(defaults) != 1:
        return None, ()
    root = Path(override or os.fsdecode(next(iter(defaults)))).resolve()
    if not (root / "tbfast").is_file():
        return None, ()
    artifacts = list(sorted(p for p in root.iterdir() if p.is_file()))
    helpers = {p.name: digest(p) for p in artifacts}
    related = {Path(path).name: value for path, value in ds["source_hashes"].items()
               if Path(path).name in {"backalign.py", "codonutil.py", "translate.py", "util.py", "pairwise_synteny_ds.py"}}
    identity = {"version": 1, "args": list(MAFFT_ARGS), "mafft": ds["mafft_version"],
                "wrapper": ds["mafft_executable_sha256"], "helpers": helpers, "backtranslation": related,
                "platform": [platform.system(), platform.machine(), *platform.libc_ver()]}
    return identity, [(str(mafft), identity["wrapper"]), *((str(p.resolve()), helpers[p.name]) for p in artifacts)]


def reuse_sources(summary, identity):
    payload = summary.read_bytes()
    archived_hash = hashlib.sha256(payload).hexdigest()
    previous = json.loads(payload)
    if identity is None or previous.get("plan", {}).get("alignment_identity") != identity:
        raise ValueError("Prior WGD alignments do not have the same verified MAFFT and backtranslation identity")
    sources, hashes = {}, {str(summary): archived_hash}
    for name, expected in previous.get("output_hashes", {}).items():
        match = ALIGNMENT_NAME.fullmatch(name)
        if match is None:
            continue
        source = (summary.parent / name).resolve()
        if digest(source) != expected:
            raise ValueError(f"Prior WGD alignment hash does not match its summary: {name}")
        sources[match[1]] = str(source)
        hashes[str(source)] = expected
    if not sources:
        raise ValueError("Prior WGD summary has no cross-block codon alignments")
    return sources, hashes


def reuse_bundle(summaries, identity):
    bundles, hashes = {}, {}
    paths = tuple(dict.fromkeys(summaries))
    for summary in paths:
        sources, inputs = reuse_sources(summary, identity)
        for code, path in sources.items():
            aligned = bundles.setdefault(code, [])
            if path not in aligned:
                aligned.append(path)
        for path, value in inputs.items():
            if hashes.setdefault(path, value) != value:
                raise ValueError(f"Prior WGD alignment source changed between summaries: {path}")
    if len(paths) <= 1:
        # Retain the single-source representation in existing plans.
        bundles = {code: sources[0] for code, sources in bundles.items()}
    return bundles, hashes


def likelihood_reuse(summaries, identity, hashes):
    """Find verified whole-table reports, keyed by exact aligned input bytes."""
    sources = {}
    for summary in dict.fromkeys(summaries):
        payload = summary.read_bytes()
        if hashlib.sha256(payload).hexdigest() != hashes[str(summary)]:
            raise ValueError(f"Prior WGD summary changed while selecting dS reuse: {summary}")
        previous = json.loads(payload)
        if previous.get("plan", {}).get("cross_ks_likelihood_identity") != identity:
            continue
        output_hashes = previous["output_hashes"]
        for name, aligned_hash in output_hashes.items():
            match = ALIGNMENT_NAME.fullmatch(name)
            if match is None:
                continue
            code = match[1]
            report_name = f"cross_ks_pair_ds.code{code}.tsv"
            if report_name not in output_hashes:
                continue
            report = (summary.parent / report_name).resolve()
            expected = output_hashes[report_name]
            if digest(report) != expected or hashes.setdefault(str(report), expected) != expected:
                raise ValueError(f"Prior WGD dS report hash does not match its summary: {report}")
            by_alignment = sources.setdefault(code, {})
            record = {"path": str(report), "sha256": expected}
            if by_alignment.setdefault(aligned_hash, record)["sha256"] != expected:
                raise ValueError(f"Conflicting verified WGD dS reports for the same alignment: {report}")
    return sources


def reuse_likelihood_report(reuse, code, pairs_file, ds_file):
    if not reuse or not (entries := reuse.get(str(code))):
        return False
    record = entries.get(digest(pairs_file))
    if record is None:
        return False
    shutil.copyfile(record["path"], ds_file)
    if digest(ds_file) != record["sha256"]:
        raise ValueError("Prior WGD dS report changed while copying")
    return True


def matching_cached_rows(paths, selected, sequences):
    if isinstance(paths, (str, Path)):
        paths = (paths,)
    result = {}
    for path in paths:
        for pair, row in matching_rows(path, selected, sequences).items():
            if result.setdefault(pair, row) != row:
                raise ValueError(f"Conflicting verified WGD alignments for {row['pair_id']}: {path}")
    return result


def matching_rows(path, selected, sequences):
    requests = {pair: info for pair, info in selected}
    rows, seen = {}, set()
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if tuple(reader.fieldnames or ()) != PAIR_COLUMNS:
            raise ValueError(f"Invalid prior WGD alignment header: {path}")
        for row in reader:
            if None in row or any(val is None for val in row.values()):
                raise ValueError(f"Invalid prior WGD alignment table: {path}")
            pair = row["target_gene"], row["query_gene"]
            if row["pair_id"] != "|".join(pair) or pair in seen:
                raise ValueError(f"Duplicate or invalid pair in prior WGD alignments: {row['pair_id']}")
            seen.add(pair)
            if pair not in requests:
                continue
            info = requests[pair]
            first, second = row["sequence_1"], row["sequence_2"]
            if (not first or len(first) != len(second) or len(first) % 3
                    or any(not CODON_ALIGNMENT.fullmatch(seq) for seq in (first, second))):
                raise ValueError(f"Malformed prior WGD codon alignment: {row['pair_id']}")
            if all(seq.replace("-", "") == sequences[species][gene][0]
                   for seq, species, gene in zip((first, second), info[:2], pair, strict=True)):
                rows[pair] = row
    return rows
