"""Run a large IQ-TREE ML search using a content-bound native checkpoint."""

import argparse
import gzip
import json
import shutil
import subprocess
import zlib
from pathlib import Path

from runtime_asset_io import atomic_json, file_sha


def read_receipt(path):
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, UnicodeError, json.JSONDecodeError):
        return None


def checkpoint_is_running(path):
    if not path.is_file():
        return False
    try:
        with gzip.open(path, "rb") as stream:
            if not stream.readline().startswith(b"--- # IQ-TREE Checkpoint"):
                return False
            # IQ-TREE refuses to rerun a completed search. Consume the stream
            # to validate its gzip trailer before handing it back to IQ-TREE.
            finished = False
            for line in stream:
                if line.startswith(b"finished:"):
                    finished = line.strip() != b"finished: false"
            return not finished
    except (EOFError, gzip.BadGzipFile, zlib.error):
        return False


def run_search(state, alignment, executable, options):
    state = state.absolute()
    alignment = alignment.resolve(strict=True)
    binary = Path(shutil.which(executable) or executable).resolve(strict=True)
    command = [str(binary), "-s", "alignment.fa", "--prefix", "search", *options]
    request = {"schema": 1, "alignment": file_sha(alignment),
               "iqtree": file_sha(binary), "argv": command}
    request_path = state / "request.json"
    complete_path = state / "complete.json"
    tree = state / "search.treefile"

    if state.is_symlink():
        raise ValueError(f"IQ-TREE checkpoint directory must not be a symlink: {state}")
    if read_receipt(request_path) == request:
        if tree.is_file() and tree.stat().st_size and read_receipt(complete_path) == file_sha(tree):
            print("IQ-TREE: reusing the completed matching search after interrupted publication.", flush=True)
            return 0
        if checkpoint_is_running(state / "search.ckp.gz"):
            print("IQ-TREE: resuming the matching native checkpoint.", flush=True)
        else:
            shutil.rmtree(state)
    elif state.exists():
        shutil.rmtree(state)

    if not state.exists():
        state.mkdir(parents=True, mode=0o700)
        atomic_json(request_path, request)
    alignment_link = state / "alignment.fa"
    alignment_link.unlink(missing_ok=True)
    alignment_link.symlink_to(alignment)
    complete_path.unlink(missing_ok=True)
    tree.unlink(missing_ok=True)
    returncode = subprocess.run(command, cwd=state, check=False).returncode
    if returncode != 0:
        return returncode if returncode > 0 else 128 - returncode
    if not tree.is_file() or not tree.stat().st_size:
        raise ValueError("IQ-TREE succeeded without a nonempty ML tree")
    atomic_json(complete_path, file_sha(tree))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--executable", default="iqtree")
    parser.add_argument("options", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.options or args.options[0] != "--":
        parser.error("pass native IQ-TREE options after --")
    try:
        return run_search(args.state, args.alignment, args.executable, args.options[1:])
    except (OSError, ValueError) as exc:
        parser.exit(1, f"IQ-TREE checkpoint error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
