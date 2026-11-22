"""Recover unambiguous missing FASTA record separators in bounded chunks."""

import gzip
import io
import os
import re
import tempfile
from contextlib import contextmanager
from pathlib import Path

BLOCK_SIZE = 1024 * 1024
_MAX_HEADER = 1024 * 1024
_DNA = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
_LINE_END = re.compile(r"[\r\n]")
_DNA_AND_LINE_BREAKS = ("".join(_DNA) + "\r\n").encode("ascii")
_NUMBERED_ID = re.compile(r"([A-Za-z_.-]*)([0-9]{1,18})([A-Za-z_.-]*)", re.ASCII)


def _successor(previous, token):
    """Match consecutive ordinals within the same numbered header family."""
    if len(previous) > 160 or len(token) > 160:
        return False
    a, b = _NUMBERED_ID.fullmatch(previous), _NUMBERED_ID.fullmatch(token)
    if not a or not b:
        return False
    prefix, number, suffix = a.groups()
    other_prefix, other_number, other_suffix = b.groups()
    return (prefix == other_prefix and suffix == other_suffix
            and len(number) == len(other_number) and int(other_number) == int(number) + 1)


def fasta_blocks(stream, repairs, *, chunk_size=BLOCK_SIZE, validate=False, limit_headers=True):
    """Yield original blocks, header spans, and offsets for evidenced insertions.

    Edits identify the decoded-stream offset, immediately after the original
    sequence. The stream's newline handling defines its offset convention.
    Permissive formatting can allow longer ordinary headers. Inline candidates
    always use the header limit.
    """
    offset = 0
    last = "\n"
    pending = ""
    previous = ""
    seen = set()
    duplicate = None
    awaiting = None
    has_sequence = False

    def check_sequence(text):
        nonlocal has_sequence
        if not text.isascii() or text.encode("ascii").translate(None, _DNA_AND_LINE_BREAKS):
            raise ValueError("Invalid DNA in FASTA after boundary repair")
        if not has_sequence and text.strip("\r\n"):
            if not previous:
                raise ValueError("Sequence precedes a FASTA ID after boundary repair")
            has_sequence = True

    while True:
        data = stream.read(chunk_size)
        end = not data
        data = pending + data
        pending = ""
        if not data:
            break
        headers, cuts = [], []
        cursor = 0
        while (marker := data.find(">", cursor)) >= 0:
            if awaiting:
                text = data[cursor:marker].lstrip()
                if not text or text[0] not in _DNA:
                    raise ValueError(f"Repaired FASTA contig {awaiting!r} has no initial DNA sequence")
                awaiting = None
            if validate:
                check_sequence(data[cursor:marker])
            flank = data[marker - 1] if marker else last
            inline = flank not in "\r\n"
            limit_header = limit_headers or inline
            newline = _LINE_END.search(data, marker + 1)
            if newline is None and not end:
                pending = data[marker:]
                if limit_header and len(pending) > _MAX_HEADER:
                    raise ValueError("Genome FASTA header exceeds 1 MiB")
                data = data[:marker]
                cursor = len(data)
                break
            stop = len(data) if newline is None else newline.start()
            if limit_header and stop - marker > _MAX_HEADER:
                raise ValueError("Genome FASTA header exceeds 1 MiB")
            header = data[marker + 1:stop].rstrip("\r")
            fields = header.split(None, 1)
            token = fields[0] if fields else ""
            if inline:
                at = offset + marker
                if (flank not in _DNA or header != token or token in seen
                        or not _successor(previous, token)):
                    raise ValueError(f"Ambiguous inline FASTA header at character {at} after {previous!r}")
                cuts.append(marker)
                repairs.append(dict(character_offset=at, previous_contig=previous, contig=token))
                awaiting = token
            if validate:
                if not token or (previous and not has_sequence):
                    raise ValueError("Missing FASTA ID or empty sequence after boundary repair")
                has_sequence = False
            if token in seen:
                duplicate = token
            seen.add(token)
            previous = token
            headers.append((marker, stop))
            cursor = stop + 1
        if validate:
            check_sequence(data[cursor:])
        if awaiting:
            text = data[cursor:].lstrip()
            if text:
                if text[0] not in _DNA:
                    raise ValueError(f"Repaired FASTA contig {awaiting!r} has no initial DNA sequence")
                awaiting = None
        if data:
            yield data, headers, cuts
            last = data[-1]
            offset += len(data)
        if end:
            break
    if awaiting:
        raise ValueError(f"Repaired FASTA contig {awaiting!r} is empty")
    if validate and not has_sequence:
        raise ValueError("Missing or empty last FASTA sequence after boundary repair")
    if (validate or repairs) and duplicate is not None:
        raise ValueError(f"Duplicate genome FASTA ID after boundary repair: {duplicate!r}")


def _repaired_chunks(stream, edits):
    for data, _, cuts in fasta_blocks(stream, edits, validate=True):
        start = 0
        for cut in cuts:
            yield data[start:cut] + "\n"
            start = cut
        yield data[start:]


@contextmanager
def _source_binary(path):
    with Path(path).open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
        raw.seek(0)
        if compressed:
            with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
                yield stream, True
        else:
            yield raw, False


def _may_need_repair(path):
    """Keep normal acquisition on a byte scan without decoding or recompressing."""
    with _source_binary(path) as (stream, _):
        previous = 10
        in_header = False
        while block := stream.read(BLOCK_SIZE):
            cursor = header = 0
            while (marker := block.find(b">", cursor)) >= 0:
                flank = block[marker - 1] if marker else previous
                if flank in (10, 13):
                    header = marker + 1
                    in_header = True
                elif (not in_header or block.find(b"\n", header, marker) >= 0
                        or block.find(b"\r", header, marker) >= 0):
                    return True
                cursor = marker + 1
            if in_header:
                in_header = block.find(b"\n", header) < 0 and block.find(b"\r", header) < 0
            previous = block[-1]
    return False


@contextmanager
def _source_text(path):
    with _source_binary(path) as (data, compressed):
        with io.TextIOWrapper(data, encoding="utf-8", newline="") as text:
            yield text, compressed


def repair_bound_genome(path):
    """Repair a previously checksummed private staging file, if necessary."""
    path = Path(path)
    if not _may_need_repair(path):
        return []
    edits = []
    with tempfile.TemporaryDirectory(prefix=".fasta-boundaries-", dir=path.parent) as directory:
        corrected = Path(directory) / "genome"
        with _source_text(path) as (source, compressed), corrected.open("wb") as output:
            if compressed:
                with gzip.GzipFile(filename="", fileobj=output, mode="wb", compresslevel=1, mtime=0) as archive:
                    with io.TextIOWrapper(archive, encoding="utf-8", newline="") as text:
                        text.writelines(_repaired_chunks(source, edits))
            else:
                with io.TextIOWrapper(output, encoding="utf-8", newline="") as text:
                    text.writelines(_repaired_chunks(source, edits))
        if not edits:
            raise ValueError("FASTA changed between boundary detection and repair")
        with corrected.open("rb") as verified:
            os.fchmod(verified.fileno(), 0o600)
            os.fsync(verified.fileno())
        corrected.replace(path)
    return edits


def repair_audit(source_sha, output_sha, edits, *, newline_mode):
    return dict(schema_version=1, kind="genegalleon-fasta-boundary-repair",
                source_sha256=source_sha, output_sha256=output_sha,
                offsets=newline_mode, insertions=edits)
