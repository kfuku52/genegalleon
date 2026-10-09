"""Read only declared coding/exon intervals into private disk scratch.

The compressed genome is scanned once with bounded fragments. Full chromosome
strings and a decompressed genome copy are never retained. The reader lives
only for one formatter invocation and does not weaken source/provenance checks.
"""
import bisect
import contextlib
import tempfile
from collections import defaultdict
from pathlib import Path

from .reference import (
    GenomeReferenceIndex,
    add_genome_reference,
    genome_fragments,
    genome_header_aliases,
)


def stat_identity(path):
    info = Path(path).stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def merged_intervals(intervals):
    result = []
    for start, end in sorted(set(intervals)):
        if start < 0 or end <= start:
            raise ValueError("Invalid requested genome interval")
        if result and start <= result[-1][1]:
            result[-1] = result[-1][0], max(end, result[-1][1])
        else:
            result.append((start, end))
    return result


class SequenceSlice:
    def __init__(self, reader, identifier):
        self.reader, self.identifier = reader, identifier

    def __len__(self):
        return self.reader.index[self.identifier]

    def __getitem__(self, key):
        if not isinstance(key, slice):
            raise TypeError("Genome interval access requires a slice")
        start, end, stride = key.indices(len(self))
        if stride != 1:
            raise ValueError("Genome interval access requires a contiguous slice")
        return self.reader.fetch(self.identifier, start, end)


class GenomeIntervals:
    def __init__(self, path, intervals, *, scratch_dir=None):
        self.path = Path(path)
        self.identity = stat_identity(self.path)
        self.index = GenomeReferenceIndex()
        self.canonical_ids = self.index.canonical_ids
        self.regions = {}
        self.starts = {}
        self.scratch = tempfile.TemporaryFile(dir=scratch_dir)
        try:
            self._scan(intervals)
        except BaseException:
            self.close()
            raise

    def _scan(self, intervals):
        requested = defaultdict(list)
        for seqid, start, end in intervals:
            requested[str(seqid)].append((int(start), int(end)))
        requested = {key: merged_intervals(value) for key, value in requested.items()}
        seen, header, token, length, selected, cursor = set(), None, "", 0, [], 0

        def finish():
            if header is not None:
                add_genome_reference(self.index, seen, header, length)

        for kind, text in genome_fragments(self.path):
            if kind == "header":
                finish()
                header, length, cursor = text, 0, 0
                token, aliases, _declared = genome_header_aliases(header)
                aliases |= {name.removeprefix("lcl|") for name in aliases}
                selected = merged_intervals(span for name in aliases for span in requested.get(name, ()))
                self.regions[token] = []
                self.starts[token] = []
            elif text:
                if header is None:
                    raise ValueError("Genome sequence precedes its FASTA header")
                if not text.isascii():
                    raise ValueError("Genome FASTA sequence must contain ASCII characters")
                end = length + len(text)
                while cursor < len(selected):
                    start, stop = selected[cursor]
                    if start >= end:
                        break
                    a, b = max(start, length), min(stop, end)
                    if a < b:
                        if start >= length:
                            self.regions[token].append((start, stop, self.scratch.tell()))
                            self.starts[token].append(start)
                        self.scratch.write(text[a - length:b - length].upper().encode("ascii"))
                    if stop > end:
                        break
                    cursor += 1
                length = end
        finish()
        if not seen:
            raise ValueError("Genome FASTA contains no records")
        if stat_identity(self.path) != self.identity:
            raise OSError("Genome changed while extracting declared intervals")
        self.scratch.flush()

    @property
    def references(self):
        return tuple(self.regions)

    @property
    def lengths(self):
        return tuple(self.index[name] for name in self.references)

    def get_reference_length(self, name):
        return self.index[name]

    def __contains__(self, name):
        return name in self.index

    def __getitem__(self, name):
        return SequenceSlice(self, name)

    def get(self, name, default=None):
        return self[name] if name in self else default

    def fetch(self, name, start, end):
        canonical = self.canonical_ids[name]
        if not 0 <= start <= end <= self.index[canonical]:
            raise ValueError("CDS/exon interval outside genome bounds: " + name)
        if start == end:
            return ""
        i = bisect.bisect_right(self.starts[canonical], start) - 1
        if i < 0 or end > self.regions[canonical][i][1]:
            raise ValueError("Genome access outside declared CDS/exon intervals: " + name)
        first, _last, offset = self.regions[canonical][i]
        self.scratch.seek(offset + start - first)
        sequence = self.scratch.read(end - start)
        if len(sequence) != end - start:
            raise OSError("Incomplete declared genome interval")
        return sequence.decode("ascii")

    def close(self):
        self.scratch.close()


@contextlib.contextmanager
def genome_intervals(path, intervals, *, scratch_dir=None):
    reader = GenomeIntervals(path, intervals, scratch_dir=scratch_dir)
    try:
        yield reader
        if stat_identity(path) != reader.identity:
            raise OSError("Genome source changed during interval use")
    finally:
        reader.close()


def cds_intervals(features_by_transcript):
    return ((row["seqid"], row["start"] - 1, row["end"])
            for rows in features_by_transcript.values() for row in rows)
