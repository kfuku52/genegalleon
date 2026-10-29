"""Replayable text filtering for whole-line comments, preserving literal '#'."""

import gzip
import io


class LineCommentStream(io.TextIOBase):
    def __init__(self, source, *, quoted=False, gff=False):
        super().__init__()
        self.owns_source = not hasattr(source, "read")
        self.source = None
        if self.owns_source:
            opener = gzip.open if str(source).lower().endswith(".gz") else open
            source = opener(source, "rt", encoding="utf-8", newline="")
        self.source = source
        self.quoted = quoted
        self.gff = gff
        self.pending = ""
        self.in_quote = self.in_fasta = False

    def _quote_state(self, line):
        # Only delimiters and quote runs need inspection. Unquoted fields can
        # contain literal quotes; doubled quotes inside a quoted field escape.
        position = 0
        while position < len(line):
            if self.in_quote:
                position = line.find('"', position)
                if position < 0:
                    return
                position += 1
                if position < len(line) and line[position] == '"':
                    position += 1
                    continue
                self.in_quote = False
            elif line[position] == '"':
                self.in_quote = True
                position += 1
                continue
            position = line.find("\t", position)
            if position < 0:
                return
            position += 1

    def _line(self):
        for line in self.source:
            if self.in_fasta:
                continue  # Consume trailers/encoding through EOF, without tabulating FASTA.
            if self.gff and (line.rstrip("\r\n") == "##FASTA" or line.startswith(">")):
                self.in_fasta = True
                continue
            if not self.in_quote and line.startswith("#"):
                # Preserve physical row numbers in pandas parse diagnostics.
                return "\n"
            if self.quoted and (self.in_quote or '"' in line):
                self._quote_state(line)
            return line
        return ""

    def read(self, size=-1):
        if size == 0:
            return ""
        pieces, length = [self.pending], len(self.pending)
        self.pending = ""
        while size < 0 or length < size:
            line = self._line()
            if not line:
                break
            pieces.append(line)
            length += len(line)
        result = "".join(pieces)
        if size >= 0:
            result, self.pending = result[:size], result[size:]
        return result

    def readline(self, size=-1):
        if size == 0:
            return ""
        if not self.pending:
            self.pending = self._line()
        end = self.pending.find("\n") + 1 or len(self.pending)
        if size >= 0:
            end = min(end, size)
        result, self.pending = self.pending[:end], self.pending[end:]
        return result

    def seek(self, offset, whence=0):
        if offset != 0 or whence != 0:
            raise io.UnsupportedOperation("Only rewinding is supported")
        self.source.seek(0)
        self.pending = ""
        self.in_quote = self.in_fasta = False
        return 0

    def close(self):
        if not self.closed and self.owns_source and self.source is not None:
            self.source.close()
        super().close()
