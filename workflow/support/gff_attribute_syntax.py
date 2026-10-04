"""Conservative attribute syntax repair at formatting, validation at consumption.

Only the numeric continuation of funannotate gene Name / mRNA product is
recoverable without guessing an attribute's meaning. Structural attributes are
never joined. Work on encoded text so existing escapes and list delimiters stay
unchanged.
"""
import gzip
import hashlib
import re
import shlex
from pathlib import Path

GFF_ATTRIBUTE_SYNTAX_VERSION = 1
KEY_VALUE = re.compile(r"^[^\s=;]+=")
NAME_CONTINUATION = re.compile(r"\d+(?:_\d+)?")
PRODUCT_CONTINUATION = re.compile(r"\d+(?:_\d+)?(?:,\s*variant\s+\d+)?")


def file_sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def validate_attributes(text):
    """Validate GFF key=value or GTF quoted attributes without decoding them."""
    if text in ("", "."):
        return
    if KEY_VALUE.match(text):
        for field in text.strip(";").split(";"):
            if field and not KEY_VALUE.match(field):
                raise ValueError("Expected a GFF key=value attribute: " + repr(field))
        return
    lexer = shlex.shlex(text, posix=True, punctuation_chars=";")
    lexer.whitespace_split = True
    lexer.commenters = ""
    tokens = iter(lexer)
    for key in tokens:
        if key == ";":
            continue
        value = next(tokens, None)
        if value is None or value == ";":
            raise ValueError("Missing GTF attribute value: " + key)
        separator = next(tokens, None)
        if separator not in (None, ";"):
            raise ValueError("Expected ';' after GTF attribute: " + key)


def normalise_attributes(text, source, feature, *, allow_bare=False):
    if allow_bare and text not in ("", ".") and not re.search(r"[\s=;]", text):
        return text  # Existing GFACS bare-ID adapter runs before output validation.
    if not KEY_VALUE.match(text):
        validate_attributes(text)
        return text
    repaired = []
    previous_key = None
    for field in text.split(";"):
        if not field:
            repaired.append(field)
            previous_key = None
        elif KEY_VALUE.match(field):
            repaired.append(field)
            previous_key = field.split("=", 1)[0]
        else:
            pattern = None
            if source.lower() == "funannotate":
                if feature.lower() == "gene" and previous_key == "Name":
                    pattern = NAME_CONTINUATION
                elif feature.lower() == "mrna" and previous_key == "product":
                    pattern = PRODUCT_CONTINUATION
            if pattern is None or pattern.fullmatch(field) is None:
                raise ValueError("Unrecoverable GFF attribute fragment: " + repr(field))
            repaired[-1] += "%3B" + field.replace(",", "%2C")
            previous_key = None  # Multiple orphan fields are ambiguous.
    result = ";".join(repaired)
    validate_attributes(result)
    return result


def normalise_line(line, path, line_number, changes):
    if not line.strip() or line.startswith("#"):
        return line
    parts = line.rstrip("\r\n").split("\t")
    if len(parts) != 9:
        raise ValueError(f"{path}:{line_number}: Expected nine GFF columns")
    before = parts[8]
    try:
        parts[8] = normalise_attributes(before, parts[1], parts[2], allow_bare=True)
    except ValueError as error:
        raise ValueError(f"{path}:{line_number}: {error}") from error
    if parts[8] == before:
        return line
    changes.append({"source_line": line_number, "before": before, "after": parts[8],
                    "reason": "escaped_funannotate_metadata_semicolon"})
    return "\t".join(parts) + line[len(line.rstrip("\r\n")):]


def validated_lines(lines, path):
    for line_number, line in enumerate(lines, 1):
        if line.startswith("##FASTA"):
            break
        if line.strip() and not line.startswith("#"):
            parts = line.rstrip("\r\n").split("\t")
            try:
                if len(parts) != 9:
                    raise ValueError("Expected nine GFF columns")
                validate_attributes(parts[8])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: {error}; regenerate with gg_input_generation") from error
        yield line


def validate_gff(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for _line in validated_lines(handle, path):
            pass


def syntax_audit(source_path, output_path, changes):
    return {"version": GFF_ATTRIBUTE_SYNTAX_VERSION, "source_sha256": file_sha256(source_path),
            "output_sha256": file_sha256(output_path), "changed_rows": len(changes), "changes": changes}
