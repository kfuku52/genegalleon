"""Disambiguate reused GFF IDs only with explicit, local parent models.

A repeated CDS ID within one parent is a multipart feature, not a collision.
No parent models or coordinates are invented. Ordered fragments and biological
exceptions retain their source identities and their existing validators.
"""

import hashlib
import json
from collections import defaultdict
from pathlib import Path

from format_species_writers import open_text

from .common import parse_gff_attributes

VERSION = 1
ANCHORS = frozenset({"gene", "pseudogene", "mrna", "rna", "transcript", "primary_transcript"})
REFERENCES = frozenset({"Parent", "Derives_from", "gene_id", "gene", "locus_tag", "transcript_id"})
MARKER = "##genegalleon-locus-identity version="


class LocusIdentityError(ValueError):
    def __init__(self, audit):
        self.audit = audit
        sample = audit["problems"][0]
        super().__init__("Ambiguous GFF locus identity: {} at line {} ({})".format(
            sample["source_id"], sample["line"], sample["reason"]))


def rows(path):
    with open_text(path, "rt", errors="replace") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            fields = line.rstrip("\r\n").split("\t")
            if line.startswith("#") or len(fields) != 9:
                continue
            attrs = parse_gff_attributes(fields[8])
            yield number, fields, attrs


def geometry(fields):
    return fields[0], fields[6], int(fields[3]), int(fields[4])


def contains(parent, child):
    return parent[:2] == child[:2] and 1 <= parent[2] <= child[2] <= child[3] <= parent[3]


def suffix(identifier, scope):
    return identifier + ".locus" + hashlib.sha256(
        json.dumps((identifier, scope), separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()[:16]


def rewrite_attributes(text, replacements):
    """Replace parsed identity tokens, preserving unrelated source attributes."""
    from urllib.parse import quote, unquote
    changed, output = {}, []
    for field in text.split(";"):
        key, separator, value = field.partition("=")
        key = key.strip()
        if not separator or key not in {"ID", *REFERENCES}:
            output.append(field)
            continue
        values = value.split(",")
        updated = [replacements.get(unquote(item), unquote(item)) for item in values]
        if updated != [unquote(item) for item in values]:
            changed[key] = value
            field = key + "=" + ",".join(quote(item, safe=":._-|@") for item in updated)
        output.append(field)
    for key, value in changed.items():
        output.append("gg_original_" + key.lower() + "=" + value)
    if "ID" in changed:
        # The original alias remains available to uniquely identified protein
        # records. An ambiguous old FASTA ID must still fail, never pick a copy.
        for index, field in enumerate(output):
            if field.partition("=")[0].strip() == "Alias":
                output[index] = field + "," + changed["ID"]
                break
        else:
            output.append("Alias=" + changed["ID"])
    return ";".join(output)


def normalise_locus_ids(source, scratch):
    anchors, ids, protected, coarse_leaves = defaultdict(dict), set(), set(), defaultdict(set)
    for number, fields, attrs in rows(source):
        values = attrs.get("ID", ())
        ids.update(values)
        if (attrs.get("part") or attrs.get("number") or attrs.get("exception")
                or any(value.lower() == "true" for value in attrs.get("is_ordered", ()))):
            protected.update(values)
            protected.update(attrs.get("Parent", ()))
        if len(values) == 1 and fields[2].lower() in ANCHORS:
            scope = (fields[2].lower(), geometry(fields), tuple(sorted(attrs.get("Parent", ()))))
            anchors[values[0]].setdefault(scope, number)
        elif len(values) == 1:
            coarse_leaves[values[0]].add((fields[2].lower(), geometry(fields)[:2],
                                        tuple(sorted(attrs.get("Parent", ())))))

    while True:
        previous = len(protected)
        for identifier, models in anchors.items():
            if any(set(scope[2]) & protected for scope in models):
                protected.add(identifier)
        for identifier, models in coarse_leaves.items():
            if any(set(scope[2]) & protected for scope in models):
                protected.add(identifier)
        if len(protected) == previous:
            break

    changes, problems = {}, []

    def problem(identifier, number, reason):
        item = dict(source_id=identifier, line=number, reason=reason)
        if item not in problems:
            problems.append(item)

    for identifier, models in anchors.items():
        if len(models) < 2 or identifier in protected:
            continue
        kinds = {scope[0] for scope in models}
        if kinds == {"gene", "mrna"}:
            continue  # Existing exports may share their gene and RNA ID.
        if len(kinds) != 1:
            problem(identifier, min(models.values()), "conflicting feature types")
            continue
        spans = sorted(scope[1] for scope in models)
        if any(left[:2] == right[:2] and right[2] <= left[3]
               for left, right in zip(spans, spans[1:], strict=False)):
            problem(identifier, min(models.values()), "overlapping alternative definitions")
            continue
        for scope, number in models.items():
            changes[(identifier, scope)] = dict(source_id=identifier,
                target_id=suffix(identifier, scope), feature_type=scope[0],
                seqid=scope[1][0], strand=scope[1][1], start=scope[1][2], end=scope[1][3], line=number)

    if not changes and not problems and not any(len(models) > 1 and identifier not in protected
                                               for identifier, models in coarse_leaves.items()):
        return Path(source), dict(version=VERSION, status="unchanged", mappings=[], problems=[],
                                 original_sources_modified=False)
    del coarse_leaves

    def parent_scope(identifier, child):
        matches = [scope for scope in anchors.get(identifier, ()) if contains(scope[1], child)]
        return matches[0] if len(matches) == 1 else None

    def bindings(fields, attrs):
        return tuple((parent, parent_scope(parent, geometry(fields))) for parent in sorted(attrs.get("Parent", ())))

    leaves, coding = defaultdict(dict), set()
    for number, fields, attrs in rows(source):
        bound = bindings(fields, attrs)
        for parent, scope in bound:
            if scope is not None:
                if fields[2].lower() == "cds":
                    coding.add((parent, scope))
                    # Require the RNA's explicit gene parent as well. This
                    # proves a separate coding locus for every renamed gene.
                    for gene in scope[2]:
                        gene_scope = parent_scope(gene, scope[1])
                        if gene_scope is not None:
                            coding.add((gene, gene_scope))
        values = attrs.get("ID", ())
        if len(values) != 1 or fields[2].lower() in ANCHORS:
            continue
        scope = (fields[2].lower(), geometry(fields)[:2], bound)
        record = leaves[values[0]].setdefault(scope, dict(line=number, start=int(fields[3]), end=int(fields[4]), blocks=set()))
        record["start"], record["end"] = min(record["start"], int(fields[3])), max(record["end"], int(fields[4]))
        record["blocks"].add((int(fields[3]), int(fields[4]), fields[7]))
    for identifier, models in leaves.items():
        if len(models) < 2 or identifier in protected:
            continue
        physical_models = {(scope[0], scope[1], tuple(sorted(record["blocks"])))
                           for scope, record in models.items()}
        if len(physical_models) == 1:
            continue  # A single physical feature may have multiple parents.
        for scope, record in models.items():
            if not scope[2] or any(parent_model is None for _parent, parent_model in scope[2]):
                problem(identifier, record["line"], "missing, nonlocal or ambiguous parent model")
                continue
            changes[(identifier, scope)] = dict(source_id=identifier,
                target_id=suffix(identifier, scope), feature_type=scope[0],
                seqid=scope[1][0], strand=scope[1][1], start=record["start"], end=record["end"], line=record["line"])

    for (identifier, scope), item in changes.items():
        if scope[0] in ANCHORS and (identifier, scope) not in coding:
            problem(identifier, item["line"], "no unambiguous local coding descendants")
        if item["target_id"] in ids:
            problem(identifier, item["line"], "generated ID collides with a source ID")
    targets = [item["target_id"] for item in changes.values()]
    if len(targets) != len(set(targets)):
        problem("generated IDs", 0, "generated IDs are not unique")

    changed_ids = {identifier for identifier, _scope in changes}

    def check_ancestry(identifier, scope, seen):
        if (identifier, scope) in seen:
            problem(identifier, anchors[identifier][scope], "cyclic parent ancestry")
            return
        for parent in scope[2]:
            parent_model = parent_scope(parent, scope[1])
            if parent_model is None:
                problem(identifier, anchors[identifier][scope], "missing, nonlocal or ambiguous ancestor model")
            else:
                check_ancestry(parent, parent_model, seen | {(identifier, scope)})

    for identifier, scope in changes:
        if scope[0] in ANCHORS:
            check_ancestry(identifier, scope, set())
        else:
            for parent, parent_model in scope[2]:
                if parent_model is not None:
                    check_ancestry(parent, parent_model, set())

    def replacements(fields, attrs, number):
        result = {}
        values = attrs.get("ID", ())
        bound = bindings(fields, attrs)
        if len(values) == 1:
            identifier = values[0]
            scope = ((fields[2].lower(), geometry(fields), tuple(sorted(attrs.get("Parent", ()))))
                     if fields[2].lower() in ANCHORS else
                     (fields[2].lower(), geometry(fields)[:2], bound))
            item = changes.get((identifier, scope))
            if item:
                result[identifier] = item["target_id"]
        affected = bool(set(values) & changed_ids or set(attrs.get("Parent", ())) & changed_ids)
        for parent, scope in bound:
            item = changes.get((parent, scope))
            if item:
                result[parent] = item["target_id"]
            if affected and scope is None:
                problem(parent, number, "child is outside a unique explicit parent model")
        for key in REFERENCES - {"Parent"}:
            for identifier in attrs.get(key, ()):
                if identifier not in changed_ids or identifier in result:
                    continue
                scope = parent_scope(identifier, geometry(fields))
                item = changes.get((identifier, scope))
                if item:
                    result[identifier] = item["target_id"]
                else:
                    problem(identifier, number, "identity reference cannot be assigned to one locus")
        return result

    for number, fields, attrs in rows(source):
        if any(set(attrs.get(key, ())) & changed_ids for key in {"ID", *REFERENCES}):
            replacements(fields, attrs, number)
    audit = dict(version=VERSION, status="blocked" if problems else "repaired" if changes else "unchanged",
                 mappings=sorted(changes.values(), key=lambda item: (item["source_id"], item["target_id"])),
                 problems=problems, original_sources_modified=False)
    digest = hashlib.sha256()
    with Path(source).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    audit["source_sha256"] = digest.hexdigest()
    if problems:
        raise LocusIdentityError(audit)
    if not changes:
        return Path(source), audit
    output = Path(scratch) / (digest.hexdigest() + ".locus.gff")
    with open_text(source, "rt", errors="replace") as handle, output.open("w") as dest:
        dest.write("##gff-version 3\n")
        dest.write(MARKER + str(VERSION) + " source_sha256=" + digest.hexdigest() + "\n")
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##gff-version"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if not line.startswith("#") and len(fields) == 9:
                attrs = parse_gff_attributes(fields[8])
                fields[8] = rewrite_attributes(fields[8], replacements(fields, attrs, number))
                line = "\t".join(fields) + "\n"
            dest.write(line)
    return output, audit
