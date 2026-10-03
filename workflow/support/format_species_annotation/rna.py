"""Extract annotated CDS from explicitly identified transcript FASTA inputs."""

import re
from collections import defaultdict

from format_species_writers import open_text

from .common import first_token, merge_coordinate_intervals, parse_gff_attributes
from .source_identity import source_annotation_path


def build_rna_coding_index(gff_path):
    features = {}
    aliases = defaultdict(set)
    rows = []
    with open_text(source_annotation_path(gff_path), "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            kind = parts[2].lower()
            attrs = parse_gff_attributes(parts[8])
            ids = attrs.get("ID", ())
            parents = tuple(attrs.get("Parent", ())) + tuple(attrs.get("Parent_Accession", ()))
            record = dict(kind=kind, parents=parents, seqid=parts[0], strand=parts[6],
                          start=int(parts[3]), end=int(parts[4]))
            if ids and kind not in ("cds", "exon"):
                identifier = ids[0]
                if identifier in features and features[identifier] != record:
                    raise ValueError("Conflicting RNA annotation feature: " + identifier)
                features[identifier] = record
                for key in ("ID", "Accession", "transcript_id", "Name"):
                    for alias in attrs.get(key, ()):
                        aliases[alias].add(identifier)
            if kind in ("cds", "exon"):
                rows.append(record)

    cache = {}

    def rna_parents(identifier):
        if identifier in cache:
            return cache[identifier]
        pending = [identifier]
        visited, result = set(), set()
        while pending:
            current = pending.pop()
            if current in visited:
                continue
            visited.add(current)
            for actual in aliases.get(current, {current}):
                record = features.get(actual)
                if record is None:
                    continue
                if record["kind"] in ("mrna", "transcript"):
                    result.add(actual)
                else:
                    pending.extend(record["parents"])
        cache[identifier] = result
        return result

    models = {identifier: dict(record, exon=[], cds=[]) for identifier, record in features.items()
              if record["kind"] in ("mrna", "transcript")}
    for record in rows:
        targets = set()
        for parent in record["parents"]:
            targets.update(rna_parents(parent))
        for target in targets:
            model = models[target]
            if (record["seqid"], record["strand"]) != (model["seqid"], model["strand"]):
                raise ValueError("RNA coding coordinates span different loci: " + target)
            model[record["kind"]].append((record["start"], record["end"]))
    return dict(models=models, aliases={key: tuple(sorted(value & models.keys()))
                                       for key, value in aliases.items() if value & models.keys()})


def extract_input_cds(task, header, sequence):
    """Return CDS, or None for a declared transcript with no coding annotation."""
    sequence = re.sub(r"\s+", "", sequence).upper()
    span = re.search(r"(?:^|\s)CDS=(\d+)-(\d+)(?:\s|$)", header)
    conversion = task.setdefault("_rna_conversion_audit", dict(trimmed=[], excluded_noncoding=[]))
    identifier = first_token(header)
    if span:
        start, end = map(int, span.groups())
        if not 1 <= start <= end <= len(sequence):
            raise ValueError("Transcript CDS span outside FASTA sequence: " + identifier)
        coding = sequence[start - 1:end]
    elif re.search(r"(?:^|\s)Type=(?:mRNA|transcript)(?:\s|$)", header, re.I):
        if not task.get("gff_path"):
            raise ValueError("Transcript FASTA requires GFF coding annotations: " + identifier)
        if "_rna_coding_index" not in task:
            task["_rna_coding_index"] = build_rna_coding_index(task["gff_path"])
        index = task["_rna_coding_index"]
        hits = index["aliases"].get(identifier, ())
        if len(hits) != 1:
            raise ValueError("Transcript FASTA has missing or ambiguous GFF RNA identity: " + identifier)
        model = index["models"][hits[0]]
        cds = merge_coordinate_intervals(model["cds"])
        if not cds:
            conversion["excluded_noncoding"].append(identifier)
            return None
        exons = merge_coordinate_intervals(model["exon"])
        if not exons or sum(end - start + 1 for start, end in exons) != len(sequence):
            raise ValueError("Transcript FASTA length disagrees with GFF exons: " + identifier)
        ordered_exons = list(reversed(exons)) if model["strand"] == "-" else exons
        slices, offset, covered = [], 0, 0
        for exon_start, exon_end in ordered_exons:
            for cds_start, cds_end in cds:
                start, end = max(exon_start, cds_start), min(exon_end, cds_end)
                if start > end:
                    continue
                relative = exon_end - end if model["strand"] == "-" else start - exon_start
                slices.append((offset + relative, offset + relative + end - start + 1))
                covered += end - start + 1
            offset += exon_end - exon_start + 1
        if covered != sum(end - start + 1 for start, end in cds):
            raise ValueError("GFF CDS lies outside transcript exons: " + identifier)
        coding = "".join(sequence[start:end] for start, end in sorted(slices))
    else:
        return sequence
    if len(coding) != len(sequence):
        conversion["trimmed"].append(dict(raw_id=identifier, transcript_length=len(sequence),
                                         cds_length=len(coding)))
    return coding
