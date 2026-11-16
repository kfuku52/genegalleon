"""Use TAIR locus identifiers in the Arabidopsis TAIR12 annotation view."""

import re
import tempfile
from collections import defaultdict
from pathlib import Path

from format_species_writers import open_text

from .common import parse_gff_attributes
from .locus_identity import VERSION, LocusIdentityError, geometry, rewrite_attributes

TAIR12_FEATURE = re.compile(
    r"(?P<role>(?:(?:gene|rna|exon|cds|id)-)?)(?P<locus>(?:TAIR12_)+(?P<agi>AT[1-5CM]G[0-9]{5}))"
    r"(?P<suffix>(?:[.-][0-9]+)*)"
)


def identity_error(identifier, line, reason):
    return LocusIdentityError({
        "version": VERSION, "status": "blocked", "mappings": [],
        "problems": [{"source_id": identifier, "line": line, "reason": reason}],
        "original_sources_modified": False,
    })


def normalise_tair12_ids(source, scratch):
    """Apply a single checked mapping to feature IDs and their references."""
    declared_ids = set()
    candidates = {}
    genes = {}
    declared_gene_tokens = defaultdict(set)
    with open_text(source, "rt", errors="replace") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            is_gene = fields[2].lower() in {"gene", "pseudogene"}
            keys = {"ID", "locus_tag", "gene_id"} if is_gene else {"ID"}
            attrs = parse_gff_attributes(fields[8], keys=keys)
            for identifier in attrs.get("ID", ()):
                declared_ids.add(identifier)
                if is_gene:
                    for token in (*attrs.get("locus_tag", ()), *attrs.get("gene_id", ())):
                        declared_gene_tokens[token].add(identifier)
                match = TAIR12_FEATURE.fullmatch(identifier)
                if match is None:
                    continue
                role, locus, agi, suffix = match.group("role", "locus", "agi", "suffix")
                target = role + agi + suffix
                candidates.setdefault(identifier, (locus, target, number))
                if is_gene and role in {"", "gene-"} and not suffix:
                    model = (identifier, geometry(fields), agi)
                    if locus in genes and genes[locus] != model:
                        raise identity_error(identifier, number, "TAIR12 ID names distinct gene loci")
                    genes[locus] = model
    if not genes:
        return Path(source), []

    replacements = {}
    owners = {}
    for identifier, (locus, target, number) in candidates.items():
        if locus not in genes:
            continue
        gene_id, _geometry, agi = genes[locus]
        if identifier == gene_id:
            target = agi
        if (target in declared_ids or target in owners
                or (identifier == gene_id and declared_gene_tokens[target] - {identifier})):
            raise identity_error(identifier, number, "TAIR12 canonical ID collision: " + target)
        owners[target] = identifier
        replacements[identifier] = target
    for locus, (_gene_id, _geometry, agi) in genes.items():
        # NCBI GFF can add a feature-role prefix to ID but not to locus_tag.
        replacements.setdefault(locus, agi)

    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", prefix="tair12-", suffix=".gff", dir=scratch, delete=False,
    ) as dest:
        output = Path(dest.name)
        with open_text(source, "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    dest.write(line)
                    dest.writelines(handle)
                    break
                if "TAIR12_" in line and not line.startswith("#"):
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) == 9:
                        fields[8] = rewrite_attributes(fields[8], replacements)
                        line = "\t".join(fields) + "\n"
                dest.write(line)

    mappings = [{"source_id": gene_id, "target_id": agi} for gene_id, _, agi in sorted(genes.values())]
    return output, mappings
