"""Independently check selected CDS against explicit source GFF gene parents.

This parser deliberately does not use the formatter's aggregation or rescue
index. Missing/ambiguous parent evidence is reported, never guessed from a
numeric transcript suffix. Supported source-identity normalization retains the
author's explicit identities when the original export's Parent labels are wrong.
"""

import re
from collections import defaultdict
from urllib.parse import unquote

from format_species_annotation.source_identity import source_annotation_path
from format_species_writers import open_text


def aliases(value):
    result = {value}
    while re.match(r"^(?:gene|rna|mrna|transcript|cds|protein)[:-]", value, re.I):
        value = re.sub(r"^[^:-]+[:-]", "", value, count=1)
        result.add(value)
    return result


class SourceGeneSelection:
    def __init__(self, gff_path):
        self.active = gff_path is not None
        self.parents = defaultdict(set)
        self.genes = set()
        self.nodes_by_alias = defaultdict(set)
        self.cache = {}
        self.output_roots = defaultdict(set)
        self.selected_roots = {}
        self.resolved_records = 0
        self.unresolved_records = 0
        if not self.active:
            return
        with open_text(source_annotation_path(gff_path), "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 9:
                    continue
                attrs = defaultdict(list)
                for item in fields[8].split(";"):
                    if "=" in item:
                        key, values = item.split("=", 1)
                        attrs[key.lower()].extend(unquote(value) for value in values.split(","))
                ids = attrs.get("id", [])
                parents = set(attrs.get("parent", []))
                kind = fields[2].lower()
                if kind in ("gene", "pseudogene"):
                    self.genes.update(ids)
                for identifier in ids:
                    self.parents[identifier].update(parents - {identifier})
                targets = parents if kind == "cds" and parents else set(ids)
                for key in ("id", "name", "alias", "accession", "protein_id", "transcript_id",
                            "orig_transcript_id", "orig_protein_id", "parent_accession", "cds"):
                    for value in attrs.get(key, []):
                        for alias in aliases(value):
                            self.nodes_by_alias[alias].update(targets)

    def roots(self, node):
        if node not in self.cache:
            found, seen, pending = set(), set(), [node]
            while pending:
                current = pending.pop()
                if current in seen:
                    continue
                seen.add(current)
                if current in self.genes:
                    found.add(current)
                else:
                    pending.extend(self.parents.get(current, ()))
            self.cache[node] = found
        return self.cache[node]

    def header_roots(self, header):
        if not self.active:
            return set()
        # Sequence/RNA identity takes precedence over a potentially erroneous
        # descriptive gene/locus tag (e.g. anonymous NCBI headers).
        primary = {header.split()[0]}
        local = re.match(r"^lcl\|.+_cds_(.+)_\d+$", header.split()[0])
        if local:
            primary.add(local.group(1))
        fields = header.split("||")
        if len(fields) >= 7 and fields[6].upper() == "CDS":
            primary.add(fields[4])
        for tag in ("protein_id", "transcript_id", "orig_protein_id", "orig_transcript_id"):
            match = re.search(r"\[" + tag + r"=([^\]]+)\]", header)
            if match:
                primary.add(match.group(1))
        for tokens in (primary, set(re.findall(r"\[gene=([^\]]+)\]", header))):
            found = set()
            for token in tokens:
                for alias in aliases(token):
                    for node in self.nodes_by_alias.get(alias, ()):
                        found.update(self.roots(node))
            if found:
                return found
        return set()

    def observe(self, header, output_id, selected):
        roots = self.header_roots(header)
        if self.active:
            if len(roots) == 1:
                self.resolved_records += 1
                self.output_roots[output_id].update(roots)
            else:
                self.unresolved_records += 1
        if selected:
            self.selected_roots[output_id] = roots

    def validate(self):
        retained = defaultdict(set)
        for output_id, roots in self.selected_roots.items():
            if len(roots) == 1:
                retained[next(iter(roots))].add(output_id)
        splits = {gene: sorted(ids) for gene, ids in retained.items() if len(ids) > 1}
        merges = {output_id: sorted(genes) for output_id, genes in self.output_roots.items() if len(genes) > 1}
        if splits or merges:
            raise ValueError("Source GFF gene ownership failed: retained_isoform_gene_groups={} "
                             "distinct_source_gene_merges={} sample={}".format(
                                 len(splits), len(merges), list(splits.items())[:3] + list(merges.items())[:3]))
        return dict(source_gene_check="explicit_parent" if self.active else "not_applicable",
                    source_gene_resolved_records=self.resolved_records,
                    source_gene_unresolved_records=self.unresolved_records,
                    source_gene_retained_isoform_groups=0, source_gene_distinct_merges=0)
