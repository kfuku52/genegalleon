"""Independently check selected CDS against explicit source GFF/GTF gene parents.

This parser deliberately does not use the formatter's aggregation or rescue
index. Missing/ambiguous parent evidence is reported, never guessed from a
numeric transcript suffix. Supported source-identity normalization retains the
author's explicit identities when the original export's Parent labels are wrong.
"""

import re
import shlex
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
    def __init__(self, gff_path, *, repair_locus_ids=False):
        self.active = gff_path is not None
        self.parents = defaultdict(set)
        self.genes = set()
        self.nodes_by_alias = defaultdict(set)
        self.cds_parts = defaultdict(list)
        self.cds_location_nodes = None
        self.cache = {}
        self.output_roots = defaultdict(set)
        self.selected_roots = {}
        self.resolved_records = 0
        self.unresolved_records = 0
        if not self.active:
            return
        declared_nodes = set()
        missing_parent_axes = defaultdict(set)
        gtf_axes = defaultdict(set)
        gtf_gene_parents = defaultdict(set)
        with open_text(source_annotation_path(gff_path, repair_locus_ids=repair_locus_ids), "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 9:
                    continue
                attrs = defaultdict(list)
                is_gtf = re.match(r"^[^=\s;]+\s+\S", fields[8].strip()) is not None
                if is_gtf:
                    lexer = shlex.shlex(fields[8], posix=True, punctuation_chars=";")
                    lexer.whitespace_split = True
                    lexer.commenters = ""
                    tokens = iter(lexer)
                    for key in tokens:
                        if key == ";":
                            continue
                        value = next(tokens, None)
                        separator = next(tokens, None)
                        if value in (None, ";") or separator not in (None, ";"):
                            raise ValueError("Invalid source GTF attributes: " + fields[8])
                        attrs[key.lower()].append(value)
                else:
                    for item in fields[8].split(";"):
                        if "=" in item:
                            key, values = item.split("=", 1)
                            key = key.strip().lower()
                            # GFF ID is opaque, even when it contains commas.
                            items = [values] if key == "id" else values.split(",")
                            attrs[key].extend(unquote(value.strip()) for value in items)
                kind = fields[2].lower()
                if is_gtf:
                    gene_ids = set(attrs.get("gene_id", ())) - {""}
                    transcript_ids = set(attrs.get("transcript_id", ())) - {""}
                    if len(gene_ids) > 1 or len(transcript_ids) > 1 or gene_ids & transcript_ids:
                        raise ValueError("Ambiguous explicit source GTF identities")
                    axis = (fields[0], fields[6])
                    for identifier in gene_ids | transcript_ids:
                        gtf_axes[identifier].add(axis)
                    for transcript in transcript_ids:
                        gtf_gene_parents[transcript].update(gene_ids)
                    # GTF explicitly repeats this hierarchy on CDS/exon rows;
                    # no gene feature or numeric transcript suffix is required.
                    self.genes.update(gene_ids)
                    for gene in gene_ids:
                        self.nodes_by_alias[gene].add(gene)
                    attrs["id"] = list(gene_ids if kind in ("gene", "pseudogene") else transcript_ids)
                    attrs["parent"] = list(gene_ids if transcript_ids else ())
                ids = attrs.get("id", [])
                parents = set(attrs.get("parent", []))
                declared_nodes.update(ids)
                if kind.endswith("rna") or "transcript" in kind:
                    for parent in parents - set(ids):
                        missing_parent_axes[parent].add((fields[0], fields[6]))
                if kind in ("gene", "pseudogene"):
                    self.genes.update(ids)
                for identifier in ids:
                    self.parents[identifier].update(parents - {identifier})
                targets = parents if kind == "cds" and parents and not is_gtf else set(ids)
                if kind == "cds":
                    for target in targets:
                        self.cds_parts[(target, fields[0], fields[6])].append(
                            (int(fields[3]), int(fields[4])))
                for key in ("id", "name", "alias", "accession", "protein_id", "transcript_id",
                            "orig_transcript_id", "orig_protein_id", "parent_accession", "cds"):
                    for value in attrs.get(key, []):
                        for alias in aliases(value):
                            self.nodes_by_alias[alias].update(targets)
        if set(gtf_gene_parents) & self.genes:
            raise ValueError("Ambiguous source GTF gene/transcript identity collision")
        if any(len(axes) != 1 for axes in gtf_axes.values()):
            raise ValueError("Conflicting axes for explicit source GTF identity")
        if any(len(genes) > 1 for genes in gtf_gene_parents.values()):
            raise ValueError("Conflicting gene_id for explicit source GTF transcript_id")
        for parent, axes in missing_parent_axes.items():
            if parent in declared_nodes:
                continue
            if len(axes) != 1:
                raise ValueError("Conflicting axes for explicit missing GFF Parent: " + parent)
            self.genes.add(parent)

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
        location_roots = self.ncbi_location_roots(header)
        gene_tags = set(re.findall(r"\[gene=([^\]]+)\]", header))
        gene_tags.update(re.findall(r"(?<!\S)gene=([^\s\[\]]+)", header))
        for tokens in (primary, gene_tags):
            found = set()
            for token in tokens:
                for alias in aliases(token):
                    for node in self.nodes_by_alias.get(alias, ()):
                        found.update(self.roots(node))
            if len(location_roots) == 1 and (not found or location_roots <= found):
                return location_roots
            if len(found) == 1 and location_roots and found.isdisjoint(location_roots):
                raise ValueError("Source CDS identifier and exact GFF location disagree: " + header.split()[0])
            if found:
                return found
        return set()

    def ncbi_location_roots(self, header):
        """Disambiguate reused protein accessions with exact source CDS coordinates.

        Geometry only selects declared GFF ancestry; it never defines a gene or
        merges nearby/overlapping models. Unsupported locations remain unresolved.
        """
        token = header.split()[0]
        if not token.startswith("lcl|") or "_cds_" not in token:
            return set()
        match = re.search(r"\[location=([^\]]+)\]", header)
        if not match:
            return set()
        location = re.sub(r"\s+", "", match[1])
        complements = 0
        while True:
            wrapper = re.fullmatch(r"(complement|join|order)\((.*)\)", location, re.I)
            if not wrapper:
                break
            complements += wrapper[1].lower() == "complement"
            location = wrapper[2]
        intervals = []
        for piece in location.split(","):
            bounds = re.fullmatch(r"[<>]?(\d+)(?:\.\.[<>]?(\d+))?", piece)
            if not bounds:
                return set()
            start, end = int(bounds[1]), int(bounds[2] or bounds[1])
            if not 1 <= start <= end:
                return set()
            intervals.append((start, end))
        if self.cds_location_nodes is None:
            self.cds_location_nodes = defaultdict(set)
            for (node, seqid, strand), parts in self.cds_parts.items():
                self.cds_location_nodes[(seqid, strand, tuple(sorted(parts)))].add(node)
        signature = (token[4:].rsplit("_cds_", 1)[0],
                     "-" if complements % 2 else "+", tuple(sorted(intervals)))
        found = set()
        for node in self.cds_location_nodes.get(signature, ()):
            found.update(self.roots(node))
        return found

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
        complete = self.active and self.resolved_records > 0 and self.unresolved_records == 0
        status = ("explicit_parent" if complete else
                  "partial_explicit_parent" if self.resolved_records else "unresolved")
        return dict(source_gene_check=status if self.active else "not_applicable",
                    source_gene_complete=complete,
                    source_gene_resolved_records=self.resolved_records,
                    source_gene_unresolved_records=self.unresolved_records,
                    source_gene_unresolved_selected_records=sum(
                        len(roots) != 1 for roots in self.selected_roots.values()),
                    source_gene_retained_isoform_groups=0, source_gene_distinct_merges=0)
