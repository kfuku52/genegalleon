"""Admission of formatted existing models to rescue protein anchors."""
from collections import defaultdict
from dataclasses import replace

from kffractbias.io import ATTRIBUTE_PRIORITY, FEATURE_PRIORITY, annotation_to_genes

try:
    from cds_model_normalisation import CdsModelNormaliser
    from pairwise_synteny import prepare_genome
except ImportError:
    from .cds_model_normalisation import CdsModelNormaliser
    from .pairwise_synteny import prepare_genome


class AnchorAdmission(CdsModelNormaliser):
    def map_annotation(self, source, fasta_aliases):
        feature, attribute = source["feature"], source["attribute"]
        if bool(feature) != bool(attribute):
            raise ValueError("GFF feature and attribute must be specified together")
        features = (feature,) if feature else FEATURE_PRIORITY
        attributes = (attribute,) if attribute else (*ATTRIBUTE_PRIORITY, "Dbxref", "Accession")
        feature_rank = {name: -i for i, name in enumerate(features)}
        attribute_rank = {name: -i for i, name in enumerate(attributes)}
        matches = defaultdict(dict)
        rows = defaultdict(list)
        for row in self.features:
            kind = row["feature"]
            if kind not in feature_rank:
                continue
            rows[kind].append(row)
            for key, value in row["attributes"].items():
                if key not in attribute_rank:
                    continue
                aliases = matches[kind, key]
                for raw in value.split(","):
                    # Match the formatter's GeneID spelling and the structural
                    # EMBL/NCBI gene prefix without replacing existing exact IDs.
                    if raw in fasta_aliases:
                        aliases[raw] = raw
                    elif raw.startswith("GeneID:") and "GeneID" + raw[7:] in fasta_aliases:
                        aliases[raw] = "GeneID" + raw[7:]
                    elif key == "Parent" and raw.startswith("gene-") and raw[5:] in fasta_aliases:
                        aliases[raw] = raw[5:]

        ranked = sorted((pair for pair in matches if matches[pair]), reverse=True,
                        key=lambda pair: (len(matches[pair]), feature_rank[pair[0]], attribute_rank[pair[1]]))
        rejected = []
        for kind, key in ranked:
            aliases = matches[kind, key]
            conflict = self.mapping_conflict(rows[kind], key, aliases)
            if conflict is not None:
                rejected.append(f"{kind}/{key}: {conflict}")
                continue
            mapping = annotation_to_genes(source["gff"], set(aliases), feature=kind, attribute=key)
            self.anchor_rows = rows[kind]
            self.anchor_aliases = {}
            for original, canonical in aliases.items():
                key = "GeneID" + original[7:] if original.startswith("GeneID:") else original
                if key != canonical:
                    self.anchor_aliases[key] = canonical
            if all(original == canonical for original, canonical in aliases.items()):
                return mapping
            return replace(mapping, genes=tuple(replace(gene, gene_id=aliases[gene.gene_id]) for gene in mapping.genes),
                           locus_by_id={aliases[identifier]: locus for identifier, locus in mapping.locus_by_id.items()})
        if rejected:
            raise ValueError("No unambiguous GFF-to-CDS mapping for " + str(source["gff"]) + " (" + "; ".join(rejected) + ")")
        raise ValueError(f"No GFF identifiers matched the CDS FASTA in {source['gff']}")

    def mapping_conflict(self, rows, attribute, aliases):
        if len(set(aliases.values())) != len(aliases):
            return "multiple GFF identifiers match the same CDS"
        positions = {}
        for row in rows:
            attr = row["attributes"]
            value = attr.get(attribute)
            if not value:
                continue
            for raw in value.split(","):
                if raw not in aliases:
                    continue
                if row["feature"] == "gene":
                    locus = attr.get("ID") or attr.get("gene_id") or raw
                else:
                    locus = attr.get("gene_id") or (raw if attribute == "Parent" else
                                                     attr.get("Parent") or attr.get("ID") or raw)
                    if row["feature"] == "CDS" and not attr.get("gene_id"):
                        # CDS parents are transcripts. Isoforms of a single
                        # gene can contribute to the same locus identifier.
                        roots = set()
                        for owner in locus.split(","):
                            matches = self.nodes.get(owner, ())
                            roots.update(r["attributes"].get("Parent", owner) for r in matches)
                            if not matches:
                                roots.add(owner)
                        if len(roots) != 1:
                            return f"{raw!r} maps to multiple parent loci"
                        [locus] = roots
                if "," in locus:
                    return f"{raw!r} maps to multiple parent loci"
                seqid, strand = row["seqid"], row["strand"]
                previous = positions.get(raw)
                if previous is not None:
                    other_seqid, other_strand, other_locus = previous
                    if (seqid != other_seqid or locus != other_locus
                            or (strand != other_strand and strand in "+-" and other_strand in "+-")):
                        return f"{raw!r} has conflicting sequence, strand or parent loci"
                    if strand == ".":
                        continue
                positions[raw] = (seqid, strand, locus)
        return None

    def __call__(self, gene, sequence, mapping):
        if self.lookup is None:
            self.configure(mapping, feature_rows=self.anchor_rows)
            for original, canonical in self.anchor_aliases.items():
                self.lookup[canonical] = self.lookup.pop(original, set())
        identifier = gene.gene_id
        alias = identifier.removeprefix(self.source["species"] + "_")
        keys = self.lookup.get(identifier, set()) | self.lookup.get(alias, set())
        return self.normalise(identifier, sequence, keys)


def prepare_rescue_genome(source, directory, side, minimum_mapping_fraction, *, required_ids=()):
    admission = AnchorAdmission(source, directory, side)
    try:
        genes, metadata = prepare_genome(source, directory, side, minimum_mapping_fraction,
                                        protein_transform=admission, annotation_mapper=admission.map_annotation)
        summary = admission.audit()
        if not genes:
            raise ValueError("No usable rescue anchors; see anchor_admission audit")
        required = set(required_ids)
        if any(row["status"] != "unchanged" for row in admission.rows if row["original_id"] in required):
            raise ValueError("New rescued model failed strict original-translation admission")
        if required - {row["original_id"] for row in admission.rows}:
            raise ValueError("New rescued model missing from annotation representatives")
        return genes, {**metadata, "anchor_admission": summary}
    finally:
        admission.close()
