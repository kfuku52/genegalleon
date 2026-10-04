"""Admission of formatted existing models to rescue protein anchors."""
import re
from collections import defaultdict
from dataclasses import replace

from kffractbias.io import ATTRIBUTE_PRIORITY, FEATURE_PRIORITY, annotation_to_genes

try:
    from cds_model_normalisation import CdsModelNormaliser, tokens
    from pairwise_synteny import prepare_genome
except ImportError:
    from .cds_model_normalisation import CdsModelNormaliser, tokens
    from .pairwise_synteny import prepare_genome


class AnchorAdmission(CdsModelNormaliser):
    def map_annotation(self, source, fasta_aliases):
        # GeneGalleon's formatter emits GeneID123 from NCBI GeneID:123 and
        # gene accessions from GWH. Resolve these own output conventions through
        # the public mapper's explicit feature/attribute interface.
        expanded = fasta_aliases | {"GeneID:" + alias[6:] for alias in fasta_aliases if re.fullmatch(r"GeneID[0-9]+", alias)}
        feature, attribute = source["feature"] or None, source["attribute"] or None
        if feature is None and attribute is None:
            attributes_order = (*ATTRIBUTE_PRIORITY, "Dbxref", "Accession")
            matches = defaultdict(set)
            for row in self.features:
                if row["feature"] not in FEATURE_PRIORITY:
                    continue
                for key in attributes_order:
                    matches[row["feature"], key].update(tokens(row["attributes"].get(key, "")) & fasta_aliases)
            if matches:
                feature, attribute = max(matches, key=lambda key: (len(matches[key]), -FEATURE_PRIORITY.index(key[0]),
                                                                  -attributes_order.index(key[1])))
                if not matches[feature, attribute]:
                    feature, attribute = None, None
        mapping = annotation_to_genes(source["gff"], expanded, feature=feature, attribute=attribute)
        def canonical(identifier):
            return identifier.replace("GeneID:", "GeneID", 1) if identifier.startswith("GeneID:") else identifier
        return replace(mapping, genes=tuple(replace(gene, gene_id=canonical(gene.gene_id)) for gene in mapping.genes),
                       locus_by_id={canonical(key): value for key, value in mapping.locus_by_id.items()})

    def __call__(self, gene, sequence, mapping):
        if self.lookup is None:
            self.configure(mapping)
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
