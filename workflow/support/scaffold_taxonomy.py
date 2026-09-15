#!/usr/bin/env python3
"""Pre-filter CDS taxonomy composition and HGT host-scaffold context.

Never infer origin from a best hit or treat an unresolved rank as compatible.
"""

import argparse
import gzip
from collections import Counter
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
METRICS = ("total_count", "compatible_count", "incompatible_count", "unresolved_count", "cds_id_count",
           "classified_fraction", "compatible_fraction", "compatible_all_fraction")
CONTEXT_COLUMNS = [f"host_scaffold_{background}{rank}_{metric}"
                   for background in ("", "background_") for rank in RANKS for metric in METRICS]
GENE_COLUMNS = ["host_scaffold_status", "host_scaffold_id", "host_scaffold_locus_id",
                "host_scaffold_count_unit", *CONTEXT_COLUMNS]
BRANCH_COLUMNS = ["host_scaffold_status", "host_scaffold_recipient_gene_count",
                  "host_scaffold_mapped_gene_count", "host_scaffold_count", *CONTEXT_COLUMNS]


def species_key(value):
    return str(value).strip().replace(" ", "_")


def composition(labels, count_units):
    counts = Counter(labels)
    total = sum(counts.values())
    compatible = counts["compatible"]
    incompatible = counts["incompatible"]
    classified = compatible + incompatible
    return dict(zip(METRICS, (total, compatible, incompatible, counts["unresolved"],
                             sum(unit == "cds_id" for unit in count_units),
                             classified / total if total else None,
                             compatible / classified if classified else None,
                             compatible / total if total else None), strict=True))


def gff_loci(path):
    """Exact GFF3 Parent / GTF gene_id relationships; no ID suffix heuristics."""
    parents, genes = {}, set()
    if not path:
        return {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed GFF record")
            attrs = {}
            for field in fields[8].split(";"):
                parts = field.strip().split("=", 1) if "=" in field else field.strip().split(None, 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip().strip('"')
            identifier = unquote(attrs.get("ID", ""))
            gene = unquote(attrs.get("gene_id", ""))
            transcript = unquote(attrs.get("transcript_id", ""))
            if fields[2].lower() in {"gene", "pseudogene"} and identifier:
                genes.add(identifier)
            if gene:
                genes.add(gene)
                if transcript:
                    parents.setdefault(transcript, set()).add(gene)
            if identifier and attrs.get("Parent"):
                parents.setdefault(identifier, set()).update(unquote(x) for x in attrs["Parent"].split(","))

    def ancestors(identifier, seen):
        if identifier in seen:
            raise ValueError("Cycle in GFF Parent relationships")
        if identifier in genes:
            return {identifier}
        return set().union(*(ancestors(p, seen | {identifier}) for p in parents.get(identifier, ())))

    result = {}
    for identifier in parents.keys() | genes:
        roots = ancestors(identifier, set())
        if len(roots) == 1:
            result[identifier] = next(iter(roots))
    return result


class RankResolver:
    def __init__(self, ncbi):
        self.ncbi = ncbi

    @lru_cache(maxsize=None)
    def ranks(self, taxid):
        if taxid <= 0:
            return {}
        try:
            lineage = self.ncbi.get_lineage(taxid)
        except ValueError:
            return {}
        ranks = self.ncbi.get_rank(lineage)
        result = {rank: tid for tid, rank in ranks.items()}
        if "domain" not in result and "superkingdom" in result:
            result["domain"] = result["superkingdom"]
        return result


def build_tables(gff_info, taxonomy, species, host_taxid, resolver, loci=None):
    """One row per input CDS identifier and rank; aggregates count loci once.

    Missing explicit locus ancestry uses the CDS identifier as its counting
    unit and is disclosed. Conflicting isoform labels make a locus unresolved.
    """
    loci = loci or {}
    required = {"gene_id", "chromosome"}
    if not required.issubset(gff_info.columns):
        raise ValueError("GFF info requires gene_id and chromosome")
    gff_info = gff_info.fillna("").drop_duplicates()
    if gff_info.gene_id.duplicated().any() or taxonomy.gene_id.duplicated().any():
        raise ValueError("Duplicate gene IDs in GFF info or CDS taxonomy")
    host = resolver.ranks(int(host_taxid))
    if not host:
        raise ValueError("Host taxid could not be resolved")
    assigned = taxonomy.set_index("gene_id").lca_taxid.to_dict()
    rows = []
    for row in gff_info.to_dict("records"):
        if not row["chromosome"] or row.get("splice_mode") == "trans-splicing":
            continue
        gene = str(row["gene_id"])
        locus = loci.get(str(row.get("gff_transcript_id", ""))) or loci.get(gene)
        rank_ids = resolver.ranks(int(assigned.get(gene, 0)))
        for rank in RANKS:
            h, t = host.get(rank), rank_ids.get(rank)
            label = "unresolved" if not h or not t else ("compatible" if h == t else "incompatible")
            rows.append(dict(species=species_key(species), gene_id=gene, scaffold=str(row["chromosome"]),
                             locus_id=locus or gene, count_unit="gff_locus" if locus else "cds_id",
                             rank=rank, host_taxid=h, label=label))
    columns = ["species", "gene_id", "scaffold", "locus_id", "count_unit", "rank", "host_taxid", "label"]
    genes = pd.DataFrame(rows, columns=columns)
    if genes.empty:
        return genes, pd.DataFrame(columns=["species", "scaffold", "rank", "host_taxid", *METRICS])
    # A locus on multiple scaffolds cannot be treated as single-scaffold evidence.
    if (genes.groupby("locus_id").scaffold.nunique() > 1).any():
        raise ValueError("GFF locus maps to multiple scaffolds")
    for _, group in genes.groupby(["scaffold", "locus_id", "rank"]):
        if group.label.nunique() != 1:
            genes.loc[group.index, "label"] = "unresolved"
    summaries = []
    for (scaffold, rank), group in genes.groupby(["scaffold", "rank"], sort=True):
        loci_group = group.drop_duplicates("locus_id")
        summaries.append(dict(species=species_key(species), scaffold=scaffold, rank=rank,
                              host_taxid=host.get(rank), **composition(loci_group.label, loci_group.count_unit)))
    return genes, pd.DataFrame(summaries)


def recipient_species(tree_path):
    if not tree_path or not Path(tree_path).is_file():
        return {}
    from Bio import Phylo
    tree = Phylo.read(tree_path, "newick")
    result = {}
    for node in tree.find_clades():
        # Bio.Phylo parses numeric internal labels as confidence values.
        label = node.name
        if label is None and node.confidence is not None:
            label = format(node.confidence, "g")
        if label is None:
            continue
        key = species_key(label)
        if key in result:
            raise ValueError(f"Ambiguous species-tree label: {label}")
        result[key] = {species_key(tip.name) for tip in node.get_terminals()}
    return result


def attach_context(branches, genes, directory, tree_path):
    """Exclude the union of all candidate loci, across OGs, from background.

    Branch pools unique recipient scaffolds, counting each locus once. Gene
    rows retain all existing candidates and describe that gene's host scaffold.
    """
    branches = pd.concat([branches.drop(columns=BRANCH_COLUMNS, errors="ignore"),
                          pd.DataFrame(None, index=branches.index, columns=BRANCH_COLUMNS, dtype=object)], axis=1)
    genes = pd.concat([genes.drop(columns=GENE_COLUMNS, errors="ignore"),
                      pd.DataFrame(None, index=genes.index, columns=GENE_COLUMNS, dtype=object)], axis=1)
    for frame in (branches, genes):
        frame["host_scaffold_status"] = "missing_scaffold_taxonomy"
    files = sorted(Path(directory).glob("*_gene_taxonomy.tsv")) if directory else []
    if not files:
        return branches, genes
    required = {"species", "gene_id", "scaffold", "locus_id", "count_unit", "rank", "label", "host_taxid"}
    candidate_keys = {(species_key(g.gene_taxon), str(g.gene_id)) for g in genes.itertuples()}
    lookup, scaffold_metrics, seen_species = {}, {}, set()
    # Genome outputs are per species: never retain all species' full gene tables
    # together. Only candidate identities and candidate-scaffold totals survive.
    for path in files:
        data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if not required.issubset(data.columns):
            raise ValueError("Invalid scaffold gene taxonomy schema")
        data["species"] = data.species.map(species_key)
        species = set(data.species)
        if species & seen_species or data.duplicated(["species", "gene_id", "rank"]).any():
            raise ValueError("Duplicate species or species/gene/rank scaffold taxonomy rows")
        seen_species.update(species)
        if not set(data.label).issubset({"compatible", "incompatible", "unresolved"}):
            raise ValueError("Invalid scaffold taxonomy label")
        identities = data.drop_duplicates(["species", "gene_id"])
        matched = {(r.species, r.gene_id): r for r in identities.itertuples()
                   if (r.species, r.gene_id) in candidate_keys}
        lookup.update(matched)
        excluded = {(r.species, r.locus_id) for r in matched.values()}
        target_scaffolds = {(r.species, r.scaffold) for r in matched.values()}
        for key, group in data.groupby(["species", "scaffold"]):
            if key not in target_scaffolds:
                continue
            metrics = {}
            for rank in RANKS:
                ranked = group.loc[group["rank"].eq(rank)].drop_duplicates("locus_id")
                if ranked.empty:
                    raise ValueError(f"Missing rank {rank} for scaffold {key}")
                for background in (False, True):
                    selected = ranked.loc[[(key[0], locus) not in excluded for locus in ranked.locus_id]] if background else ranked
                    prefix = f"host_scaffold_{'background_' if background else ''}{rank}_"
                    metrics.update({prefix + k: v for k, v in composition(selected.label, selected.count_unit).items()})
            scaffold_metrics[key] = metrics
    for index, row in genes.iterrows():
        record = lookup.get((species_key(row.gene_taxon), str(row.gene_id)))
        if record is None:
            genes.at[index, "host_scaffold_status"] = "gene_not_mapped"
            continue
        for key, value in scaffold_metrics[(record.species, record.scaffold)].items():
            genes.at[index, key] = value
        for key, value in {"status": "measured", "id": record.scaffold, "locus_id": record.locus_id,
                           "count_unit": record.count_unit}.items():
            genes.at[index, "host_scaffold_" + key] = value
    recipients = recipient_species(tree_path)
    gene_groups = {key: group for key, group in genes.groupby("orthogroup")}
    for index, branch in branches.iterrows():
        transfer = str(branch.generax_transfer).split("@")
        taxa = recipients.get(species_key(transfer[2])) if len(transfer) == 3 and transfer[0] == "Y" else None
        if taxa is None:
            branches.at[index, "host_scaffold_status"] = "recipient_unresolved"
            continue
        group = gene_groups.get(branch.orthogroup, genes.iloc[:0])
        ids = {g.strip() for g in str(branch.candidate_genes).split(";")}
        selected = group.loc[group.gene_id.isin(ids) & group.gene_taxon.map(species_key).isin(taxa)]
        mapped = [lookup[(species_key(g.gene_taxon), str(g.gene_id))] for g in selected.itertuples()
                  if (species_key(g.gene_taxon), str(g.gene_id)) in lookup]
        keys = {(r.species, r.scaffold) for r in mapped}
        branches.at[index, "host_scaffold_recipient_gene_count"] = len(selected)
        branches.at[index, "host_scaffold_mapped_gene_count"] = len(mapped)
        branches.at[index, "host_scaffold_count"] = len(keys)
        branches.at[index, "host_scaffold_status"] = ("measured" if len(mapped) == len(selected) else "partial") if keys else "no_recipient_scaffold"
        if not keys:
            continue
        for background in ("", "background_"):
            for rank in RANKS:
                prefix = f"host_scaffold_{background}{rank}_"
                # Pool counts over unique scaffolds, then recompute fractions.
                totals = {name: sum(scaffold_metrics[key][prefix + name] for key in keys) for name in METRICS[:5]}
                n, c, i = totals["total_count"], totals["compatible_count"], totals["incompatible_count"]
                totals.update(classified_fraction=(c + i) / n if n else None,
                              compatible_fraction=c / (c + i) if c + i else None,
                              compatible_all_fraction=c / n if n else None)
                for name, value in totals.items():
                    branches.at[index, prefix + name] = value
    return branches, genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff-info", required=True)
    parser.add_argument("--gff", default="")
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--host-taxid", type=int)
    parser.add_argument("--taxonomy-dbfile", required=True)
    parser.add_argument("--gene-out", required=True)
    parser.add_argument("--scaffold-out", required=True)
    args = parser.parse_args()
    from ete4 import NCBITaxa
    if not Path(args.taxonomy_dbfile).is_file():
        raise FileNotFoundError(args.taxonomy_dbfile)
    ncbi = NCBITaxa(dbfile=args.taxonomy_dbfile)
    taxid = args.host_taxid
    if taxid is None:
        name = args.species.replace("_", " ")
        ids = ncbi.get_name_translator([name]).get(name, [])
        if len(ids) != 1:
            raise ValueError(f"Host name must resolve uniquely; supply --host-taxid: {name}")
        taxid = ids[0]
    taxonomy = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str, keep_default_na=False)
    if taxonomy.shape[1] != 9:
        raise ValueError("Expected 9-column MMseqs2 CDS taxonomy output")
    taxonomy = taxonomy.iloc[:, :2].set_axis(["gene_id", "lca_taxid"], axis=1)
    taxonomy.lca_taxid = pd.to_numeric(taxonomy.lca_taxid, errors="raise").astype("int64")
    genes, scaffolds = build_tables(pd.read_csv(args.gff_info, sep="\t", dtype=str, keep_default_na=False),
                                    taxonomy, args.species, taxid, RankResolver(ncbi), gff_loci(args.gff))
    for frame, path in ((genes, args.gene_out), (scaffolds, args.scaffold_out)):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main()
