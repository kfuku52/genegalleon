"""Explicit gene exclusions in a copied native orthogroup catalog.

The application supplies the exclusion list; this module never infers organelle
identity from a gene name. Original tables and every excluded membership remain
in the derived catalog. Unchanged families retain their verified annotations.
"""

import shutil
from pathlib import Path

import pandas


def read_excluded_genes(path):
    genes = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        gene = line.strip()
        if not gene:
            continue
        if any(char.isspace() or char == "," for char in gene):
            raise ValueError("Gene exclusions require one exact identifier per line")
        genes.append(gene)
    if not genes or len(genes) != len(set(genes)):
        raise ValueError("Gene exclusion list is empty or contains duplicate identifiers")
    return frozenset(genes)


def filter_members(members, counts, excluded):
    species = [name for name in members.columns if name != "Orthogroup"]
    if (members["Orthogroup"].duplicated().any()
            or counts["Orthogroup"].duplicated().any()
            or members["Orthogroup"].tolist() != counts["Orthogroup"].tolist()
            or not set(species + ["Total"]).issubset(counts.columns)):
        raise ValueError("Gene exclusion requires aligned unique native families and species counts")
    filtered = members.copy()
    updated = counts.copy()
    found, changes, changed = set(), [], set()
    member_rows = members[species].itertuples(index=False, name=None)
    count_rows = counts[species + ["Total"]].itertuples(index=False, name=None)
    for position, (cells, expected) in enumerate(zip(member_rows, count_rows, strict=True)):
        family = members["Orthogroup"].iloc[position]
        parsed = [([] if pandas.isna(cell) else [g.strip() for g in str(cell).split(",") if g.strip()]) for cell in cells]
        all_genes = [gene for cell in parsed for gene in cell]
        if (any(len(cell) != n for cell, n in zip(parsed, expected[:-1], strict=True))
                or len(all_genes) != expected[-1] or len(set(all_genes)) != len(all_genes)):
            raise ValueError("Native membership does not match its exact species counts: " + str(family))
        for column, genes in zip(species, parsed, strict=True):
            removed = [gene for gene in genes if gene in excluded]
            if not removed:
                continue
            changed.add(family)
            found.update(removed)
            kept = [gene for gene in genes if gene not in excluded]
            filtered.at[filtered.index[position], column] = ", ".join(kept)
            updated.at[updated.index[position], column] = len(kept)
            changes.extend({"Orthogroup": family, "Species": column, "GeneID": gene} for gene in removed)
        if family in changed:
            updated.at[updated.index[position], "Total"] = sum(updated.iloc[position][species])
    if found != set(excluded):
        raise ValueError("Excluded identifiers are absent from the native catalog: " + ", ".join(sorted(set(excluded) - found)))
    return filtered, updated, changed, pandas.DataFrame(changes, columns=["Orthogroup", "Species", "GeneID"])


def prepare_annotation(args, *, load_table, reuse_cache, protein_lengths, choose_representatives, annotate):
    directory = Path(args.dir_orthofinder_og)
    sources = [directory / "Orthogroups.tsv", directory / "Orthogroups.GeneCount.tsv"]
    backups = [directory / "Orthogroups.before_gene_exclusion.tsv",
               directory / "Orthogroups.GeneCount.before_gene_exclusion.tsv"]
    if any(path.is_symlink() or not path.is_file() for path in sources):
        raise ValueError("Gene exclusion requires regular copied native catalog tables")
    if any(path.exists() for path in backups):
        raise ValueError("Gene exclusions require a fresh copied catalog; preserved originals already exist")
    members = load_table(sources[0], gene_membership=True)
    counts = load_table(sources[1])
    filtered, updated, changed, report = filter_members(members, counts, read_excluded_genes(args.exclude_gene_ids))
    scope = getattr(args, "selection_scope", "")
    if scope:
        requested = load_table(scope)["Orthogroup"].tolist()
        if len(requested) != len(set(requested)) or not set(requested).issubset(set(counts["Orthogroup"])):
            raise ValueError("Selection scope must name unique existing native families")
    cache = getattr(args, "annotation_cache", "")
    if cache:
        cached = load_table(reuse_cache(args, cache))
        result = cached.copy()
        for column in counts.columns:
            result[column] = updated[column]
        annotation_columns = [column for column in result if column.startswith(("geneid_", "besthit_"))]
        result.loc[result["Orthogroup"].isin(changed), annotation_columns] = ""
        reannotate = updated["Orthogroup"].isin(changed) & (updated["Total"] > 0)
    else:
        result = updated.copy()
        reannotate = updated["Total"] > 0
    if reannotate.any():
        lengths = protein_lengths(args)
        if lengths.empty:
            raise ValueError("No protein lengths available for the filtered native families")
        representatives = choose_representatives(
            filtered.loc[reannotate].copy(), updated.loc[reannotate].copy(), lengths, args)
        replacements = annotate(representatives, args).set_index("Orthogroup")
        result = result.set_index("Orthogroup")
        for column in replacements:
            if column not in result:
                result[column] = ""
            result.loc[replacements.index, column] = replacements[column]
        result = result.reset_index()
    for source, backup in zip(sources, backups, strict=True):
        shutil.copy2(source, backup)
    filtered.to_csv(sources[0], sep="\t", index=False)
    updated.to_csv(sources[1], sep="\t", index=False)
    report.to_csv(directory / "Orthogroups.GeneExclusion.tsv", sep="\t", index=False)
    output = directory / "Orthogroups.GeneCount.annotated.tsv"
    result.to_csv(output, sep="\t", index=False)
    return str(output)
