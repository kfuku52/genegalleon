"""Directory catalogs for retained workflow results; inputs and live caches stay raw."""

ORTHOFINDER_DIRECTORIES = tuple(
    prefix + name
    for prefix in ("", "core/")
    for name in (
        "Orthogroup_Sequences", "MultipleSequenceAlignments", "Gene_Trees",
        "Resolved_Gene_Trees", "Single_Copy_Orthologue_Sequences", "Orthologues",
        "Putative_Xenologs", "Gene_Duplication_Events",
    )
)
DIRECTORY_SCOPES = {
    # One immutable submission batch per root; scheduler writers do not take
    # GeneGalleon locks, so conversion additionally requires explicit quiescence.
    "scheduler_logs": ("scheduler_output",),
    "genome_evolution": (
        "busco_cds_fasta", "busco_mafft", "busco_trimal", "busco_iqtree_dna", "busco_iqtree_pep",
        "busco_reconciliation_dna", "busco_reconciliation_pep", "busco_rooted_txt_dna",
        "busco_rooted_txt_pep", "busco_rooted_nwk_dna", "busco_rooted_nwk_pep",
    ),
    "hgt": ("tree_plot", "tree_plot_input", "artifact_provenance", "trait_focus"),
    "analyses": (
        "inputs/families", "rendered_global/families", "rendered_gall/families",
        "package/provenance", "preparation",
    ),
}
ALL_DIRECTORIES = tuple(dict.fromkeys(ORTHOFINDER_DIRECTORIES + tuple(
    name for names in DIRECTORY_SCOPES.values() for name in names
)))


if __name__ == "__main__":
    # Exit 0 only for a potential managed directory. The owning provenance
    # command still verifies the marker and every declared content hash.
    import argparse
    import sys
    from pathlib import PurePath

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Potential managed output directory (exit 0 when recognized, otherwise 1).")
    args = parser.parse_args()
    sys.exit(0 if PurePath(args.path).name in {PurePath(name).name for name in ALL_DIRECTORIES} else 1)
