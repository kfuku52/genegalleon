"""Directory catalogs for retained workflow results; inputs and live caches stay raw."""

from pathlib import PurePath


def is_transcriptome_quant_species(path):
    """The per-species output of the locked, quantification-only workflow."""
    path = PurePath(path)
    return path.parent.parts[-3:] == (
        "output", "transcriptome_quantification", "amalgkit_quant",
    )


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
    "cafe_plots": ("each_family_plot",),
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

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Potential managed output directory (exit 0 when recognized, otherwise 1).")
    args = parser.parse_args()
    recognized = (is_transcriptome_quant_species(args.path)
                  or PurePath(args.path).name in {PurePath(name).name for name in ALL_DIRECTORIES})
    sys.exit(0 if recognized else 1)
