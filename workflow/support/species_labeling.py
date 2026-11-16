"""Compatibility exports for NWKIT extended taxonomic labels."""

from nwkit.taxonomic_labels import (
    SPECIES_LABEL_TERMINAL_SUFFIXES as SPECIES_LABEL_TERMINAL_SUFFIXES,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_DISPLAY_RANKS as TAXONOMIC_DISPLAY_RANKS,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_GENUS_ONLY_PLACEHOLDERS as TAXONOMIC_GENUS_ONLY_PLACEHOLDERS,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_HYBRID_CONNECTORS as TAXONOMIC_HYBRID_CONNECTORS,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_INFRASPECIFIC_RANK_ALIASES as TAXONOMIC_INFRASPECIFIC_RANK_ALIASES,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_INFRASPECIFIC_RANKS as TAXONOMIC_INFRASPECIFIC_RANKS,
)
from nwkit.taxonomic_labels import (
    TAXONOMIC_PROXIMITY_QUALIFIERS as TAXONOMIC_PROXIMITY_QUALIFIERS,
)
from nwkit.taxonomic_labels import (
    base_species_label as base_species_label,
)
from nwkit.taxonomic_labels import (
    canonical_taxonomic_token as canonical_taxonomic_token,
)
from nwkit.taxonomic_labels import (
    display_rank_token as display_rank_token,
)
from nwkit.taxonomic_labels import (
    extract_species_label as extract_species_label,
)
from nwkit.taxonomic_labels import (
    is_hybrid_binomial_connector as is_hybrid_binomial_connector,
)
from nwkit.taxonomic_labels import (
    is_hybrid_connector_token as is_hybrid_connector_token,
)
from nwkit.taxonomic_labels import (
    is_species_rank_or_qualifier_token as is_species_rank_or_qualifier_token,
)
from nwkit.taxonomic_labels import (
    matches_species_label as matches_species_label,
)
from nwkit.taxonomic_labels import (
    normalize_taxonomic_name_text as normalize_taxonomic_name_text,
)
from nwkit.taxonomic_labels import (
    scientific_name_from_label as scientific_name_from_label,
)
from nwkit.taxonomic_labels import (
    species_key_token as species_key_token,
)
from nwkit.taxonomic_labels import (
    species_label_from_taxonomic_text as species_label_from_taxonomic_text,
)
from nwkit.taxonomic_labels import (
    species_label_prefix_token as species_label_prefix_token,
)
from nwkit.taxonomic_labels import (
    species_prefix_token_count as species_prefix_token_count,
)
from nwkit.taxonomic_labels import (
    strip_species_label as strip_species_label,
)
from nwkit.taxonomic_labels import (
    strip_species_label_terminal_suffixes as strip_species_label_terminal_suffixes,
)
from nwkit.taxonomic_labels import (
    taxonomic_token_key as taxonomic_token_key,
)
from nwkit.taxonomic_labels import (
    tokenize_taxonomic_name as tokenize_taxonomic_name,
)


def assembly_only_species(assembly_dir, gff_dir, genome_dir):
    """Exclude known assemblies from genomic rescue only when both genomic inputs are absent."""
    from pathlib import Path

    if not assembly_dir:
        return set()
    assemblies = Path(assembly_dir) / 'assembled_transcripts_with_isoforms'
    names = {path.name.removesuffix('.gz').removesuffix('_isoform.fa')
             for path in assemblies.iterdir() if path.is_file()
             and path.name.removesuffix('.gz').endswith('_isoform.fa')} if assemblies.is_dir() else set()
    for directory, suffixes in ((gff_dir, ('.gff', '.gff3')), (genome_dir, ('.fa', '.fas', '.fasta', '.fna'))):
        if directory and Path(directory).is_dir():
            names.difference_update(extract_species_label(path.name) for path in Path(directory).iterdir()
                                    if path.is_file() and path.name.removesuffix('.gz').endswith(suffixes))
    return names
