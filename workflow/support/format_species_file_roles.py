"""Filename role markers shared by discovery, manifests and provider selection."""

import re
from pathlib import Path

CDS_MARKER = r"(?:annotated)?cdss?"
TRANSCRIPT_MARKER = r"(?:transcripts?|mrnas?|cdnas?)"
PROTEIN_MARKER = r"(?:(?:annotated)?proteins?|peps?|peptides?|polypeptides?|translated[._-]cds)"
NONCODING_MARKER = r"(?:(?:t|r|nc)?rnas?|rna[._-]from[._-]genomic)"
REPEAT_MARKER = r"(?:repeats?|repeatmask(?:er)?|te|teanno|transpos(?:on)?)"
ASSEMBLY_MARKER = (
    r"(?:genomes?|genomic|assembly|dna|toplevel|pseudomolecules?|pseudomol|pmol|"
    r"(?:chr|chrom|chromosome)(?:[0-9][a-z0-9]*|[xymcp]|all)?|"
    r"scaffolds?(?:[0-9][a-z0-9]*)?|asm[0-9][a-z0-9]*|"
    r"(?:soft|hard|repeat)?masked)"
)
HAPLOTYPE_MARKER = r"(?:hap(?:[0-9][0-9a-z]*)?|haplotype[0-9a-z]*)"


def filename_role_text(name, species_key=""):
    """Exclude directories and an explicitly known species prefix from roles."""
    text = Path(str(name or "")).name.lower()
    species = str(species_key or "").strip().lower()
    if species:
        prefix = r"[._ -]+".join(re.escape(token) for token in re.split(r"[_\s]+", species))
        text = re.sub(r"^" + prefix + r"(?=[._-]|$)", "", text, count=1).lstrip("._-")
    return text


def has_filename_marker(name, pattern, species_key=""):
    return re.search(
        r"(?:^|[._-])(?:" + pattern + r")(?=[._-]|$)",
        filename_role_text(name, species_key),
    ) is not None


def has_non_cds_role(name, species_key=""):
    """Keep protein, noncoding RNA and repeat FASTAs out of CDS fallbacks too."""
    return has_filename_marker(
        name, "|".join((PROTEIN_MARKER, NONCODING_MARKER, REPEAT_MARKER)), species_key
    )
