"""Delimited assembly identifiers and coherent NCBI annotation bundles."""

import re

NCBI_ACCESSION_SEARCH = re.compile(
    r"(?<![A-Za-z0-9])(GC[AF]_[0-9]+\.[0-9]+)(?![A-Za-z0-9]|\.[0-9])", re.IGNORECASE)
GWH_ACCESSION_SEARCH = re.compile(
    r"(?<![A-Za-z0-9])(GWH[A-Z0-9]+(?:\.[0-9]+)?)(?![A-Za-z0-9]|\.[0-9])", re.IGNORECASE)


def gwh_accession_matches(text, accession):
    tokens = {token.upper() for token in GWH_ACCESSION_SEARCH.findall(str(text))}
    requested = str(accession).upper()
    if "." in requested:
        return tokens == {requested}
    return len(tokens) == 1 and next(iter(tokens)).split(".", 1)[0] == requested


def select_ncbi_bundle(groups):
    """Index each file once, keeping accession versions distinct.

    Preserve the established completeness score and accession tie-break. An
    untagged directory retains legacy discovery; tagged files cannot recruit
    another assembly version or an untagged companion into their bundle.
    """
    bundles = {}
    for index, paths in enumerate(groups):
        for path in paths:
            accessions = {token.upper() for token in NCBI_ACCESSION_SEARCH.findall(path.name)}
            if len(accessions) > 1:
                raise ValueError("multiple assembly accessions in filename: " + path.name)
            if accessions:
                accession = accessions.pop()
                bundles.setdefault(accession, [[], [], [], []])[index].append(path)
    if not bundles:
        return "", groups, ()
    coherent = []
    for accession, bundle in bundles.items():
        cds, gff, gbff, genome = bundle
        if gbff or gff and (cds or genome):
            coherent.append((-(bool(cds) + bool(gff) + bool(genome)), accession))
    if not coherent and len(bundles) > 1:
        raise ValueError("annotation files span incompatible assembly accessions: " + ",".join(sorted(bundles)))
    selected = min(coherent)[1] if coherent else next(iter(bundles))
    return selected, bundles[selected], tuple(sorted(bundles))
