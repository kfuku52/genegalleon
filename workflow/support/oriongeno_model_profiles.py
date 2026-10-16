"""Native checkpoint directory names published by OrionGeno."""

PROFILES = frozenset({'plants', 'mammals', 'birds', 'fish', 'other_vertebrates',
                      'arthropods', 'other_invertebrates', 'fungi'})

# Match the published checkpoint groups. The upstream benchmark identifies the
# fish checkpoint as Actinopterygii; other vertebrates use the broader model.
ANCHORS = {'plants': (33090, 'Viridiplantae'), 'mammals': (40674, 'Mammalia'),
           'birds': (8782, 'Aves'), 'fish': (7898, 'Actinopterygii'),
           'other_vertebrates': (7742, 'Vertebrata'), 'arthropods': (6656, 'Arthropoda'),
           'other_invertebrates': (33208, 'Metazoa'), 'fungi': (4751, 'Fungi')}
ANCESTOR_PROFILES = {'mammals': {'other_vertebrates', 'other_invertebrates'},
                     'birds': {'other_vertebrates', 'other_invertebrates'},
                     'fish': {'other_vertebrates', 'other_invertebrates'},
                     'other_vertebrates': {'other_invertebrates'},
                     'arthropods': {'other_invertebrates'}}


def profile_for_lineage(lineage):
    """Choose the most specific checkpoint, excluding its broader remainder groups."""
    by_taxid = {taxid: (profile, name) for profile, (taxid, name) in ANCHORS.items()}
    matches = []
    for node in lineage:  # The verified NCBI track is ordered root to species.
        if node['taxid'] in by_taxid:
            profile, name = by_taxid[node['taxid']]
            if node['spname'] != name:
                raise ValueError('NCBI model anchor name contradicts its TaxID')
            matches.append(profile)
    if not matches:
        raise ValueError('No supported OrionGeno model in NCBI ancestors')
    profile = matches[-1]
    if set(matches[:-1]) != ANCESTOR_PROFILES.get(profile, set()):
        raise ValueError('Conflicting OrionGeno clades in NCBI ancestors')
    return profile, ANCHORS[profile][0]
