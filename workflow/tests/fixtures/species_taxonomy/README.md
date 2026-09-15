# Species taxonomy fixture

`species.tsv` lists eight real plant species without supplying TaxIDs, so the
tests and example exercise automatic name resolution.

`lineages.json` is a small, fixed subset extracted from GeneGalleon's existing
NCBI/ETE cache on 2026-09-15. Tests build a temporary SQLite database from this
file and never download taxonomy data. Production does not use this fixture.

`species_tree.nwk` is an illustrative topology with artificial branch lengths
and numeric support labels for preservation checks. It is not a phylogenetic
estimate or a source of biological evidence.

See `docs/species-taxonomy.md` for runnable plotting examples using the actual
shared taxonomy database.
