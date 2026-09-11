# IQ-TREE zero-frequency ancestral-state failure

Observed with IQ-TREE 3.1.4 in the ARM64 runtime described in the parent record.
The code-2 alignment has no internal stops, but lacks TGG observations.
Run from this directory with a writable working directory:

```bash
iqtree -s csubst.fasta -te csubst.nwk -m GY+F+R4 -T AUTO \
  --threads-max 2 --seqtype CODON2 --prefix reproduction \
  --ancestral --rate --seed 12345 --redo
```

The observed run exits successfully, logs numerical underflow, and writes NaN
ancestral probabilities. `state-excerpt.tsv` retains the first affected row;
`csubst.log` is the original log. This is an unresolved upstream issue, not
a GeneGalleon workaround or a passing validation case.
