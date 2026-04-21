# Gene-Family Outputs and Progress Monitoring

This page expands the short stage references for
`gg_gene_evolution_entrypoint.sh`, especially for users familiar with the old
`gfe_geneFamilyPhylogeny` workflow.

## Query2family input recap

In `mode_gene_evolution=query2family`:

- place one file per family under `workspace/input/query_gene/`,
- the file basename becomes the family/task ID,
- the scheduler array size is the number of files in that directory.

Each input file should contain query genes for one gene family only. If genes
from different families are combined in one file, GeneGalleon treats them as one
family-level task, which can artificially merge unrelated family phylogenies into
one unnatural tree and PDF output.

Accepted query-file forms are:

- plain gene ID list,
- in-frame CDS FASTA,
- protein FASTA.

See [Input Conventions](input-conventions.md) for concrete examples.

## Main per-family outputs

If the input file is `workspace/input/query_gene/2_WOX`, the main outputs are
typically written under `workspace/output/query2family/` as:

- `query_gene/2_WOX`
- `query_aa_fasta/2_WOX_query.aa.fa.gz`
- `cds_fasta/2_WOX_cds.fa.gz`
- `mafft/`, `clipkit/`, and/or `trimal/` alignment products depending on the enabled cleaning steps
- `iqtree_tree/2_WOX_iqtree.nwk`
- `rooted_tree/2_WOX_root.nwk`
- `dated_tree/2_WOX_dated.nwk` and `dated_tree_log/2_WOX_dated.log.txt` when `run_tree_dating=1`
- `stat_branch/2_WOX_stat.branch.tsv`
- `stat_tree/2_WOX_stat.tree.tsv`
- `tree_plot/2_WOX_tree_plot.pdf`

In `mode_gene_evolution=orthogroup`, the same style of per-family outputs is
written under `workspace/output/orthogroup/`.

## How to read `tree_plot`

`tree_plot` is generated from `stat_branch/*.tsv` and summarizes many
gene-family attributes around the inferred tree. In the current default
GeneGalleon configuration, the panel order is:

1. gene tree
2. expression heatmap
3. expression point plot
4. cluster-membership panel
5. local synteny panel
6. tip labels
7. signal peptide summary
8. transmembrane-domain summary
9. intron-count panel
10. protein-domain panel from RPS-BLAST
11. alignment panel
12. promoter-motif panel from FIMO
13. MEME motif summary
14. ortholog-context panel

Practical interpretation:

- the leftmost tree is the anchor and carries the branch length, support, and
  branch-color settings chosen by `treevis_*` parameters,
- expression, cluster, and synteny panels help interpret lineage-specific
  expansions or expression shifts,
- structural panels summarize domains, introns, signal peptides, and membrane
  predictions on the same row order as the tree tips,
- if an upstream analysis was disabled, or the corresponding inputs were not
  available, the associated panel may be blank or minimally populated.

`stat.branch.tsv` is the master table that collects per-branch and per-tip
annotations for plotting. `stat.tree.tsv` is the paired tree-level summary.

## Orthogroup-scale progress summaries

GeneGalleon provides a built-in summary wrapper for orthogroup and
transcriptome-scale runs:

```bash
cd workflow
bash gg_progress_summary_entrypoint.sh
```

This writes:

- `workspace/orthogroup_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

For orthogroup runs, `orthogroup_summary.tsv` is useful because:

- it adds `GG_ARRAY_TASK_ID`, which is the row index to resubmit,
- it appends AMAS-derived alignment statistics such as
  `Parsimony_informative_sites_clean`,
- it adds one `1/0` completion column per visible output subdirectory under
  `workspace/output/orthogroup/`.

Practical reading tips:

- rows with `Parsimony_informative_sites_clean == 0` cannot produce normal
  IQ-TREE-based downstream outputs,
- rows with `0` in late-stage columns such as `stat_tree`, `stat_branch`, or
  `tree_plot` are the first candidates to inspect or rerun.

Transcriptome summary rows are species based rather than orthogroup based.
The wrapper also prints incomplete species IDs inferred from missing
`amalgkit_getfastq/<species>/*.safely_removed` markers.

## Query2family completion audit

There is currently no dedicated `query2family_summary.tsv`.
For large query2family runs, compare the input basenames with a downstream
completion marker such as `tree_plot` or `stat_branch`.

Example using `tree_plot` as the completion marker:

```bash
python - <<'PY'
from pathlib import Path

query_dir = Path("workspace/input/query_gene")
plot_dir = Path("workspace/output/query2family/tree_plot")

queries = sorted(p.name for p in query_dir.iterdir() if p.is_file())
finished = {p.name[:-len("_tree_plot.pdf")] for p in plot_dir.glob("*_tree_plot.pdf")}

for task_id, family in enumerate(queries, start=1):
    if family not in finished:
        print(f"{task_id}\t{family}")
PY
```

If `run_tree_plot=0`, switch the marker to `stat_branch/*_stat.branch.tsv`
instead of `tree_plot/*_tree_plot.pdf`.

## Rerunning incomplete tasks

Once you know the task IDs:

Slurm:

```bash
sbatch --array=17,42,105 workflow/gg_gene_evolution_entrypoint.sh
```

UGE:

```bash
qsub -t 17,42,105 workflow/gg_gene_evolution_entrypoint.sh
```

Local rerun of one task:

```bash
GG_ARRAY_TASK_ID=17 GG_TASK_CPUS=4 bash workflow/gg_gene_evolution_entrypoint.sh
```

If your site prefers one submission per task rather than an explicit list-style
array, submit a short loop around the same wrapper.
