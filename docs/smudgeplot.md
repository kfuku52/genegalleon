# Read k-mer structure with Smudgeplot

Set `run_smudgeplot=1` in `workflow/gg_genome_annotation_entrypoint.sh` to run
[FastK](https://github.com/thegenemyers/FASTK) and
[Smudgeplot](https://github.com/KamilSJaron/smudgeplot) on each species' FASTQ
files in `workspace/input/species_dnaseq/<species>/`. Plain and gzip-compressed
FASTQ and symlinked inputs are supported. Other annotation stages remain
independently configurable. Rebuild the runtime when updating from an image
without these tools.

Discovery ignores hidden files and directories, resolves symlinks, and rejects
empty files, duplicate files (including hardlinks), broken FASTQ links and
repeated or cyclic directories. Only the selected FASTQ files are hashed for
provenance. Unrelated files do not invalidate a completed analysis.

Canonical input, output and database paths must contain no whitespace or shell
metacharacters. The current upstream native programs interpolate these paths
into shell commands; GeneGalleon rejects unsafe paths before counting.

The stage publishes
`workspace/output/species_dnaseq_smudgeplot/<species>_smudgeplot.zip`, containing
linear and log-density PDF plots, the FastK histogram, the k-mer pair table,
Smudgeplot's inference report, commands and input metadata. The large FastK
database is temporary. Input-aware provenance governs reruns and atomic ZIP
publication prevents incomplete archives from replacing prior results.

Parameters:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `smudgeplot_kmer_length` | 21 | k-mer length |
| `smudgeplot_lower_count` | `auto` | Smudgeplot histogram cutoff, or an explicit positive count |
| `smudgeplot_aggregation_distance` | 2 | Smudgeplot's local aggregation neighbourhood |

Inspect the histogram together with GenomeScope before interpreting the plot.
Smudgeplot's automatic cutoff is at least 10, which can remove genuine
heterozygous k-mers in low-coverage datasets. Set the cutoff explicitly when
needed and assess sensitivity to nearby values. The inferred haplotype structure
supports ploidy assessment; limited coverage, repeats and divergent haplotypes
can prevent a decisive inference.

For a standalone run in the GeneGalleon runtime:

```bash
python workflow/support/run_smudgeplot.py \
  --reads /data/species.hifi.fastq.gz \
  --output-dir /results/species.smudgeplot \
  --kmer-length 21 --lower-count 4 --threads 4 --memory-gb 12
```

The output directory must be empty. Add `--database-dir /scratch/species.fastk`
to retain the database for cutoff sensitivity analyses; this directory must also
be empty. The example cutoff is illustrative and must be chosen for the data.

Smudgeplot runs in an isolated Python environment because its weighted
percentile computation requires NumPy >=2.0. This preserves the main runtime's
NumPy constraints for other scientific tools. Exact upstream revisions are
recorded in `/opt/pg/logs/source_revisions.tsv`.

### Upstream issues found during integration

Smudgeplot's `src/lib/PloidyPlot.c` passes the database and temporary paths to
`Logex`, `Symmex` and `Fastrm` through unquoted shell commands. A FastK database
under `database with spaces/` is counted successfully but `smudgeplot hetmers`
fails in `Symmex`. Fixing path handling belongs in the native dependency;
GeneGalleon's preflight restriction can be removed once that behaviour is fixed
and tested upstream.

The same source uses `mktemp("._SPAIR.XXXX")`, with four trailing `X` characters
instead of the six required by `mktemp`. On Linux this produces an empty root,
so conditioning files are named `.trim.ktab` and `.symx.ktab` in the working
directory, and failures can leave them behind. This upstream defect remains
unpatched here. Inspect failed output directories before retaining or sharing
them; completed stage archives are published only after successful analysis.

An AddressSanitizer build (`make CFLAGS='-O1 -g -fsanitize=address
-fno-omit-frame-pointer'`) also detects a heap-buffer-overflow in the bundled
`src/lib/libfastk.c:Current_Entry`: it allocates `S->pbyte` but writes the prefix
and suffix together. The current FastK library allocates `S->tbyte` for this
entry. Smudgeplot must update its bundled library upstream; GeneGalleon does not
overwrite vendored dependency sources during the build.

In a synthetic diploid dataset with about 12-fold haploid coverage, automatic
cutoff 10 caused `smudgeplot all` to fail during smudge classification with a
pandas `.str` accessor error. The same data worked at explicit cutoff 4.
GeneGalleon preserves the histogram and command log and reports the failure;
it does not lower the cutoff automatically. The runtime checks separately cover
manual cutoff at about 12-fold coverage and automatic cutoff at about 24-fold
coverage.
