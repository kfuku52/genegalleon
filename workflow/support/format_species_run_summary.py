"""Final statistics, reporting, and status for species formatting runs."""

import json
from pathlib import Path


def finalize_run(
    *,
    args,
    stats,
    output_cds_dir,
    output_gff_dir,
    output_genome_dir,
    species_summary_path,
    error_count,
    failed_format_tasks,
):
    if args.stats_output != "":
        stats_path = Path(args.stats_output).expanduser().resolve()
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "wt", encoding="utf-8") as handle:
            json.dump(stats, handle, ensure_ascii=True, indent=2, sort_keys=True)

    print(
        "Finished. species processed: {}, CDS aggregated away: {}, CDS before/after: {}/{}, first CDS sequence: {}, output CDS dir: {}, output GFF dir: {}, output genome dir: {}, species summary: {}, dry_run={}".format(
            stats["species_processed"],
            stats["aggregated_cds_removed"],
            stats["cds_sequences_before"],
            stats["cds_sequences_after"],
            stats["cds_first_sequence_name"] if stats["cds_first_sequence_name"] != "" else "NA",
            output_cds_dir,
            output_gff_dir,
            output_genome_dir,
            species_summary_path,
            int(args.dry_run),
        )
    )
    if error_count > 0 or failed_format_tasks > 0:
        return 2
    return 0
