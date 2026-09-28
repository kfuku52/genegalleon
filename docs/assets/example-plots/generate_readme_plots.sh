#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
support_dir="${repo_root}/workflow/support"
family_dir="${repo_root}/workspace/output/query2family"
species_tree="${family_dir}/parameters/undated_species_tree.pruned.nwk"

if [[ "${1:-}" == "--in-runtime" ]]; then
  stage=$2
  cd "${stage}"
  mkdir -p Rlib
  R CMD INSTALL --library="${stage}/Rlib" "${support_dir}/treevis"
  export R_LIBS_USER="${stage}/Rlib"

  gzip -cd "${family_dir}/mafft/AHA_cds.aln.fa.gz" > AHA_untrimmed.fa
  gzip -cd "${family_dir}/clipkit/AHA_cds.clipkit.fa.gz" > AHA_trimmed.fa

  # AHA query proteins exceed 300 aa, so the workflow's auto E-value is 0.01.
  python "${support_dir}/synteny_neighbors.py" \
    --focal_cds_fasta "${family_dir}/cds_fasta/AHA_cds.fa.gz" \
    --dir_sp_cds "${repo_root}/workspace/input/species_cds" \
    --dir_sp_gff "${repo_root}/workspace/input/species_gff" \
    --cache_dir "${stage}/synteny_cache" \
    --lock_dir "${stage}/synteny_locks" \
    --gff2genestat_script "${support_dir}/gff2genestat.py" \
    --input_sequence_mode cds --window 20 --evalue 0.01 \
    --genetic_code 1 --threads 2 --outfile AHA_synteny.tsv

  python "${support_dir}/gff2genestat.py" \
    --validate-cds-length --structure-policy report --phase-policy report \
    --dir_gff "${repo_root}/workspace/input/species_gff" \
    --feature CDS --multiple_hits longest \
    --seqfile "${family_dir}/cds_fasta/AHA_cds.fa.gz" \
    --ncpu 2 --outfile AHA_gff.tsv

  python "${support_dir}/orthogroup_statistics.py" \
    --species_tree "${species_tree}" \
    --untrimmed_aln AHA_untrimmed.fa --trimmed_aln AHA_trimmed.fa \
    --unrooted_tree "${family_dir}/iqtree_tree/AHA_iqtree.nwk" \
    --rooted_tree "${family_dir}/rooted_tree/AHA_root.nwk" \
    --rooting_log "${family_dir}/rooted_tree_log/AHA_root.txt" \
    --expression "${family_dir}/character_expression/AHA_expression.tsv" \
    --character_gff AHA_gff.tsv \
    --rpsblast "${family_dir}/rpsblast/AHA_rpsblast.tsv" \
    --synteny AHA_synteny.tsv \
    --clade_ortholog_prefix Arabidopsis_thaliana_ --ncpu 2

  python "${support_dir}/annotate_stat_branch_query_markers.py" \
    --stat_branch orthogroup.branch.tsv \
    --query_gene "${repo_root}/workspace/input/query_gene/AHA" \
    --query_aa_fasta "${family_dir}/query_aa_fasta/AHA_query.aa.fa.gz" \
    --query_blast "${family_dir}/query_blast/AHA_query_blast.tsv" \
    --outfile AHA_stat.branch.tsv --min_query_blast_coverage 0.25

  python - <<'PY'
import csv
from pathlib import Path

source = Path("AHA_stat.branch.tsv")
with source.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = list(reader.fieldnames or [])
    rows = list(reader)
fields.append("dup_conf_label")
for row in rows:
    score = row.get("dup_conf_score", "")
    row["dup_conf_label"] = f"{float(score):.2f}" if row.get("so_event") == "D" and score else ""
with source.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY

  Rscript "${support_dir}/stat_branch2tree_plot.r" \
    --stat_branch=AHA_stat.branch.tsv \
    --panel_widths_mm=tree:60 \
    --panel1=tree,bl_rooted,dup_conf_label,species,L \
    --panel2=heatmap,no,abs,_,expression_ \
    --panel3=cluster_membership,100000 \
    --panel4=synteny,AHA_synteny.tsv,5 \
    --panel5=tiplabel \
    --panel6=categorical,query_marker,Query,- \
    --panel7=intron_number \
    --panel8=gene_structure,compressed,23,AHA_untrimmed.fa \
    --panel9="domain,${family_dir}/rpsblast/AHA_rpsblast.tsv" \
    --panel10=alignment,AHA_trimmed.fa,AHA_untrimmed.fa \
    --panel11=ortholog,Arabidopsis_thaliana_,no \
    --panel12=synteny_similarity,AHA_synteny.tsv,20,15 \
    --panel13=sequence_similarity,AHA_trimmed.fa,cds,15,1 \
    --show_branch_id=no --event_method=species_overlap \
    --species_color_table=PLACEHOLDER \
    --pie_chart_value_transformation=identity --long_branch_display=auto

  python "${support_dir}/gene_family_presence_absence.py" \
    --mode query2family --dir_gene_family "${family_dir}" \
    --dir_query_gene "${repo_root}/workspace/input/query_gene" \
    --species_tree "${species_tree}" \
    --out_presence query2family_presence_absence.tsv \
    --out_copy_number query2family_copy_number.tsv \
    --out_long query2family_presence_absence.long.tsv \
    --out_plot_presence query2family_presence_absence.plot.tsv \
    --out_plot_copy_number query2family_copy_number.plot.tsv \
    --out_plot_long query2family_presence_absence.plot.long.tsv \
    --out_selection query2family_presence_absence.plot_selection.tsv \
    --include_incomplete 0 --max_families all

  # Overlay the newly calculated AHA neighborhoods without changing curated
  # query2family outputs. The ortholog summarizer reads physical store files.
  mkdir -p staged_family/{cds_fasta,stat_branch,synteny}
  for family_id in AHA STRICTCHK YABBY; do
    ln "${family_dir}/cds_fasta/${family_id}_cds.fa.gz" "staged_family/cds_fasta/${family_id}_cds.fa.gz"
    ln "${family_dir}/stat_branch/${family_id}_stat.branch.tsv" "staged_family/stat_branch/${family_id}_stat.branch.tsv"
    if [[ "${family_id}" == AHA ]]; then
      ln AHA_synteny.tsv staged_family/synteny/AHA_synteny.tsv
    else
      ln "${family_dir}/synteny/${family_id}_synteny.tsv" "staged_family/synteny/${family_id}_synteny.tsv"
    fi
  done

  python "${support_dir}/query_gene_orthologs.py" \
    --basis reference_species \
    --dir_gene_family staged_family \
    --dir_query_gene "${repo_root}/workspace/input/query_gene" \
    --family_file query2family_presence_absence.plot_selection.tsv \
    --reference_species Arabidopsis_thaliana \
    --out_columns query2family_reference_gene_orthologs.columns.tsv \
    --out_glyphs query2family_reference_gene_orthologs.glyphs.tsv \
    --out_tree query2family_reference_gene_orthologs.tree.tsv \
    --out_synteny query2family_reference_gene_orthologs.synteny.tsv \
    --out_ufboot query2family_reference_gene_orthologs.ufboot.tsv

  Rscript "${support_dir}/plot_query2family_presence_absence.R" \
    "--species_tree=${species_tree}" \
    --long_table=query2family_presence_absence.plot.long.tsv \
    --ortholog_column_table=query2family_reference_gene_orthologs.columns.tsv \
    --ortholog_glyph_table=query2family_reference_gene_orthologs.glyphs.tsv \
    --ortholog_tree_table=query2family_reference_gene_orthologs.tree.tsv \
    --ortholog_synteny_table=query2family_reference_gene_orthologs.synteny.tsv \
    --ortholog_ufboot_table=query2family_reference_gene_orthologs.ufboot.tsv \
    "--species_mapping_tree=${species_tree}" \
    --ortholog_basis=reference_species \
    --reference_species=Arabidopsis_thaliana \
    --evidence_layout=band --value=presence --width=7.2 \
    --out_pdf=query2family_reference_gene_orthologs.pdf

  exit 0
fi

cd "${repo_root}"
stage="$(mktemp -d "${repo_root}/workspace/output/readme-plots.XXXXXX")"
trap 'rm -rf -- "${stage}"' EXIT
bash workflow/tests/run_in_runtime.sh bash "${script_dir}/generate_readme_plots.sh" --in-runtime "${stage}"

pdftoppm_bin=""
IFS=: read -r -a path_entries <<< "${PATH}"
for path_entry in "${path_entries[@]}"; do
  candidate="${path_entry}/pdftoppm"
  if [[ -x "${candidate}" ]] && "${candidate}" -v >/dev/null 2>&1; then
    pdftoppm_bin="${candidate}"
    break
  fi
done
if [[ -z "${pdftoppm_bin}" ]]; then
  echo "A working pdftoppm executable is required to render the README PNGs." >&2
  exit 1
fi
"${pdftoppm_bin}" -f 1 -singlefile -r 180 -png \
  "${stage}/stat_branch2tree_plot.pdf" "${script_dir}/readme-aha-tree-plot"
"${pdftoppm_bin}" -f 1 -singlefile -r 240 -png \
  "${stage}/query2family_reference_gene_orthologs.pdf" "${script_dir}/readme-test-family-presence-absence"
