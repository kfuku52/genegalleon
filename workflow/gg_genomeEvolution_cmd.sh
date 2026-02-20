#!/usr/bin/env bash

#run a busco analysis for all species in species_cds and extract shared complete BUSCO genes among all species.
#then make cds fasta files including all species' cds sequenes for each BUSCO genes and make alignment
#then concatenate all fasta files for further downstream analysis -> task = trimal in gg_species_tree.

### Start: Modify this block to tailor your analysis ###

# Workflow
# Values should be boolean (0 or 1).
run_busco=1
run_get_busco_summary=1
# Polyploidization history estimation
run_busco_getfasta=1
run_busco_mafft=1
run_busco_trimal=1
run_busco_iqtree_dna=1
run_busco_iqtree_pep=1
run_busco_notung_root_dna=1
run_busco_notung_root_pep=1
run_busco_root_dna=1
run_busco_root_pep=1
run_busco_grampa_dna=1
run_busco_grampa_pep=1
run_orthogroup_grampa=1 # Requires outputs of run_tree_root=1 of gg_geneFamilyPhylogeny with mode_orthogroup=1
# Gene family evolution
run_cafe=0
run_go_enrichment=0

# Other parameters
genetic_code=1 # Integer. See here https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
busco_lineage="embryophyta_odb12" # See here for available datasets: https://busco-data.ezlab.org/v5/data/lineages/
nucleotide_model='GTR+R4' # IQ-TREE's model
protein_model='LG+R4' # IQ-TREE's model
grampa_h1='' # '' for exploratory analysis. If an internal node should be specified, run grampa once with grampa_h1='' to obtain internal node IDs.
min_gene_orthogroup_grampa=5 # Lower threshold for the orthogroup size to be included in the orthogroup-based grampa analysis.
max_gene_orthogroup_grampa=50 # Upper threshold for the orthogroup size to be included in the orthogroup-based grampa analysis.
notung_jar="/usr/local/bin/Notung.jar"
max_size_differential_cafe=9999999 # Maximum allowed size differential among species for each orthogroup to be included in the CAFE5 analysis.
n_gamma_cats_cafe=4 # Number of gamma rate categories to use in CAFE5.
target_branch_go="<25>" # Target branch for GO enrichment analysis (species name or branch ID; see workspace/genome_evolution/cafe/summary_plot/branch_id.pdf)
change_direction_go="increase" # Direction of gene change to focus on ("increase" or "decrease") in GO enrichment analysis.
species_go_annotation="Arabidopsis_thaliana" # Species used for GO annotation
go_category="BP,MF,CC" # GO categories to analyze (BP: Biological Process, MF: Molecular Function, CC: Cellular Component)

### End: Modify this block to tailor your analysis ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

copy_busco_tables() {
  local busco_root_dir="$1"
  local lineage="$2"
  local file_full="$3"
  local file_short="$4"
  local run_dir="${busco_root_dir}/run_${lineage}"
  local full_src=""
  local short_src=""
  local -a short_candidates=()

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_src="${run_dir}/full_table.tsv"
  elif [[ -s "${run_dir}/full_table.tsv.gz" ]]; then
    gzip -cd "${run_dir}/full_table.tsv.gz" > "./tmp.busco.full_table.tsv"
    full_src="./tmp.busco.full_table.tsv"
  fi

  if [[ -s "${run_dir}/short_summary.txt" ]]; then
    short_src="${run_dir}/short_summary.txt"
  else
    mapfile -t short_candidates < <(find "${run_dir}" -maxdepth 1 -type f -name "short_summary*.txt" | sort)
    if [[ ${#short_candidates[@]} -gt 0 ]]; then
      short_src="${short_candidates[0]}"
    fi
  fi

  if [[ -z "${full_src}" || -z "${short_src}" ]]; then
    echo "BUSCO outputs were not found under ${run_dir}."
    echo "Expected full_table.tsv(.gz) and short_summary*.txt."
    rm -f "./tmp.busco.full_table.tsv"
    return 1
  fi

  cp_out "${full_src}" "${file_full}"
  cp_out "${short_src}" "${file_short}"
  rm -f "./tmp.busco.full_table.tsv"
}

dir_sp_cds="${dir_pg_input}/species_cds"
dir_og="${dir_pg_output}/orthogroup"
dir_og_rooted_tree="${dir_og}/rooted_tree"
dir_species_tree="${dir_pg_output}/species_tree"
dir_genome_evolution="${dir_pg_output}/genome_evolution"
dir_species_busco_full="${dir_pg_output}/species_cds_busco_full"
dir_species_busco_short="${dir_pg_output}/species_cds_busco_short"
dir_busco_summary_table="${dir_genome_evolution}/busco_summary_table"
dir_busco_fasta="${dir_genome_evolution}/busco_cds_fasta"
dir_busco_mafft="${dir_genome_evolution}/busco_mafft"
dir_busco_trimal="${dir_genome_evolution}/busco_trimal"
dir_busco_iqtree_dna="${dir_genome_evolution}/busco_iqtree_dna"
dir_busco_notung_dna="${dir_genome_evolution}/busco_notung_dna"
dir_busco_rooted_nwk_dna="${dir_genome_evolution}/busco_rooted_nwk_dna"
dir_busco_rooted_txt_dna="${dir_genome_evolution}/busco_rooted_txt_dna"
dir_busco_grampa_dna="${dir_genome_evolution}/grampa_busco_dna"
dir_busco_iqtree_pep="${dir_genome_evolution}/busco_iqtree_pep"
dir_busco_notung_pep="${dir_genome_evolution}/busco_notung_pep"
dir_busco_rooted_nwk_pep="${dir_genome_evolution}/busco_rooted_nwk_pep"
dir_busco_rooted_txt_pep="${dir_genome_evolution}/busco_rooted_txt_pep"
dir_busco_grampa_pep="${dir_genome_evolution}/grampa_busco_pep"
dir_orthogroup_grampa="${dir_genome_evolution}/grampa_orthogroup"
dir_cafe="${dir_genome_evolution}/cafe"
dir_cafe_orthogroup_selection="${dir_cafe}/orthogroup_selection"
dir_cafe_output="${dir_cafe}/cafe_output"
dir_cafe_summary_plot="${dir_cafe}/summary_plot"
dir_cafe_each_family_plot="${dir_cafe}/each_family_plot"
dir_species_cds_annotation="${dir_pg_output}/species_cds_annotation"
dir_go_enrichment="${dir_cafe}/go_enrichment"
dir_tmp="${dir_genome_evolution}/tmp"

check_species_cds ${dir_pg}
check_if_species_files_unique ${dir_sp_cds}
memory_notung=${MEM_PER_SLOT}

file_genecount="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv"
file_dated_species_tree="${dir_species_tree}/dated_species_tree.nwk"
file_busco_summary_table="${dir_busco_summary_table}/busco_summary.tsv"
file_busco_grampa_dna="${dir_busco_grampa_dna}/grampa_summary.tsv"
file_busco_grampa_pep="${dir_busco_grampa_pep}/grampa_summary.tsv"
file_orthogroup_grampa="${dir_orthogroup_grampa}/grampa_summary.tsv"
file_gene_id="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.selected.tsv"
file_go_annotation="${dir_species_cds_annotation}/${species_go_annotation}_annotation.tsv"

if [[ ${run_busco_notung_root_dna} -eq 1 || ${run_busco_notung_root_pep} -eq 1 || ${run_busco_root_dna} -eq 1 || ${run_busco_root_pep} -eq 1 || ${run_busco_grampa_dna} -eq 1 || ${run_busco_grampa_pep} -eq 1 ]]; then
  if [[ ! -s "${notung_jar}" ]]; then
    echo "Notung jar was not found: ${notung_jar}"
    echo "Disabling NOTUNG-dependent BUSCO tasks: run_busco_notung_root_dna, run_busco_notung_root_pep, run_busco_root_dna, run_busco_root_pep, run_busco_grampa_dna, run_busco_grampa_pep"
    run_busco_notung_root_dna=0
    run_busco_notung_root_pep=0
    run_busco_root_dna=0
    run_busco_root_pep=0
    run_busco_grampa_dna=0
    run_busco_grampa_pep=0
  fi
fi

if [[ ${run_orthogroup_grampa} -eq 1 ]]; then
  if [[ ! -s "${file_genecount}" ]]; then
    echo "Disabling run_orthogroup_grampa because required file is missing: ${file_genecount}"
    run_orthogroup_grampa=0
  elif [[ ! -d "${dir_og_rooted_tree}" ]]; then
    echo "Disabling run_orthogroup_grampa because required directory is missing: ${dir_og_rooted_tree}"
    run_orthogroup_grampa=0
  fi
fi

ensure_dir "${dir_tmp}"
cd ${dir_tmp}

enable_all_run_flags_for_debug_mode

task="BUSCO analysis of species-wise CDS files"
if [[ ${run_busco} -eq 1 ]]; then
  ensure_dir "${dir_species_busco_full}"
  ensure_dir "${dir_species_busco_short}"
  species_cds_fasta=()
  mapfile -t species_cds_fasta < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  echo "Number of CDS files for BUSCO: ${#species_cds_fasta[@]}"
  if [[ ${#species_cds_fasta[@]} -eq 0 ]]; then
    echo "No CDS file found. Exiting."
    exit 1
  fi
  for cds_full in "${species_cds_fasta[@]}"; do
    cds=$(basename "${cds_full}")
    sp_ub=$(gg_species_name_from_path "${cds}")
    file_sp_busco_full="${dir_species_busco_full}/${sp_ub}.busco.full.tsv"
    file_sp_busco_short="${dir_species_busco_short}/${sp_ub}.busco.short.txt"
    if [[ ! -s ${file_sp_busco_full} || ! -s ${file_sp_busco_short} ]]; then
      gg_step_start "${task}: ${cds}"

      if [[ "${cds}" == *gz ]]; then # Do not quote like "*gz"
        echo "Decompressing gzipped CDS fasta for BUSCO: ${cds}"
        seqkit seq --threads ${NSLOTS} "${dir_sp_cds}/${cds}" --out-file "tmp.busco_input.cds.fasta"
      else
        cp_out ${dir_sp_cds}/${cds} ./tmp.busco_input.cds.fasta
      fi

      dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}")
      if [[ $? -ne 0 ]]; then
        echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
        exit 1
      fi
      dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

      busco \
      --in "tmp.busco_input.cds.fasta" \
      --mode "transcriptome" \
      --out "busco_tmp" \
      --cpu ${NSLOTS} \
      --force \
      --evalue 1e-03 \
      --limit 20 \
      --lineage_dataset ${dir_busco_lineage} \
      --download_path ${dir_busco_db} \
      --offline

      if [[ $? -eq 0 ]]; then
        if copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
          rm -r './busco_tmp'
        else
          echo "Failed to locate normalized BUSCO outputs for ${sp_ub}. Exiting."
          exit 1
        fi
      fi
    else
      echo "Skipped BUSCO: ${cds}"
    fi
  done
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="Collecting IDs of common BUSCO genes"
if [[ ! -s ${file_busco_summary_table} && ${run_get_busco_summary} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_busco_summary_table}"

  python ${dir_myscript}/collect_common_BUSCO_genes.py \
  --busco_outdir ${dir_species_busco_full} \
  --ncpu ${NSLOTS} \
  --outfile "tmp.busco_summary_table.tsv"
  mv_out "tmp.busco_summary_table.tsv" ${file_busco_summary_table}

  num_busco_ids=$(( $(wc -l < "${file_busco_summary_table}") - 1 ))
  echo "Number of BUSCO genes: ${num_busco_ids}"
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="Generating fasta files for individual single-copy genes"
if [[ ${run_busco_getfasta} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_fasta}"
  busco_rows=()
  mapfile -t busco_rows < <(tail -n +2 "${file_busco_summary_table}")
  num_busco_ids=${#busco_rows[@]}

  generate_single_copy_fasta() {
    local busco_idx=$1
    local busco_row="${busco_rows[${busco_idx}]}"
    local busco_id
    local outfile2
    local cols=()
    local genes=()
    local genes2=()
    IFS=$'\t' read -r -a cols <<< "${busco_row}"
    busco_id="${cols[0]:-}"
    if [[ -z "${busco_id}" ]]; then
      return 0
    fi
    if [[ ${#cols[@]} -gt 3 ]]; then
      genes=( "${cols[@]:3}" )
    fi
    outfile2="${dir_busco_fasta}/${busco_id}.busco.cds.fasta"
    if [[ -s ${outfile2} ]]; then
      return 0
    fi
    echo "busco_id: ${busco_id}"
    if [[ ! -s ${outfile2} ]]; then
      mapfile -t genes2 < <(gg_busco_gene_tokens "split_duplicated" "${genes[@]}")
      if [[ ${#genes2[@]} -eq 0 ]]; then
        echo "Skipping. ${busco_id} has no genes in the selected species."
        return 0
      fi
      gg_seqkit_grep_by_patterns_from_infile_list 1 "species_cds_fasta_list.txt" "${genes2[@]}" \
      | gg_prepare_cds_fasta_stream 1 \
      > ${outfile2}
      if [[ ! -s ${outfile2} ]]; then
        echo "File is empty. Removing: ${outfile2}"
        rm ${outfile2}
      fi
    fi
  }

  gg_find_fasta_files "${dir_sp_cds}" 1 > species_cds_fasta_list.txt
  for (( busco_idx=0; busco_idx<num_busco_ids; busco_idx++ )); do
    wait_until_jobn_le ${NSLOTS}
    generate_single_copy_fasta "${busco_idx}" &
  done
  wait
  rm tmp.*
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment of duplicate-containing BUSCO genes"
if [[ ${run_busco_mafft} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_mafft}"

  run_mafft() {
    infile=$1
    infile_base=${infile%%.*}
    outfile=${dir_busco_mafft}/${infile_base}.busco.cds.aln.fasta
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"

    cdskit mask \
    --seqfile ${dir_busco_fasta}/${infile} \
    --outfile tmp.${infile_base}.cds.fasta

    num_seq=$(gg_count_fasta_records "tmp.${infile_base}.cds.fasta")
    if [[ ${num_seq} -lt 2 ]]; then
      echo "Skipped MAFFT/backalign because fewer than 2 sequences were found: ${infile}"
      cp_out "tmp.${infile_base}.cds.fasta" "${outfile}"
      rm tmp.${infile_base}*
      return 0
    fi

    seqkit translate \
    --allow-unknown-codon \
    --transl-table ${genetic_code} \
    --threads 1 \
    tmp.${infile_base}.cds.fasta \
    > tmp.${infile_base}.pep.fasta

    mafft \
    --auto \
    --amino \
    --thread 1 \
    tmp.${infile_base}.pep.fasta \
    > tmp.${infile_base}.pep.aln.fasta

    cdskit backalign \
    --seqfile tmp.${infile_base}.cds.fasta \
    --aa_aln tmp.${infile_base}.pep.aln.fasta \
    --codontable ${genetic_code} \
    --outfile tmp.${infile_base}.cds.aln.fasta

    if [[ -s tmp.${infile_base}.cds.aln.fasta ]]; then
      cp_out tmp.${infile_base}.cds.aln.fasta ${outfile}
    fi
    rm tmp.${infile_base}*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_fasta}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_mafft ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="TrimAl of duplicate-containing BUSCO genes"
if [[ ${run_busco_trimal} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_trimal}"

  run_trimal() {
    infile=$1
    infile_base=${infile%%.*}
    outfile="${dir_busco_trimal}/${infile_base}.busco.trimal.fasta"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    seqkit seq --remove-gaps --threads 1 ${dir_busco_mafft}/${infile} > tmp.${infile_base}.degap.fasta
    seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_busco_mafft}/${infile} > tmp.${infile_base}.pep.fasta

    trimal \
    -in tmp.${infile_base}.pep.fasta \
    -backtrans tmp.${infile_base}.degap.fasta \
    -out tmp.${infile_base}.trimal.fasta \
    -ignorestopcodon \
    -automated1

    if [[ -s tmp.${infile_base}.trimal.fasta ]]; then
      mv_out tmp.${infile_base}.trimal.fasta ${outfile}
    fi
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_mafft}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_trimal ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_iqtree_dna} -eq 1 ]]; then
    gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_dna}"

  busco_iqtree_dna() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    cp_out "${indir}/${infile}" "./tmp.${infile_base}.input.fasta"

    iqtree \
    -s "./tmp.${infile_base}.input.fasta" \
    -m ${nucleotide_model} \
    -T 1 \
    --prefix tmp.${infile_base} \
    --seed 12345 \
    --redo

    mv_out tmp.${infile_base}.treefile ${outfile}
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_iqtree_dna ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO protein trees"
if [[ ${run_busco_iqtree_pep} -eq 1 ]]; then
    gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_pep}"

  busco_iqtree_pep() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    seqkit translate --transl-table ${genetic_code} --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"

    iqtree \
    -s "./tmp.${infile_base}.input.fasta" \
    -m ${protein_model} \
    -st AA \
    -T 1 \
    --prefix tmp.${infile_base} \
    --seed 12345 \
    --redo

    mv_out tmp.${infile_base}.treefile ${outfile}
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_iqtree_pep ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_notung () {
  infile=$1
  indir=$2
  outdir=$3
  busco_id=${infile%%.*}
  outfile="${outdir}/${busco_id}.busco.notung.root.zip"
  if [[ -s "${outfile}" ]]; then
    return 0
  fi
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -r "./${busco_id}.notung.root"
  fi
  java -jar -Xmx${memory_notung}g ${notung_jar} \
  -s "${file_dated_species_tree}" \
  -g "${indir}/${infile}" \
  --root \
  --infertransfers "false" \
  --treeoutput newick \
  --log \
  --treestats \
  --events \
  --parsable \
  --speciestag prefix \
  --allopt \
  --maxtrees 1000 \
  --nolosses \
  --outputdir "./${busco_id}.notung.root"
  if [[ -e "${busco_id}.notung.root/${busco_id}.busco.nwk.rooting.0" ]]; then
    zip -rq "${busco_id}.busco.notung.root.zip" "${busco_id}.notung.root"
    mv_out "${busco_id}.busco.notung.root.zip" "${outfile}"
    rm -r "./${busco_id}.notung.root"
  fi
}

task="NOTUNG rooting of duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_notung_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${NSLOTS}/2))
    busco_notung ${infile} ${dir_busco_iqtree_dna} ${dir_busco_notung_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="NOTUNG rooting of duplicate-containing BUSCO protein trees"
if [[ ${run_busco_notung_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${NSLOTS}/2))
    busco_notung ${infile} ${dir_busco_iqtree_pep} ${dir_busco_notung_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_species_tree_assisted_gene_tree_rooting () {
  infile=$1
  indir=$2
  intreedir=$3
  outdir_txt=$4
  outdir_nwk=$5
  busco_id=${infile%%.*}
  intree="${intreedir}/${busco_id}.busco.nwk"
  outfile_txt="${outdir_txt}/${busco_id}.busco.root.txt"
  outfile_nwk="${outdir_nwk}/${busco_id}.busco.root.nwk"
  if [[ -s ${outfile_txt} && -s ${outfile_nwk} ]]; then
    return 0
  fi
  echo "Start NOTUNG root: ${busco_id}"
  if [[ -e ./${busco_id}.notung.root ]]; then
    rm -r ./${busco_id}.notung.root
  fi
  cp_out ${indir}/${infile} .
  unzip -q ${infile}

  Rscript ${dir_myscript}/gene_tree_rooting.r \
  ./${busco_id}.notung.root/ \
  ${intree} \
  ${busco_id}.root.nwk \
  ${NSLOTS} \
  2>&1 | tee ${busco_id}.root.txt

  if [[ -s ${busco_id}.root.nwk ]]; then
    mv_out ${busco_id}.root.txt ${outfile_txt}
    mv_out ${busco_id}.root.nwk ${outfile_nwk}
    rm ${infile}
    rm -r ${busco_id}.notung.root
  fi
}

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_root_dna} -eq 1 ]]; then
    gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_species_tree_assisted_gene_tree_rooting ${infile} ${dir_busco_notung_dna} ${dir_busco_iqtree_dna} ${dir_busco_rooted_txt_dna} ${dir_busco_rooted_nwk_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO protein trees"
if [[ ${run_busco_root_pep} -eq 1 ]]; then
    gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_species_tree_assisted_gene_tree_rooting ${infile} ${dir_busco_notung_pep} ${dir_busco_iqtree_pep} ${dir_busco_rooted_txt_pep} ${dir_busco_rooted_nwk_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_grampa () {
  indir=$1
  outdir=$2
  outfile=$3
  ensure_dir "${outdir}"
  local nwk_files=()
  shopt -s nullglob
  nwk_files=( "${indir}"/*.nwk )
  shopt -u nullglob
  if [[ ${#nwk_files[@]} -eq 0 ]]; then
    echo "Skipping Grampa because no rooted gene trees were found in: ${indir}"
    return 0
  fi
  if [[ -e "./grampa_out" ]]; then
    rm -r "./grampa_out"
  fi
  nwkit drop --infile ${file_dated_species_tree} --target 'intnode' --name 'yes' \
  | sed -e "s/_/-/g" \
  > "grampa_input_species_tree.nwk"

  cat "${nwk_files[@]}" \
  | sed -E "s/([(,])([^_)(,:]+)_([^_)(,:]+)_([^)(,:]+)([)(,:])/\1\4|||\2\-\3\5/g" \
  | sed -e "s/_/-/g" -e "s/|||/_/g" \
  > "grampa_input_gene_trees.nwk"

  printf "%s\n" "${nwk_files[@]##*/}" > "busco_genetree_filenames.txt"

  if [[ ${grampa_h1} == "" ]]; then
    h1_param=""
  else
    local grampa_h1_normalized=${grampa_h1//_/-}
    grampa_h1_normalized=${grampa_h1_normalized//[[:space:]]/-}
    h1_param="-h1 ${grampa_h1_normalized}"
  fi

  grampa.py \
  -s "grampa_input_species_tree.nwk" \
  -g "grampa_input_gene_trees.nwk" \
  -o "./grampa_out" \
  -p ${NSLOTS} \
  -v -1 \
  --maps \
  ${h1_param}

  if [[ -s "./grampa_out/grampa_trees_filtered.txt" ]]; then
    # Filtered gene trees are not analyzed by Grampa, but it does not disturb gene tree IDs.
    # For example, if the 136th gene tree is filtered out, it is replaced with a placeholder text in grampa_trees_filtered.txt.
    # And GT-136 does not appear in the Grampa outputs. Still, the 137th gene tree is labeled correctly as GT-137.
    sed -e "s/$/;/" "./grampa_out/grampa_trees_filtered.txt" > "grampa_input_gene_trees_filtered.nwk"
  fi
  python ${dir_myscript}/parse_grampa.py \
  --grampa_det "./grampa_out/grampa_det.txt" \
  --grampa_out "./grampa_out/grampa_out.txt" \
  --gene_trees "./grampa_input_gene_trees.nwk" \
  --species_tree "./grampa_input_species_tree.nwk" \
  --ncpu ${NSLOTS} \
  --sorted_gene_tree_file_names "./busco_genetree_filenames.txt"

  if [[ -s "grampa_out/grampa_checknums.txt" && -s "grampa_out/grampa_det.txt" && -s "grampa_out/grampa_out.txt" && -s "grampa_summary.tsv" ]]; then
    awk -F'\t' '/^The MUL-tree with the minimum parsimony score/ {print $NF}' "./grampa_out/grampa_out.txt" > "${outdir}/best_mul_tree.nwk"
    mv_out "./grampa_out/grampa_checknums.txt" "${outdir}"
    mv_out "./grampa_out/grampa_det.txt" "${outdir}"
    mv_out "./grampa_out/grampa_out.txt" "${outdir}"
    if [[ -s "./grampa_input_gene_trees_filtered.nwk" ]]; then
      mv_out "./grampa_input_gene_trees_filtered.nwk" "${outdir}"
    fi
    mv_out "./grampa_input_species_tree.nwk" "${outdir}"
    mv_out "./grampa_input_gene_trees.nwk" "${outdir}"
    mv_out "./grampa_summary.tsv" "${outfile}"
    rm -r "./grampa_out"
  else
    echo "Grampa output files are missing. Ending Grampa."
  fi
}

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO DNA trees"
if [[ ! -s ${file_busco_grampa_dna} && ${run_busco_grampa_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa ${dir_busco_rooted_nwk_dna} ${dir_busco_grampa_dna} ${file_busco_grampa_dna}
else
  gg_step_skip "${task}"
fi

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO protein trees"
if [[ ! -s ${file_busco_grampa_pep} && ${run_busco_grampa_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa ${dir_busco_rooted_nwk_pep} ${dir_busco_grampa_pep} ${file_busco_grampa_pep}
else
  gg_step_skip "${task}"
fi

task="Orthogroup-based Grampa analysis for polyploidization history"
if [[ ! -s ${file_orthogroup_grampa} && ${run_orthogroup_grampa} -eq 1 ]]; then
  gg_step_start "${task}"

  og_ids=()
  mapfile -t og_ids < <(python -c 'import pandas,sys; d=pandas.read_csv(sys.argv[1],sep="\t",header=0); ids=d.loc[(d["Total"]>=int(sys.argv[2]))&(d["Total"]<=int(sys.argv[3])),"Orthogroup"].astype(str).tolist(); print("\n".join(ids))' ${file_genecount} ${min_gene_orthogroup_grampa} ${max_gene_orthogroup_grampa})
  file_names=()
  mapfile -t file_names < <(gg_find_file_basenames "${dir_og_rooted_tree}")
  echo "Number of files in ${dir_og_rooted_tree}: ${#file_names[@]}"
  echo "Number of selected orthogroups with ${min_gene_orthogroup_grampa}<=gene number<=${max_gene_orthogroup_grampa}: ${#og_ids[@]}"
  if [[ -e ./tmp.orthogroup_grampa_indir ]]; then
    rm -r ./tmp.orthogroup_grampa_indir
  fi
  mkdir ./tmp.orthogroup_grampa_indir
  for file_name in "${file_names[@]}"; do
    for og_id in "${og_ids[@]}"; do
      if [[ ${file_name} == ${og_id}* ]]; then
        cp_out ${dir_og_rooted_tree}/${file_name} ./tmp.orthogroup_grampa_indir
        mapfile -t og_ids < <(printf "%s\n" "${og_ids[@]}" | grep -v -Fx "${og_id}" || true)
        break
      fi
    done
  done

  busco_grampa ./tmp.orthogroup_grampa_indir ${dir_orthogroup_grampa} ${file_orthogroup_grampa}
else
  gg_step_skip "${task}"
fi

task='CAFE analysis'
disable_if_no_input_file "run_cafe" ${file_genecount} ${file_dated_species_tree}
if [[ (! -s "${dir_cafe_summary_plot}/summary_all.pdf" || ! -s "${dir_cafe_summary_plot}/summary_significant.pdf") && ${run_cafe} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_cafe}"

  if [[ ! -s "${dir_cafe_output}/Gamma_asr.tre" || ! -s "${dir_cafe_output}/Gamma_count.tab" || ! -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d ${dir_cafe_orthogroup_selection} ]]; then
      rm -r ${dir_cafe_orthogroup_selection}
    fi
    python ${dir_myscript}/cafe_orthogroup_selection.py \
    --genecount ${file_genecount} \
    --dated_species_tree ${file_dated_species_tree} \
    --output_dir ${dir_cafe_orthogroup_selection} \
    --max_size_differential ${max_size_differential_cafe}

    if [[ -d ${dir_cafe_output} ]]; then
      rm -r ${dir_cafe_output}
    fi
    cafe5 \
    --infile ${dir_cafe_orthogroup_selection}/cafe_input.tsv \
    --tree ${file_dated_species_tree} \
    --n_gamma_cats ${n_gamma_cats_cafe} \
    --pvalue 0.05 \
    --cores ${NSLOTS} \
    --output_prefix ${dir_cafe_output}
  else
    echo "CAFE output files already exist. Skipping CAFE run."
  fi

  if [[ -s "${dir_cafe_output}/Gamma_asr.tre" && -s "${dir_cafe_output}/Gamma_count.tab" && -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d ${dir_cafe_each_family_plot} ]]; then
      rm -r ${dir_cafe_each_family_plot}
    fi
    Rscript ${dir_myscript}/cafe_plot_each_family.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_output}/Gamma_count.tab \
    ${dir_cafe_output}/Gamma_change.tab \
    ${dir_cafe_each_family_plot} \
    ${NSLOTS}
    status=$?
    if [[ ${status} -ne 0 ]]; then
      echo "Error in Rscript cafe_plot_each_family.r. Exiting."
      exit 1
    fi

    if [[ -d ${dir_cafe_summary_plot} ]]; then
      rm -r ${dir_cafe_summary_plot}
    fi
    Rscript ${dir_myscript}/cafe_plot_summary.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_output}/Gamma_change.tab \
    ${dir_cafe_summary_plot}

    Rscript ${dir_myscript}/cafe_plot_branch_id.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_summary_plot}
  else
    echo "CAFE did not finish successfully. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="GO enrichment analysis"
disable_if_no_input_file "run_go_enrichment" "${dir_cafe_output}/Gamma_change.tab" "${dir_cafe_output}/Gamma_branch_probabilities.tab" ${file_gene_id} ${file_go_annotation}
if [[(! -s "${dir_go_enrichment}/enrichment_significant_${change_direction_go}_${target_branch_go}_significant_go.tsv") && ${run_go_enrichment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_go_enrichment}"
  Rscript ${dir_myscript}/cafe_go_enrichment.r \
  "${dir_cafe_output}/Gamma_change.tab" \
  "${dir_cafe_output}/Gamma_branch_probabilities.tab" \
  "${file_gene_id}" \
  "${file_go_annotation}" \
  "${dir_go_enrichment}" \
  "${target_branch_go}" \
  "${change_direction_go}" \
  "${go_category}"
  status=$?
  if [[ ${status} -ne 0 ]]; then
    echo "Error in Rscript cafe_go_enrichment.r. Exiting."
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  fi
fi

echo $(date): end
