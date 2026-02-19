#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

mode_sraid=1 # Need input at workspace/input/species_sra_list. File name should be GENUS_SPECIES.txt. One SRA or BioProject ID per line.
mode_fastq=0 # Need input at workspace/input/species_rnaseq. Multiple fastq can be placed under workspace/input/species_rnaseq/GENUS_SPECIES.
mode_metadata=0 # Need input at workspace/input/transcriptome_assembly/amalgkit_metadata. File name should be 'GENUS_SPECIES.metadata.tsv'. The 'scientific_name' column should be set to 'GENUS SPECIES'.

run_amalgkit_metadata_or_integrate=1 # Metadata retrieval.
run_amalgkit_getfastq=1 # fastq generation from NCBI SRA.
run_rRNA_contamination_report=1 # Runs amalgkit getfastq prepared .fastq files against rRNA database to determine rRNA contamination.
run_assembly=1 # Transcriptome assembly with Trinity or rnaSPAdes.
run_longestcds=1 # Longest CDS extraction.
run_longestcds_fx2tab=1 # Compute sequence stats of the longest CDS sequences. https://bioinf.shenwei.me/seqkit/usage/#fx2tab-tab2fx
run_longestcds_mmseqs2taxonomy=0 # Taxonomic assignment of longest CDS sequences. 16CPUS * 8GRAM should work. https://github.com/timkahlke/BASTA
run_longestcds_contamination_removal=0 # Removal of taxonomically non-compatible contaminated sequences using the result of MMseqs2 Taxonomy.
run_busco1=1 # BUSCO for transcriptome assembly that contains isoforms.
run_busco2=1 # BUSCO for longest CDS
run_busco3=0 # BUSCO for contamination-removed longest CDS
run_assembly_stat=1 # seqkit stat for assembly and extracted CDS files.
run_amalgkit_quant=1 # Expression level quantification against the transcripts with longest CDS or species_cds.
run_amalgkit_merge=1 # Merging expression levels across samples. This command will run for each species per batch job, unlike the intended use of `amalgkit merge`, which can process all species at once.
run_multispecies_summary=1 # Generate plots and tables that summarize annotation info across species.

remove_amalgkit_fastq_after_completion=1 # To save disk space after getting all output files.
max_assembly_input_fastq_size="30,000,000,000" # The maximum input total read length (bp) for transcriptome assembly. If the input file size is larger than this, reads will be randomly sub-sampled. This is to avoid memory overflow in transcriptome assembly. 30,000,000,000 bp = 30 GB.
assembly_method="Trinity" # {Trinity,rnaSPAdes}
protocol_rna_seq="mixed" # {same,mixed} Option for rnaSPAdes. If "mixed", top 9 fastq files by size are selected, as rnaSPAdes supports fewer than 10 libraries with mixed protocols. https://ablab.github.io/spades/rna.html
busco_lineage="embryophyta_odb12" # See here for available datasets: https://busco-data.ezlab.org/v5/data/lineages/
contamination_removal_rank="phylum" # {phylum, class, order, family, genus, species} The rank of taxonomic assignment to be used for contamination removal.
kallisto_reference="longest_cds" # {species_cds, longest_transcript, longest_cds, contamination_removed_longest_cds}
orf_aggregation_level="i" # {c,g,i,p}, {aggressive <-> moderate}, Which clustering level ORFs are aggregated at, assuming Trinity gene naming like "Drosophyllum_lusitanicum_DN269437−c0−g1−i1.p1". "c" is not supported for rnaSPAdes naming like "Drosophyllum_lusitanicum_g1−i1.p1".
assembly_cpu_offset=0 # $((NSLOTS/4)) # In NIG, Trinity jobs are often killed if CPU or RAM is fully used.
assembly_ram_offset=4 # GB In NIG, Trinity jobs are often killed if CPU or RAM is fully used.

### End: Modify this block to tailor your analysis ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

# Setting modes
if enable_all_run_flags_for_debug_mode; then
  echo "gg debug mode: max_assembly_input_fastq_size is set to 30,000 bp."
  max_assembly_input_fastq_size="30,000,000"
fi

dir_transcriptome_assembly_input="${dir_pg_input}/transcriptome_assembly"
dir_transcriptome_assembly_output="${dir_pg_output}/transcriptome_assembly"
dir_input_fastq="${dir_pg_input}/species_rnaseq"
dir_input_sra_list="${dir_pg_input}/species_sra_list"
dir_tmp_main="${dir_transcriptome_assembly_output}/tmp"
dir_amalgkit_metadata="${dir_transcriptome_assembly_input}/amalgkit_metadata"
dir_amalgkit_getfastq="${dir_transcriptome_assembly_output}/amalgkit_getfastq"
dir_rRNA_contamination_report="${dir_transcriptome_assembly_output}/amalgkit_metadata_rRNA_mapping_rate"
dir_assembly="${dir_transcriptome_assembly_output}/assembled_transcripts_with_isoforms"
dir_longestcds="${dir_transcriptome_assembly_output}/longest_cds"
dir_longestcds_transcript="${dir_transcriptome_assembly_output}/longest_cds_transcript"
dir_longestcds_fx2tab="${dir_transcriptome_assembly_output}/longest_cds_fx2tab"
dir_longestcds_mmseqs2taxonomy="${dir_transcriptome_assembly_output}/longest_cds_mmseqs2taxonomy"
dir_longestcds_contamination_removal_fasta="${dir_transcriptome_assembly_output}/longest_cds_contamination_removal_fasta"
dir_longestcds_contamination_removal_tsv="${dir_transcriptome_assembly_output}/longest_cds_contamination_removal_tsv"
dir_assembly_stat="${dir_transcriptome_assembly_output}/assembly_stat"
dir_busco_full1="${dir_transcriptome_assembly_output}/busco_full_cdna_isoforms"
dir_busco_short1="${dir_transcriptome_assembly_output}/busco_short_cdna_isoforms"
dir_busco_full2="${dir_transcriptome_assembly_output}/busco_full_longest_cds"
dir_busco_short2="${dir_transcriptome_assembly_output}/busco_short_longest_cds"
dir_busco_full3="${dir_transcriptome_assembly_output}/busco_full_longest_cds_contamination_removal"
dir_busco_short3="${dir_transcriptome_assembly_output}/busco_short_longest_cds_contamination_removal"
dir_amalgkit_quant="${dir_transcriptome_assembly_output}/amalgkit_quant"
dir_amalgkit_merge="${dir_transcriptome_assembly_output}/amalgkit_merge"
dir_multispecies_summary="${dir_transcriptome_assembly_output}/multispecies_summary"
dir_mmseqs2_db="${dir_pg_output}/db_mmseqs2"

# Backward compatibility for legacy input layout:
# - workspace/input/transcriptome_assembly/input_fastq
# - workspace/input/transcriptome_assembly/input_sra_list
if [[ ! -d "${dir_input_fastq}" && -d "${dir_transcriptome_assembly_input}/input_fastq" ]]; then
  dir_input_fastq="${dir_transcriptome_assembly_input}/input_fastq"
fi
if [[ ! -d "${dir_input_sra_list}" && -d "${dir_transcriptome_assembly_input}/input_sra_list" ]]; then
  dir_input_sra_list="${dir_transcriptome_assembly_input}/input_sra_list"
fi

if [[ ${mode_fastq} -eq 1 && ! -d "${dir_input_fastq}" ]]; then
  echo "Missing mode_fastq input directory: ${dir_input_fastq}"
  mode_fastq=0
fi
if [[ ${mode_sraid} -eq 1 ]]; then
  if [[ ! -d "${dir_input_sra_list}" ]]; then
    echo "Missing mode_sraid input directory: ${dir_input_sra_list}"
    mode_sraid=0
  else
    shopt -s nullglob
    sra_mode_files=( "${dir_input_sra_list}"/* )
    shopt -u nullglob
    if [[ ${#sra_mode_files[@]} -eq 0 ]]; then
      echo "Input directory is empty for mode_sraid: ${dir_input_sra_list}"
      mode_sraid=0
    fi
  fi
fi
if [[ ${mode_metadata} -eq 1 ]]; then
  if [[ ! -d "${dir_amalgkit_metadata}" ]]; then
    echo "Missing mode_metadata input directory: ${dir_amalgkit_metadata}"
    mode_metadata=0
  else
    shopt -s nullglob
    metadata_mode_files=( "${dir_amalgkit_metadata}"/* )
    shopt -u nullglob
    if [[ ${#metadata_mode_files[@]} -eq 0 ]]; then
      echo "Input directory is empty for mode_metadata: ${dir_amalgkit_metadata}"
      mode_metadata=0
    fi
  fi
fi
if [[ ${mode_fastq} -eq 0 && ${mode_sraid} -eq 0 && ${mode_metadata} -eq 0 ]]; then
  echo "No valid transcriptome input mode is available. Skipping transcriptome assembly workflow."
  exit 0
fi

if [[ ${mode_fastq} -eq 1 ]]; then
  echo 'Mode: fastq input'
  if [[ ! -d "${dir_input_fastq}" ]]; then
    echo "Input directory does not exist: ${dir_input_fastq}"
    exit 1
  fi
  shopt -s nullglob
  dirs=( "${dir_input_fastq}"/* )
  shopt -u nullglob
  id=$((SGE_TASK_ID-1))
  dir_species_fastq="${dirs[${id}]}"
  sp_ub="$(basename "${dir_species_fastq}")"
  shopt -s nullglob
  files_fastq=( "${dir_species_fastq}"/* )
  shopt -u nullglob
  echo "Input fastq directory: ${dir_species_fastq}"
  echo "Species: ${sp_ub}"
  echo "Input fastq files: ${files_fastq[@]}"
  if [[ ! -d "${dir_species_fastq}" || ${#files_fastq[@]} -eq 0 ]]; then
    echo "Input fastq file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
elif [[ ${mode_sraid} -eq 1 ]]; then
  echo 'Mode: sraid input'
  if [[ ! -d "${dir_input_sra_list}" ]]; then
    echo "Input directory does not exist: ${dir_input_sra_list}"
    exit 1
  fi
  shopt -s nullglob
  files=( "${dir_input_sra_list}"/* )
  shopt -u nullglob
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "Input directory is empty: ${dir_input_sra_list}"
    exit 1
  fi
  id=$((SGE_TASK_ID-1))
  file_input_sra_list="${files[${id}]}"
  if [[ ! -f "${file_input_sra_list}" ]]; then
    echo "Input SRA list file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
  sra_ids=( $(cat "${file_input_sra_list}") )
  sp_ub="$(gg_species_name_from_path_or_dot "${file_input_sra_list}")"
  echo "Input SRA list file: ${file_input_sra_list}"
  echo "Species: ${sp_ub}"
  echo "Number of input SRA IDs: ${#sra_ids[@]}"
  echo "Input SRA IDs: ${sra_ids[@]}"
  if [[ ${#sra_ids[@]} -eq 0 ]]; then
    echo "SRA IDs not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
elif [[ ${mode_metadata} -eq 1 ]]; then
  echo 'Mode: metadata input'
  if [[ ! -d "${dir_amalgkit_metadata}" ]]; then
    echo "Input directory does not exist: ${dir_amalgkit_metadata}"
    exit 1
  fi
  shopt -s nullglob
  files=( "${dir_amalgkit_metadata}"/* )
  shopt -u nullglob
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "Input directory is empty: ${dir_amalgkit_metadata}"
    exit 1
  fi
  id=$((SGE_TASK_ID-1))
  file_metadata="${files[${id}]}"
  sp_ub="$(gg_species_name_from_path_or_dot "${file_metadata}")"
  echo "Input metadata file: ${file_metadata}"
  echo "Species: ${sp_ub}"
  if [[ ! -f "${file_metadata}" ]]; then
    echo "Input metadata file not found, probably due to the out-of-range specification of array jobs. Exiting."
    exit 1
  fi
fi

dir_tmp="${dir_tmp_main}/${SGE_TASK_ID}_${sp_ub}"
dir_amalgkit_getfastq_sp="${dir_amalgkit_getfastq}/${sp_ub}"
file_amalgkit_metadata="${dir_amalgkit_metadata}/${sp_ub}.metadata.tsv"
file_amalgkit_getfastq_safely_removed_flag=${dir_amalgkit_getfastq}/${sp_ub}.safely_removed.txt
file_rRNA_contamination_report="${dir_rRNA_contamination_report}/${sp_ub}_rRNA_mapping_rate.tsv"
file_isoform="${dir_assembly}/${sp_ub}.isoform.fa.gz"
file_longestcds="${dir_longestcds}/${sp_ub}_longestCDS.fa.gz"
file_longestcds_transcript="${dir_longestcds_transcript}/${sp_ub}_longestCDS.transcript.fa.gz"
file_longestcds_fx2tab="${dir_longestcds_fx2tab}/${sp_ub}_longestCDS.fx2tab_cds.tsv"
file_longestcds_mmseqs2taxonomy="${dir_longestcds_mmseqs2taxonomy}/${sp_ub}_longestCDS.mmseqs2taxonomy.tsv"
file_longestcds_contamination_removal_fasta="${dir_longestcds_contamination_removal_fasta}/${sp_ub}_longestCDS_contamination_removal.fa.gz"
file_longestcds_contamination_removal_tsv="${dir_longestcds_contamination_removal_tsv}/${sp_ub}_longestCDS_contamination_removal.tsv"
file_assembly_stat="${dir_assembly_stat}/${sp_ub}_assembly_stat.tsv"
file_busco_full1="${dir_busco_full1}/${sp_ub}.busco.full.tsv"
file_busco_short1="${dir_busco_short1}/${sp_ub}.busco.short.txt"
file_busco_full2="${dir_busco_full2}/${sp_ub}.busco.full.tsv"
file_busco_short2="${dir_busco_short2}/${sp_ub}.busco.short.txt"
file_busco_full3="${dir_busco_full3}/${sp_ub}.busco.full.tsv"
file_busco_short3="${dir_busco_short3}/${sp_ub}.busco.short.txt"
file_amalgkit_merge_efflen="${dir_amalgkit_merge}/${sp_ub}/${sp_ub}_eff_length.tsv"
file_amalgkit_merge_count="${dir_amalgkit_merge}/${sp_ub}/${sp_ub}_est_counts.tsv"
file_amalgkit_merge_tpm="${dir_amalgkit_merge}/${sp_ub}/${sp_ub}_tpm.tsv"
file_amalgkit_merge_metadata="${dir_amalgkit_merge}/${sp_ub}/${sp_ub}_metadata.tsv"
file_multispecies_summary="${dir_multispecies_summary}/assembly_stat_summary.pdf"

ensure_dir "${dir_tmp}"
cd ${dir_tmp}

task="amalgkit metadata/integrate"
if [[ ! -s ${file_amalgkit_metadata} && ${run_amalgkit_metadata_or_integrate} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  if [[ -e "./metadata" ]]; then
    rm -r "./metadata"
  fi

  if [[ ${mode_sraid} -eq 1 ]]; then
    search_string=$(cat ${file_input_sra_list} | tr -s "\n" | sed ':a;N;$!ba;s/\n/ OR /g' | sed -e "s/ OR $//" -e "s/^/(/" -e "s/$/)/")
    search_string=$(echo "${search_string}" | sed -e 's/$/ AND "Illumina"[Platform] AND ("RNA-seq"[Strategy] OR "EST"[Strategy])/')
    echo "Entrez search string: ${search_string}"

    amalgkit metadata \
    --out_dir "./" \
    --search_string "${search_string}"

    sp_space=$(echo ${sp_ub} | sed -e "s/_/ /g")
    {
        head -n 1 "./metadata/metadata.tsv";
        cat "./metadata/metadata.tsv" | grep -e "${sp_space}";
    } | sed -e "s/\t\t\tno\t/\tyes\tyes\tno\t/g" > "./metadata.tsv"

  elif [[ ${mode_fastq} -eq 1 ]]; then
    amalgkit integrate \
    --out_dir "./" \
    --fastq_dir ${dir_species_fastq} \
    --threads ${NSLOTS} \
    --remove_tmp yes

    mv ./metadata_private_fastq.tsv ./metadata/metadata.tsv
    python -c "import pandas as pd; d=pd.read_csv('./metadata/metadata.tsv',sep='\t',header=0); d.loc[:,'scientific_name']='${sp_ub}'; d.to_csv('./metadata.tsv',sep='\t',index=False)"
  fi

  if [[ -s "./metadata.tsv" ]]; then
    mv_out "./metadata.tsv" ${file_amalgkit_metadata}
  else
    echo "metadata.tsv not found. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task="amalgkit getfastq"
amalgkit_fastq_files=( $(find ${dir_amalgkit_getfastq_sp} -name "*.amalgkit.fastq.gz") )
echo "Number of amalgkit getfastq fastq files: ${#amalgkit_fastq_files[@]}"
echo "is_fastq_requiring_downstream_analysis_done: $(is_fastq_requiring_downstream_analysis_done)"
if [[ ( ${#amalgkit_fastq_files[@]} -eq 0 && ${run_amalgkit_getfastq} -eq 1 ) && $(is_fastq_requiring_downstream_analysis_done) -eq 0 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_dir "${dir_amalgkit_getfastq_sp}"

  if [[ -e ${file_amalgkit_getfastq_safely_removed_flag} ]]; then
    rm ${file_amalgkit_getfastq_safely_removed_flag}
  fi

  amalgkit getfastq \
  --out_dir ${dir_tmp} \
  --metadata ${file_amalgkit_metadata} \
  --threads ${NSLOTS} \
  --remove_sra yes \
  --remove_tmp yes \
  --read_name 'trinity' \
  --aws yes \
  --ncbi yes \
  --redo no
  status_amalgkit=$?
  if [[ ${status_amalgkit} -eq 0 ]]; then
    echo "amalgkit getfastq safely finished."
    mv ${dir_tmp}/getfastq/* ${dir_amalgkit_getfastq_sp}
    rm -r ${dir_tmp}/getfastq
  else
    echo "amalgkit getfastq did not safely finish. Exiting."
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='rRNA contamination report'
if [[ ( ! -s ${file_rRNA_contamination_report} ) && ${run_rRNA_contamination_report} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  silva_rrna_ref=$(ensure_silva_rrna_ref_db "${dir_pg}")
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare SILVA rRNA reference. Exiting."
    exit 1
  fi

  recreate_dir "./fasta"
  recreate_dir "./index"
  recreate_dir "./quant"
  recreate_dir "./merge"
  if [[ -e "./getfastq" ]]; then
    rm -r "./getfastq"
  fi
  ln -s ${dir_amalgkit_getfastq_sp} "./getfastq"
  path_reference_fasta_link="./fasta/${sp_ub}_dummy_rRNA_for_kallisto_index.fa.gz"
  ln -s ${silva_rrna_ref} ${path_reference_fasta_link}
  kallisto index --make-unique -i "./index/${sp_ub}.idx" ${path_reference_fasta_link}; status_kallisto=$?
  if [[ ${status_kallisto} -ne 0 ]]; then
    echo "kallisto index failed with exit code ${status_kallisto}"
    exit 1
  fi

  amalgkit quant \
  --out_dir "./" \
  --threads ${NSLOTS} \
  --metadata "${file_amalgkit_metadata}" \
  --clean_fastq no \
  --fasta_dir "./fasta" \
  --build_index no
  exit_code_amalgkit_quant=$?

  amalgkit merge \
  --out_dir "./" \
  --metadata "${file_amalgkit_metadata}"
  exit_code_amalgkit_merge=$?

  if [[ ${exit_code_amalgkit_quant} -ne 0 || ${exit_code_amalgkit_merge} -ne 0 ]]; then
    echo "amalgkit quant for rRNA contamination rate estimation failed with exit code ${exit_code_amalgkit_quant}"
    echo "If you got an error 'FileNotFoundError: Could not find index file.', the specified species name may not be the same as the one in the metadata.tsv file."
    echo "In that case, please correct the species name in the metadata.tsv file (or in your input files/directories) and run this pipeline again."
    echo "Exiting."
    exit 1
  else
    echo "amalgkit quant for rRNA contamination rate estimation finished successfully"
    mv_out ./merge/metadata.tsv ${file_rRNA_contamination_report}
    rm ${path_reference_fasta_link}
    rm -r "./fasta"
    rm "./getfastq" # Do not put -r, otherwise the original getfastq files will be deleted.
  fi
else
  echo "$(date): Skipped: ${task}"
fi

task='De novo transcriptome assembly'
if [[ ! -s ${file_isoform} && ${run_assembly} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_parent_dir "${file_isoform}"

  mapfile -t files_right < <(find "${dir_amalgkit_getfastq_sp}" -name "*_2.amalgkit.fastq.gz" | sort)
  if [[ ${#files_right[@]} -eq 0 ]]; then
    echo "Paired-end samples were not detected. SE reads will be used for transcriptome assembly."
    lib_layout='single'
  else
    echo "Paired-end samples were detected. PE reads will be used for transcriptome assembly and SE reads, if any, will not be used."
    lib_layout='paired'
  fi

  selected_fastq_dir="./tmp_selected_fastq"
  if [[ -e ${selected_fastq_dir} ]]; then
    rm -r ${selected_fastq_dir}
  fi
  mkdir ${selected_fastq_dir}
  if [[ ${assembly_method} == "rnaSPAdes" && ${protocol_rna_seq} == "mixed" ]]; then
    if [[ ${lib_layout} == 'single' ]]; then
      mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -name "*.amalgkit.fastq.gz")
    elif [[ ${lib_layout} == 'paired' ]]; then
      mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -name "*_1.amalgkit.fastq.gz")
    fi
    n_pair=${#candidate_fastq_files[@]}
    if [[ ${n_pair} -ge 10 ]]; then
      echo "Selecting top 9 ${lib_layout}-end fastq files by size, as rnaSPAdes accepts fewer than 10 libraries with different protocols."
      echo "For details, please refer to the rnaSPAdes manual: https://ablab.github.io/spades/rna.html"
      echo "Selected fastq files:"
      if [[ ${lib_layout} == 'single' ]]; then
        find "${dir_amalgkit_getfastq_sp}" -name "*.amalgkit.fastq.gz" -exec du -b {} + \
        | sort -nr | head -n 9 | cut -f2 | while read -r file; do
          cp ${file} ${selected_fastq_dir}; echo ${file}
        done
      elif [[ ${lib_layout} == 'paired' ]]; then
        find "${dir_amalgkit_getfastq_sp}" -name "*_1.amalgkit.fastq.gz" -exec du -b {} + \
        | sort -nr | head -n 9 | cut -f2 | while read -r file1; do
          file2="${file1/_1.amalgkit.fastq.gz/_2.amalgkit.fastq.gz}"
          cp ${file1} ${selected_fastq_dir}; echo ${file1}
          cp ${file2} ${selected_fastq_dir}; echo ${file2}
        done
      fi
    else
      echo "All ${lib_layout}-end fastq files will be used for transcriptome assembly."
      selected_fastq_dir=${dir_amalgkit_getfastq_sp}
    fi
  else
    echo "All ${lib_layout}-end fastq files will be used for transcriptome assembly."
    selected_fastq_dir=${dir_amalgkit_getfastq_sp}
  fi

  if [[ ${lib_layout} == 'single' ]]; then
    total_fastq_len=$(get_total_fastq_len ${selected_fastq_dir} "*.amalgkit.fastq.gz")
  elif [[ ${lib_layout} == 'paired' ]]; then
    total_fastq_len1=$(get_total_fastq_len ${selected_fastq_dir} "*_1.amalgkit.fastq.gz")
    total_fastq_len2=$(get_total_fastq_len ${selected_fastq_dir} "*_2.amalgkit.fastq.gz")
    total_fastq_len=$((${total_fastq_len1}+${total_fastq_len2}))
  fi
  max_assembly_input_fastq_size=$(echo ${max_assembly_input_fastq_size} | sed -e "s/,//g")
  if [[ ${total_fastq_len} -gt ${max_assembly_input_fastq_size} ]]; then
    echo "Total ${lib_layout} fastq length is ${total_fastq_len} bp, which is greater than ${max_assembly_input_fastq_size} bp."
    echo "Only ${max_assembly_input_fastq_size} bp will be used."
    assembly_input_fastq_dir="./tmp_assembly_input_fastq"
    if [[ -e ${assembly_input_fastq_dir} ]]; then
      rm -r ${assembly_input_fastq_dir}
    fi
    mkdir ${assembly_input_fastq_dir}
    proportion=$(echo "${max_assembly_input_fastq_size} ${total_fastq_len}" | awk '{printf "%.3f\n", $1/$2}')
    echo "Proportion of fastq reads to be used: ${proportion} (${max_assembly_input_fastq_size}/${total_fastq_len})"
    files=()
    if [[ ${lib_layout} == 'single' ]]; then
      mapfile -t files < <(find "${selected_fastq_dir}" -name "*.amalgkit.fastq.gz" | sort)
    elif [[ ${lib_layout} == 'paired' ]]; then
      mapfile -t files1 < <(find "${selected_fastq_dir}" -name "*_1.amalgkit.fastq.gz" | sort)
      mapfile -t files2 < <(find "${selected_fastq_dir}" -name "*_2.amalgkit.fastq.gz" | sort)
      files=( "${files1[@]}" "${files2[@]}" )
    fi
    for file in "${files[@]}"; do
      seqkit sample --proportion ${proportion} --rand-seed 11 --out-file "${assembly_input_fastq_dir}/$(basename "${file}")" "${file}"
    done
    echo "Total fastq length of the subsampled fastq files: $(get_total_fastq_len ${assembly_input_fastq_dir} "*.amalgkit.fastq.gz") bp"
  else
    echo "Total fastq length is ${total_fastq_len} bp, which is less than ${max_assembly_input_fastq_size} bp. All fastq reads will be used."
    assembly_input_fastq_dir=${selected_fastq_dir}
  fi

  nslots_assembly=$((${NSLOTS}-${assembly_cpu_offset}))
  memory_assembly=$((${MEM_PER_HOST}-${assembly_ram_offset}))
  bflyHeapSpaceMax=${MEM_PER_SLOT} # GB
  echo "Number of offset CPUs for transcriptome assembly is ${assembly_cpu_offset}."
  echo "${NSLOTS} CPUs and ${MEM_PER_HOST}G RAM are allocated to this job."
  echo "${nslots_assembly} CPUs and ${memory_assembly}G RAM are used for transcriptome assembly."

  files_single=()
  files_left=()
  files_right=()
  mapfile -t files_single < <(find "${assembly_input_fastq_dir}" -name "*.amalgkit.fastq.gz" | sort)
  mapfile -t files_left < <(find "${assembly_input_fastq_dir}" -name "*_1.amalgkit.fastq.gz" | sort)
  mapfile -t files_right < <(find "${assembly_input_fastq_dir}" -name "*_2.amalgkit.fastq.gz" | sort)

  if [[ ${lib_layout} == 'single' && ${#files_single[@]} -eq 0 ]]; then
    echo "No single-end fastq files were found for transcriptome assembly. Exiting."
    exit 1
  fi
  if [[ ${lib_layout} == 'paired' ]]; then
    if [[ ${#files_left[@]} -eq 0 || ${#files_right[@]} -eq 0 || ${#files_left[@]} -ne ${#files_right[@]} ]]; then
      echo "Paired-end input files were not detected correctly in ${assembly_input_fastq_dir}."
      echo "Detected left/right counts: ${#files_left[@]}/${#files_right[@]}. Exiting."
      exit 1
    fi
  fi

  if [[ ${assembly_method} == 'Trinity' ]]; then
    if [[ ${lib_layout} == 'single' ]]; then
      in_single="$(IFS=","; echo "${files_single[*]}")"
      trinity_input="--single ${in_single}"
    elif [[ ${lib_layout} == 'paired' ]]; then
      in_left="$(IFS=","; echo "${files_left[*]}")"
      in_right="$(IFS=","; echo "${files_right[*]}")"
      trinity_input="--left ${in_left} --right ${in_right}"
    fi
    Trinity \
    --seqType fq \
    --CPU ${nslots_assembly} \
    --max_memory ${memory_assembly}G \
    --min_contig_length 200 \
    --output trinity \
    --full_cleanup \
    --NO_SEQTK \
    --bflyHeapSpaceMax ${bflyHeapSpaceMax}G \
    --bflyGCThreads 1 \
    --bflyCalculateCPU \
    ${trinity_input}
    # For --NO_SEQTK, see https://github.com/trinityrnaseq/trinityrnaseq/issues/787
    if [[ -s "${dir_tmp}/trinity.Trinity.fasta" ]]; then
      seqkit seq --threads ${NSLOTS} "${dir_tmp}/trinity.Trinity.fasta" --out-file "${file_isoform}"
    fi
  elif [[ ${assembly_method} == 'rnaSPAdes' ]]; then
    if [[ ${protocol_rna_seq} == "same" ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        rnaspades_input=$(for i in "${!files_single[@]}"; do echo -n " --s1 ${files_single[i]}"; done)
      elif [[ ${lib_layout} == 'paired' ]]; then
        rnaspades_input=$(for i in "${!files_left[@]}"; do echo -n " --pe1-1 ${files_left[i]} --pe1-2 ${files_right[i]}"; done)
      fi
    elif [[ ${protocol_rna_seq} == "mixed" ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        rnaspades_input=$(for i in "${!files_single[@]}"; do j=$((i + 1)); echo -n " --s${j} ${files_single[i]}"; done)
      elif [[ ${lib_layout} == 'paired' ]]; then
        rnaspades_input=$(for i in "${!files_left[@]}"; do j=$((i + 1)); echo -n " --pe${j}-1 ${files_left[i]} --pe${j}-2 ${files_right[i]}"; done)
      fi
    else
      echo "Invalid value for 'protocol_rna_seq'. Please specify either 'same' or 'mixed'."
      echo "Exiting."
      exit 1
    fi
    if [[ -d "${dir_tmp}/rnaspades_output" ]]; then
      rm -r "${dir_tmp}/rnaspades_output"
    fi
    rnaspades.py \
    --threads ${nslots_assembly} \
    --memory ${memory_assembly} \
    -o rnaspades_output \
    ${rnaspades_input}
    if [[ -s "${dir_tmp}/rnaspades_output/transcripts.fasta" ]]; then
      seqkit seq --threads ${NSLOTS} "${dir_tmp}/rnaspades_output/transcripts.fasta" --out-file "${file_isoform}"
    fi
  else
    echo "Invalid value for 'assembly_method'. Please specify either 'Trinity' or 'rnaSPAdes'."
    echo "Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='Longest CDS extraction'
disable_if_no_input_file "run_longestcds" ${file_isoform}
if [[ ( ! -s ${file_longestcds} || ! -s ${file_longestcds_transcript} ) && ${run_longestcds} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_parent_dir "${file_longestcds}"
  ensure_parent_dir "${file_longestcds_transcript}"

  if [[ ${orf_aggregation_level} = "p" ]]; then
    aggregate_expression="\.${orf_aggregation_level}[0-9].*"
  else
    aggregate_expression="\-${orf_aggregation_level}[0-9].*"
  fi
  if [[ ${orf_aggregation_level} == "c" && ${assembly_method} == "rnaSPAdes" ]]; then
    echo "The aggregation level 'c' is not supported for rnaSPAdes assemblies. Please set 'orf_aggregation_level' to either 'i' or 'g'."
    echo "Exiting"
    exit 1
  fi

  echo "scientific_name: ${sp_ub}"
  if [[ ${assembly_method} == 'Trinity' ]]; then
    seqkit sort --by-length --reverse --threads ${NSLOTS} "${file_isoform}" \
    | sed -e "s/_/-/g" -e "s/^>[[:space:]]*/>${sp_ub}_/" -e "s/TRINITY-//" -e "s/[[:space:]].*//" \
    > "${sp_ub}.tmp.renamed.fa"
  elif [[ ${assembly_method} == 'rnaSPAdes' ]]; then
    seqkit sort --by-length --reverse --threads ${NSLOTS} "${file_isoform}" \
    | sed -e "s/^>.*_g\([0-9]\+\)_i\([0-9]\+\)$/>${sp_ub}_g\1-i\2/" \
    > "${sp_ub}.tmp.renamed.fa"
  else
    echo "Invalid value for 'assembly_method'. Please specify either 'Trinity' or 'rnaSPAdes'."
    echo "Exiting."
    exit 1
  fi

  TransDecoder.LongOrfs -G universal -m 50 -t "${sp_ub}.tmp.renamed.fa"

  seqkit rmdup --by-seq "${sp_ub}.tmp.renamed.fa.transdecoder_dir/longest_orfs.cds" \
  | cdskit pad \
  | seqkit sort --by-name --threads ${NSLOTS} \
  | cdskit aggregate -x "${aggregate_expression}" \
  > "${sp_ub}.tmp.aggregated.fa"

  grep -e "^>" "${sp_ub}.tmp.aggregated.fa" | sed -e "s/>//" -e "s/[[:space:]].*//" -e "s/\.p[0-9].*//" \
  > "${sp_ub}.tmp.aggregated.transcript_id.txt"

  if [[ -s "${sp_ub}.tmp.aggregated.fa" && -s "${sp_ub}.tmp.aggregated.transcript_id.txt" ]]; then
    echo "Output file detected for the task: ${task}"
    seqkit grep \
    --threads ${NSLOTS} \
    --pattern-file "${sp_ub}.tmp.aggregated.transcript_id.txt" \
    "${sp_ub}.tmp.renamed.fa" \
    --out-file "${file_longestcds_transcript}"
    sed -e "s/[[:space:]].*//" -e "s/${aggregate_expression}//" "${sp_ub}.tmp.aggregated.fa" \
    | seqkit seq --threads ${NSLOTS} --out-file "${file_longestcds}"
  else
    echo "Output file not detected for the task: ${task}"
  fi
  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task="seqkit fx2tab for the longest CDS sequences"
disable_if_no_input_file "run_longestcds_fx2tab" ${file_longestcds}
if [[ ! -s ${file_longestcds_fx2tab} && ${run_longestcds_fx2tab} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  seqkit fx2tab \
  --threads ${NSLOTS} \
  --length \
  --name \
  --gc \
  --gc-skew \
  --header-line \
  --only-id \
  ${file_longestcds} \
  > "tmp.cds_length.tsv"

  if [[ -s "tmp.cds_length.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "tmp.cds_length.tsv" ${file_longestcds_fx2tab}
  fi
else
  echo "$(date): Skipped: ${task}"
fi

task="MMseqs2 Taxonomy of the CDS sequences"
disable_if_no_input_file "run_longestcds_mmseqs2taxonomy" ${file_longestcds}
if [[ ! -s ${file_longestcds_mmseqs2taxonomy} && ${run_longestcds_mmseqs2taxonomy} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  ensure_mmseqs_uniref90_db "${dir_mmseqs2_db}" "${NSLOTS}"
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare MMseqs2 UniRef90 DB. Exiting."
    exit 1
  fi

  if [[ ! -s "queryDB" ]]; then
    mmseqs createdb ${file_longestcds} queryDB
  fi

  if [[ ! -e "tmp_mmseqs2" ]]; then
    mkdir -p "tmp_mmseqs2"
  fi

  mmseqs taxonomy "queryDB" "${dir_mmseqs2_db}/UniRef90_DB" "output_prefix" "tmp" \
  --split-mode 2 \
  --split-memory-limit $((${MEM_PER_HOST}*3/4))G \
  --majority 0.5 \
  --lca-mode 3 \
  --vote-mode 1 \
  --tax-lineage 2 \
  --orf-filter 1 \
  --threads ${NSLOTS}

  mmseqs createtsv "queryDB" "output_prefix" "result.tsv" --threads ${NSLOTS}

  if [[ -s "result.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    mv_out "result.tsv" ${file_longestcds_mmseqs2taxonomy}
    rm queryDB*
    rm output_prefix.*
    rm -r "tmp_mmseqs2"
  fi
else
  echo "$(date): Skipped: ${task}"
fi

task="Contaminated sequence removal from the CDS sequences"
disable_if_no_input_file "run_longestcds_contamination_removal" ${file_longestcds} ${file_longestcds_fx2tab} ${file_longestcds_mmseqs2taxonomy}
if [[ ( ! -s ${file_longestcds_contamination_removal_fasta} || ! -s ${file_longestcds_contamination_removal_tsv} ) && ${run_longestcds_contamination_removal} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_parent_dir "${file_longestcds_contamination_removal_fasta}"

  if ! ensure_ete_taxonomy_db "${dir_pg}"; then
    echo "Failed to prepare ETE taxonomy DB. Exiting."
    exit 1
  fi

  python ${dir_myscript}/remove_contaminated_sequences.py \
  --fasta_file ${file_longestcds} \
  --mmseqs2taxonomy_tsv ${file_longestcds_mmseqs2taxonomy} \
  --fx2tab_tsv ${file_longestcds_fx2tab} \
  --species_name ${sp_ub} \
  --rank ${contamination_removal_rank} \
  --ncpu ${NSLOTS} \
  --rename_seq "no" \
  --verbose "no"

  if [[ -s "clean_sequences.fa" && -s "lineage_compatibility.tsv" ]]; then
    echo "Output file detected for the task: ${task}"
    seqkit seq --threads ${NSLOTS} "clean_sequences.fa" --out-file ${file_longestcds_contamination_removal_fasta}
    rm "clean_sequences.fa"
    mv_out "lineage_compatibility.tsv" ${file_longestcds_contamination_removal_tsv}
  fi
else
  echo "$(date): Skipped: ${task}"
fi

task='BUSCO for cDNA isoforms (isoform.fasta)'
disable_if_no_input_file "run_busco1" ${file_isoform}
if [[ ( ! -s ${file_busco_full1} || ! -s ${file_busco_short1} ) && ${run_busco1} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  seqkit seq --threads ${NSLOTS} ${file_isoform} --out-file "busco_infile_cdna.fa"

  dir_busco_db="/usr/local/db/busco_downloads"
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  if [[ ! -e ${dir_busco_lineage} ]]; then
    dir_busco_db=${dir_pg_db}/busco_downloads
    dir_busco_lineage=${dir_busco_db}/lineages/${busco_lineage}

    mkdir -p ${dir_busco_db}
    flock ${dir_busco_db}/${busco_lineage}.lock -c "
      if [ ! -e '${dir_busco_db}/lineages/${busco_lineage}' ]; then
        echo 'Starting BUSCO dataset download.'
        busco --download ${busco_lineage}
        mv busco_downloads/* ${dir_busco_db}
        rm -r busco_downloads
        echo 'BUSCO dataset download has been finished.'
      fi
    "
  fi

  busco \
  --in "busco_infile_cdna.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu ${NSLOTS} \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset ${dir_busco_lineage} \
  --download_path ${dir_busco_db} \
  --offline

  cp_out ./busco_tmp/run_${busco_lineage}/full_table.tsv ${file_busco_full1}
  cp_out ./busco_tmp/run_${busco_lineage}/short_summary.txt ${file_busco_short1}
  rm -r './busco_tmp'
  rm "busco_infile_cdna.fa"

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='BUSCO for longest CDS'
disable_if_no_input_file "run_busco2" ${file_longestcds}
if [[ ( ! -s ${file_busco_full2} || ! -s ${file_busco_short2} ) && ${run_busco2} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  seqkit seq --threads ${NSLOTS} ${file_longestcds} --out-file "busco_infile_cds.fa"

  dir_busco_db="/usr/local/db/busco_downloads"
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  if [[ ! -e ${dir_busco_lineage} ]]; then
    dir_busco_db=${dir_pg_db}/busco_downloads
    dir_busco_lineage=${dir_busco_db}/lineages/${busco_lineage}

    mkdir -p ${dir_busco_db}
    flock ${dir_busco_db}/${busco_lineage}.lock -c "
      if [ ! -e '${dir_busco_db}/lineages/${busco_lineage}' ]; then
        echo 'Starting BUSCO dataset download.'
        busco --download ${busco_lineage}
        mv busco_downloads/* ${dir_busco_db}
        rm -r busco_downloads
        echo 'BUSCO dataset download has been finished.'
      fi
    "
  fi

  busco \
  --in "busco_infile_cds.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu ${NSLOTS} \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset ${dir_busco_lineage} \
  --download_path ${dir_busco_db} \
  --offline

  cp_out ./busco_tmp/run_${busco_lineage}/full_table.tsv ${file_busco_full2}
  cp_out ./busco_tmp/run_${busco_lineage}/short_summary.txt ${file_busco_short2}
  rm -r './busco_tmp'
  rm "busco_infile_cds.fa"

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='BUSCO for contamination-removed longest CDS'
disable_if_no_input_file "run_busco3" ${file_longestcds_contamination_removal_fasta}
if [[ ( ! -s ${file_busco_full3} || ! -s ${file_busco_short3} ) && ${run_busco3} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  seqkit seq --threads ${NSLOTS} ${file_longestcds_contamination_removal_fasta} --out-file "busco_infile_cds.fa"

  dir_busco_db="/usr/local/db/busco_downloads"
  dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

  if [[ ! -e ${dir_busco_lineage} ]]; then
    dir_busco_db=${dir_pg_db}/busco_downloads
    dir_busco_lineage=${dir_busco_db}/lineages/${busco_lineage}

    mkdir -p ${dir_busco_db}
    flock ${dir_busco_db}/${busco_lineage}.lock -c "
      if [ ! -e '${dir_busco_db}/lineages/${busco_lineage}' ]; then
        echo 'Starting BUSCO dataset download.'
        busco --download ${busco_lineage}
        mv busco_downloads/* ${dir_busco_db}
        rm -r busco_downloads
        echo 'BUSCO dataset download has been finished.'
      fi
    "
  fi

  busco \
  --in "busco_infile_cds.fa" \
  --mode transcriptome \
  --out "busco_tmp" \
  --cpu ${NSLOTS} \
  --force \
  --evalue 1e-03 \
  --limit 20 \
  --lineage_dataset ${dir_busco_lineage} \
  --download_path ${dir_busco_db} \
  --offline

  cp_out ./busco_tmp/run_${busco_lineage}/full_table.tsv ${file_busco_full3}
  cp_out ./busco_tmp/run_${busco_lineage}/short_summary.txt ${file_busco_short3}
  rm -r './busco_tmp'
  rm "busco_infile_cds.fa"

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='Assembly statistics'
disable_if_no_input_file "run_assembly_stat" ${file_isoform}
if [[ ${run_longestcds} -eq 1 ]]; then disable_if_no_input_file "run_assembly_stat" ${file_longestcds}; fi
if [[ ${run_longestcds_contamination_removal} -eq 1 ]]; then disable_if_no_input_file "run_assembly_stat" ${file_longestcds_contamination_removal_fasta}; fi
if [[ ! -s ${file_assembly_stat} && ${run_assembly_stat} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  input_files=${file_isoform}
  if [[ -s ${file_longestcds} ]]; then
    input_files="${input_files} ${file_longestcds}"
  fi
  if [[ -s ${file_longestcds_contamination_removal_fasta} ]]; then
    input_files="${input_files} ${file_longestcds_contamination_removal_fasta}"
  fi

  seqkit stats \
  --all \
  --tabular \
  --threads ${NSLOTS} \
  --out-file assembly_stat.tsv \
  ${input_files}

  if [[ -s assembly_stat.tsv ]]; then
    echo "The task ${task} was successful."
    mv_out assembly_stat.tsv ${file_assembly_stat}
  fi
  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='amalgkit quant'
disable_if_no_input_file "run_amalgkit_quant" ${file_amalgkit_metadata}
if [[ ( ! -s ${file_amalgkit_merge_efflen} || ! -s ${file_amalgkit_merge_count} || ! -s ${file_amalgkit_merge_tpm} ) && ${run_amalgkit_quant} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_dir "${dir_amalgkit_quant}/${sp_ub}"

  if [[ -e "./getfastq" ]]; then
    rm -r "./getfastq"
  fi
  ln -s ${dir_amalgkit_getfastq_sp} "./getfastq"
  recreate_dir "./fasta"
  recreate_dir "./index"
  recreate_dir "./quant"

  path_reference_fasta_link="./fasta/${sp_ub}_for_kallisto_index.fasta"
  if [[ ${kallisto_reference} == 'species_cds' ]]; then
    path_kallisto_reference_fasta=$(ls ${dir_pg_input}/species_cds/${sp_ub}_*)
  elif [[ ${kallisto_reference} == 'longest_transcript' ]]; then
    path_kallisto_reference_fasta=${file_longestcds_transcript}
  elif [[ ${kallisto_reference} == 'longest_cds' ]]; then
    path_kallisto_reference_fasta=${file_longestcds}
  elif [[ ${kallisto_reference} == 'contamination_removed_longest_cds' ]]; then
    path_kallisto_reference_fasta=${file_longestcds_contamination_removal_fasta}
  else
    echo "Please check the input parameter. kallisto_reference must not be: ${kallisto_reference}"
    exit 1
  fi

  echo "kallisto reference = ${kallisto_reference}: ${path_kallisto_reference_fasta}"
  if [[ -e ${path_kallisto_reference_fasta} ]]; then
    ln -s ${path_kallisto_reference_fasta} ${path_reference_fasta_link}
  else
    echo "kallisto reference fasta file was not found in: ${path_kallisto_reference_fasta}"
    exit 1
  fi

  amalgkit quant \
   --out_dir "./" \
  --threads ${NSLOTS} \
  --metadata "${file_amalgkit_metadata}" \
  --clean_fastq no \
  --fasta_dir "./fasta" \
  --build_index yes
  exit_code_amalgkit_quant=$?

  if [[ ${exit_code_amalgkit_quant} -ne 0 ]]; then
    echo "amalgkit quant failed with exit code ${exit_code_amalgkit_quant}"
    exit 1
  else
    echo "amalgkit quant finished successfully"
    if [[ ! -e "${dir_amalgkit_quant}/${sp_ub}" ]]; then
      mkdir -p "${dir_amalgkit_quant}/${sp_ub}"
    fi
    mv ./quant/* "${dir_amalgkit_quant}/${sp_ub}"
    rm -r "./quant"
    rm "./getfastq" # Do not put -r, otherwise the original getfastq files will be deleted.
  fi

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='amalgkit merge'
disable_if_no_input_file "run_amalgkit_merge" ${file_amalgkit_metadata}
if [[ ( ! -s ${file_amalgkit_merge_efflen} || ! -s ${file_amalgkit_merge_count} || ! -s ${file_amalgkit_merge_tpm} || ! -s ${file_amalgkit_merge_metadata} ) && ${run_amalgkit_merge} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)

  if [[ -e ./quant ]]; then
    rm -r ./quant
  fi
  ln -s ${dir_amalgkit_quant}/${sp_ub} ./quant

  amalgkit merge \
  --out_dir "./" \
  --metadata "${file_amalgkit_metadata}"

  #sp_metadata=$(python -c "import pandas; d=pandas.read_csv('${file_amalgkit_metadata}',sep='\t',header=0); print(d.at[0,'scientific_name'].replace(' ','_'))")
  if [[ -s "./merge/${sp_ub}/${sp_ub}_eff_length.tsv" ]]; then
    echo "Copying amalgkit merge outputs from: ./merge/${sp_ub}"
    mv "./merge/${sp_ub}" "${dir_amalgkit_merge}"
    mv_out "./merge/metadata.tsv" "${file_amalgkit_merge_metadata}"
    rm -r "./merge"
    rm "./quant" # Do not put -r, otherwise the original quant files will be deleted.
  else
    echo "amalgkit merge outputs were not found in: ./merge/${sp_ub}"
  fi
  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

task='Multispecies summary'
is_output_older_than_inputs "^file_" ${file_multispecies_summary}; summary_flag=$?
if [[ ${run_multispecies_summary} -eq 1 && ${summary_flag} -eq 1 ]]; then
  echo "$(date): Start: ${task}" | tee >(cat >&2)
  ensure_dir "${dir_multispecies_summary}"
  cd ${dir_multispecies_summary}

  python ${dir_myscript}/collect_common_BUSCO_genes.py \
  --busco_outdir ${dir_busco_full2} \
  --ncpu ${NSLOTS} \
  --outfile busco_table.tsv

  for dir_busco in ${dir_busco_full1} ${dir_busco_full2} ${dir_busco_full3}; do
    if [[ -z "$(find "${dir_busco}" -mindepth 1 -print -quit)" ]]; then
      echo "Skipping. No BUSCO output was found in: ${dir_busco}"
      continue
    fi
    Rscript ${dir_myscript}/annotation_summary.r \
    --dir_species_cds_busco="${dir_busco}" \
    --tree_annotation_dir="${dir_myscript}/tree_annotation" \
    --min_og_species='auto'
    mv "annotation_summary.tsv" "$(basename ${dir_busco}).tsv"
    mv "busco_cds.svg" "$(basename ${dir_busco}).svg"
    mv "busco_cds.pdf" "$(basename ${dir_busco}).pdf"
  done

  Rscript ${dir_myscript}/multispecies_transcriptome_summary.r \
  --dir_assembly_stat=${dir_assembly_stat} \
  --dir_amalgkit_metadata=${dir_amalgkit_metadata} \
  --dir_amalgkit_merge=${dir_amalgkit_merge} \
  --dir_busco_isoform=${dir_busco_full1} \
  --dir_busco_longest_cds=${dir_busco_full2} \
  --dir_myscript=${dir_myscript}

  if [[ -e "Rplots.pdf" ]]; then
    rm "Rplots.pdf"
  fi
  cd ${dir_tmp}

  echo "$(date): End: ${task}"
else
  echo "$(date): Skipped: ${task}"
fi

if [[ ${remove_amalgkit_fastq_after_completion} -eq 1 && $(is_fastq_requiring_downstream_analysis_done) -eq 1 ]]; then
  echo "remove_amalgkit_fastq_after_completion=1: All necessary output files were detected. amalgkit getfastq outputs will be removed."
  if [[ -e ${dir_amalgkit_getfastq_sp} ]]; then
    rm -r ${dir_amalgkit_getfastq_sp}
    ensure_parent_dir "${file_amalgkit_getfastq_safely_removed_flag}"
    echo "Fastq files for this species have been safely removed." > ${file_amalgkit_getfastq_safely_removed_flag}
  fi
else
  echo "fastp fastq files will not be removed."
fi

remove_empty_subdirs ${dir_transcriptome_assembly_output}
if [[ ${delete_tmp_dir} -eq 1 && $(is_fastq_requiring_downstream_analysis_done) -eq 1 ]]; then
  echo "delete_tmp_dir=1: All necessary output files were detected. Deleting ${dir_tmp}"
  rm -r ${dir_tmp}
else
  echo "Tmp directory will not be deleted: ${dir_tmp}"
fi

echo "$(date): Exiting Singularity environment"
