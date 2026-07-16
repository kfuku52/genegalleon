# shellcheck shell=bash
# Sourced by gg_transcriptome_generation_core.sh.

task='De novo transcriptome assembly'
if [[ ! -s "${file_isoform}" && ${run_assembly} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_isoform}"
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    echo "amalgkit getfastq output directory not found: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi
  if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2> /dev/null)" ]]; then
    echo "No amalgkit getfastq FASTQ files were found in: ${dir_amalgkit_getfastq_sp}. Exiting."
    exit 1
  fi

  assembly_cpus=$((${GG_TASK_CPUS} - ${assembly_cpu_offset}))
  if ! gg_is_nonnegative_integer "${assembly_ram_offset:-0}"; then
    echo "Invalid assembly_ram_offset=${assembly_ram_offset:-unset}; using 0."
    assembly_ram_offset=0
  fi
  assembly_ram_reserved_gb=${assembly_ram_offset}
  if [[ "${assembly_ram_reserved_gb}" -lt "${GG_MEM_TOOL_RESERVE_GB}" ]]; then
    assembly_ram_reserved_gb=${GG_MEM_TOOL_RESERVE_GB}
  fi
  assembly_mem_gb=$((${GG_MEM_TOTAL_GB} - ${assembly_ram_reserved_gb}))
  if [[ ${assembly_cpus} -lt 1 ]]; then
    echo "Adjusted assembly_cpus from ${assembly_cpus} to 1 because it must be >=1."
    assembly_cpus=1
  fi
  if [[ ${assembly_mem_gb} -lt 1 ]]; then
    echo "Adjusted assembly_mem_gb from ${assembly_mem_gb}G to 1G because it must be >=1G."
    assembly_mem_gb=1
  fi
  bflyHeapSpaceMax=$(gg_memory_fraction_gb "${assembly_mem_gb}" 1 "${assembly_cpus}") # GB
  echo "Number of offset CPUs for transcriptome assembly is ${assembly_cpu_offset}."
  echo "Number of reserved GB for transcriptome assembly is ${assembly_ram_reserved_gb}."
  echo "${GG_TASK_CPUS} CPUs and ${GG_MEM_TOTAL_GB}G RAM are allocated to this job."
  echo "${assembly_cpus} CPUs and ${assembly_mem_gb}G RAM are used for transcriptome assembly."

  if [[ "${effective_assembly_method}" == 'rna-bloom2' ]]; then
    files_long=()
    load_classified_getfastq_files "${file_amalgkit_read_technology}"
    files_long=("${classified_long_fastq_files[@]}")
    if [[ ${#files_long[@]} -eq 0 ]]; then
      echo "No long-read FASTQ files were detected for RNA-Bloom2 assembly. Exiting."
      exit 1
    fi

    total_fastq_len=$(get_total_fastq_len_from_files "${files_long[@]}")
    max_assembly_input_fastq_size="${max_assembly_input_fastq_size//,/}"
    if [[ ${total_fastq_len} -gt ${max_assembly_input_fastq_size} ]]; then
      echo "Total long-read fastq length is ${total_fastq_len} bp, which is greater than ${max_assembly_input_fastq_size} bp."
      echo "Only ${max_assembly_input_fastq_size} bp will be used."
      assembly_input_fastq_dir="./tmp_assembly_input_fastq"
      if [[ -e "${assembly_input_fastq_dir}" ]]; then
        rm -rf -- "${assembly_input_fastq_dir}"
      fi
      mkdir -p "${assembly_input_fastq_dir}"
      proportion=$(awk -v max="${max_assembly_input_fastq_size}" -v total="${total_fastq_len}" 'BEGIN {printf "%.3f\n", max/total}')
      echo "Proportion of long-read fastq reads to be used: ${proportion} (${max_assembly_input_fastq_size}/${total_fastq_len})"
      for file in "${files_long[@]}"; do
        seqkit sample --proportion "${proportion}" --rand-seed 11 --out-file "${assembly_input_fastq_dir}/$(basename "${file}")" "${file}"
      done
      mapfile -t files_long < <(find "${assembly_input_fastq_dir}" -type f -name "*.amalgkit.fastq.gz" | sort)
      echo "Total fastq length of the subsampled long-read fastq files: $(get_total_fastq_len_from_files "${files_long[@]}") bp"
    else
      echo "Total long-read fastq length is ${total_fastq_len} bp, which is less than ${max_assembly_input_fastq_size} bp. All long-read fastq reads will be used."
    fi

    rnabloom_input_args=(-long)
    rnabloom_input_args+=("${files_long[@]}")
    rnabloom_extra_args=()
    if [[ ${detected_has_pacbio} -eq 1 ]]; then
      rnabloom_extra_args+=(-lrpb)
    fi
    if [[ ${detected_has_ont_direct_rna} -eq 1 ]]; then
      rnabloom_extra_args+=(-stranded)
    fi
    if [[ -d "${dir_tmp}/rnabloom_output" ]]; then
      rm -rf -- "${dir_tmp}/rnabloom_output"
    fi
    rnabloom \
      "${rnabloom_extra_args[@]}" \
      "${rnabloom_input_args[@]}" \
      -t "${assembly_cpus}" \
      -mem "${assembly_mem_gb}" \
      -outdir rnabloom_output
    if [[ -s "${dir_tmp}/rnabloom_output/rnabloom.transcripts.fa" ]]; then
      seqkit seq --threads "${GG_TASK_CPUS}" "${dir_tmp}/rnabloom_output/rnabloom.transcripts.fa" --out-file "tmp.isoform.fa.gz"
      mv_out "tmp.isoform.fa.gz" "${file_isoform}"
    else
      echo "RNA-Bloom2 did not produce rnabloom.transcripts.fa in: ${dir_tmp}/rnabloom_output. Exiting."
      exit 1
    fi
  else
    mapfile -t files_right < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*_2.amalgkit.fastq.gz" | sort)
    if [[ ${#files_right[@]} -eq 0 ]]; then
      echo "Paired-end samples were not detected. SE reads will be used for transcriptome assembly."
      lib_layout='single'
    else
      echo "Paired-end samples were detected. PE reads will be used for transcriptome assembly and SE reads, if any, will not be used."
      lib_layout='paired'
    fi

    selected_fastq_dir="./tmp_selected_fastq"
    if [[ -e "${selected_fastq_dir}" ]]; then
      rm -rf -- "${selected_fastq_dir}"
    fi
    mkdir -p "${selected_fastq_dir}"
    if [[ "${effective_assembly_method}" == "rnaspades" && ${protocol_rna_seq} == "mixed" ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" | sort)
      elif [[ ${lib_layout} == 'paired' ]]; then
        mapfile -t candidate_fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
      fi
      n_pair=${#candidate_fastq_files[@]}
      if [[ ${n_pair} -ge 10 ]]; then
        echo "Selecting top 9 ${lib_layout}-end fastq files by size to fit SPAdes/rnaSPAdes numbered library input slots (--s1..--s9 or --pe1..--pe9)."
        echo "For details, please refer to the rnaSPAdes manual: https://ablab.github.io/spades/rna.html"
        echo "Selected fastq files:"
        if [[ ${lib_layout} == 'single' ]]; then
          find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -exec du -b {} + |
            sort -nr | head -n 9 | cut -f2 | while read -r file; do
            cp_out "${file}" "${selected_fastq_dir}"
            echo "${file}"
          done
        elif [[ ${lib_layout} == 'paired' ]]; then
          find "${dir_amalgkit_getfastq_sp}" -type f -name "*_1.amalgkit.fastq.gz" -exec du -b {} + |
            sort -nr | head -n 9 | cut -f2 | while read -r file1; do
            file2="${file1/_1.amalgkit.fastq.gz/_2.amalgkit.fastq.gz}"
            cp_out "${file1}" "${selected_fastq_dir}"
            echo "${file1}"
            cp_out "${file2}" "${selected_fastq_dir}"
            echo "${file2}"
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
      total_fastq_len=$(get_total_fastq_len "${selected_fastq_dir}" "*.amalgkit.fastq.gz")
    elif [[ ${lib_layout} == 'paired' ]]; then
      total_fastq_len1=$(get_total_fastq_len "${selected_fastq_dir}" "*_1.amalgkit.fastq.gz")
      total_fastq_len2=$(get_total_fastq_len "${selected_fastq_dir}" "*_2.amalgkit.fastq.gz")
      total_fastq_len=$((${total_fastq_len1} + ${total_fastq_len2}))
    fi
    max_assembly_input_fastq_size="${max_assembly_input_fastq_size//,/}"
    if [[ ${total_fastq_len} -gt ${max_assembly_input_fastq_size} ]]; then
      echo "Total ${lib_layout} fastq length is ${total_fastq_len} bp, which is greater than ${max_assembly_input_fastq_size} bp."
      echo "Only ${max_assembly_input_fastq_size} bp will be used."
      assembly_input_fastq_dir="./tmp_assembly_input_fastq"
      if [[ -e "${assembly_input_fastq_dir}" ]]; then
        rm -rf -- "${assembly_input_fastq_dir}"
      fi
      mkdir -p "${assembly_input_fastq_dir}"
      proportion=$(awk -v max="${max_assembly_input_fastq_size}" -v total="${total_fastq_len}" 'BEGIN {printf "%.3f\n", max/total}')
      echo "Proportion of fastq reads to be used: ${proportion} (${max_assembly_input_fastq_size}/${total_fastq_len})"
      files=()
      if [[ ${lib_layout} == 'single' ]]; then
        mapfile -t files < <(find "${selected_fastq_dir}" -type f -name "*.amalgkit.fastq.gz" | sort)
      elif [[ ${lib_layout} == 'paired' ]]; then
        mapfile -t files1 < <(find "${selected_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
        mapfile -t files2 < <(find "${selected_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz" | sort)
        files=("${files1[@]}" "${files2[@]}")
      fi
      for file in "${files[@]}"; do
        seqkit sample --proportion "${proportion}" --rand-seed 11 --out-file "${assembly_input_fastq_dir}/$(basename "${file}")" "${file}"
      done
      echo "Total fastq length of the subsampled fastq files: $(get_total_fastq_len "${assembly_input_fastq_dir}" "*.amalgkit.fastq.gz") bp"
    else
      echo "Total fastq length is ${total_fastq_len} bp, which is less than ${max_assembly_input_fastq_size} bp. All fastq reads will be used."
      assembly_input_fastq_dir=${selected_fastq_dir}
    fi

    files_single=()
    files_left=()
    files_right=()
    mapfile -t files_single < <(find "${assembly_input_fastq_dir}" -type f -name "*.amalgkit.fastq.gz" | sort)
    mapfile -t files_left < <(find "${assembly_input_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz" | sort)
    mapfile -t files_right < <(find "${assembly_input_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz" | sort)

    if [[ ${lib_layout} == 'single' && ${#files_single[@]} -eq 0 ]]; then
      echo "No single-end fastq files were found for transcriptome assembly. Exiting."
      exit 1
    fi
    if [[ ${lib_layout} == 'paired' ]]; then
      if [[ ${#files_left[@]} -eq 0 || ${#files_right[@]} -eq 0 ]]; then
        echo "Paired-end input files were not detected correctly in ${assembly_input_fastq_dir}."
        echo "Detected left/right counts: ${#files_left[@]}/${#files_right[@]}. Exiting."
        exit 1
      fi
      if ! filter_valid_paired_fastq_files "${dir_tmp}/paired_fastq_validation.tsv"; then
        exit 1
      fi
    fi

    if [[ "${effective_assembly_method}" == 'trinity' ]]; then
      if [[ ${lib_layout} == 'single' ]]; then
        in_single="$(
          IFS=","
          echo "${files_single[*]}"
        )"
      elif [[ ${lib_layout} == 'paired' ]]; then
        in_left="$(
          IFS=","
          echo "${files_left[*]}"
        )"
        in_right="$(
          IFS=","
          echo "${files_right[*]}"
        )"
      fi
      if [[ ${lib_layout} == 'single' ]]; then
        Trinity \
          --seqType fq \
          --CPU "${assembly_cpus}" \
          --max_memory "${assembly_mem_gb}G" \
          --min_contig_length 200 \
          --output trinity \
          --full_cleanup \
          --NO_SEQTK \
          --bflyHeapSpaceMax "${bflyHeapSpaceMax}G" \
          --bflyGCThreads 1 \
          --bflyCalculateCPU \
          --single "${in_single}"
      else
        Trinity \
          --seqType fq \
          --CPU "${assembly_cpus}" \
          --max_memory "${assembly_mem_gb}G" \
          --min_contig_length 200 \
          --output trinity \
          --full_cleanup \
          --NO_SEQTK \
          --bflyHeapSpaceMax "${bflyHeapSpaceMax}G" \
          --bflyGCThreads 1 \
          --bflyCalculateCPU \
          --left "${in_left}" \
          --right "${in_right}"
      fi
      # For --NO_SEQTK, see https://github.com/trinityrnaseq/trinityrnaseq/issues/787
      if [[ -s "${dir_tmp}/trinity.Trinity.fasta" ]]; then
        seqkit seq --threads "${GG_TASK_CPUS}" "${dir_tmp}/trinity.Trinity.fasta" --out-file "tmp.isoform.fa.gz"
        mv_out "tmp.isoform.fa.gz" "${file_isoform}"
      fi
    elif [[ "${effective_assembly_method}" == 'rnaspades' ]]; then
      rnaspades_transcript_fasta=""
      rnaspades_input_args=()
      if [[ ${protocol_rna_seq} == "same" ]]; then
        if [[ ${lib_layout} == 'single' ]]; then
          for i in "${!files_single[@]}"; do
            rnaspades_input_args+=(--s1 "${files_single[i]}")
          done
        elif [[ ${lib_layout} == 'paired' ]]; then
          for i in "${!files_left[@]}"; do
            rnaspades_input_args+=(--pe1-1 "${files_left[i]}" --pe1-2 "${files_right[i]}")
          done
        fi
      elif [[ ${protocol_rna_seq} == "mixed" ]]; then
        if [[ ${lib_layout} == 'single' ]]; then
          for i in "${!files_single[@]}"; do
            j=$((i + 1))
            rnaspades_input_args+=("--s${j}" "${files_single[i]}")
          done
        elif [[ ${lib_layout} == 'paired' ]]; then
          for i in "${!files_left[@]}"; do
            j=$((i + 1))
            rnaspades_input_args+=("--pe${j}-1" "${files_left[i]}" "--pe${j}-2" "${files_right[i]}")
          done
        fi
      else
        echo "Invalid value for 'protocol_rna_seq'. Please specify either 'same' or 'mixed'."
        echo "Exiting."
        exit 1
      fi
      if [[ -d "${dir_tmp}/rnaspades_output" ]]; then
        rm -rf -- "${dir_tmp}/rnaspades_output"
      fi
      OMP_NUM_THREADS="${assembly_cpus}" \
      OMP_THREAD_LIMIT="${assembly_cpus}" \
      rnaspades.py \
        --threads "${assembly_cpus}" \
        --memory "${assembly_mem_gb}" \
        -o rnaspades_output \
        "${rnaspades_input_args[@]}"
      if rnaspades_transcript_fasta=$(resolve_rnaspades_transcript_fasta "${dir_tmp}/rnaspades_output"); then
        echo "Using rnaSPAdes transcript fasta: ${rnaspades_transcript_fasta}"
        seqkit seq --threads "${GG_TASK_CPUS}" "${rnaspades_transcript_fasta}" --out-file "tmp.isoform.fa.gz"
        if [[ ! -s "tmp.isoform.fa.gz" ]]; then
          echo "Failed to generate tmp.isoform.fa.gz from rnaSPAdes output: ${rnaspades_transcript_fasta}. Exiting."
          exit 1
        fi
        mv_out "tmp.isoform.fa.gz" "${file_isoform}"
      else
        echo "rnaSPAdes did not produce a supported transcript fasta under: ${dir_tmp}/rnaspades_output."
        echo "Checked: transcripts.fasta, soft_filtered_transcripts.fasta, hard_filtered_transcripts.fasta"
        exit 1
      fi
    fi
  fi

  if [[ ! -s "${file_isoform}" ]]; then
    echo "Transcriptome assembly did not generate: ${file_isoform}. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task='Corset clustering of long-read transcripts'
if [[ "${effective_assembly_method}" == 'rna-bloom2' && -s "${file_isoform}" && ! -s "${file_corset_clusters}" && ${run_longestcds} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_corset_clusters}"
  ensure_parent_dir "${file_corset_counts}"
  load_classified_getfastq_files "${file_amalgkit_read_technology}"
  if [[ ${#classified_long_fastq_files[@]} -eq 0 ]]; then
    echo "No long-read FASTQ files were detected for Corset clustering. Exiting."
    exit 1
  fi

  recreate_dir "./corset_work"
  cd "./corset_work" || exit 1
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_isoform}" --out-file "assembly.fa"

  corset_bams=()
  corset_names=()
  corset_groups=()
  corset_minimap2_preset="map-ont"
  if [[ ${detected_has_pacbio} -eq 1 ]]; then
    corset_minimap2_preset="map-pb"
  fi
  for long_fastq in "${classified_long_fastq_files[@]}"; do
    sample_name="$(basename "${long_fastq}")"
    sample_name="${sample_name%.amalgkit.fastq.gz}"
    bam_file="${sample_name}.bam"
    minimap2 \
      -a \
      -x "${corset_minimap2_preset}" \
      -N 50 \
      --secondary=yes \
      -t "${GG_TASK_CPUS}" \
      "assembly.fa" \
      "${long_fastq}" \
      | samtools view -@ "${GG_TASK_CPUS}" -b -o "${bam_file}" -
    corset_bams+=( "${bam_file}" )
    corset_names+=( "${sample_name}" )
    corset_groups+=( "1" )
  done

  if [[ ${#corset_bams[@]} -eq 0 ]]; then
    echo "No BAM files were generated for Corset clustering. Exiting."
    exit 1
  fi

  corset \
    -D 99999999999 \
    -g "$(csv_join_from_array "${corset_groups[@]}")" \
    -n "$(csv_join_from_array "${corset_names[@]}")" \
    "${corset_bams[@]}"

  if [[ -s "clusters.txt" ]]; then
    mv_out "clusters.txt" "${file_corset_clusters}"
  else
    echo "Corset did not produce clusters.txt. Exiting."
    exit 1
  fi
  if [[ -s "counts.txt" ]]; then
    mv_out "counts.txt" "${file_corset_counts}"
  fi
  cd "${dir_tmp}" || exit 1
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi
