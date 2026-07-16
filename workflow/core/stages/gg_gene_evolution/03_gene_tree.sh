# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="IQ-TREE"
disable_if_no_input_file "run_iqtree" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_iqtree_tree}" && ${run_iqtree} -eq 1 ]]; then
  gg_step_start "${task}"
  num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_seq} -ge 4 ]]; then
    if [[ ${run_generax} -eq 1 ]]; then
      other_iqtree_params=()
      file_tree="${og_id}.treefile"
      use_ufboot=0
      echo "run_generax=1: disabling UFBOOT in the initial IQ-TREE run. Support will be computed on the GeneRax topology."
    else
      other_iqtree_params=(--ufboot 1000 --bnni)
      file_tree="${og_id}.contree"
      use_ufboot=1
    fi
  else
    other_iqtree_params=()
    file_tree="${og_id}.treefile"
    use_ufboot=0
  fi
  if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
    if [[ ${use_ufboot} -eq 1 ]]; then
      echo "Disabling IQ-TREE UFBOOT because fast mode is enabled for large alignments (${num_seq} > ${iqtree_fast_mode_gt})."
      other_iqtree_params=()
      file_tree="${og_id}.treefile"
      use_ufboot=0
    fi
    other_iqtree_params+=(--fast)
  fi

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    assert_gene_evolution_aa_model_for_protein_mode "${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
  elif [[ ${run_generax} -eq 1 ]]; then
    if gene_evolution_model_is_aa "${generax_model}"; then
      echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
        sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
          > iqtree_input.fa
    else
      echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
      seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
    fi
  else
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
  fi
  echo "IQ-TREE starting..."

  iqtree_model_string=${generax_model}
  build_iqtree_mem_args

  iqtree \
    -s iqtree_input.fa \
    -m "${iqtree_model_string}" \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
    --prefix "${og_id}" \
    "${IQTREE_MEM_ARGS[@]}" \
    --seed 12345 \
    --redo \
    "${other_iqtree_params[@]}"

  cp_out "${file_tree}" "${file_og_iqtree_tree}"
else
  gg_step_skip "${task}"
fi

task="Gene tree rooting"
disable_if_no_input_file "run_tree_root" "${file_og_unrooted_tree_analysis}"
if [[ (! -s "${file_og_rooted_tree}" || ! -s "${file_og_rooted_log}") && ${run_tree_root} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${tree_rooting_method}" == "notung" ]]; then
    if [[ ! -s "${species_tree_pruned}" ]]; then
      echo "tree_rooting_method=notung requires species tree: ${species_tree_pruned}"
      exit 1
    fi
    if [[ -e "./${og_id}.notung.root" ]]; then
      rm -rf -- "./${og_id}.notung.root"
    fi

    echo "memory_notung: ${memory_notung}"
    java -jar -Xmx${memory_notung}g "${notung_jar}" \
      -s "${species_tree_pruned}" \
      -g "${file_og_unrooted_tree_analysis}" \
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
      --outputdir "./${og_id}.notung.root"

    rooted_candidates=()
    mapfile -t rooted_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "${og_id}.iqtree.nwk.rooting.*" | sort -V)
    selected_rooted_tree=""
    for candidate in "${rooted_candidates[@]}"; do
      if [[ "${candidate}" =~ \.rooting\.[0-9]+$ ]]; then
        selected_rooted_tree="${candidate}"
        break
      fi
    done
    if [[ -z "${selected_rooted_tree}" ]]; then
      echo "NOTUNG did not generate rooted-tree candidates in ./${og_id}.notung.root"
      exit 1
    fi

    nwkit label --target intnode --force yes --infile "${selected_rooted_tree}" --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" "${file_og_rooted_tree}"
    {
      echo "tree_rooting_method=notung"
      echo "selected_rooting=${selected_rooted_tree}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" "${file_og_rooted_log}"
  else
    nwkit_root_method="${tree_rooting_method}"
    if [[ "${nwkit_root_method}" == "md" ]]; then
      nwkit_root_method="mv"
    fi
    nwkit_root_args=(root --method "${nwkit_root_method}" --infile "${file_og_unrooted_tree_analysis}")
    if [[ "${nwkit_root_method}" == "taxonomy" ]]; then
      nwkit_root_args+=(--species-parser "${species_label_parser}")
      if [[ -n "${species_label_regex}" ]]; then
        nwkit_root_args+=(--species-regex "${species_label_regex}")
      fi
      if [[ -n "${species_label_map_tsv}" ]]; then
        nwkit_root_args+=(--species-map-tsv "${species_label_map_tsv}")
      fi
    fi
    nwkit "${nwkit_root_args[@]}" |
      nwkit label --target intnode --force yes --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" "${file_og_rooted_tree}"
    {
      echo "tree_rooting_method=${tree_rooting_method}"
      echo "nwkit_method=${nwkit_root_method}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" "${file_og_rooted_log}"
  fi
else
  gg_step_skip "${task}"
fi

task="Orthogroup extraction with NWKIT"
run_orthogroup_extraction_original="${run_orthogroup_extraction}" # This variable may be disabled by disable_if_no_input_file but the original value is necessary to properly update file_og_*_analysis
disable_if_no_input_file "run_orthogroup_extraction" "${file_query_gene:-}" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
if [[ (! -s "${file_og_orthogroup_extraction_nwk}" || ! -s "${file_og_orthogroup_extraction_rooted_nwk}" || ! -s "${file_og_orthogroup_extraction_fasta}") && ${run_orthogroup_extraction} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
    echo "Fasta format was detected. Query IDs absent from the tree will be replaced by their best tree-backed query BLAST hit."
  else
    echo "Gene IDs were detected."
    cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
  fi

  run_nwkit_subtree() {
    local infile=$1
    echo "Running nwkit subtree for ${infile}"
    local info_txt
    info_txt=$(nwkit subtree --infile "${infile}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup_conf_score_threshold 0 2> /dev/null | nwkit info 2> /dev/null)
    local num_leaf
    num_leaf=$(awk -F': *' '/Number of leaves/ {print $2; exit}' <<< "${info_txt}")
    printf '%s\t%s\n' "${num_leaf}" "${infile}" >> tmp_num_leaf.tsv
  }

  subtree_infiles=()
  if [[ "${tree_rooting_method}" == "notung" && -d "./${og_id}.notung.root" ]]; then
    mapfile -t subtree_infiles < <(
      find "./${og_id}.notung.root" -maxdepth 1 -type f |
        awk -v og="${og_id}" '$0 ~ (og "\\.iqtree\\.nwk\\.rooting\\.[0-9]+$") {print}' |
        sort -V
    )
  fi
  if [[ ${#subtree_infiles[@]} -eq 0 ]]; then
    if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
      subtree_infiles=("${file_og_rooted_tree_analysis}")
    else
      echo "No rooted tree is available for orthogroup extraction."
      exit 1
    fi
  fi
  comma_separated_genes=$(
    python - "${file_query_gene}" "${file_og_query_blast:-}" "${subtree_infiles[0]}" << 'PY'
import csv
import math
import os
import sys

query_gene_path = sys.argv[1]
query_blast_path = sys.argv[2]
tree_path = sys.argv[3]


def normalize_id(value):
    return value.strip().replace("−", "-")


def parse_query_ids(path):
    query_ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first_char = handle.read(1)
        handle.seek(0)
        if first_char == ">":
            for line in handle:
                if not line.startswith(">"):
                    continue
                query_id = normalize_id(line[1:].split()[0])
                if query_id:
                    query_ids.append(query_id)
        else:
            for line in handle:
                query_id = normalize_id(line)
                if query_id:
                    query_ids.append(query_id)
    return query_ids


def parse_newick_leaves(path):
    text = open(path, "r", encoding="utf-8", errors="replace").read()
    leaves = set()
    index = 0
    while index < len(text):
        if text[index] not in "(,":
            index += 1
            continue
        index += 1
        while index < len(text) and text[index].isspace():
            index += 1
        if index >= len(text) or text[index] == ")":
            continue
        if text[index] == "'":
            index += 1
            token = []
            while index < len(text):
                char = text[index]
                if char == "'":
                    if index + 1 < len(text) and text[index + 1] == "'":
                        token.append("'")
                        index += 2
                        continue
                    index += 1
                    break
                token.append(char)
                index += 1
            label = "".join(token)
        else:
            start = index
            while index < len(text) and text[index] not in ":,();":
                index += 1
            label = text[start:index]
        label = normalize_id(label)
        if label:
            leaves.add(label)
    return leaves


def to_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def best_tree_backed_hits(path, leaves):
    best = {}
    if not path or not os.path.exists(path):
        return best
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_index, row in enumerate(reader):
            qacc = normalize_id(row.get("qacc", ""))
            sacc = normalize_id(row.get("sacc", ""))
            if not qacc or not sacc or sacc not in leaves:
                continue
            evalue = to_float(row.get("min_evalue") or row.get("evalue"), math.inf)
            bitscore = to_float(row.get("max_bitscore") or row.get("bitscore"), -math.inf)
            rank = (evalue, -bitscore, row_index)
            current = best.get(qacc)
            if current is None or rank < current["rank"]:
                best[qacc] = {
                    "sacc": sacc,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "rank": rank,
                }
    return best


query_ids = parse_query_ids(query_gene_path)
tree_leaves = parse_newick_leaves(tree_path)
best_hits = best_tree_backed_hits(query_blast_path, tree_leaves)
resolved = []
seen = set()

for query_id in query_ids:
    if query_id in tree_leaves:
        seed = query_id
    elif query_id in best_hits:
        hit = best_hits[query_id]
        seed = hit["sacc"]
        print(
            "Orthogroup extraction seed fallback: "
            f"query {query_id} was not found in the tree; using best tree-backed "
            f"query BLAST hit {seed} (evalue={hit['evalue']}, bitscore={hit['bitscore']}).",
            file=sys.stderr,
        )
    else:
        print(
            "Warning: Orthogroup extraction seed skipped because neither the query "
            f"nor a tree-backed query BLAST hit was found: {query_id}",
            file=sys.stderr,
        )
        continue
    if seed in seen:
        print(
            f"Orthogroup extraction seed duplicate skipped: {seed} "
            f"(from query {query_id}).",
            file=sys.stderr,
        )
        continue
    seen.add(seed)
    resolved.append(seed)

if not resolved:
    print(
        "No orthogroup extraction seed genes were found in the rooted tree. Exiting.",
        file=sys.stderr,
    )
    sys.exit(1)

print(",".join(resolved))
PY
  )
  echo "Seed genes for orthogroup extraction: ${comma_separated_genes}"

  printf 'num_leaf\tfile\n' > tmp_num_leaf.tsv
  for subtree_infile in "${subtree_infiles[@]}"; do
    wait_until_jobn_le ${GG_TASK_CPUS}
    run_nwkit_subtree "${subtree_infile}"
  done

  if ! IFS=$'\t' read -r min_leaf_num min_leaf_file max_leaf_num max_leaf_file < <(
    awk -F'\t' '
      NR==1 {next}
      NR==2 {min=$1; max=$1; min_file=$2; max_file=$2; next}
      {
        if ($1 < min) {min=$1; min_file=$2}
        if ($1 > max) {max=$1; max_file=$2}
      }
      END {
        if (NR < 2) exit 1
        printf "%s\t%s\t%s\t%s\n", min, min_file, max, max_file
      }
    ' tmp_num_leaf.tsv
  ); then
    echo "Failed to parse tmp_num_leaf.tsv."
    exit 1
  fi
  echo "Minimum number of orthogroup subtree leaves after checking all rooting positions: ${min_leaf_num} in ${min_leaf_file} (will be used for orthogroup extraction)"
  echo "Maximum number of orthogroup subtree leaves after checking all rooting positions: ${max_leaf_num} in ${max_leaf_file} (shown just as a reference)"

  nwkit subtree --infile "${min_leaf_file}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup_conf_score_threshold 0 \
    --outfile "${og_id}.orthogroup_seed.tmp.nwk"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file tmp.trimmed.input.fasta
  nwkit intersection \
    --infile "${og_id}.orthogroup_seed.tmp.nwk" \
    --outfile /dev/null \
    --seqin tmp.trimmed.input.fasta \
    --seqout tmp.fasta \
    --match "complete"
  rm -f -- tmp.trimmed.input.fasta

  # Preserve IQ-TREE support values in the extracted unrooted tree.
  nwkit intersection \
    --infile "${file_og_iqtree_tree}" \
    --outfile "${og_id}.orthogroup_extraction.tmp.nwk" \
    --seqin tmp.fasta \
    --seqout /dev/null \
    --match "complete"
  mv_out "${og_id}.orthogroup_extraction.tmp.nwk" "${file_og_orthogroup_extraction_nwk}"
  mv_out "${og_id}.orthogroup_seed.tmp.nwk" "${file_og_orthogroup_extraction_rooted_nwk}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" tmp.fasta --out-file "${og_id}.orthogroup_extraction.tmp.fasta"
  else
    cdskit hammer --nail 4 -s tmp.fasta -o "${og_id}.orthogroup_extraction.tmp.fasta"
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.orthogroup_extraction.tmp.fasta" --out-file "${og_id}.orthogroup_extraction.out.fa.gz"
  mv_out "${og_id}.orthogroup_extraction.out.fa.gz" "${file_og_orthogroup_extraction_fasta}"
  rm -f -- "${og_id}.orthogroup_extraction.tmp.fasta"
else
  gg_step_skip "${task}"
fi
if [[ ${run_orthogroup_extraction_original} -eq 1 && -s "${file_og_orthogroup_extraction_nwk}" && -s "${file_og_orthogroup_extraction_rooted_nwk}" && -s "${file_og_orthogroup_extraction_fasta}" ]]; then
  set_analysis_file unrooted_tree "${file_og_orthogroup_extraction_nwk}"
  set_analysis_file rooted_tree "${file_og_orthogroup_extraction_rooted_nwk}"
  set_analysis_file trimmed_aln "${file_og_orthogroup_extraction_fasta}"
fi

task="GeneRax"
disable_if_no_input_file "run_generax" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${species_tree_pruned}"
if [[ ! -s "${file_og_generax_nhx}" && ${run_generax} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    assert_gene_evolution_aa_model_for_protein_mode "${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
  elif gene_evolution_model_is_aa "${generax_model}"; then
    echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
    seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > generax_input_alignment.fas
  else
    echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
  fi

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_unrooted_tree_analysis}" \
    --outfile generax_input_gene_tree.nwk

  #avoid multifurcating tree
  R -q -e "library(ape); t=read.tree(\"generax_input_gene_tree.nwk\"); t=multi2di(t,random=FALSE); write.tree(t, \"generax_input_gene_tree_bi.nwk\")"

  generate_generax_mapfile() {
    # https://github.com/BenoitMorel/GeneRax/wiki/Gene-to-species-mapping
    my_aln_file=$1
    awk '/^>/ {sub(/^>/, "", $0); print}' "${my_aln_file}" > tmp.gene_names.txt
    while IFS= read -r gene_name; do gg_species_name_from_path "${gene_name}"; done < tmp.gene_names.txt > tmp.species_names.txt
    paste tmp.gene_names.txt tmp.species_names.txt > generax_map.txt
    rm -f -- tmp.gene_names.txt tmp.species_names.txt
  }
  generate_generax_mapfile generax_input_alignment.fas

  printf '%s\n' \
    '[FAMILIES]' \
    '- family_1' \
    'starting_gene_tree = generax_input_gene_tree_bi.nwk' \
    'alignment = generax_input_alignment.fas' \
    'mapping = generax_map.txt' \
    "subst_model = ${generax_model}" \
    > generax_families.txt

  # GeneRax runs within one scheduler task; keep OpenMPI from probing ssh/rsh inside containers.
  mpiexec_args=(mpiexec -oversubscribe -np "${GG_TASK_CPUS}")
  mpi_env_args=(env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)
  if [[ "$(id -u)" -eq 0 ]]; then
    mpiexec_args+=(--allow-run-as-root)
  fi
  "${mpi_env_args[@]}" "${mpiexec_args[@]}" generax \
    --species-tree "${species_tree_pruned}" \
    --families generax_families.txt \
    --strategy "SPR" \
    --rec-model "${generax_rec_model}" \
    --prefix "generax_${og_id}" \
    --per-family-rates \
    --skip-family-filtering \
    --mad-rooting \
    --seed 12345 < /dev/null

  echo "GeneRax exit code = $?"

  generax_out_sptree="./generax_${og_id}/species_trees/starting_species_tree.newick" # generax v2.0
  if [[ -s "${generax_out_sptree}" ]]; then
    lock_file="${species_tree_generax}.lock"
    (
      if ! gg_shared_lock_acquire "${lock_file}" "GeneRax species tree copy"; then
        exit 1
      fi
      gg_shared_lock_start_heartbeat "${lock_file}"
      heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
      cleanup_generax_tree_lock() {
        gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
        gg_shared_lock_release "${lock_file}"
      }
      trap cleanup_generax_tree_lock EXIT
      if [[ ! -s "${species_tree_generax}" ]]; then
        echo "copying GeneRax output species tree (first writer only)."
        cp_out "${generax_out_sptree}" "${species_tree_generax}"
      fi
    ) || exit 1
  elif [[ ! -s "${species_tree_generax}" ]]; then
    echo "GeneRax species tree file was not found yet: ${generax_out_sptree}"
  fi
  echo "copying GeneRax output gene tree."
  reconciled_base="./generax_${og_id}/reconciliations/family_1_reconciliated"
  reconciled_xml="${reconciled_base}.xml"
  reconciled_nhx="${reconciled_base}.nhx"
  if [[ -e "${reconciled_nhx}" ]]; then
    echo "GeneRax outfile was found. Copying."
    nwkit nhx2nwk --infile "${reconciled_nhx}" --outfile "${og_id}.generax.tmp.nwk"
    mv_out "${og_id}.generax.tmp.nwk" "${file_og_generax_nwk}"
    cp_out "${reconciled_xml}" "${file_og_generax_xml}"
    cp_out "${reconciled_nhx}" "${file_og_generax_nhx}"
  else
    echo "GeneRax outfile was not found. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE UFBOOT on GeneRax topology"
if [[ ${run_generax} -eq 1 ]]; then
  if [[ ! -s "${file_og_generax_nwk}" ]]; then
    echo "Skipped: ${task}. Missing GeneRax output tree: ${file_og_generax_nwk}"
  elif [[ ! -s "${file_og_trimmed_aln_analysis}" ]]; then
    echo "Skipped: ${task}. Missing alignment: ${file_og_trimmed_aln_analysis}"
  elif [[ ! -s "${file_og_iqtree_generax_ufboot}" || "${file_og_generax_nwk}" -nt "${file_og_iqtree_generax_ufboot}" || "${file_og_trimmed_aln_analysis}" -nt "${file_og_iqtree_generax_ufboot}" ]]; then
    gg_step_start "${task}"
    num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
    if [[ ${num_seq} -lt 4 ]]; then
      echo "UFBOOT requires >=4 sequences. Using the GeneRax topology without bootstrap support."
      nwkit drop --target root --length yes --infile "${file_og_generax_nwk}" --outfile "${og_id}.generax_ufboot.tmp.nwk"
      mv_out "${og_id}.generax_ufboot.tmp.nwk" "${file_og_iqtree_generax_ufboot}"
    else
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        assert_gene_evolution_aa_model_for_protein_mode "${task}"
        seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.generax_ufboot.input.fa"
      elif gene_evolution_model_is_aa "${generax_model}"; then
        echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
        seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
          sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
            > "${og_id}.generax_ufboot.input.fa"
      else
        echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
        seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.generax_ufboot.input.fa"
      fi

      nwkit drop --target intnode --support yes --name yes --infile "${file_og_generax_nwk}" |
        nwkit drop --target all --length yes --outformat 9 --outfile "${og_id}.generax_ufboot.constraint.nwk"

      other_iqtree_params=(--ufboot 1000 --bnni)
      if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
        echo "Skipping IQ-TREE --fast in UFBOOT-on-GeneRax mode because the options are incompatible."
      fi

      build_iqtree_mem_args
      iqtree \
        -s "${og_id}.generax_ufboot.input.fa" \
        -g "${og_id}.generax_ufboot.constraint.nwk" \
        -m "${generax_model}" \
        -T AUTO \
        --threads-max "${GG_TASK_CPUS}" \
        --prefix "${og_id}.generax_ufboot" \
        "${IQTREE_MEM_ARGS[@]}" \
        --seed 12345 \
        --redo \
        "${other_iqtree_params[@]}"

      if [[ -s "${og_id}.generax_ufboot.contree" ]]; then
        cp_out "${og_id}.generax_ufboot.contree" "${file_og_iqtree_generax_ufboot}"
      elif [[ -s "${og_id}.generax_ufboot.treefile" ]]; then
        cp_out "${og_id}.generax_ufboot.treefile" "${file_og_iqtree_generax_ufboot}"
      else
        echo "IQ-TREE UFBOOT on GeneRax topology failed to generate a tree file."
        exit 1
      fi
      rm -f -- "${og_id}.generax_ufboot.input.fa" "${og_id}.generax_ufboot.constraint.nwk"
    fi
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi
if [[ ${run_generax} -eq 1 && -s "${file_og_iqtree_generax_ufboot}" ]]; then
  set_analysis_file unrooted_tree "${file_og_iqtree_generax_ufboot}"
fi

task="NOTUNG reconciliation"
disable_if_no_input_file "run_notung_reconcil" "${file_og_rooted_tree}" "${species_tree_pruned}"
if [[ ! -s "${file_og_notung_reconcil}" && ${run_notung_reconcil} -eq 1 ]]; then
  gg_step_start "${task}"

  echo "memory_notung: ${memory_notung}"

  if [[ -s "./${og_id}.root.nwk" ]]; then
    rm -f -- "${og_id}.root.nwk"
  fi
  if [[ -e "./${og_id}.notung_reconcile" ]]; then
    rm -rf -- "${og_id}.notung_reconcile"
  fi

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_rooted_tree}" \
    --outfile "${og_id}.root.nwk"

  java -jar -Xmx${memory_notung}g "${notung_jar}" \
    -s "${species_tree_pruned}" \
    -g "${og_id}.root.nwk" \
    --reconcile \
    --infertransfers "false" \
    --treeoutput newick \
    --log \
    --treestats \
    --events \
    --parsable \
    --speciestag prefix \
    --maxtrees 1 \
    --nolosses \
    --outputdir ./${og_id}.notung_reconcile

  if [[ -s "${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.parsable.txt" || -s "${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.0.parsable.txt" ]]; then
    zip -rq "${og_id}.notung_reconcile.zip" "${og_id}.notung_reconcile"
    cp_out "${og_id}.notung_reconcile.zip" "${file_og_notung_reconcil}"
  fi
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided divergence time estimation"
disable_if_no_input_file "run_tree_dating" "${species_tree_pruned}" "${file_og_unrooted_tree_analysis}"
if [[ (! -s "${file_og_dated_tree}" || ! -s "${file_og_dated_tree_log}") && ${run_tree_dating} -eq 1 ]]; then
  gg_step_start "${task}"
  radte_args=()

  if [[ ${run_generax} -eq 1 ]]; then
    radte_args+=("--species_tree=${species_tree_generax}")
    radte_args+=("--generax_nhx=${file_og_generax_nhx}")
  else
    if [[ -e "./${og_id}.notung_reconcile" ]]; then
      rm -rf -- "./${og_id}.notung_reconcile"
    fi
    cp_out "${file_og_notung_reconcil}" .
    unzip -q "$(basename "${file_og_notung_reconcil}")"
    if [[ -s ./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.0 ]]; then
      cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0 ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled
      cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0.parsable.txt ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.parsable.txt
    fi
    radte_args+=("--species_tree=${species_tree_pruned}")
    radte_args+=("--gene_tree=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled")
    radte_args+=("--notung_parsable=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.parsable.txt")
  fi
  radte_args+=("--species-parser=${species_label_parser}")
  if [[ -n "${species_label_regex}" ]]; then
    radte_args+=("--species-regex=${species_label_regex}")
  fi
  if [[ -n "${species_label_map_tsv}" ]]; then
    radte_args+=("--species-map-tsv=${species_label_map_tsv}")
  fi

  radte.r \
    "${radte_args[@]}" \
    --max_age="${radte_max_age}" \
    --chronos_lambda=1 \
    --chronos_model=discrete \
    --pad_short_edge=0.001 \
    2>&1 | tee radte.log

  constrained_node=$(awk -F': *' '/^Calibrated nodes:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' radte.log)
  echo "${constrained_node}" > "${og_id}.dated.log.txt"

  if grep -q ":-" radte_gene_tree_output.nwk; then
    contain_negative_bl=1
  else
    contain_negative_bl=0
  fi
  if [[ ${contain_negative_bl} -eq 1 ]]; then
    echo "Dated tree has negative branch length. Deleting output files depending on the tree file."
    for key in l1ou pem scm dated stat tree_plot; do
      files=()
      mapfile -t files < <(compgen -A variable "file_og_${key}")
      for f in "${files[@]}"; do
        target_file="${!f}"
        if [[ -e "${target_file}" ]]; then
          echo "deleting: ${target_file}"
          rm -f -- "${target_file}"
        fi
      done
    done
  else
    echo "Dated tree has no negative branch length. Continue."
    cp_out radte_calibrated_nodes.txt "${file_og_dated_tree_log}"
    cp_out radte_gene_tree_output.nwk "${file_og_dated_tree}"
  fi
else
  gg_step_skip "${task}"
fi
