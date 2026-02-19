#!/bin/bash

# Check if input files, output files, and max length are provided
if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <input.fasta> <output.fasta> <input.newick> <output.newick> <max_length>"
  exit 1
fi

input_fasta=$1
output_fasta=$2
input_newick=$3
output_newick=$4
max_length=$5

# Create temporary files
temp_fasta=$(mktemp)
temp_newick=$(mktemp)

# Ensure long_names.txt exists even if it's empty
touch long_names.txt

# Step 1 & 2: Extract long names and generate shortened names from FASTA file
awk -v max_length="$max_length" '/^>/ {
    seq_name = substr($0, 2);
    if (length($0)-1 > max_length) {
        print seq_name > "long_names.txt";
        print "Detected long sequence name: " seq_name | "cat 1>&2";
    }
}' "$input_fasta"

# Create a file with old name and new name pairs
awk -v max_length="$max_length" '{
    new_name = substr($0, 1, max_length - length(NR)) "_" NR;
    print $0, new_name;
}' long_names.txt > name_pairs.txt

# Function to replace names in a file
replace_names() {
    local input_file=$1
    local output_file=$2
    local temp_file=$(mktemp)

    cp "$input_file" "$temp_file"
    while read old_name new_name; do
        # Check if file is Newick format (based on file extension)
        if [[ "$input_file" == *.nwk ]]; then
            # Use Perl for more sophisticated regex in Newick file
            perl -i -pe "s/(?<=\(|,)$old_name(?=:|\)|,)/$new_name/g" "$temp_file"
        else
            # Use Perl for precise replacement in FASTA file, ensuring whole-word matching
            perl -i -pe "s/^>$old_name(?=\$| )/>$new_name/g" "$temp_file"
        fi
    done < name_pairs.txt

    if [[ "$input_file" == "$output_file" ]]; then
        mv "$temp_file" "$output_file"
    else
        cp "$temp_file" "$output_file"
        rm "$temp_file"
    fi
}

# Step 3: Replace in FASTA file
replace_names "$input_fasta" "$output_fasta"

# Step 4: Replace in Newick file
replace_names "$input_newick" "$output_newick"

# Optionally, clean up other temporary files
[ -e long_names.txt ] && rm long_names.txt
[ -e name_pairs.txt ] && rm name_pairs.txt
