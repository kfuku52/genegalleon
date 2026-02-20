#!/usr/bin/env python
# coding: utf-8

import argparse
import math
import re

import pandas
from kftools.kfphylo import *
from kftools.kfseq import *


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--iqtree', metavar='PATH', type=str, required=True, help='.iqtree output from IQ-TREE')
    parser.add_argument('--log', metavar='PATH', type=str, required=True, help='.log output from IQ-TREE')
    parser.add_argument('--state', metavar='PATH', type=str, required=True, help='.state output from IQ-TREE')
    parser.add_argument('--alignment', metavar='PATH', type=str, required=True, help='input alignment for IQ-TREE')
    parser.add_argument('--treefile', metavar='PATH', type=str, required=True, help='.treefile output from IQ-TREE')
    parser.add_argument('--rooted_tree', metavar='PATH', type=str, required=True, help='A rooted newick tree')
    parser.add_argument('--genetic_code', metavar='INT', type=int, required=True, help='NCBI genetic code')
    return parser


VALID_NUCLEOTIDES = ('A', 'C', 'G', 'T')


def parse_codon_freqs_from_iqtree(iqtree_path):
    # IQ-TREE2 writes codon pi(...) values, but IQ-TREE3 may not.
    pattern = re.compile(r'pi\(([ACGT]{3})\)\s*=\s*([0-9eE+\-\.]+)')
    codon_freqs = {}
    with open(iqtree_path, 'r') as file:
        for line in file:
            for codon, value in pattern.findall(line):
                try:
                    codon_freqs[codon] = float(value)
                except ValueError:
                    continue
    return codon_freqs


def safe_codon_freqs_to_nuc_freqs(codon_freqs, model):
    if codon_freqs is None:
        return None
    cleaned = {}
    for codon, value in codon_freqs.items():
        if re.fullmatch(r'[ACGT]{3}', str(codon)) is None:
            continue
        try:
            freq = float(value)
        except (TypeError, ValueError):
            continue
        if (not math.isfinite(freq)) or (freq < 0):
            continue
        cleaned[codon] = freq
    if (len(cleaned) == 0) or (sum(cleaned.values()) <= 0):
        return None
    try:
        return codon2nuc_freqs(codon_freqs=cleaned, model=model)
    except ZeroDivisionError:
        return None


def normalize_nuc_counts(counts):
    total = sum(counts.values())
    if total <= 0:
        return {nuc: 0.25 for nuc in VALID_NUCLEOTIDES}
    return {nuc: counts[nuc] / total for nuc in VALID_NUCLEOTIDES}


def read_fasta_sequences(alignment_file):
    sequences = {}
    current_header = None
    chunks = []
    with open(alignment_file, 'r') as file:
        for raw_line in file:
            line = raw_line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(chunks)
                current_header = line[1:].split()[0]
                chunks = []
            else:
                if current_header is None:
                    continue
                chunks.append(line)
    if current_header is not None:
        sequences[current_header] = ''.join(chunks)
    return sequences


def normalize_label(label):
    return re.sub(r'[^A-Za-z0-9]', '_', str(label))


def alignment_subset_nuc_freqs(alignment_file, model, leaf_names=None):
    if 'F3X4' not in model:
        raise ValueError("alignment_subset_nuc_freqs supports only F3X4 models.")

    all_sequences = read_fasta_sequences(alignment_file)
    if len(all_sequences) == 0:
        raise ValueError(f'No sequences were found in alignment: {alignment_file}')

    if leaf_names is None:
        selected_headers = list(all_sequences.keys())
    else:
        leaf_names = [str(name) for name in leaf_names if str(name) != '']
        leaf_set = set(leaf_names)
        leaf_norm_set = {normalize_label(name) for name in leaf_names}
        selected_headers = []
        for header in all_sequences.keys():
            if header in leaf_set:
                selected_headers.append(header)
                continue
            if any(header.startswith(name) for name in leaf_names):
                selected_headers.append(header)
                continue
            if normalize_label(header) in leaf_norm_set:
                selected_headers.append(header)
    if len(selected_headers) == 0:
        raise ValueError('No matching sequences were found in alignment for the target leaf set.')

    counts = [{nuc: 0.0 for nuc in VALID_NUCLEOTIDES} for _ in range(3)]
    for header in selected_headers:
        seq = all_sequences[header].upper().replace('U', 'T')
        usable_len = len(seq) - (len(seq) % 3)
        for i in range(0, usable_len, 3):
            codon = seq[i:(i + 3)]
            if any(base not in VALID_NUCLEOTIDES for base in codon):
                continue
            for pos, base in enumerate(codon):
                counts[pos][base] += 1.0
    return [normalize_nuc_counts(c) for c in counts]


def run(args):
    print('args:', vars(args))

    with open(args.log, 'r') as file:
        text_log = file.read().split('\n')
    line_omega = [t for t in text_log if 'Nonsynonymous/synonymous ratio' in t][0]
    line_kappa = [t for t in text_log if 'Transition/transversion ratio' in t][0]
    line_alpha = [t for t in text_log if 'Gamma shape alpha' in t][0]
    line_command = [t for t in text_log if 'Command: ' in t][0]
    value_omega = float(re.sub(r'.*: ', '', line_omega))
    value_kappa = float(re.sub(r'.*: ', '', line_kappa))
    value_alpha = float(re.sub(r'.*: ', '', line_alpha))
    value_model = re.sub(r' .*', '', re.sub(r'.*-m ', '', line_command))
    print('omega =', value_omega)
    print('kappa =', value_kappa)
    print('alpha =', value_alpha)
    print('model =', value_model)

    codon_freqs = parse_codon_freqs_from_iqtree(args.iqtree)
    nuc_freqs = safe_codon_freqs_to_nuc_freqs(codon_freqs=codon_freqs, model=value_model)
    if nuc_freqs is None:
        print('No valid codon pi(...) frequencies were found. Using alignment-derived nucleotide frequencies.')
        nuc_freqs = alignment_subset_nuc_freqs(alignment_file=args.alignment, model=value_model)
    print('equilibrium nucleotide frequency:', nuc_freqs)

    thetas = nuc_freq2theta(nuc_freqs=nuc_freqs)
    print('equilibrium theta:', thetas)

    iqtree_tree = transfer_root(tree_to=args.treefile, tree_from=args.rooted_tree)

    anc_state = pandas.read_csv(args.state, sep='\t', comment='#')
    subroot_nodes = list(iqtree_tree.get_children())

    subroot_states = {}
    subroot_codon_freqs = {}
    subroot_nuc_freqs = {}
    subroot_thetas = {}
    for subroot_node in subroot_nodes:
        snn = subroot_node.name
        subroot_state = anc_state.loc[(anc_state['Node'] == snn), :]
        subroot_states[snn] = subroot_state

        subroot_nuc_freq = None
        if subroot_state.size != 0:
            codon_columns = subroot_state.columns[subroot_state.columns.str.startswith('p_')]
            subroot_codon_freq = subroot_state.loc[:, codon_columns].mean(axis=0)
            subroot_codon_freq.index = codon_columns.str.replace('p_', '')
            subroot_codon_freq = subroot_codon_freq.to_dict()
            subroot_codon_freqs[snn] = subroot_codon_freq
            subroot_nuc_freq = safe_codon_freqs_to_nuc_freqs(codon_freqs=subroot_codon_freq, model=value_model)
        if subroot_nuc_freq is None:
            print('Using alignment-derived subroot nucleotide frequencies for:', snn)
            subroot_nuc_freq = alignment_subset_nuc_freqs(
                alignment_file=args.alignment,
                model=value_model,
                leaf_names=list(subroot_node.leaf_names()),
            )
        subroot_nuc_freqs[snn] = subroot_nuc_freq
        subroot_thetas[snn] = nuc_freq2theta(nuc_freqs=subroot_nuc_freqs[snn])

    print('subroot nucleotide frequency:', subroot_nuc_freqs)
    root_thetas = weighted_mean_root_thetas(subroot_thetas, iqtree_tree, model=value_model)
    print('root theta:', root_thetas)

    num_gamma_category = re.sub(r'.*\+G', '', value_model)

    out = ''
    out += 'alphabet=Codon(letter=DNA)' + '\n'
    out += 'genetic_code=' + str(args.genetic_code) + '\n'
    out += 'input.sequence.file=$(SEQ)' + '\n'
    out += 'input.sequence.format=Fasta' + '\n'
    out += 'input.sequence.remove_stop_codons=yes' + '\n'
    out += 'input.tree.file=$(TREE)' + '\n'
    out += 'input.tree.format=Newick' + '\n'
    out += 'model=YN98(kappa=' + str(value_kappa) + ',omega=' + str(value_omega) + ',initFreqs=observed)' + '\n'
    out += 'rate_distribution=Gamma(n=' + num_gamma_category + ',alpha=' + str(value_alpha) + ')' + '\n'
    out += 'map.type=DnDs' + '\n'
    out += 'output.counts=PerBranch(prefix=$(OUT).)' + '\n'
    out += 'output.tree_with_id.file=$(OUT).with_id.nwk' + '\n'
    print(out)

    with open('iqtree2mapnh.params', 'w') as file:
        file.write(out)
    iqtree_tree.write(outfile='iqtree2mapnh.nwk', parser=5)
    print('iqtree2mapnh, done!')


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
