import argparse
import os
import re
import shlex
import subprocess
import sys

import numpy
import pandas


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_genome', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--geneinfo_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--seqkit_exe', metavar='PATH', default='seqkit', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='promoter_sequence.fasta', type=str, help='')
    parser.add_argument('--promoter_bp', metavar='INT', default=2000, type=int, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of seqkit threads.')
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    return parser


def print_edge_seq(df, i, mode):
    txt = 'Located at the edge of the reference chromosome/scaffold. '
    if mode == 'all':
        txt += 'No promoter sequence is available: {}, '
    elif mode == 'part':
        txt += 'The promoter sequence is truncated: {}, '
    txt += 'start = {:,}, end = {:,}, strand = {}'
    print(txt.format(df.at[i, 'gene_id'], df.at[i, 'start'], df.at[i, 'end'], df.at[i, 'strand']))


def get_genome_file(dir_genome, sp):
    genome_files = os.listdir(dir_genome)
    genome_files = [file for file in genome_files if file.startswith(sp)]
    if len(genome_files) == 0:
        print('Genome fasta not found. Skipping {}'.format(sp))
        return None
    if len(genome_files) == 1:
        print('Genome fasta for {} is {}'.format(sp, genome_files[0]))
        return os.path.join(dir_genome, genome_files[0])
    sys.stderr.write('Multiple genome fasta files found. Skipping {}. {}\n'.format(sp, ', '.join(genome_files)))
    return None


def run_pipeline(commands):
    processes = []
    prev_stdout = None
    for command in commands:
        proc = subprocess.Popen(command, stdin=prev_stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if prev_stdout is not None:
            prev_stdout.close()
        processes.append(proc)
        prev_stdout = proc.stdout
    stdout_final, stderr_final = processes[-1].communicate()
    stderr_texts = []
    return_codes = []
    for proc in processes[:-1]:
        proc.wait()
        return_codes.append(proc.returncode)
        stderr_texts.append(proc.stderr.read().decode('utf8'))
    return_codes.append(processes[-1].returncode)
    stderr_texts.append(stderr_final.decode('utf8'))
    return stdout_final, return_codes, stderr_texts


def run_seqkit(args, df, i, genome_file):
    strand = df.at[i, 'strand']
    if strand == '+':
        region_end = df.at[i, 'start'] - 1
        if region_end <= 0:
            print_edge_seq(df, i, mode='all')
            return ''
        region_start = region_end - args.promoter_bp + 1
        if region_start <= 0:
            print_edge_seq(df, i, mode='part')
            region_start = 1
    elif strand == '-':
        region_start = df.at[i, 'end'] + 1
        region_end = region_start + args.promoter_bp - 1
    else:
        sys.stderr.write('Invalid strand value for {}: {}. Skipping.\n'.format(df.at[i, 'gene_id'], strand))
        return ''

    seqkit_region = str(region_start) + ':' + str(region_end)
    command1 = [args.seqkit_exe, 'grep', '--pattern', df.at[i, 'chromosome'], '--threads', str(args.ncpu), genome_file]
    command2 = [args.seqkit_exe, 'subseq', '--region', seqkit_region, '--threads', str(args.ncpu)]
    command3 = [args.seqkit_exe, 'seq', '--reverse', '--complement', '--seq-type', 'DNA']
    cmd_to_str = lambda cmd: ' '.join([shlex.quote(str(x)) for x in cmd])
    if strand == '+':
        commands = [command1, command2]
    else:
        commands = [command1, command2, command3]
    print('Command:', ' | '.join(cmd_to_str(cmd) for cmd in commands))
    stdout_bytes, return_codes, stderr_texts = run_pipeline(commands)
    if any(code != 0 for code in return_codes):
        sys.stderr.write("seqkit did not finish safely.\n")
        for stderr_text in stderr_texts:
            if stderr_text != '':
                sys.stderr.write(stderr_text)
                sys.stderr.write('\n')
    fasta = stdout_bytes.decode('utf8')
    fasta = re.sub('>.*', '>' + df.at[i, 'gene_id'], fasta)
    return fasta


def sanitize_fasta(fasta_records):
    fasta_records = [f for f in fasta_records if f != '']
    new_fasta_records = list()
    for fasta_record in fasta_records:
        seq = re.sub('\n', '', re.sub('>.*', '', fasta_record))
        if len(seq) > 0:
            new_fasta_records.append(fasta_record)
        else:
            print('Empty sequence. Removed:')
            print(fasta_record)
    return new_fasta_records


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    print('get_promoter_fasta.py started.')

    df = pandas.read_csv(args.geneinfo_tsv, sep='\t', header=0)
    species_key = df.loc[:, 'gene_id'].astype(str).str.extract(r'^([^_]+_[^_]+)', expand=False)
    df.loc[:, '_species_key'] = species_key.fillna('')
    spp = sorted([sp for sp in df.loc[:, '_species_key'].unique() if sp != ''])
    species_indices = df.groupby('_species_key', sort=False).indices
    fasta_records = list()
    for sp in spp:
        print('\nHandling the reference genome of {}'.format(sp))
        sp_ind = species_indices.get(sp, [])
        genome_file = get_genome_file(dir_genome=args.dir_genome, sp=sp)
        if genome_file is None:
            continue
        for i in sp_ind:
            fasta_records.append(run_seqkit(args, df, i, genome_file))

    new_fasta_records = sanitize_fasta(fasta_records)
    print('Number of rows in input tsv: {:,}'.format(df.shape[0]))
    print('Number of sequences in output fasta: {:,}'.format(len(new_fasta_records)))

    with open(args.outfile, 'w') as f:
        f.write('\n'.join(new_fasta_records))
    print('get_promoter_fasta.py is completed.')


if __name__ == '__main__':
    main()
