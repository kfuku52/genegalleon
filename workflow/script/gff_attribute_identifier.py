#! /usr/bin/env python

import argparse
import datetime
import gzip
import pandas
import re
import sys
import time

pandas.options.mode.chained_assignment = None

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--gff_file', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--gff_feature', metavar='STR', default='mRNA', type=str, help='')
    return parser

def read_fasta_seqname(file_path):
    raw_seqnames = []
    if file_path.endswith('.gz'):
        with gzip.open(file_path,'r') as f:
            for line in f.readlines():
                line = line.decode("utf-8")
                if ">" in line:
                    raw_seqnames.append(line)
    else:
        with open(file_path,'r') as f:
            for line in f.readlines():
                if ">" in line:
                    raw_seqnames.append(line)
    seqnames = pandas.Series(raw_seqnames)
    seqnames = seqnames.str.replace('^>','',regex=True)
    seqnames = seqnames.str.replace('\n$','',regex=True)
    return seqnames

def get_gene_names(seq_names):
    gene_names = seq_names.str.replace('_', 'PLACEHOLDERTEXT', 2)
    gene_names = gene_names.str.replace(r'.*PLACEHOLDERTEXT', '', regex=True)
    return gene_names

def filename2sciname(file_name):
    tmp = re.sub('.*/', '', file_name)
    tmp = tmp.replace('_','|',1)
    tmp = re.sub(r'[-_.,].*', '', tmp)
    sci_name = tmp.replace('|', '_')
    return sci_name

def identify_attribute_key(args):
    seq_names = read_fasta_seqname(file_path=args.fasta_file)
    sci_name = filename2sciname(file_name=args.fasta_file)
    seq_names = seq_names.str.replace('^'+sci_name+'_', '', regex=True).sort_values().reset_index(drop=True)
    gff_cols = ['sequence','source','feature','start','end','score','strand','phase','attributes']
    gff = pandas.read_csv(args.gff_file, sep='\t', header=None, comment='#', low_memory=False, quoting=3)
    gff.columns = gff_cols
    gff = gff.loc[(gff['feature']==args.gff_feature),:].reset_index(drop=True)
    for i in gff.index:
        attrs = {re.sub('=.*', '', item): re.sub('.*=', '', item) for item in gff.loc[i,'attributes'].split(';')}
        for key, value in attrs.items():
            if (value == seq_names).any():
                return key
    return None


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    start_time = time.time()
    sys.stderr.write('gff_attribute_identifier.py started: {}\n'.format(datetime.datetime.now()))
    attr_key = identify_attribute_key(args)
    if attr_key is not None:
        print(attr_key)
        txt = 'gff_attribute_identifier.py completed in {:,} secs: {}\n'
        sys.stderr.write(txt.format(int(time.time()-start_time), datetime.datetime.now()))
        sys.exit(0)
    txt = 'gff_attribute_identifier.py failed. Time elapsed: {:,} secs. {}\n'
    sys.stderr.write(txt.format(int(time.time()-start_time), datetime.datetime.now()))
    sys.exit(1)


if __name__ == '__main__':
    main()
