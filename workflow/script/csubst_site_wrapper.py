#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import glob
import gzip
import matplotlib
import numpy
import os
import pandas
import re
import shutil
import sqlalchemy
import subprocess
import textwrap
import zipfile
from concurrent.futures import ProcessPoolExecutor

try:
    from distutils.util import strtobool
except ModuleNotFoundError:
    def strtobool(val):
        val = str(val).strip().lower()
        if val in ('y', 'yes', 't', 'true', 'on', '1'):
            return 1
        if val in ('n', 'no', 'f', 'false', 'off', '0'):
            return 0
        raise ValueError("invalid truth value {!r}".format(val))
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

pandas.set_option("display.max_columns", None)

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none' # none, path, or svgfont

def plot_scatter(cb, xcol, ycol, outbase, xmax=10, ymax=10, nbin=50, polygon_xmin=0, polygon_ymin=0):
    axis_margin_rate = 1.02
    cb_plot = cb.copy()
    cb_plot[xcol] = cb_plot[xcol].clip(0, xmax)
    cb_plot[ycol] = cb_plot[ycol].clip(0, ymax)
    wr = (0.8, 0.2)
    hr = (0.2, 0.8)
    fig, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, figsize=(3.2, 3.2), sharey=False, sharex=False,
                                           width_ratios=wr, height_ratios=hr)
    axes = axes.flat
    axes[1].axis('off')
    ax = axes[0]
    ax.axvline(x=min_OCNany2spe, color='red', alpha=0.3)
    ax.hist(x=cb_plot[xcol], bins=numpy.arange(0, 1.001 + (1 / nbin), 1 / nbin) * xmax, color='darkblue')
    ax.set_xlim(0, xmax * axis_margin_rate)
    ax.set_xticks([])
    ax.set_yscale('log')
    ax.set_ylabel('Count')
    ax = axes[3]
    ax.axhline(y=min_omegaCany2spe, color='red', alpha=0.3)
    ax.hist(x=cb_plot[ycol], bins=numpy.arange(0, 1.001 + (1 / nbin), 1 / nbin) * ymax, color='darkblue',
            orientation='horizontal')
    ax.set_ylim(0, ymax * axis_margin_rate)
    ax.set_yticks([])
    ax.set_xscale('log')
    ax.set_xlabel('Count')
    ax = axes[2]
    plot_range = [0, xmax * axis_margin_rate, 0, ymax * axis_margin_rate]
    polygon_coord = [[polygon_xmin, polygon_ymin], [xmax, polygon_ymin], [xmax, ymax],
                     [polygon_xmin, ymax]]
    ax.add_patch(matplotlib.patches.Polygon(polygon_coord, color='red', alpha=0.3))
    ax.hexbin(
        cb_plot[xcol].values,
        cb_plot[ycol].values,
        gridsize=nbin,
        extent=plot_range,
        cmap='jet',
        mincnt=1,
    )
    ax.set_xlim(0, xmax * axis_margin_rate)
    ax.set_ylim(0, ymax * axis_margin_rate)
    fig.tight_layout(pad=0.25, w_pad=1.5, h_pad=0.5)
    for ext in ['svg', 'pdf']:
        outpath = os.path.join(outbase + '.' + ext)
        fig.savefig(outpath, format=ext) 
    return None   

def combine_pdfs(pdf_file_paths, output_path):
    """
    Combines multiple PDF files into a single PDF file using pypdf's PdfWriter.
    """
    # Check that all PDF file paths exist
    for pdf_file_path in pdf_file_paths:
        if not os.path.exists(pdf_file_path):
            raise Exception(f"PDF file not found: {pdf_file_path}")
            return None

    writer = PdfWriter()

    # Read and append each PDF
    for pdf_path in pdf_file_paths:
        reader = PdfReader(pdf_path)
        for page in reader.pages:
            writer.add_page(page)

    # Write out the combined PDF
    with open(output_path, "wb") as f_out:
        writer.write(f_out)

    print(f"Successfully created combined PDF: {output_path}")
    return None

def extract_pdb_id(indir):
    pattern = re.compile(r'^csubst_site\.([^.]+)\.fa$')
    if not os.path.isdir(indir):
        print(f"PDB ID search directory not found: {indir}", flush=True)
        return None
    for filename in os.listdir(indir):
        match = pattern.search(filename)
        if match:
            pdb_id = match.group(1)
            print(f"Found file: {filename} --> Extracted PDB ID: {pdb_id}", flush=True)
            return pdb_id
    print(f"PDB ID not found in: {indir}", flush=True)
    return None

def get_annotation_text(og, arity, branch_id_str, trait, min_OCNany2spe, min_omegaCany2spe, min_OCNCoD, besthit_values):
    if besthit_values is None:
        besthit_values = ["NA", "NA", "NA", "NA", "NA"]
    if len(besthit_values) < 5:
        besthit_values = list(besthit_values) + (["NA"] * (5 - len(besthit_values)))
    annotation_text = f"""Orthogroup: {og}

K: {arity}

Branch ID: {branch_id_str}

Trait: {trait}

Minimum OCNany2spe: {min_OCNany2spe}

Minimum omegaCany2spe: {min_omegaCany2spe}

Minimum OCNCoD: {min_OCNCoD}

Time: {datetime.datetime.now()}

Annotation of the gene as being in the 5th, 25th, 50th, 75th, or 95th percentile for length within its orthogroup:

{besthit_values[0]}

{besthit_values[1]}

{besthit_values[2]}

{besthit_values[3]}

{besthit_values[4]}
"""
    return annotation_text

def process_index(og, branch_id_str, dir_out, dir_og, file_trait_color, ncpu, annotation_text):
    previous_cwd = os.getcwd()
    dir_iqtree_anc = os.path.join(dir_og, "iqtree.anc")
    dir_tree_plot = os.path.join(dir_og, "tree_plot")
    dir_stat_branch = os.path.join(dir_og, "stat.branch")
    dir_out_og = os.path.join(dir_out, og+'_'+branch_id_str.replace(',', '_'))
    print('{}, --branch_ids {}: wd: {}'.format(og, branch_id_str, dir_out_og), flush=True)
    if not os.path.exists(dir_out_og):
        os.makedirs(dir_out_og)
    os.chdir(dir_out_og)
    iqtree_tree_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.treefile')
    iqtree_state_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.state')
    iqtree_rate_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.rate')
    iqtree_iqtree_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.iqtree')
    iqtree_log_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.log')
    iqtree_ckp_file = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.ckp.gz')
    file_summary = os.path.join(dir_out_og, f'summary.{og}_branch_id{branch_id_str}.pdf')
    if os.path.exists(file_summary):
        print(f'Skipped. Outfile already exists: {file_summary}', flush=True)
        return og, None
    try:
        # Find csubst_site output tsv: csubst_site.tsv or csubst_site.<pdb_id>.tsv
        dir_csubst_site = os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str)
        file_csubst_out_candidates = glob.glob(os.path.join(dir_csubst_site, 'csubst_site*.tsv'))
        # Exclude known fixed filenames
        exclude_files = {'csubst_site.state_N.tsv', 'csubst_site.state_S.tsv'}
        file_csubst_out_list = [f for f in file_csubst_out_candidates if os.path.basename(f) not in exclude_files]
        file_csubst_out = file_csubst_out_list[0] if file_csubst_out_list else None
        if file_csubst_out is not None and os.path.exists(file_csubst_out):
            print(f'Skipped csubst site. Outfile already exists: {file_csubst_out}', flush=True)
        else:
            print(f'Running csubst site. Output file not found: {file_csubst_out}', flush=True)
            path_iqtree_zip = os.path.join(dir_iqtree_anc, og+'.iqtree.anc.zip')
            with zipfile.ZipFile(path_iqtree_zip, "r") as zip_ref:
                zip_ref.extractall(dir_out_og)
            cmd = ['csubst', 'site']
            cmd += ['--alignment_file', os.path.join(og+'.iqtree.anc', 'csubst.fasta')]
            cmd += ['--rooted_tree_file', os.path.join(og+'.iqtree.anc', 'csubst.nwk')]
            cmd += ['--branch_id', branch_id_str]
            cmd += ['--ml_anc', 'no']
            cmd += ['--threads', str(max(1, int(ncpu)))]
            cmd += ['--iqtree_treefile', iqtree_tree_file]
            cmd += ['--iqtree_state', iqtree_state_file]
            cmd += ['--iqtree_rate', iqtree_rate_file]
            cmd += ['--iqtree_iqtree', iqtree_iqtree_file]
            cmd += ['--iqtree_log', iqtree_log_file]
            cmd += ['--mafft_exe', 'mafft']
            cmd += ['--pdb', 'besthit']
            print('COMMAND: {}'.format(' '.join(cmd)), flush=True)
            subprocess.run(cmd, check=True)
        print(f'{datetime.datetime.now()}: csubst site done: {og}', flush=True)
        file_tree_plot = os.path.join(dir_out_og, og+'.tree_plot.pdf')
        if os.path.exists(file_tree_plot):
            print(f'Tree plot skipped: outfile already exists: {file_tree_plot}', flush=True)
        else:
            print(f'Running stat_branch2tree_plot to generate a tree pdf file: {file_tree_plot}', flush=True)
            run_stat_branch2tree_plot(
                og,
                branch_id_str,
                file_trait_color=file_trait_color,
                dir_out_og=dir_out_og,
                dir_og=dir_og,
                ncpu=ncpu,
            )
            print(f'{datetime.datetime.now()}: stat_branch2tree_plot done: {og}', flush=True)
        if os.path.exists('annotation_text.pdf'):
            print(f'Annotation text generation skipped: outfile already exists: annotation_text.pdf for {og}', flush=True)
        else:
            print(f'Generating annotation text: annotation_text.pdf for {og}', flush=True)
            create_pdf(annotation_text, 'annotation_text.pdf')
        pdb_id = extract_pdb_id(indir=os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str))
        path_pymol_pdf = os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.{pdb_id}.pymol.pdf') if pdb_id is not None else None
        if pdb_id is None:
            pdf_file_paths = [
                os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.pdf'),
                file_tree_plot,
                'annotation_text.pdf',
            ]
        elif not os.path.exists(path_pymol_pdf):
            pdf_file_paths = [
                os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.{pdb_id}.pdf'),
                file_tree_plot,
                'annotation_text.pdf',
            ]
        else:
            pdf_file_paths = [
                os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.{pdb_id}.pdf'),
                os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.{pdb_id}.pymol.pdf'),
                file_tree_plot,
                'annotation_text.pdf',
            ]
        combine_pdfs(pdf_file_paths, output_path=file_summary)
        remove_files = [
            'tmp.csubst.sub_tensor.S.mmap',
            'tmp.csubst.sub_tensor.N.mmap',
            'Rplots.pdf',
            'species_color.tsv',
            'annotation_text.pdf',
            file_tree_plot,
            iqtree_tree_file,
            iqtree_state_file,
            iqtree_rate_file,
            iqtree_iqtree_file,
            iqtree_log_file,
            iqtree_ckp_file,
        ] + glob.glob('*.cif')
        if os.path.exists(file_summary):
            for remove_file in remove_files:
                if os.path.exists(remove_file): os.remove(remove_file)
        return og,None
    except Exception as e:
        print(f'process_index: Error in {og}: {e}', flush=True)
        return og,e
    finally:
        os.chdir(previous_cwd)

def skip_lower_order(cb_passed, trait, already_analyzed_in_greater_K):
    num_before_filtering = cb_passed.shape[0]
    bid_cols = cb_passed.columns[cb_passed.columns.str.startswith('branch_id_')]
    is_already_analyzed = numpy.zeros(cb_passed.shape[0], dtype=bool)
    for i in cb_passed.index:
        flag_analyzed = False
        og = cb_passed.at[i,'orthogroup']
        branch_id_set = set(cb_passed.loc[i, bid_cols].tolist())
        if not og in already_analyzed_in_greater_K[trait].keys():
            flag_analyzed = False
            already_analyzed_in_greater_K[trait][og] = list()
        else:
            for already_analyzed_set in already_analyzed_in_greater_K[trait][og]:
                if branch_id_set.issubset(already_analyzed_set):
                    flag_analyzed = True
                    break
        if flag_analyzed:
            is_already_analyzed[i] = True
            print(f'Skipped. Subset of already analyzed higher order convergence: {og} {branch_id_set}', flush=True)
        else:
            already_analyzed_in_greater_K[trait][og].append(branch_id_set)
    cb_passed = cb_passed.loc[~is_already_analyzed,:].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity}: Skipped branch combinations that are subsets of already analyzed higher-order convergence: {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed, already_analyzed_in_greater_K

def generate_trait_colors(file_trait, trait_names):
    df_trait = pandas.read_csv(file_trait, sep='\t', header=0, index_col=None)
    for trait_name in trait_names:
        df_trait_color = df_trait.loc[:,['species', trait_name]]
        is_foreground = (df_trait_color[trait_name]==1)
        df_trait_color.columns = ['species', 'color']
        df_trait_color['color'] = 'black'
        df_trait_color.loc[is_foreground,'color'] = 'firebrick'
        df_trait_color.to_csv(f'trait_{trait_name}.color.tsv', sep='\t', index=False)
    return None

def resolve_existing_path(candidates):
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return None

def get_alignment_for_tree_plot(dir_og, og, dir_out_og):
    alignment_candidates = [
        os.path.join(dir_og, 'clipkit', og + '.cds.clipkit.fa.gz'),
        os.path.join(dir_og, 'clipkit', og + '.cds.clipkit.fasta'),
        os.path.join(dir_og, 'clipkit', og + '.cds.clipkit.fa'),
    ]
    alignment_path = resolve_existing_path(alignment_candidates)
    if alignment_path is None:
        raise FileNotFoundError(
            f'Alignment file was not found for tree plotting. Checked: {alignment_candidates}'
        )
    if alignment_path.endswith('.gz'):
        plain_path = os.path.join(dir_out_og, og + '.cds.clipkit.plot.fasta')
        if (not os.path.exists(plain_path)) or (os.path.getmtime(plain_path) < os.path.getmtime(alignment_path)):
            with gzip.open(alignment_path, 'rt') as fin, open(plain_path, 'w') as fout:
                shutil.copyfileobj(fin, fout)
        return plain_path
    return alignment_path

def run_stat_branch2tree_plot(og, branch_id_str, file_trait_color, dir_out_og, dir_og, ncpu=1):
    dir_myscript = os.path.realpath(os.path.dirname(__file__))
    dir_tree_annotation = os.path.join(dir_myscript, 'tree_annotation')
    file_stat_branch = os.path.join(dir_og, 'stat.branch', og+'.stat.branch.tsv')
    file_og_rpsblast = os.path.join(dir_og, 'rpsblast', og+'.rpsblast.tsv')
    file_og_alignment = get_alignment_for_tree_plot(dir_og=dir_og, og=og, dir_out_og=dir_out_og)
    file_csubst_input_fasta = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.fasta')
    pdb_id = extract_pdb_id(indir=os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str))
    if pdb_id is None:
        file_csubst_site_tsv = os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.tsv')
    else:
        file_csubst_site_tsv = os.path.join(dir_out_og, 'csubst_site.branch_id'+branch_id_str, f'csubst_site.{pdb_id}.tsv')
    df_csubst_site = pandas.read_csv(file_csubst_site_tsv, sep='\t', header=0, index_col=None)
    convergent_sites = df_csubst_site.loc[(df_csubst_site['OCNany2spe'] > 0.5), 'codon_site_alignment'].tolist()
    convergent_site_str = ':'.join([ str(cs) for cs in convergent_sites ])
    print(f'Convergent sites extracted from {file_csubst_site_tsv}: {convergent_site_str}', flush=True)
    file_tree_plot_out = og+'.tree_plot.pdf'
    if os.path.exists(file_tree_plot_out):
        print(f'Tree plot skipped: outfile already exists: {file_tree_plot_out}', flush=True)
        return None
    cmd = ['Rscript', os.path.join(dir_myscript, 'stat_branch2tree_plot.r')]
    cmd.append('--stat_branch='+file_stat_branch)
    cmd.append('--tree_annotation_dir='+dir_tree_annotation)
    cmd.append('--max_delta_intron_present=-0.5')
    cmd.append('--width=7.2')
    cmd.append('--rel_widths=')
    cmd.append('--panel1=tree,bl_rooted,support_unrooted,species,L')
    cmd.append('--panel2=heatmap,no,abs,_,expression_')
    #cmd.append('--panel3=pointplot,no,rel,_,expression_')
    cmd.append('--panel4=cluster_membership,100000')
    cmd.append('--panel5=tiplabel')
    cmd.append('--panel6=signal_peptide')
    cmd.append('--panel7=transmembrane_domain')
    cmd.append('--panel8=intron_number')
    cmd.append('--panel9=domain,'+file_og_rpsblast)
    cmd.append(f'--panel10=amino_acid_site,1,{convergent_site_str},{file_csubst_input_fasta}')
    cmd.append('--panel11=alignment,'+file_og_alignment)
    cmd.append('--panel12=fimo,2000,0.05')
    cmd.append('--show_branch_id=yes')
    cmd.append('--event_method=species_overlap')
    cmd.append('--species_color_table='+file_trait_color)
    cmd.append('--branch_combination='+branch_id_str)
    cmd.append('--ncpu='+str(max(1, int(ncpu))))
    print('COMMAND: {}'.format(' '.join(cmd)), flush=True)
    subprocess.run(cmd, check=True)
    shutil.move('stat_branch2tree_plot.pdf', file_tree_plot_out)
    return None

def find_file_trait_color(trait):
    file_name = f'trait_{trait}.color.tsv'
    file_trait_color = os.path.join(os.getcwd(), file_name)
    return file_trait_color

def filter_max_per_og(cb_passed, arity, trait):
    sampled_indices = []
    num_before_filtering = cb_passed.shape[0]
    for og in cb_passed['orthogroup'].drop_duplicates():
        og_indices = cb_passed[cb_passed['orthogroup'] == og].index
        if (og_indices.shape[0] > args.max_per_og)&(args.max_per_og > 0):
            print(f'K = {arity}, {og}: Number of branch combinations ({og_indices.shape[0]}) exceeded --max_per_og. Only {args.max_per_og} combinations will be analyzed.')
            sampled_indice_positions = numpy.floor(numpy.linspace(0, og_indices.shape[0]-1, args.max_per_og)).astype(int)
            sampled_indices += og_indices[sampled_indice_positions].tolist()
        else:
            sampled_indices += og_indices.tolist()
    cb_passed = cb_passed.loc[sampled_indices,:].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity}: Skipped branch combinations that are not top {args.max_per_og} within individual orthogroups: {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed

def filter_fg_stem_ratio(cb_passed, arity, args, no_trait_name):
    min_fg_stem_num = int(numpy.round(args.min_fg_stem_ratio * arity, decimals=0))
    num_before_filtering = cb_passed.shape[0]
    if no_trait_name:
        col = 'branch_num_fg_stem'
    else:
        col = 'branch_num_fg_stem_'+trait
    cb_passed = cb_passed.loc[cb_passed[col] >= min_fg_stem_num,:].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity}: Skipped branch combinations with less than {min_fg_stem_num} foreground stem branches (prop = {args.min_fg_stem_ratio}): {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed

def write_annotated_table(cb_passed, dir_of, dir_out, out_name):
    path_annot = os.path.join(dir_of, 'Orthogroups', 'Orthogroups.GeneCount.annotated.tsv')
    annot = pandas.read_csv(path_annot, sep='\t', header=0, low_memory=False)
    annot.columns = annot.columns.str.replace('Orthogroup', 'orthogroup')
    cols = ['orthogroup'] + annot.columns[annot.columns.str.startswith('besthit_')].tolist()
    annot = annot.loc[:,cols]
    cb_passed = pandas.merge(cb_passed, annot, how='left', on='orthogroup', sort=False)
    del annot
    outfile = os.path.join(dir_out, out_name+'.tsv')
    print(f'Writing the table to: {outfile}', flush=True)
    cb_passed.to_csv(outfile, sep='\t', index=False)
    return cb_passed

def filter_max_per_K_1st(cb_passed, arity, trait):
    print(f'Number of branch combinations after filtering exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}.')
    print(f'Selecting one branch combination per orthogroup with the highest OCN.', flush=True)
    num_before_filtering = cb_passed.shape[0]
    cb_passed = cb_passed.sort_values(by=['orthogroup','OCNany2spe'], ascending=[True, False])
    cb_passed = cb_passed.drop_duplicates(subset=['orthogroup'], keep='first').reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity} for {trait}: Skipped branch combinations with non-highest OCN within orthogroups: {num_before_filtering} -> {num_after_filtering}', flush=True)
    outfile = os.path.join(dir_out, out_name + '.txt')
    print(f'Writing a placeholder file to: {outfile}', flush=True)
    with open(outfile, 'w') as f:
        f.write(f'One branch combination per orthogroup was analyzed because the number of branch combinations after all filtering still exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}.\n')
    return cb_passed    

def filter_max_per_K_2nd(cb_passed, arity, trait):
    print(f'Number of branch combinations after one-per-orthogroup filtering still exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}.')
    print(f'Only the top {args.max_per_K} high-OCN branch combinations will be analyzed for K = {arity} of {trait}', flush=True)
    num_before_filtering = cb_passed.shape[0]
    cb_passed = cb_passed.sort_values(by='OCNany2spe', ascending=False).iloc[:args.max_per_K,:].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity} for {trait}: Skipped branch combinations that are not top {args.max_per_K} at K = {arity}: {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed

def create_pdf(text, filename):
    # Constants
    INCH = 72  # 1 inch = 72 points
    PAGE_WIDTH_INCH = 7.2  # Standard letter width
    PAGE_HEIGHT_INCH = 9.7  # Standard letter height
    PAGE_WIDTH = PAGE_WIDTH_INCH * INCH  # 612 points
    PAGE_HEIGHT = PAGE_HEIGHT_INCH * INCH  # 792 points
    TOP_MARGIN = INCH  # 1 inch
    BOTTOM_MARGIN = INCH  # 1 inch
    LEFT_MARGIN = INCH  # 1 inch
    RIGHT_MARGIN = INCH  # 1 inch
    LINE_HEIGHT = 10  # points
    FONT_SIZE = 8
    FONT_NAME = "Helvetica"
    
    # Register a TrueType font if needed (optional)
    # pdfmetrics.registerFont(TTFont('Helvetica', 'Helvetica.ttf'))
    
    # Create a canvas with standard letter size
    c = canvas.Canvas(filename, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))
    
    # Set font
    c.setFont(FONT_NAME, FONT_SIZE)
    
    # Calculate the maximum width available for text
    max_width = PAGE_WIDTH - LEFT_MARGIN - RIGHT_MARGIN
    
    # Approximate average character width (depends on font and size)
    # Helvetica has an average width of about 0.5 * font size
    # This is a simplification; for precise wrapping, use string width calculations
    avg_char_width = FONT_SIZE * 0.5
    max_chars_per_line = int(max_width / avg_char_width)
    
    # Split the text into paragraphs
    paragraphs = text.split('\n\n')
    
    # Starting position
    x = LEFT_MARGIN
    y = PAGE_HEIGHT - TOP_MARGIN
    
    for para in paragraphs:
        # Use textwrap to wrap each paragraph
        wrapped_lines = textwrap.wrap(para, width=max_chars_per_line)
        for line in wrapped_lines:
            if y < BOTTOM_MARGIN:
                c.showPage()
                c.setFont(FONT_NAME, FONT_SIZE)
                y = PAGE_HEIGHT - TOP_MARGIN
            c.drawString(x, y, line)
            y -= LINE_HEIGHT
        # Add extra space after each paragraph
        y -= LINE_HEIGHT / 2
    
    # Save the PDF
    c.save()

if __name__ == '__main__':
    print(f'{datetime.datetime.now()}: {__file__} started.', flush=True)

    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_orthogroup', metavar='PATH', default='../workspace/output/orthogroup', type=str, help='')
    parser.add_argument('--dir_orthofinder', metavar='PATH', default='../workspace/output/orthofinder', type=str, help='')
    parser.add_argument('--dir_out', metavar='PATH', default='.', type=str, help='')
    parser.add_argument('--file_trait', metavar='PATH', default='.', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=2, type=int, help='Number of CPUs to use.')
    parser.add_argument('--arity_range', metavar='STR', default='3-10', type=str,
                        help='Hyphen-separated range of arity (K).')
    parser.add_argument('--trait', metavar='STR', default='all', type=str,
                        help='Comma-separated list of traits for analysis. "all" for all traits.')    
    parser.add_argument('--skip_lower_order', metavar='yes|no', default='yes', type=strtobool,
                        help='Skip if the branch combination is a subset of already analyzed higher-order combinations.')
    parser.add_argument('--min_fg_stem_ratio', metavar='FLOAT', default=0.5, type=float, help='Minimum proportion of foreground stem branches in branch combinations. For exampmle, 3 or more branches should be foreground stems at K = 6 if set to 0.5.')
    parser.add_argument('--max_per_og', metavar='INT', default=0, type=int, help='The maximum number of branch combinations to analyze per orthogroup per K.')
    parser.add_argument('--min_OCNany2spe', metavar='FLOAT', default=5.0, type=float, help='')
    parser.add_argument('--min_omegaCany2spe', metavar='FLOAT', default=5.0, type=float, help='')
    parser.add_argument('--min_OCNCoD', metavar='FLOAT', default=5.0, type=float, help='')
    parser.add_argument('--max_per_K', metavar='INT', default=300, type=int,
                        help='Maximum number of branch combinations to analyze per K per trait. '
                             'If exceeded, only one branch combinations per orthogroup will be analyzed. '
                             'If still exceeded, the K will be skipped.')
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))

    if not 'cwd' in locals():
        cwd = os.getcwd()
    print('Working at: {}'.format(cwd))
    min_OCNany2spe = args.min_OCNany2spe
    min_omegaCany2spe = args.min_omegaCany2spe
    min_OCNCoD = args.min_OCNCoD
    dir_og = os.path.realpath(args.dir_orthogroup)
    dir_of = os.path.realpath(args.dir_orthofinder)
    db_file = os.path.join(dir_og, 'gg_orthogroup.db')
    arity_min = int(args.arity_range.split('-')[0])
    arity_max = int(args.arity_range.split('-')[1])
    arity_range = numpy.arange(arity_min, arity_max+1)[::-1]
    print('Arity range: {}'.format(' '.join([ str(a) for a in arity_range ])), flush=True)
    conn = sqlalchemy.create_engine("sqlite:///" + db_file)
    inspector = sqlalchemy.inspect(conn)
    table_names = inspector.get_table_names()
    print(f"Tables in the database: {' '.join(table_names)}", flush=True)
    column_names = dict()
    for table in table_names:
        print('')
        query = sqlalchemy.text("PRAGMA TABLE_INFO(" + table + ")")
        columns = pandas.read_sql_query(sql=query, con=conn.connect(), index_col=None, coerce_float=True, chunksize=None)
        column_names[table] = columns['name'].tolist()
        print(f"Column names of {table}: {' '.join(column_names[table])}")
    trait_names = args.trait.split(',')
    if trait_names[0] == 'all':
        for arity in arity_range:
            if 'cb' + str(arity) in column_names.keys():
                colnames = column_names['cb' + str(arity)]
                trait_names = [ colname.replace('is_fg_', '') for colname in colnames if colname.startswith('is_fg_') ]
                break
    if len(trait_names) == 0:
        no_trait_name = True
        df_trait = pandas.read_csv(args.file_trait, sep='\t', header=0, index_col=None)
        print(f'No trait name found in cb tables. Analyzing only the first trait in {args.file_trait}: {df_trait.columns[1]}', flush=True)
        trait_names = [df_trait.columns[1],]
    else:
        no_trait_name = False
        print(f'Specified or detected traits: {",".join(trait_names)}', flush=True)
    generate_trait_colors(file_trait=args.file_trait, trait_names=trait_names)
    already_analyzed_in_greater_K = dict()
    for trait in trait_names:
        already_analyzed_in_greater_K[trait] = dict() # Initialize
    for arity in arity_range:
        print('')
        if 'cb'+str(arity) not in table_names:
            print(f'cb{arity} is not in the database. Skipping.', flush=True)
            continue
        sql_cb = sqlalchemy.text(f'SELECT * from cb{arity}')
        cb = pandas.read_sql_query(sql=sql_cb, con=conn.connect(), index_col=None, coerce_float=True)
        sum_ocn = cb.loc[:, 'OCNany2spe'].sum()
        sum_ecn = cb.loc[:, 'ECNany2spe'].sum()
        sum_ocs = cb.loc[:, 'OCSany2spe'].sum()
        sum_ecs = cb.loc[:, 'ECSany2spe'].sum()
        print(f'Total observed/expected nonsynonymous convergence at K = {arity}: {int(sum_ocn):,} / {int(sum_ecn):,}')
        print(f'Total observed/expected synonymous convergence at K = {arity}: {int(sum_ocs):,} / {int(sum_ecs):,}')
        for trait in trait_names:
            flag_zip = True
            out_name = f'csubst_site_{trait}_K{arity}_minOCN{min_OCNany2spe}_minomegaC{min_omegaCany2spe}_minOCNCoD{min_OCNCoD}'
            print(f'{datetime.datetime.now()}: Started processing: {out_name}', flush=True)
            file_trait_color = find_file_trait_color(trait)
            dir_out = os.path.realpath(os.path.join(args.dir_out, out_name))
            out_zip = os.path.join(cwd, dir_out+'.zip')
            if (not os.path.exists(dir_out)) & (not os.path.exists(out_zip)):
                os.makedirs(dir_out)
            conditions = True
            conditions &= (cb['OCNany2spe']>=min_OCNany2spe)
            conditions &= (cb['omegaCany2spe']>=min_omegaCany2spe)
            conditions &= (cb['OCNCoD']>=min_OCNCoD)
            if 'is_fg' in cb.columns:
                conditions &= (cb.loc[:,'is_fg'] == 'Y')
            else:
                conditions &= (cb.loc[:,'is_fg_'+trait] == 'Y')
            print(f'Number of branch combinations for {trait} with OCNany2spe>={min_OCNany2spe}, '
                  f'omegaCany2spe>={min_omegaCany2spe}, and OCNCoD>={min_OCNCoD} at K={arity}: {conditions.sum()} '
                  f'(before removing branch combinations already analyzed in greater K)', flush=True)
            if conditions.sum()==0:
                outfile = os.path.join(dir_out, out_name+'.txt')
                print('No branch combinations to analyze.', flush=True)
                if os.path.exists(dir_out) & (not os.path.exists(out_zip)):
                    print(f'Writing a placeholder file to: {outfile}', flush=True)
                    with open(outfile, 'w') as f:
                        f.write('No analyzable branch combinations.\n')
                    shutil.make_archive(dir_out, 'zip', dir_out)
                    shutil.rmtree(dir_out)
                continue
            cb_passed = cb.loc[conditions,:].sort_values(by='orthogroup').reset_index(drop=True)                
            if os.path.exists(out_zip):
                print(f'Skipped arity = {arity}. Output zip file already exists: {out_zip}', flush=True)
                continue
            plot_scatter(cb=cb_passed, xcol='OCNany2spe', ycol='omegaCany2spe',
                         outbase=os.path.join(dir_out, out_name + '_OCNany2spe-omegaC'),
                         polygon_xmin=min_OCNany2spe, polygon_ymin=min_omegaCany2spe)
            plot_scatter(cb=cb_passed, xcol='OCNany2spe', ycol='OCNCoD',
                         outbase=os.path.join(dir_out, out_name + '_OCNany2spe-OCNCoD'),
                         polygon_xmin=min_OCNany2spe, polygon_ymin=min_OCNCoD)
            if args.skip_lower_order:
                out = skip_lower_order(cb_passed, trait, already_analyzed_in_greater_K)
                cb_passed, already_analyzed_in_greater_K = out
                del out
            if args.min_fg_stem_ratio > 0:
                cb_passed = filter_fg_stem_ratio(cb_passed, arity, args, no_trait_name)
            cb_passed = write_annotated_table(cb_passed, dir_of, dir_out, out_name)
            cb_passed = filter_max_per_og(cb_passed, arity, trait)            
            if args.max_per_K < cb_passed.shape[0]:
                cb_passed = filter_max_per_K_1st(cb_passed, arity, trait)
            if args.max_per_K < cb_passed.shape[0]:
                cb_passed = filter_max_per_K_2nd(cb_passed, arity, trait)
            print('Starting analyzing individual branch combinations.', flush=True)
            branch_id_cols = cb_passed.columns[cb_passed.columns.str.startswith('branch_id_')]
            og_list = cb_passed['orthogroup'].tolist()
            branch_id_str_list = [ ','.join([ str(bid) for bid in cb_passed.loc[i,branch_id_cols].tolist() ]) for i in range(cb_passed.shape[0]) ]
            besthit_cols = ['besthit_0.05', 'besthit_0.25', 'besthit_0.5', 'besthit_0.75', 'besthit_0.95']
            annotation_text_list = []
            for i in range(cb_passed.shape[0]):
                besthit_values = []
                for col in besthit_cols:
                    if (col in cb_passed.columns) and pandas.notna(cb_passed.at[i, col]):
                        besthit_values.append(str(cb_passed.at[i, col]))
                    else:
                        besthit_values.append('NA')
                annotation_text_list.append(
                    get_annotation_text(
                        og=cb_passed.at[i, 'orthogroup'],
                        arity=arity,
                        branch_id_str=branch_id_str_list[i],
                        trait=trait,
                        min_OCNany2spe=min_OCNany2spe,
                        min_omegaCany2spe=min_omegaCany2spe,
                        min_OCNCoD=min_OCNCoD,
                        besthit_values=besthit_values,
                    )
                )
            dir_out_list = [ dir_out for i in range(cb_passed.shape[0]) ]
            dir_og_list = [dir_og,] * cb_passed.shape[0]
            file_trait_color_list = [file_trait_color] * cb_passed.shape[0]
            ncpu_list = [args.ncpu] * cb_passed.shape[0]
            with ProcessPoolExecutor(max_workers=args.ncpu) as executor:
                results = executor.map(
                    process_index,
                    og_list,
                    branch_id_str_list,
                    dir_out_list,
                    dir_og_list,
                    file_trait_color_list,
                    ncpu_list,
                    annotation_text_list,
                )
                for og,result in results:
                    if isinstance(result, Exception):
                        print(f"Error in {og}: {result}")
                        flag_zip = False
            os.chdir(cwd)
            if flag_zip:
                print(f'No error detected. Zipping {dir_out}', flush=True)
                shutil.make_archive(dir_out, 'zip', dir_out)
                shutil.rmtree(dir_out)
            else:
                print(f'Error detected. Skipping zipping {dir_out}', flush=True)
    conn.dispose()
    print(f'{datetime.datetime.now()}: {__file__} completed!', flush=True)
