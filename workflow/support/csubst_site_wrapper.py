#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import glob
import gzip
import numpy
import os
import pandas
import re
import shutil
import subprocess
import textwrap
import zipfile
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
try:
    import sqlalchemy
except ImportError:
    sqlalchemy = None

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
try:
    from pypdf import PdfReader, PdfWriter
except ImportError:
    PdfReader = None
    PdfWriter = None
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
except ImportError:
    canvas = None
    letter = None
    pdfmetrics = None
    TTFont = None

pandas.set_option("display.max_columns", None)


def _get_matplotlib():
    import matplotlib

    matplotlib.rcParams['font.size'] = 8
    matplotlib.rcParams['font.family'] = 'Helvetica'
    matplotlib.rcParams['svg.fonttype'] = 'none' # none, path, or svgfont
    return matplotlib


def require_sqlalchemy():
    if sqlalchemy is None:
        raise ImportError("sqlalchemy is required to read orthogroup databases.")


def require_pdf_dependencies():
    if PdfReader is None or PdfWriter is None:
        raise ImportError("pypdf is required to generate combined PDF reports.")
    if canvas is None or letter is None or pdfmetrics is None or TTFont is None:
        raise ImportError("reportlab is required to generate annotation PDFs.")

def plot_scatter(cb, xcol, ycol, outbase, xmax=10, ymax=10, nbin=50, polygon_xmin=0, polygon_ymin=0):
    matplotlib = _get_matplotlib()
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
    require_pdf_dependencies()
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
    pattern = re.compile(r'^csubst_sites?\.([^.]+)\.fa$')
    if not os.path.isdir(indir):
        print(f"PDB ID search directory not found: {indir}", flush=True)
        return None
    for filename in sorted(f for f in os.listdir(indir) if not f.startswith('.')):
        match = pattern.search(filename)
        if match:
            pdb_id = match.group(1)
            print(f"Found file: {filename} --> Extracted PDB ID: {pdb_id}", flush=True)
            return pdb_id
    print(f"PDB ID not found in: {indir}", flush=True)
    return None


def get_site_dir_candidates(dir_out_og, branch_id_str):
    return [
        os.path.join(dir_out_og, 'csubst_sites.branch_id' + branch_id_str),
        os.path.join(dir_out_og, 'csubst_site.branch_id' + branch_id_str),
    ]


def resolve_site_output_dir(dir_out_og, branch_id_str):
    candidates = get_site_dir_candidates(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
    for candidate in candidates:
        if os.path.isdir(candidate):
            return candidate
    return candidates[0]


def load_site_output_manifest(site_dir):
    manifest_candidates = [
        os.path.join(site_dir, 'csubst_sites.outputs.tsv'),
        os.path.join(site_dir, 'csubst_site.outputs.tsv'),
    ]
    manifest_path = resolve_existing_path(manifest_candidates)
    if manifest_path is None:
        return None
    return pandas.read_csv(manifest_path, sep='\t', header=0, index_col=None)


def resolve_manifest_output_path(site_dir, row):
    output_path = str(row.get('output_path', '')).strip()
    if output_path not in ('', 'nan'):
        if os.path.isabs(output_path):
            return output_path
        return os.path.join(site_dir, output_path)
    output_file = str(row.get('output_file', '')).strip()
    if output_file not in ('', 'nan'):
        return os.path.join(site_dir, output_file)
    return None


def find_site_output_from_manifest(manifest_df, site_dir, output_kinds):
    if manifest_df is None:
        return None
    if 'output_kind' not in manifest_df.columns:
        return None
    if isinstance(output_kinds, str):
        output_kinds = [output_kinds]
    manifest_df = manifest_df.copy()
    if 'file_exists' in manifest_df.columns:
        is_existing = manifest_df['file_exists'].astype(str).str.upper().isin({'Y', 'YES', 'TRUE', '1'})
        manifest_df = manifest_df.loc[is_existing, :].reset_index(drop=True)
    for output_kind in output_kinds:
        rows = manifest_df.loc[manifest_df['output_kind'].astype(str) == str(output_kind), :]
        for _, row in rows.iterrows():
            output_path = resolve_manifest_output_path(site_dir=site_dir, row=row)
            if (output_path is not None) and os.path.exists(output_path):
                return output_path
    return None


def find_first_matching_file(site_dir, glob_patterns, exclude_predicate=None):
    for glob_pattern in glob_patterns:
        for candidate in sorted(glob.glob(os.path.join(site_dir, glob_pattern))):
            if exclude_predicate is not None and exclude_predicate(candidate):
                continue
            if os.path.isfile(candidate):
                return candidate
    return None


def is_auxiliary_site_tsv(path):
    basename = os.path.basename(path)
    return (
        basename.endswith('.outputs.tsv')
        or basename.endswith('.state_N.tsv')
        or basename.endswith('.state_S.tsv')
        or '.tree_site.' in basename
    )


def is_auxiliary_site_pdf(path):
    basename = os.path.basename(path)
    return (
        basename.endswith('.pymol.pdf')
        or basename.endswith('.state.pdf')
        or '.tree_site.' in basename
    )


def resolve_site_artifacts(dir_out_og, branch_id_str):
    site_dir = resolve_site_output_dir(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
    manifest_df = load_site_output_manifest(site_dir=site_dir)
    pdb_id = extract_pdb_id(indir=site_dir)

    site_table_candidates = []
    site_summary_pdf_candidates = []
    pymol_pdf_candidates = []
    if pdb_id is not None:
        site_table_candidates.extend([
            os.path.join(site_dir, f'csubst_sites.{pdb_id}.tsv'),
            os.path.join(site_dir, f'csubst_site.{pdb_id}.tsv'),
        ])
        site_summary_pdf_candidates.extend([
            os.path.join(site_dir, f'csubst_sites.{pdb_id}.pdf'),
            os.path.join(site_dir, f'csubst_site.{pdb_id}.pdf'),
        ])
        pymol_pdf_candidates.extend([
            os.path.join(site_dir, f'csubst_sites.{pdb_id}.pymol.pdf'),
            os.path.join(site_dir, f'csubst_site.{pdb_id}.pymol.pdf'),
        ])
    site_table_candidates.extend([
        os.path.join(site_dir, 'csubst_sites.tsv'),
        os.path.join(site_dir, 'csubst_site.tsv'),
    ])
    site_summary_pdf_candidates.extend([
        os.path.join(site_dir, 'csubst_sites.pdf'),
        os.path.join(site_dir, 'csubst_site.pdf'),
    ])

    site_table_tsv = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=['site_table_tsv'],
    )
    if site_table_tsv is None:
        site_table_tsv = resolve_existing_path(site_table_candidates)
    if site_table_tsv is None:
        site_table_tsv = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=['csubst_sites*.tsv', 'csubst_site*.tsv'],
            exclude_predicate=is_auxiliary_site_tsv,
        )

    site_summary_pdf = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=['site_summary_pdf'],
    )
    if site_summary_pdf is None:
        site_summary_pdf = resolve_existing_path(site_summary_pdf_candidates)
    if site_summary_pdf is None:
        site_summary_pdf = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=['csubst_sites*.pdf', 'csubst_site*.pdf'],
            exclude_predicate=is_auxiliary_site_pdf,
        )

    pymol_summary_pdf = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=['pymol_summary_pdf'],
    )
    if pymol_summary_pdf is None:
        pymol_summary_pdf = resolve_existing_path(pymol_pdf_candidates)
    if pymol_summary_pdf is None:
        pymol_summary_pdf = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=['csubst_sites*.pymol.pdf', 'csubst_site*.pymol.pdf'],
        )

    return {
        'site_dir': site_dir,
        'manifest_df': manifest_df,
        'pdb_id': pdb_id,
        'site_table_tsv': site_table_tsv,
        'site_summary_pdf': site_summary_pdf,
        'pymol_summary_pdf': pymol_summary_pdf,
    }

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
    dir_iqtree_anc = os.path.join(dir_og, "iqtree_anc")
    dir_tree_plot = os.path.join(dir_og, "tree_plot")
    dir_stat_branch = os.path.join(dir_og, "stat_branch")
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
        artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
        file_csubst_out = artifacts['site_table_tsv']
        if file_csubst_out is not None and os.path.exists(file_csubst_out):
            print(f'Skipped csubst sites. Outfile already exists: {file_csubst_out}', flush=True)
        else:
            print(f'Running csubst sites. Output file not found: {file_csubst_out}', flush=True)
            path_iqtree_zip = os.path.join(dir_iqtree_anc, og+'.iqtree.anc.zip')
            with zipfile.ZipFile(path_iqtree_zip, "r") as zip_ref:
                zip_ref.extractall(dir_out_og)
            cmd = ['csubst', 'sites']
            cmd += ['--alignment_file', os.path.join(og+'.iqtree.anc', 'csubst.fasta')]
            cmd += ['--rooted_tree_file', os.path.join(og+'.iqtree.anc', 'csubst.nwk')]
            cmd += ['--branch_id', branch_id_str]
            cmd += ['--threads', str(max(1, int(ncpu)))]
            cmd += ['--iqtree_treefile', iqtree_tree_file]
            cmd += ['--iqtree_state', iqtree_state_file]
            cmd += ['--iqtree_rate', iqtree_rate_file]
            cmd += ['--iqtree_iqtree', iqtree_iqtree_file]
            cmd += ['--iqtree_log', iqtree_log_file]
            cmd += ['--pdb', 'besthit']
            print('COMMAND: {}'.format(' '.join(cmd)), flush=True)
            subprocess.run(cmd, check=True)
        print(f'{datetime.datetime.now()}: csubst sites done: {og}', flush=True)
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
        artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
        if artifacts['site_summary_pdf'] is None:
            raise FileNotFoundError(
                f'CSUBST site-summary PDF was not found in {artifacts["site_dir"]}.'
            )
        pdf_file_paths = [artifacts['site_summary_pdf']]
        if artifacts['pymol_summary_pdf'] is not None and os.path.exists(artifacts['pymol_summary_pdf']):
            pdf_file_paths.append(artifacts['pymol_summary_pdf'])
        pdf_file_paths.extend([file_tree_plot, 'annotation_text.pdf'])
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

def skip_lower_order(cb_passed, arity, trait, already_analyzed_in_greater_K):
    num_before_filtering = cb_passed.shape[0]
    bid_cols = cb_passed.columns[cb_passed.columns.str.startswith('branch_id_')].tolist()
    keep_indices = []
    analyzed_by_og = already_analyzed_in_greater_K[trait]
    for row in cb_passed.loc[:, ['orthogroup'] + bid_cols].itertuples(index=True, name=None):
        row_index = row[0]
        og = row[1]
        branch_id_set = frozenset(row[2:])
        analyzed_sets = analyzed_by_og.setdefault(og, [])
        if any(branch_id_set.issubset(analyzed_set) for analyzed_set in analyzed_sets):
            print(f'Skipped. Subset of already analyzed higher order convergence: {og} {branch_id_set}', flush=True)
            continue
        analyzed_sets.append(branch_id_set)
        keep_indices.append(row_index)
    cb_passed = cb_passed.loc[keep_indices, :].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity}: Skipped branch combinations that are subsets of already analyzed higher-order convergence: {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed, already_analyzed_in_greater_K

def generate_trait_colors(df_trait, trait_names):
    for trait_name in trait_names:
        df_trait_color = df_trait.loc[:,['species', trait_name]]
        is_foreground = (df_trait_color[trait_name]==1)
        df_trait_color.columns = ['species', 'color']
        df_trait_color['color'] = 'black'
        df_trait_color.loc[is_foreground,'color'] = 'firebrick'
        df_trait_color.to_csv(f'trait_{trait_name}.color.tsv', sep='\t', index=False)
    return None


def load_annotation_besthits(dir_of):
    path_annot = os.path.join(dir_of, 'Orthogroups', 'Orthogroups.GeneCount.annotated.tsv')
    annot = pandas.read_csv(path_annot, sep='\t', header=0, low_memory=False)
    annot.columns = annot.columns.str.replace('Orthogroup', 'orthogroup')
    cols = ['orthogroup'] + annot.columns[annot.columns.str.startswith('besthit_')].tolist()
    return annot.loc[:, cols]


def get_cb_required_columns(cb_columns, trait_names):
    required = {
        'orthogroup',
        'OCNany2spe',
        'ECNany2spe',
        'OCSany2spe',
        'ECSany2spe',
        'omegaCany2spe',
        'OCNCoD',
    }
    required.update(col for col in cb_columns if col.startswith('branch_id_'))
    if 'is_fg' in cb_columns:
        required.add('is_fg')
    else:
        required.update(f'is_fg_{trait}' for trait in trait_names if f'is_fg_{trait}' in cb_columns)
    if 'branch_num_fg_stem' in cb_columns:
        required.add('branch_num_fg_stem')
    else:
        required.update(
            f'branch_num_fg_stem_{trait}'
            for trait in trait_names
            if f'branch_num_fg_stem_{trait}' in cb_columns
        )
    return [col for col in cb_columns if col in required]

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
    dir_treevis = os.path.join(dir_myscript, 'treevis')
    file_stat_branch = os.path.join(dir_og, 'stat_branch', og+'.stat.branch.tsv')
    file_og_rpsblast = os.path.join(dir_og, 'rpsblast', og+'.rpsblast.tsv')
    file_og_alignment = get_alignment_for_tree_plot(dir_og=dir_og, og=og, dir_out_og=dir_out_og)
    file_csubst_input_fasta = os.path.join(dir_out_og, og+'.iqtree.anc', 'csubst.fasta')
    artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
    file_csubst_site_tsv = artifacts['site_table_tsv']
    if file_csubst_site_tsv is None:
        raise FileNotFoundError(
            f'CSUBST site table was not found in {artifacts["site_dir"]}.'
        )
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
    cmd.append('--treevis_dir='+dir_treevis)
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

def filter_max_per_og(cb_passed, arity, args):
    sampled_indices = []
    num_before_filtering = cb_passed.shape[0]
    if args.max_per_og <= 0:
        print(f'K = {arity}: Skipped branch combinations that are not top {args.max_per_og} within individual orthogroups: {num_before_filtering} -> {num_before_filtering}', flush=True)
        return cb_passed
    for og, og_indices in cb_passed.groupby('orthogroup', sort=False).indices.items():
        og_indices = numpy.asarray(og_indices, dtype=int)
        if (og_indices.shape[0] > args.max_per_og) and (args.max_per_og > 0):
            print(f'K = {arity}, {og}: Number of branch combinations ({og_indices.shape[0]}) exceeded --max_per_og. Only {args.max_per_og} combinations will be analyzed.')
            sampled_indice_positions = numpy.floor(numpy.linspace(0, og_indices.shape[0]-1, args.max_per_og)).astype(int)
            sampled_indices += og_indices[sampled_indice_positions].tolist()
        else:
            sampled_indices += og_indices.tolist()
    cb_passed = cb_passed.loc[sampled_indices,:].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(f'K = {arity}: Skipped branch combinations that are not top {args.max_per_og} within individual orthogroups: {num_before_filtering} -> {num_after_filtering}', flush=True)
    return cb_passed

def filter_fg_stem_ratio(cb_passed, arity, trait, args, no_trait_name):
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

def write_annotated_table(cb_passed, annot, dir_out, out_name):
    cb_passed = pandas.merge(cb_passed, annot, how='left', on='orthogroup', sort=False)
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
    require_pdf_dependencies()
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


def resolve_workspace_dir(workspace_dir_arg):
    if workspace_dir_arg not in ('', 'auto'):
        return os.path.realpath(workspace_dir_arg)

    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, 'input')) or os.path.isdir(os.path.join(cwd, 'output')):
        return os.path.realpath(cwd)

    workspace_candidate = os.path.join(cwd, 'workspace')
    if os.path.isdir(os.path.join(workspace_candidate, 'input')) or os.path.isdir(os.path.join(workspace_candidate, 'output')):
        return os.path.realpath(workspace_candidate)

    return os.path.realpath(workspace_candidate)


def resolve_path_arg(path_arg, workspace_dir, *relative_parts):
    if path_arg not in ('', 'auto'):
        return os.path.realpath(path_arg)
    return os.path.realpath(os.path.join(workspace_dir, *relative_parts))

if __name__ == '__main__':
    print(f'{datetime.datetime.now()}: {__file__} started.', flush=True)

    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace-dir', metavar='PATH', default='auto', type=str, help='')
    parser.add_argument('--dir_orthogroup', metavar='PATH', default='auto', type=str, help='')
    parser.add_argument('--dir_orthofinder', metavar='PATH', default='auto', type=str, help='')
    parser.add_argument('--dir_out', metavar='PATH', default='auto', type=str, help='')
    parser.add_argument('--file_trait', metavar='PATH', default='auto', type=str, help='')
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
    workspace_dir = resolve_workspace_dir(args.workspace_dir)
    dir_og = resolve_path_arg(args.dir_orthogroup, workspace_dir, 'output', 'orthogroup')
    dir_of = resolve_path_arg(args.dir_orthofinder, workspace_dir, 'output', 'orthofinder')
    args.dir_out = resolve_path_arg(args.dir_out, workspace_dir, 'output', 'csubst_site')
    args.file_trait = resolve_path_arg(args.file_trait, workspace_dir, 'input', 'species_trait', 'species_trait.tsv')
    db_file = os.path.join(dir_og, 'gg_orthogroup.db')
    arity_min = int(args.arity_range.split('-')[0])
    arity_max = int(args.arity_range.split('-')[1])
    arity_range = numpy.arange(arity_min, arity_max+1)[::-1]
    print('Arity range: {}'.format(' '.join([ str(a) for a in arity_range ])), flush=True)
    require_sqlalchemy()
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
    df_trait = pandas.read_csv(args.file_trait, sep='\t', header=0, index_col=None)
    trait_names = [t.strip() for t in args.trait.split(',') if t.strip()]
    no_trait_name = False
    if (len(trait_names) == 1) and (trait_names[0] == 'all'):
        detected_trait_names = []
        for arity in arity_range:
            if 'cb' + str(arity) in column_names.keys():
                colnames = column_names['cb' + str(arity)]
                detected_trait_names = [colname.replace('is_fg_', '') for colname in colnames if colname.startswith('is_fg_')]
                if len(detected_trait_names) > 0:
                    break
        if len(detected_trait_names) > 0:
            trait_names = detected_trait_names
            no_trait_name = False
        else:
            # Fallback when cb* tables are not generated yet.
            trait_names = [c for c in df_trait.columns.tolist() if c != 'species']
            no_trait_name = True
            print(
                f'No trait name found in cb tables. Using all traits from {args.file_trait}: {",".join(trait_names)}',
                flush=True
            )
    else:
        missing_trait_cols = [t for t in trait_names if t not in df_trait.columns]
        if len(missing_trait_cols) > 0:
            raise ValueError(f'Trait columns not found in {args.file_trait}: {missing_trait_cols}')
        no_trait_name = False
    if len(trait_names) == 0:
        raise ValueError(f'No analyzable trait columns were detected in {args.file_trait}.')
    print(f'Specified or detected traits: {",".join(trait_names)}', flush=True)
    generate_trait_colors(df_trait=df_trait, trait_names=trait_names)
    annot_besthits = load_annotation_besthits(dir_of)
    already_analyzed_in_greater_K = dict()
    for trait in trait_names:
        already_analyzed_in_greater_K[trait] = dict() # Initialize
    for arity in arity_range:
        print('')
        if 'cb'+str(arity) not in table_names:
            print(f'cb{arity} is not in the database. Skipping.', flush=True)
            continue
        cb_table = 'cb'+str(arity)
        cb_required_cols = get_cb_required_columns(column_names[cb_table], trait_names)
        select_cols_sql = ', '.join('"{}"'.format(col) for col in cb_required_cols)
        sql_cb = sqlalchemy.text(f'SELECT {select_cols_sql} from {cb_table}')
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
            if (not os.path.exists(dir_out)) and (not os.path.exists(out_zip)):
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
                if os.path.exists(dir_out) and (not os.path.exists(out_zip)):
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
                out = skip_lower_order(cb_passed, arity, trait, already_analyzed_in_greater_K)
                cb_passed, already_analyzed_in_greater_K = out
                del out
            if args.min_fg_stem_ratio > 0:
                cb_passed = filter_fg_stem_ratio(cb_passed, arity, trait, args, no_trait_name)
            cb_passed = write_annotated_table(cb_passed, annot_besthits, dir_out, out_name)
            cb_passed = filter_max_per_og(cb_passed, arity, args)
            if args.max_per_K < cb_passed.shape[0]:
                cb_passed = filter_max_per_K_1st(cb_passed, arity, trait)
            if args.max_per_K < cb_passed.shape[0]:
                cb_passed = filter_max_per_K_2nd(cb_passed, arity, trait)
            print('Starting analyzing individual branch combinations.', flush=True)
            branch_id_cols = cb_passed.columns[cb_passed.columns.str.startswith('branch_id_')].tolist()
            og_list = cb_passed['orthogroup'].tolist()
            branch_id_str_list = [
                ','.join(str(bid) for bid in branch_ids)
                for branch_ids in cb_passed.loc[:, branch_id_cols].itertuples(index=False, name=None)
            ]
            besthit_cols = ['besthit_0.05', 'besthit_0.25', 'besthit_0.5', 'besthit_0.75', 'besthit_0.95']
            besthit_frame = cb_passed.reindex(columns=besthit_cols)
            annotation_text_list = []
            for og, branch_id_str, besthit_values_raw in zip(
                og_list,
                branch_id_str_list,
                besthit_frame.itertuples(index=False, name=None),
            ):
                besthit_values = [
                    str(value) if pandas.notna(value) else 'NA'
                    for value in besthit_values_raw
                ]
                annotation_text_list.append(
                    get_annotation_text(
                        og=og,
                        arity=arity,
                        branch_id_str=branch_id_str,
                        trait=trait,
                        min_OCNany2spe=min_OCNany2spe,
                        min_omegaCany2spe=min_omegaCany2spe,
                        min_OCNCoD=min_OCNCoD,
                        besthit_values=besthit_values,
                    )
                )
            with ProcessPoolExecutor(max_workers=args.ncpu) as executor:
                results = executor.map(
                    process_index,
                    og_list,
                    branch_id_str_list,
                    repeat(dir_out),
                    repeat(dir_og),
                    repeat(file_trait_color),
                    repeat(args.ncpu),
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
