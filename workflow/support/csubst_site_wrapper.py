#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import fcntl
import glob
import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import uuid
import zipfile
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import numpy
import pandas

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gene_family_output_store import GeneFamilyOutputStore
from safe_zip_extract import extract_expected_prefix

CSUBST_SITE_ARCHIVE_COMMENT = b"GeneGalleon csubst-site archive v1\n"

try:
    import sqlalchemy
except ImportError:
    sqlalchemy = None


def _fsync_directory(path):
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def csubst_site_archive_is_complete(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        return False
    try:
        with zipfile.ZipFile(path, "r") as archive:
            if archive.comment != CSUBST_SITE_ARCHIVE_COMMENT:
                return False
            if not any(not info.is_dir() for info in archive.infolist()):
                return False
            return archive.testzip() is None
    except (OSError, RuntimeError, zipfile.BadZipFile):
        return False


def quarantine_unverified_csubst_site_archive(path):
    path = Path(path)
    if not path.exists():
        return None
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Invalid CSUBST site archive path: {path}")
    quarantine = path.with_name(f".{path.name}.unverified.{os.getpid()}.{uuid.uuid4().hex}")
    os.replace(path, quarantine)
    _fsync_directory(path.parent)
    return quarantine


def cleanup_quarantined_csubst_site_archives(path):
    path = Path(path)
    removed = False
    for candidate in path.parent.glob(f".{path.name}.unverified.*"):
        if candidate.is_file() and not candidate.is_symlink():
            candidate.unlink()
            removed = True
    if removed:
        _fsync_directory(path.parent)


def create_csubst_site_archive(directory):
    directory = Path(directory).resolve()
    if directory.is_symlink() or not directory.is_dir():
        raise RuntimeError(f"CSUBST site archive source is not a directory: {directory}")
    final_path = Path(f"{directory}.zip")
    if final_path.is_symlink():
        raise RuntimeError(f"Symlinked CSUBST site ZIP is unsupported: {final_path}")
    temporary_base = final_path.parent / (f".{final_path.stem}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    temporary_zip = Path(f"{temporary_base}.zip")
    try:
        shutil.make_archive(
            str(temporary_base),
            "zip",
            root_dir=str(directory),
        )
        with zipfile.ZipFile(temporary_zip, "a") as archive:
            archive.comment = CSUBST_SITE_ARCHIVE_COMMENT
        if not csubst_site_archive_is_complete(temporary_zip):
            raise RuntimeError(f"New CSUBST site archive failed verification: {temporary_zip}")
        with temporary_zip.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary_zip, final_path)
        _fsync_directory(final_path.parent)
        return final_path
    finally:
        temporary_zip.unlink(missing_ok=True)


class LockedMaterializationDirectory:
    def __init__(self, parent, og):
        self.parent = Path(parent).resolve()
        self.root = self.parent / ".gg_materialized"
        # One stable lock outside the per-analysis directory avoids including a
        # lock file in the result ZIP. It must not be unlinked: unlinking a
        # flocked file lets a concurrent opener lock a different inode.
        self.global_lock_path = self.parent.parent / ".gg_materialized.lock"
        self.descriptor = None
        self.name = ""
        self.parent.mkdir(parents=True, exist_ok=True)
        global_descriptor = _open_materialization_lock(self.global_lock_path)
        try:
            fcntl.flock(global_descriptor, fcntl.LOCK_EX)
            self.root.mkdir(parents=True, exist_ok=True)
            for candidate in sorted(self.root.iterdir()):
                if not candidate.is_dir() or candidate.is_symlink():
                    continue
                lock_path = candidate / ".run.lock"
                stale_descriptor = _open_materialization_lock(lock_path)
                try:
                    try:
                        fcntl.flock(
                            stale_descriptor,
                            fcntl.LOCK_EX | fcntl.LOCK_NB,
                        )
                    except BlockingIOError:
                        continue
                    shutil.rmtree(candidate)
                finally:
                    os.close(stale_descriptor)
            run_dir = Path(
                tempfile.mkdtemp(
                    prefix=f"{og}.",
                    dir=self.root,
                )
            )
            lock_path = run_dir / ".run.lock"
            self.descriptor = _open_materialization_lock(lock_path)
            fcntl.flock(self.descriptor, fcntl.LOCK_EX)
            self.name = str(run_dir)
        finally:
            fcntl.flock(global_descriptor, fcntl.LOCK_UN)
            os.close(global_descriptor)

    def cleanup(self):
        if self.descriptor is None:
            return
        global_descriptor = _open_materialization_lock(self.global_lock_path)
        try:
            fcntl.flock(global_descriptor, fcntl.LOCK_EX)
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None
            if self.name:
                shutil.rmtree(self.name, ignore_errors=True)
            try:
                self.root.rmdir()
            except OSError:
                pass
        finally:
            fcntl.flock(global_descriptor, fcntl.LOCK_UN)
            os.close(global_descriptor)


def _open_materialization_lock(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_RDWR | os.O_CREAT
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    return os.open(path, flags, 0o600)


def cleanup_materialization_metadata(parent):
    parent = Path(parent).resolve()
    root = parent / ".gg_materialized"
    try:
        root.rmdir()
    except FileNotFoundError:
        pass
    except OSError:
        return


try:
    from distutils.util import strtobool
except ModuleNotFoundError:

    def strtobool(val):
        val = str(val).strip().lower()
        if val in ("y", "yes", "t", "true", "on", "1"):
            return 1
        if val in ("n", "no", "f", "false", "off", "0"):
            return 0
        raise ValueError("invalid truth value {!r}".format(val))


try:
    from pypdf import PdfReader, PdfWriter
except ImportError:
    PdfReader = None
    PdfWriter = None
try:
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas
except ImportError:
    canvas = None
    letter = None
    pdfmetrics = None
    TTFont = None

CSUBST_NONSYN_RECODE_CHOICES = (
    "no",
    "3di20",
    "dayhoff6",
    "sr6",
    "kgb6",
    "sr4",
    "dayhoff9",
    "dayhoff12",
    "dayhoff15",
    "dayhoff18",
    "srchisq6",
    "kgbauto6",
)
CSUBST_SITES_OUTDIR = "csubst_sites"
CSUBST_SITES_OUTPUT_PREFIX = "csubst"
NSY_ALIGNMENT_SYMBOLS = tuple("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz")
STANDARD_GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def normalize_csubst_nonsyn_recode(value):
    normalized = str(value or "no").strip().lower()
    if normalized not in CSUBST_NONSYN_RECODE_CHOICES:
        allowed = ", ".join(CSUBST_NONSYN_RECODE_CHOICES)
        raise ValueError(f"Invalid csubst_nonsyn_recode: {value}. Expected one of: {allowed}")
    return normalized


def csubst_nonsyn_recode_output_suffix(value):
    recode = normalize_csubst_nonsyn_recode(value)
    if recode == "no":
        return ""
    return f"_nonsynRecode-{recode}"


def build_csubst_sites_command(
    iqtree_anc_rel_dir,
    iqtree_anc_dir,
    branch_id_str,
    ncpu,
    csubst_nonsyn_recode,
    pdb="besthit",
):
    recode = normalize_csubst_nonsyn_recode(csubst_nonsyn_recode)
    cmd = ["csubst", "sites"]
    cmd += ["--outdir", CSUBST_SITES_OUTDIR]
    cmd += ["--output_prefix", CSUBST_SITES_OUTPUT_PREFIX]
    cmd += ["--alignment_file", os.path.join(iqtree_anc_rel_dir, "csubst.fasta")]
    cmd += ["--rooted_tree_file", os.path.join(iqtree_anc_rel_dir, "csubst.nwk")]
    cmd += ["--branch_id", branch_id_str]
    cmd += ["--threads", str(max(1, int(ncpu)))]
    cmd += ["--iqtree_treefile", os.path.join(iqtree_anc_dir, "csubst.treefile")]
    cmd += ["--iqtree_state", os.path.join(iqtree_anc_dir, "csubst.state")]
    cmd += ["--iqtree_rate", os.path.join(iqtree_anc_dir, "csubst.rate")]
    cmd += ["--iqtree_iqtree", os.path.join(iqtree_anc_dir, "csubst.iqtree")]
    cmd += ["--iqtree_log", os.path.join(iqtree_anc_dir, "csubst.log")]
    if recode != "no":
        cmd += ["--nonsyn_recode", recode]
    if str(pdb or "none").strip().lower() != "none":
        cmd += ["--pdb", str(pdb)]
    return cmd


def shell_join_command(cmd):
    return " ".join(shlex.quote(str(token)) for token in cmd)


pandas.set_option("display.max_columns", None)


def _get_matplotlib():
    import matplotlib
    import matplotlib.patches
    import matplotlib.pyplot

    matplotlib.rcParams["font.size"] = 8
    matplotlib.rcParams["font.family"] = "Helvetica"
    matplotlib.rcParams["svg.fonttype"] = "none"  # none, path, or svgfont
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
    fig, axes = matplotlib.pyplot.subplots(
        nrows=2, ncols=2, figsize=(3.2, 3.2), sharey=False, sharex=False, width_ratios=wr, height_ratios=hr
    )
    axes = axes.flat
    axes[1].axis("off")
    ax = axes[0]
    ax.axvline(x=min_OCNany2spe, color="red", alpha=0.3)
    ax.hist(x=cb_plot[xcol], bins=numpy.arange(0, 1.001 + (1 / nbin), 1 / nbin) * xmax, color="darkblue")
    ax.set_xlim(0, xmax * axis_margin_rate)
    ax.set_xticks([])
    ax.set_yscale("log")
    ax.set_ylabel("Count")
    ax = axes[3]
    ax.axhline(y=min_omegaCany2spe, color="red", alpha=0.3)
    ax.hist(
        x=cb_plot[ycol],
        bins=numpy.arange(0, 1.001 + (1 / nbin), 1 / nbin) * ymax,
        color="darkblue",
        orientation="horizontal",
    )
    ax.set_ylim(0, ymax * axis_margin_rate)
    ax.set_yticks([])
    ax.set_xscale("log")
    ax.set_xlabel("Count")
    ax = axes[2]
    plot_range = [0, xmax * axis_margin_rate, 0, ymax * axis_margin_rate]
    polygon_coord = [[polygon_xmin, polygon_ymin], [xmax, polygon_ymin], [xmax, ymax], [polygon_xmin, ymax]]
    ax.add_patch(matplotlib.patches.Polygon(polygon_coord, color="red", alpha=0.3))
    ax.hexbin(
        cb_plot[xcol].values,
        cb_plot[ycol].values,
        gridsize=nbin,
        extent=plot_range,
        cmap="jet",
        mincnt=1,
    )
    ax.set_xlim(0, xmax * axis_margin_rate)
    ax.set_ylim(0, ymax * axis_margin_rate)
    fig.tight_layout(pad=0.25, w_pad=1.5, h_pad=0.5)
    for ext in ["svg", "pdf"]:
        outpath = os.path.join(outbase + "." + ext)
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
    pattern = re.compile(r"^csubst\.([^.]+)\.fa$")
    if not os.path.isdir(indir):
        print(f"PDB ID search directory not found: {indir}", flush=True)
        return None
    for filename in sorted(f for f in os.listdir(indir) if not f.startswith(".")):
        match = pattern.search(filename)
        if match:
            pdb_id = match.group(1)
            print(f"Found file: {filename} --> Extracted PDB ID: {pdb_id}", flush=True)
            return pdb_id
    print(f"PDB ID not found in: {indir}", flush=True)
    return None


def resolve_site_output_dir(dir_out_og, branch_id_str):
    return os.path.join(
        dir_out_og,
        CSUBST_SITES_OUTDIR,
        CSUBST_SITES_OUTPUT_PREFIX + ".branch_id" + branch_id_str,
    )


def load_site_output_manifest(site_dir):
    manifest_path = os.path.join(site_dir, CSUBST_SITES_OUTPUT_PREFIX + ".outputs.tsv")
    if not os.path.isfile(manifest_path):
        return None
    return pandas.read_csv(manifest_path, sep="\t", header=0, index_col=None)


def resolve_manifest_output_path(site_dir, row):
    output_path = str(row.get("output_path", "")).strip()
    if output_path not in ("", "nan"):
        if os.path.isabs(output_path):
            return output_path
        return os.path.join(site_dir, output_path)
    output_file = str(row.get("output_file", "")).strip()
    if output_file not in ("", "nan"):
        return os.path.join(site_dir, output_file)
    return None


def find_site_output_from_manifest(manifest_df, site_dir, output_kinds):
    if manifest_df is None:
        return None
    if "output_kind" not in manifest_df.columns:
        return None
    if isinstance(output_kinds, str):
        output_kinds = [output_kinds]
    manifest_df = manifest_df.copy()
    if "file_exists" in manifest_df.columns:
        is_existing = manifest_df["file_exists"].astype(str).str.upper().isin({"Y", "YES", "TRUE", "1"})
        manifest_df = manifest_df.loc[is_existing, :].reset_index(drop=True)
    for output_kind in output_kinds:
        rows = manifest_df.loc[manifest_df["output_kind"].astype(str) == str(output_kind), :]
        for _, row in rows.iterrows():
            output_path = resolve_manifest_output_path(site_dir=site_dir, row=row)
            # gg-cache-guard: audited - this resolves a manifest entry; it does not reuse a computation.
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
        basename.endswith(".outputs.tsv")
        or basename.endswith(".state_N.tsv")
        or basename.endswith(".state_S.tsv")
        or ".tree_site." in basename
    )


def is_auxiliary_site_pdf(path):
    basename = os.path.basename(path)
    return basename.endswith(".pymol.pdf") or basename.endswith(".state.pdf") or ".tree_site." in basename


def resolve_site_artifacts(dir_out_og, branch_id_str):
    site_dir = resolve_site_output_dir(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
    manifest_df = load_site_output_manifest(site_dir=site_dir)
    pdb_id = extract_pdb_id(indir=site_dir)

    site_table_candidates = []
    site_summary_pdf_candidates = []
    pymol_pdf_candidates = []
    if pdb_id is not None:
        site_table_candidates.extend(
            [
                os.path.join(site_dir, f"{CSUBST_SITES_OUTPUT_PREFIX}.{pdb_id}.tsv"),
            ]
        )
        site_summary_pdf_candidates.extend(
            [
                os.path.join(site_dir, f"{CSUBST_SITES_OUTPUT_PREFIX}.{pdb_id}.pdf"),
            ]
        )
        pymol_pdf_candidates.extend(
            [
                os.path.join(site_dir, f"{CSUBST_SITES_OUTPUT_PREFIX}.{pdb_id}.pymol.pdf"),
            ]
        )
    site_table_candidates.extend(
        [
            os.path.join(site_dir, CSUBST_SITES_OUTPUT_PREFIX + ".tsv"),
        ]
    )
    site_summary_pdf_candidates.extend(
        [
            os.path.join(site_dir, CSUBST_SITES_OUTPUT_PREFIX + ".pdf"),
        ]
    )

    site_table_tsv = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=["site_table_tsv"],
    )
    if site_table_tsv is None:
        site_table_tsv = resolve_existing_path(site_table_candidates)
    if site_table_tsv is None:
        site_table_tsv = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=[CSUBST_SITES_OUTPUT_PREFIX + "*.tsv"],
            exclude_predicate=is_auxiliary_site_tsv,
        )

    site_summary_pdf = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=["site_summary_pdf"],
    )
    if site_summary_pdf is None:
        site_summary_pdf = resolve_existing_path(site_summary_pdf_candidates)
    if site_summary_pdf is None:
        site_summary_pdf = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=[CSUBST_SITES_OUTPUT_PREFIX + "*.pdf"],
            exclude_predicate=is_auxiliary_site_pdf,
        )

    pymol_summary_pdf = find_site_output_from_manifest(
        manifest_df=manifest_df,
        site_dir=site_dir,
        output_kinds=["pymol_summary_pdf"],
    )
    if pymol_summary_pdf is None:
        pymol_summary_pdf = resolve_existing_path(pymol_pdf_candidates)
    if pymol_summary_pdf is None:
        pymol_summary_pdf = find_first_matching_file(
            site_dir=site_dir,
            glob_patterns=[CSUBST_SITES_OUTPUT_PREFIX + "*.pymol.pdf"],
        )

    return {
        "site_dir": site_dir,
        "manifest_df": manifest_df,
        "pdb_id": pdb_id,
        "site_table_tsv": site_table_tsv,
        "site_summary_pdf": site_summary_pdf,
        "pymol_summary_pdf": pymol_summary_pdf,
    }


def format_annotation_value(value):
    if value is None:
        return "NA"
    try:
        if pandas.isna(value):
            return "NA"
    except TypeError:
        pass
    return str(value)


def normalize_observed_branch_stats(observed_values):
    keys = [
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "omegaCany2any",
        "omegaCdif2spe",
        "OCNCoD",
    ]
    observed = {key: "NA" for key in keys}
    if isinstance(observed_values, dict):
        for key in keys:
            observed[key] = format_annotation_value(observed_values.get(key))
    elif observed_values is not None:
        for key, value in zip(keys, observed_values, strict=False):
            observed[key] = format_annotation_value(value)
    return observed


def append_csubst_sites_command(annotation_text, cmd):
    if not cmd:
        return annotation_text
    return f"""{annotation_text}

CSUBST sites command

{shell_join_command(cmd)}
"""


def get_annotation_text(
    og, arity, branch_id_str, trait, min_OCNany2spe, min_omegaCany2spe, min_OCNCoD, besthit_values, observed_values=None
):
    if besthit_values is None:
        besthit_values = ["NA", "NA", "NA", "NA", "NA"]
    if len(besthit_values) < 5:
        besthit_values = list(besthit_values) + (["NA"] * (5 - len(besthit_values)))
    observed = normalize_observed_branch_stats(observed_values)
    annotation_text = f"""Orthogroup: {og}

K: {arity}

Branch ID: {branch_id_str}

Trait: {trait}

Threshold

OCNany2spe: {min_OCNany2spe}

omegaCany2spe: {min_omegaCany2spe}

OCNCoD: {min_OCNCoD}

Observed

OCNany2spe: {observed["OCNany2spe"]}

ECNany2spe: {observed["ECNany2spe"]}

OCSany2spe: {observed["OCSany2spe"]}

ECSany2spe: {observed["ECSany2spe"]}

omegaCany2spe: {observed["omegaCany2spe"]}

omegaCany2any: {observed["omegaCany2any"]}

omegaCdif2spe: {observed["omegaCdif2spe"]}

OCNCoD: {observed["OCNCoD"]}

Time: {datetime.datetime.now()}

Annotation of the gene as being in the 5th, 25th, 50th, 75th, or 95th percentile for length within its orthogroup:

{besthit_values[0]}

{besthit_values[1]}

{besthit_values[2]}

{besthit_values[3]}

{besthit_values[4]}
"""
    return annotation_text


def process_index(og, branch_id_str, dir_out, dir_og, file_trait_color, ncpu, csubst_nonsyn_recode, annotation_text):
    previous_cwd = os.getcwd()
    dir_out_og = os.path.join(dir_out, og + "_" + branch_id_str.replace(",", "_"))
    print("{}, --branch_ids {}: wd: {}".format(og, branch_id_str, dir_out_og), flush=True)
    if not os.path.exists(dir_out_og):
        os.makedirs(dir_out_og)
    os.chdir(dir_out_og)
    iqtree_anc_dir = get_iqtree_anc_dir(dir_out_og=dir_out_og, og=og)
    iqtree_anc_rel_dir = os.path.basename(iqtree_anc_dir)
    csubst_sites_cmd = build_csubst_sites_command(
        iqtree_anc_rel_dir=iqtree_anc_rel_dir,
        iqtree_anc_dir=iqtree_anc_dir,
        branch_id_str=branch_id_str,
        ncpu=ncpu,
        csubst_nonsyn_recode=csubst_nonsyn_recode,
    )
    iqtree_tree_file = os.path.join(iqtree_anc_dir, "csubst.treefile")
    iqtree_state_file = os.path.join(iqtree_anc_dir, "csubst.state")
    iqtree_rate_file = os.path.join(iqtree_anc_dir, "csubst.rate")
    iqtree_iqtree_file = os.path.join(iqtree_anc_dir, "csubst.iqtree")
    iqtree_log_file = os.path.join(iqtree_anc_dir, "csubst.log")
    iqtree_ckp_file = os.path.join(iqtree_anc_dir, "csubst.ckp.gz")
    file_summary = os.path.join(dir_out_og, f"summary.{og}_branch_id{branch_id_str}.pdf")
    # gg-cache-guard: audited - the outer csubst_site artifact contract removes this directory on rebuild.
    if os.path.exists(file_summary):
        print(f"Skipped. Outfile already exists: {file_summary}", flush=True)
        os.chdir(previous_cwd)
        return og, None
    materialized_input_context = None
    try:
        effective_dir_og = dir_og
        if any(os.path.isdir(os.path.join(dir_og, name)) for name in (".gg_store", ".gg_archives")):
            materialized_input_context = LockedMaterializationDirectory(
                dir_out,
                og,
            )
            effective_dir_og = materialized_input_context.name
            materialize_csubst_site_inputs(
                dir_og=dir_og,
                og=og,
                destination_root=effective_dir_og,
            )
        artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
        file_csubst_out = artifacts["site_table_tsv"]
        if file_csubst_out is not None and os.path.exists(file_csubst_out):
            print(f"Skipped csubst sites. Outfile already exists: {file_csubst_out}", flush=True)
        else:
            print(f"Running csubst sites. Output file not found: {file_csubst_out}", flush=True)
            path_iqtree_zip = get_iqtree_anc_zip_path(dir_og=effective_dir_og, og=og)
            extract_expected_prefix(
                path_iqtree_zip,
                dir_out_og,
                f"{og}.iqtree.anc",
            )
            cmd = csubst_sites_cmd
            print("COMMAND: {}".format(" ".join(cmd)), flush=True)
            subprocess.run(cmd, check=True)
        print(f"{datetime.datetime.now()}: csubst sites done: {og}", flush=True)
        file_tree_plot = os.path.join(dir_out_og, og + ".tree_plot.pdf")
        if os.path.exists(file_tree_plot):
            print(f"Tree plot skipped: outfile already exists: {file_tree_plot}", flush=True)
        else:
            print(f"Running stat_branch2tree_plot to generate a tree pdf file: {file_tree_plot}", flush=True)
            run_stat_branch2tree_plot(
                og,
                branch_id_str,
                file_trait_color=file_trait_color,
                dir_out_og=dir_out_og,
                dir_og=effective_dir_og,
                ncpu=ncpu,
                csubst_nonsyn_recode=csubst_nonsyn_recode,
            )
            print(f"{datetime.datetime.now()}: stat_branch2tree_plot done: {og}", flush=True)
        if os.path.exists("annotation_text.pdf"):
            print(
                f"Annotation text generation skipped: outfile already exists: annotation_text.pdf for {og}", flush=True
            )
        else:
            print(f"Generating annotation text: annotation_text.pdf for {og}", flush=True)
            create_pdf(append_csubst_sites_command(annotation_text, csubst_sites_cmd), "annotation_text.pdf")
        artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
        if artifacts["site_summary_pdf"] is None:
            raise FileNotFoundError(f"CSUBST site-summary PDF was not found in {artifacts['site_dir']}.")
        pdf_file_paths = [artifacts["site_summary_pdf"]]
        if artifacts["pymol_summary_pdf"] is not None and os.path.exists(artifacts["pymol_summary_pdf"]):
            pdf_file_paths.append(artifacts["pymol_summary_pdf"])
        pdf_file_paths.extend([file_tree_plot, "annotation_text.pdf"])
        combine_pdfs(pdf_file_paths, output_path=file_summary)
        remove_files = [
            "tmp.csubst.sub_tensor.S.mmap",
            "tmp.csubst.sub_tensor.N.mmap",
            "Rplots.pdf",
            "species_color.tsv",
            "annotation_text.pdf",
            file_tree_plot,
            iqtree_tree_file,
            iqtree_state_file,
            iqtree_rate_file,
            iqtree_iqtree_file,
            iqtree_log_file,
            iqtree_ckp_file,
        ] + glob.glob("*.cif")
        if os.path.exists(file_summary):
            for remove_file in remove_files:
                if os.path.exists(remove_file):
                    os.remove(remove_file)
        return og, None
    except Exception as e:
        print(f"process_index: Error in {og}: {e}", flush=True)
        return og, e
    finally:
        os.chdir(previous_cwd)
        if materialized_input_context is not None:
            materialized_input_context.cleanup()


def process_og_batch(
    og,
    branch_id_strs,
    dir_out,
    dir_og,
    file_trait_color,
    ncpu,
    csubst_nonsyn_recode,
    annotation_texts,
):
    materialized_input_context = None
    effective_dir_og = dir_og
    try:
        if any(os.path.isdir(os.path.join(dir_og, name)) for name in (".gg_store", ".gg_archives")):
            materialized_input_context = LockedMaterializationDirectory(
                dir_out,
                og,
            )
            effective_dir_og = materialized_input_context.name
            materialize_csubst_site_inputs(
                dir_og=dir_og,
                og=og,
                destination_root=effective_dir_og,
            )
        return [
            process_index(
                og,
                branch_id_str,
                dir_out,
                effective_dir_og,
                file_trait_color,
                ncpu,
                csubst_nonsyn_recode,
                annotation_text,
            )
            for branch_id_str, annotation_text in zip(
                branch_id_strs,
                annotation_texts,
                strict=True,
            )
        ]
    finally:
        if materialized_input_context is not None:
            materialized_input_context.cleanup()


def skip_lower_order(cb_passed, arity, trait, already_analyzed_in_greater_K):
    num_before_filtering = cb_passed.shape[0]
    bid_cols = cb_passed.columns[cb_passed.columns.str.startswith("branch_id_")].tolist()
    keep_indices = []
    analyzed_by_og = already_analyzed_in_greater_K[trait]
    for row in cb_passed.loc[:, ["orthogroup"] + bid_cols].itertuples(index=True, name=None):
        row_index = row[0]
        og = row[1]
        branch_id_set = frozenset(row[2:])
        analyzed_sets = analyzed_by_og.setdefault(og, [])
        if any(branch_id_set.issubset(analyzed_set) for analyzed_set in analyzed_sets):
            print(f"Skipped. Subset of already analyzed higher order convergence: {og} {branch_id_set}", flush=True)
            continue
        analyzed_sets.append(branch_id_set)
        keep_indices.append(row_index)
    cb_passed = cb_passed.loc[keep_indices, :].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(
        f"K = {arity}: Skipped branch combinations that are subsets of already analyzed higher-order convergence: {num_before_filtering} -> {num_after_filtering}",
        flush=True,
    )
    return cb_passed, already_analyzed_in_greater_K


def is_foreground_trait_value(value):
    if pandas.isna(value):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip()
    if text == "" or text.lower() in {"na", "nan", "none", "null"}:
        return False
    try:
        return float(text) != 0
    except ValueError:
        return True


def generate_trait_colors(df_trait, trait_names):
    for trait_name in trait_names:
        df_trait_color = df_trait.loc[:, ["species", trait_name]]
        is_foreground = df_trait_color[trait_name].map(is_foreground_trait_value)
        df_trait_color.columns = ["species", "color"]
        df_trait_color["color"] = "black"
        df_trait_color.loc[is_foreground, "color"] = "firebrick"
        df_trait_color.to_csv(f"trait_{trait_name}.color.tsv", sep="\t", index=False)
    return None


def load_annotation_besthits(dir_of):
    path_annot = os.path.join(dir_of, "Orthogroups_filtered", "Orthogroups.GeneCount.annotated.tsv")
    annot = pandas.read_csv(path_annot, sep="\t", header=0, low_memory=False)
    annot.columns = annot.columns.str.replace("Orthogroup", "orthogroup")
    cols = ["orthogroup"] + annot.columns[annot.columns.str.startswith("besthit_")].tolist()
    return annot.loc[:, cols]


def get_cb_required_columns(cb_columns, trait_names):
    required = {
        "orthogroup",
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "OCNCoD",
    }
    required.update(col for col in ["omegaCany2any", "omegaCdif2spe"] if col in cb_columns)
    required.update(col for col in cb_columns if col.startswith("branch_id_"))
    if "is_fg" in cb_columns:
        required.add("is_fg")
    else:
        required.update(f"is_fg_{trait}" for trait in trait_names if f"is_fg_{trait}" in cb_columns)
    if "branch_num_fg_stem" in cb_columns:
        required.add("branch_num_fg_stem")
    else:
        required.update(
            f"branch_num_fg_stem_{trait}" for trait in trait_names if f"branch_num_fg_stem_{trait}" in cb_columns
        )
    return [col for col in cb_columns if col in required]


def resolve_existing_path(candidates):
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return None


def get_iqtree_anc_zip_path(dir_og, og):
    return os.path.join(dir_og, "iqtree_anc", og + "_iqtree.anc.zip")


def get_iqtree_anc_dir(dir_out_og, og):
    return os.path.join(dir_out_og, og + ".iqtree.anc")


def get_stat_branch_path(dir_og, og):
    return os.path.join(dir_og, "stat_branch", og + "_stat.branch.tsv")


def get_rpsblast_path(dir_og, og):
    return os.path.join(dir_og, "rpsblast", og + "_rpsblast.tsv")


def get_csubst_site_input_candidates(og):
    return {
        "iqtree_anc": [og + "_iqtree.anc.zip"],
        "stat_branch": [og + "_stat.branch.tsv"],
        "rpsblast": [og + "_rpsblast.tsv"],
        "clipkit": [
            og + "_cds.clipkit.fa.gz",
            og + "_cds.clipkit.fasta",
            og + "_cds.clipkit.fa",
        ],
        "mafft": [
            og + "_cds.aln.fa.gz",
            og + "_cds.aln.fasta",
            og + "_cds.aln.fa",
        ],
        "orthogroup_extraction_fasta": [
            og + "_orthogroup_extraction.fa.gz",
            og + "_orthogroup_extraction.fasta",
            og + "_orthogroup_extraction.fa",
        ],
        "cds_fasta": [
            og + "_cds.fa.gz",
            og + "_cds.fasta",
            og + "_cds.fa",
        ],
    }


def _parse_branch_id_list(value):
    tokens = [token.strip() for token in str(value).split(",") if token.strip()]
    if len(tokens) == 0:
        raise ValueError("Branch ID list is empty.")
    branch_ids = []
    for token in tokens:
        if re.fullmatch(r"[0-9]+", token) is None:
            raise ValueError(f"Invalid branch ID: {token}")
        branch_ids.append(int(token))
    return sorted(set(branch_ids))


def get_csubst_branch_clade_signatures(iqtree_anc_dir):
    """Return CSUBST branch IDs keyed to stable descendant-tip signatures."""
    try:
        from csubst import ete as csubst_ete
        from csubst import parser_misc as csubst_parser_misc
        from csubst import tree as csubst_tree
    except ImportError as error:
        raise RuntimeError("The csubst Python package is required to validate CSUBST branch IDs.") from error

    context = {
        "rooted_tree_file": os.path.join(iqtree_anc_dir, "csubst.nwk"),
        "iqtree_treefile": os.path.join(iqtree_anc_dir, "csubst.treefile"),
    }
    for path in context.values():
        if not os.path.isfile(path):
            raise FileNotFoundError(f"CSUBST tree input was not found: {path}")
    context = csubst_tree.read_treefile(context)
    context = csubst_parser_misc.annotate_tree(context)
    signatures = {}
    for node in context["tree"].traverse():
        branch_id = int(csubst_ete.get_prop(node, "numerical_label"))
        signatures[branch_id] = frozenset(str(name) for name in csubst_ete.get_leaf_names(node))
    return signatures


def _parse_gene_labels(value):
    if pandas.isna(value):
        return frozenset()
    return frozenset(label.strip() for label in str(value).split(";") if label.strip())


def get_stat_branch_clade_signatures(file_stat_branch):
    """Return stat_branch IDs keyed to descendant tips in its plotted topology."""
    frame = pandas.read_csv(file_stat_branch, sep="\t", header=0, low_memory=False)
    required = {"branch_id", "node_name"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"stat_branch is missing columns required for branch mapping: {', '.join(missing)}")
    branch_ids = pandas.to_numeric(frame["branch_id"], errors="coerce")
    if branch_ids.isna().any() or not numpy.equal(branch_ids, numpy.floor(branch_ids)).all():
        raise ValueError(f"stat_branch contains non-integer branch IDs: {file_stat_branch}")
    frame = frame.copy()
    frame["branch_id"] = branch_ids.astype(int)
    if frame["branch_id"].duplicated().any():
        raise ValueError(f"stat_branch contains duplicate branch IDs: {file_stat_branch}")

    if "gene_labels" in frame.columns:
        signatures = {
            int(row.branch_id): _parse_gene_labels(row.gene_labels)
            for row in frame.loc[:, ["branch_id", "gene_labels"]].itertuples(index=False)
        }
        if all(len(signature) > 0 for signature in signatures.values()):
            return signatures

    child_columns = [column for column in ("child1", "child2") if column in frame.columns]
    if len(child_columns) != 2:
        raise ValueError(
            f"stat_branch needs complete gene_labels or child1/child2 columns for branch mapping: {file_stat_branch}"
        )
    rows = frame.set_index("branch_id", drop=False)
    signatures = {}
    active = set()

    def descendant_tips(branch_id):
        branch_id = int(branch_id)
        if branch_id in signatures:
            return signatures[branch_id]
        if branch_id in active:
            raise ValueError(f"stat_branch topology contains a cycle at branch {branch_id}.")
        if branch_id not in rows.index:
            raise ValueError(f"stat_branch references missing child branch {branch_id}.")
        active.add(branch_id)
        row = rows.loc[branch_id]
        children = []
        for column in child_columns:
            child = pandas.to_numeric(pandas.Series([row[column]]), errors="coerce").iloc[0]
            if pandas.isna(child) or int(child) < 0:
                continue
            children.append(int(child))
        if len(children) == 0:
            node_name = str(row["node_name"]).strip()
            if node_name == "" or node_name.lower() == "nan":
                raise ValueError(f"stat_branch leaf {branch_id} has no node_name.")
            signature = frozenset([node_name])
        else:
            signature = frozenset().union(*(descendant_tips(child) for child in children))
        active.remove(branch_id)
        signatures[branch_id] = signature
        return signature

    for branch_id in frame["branch_id"].tolist():
        descendant_tips(branch_id)
    return signatures


def validate_all_csubst_stat_branch_identity(file_stat_branch, iqtree_anc_dir):
    """Require every CSUBST/stat_branch ID to describe the same rooted clade."""
    csubst_signatures = get_csubst_branch_clade_signatures(iqtree_anc_dir)
    stat_signatures = get_stat_branch_clade_signatures(file_stat_branch)
    csubst_ids = set(csubst_signatures)
    stat_ids = set(stat_signatures)
    mismatched_ids = sorted(
        branch_id
        for branch_id in csubst_ids.intersection(stat_ids)
        if csubst_signatures[branch_id] != stat_signatures[branch_id]
    )
    only_csubst = sorted(csubst_ids.difference(stat_ids))
    only_stat = sorted(stat_ids.difference(csubst_ids))
    if mismatched_ids or only_csubst or only_stat:
        details = []
        for branch_id in mismatched_ids[:5]:
            csubst_preview = ",".join(sorted(csubst_signatures[branch_id])[:3])
            stat_preview = ",".join(sorted(stat_signatures[branch_id])[:3])
            details.append(f"{branch_id}: CSUBST=[{csubst_preview}] stat_branch=[{stat_preview}]")
        if only_csubst:
            details.append(f"IDs only in CSUBST={only_csubst[:5]}")
        if only_stat:
            details.append(f"IDs only in stat_branch={only_stat[:5]}")
        raise RuntimeError(
            "CSUBST and stat_branch branch IDs were generated from inconsistent rooted trees. "
            "Refusing to remap IDs because that would hide stale or mixed-provenance artifacts. "
            "Rebuild iqtree_anc and its dependent CSUBST outputs from the same rooted tree used "
            "to generate stat_branch. First differences: " + "; ".join(details)
        )
    print(
        f"Validated identical CSUBST/stat_branch branch IDs for {len(csubst_ids)} rooted clades.",
        flush=True,
    )


def validate_csubst_stat_branch_identity(branch_id_str, file_stat_branch, iqtree_anc_dir):
    """Require requested and complete CSUBST/stat_branch IDs to match rooted clades."""
    csubst_branch_ids = _parse_branch_id_list(branch_id_str)
    csubst_signatures = get_csubst_branch_clade_signatures(iqtree_anc_dir)
    for csubst_branch_id in csubst_branch_ids:
        if csubst_branch_id not in csubst_signatures:
            raise ValueError(f"CSUBST branch ID was not found in its tree: {csubst_branch_id}")
    validate_all_csubst_stat_branch_identity(file_stat_branch, iqtree_anc_dir)


def materialize_csubst_site_inputs(dir_og, og, destination_root):
    store = GeneFamilyOutputStore(dir_og, family_filter=og)
    candidates_by_subdir = get_csubst_site_input_candidates(og)
    materialized = []
    for subdir, candidate_names in candidates_by_subdir.items():
        for name in candidate_names:
            if store.artifact(subdir, name) is None:
                continue
            materialized.append(
                store.materialize(
                    subdir,
                    name,
                    destination_root=destination_root,
                )
            )
    return materialized


def materialize_tree_plot_alignment(alignment_path, plain_path):
    if alignment_path.endswith(".gz"):
        if (not os.path.exists(plain_path)) or (os.path.getmtime(plain_path) < os.path.getmtime(alignment_path)):
            with gzip.open(alignment_path, "rt") as fin, open(plain_path, "w") as fout:
                shutil.copyfileobj(fin, fout)
        return plain_path
    return alignment_path


def open_text_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta_records(path):
    records = []
    name = None
    seq_lines = []
    with open_text_maybe_gzip(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "":
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_lines)))
                name = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if name is not None:
        records.append((name, "".join(seq_lines)))
    return records


def write_fasta_records(records, path):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def translate_codon_to_aa(codon):
    codon = str(codon).upper().replace("U", "T")
    if len(codon) != 3 or any(base not in {"A", "C", "G", "T"} for base in codon):
        return "-"
    aa = STANDARD_GENETIC_CODE.get(codon, "-")
    if aa == "*":
        return "-"
    return aa


def load_nonsyn_recoding_symbol_map(recoding_table_path):
    recoding_df = pandas.read_csv(recoding_table_path, sep="\t", header=0, index_col=None, dtype=str)
    required_columns = {"amino_acid", "state_id", "state_label"}
    missing_columns = required_columns.difference(recoding_df.columns)
    if len(missing_columns) > 0:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Nonsynonymous recoding table is missing required column(s): {missing_text}")
    recoding_df = recoding_df.dropna(subset=["amino_acid", "state_id", "state_label"]).copy()
    recoding_df["state_id_sort"] = pandas.to_numeric(recoding_df["state_id"], errors="coerce")
    states = recoding_df.loc[:, ["state_id_sort", "state_id", "state_label"]].drop_duplicates()
    if states["state_id_sort"].notna().all():
        states = states.sort_values(["state_id_sort", "state_label"], kind="stable")
    else:
        states = states.drop_duplicates(subset=["state_label"], keep="first")
    state_labels = states["state_label"].astype(str).tolist()
    if len(state_labels) > len(NSY_ALIGNMENT_SYMBOLS):
        raise ValueError(f"Too many recoded nonsynonymous states to plot: {len(state_labels)}")
    state_to_symbol = {
        state_label: NSY_ALIGNMENT_SYMBOLS[state_index] for state_index, state_label in enumerate(state_labels)
    }
    aa_to_symbol = {}
    for _, row in recoding_df.iterrows():
        aa = str(row["amino_acid"]).strip().upper()
        state_label = str(row["state_label"])
        if len(aa) == 1:
            aa_to_symbol[aa] = state_to_symbol[state_label]
    return aa_to_symbol


def recode_codon_alignment_records(records, aa_to_symbol):
    recoded_records = []
    for name, seq in records:
        seq = str(seq).upper().replace("U", "T")
        if (len(seq) % 3) != 0:
            raise ValueError(f"CDS alignment sequence length is not a multiple of 3: {name}")
        recoded_seq = []
        for start in range(0, len(seq), 3):
            aa = translate_codon_to_aa(seq[start : start + 3])
            recoded_seq.append(aa_to_symbol.get(aa, "-"))
        recoded_records.append((name, "".join(recoded_seq)))
    return recoded_records


def write_recoded_site_alignment(codon_alignment_path, recoding_table_path, output_path):
    aa_to_symbol = load_nonsyn_recoding_symbol_map(recoding_table_path)
    records = read_fasta_records(codon_alignment_path)
    if len(records) == 0:
        raise ValueError(f"CDS alignment FASTA contains no records: {codon_alignment_path}")
    recoded_records = recode_codon_alignment_records(records=records, aa_to_symbol=aa_to_symbol)
    write_fasta_records(records=recoded_records, path=output_path)
    return output_path


def resolve_nonsyn_recoding_table(dir_out_og, site_dir=None):
    candidates = [os.path.join(dir_out_og, "csubst_nonsyn_recoding.tsv")]
    if site_dir is not None:
        candidates.append(os.path.join(site_dir, "csubst_nonsyn_recoding.tsv"))
    return resolve_existing_path(candidates)


def prepare_recoded_site_alignment(dir_out_og, og, site_dir, codon_alignment_path, csubst_nonsyn_recode):
    recode = normalize_csubst_nonsyn_recode(csubst_nonsyn_recode)
    if recode == "no":
        return None
    if recode == "3di20":
        print("Recoded site panel skipped: 3di20 states cannot be derived from the 20-aa CDS alignment.", flush=True)
        return None
    recoding_table_path = resolve_nonsyn_recoding_table(dir_out_og=dir_out_og, site_dir=site_dir)
    if recoding_table_path is None:
        print(
            f"Recoded site panel skipped: csubst_nonsyn_recoding.tsv was not found for --nonsyn_recode {recode}.",
            flush=True,
        )
        return None
    output_path = os.path.join(dir_out_og, f"{og}_csubst_sites.{recode}.plot.fasta")
    write_recoded_site_alignment(
        codon_alignment_path=codon_alignment_path,
        recoding_table_path=recoding_table_path,
        output_path=output_path,
    )
    print(f"Recoded site alignment written: {output_path}", flush=True)
    return output_path


def get_alignment_for_tree_plot(dir_og, og, dir_out_og):
    alignment_candidates = [
        os.path.join(dir_og, "clipkit", og + "_cds.clipkit.fa.gz"),
        os.path.join(dir_og, "clipkit", og + "_cds.clipkit.fasta"),
        os.path.join(dir_og, "clipkit", og + "_cds.clipkit.fa"),
    ]
    alignment_path = resolve_existing_path(alignment_candidates)
    if alignment_path is None:
        raise FileNotFoundError(f"Alignment file was not found for tree plotting. Checked: {alignment_candidates}")
    return materialize_tree_plot_alignment(
        alignment_path,
        os.path.join(dir_out_og, og + "_cds.clipkit.plot.fasta"),
    )


def get_untrimmed_alignment_for_tree_plot(dir_og, og, dir_out_og):
    alignment_candidates = [
        os.path.join(dir_og, "mafft", og + "_cds.aln.fa.gz"),
        os.path.join(dir_og, "mafft", og + "_cds.aln.fasta"),
        os.path.join(dir_og, "mafft", og + "_cds.aln.fa"),
        os.path.join(dir_og, "orthogroup_extraction_fasta", og + "_orthogroup_extraction.fa.gz"),
        os.path.join(dir_og, "orthogroup_extraction_fasta", og + "_orthogroup_extraction.fasta"),
        os.path.join(dir_og, "orthogroup_extraction_fasta", og + "_orthogroup_extraction.fa"),
        os.path.join(dir_og, "cds_fasta", og + "_cds.fa.gz"),
        os.path.join(dir_og, "cds_fasta", og + "_cds.fasta"),
        os.path.join(dir_og, "cds_fasta", og + "_cds.fa"),
    ]
    alignment_path = resolve_existing_path(alignment_candidates)
    if alignment_path is None:
        print(
            f"Untrimmed alignment file was not found for tree plotting. "
            f"Domain coloring in the alignment panel will be disabled. Checked: {alignment_candidates}",
            flush=True,
        )
        return None
    return materialize_tree_plot_alignment(
        alignment_path,
        os.path.join(dir_out_og, og + "_cds.untrimmed.plot.fasta"),
    )


def build_alignment_panel_spec(trimmed_alignment, untrimmed_alignment=None):
    if untrimmed_alignment is not None:
        return "alignment," + trimmed_alignment + "," + untrimmed_alignment
    return "alignment," + trimmed_alignment


def build_alignment_panel_arg(trimmed_alignment, untrimmed_alignment=None):
    return "--panel11=" + build_alignment_panel_spec(trimmed_alignment, untrimmed_alignment)


def sanitize_panel_name(value):
    sanitized = re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")
    if sanitized == "":
        return "state"
    return sanitized


def build_site_state_panel_spec(recode, convergent_site_str, recoded_site_alignment):
    recode = normalize_csubst_nonsyn_recode(recode)
    qname = "site_state_" + sanitize_panel_name(recode)
    return f"site_state,{qname},{convergent_site_str},{recoded_site_alignment},Recoded state ({recode})"


def build_tree_plot_panel_args(
    file_og_rpsblast,
    file_csubst_input_fasta,
    convergent_site_str,
    file_og_alignment,
    file_og_untrimmed_alignment=None,
    recoded_site_alignment=None,
    csubst_nonsyn_recode="no",
):
    panel_args = [
        "--panel1=tree,bl_rooted,support_unrooted,species,L",
        "--panel2=heatmap,no,abs,_,expression_",
        "--panel4=cluster_membership,100000",
        "--panel5=tiplabel",
        "--panel6=signal_peptide",
        "--panel7=transmembrane_domain",
        "--panel8=intron_number",
        "--panel9=domain," + file_og_rpsblast,
        f"--panel10=amino_acid_site,1,{convergent_site_str},{file_csubst_input_fasta}",
    ]
    panel_index = 11
    if recoded_site_alignment is not None:
        panel_args.append(
            f"--panel{panel_index}="
            + build_site_state_panel_spec(
                recode=csubst_nonsyn_recode,
                convergent_site_str=convergent_site_str,
                recoded_site_alignment=recoded_site_alignment,
            )
        )
        panel_index += 1
    panel_args.append(
        f"--panel{panel_index}=" + build_alignment_panel_spec(file_og_alignment, file_og_untrimmed_alignment)
    )
    panel_index += 1
    panel_args.append(f"--panel{panel_index}=fimo,2000,0.05")
    return panel_args


def run_stat_branch2tree_plot(
    og,
    branch_id_str,
    file_trait_color,
    dir_out_og,
    dir_og,
    ncpu=1,
    csubst_nonsyn_recode="no",
    convergent_sites=None,
    file_tree_plot_out=None,
):
    dir_myscript = os.path.realpath(os.path.dirname(__file__))
    file_stat_branch = get_stat_branch_path(dir_og=dir_og, og=og)
    file_og_rpsblast = get_rpsblast_path(dir_og=dir_og, og=og)
    file_og_alignment = get_alignment_for_tree_plot(dir_og=dir_og, og=og, dir_out_og=dir_out_og)
    file_og_untrimmed_alignment = get_untrimmed_alignment_for_tree_plot(dir_og=dir_og, og=og, dir_out_og=dir_out_og)
    iqtree_anc_dir = get_iqtree_anc_dir(dir_out_og=dir_out_og, og=og)
    file_csubst_input_fasta = os.path.join(iqtree_anc_dir, "csubst.fasta")
    artifacts = resolve_site_artifacts(dir_out_og=dir_out_og, branch_id_str=branch_id_str)
    file_csubst_site_tsv = artifacts["site_table_tsv"]
    if file_csubst_site_tsv is None:
        raise FileNotFoundError(f"CSUBST site table was not found in {artifacts['site_dir']}.")
    explicit_convergent_sites = convergent_sites is not None
    if convergent_sites is None:
        df_csubst_site = pandas.read_csv(file_csubst_site_tsv, sep="\t", header=0, index_col=None)
        convergent_sites = df_csubst_site.loc[(df_csubst_site["OCNany2spe"] > 0.5), "codon_site_alignment"].tolist()
    else:
        convergent_sites = list(convergent_sites)
    if explicit_convergent_sites and len(convergent_sites) == 0:
        raise ValueError(f"No convergent sites were selected for tree plotting: {og}")
    convergent_site_str = ":".join([str(cs) for cs in convergent_sites])
    print(f"Convergent sites extracted from {file_csubst_site_tsv}: {convergent_site_str}", flush=True)
    recoded_site_alignment = prepare_recoded_site_alignment(
        dir_out_og=dir_out_og,
        og=og,
        site_dir=artifacts["site_dir"],
        codon_alignment_path=file_csubst_input_fasta,
        csubst_nonsyn_recode=csubst_nonsyn_recode,
    )
    if file_tree_plot_out is None:
        file_tree_plot_out = og + ".tree_plot.pdf"
    validate_csubst_stat_branch_identity(
        branch_id_str=branch_id_str,
        file_stat_branch=file_stat_branch,
        iqtree_anc_dir=iqtree_anc_dir,
    )
    # gg-cache-guard: audited - the outer csubst_site artifact contract removes this directory on rebuild.
    if os.path.exists(file_tree_plot_out):
        print(f"Tree plot skipped: outfile already exists: {file_tree_plot_out}", flush=True)
        return None
    cmd = ["Rscript", os.path.join(dir_myscript, "stat_branch2tree_plot.r")]
    cmd.append("--stat_branch=" + file_stat_branch)
    cmd.append("--max_delta_intron_present=-0.5")
    cmd.append("--width=7.2")
    cmd.append("--rel_widths=")
    cmd.extend(
        build_tree_plot_panel_args(
            file_og_rpsblast=file_og_rpsblast,
            file_csubst_input_fasta=file_csubst_input_fasta,
            convergent_site_str=convergent_site_str,
            file_og_alignment=file_og_alignment,
            file_og_untrimmed_alignment=file_og_untrimmed_alignment,
            recoded_site_alignment=recoded_site_alignment,
            csubst_nonsyn_recode=csubst_nonsyn_recode,
        )
    )
    cmd.append("--show_branch_id=yes")
    cmd.append("--event_method=species_overlap")
    cmd.append("--species_color_table=" + file_trait_color)
    cmd.append("--branch_combination=" + branch_id_str)
    cmd.append("--ncpu=" + str(max(1, int(ncpu))))
    print("COMMAND: {}".format(" ".join(cmd)), flush=True)
    subprocess.run(cmd, check=True)
    shutil.move("stat_branch2tree_plot.pdf", file_tree_plot_out)
    return None


def find_file_trait_color(trait):
    file_name = f"trait_{trait}.color.tsv"
    file_trait_color = os.path.join(os.getcwd(), file_name)
    return file_trait_color


def filter_max_per_og(cb_passed, arity, args):
    sampled_indices = []
    num_before_filtering = cb_passed.shape[0]
    if args.max_per_og <= 0:
        print(
            f"K = {arity}: Skipped branch combinations that are not top {args.max_per_og} within individual orthogroups: {num_before_filtering} -> {num_before_filtering}",
            flush=True,
        )
        return cb_passed
    for og, og_indices in cb_passed.groupby("orthogroup", sort=False).indices.items():
        og_indices = numpy.asarray(og_indices, dtype=int)
        if (og_indices.shape[0] > args.max_per_og) and (args.max_per_og > 0):
            print(
                f"K = {arity}, {og}: Number of branch combinations ({og_indices.shape[0]}) exceeded --max_per_og. Only {args.max_per_og} combinations will be analyzed."
            )
            sampled_indice_positions = numpy.floor(numpy.linspace(0, og_indices.shape[0] - 1, args.max_per_og)).astype(
                int
            )
            sampled_indices += og_indices[sampled_indice_positions].tolist()
        else:
            sampled_indices += og_indices.tolist()
    cb_passed = cb_passed.loc[sampled_indices, :].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(
        f"K = {arity}: Skipped branch combinations that are not top {args.max_per_og} within individual orthogroups: {num_before_filtering} -> {num_after_filtering}",
        flush=True,
    )
    return cb_passed


def filter_fg_stem_ratio(cb_passed, arity, trait, args, no_trait_name):
    min_fg_stem_num = int(numpy.round(args.min_fg_stem_ratio * arity, decimals=0))
    num_before_filtering = cb_passed.shape[0]
    if no_trait_name:
        col = "branch_num_fg_stem"
    else:
        col = "branch_num_fg_stem_" + trait
    cb_passed = cb_passed.loc[cb_passed[col] >= min_fg_stem_num, :].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(
        f"K = {arity}: Skipped branch combinations with less than {min_fg_stem_num} foreground stem branches (prop = {args.min_fg_stem_ratio}): {num_before_filtering} -> {num_after_filtering}",
        flush=True,
    )
    return cb_passed


def write_annotated_table(cb_passed, annot, dir_out, out_name):
    cb_passed = pandas.merge(cb_passed, annot, how="left", on="orthogroup", sort=False)
    outfile = os.path.join(dir_out, out_name + ".tsv")
    print(f"Writing the table to: {outfile}", flush=True)
    cb_passed.to_csv(outfile, sep="\t", index=False)
    return cb_passed


def filter_max_per_K_1st(cb_passed, arity, trait):
    print(
        f"Number of branch combinations after filtering exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}."
    )
    print("Selecting one branch combination per orthogroup with the highest OCN.", flush=True)
    num_before_filtering = cb_passed.shape[0]
    cb_passed = cb_passed.sort_values(by=["orthogroup", "OCNany2spe"], ascending=[True, False])
    cb_passed = cb_passed.drop_duplicates(subset=["orthogroup"], keep="first").reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(
        f"K = {arity} for {trait}: Skipped branch combinations with non-highest OCN within orthogroups: {num_before_filtering} -> {num_after_filtering}",
        flush=True,
    )
    outfile = os.path.join(dir_out, out_name + ".txt")
    print(f"Writing a placeholder file to: {outfile}", flush=True)
    with open(outfile, "w") as f:
        f.write(
            f"One branch combination per orthogroup was analyzed because the number of branch combinations after all filtering still exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}.\n"
        )
    return cb_passed


def filter_max_per_K_2nd(cb_passed, arity, trait):
    print(
        f"Number of branch combinations after one-per-orthogroup filtering still exceeded the maximum: {cb_passed.shape[0]} > {args.max_per_K}."
    )
    print(
        f"Only the top {args.max_per_K} high-OCN branch combinations will be analyzed for K = {arity} of {trait}",
        flush=True,
    )
    num_before_filtering = cb_passed.shape[0]
    cb_passed = cb_passed.sort_values(by="OCNany2spe", ascending=False).iloc[: args.max_per_K, :].reset_index(drop=True)
    num_after_filtering = cb_passed.shape[0]
    print(
        f"K = {arity} for {trait}: Skipped branch combinations that are not top {args.max_per_K} at K = {arity}: {num_before_filtering} -> {num_after_filtering}",
        flush=True,
    )
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
    paragraphs = text.split("\n\n")

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
    if workspace_dir_arg not in ("", "auto"):
        return os.path.realpath(workspace_dir_arg)

    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "input")) or os.path.isdir(os.path.join(cwd, "output")):
        return os.path.realpath(cwd)

    workspace_candidate = os.path.join(cwd, "workspace")
    if os.path.isdir(os.path.join(workspace_candidate, "input")) or os.path.isdir(
        os.path.join(workspace_candidate, "output")
    ):
        return os.path.realpath(workspace_candidate)

    return os.path.realpath(workspace_candidate)


def resolve_path_arg(path_arg, workspace_dir, *relative_parts):
    if path_arg not in ("", "auto"):
        return os.path.realpath(path_arg)
    return os.path.realpath(os.path.join(workspace_dir, *relative_parts))


def raise_on_processing_failures(failures):
    if not failures:
        return
    details = "; ".join(f"{og}: {error}" for og, error in failures)
    raise RuntimeError(f"CSUBST site processing failed for {len(failures)} job(s): {details}")


if __name__ == "__main__":
    print(f"{datetime.datetime.now()}: {__file__} started.", flush=True)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-dir", metavar="PATH", default="auto", type=str, help="Path used by --workspace-dir."
    )
    parser.add_argument(
        "--dir_orthogroup", metavar="PATH", default="auto", type=str, help="Path used by --dir_orthogroup."
    )
    parser.add_argument(
        "--dir_orthofinder", metavar="PATH", default="auto", type=str, help="Path used by --dir_orthofinder."
    )
    parser.add_argument("--dir_out", metavar="PATH", default="auto", type=str, help="Path used by --dir_out.")
    parser.add_argument("--file_trait", metavar="PATH", default="auto", type=str, help="Path used by --file_trait.")
    parser.add_argument("--ncpu", metavar="INT", default=2, type=int, help="Number of CPUs to use.")
    parser.add_argument(
        "--arity_range", metavar="STR", default="3-10", type=str, help="Hyphen-separated range of arity (K)."
    )
    parser.add_argument(
        "--trait",
        metavar="STR",
        default="all",
        type=str,
        help='Comma-separated list of traits for analysis. "all" for all traits.',
    )
    parser.add_argument(
        "--skip_lower_order",
        metavar="yes|no",
        default="yes",
        type=strtobool,
        help="Skip if the branch combination is a subset of already analyzed higher-order combinations.",
    )
    parser.add_argument(
        "--min_fg_stem_ratio",
        metavar="FLOAT",
        default=0.5,
        type=float,
        help="Minimum proportion of foreground stem branches in branch combinations. For exampmle, 3 or more branches should be foreground stems at K = 6 if set to 0.5.",
    )
    parser.add_argument(
        "--max_per_og",
        metavar="INT",
        default=0,
        type=int,
        help="The maximum number of branch combinations to analyze per orthogroup per K.",
    )
    parser.add_argument(
        "--min_OCNany2spe", metavar="FLOAT", default=5.0, type=float, help="Numeric value for --min_OCNany2spe."
    )
    parser.add_argument(
        "--min_omegaCany2spe", metavar="FLOAT", default=5.0, type=float, help="Numeric value for --min_omegaCany2spe."
    )
    parser.add_argument(
        "--min_OCNCoD", metavar="FLOAT", default=5.0, type=float, help="Numeric value for --min_OCNCoD."
    )
    parser.add_argument(
        "--max_per_K",
        metavar="INT",
        default=300,
        type=int,
        help="Maximum number of branch combinations to analyze per K per trait. "
        "If exceeded, only one branch combinations per orthogroup will be analyzed. "
        "If still exceeded, the K will be skipped.",
    )
    parser.add_argument(
        "--csubst_nonsyn_recode",
        metavar="STR",
        default="no",
        type=normalize_csubst_nonsyn_recode,
        help="CSUBST nonsynonymous-state recoding scheme used for csubst sites.",
    )
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))

    if "cwd" not in locals():
        cwd = os.getcwd()
    print("Working at: {}".format(cwd))
    min_OCNany2spe = args.min_OCNany2spe
    min_omegaCany2spe = args.min_omegaCany2spe
    min_OCNCoD = args.min_OCNCoD
    workspace_dir = resolve_workspace_dir(args.workspace_dir)
    dir_og = resolve_path_arg(args.dir_orthogroup, workspace_dir, "output", "orthogroup")
    dir_of = resolve_path_arg(args.dir_orthofinder, workspace_dir, "output", "orthofinder")
    args.dir_out = resolve_path_arg(args.dir_out, workspace_dir, "output", "csubst_site")
    args.file_trait = resolve_path_arg(args.file_trait, workspace_dir, "input", "species_trait", "species_trait.tsv")
    db_file = os.path.join(dir_og, "gg_orthogroup.db")
    arity_min = int(args.arity_range.split("-")[0])
    arity_max = int(args.arity_range.split("-")[1])
    arity_range = numpy.arange(arity_min, arity_max + 1)[::-1]
    print("Arity range: {}".format(" ".join([str(a) for a in arity_range])), flush=True)
    require_sqlalchemy()
    conn = sqlalchemy.create_engine("sqlite:///" + db_file)
    inspector = sqlalchemy.inspect(conn)
    table_names = inspector.get_table_names()
    print(f"Tables in the database: {' '.join(table_names)}", flush=True)
    column_names = dict()
    for table in table_names:
        print("")
        query = sqlalchemy.text("PRAGMA TABLE_INFO(" + table + ")")
        columns = pandas.read_sql_query(
            sql=query, con=conn.connect(), index_col=None, coerce_float=True, chunksize=None
        )
        column_names[table] = columns["name"].tolist()
        print(f"Column names of {table}: {' '.join(column_names[table])}")
    df_trait = pandas.read_csv(args.file_trait, sep="\t", header=0, index_col=None)
    trait_names = [t.strip() for t in args.trait.split(",") if t.strip()]
    no_trait_name = False
    if (len(trait_names) == 1) and (trait_names[0] == "all"):
        detected_trait_names = []
        for arity in arity_range:
            if "cb" + str(arity) in column_names.keys():
                colnames = column_names["cb" + str(arity)]
                detected_trait_names = [
                    colname.replace("is_fg_", "") for colname in colnames if colname.startswith("is_fg_")
                ]
                if len(detected_trait_names) > 0:
                    break
        if len(detected_trait_names) > 0:
            trait_names = detected_trait_names
            no_trait_name = False
        else:
            # Fallback when cb* tables are not generated yet.
            trait_names = [c for c in df_trait.columns.tolist() if c != "species"]
            no_trait_name = True
            print(
                f"No trait name found in cb tables. Using all traits from {args.file_trait}: {','.join(trait_names)}",
                flush=True,
            )
    else:
        missing_trait_cols = [t for t in trait_names if t not in df_trait.columns]
        if len(missing_trait_cols) > 0:
            raise ValueError(f"Trait columns not found in {args.file_trait}: {missing_trait_cols}")
        no_trait_name = False
    if len(trait_names) == 0:
        raise ValueError(f"No analyzable trait columns were detected in {args.file_trait}.")
    print(f"Specified or detected traits: {','.join(trait_names)}", flush=True)
    generate_trait_colors(df_trait=df_trait, trait_names=trait_names)
    annot_besthits = load_annotation_besthits(dir_of)
    already_analyzed_in_greater_K = dict()
    processing_failures = []
    for trait in trait_names:
        already_analyzed_in_greater_K[trait] = dict()  # Initialize
    for arity in arity_range:
        print("")
        if "cb" + str(arity) not in table_names:
            print(f"cb{arity} is not in the database. Skipping.", flush=True)
            continue
        cb_table = "cb" + str(arity)
        cb_required_cols = get_cb_required_columns(column_names[cb_table], trait_names)
        select_cols_sql = ", ".join('"{}"'.format(col) for col in cb_required_cols)
        sql_cb = sqlalchemy.text(f"SELECT {select_cols_sql} from {cb_table}")
        cb = pandas.read_sql_query(sql=sql_cb, con=conn.connect(), index_col=None, coerce_float=True)
        sum_ocn = cb.loc[:, "OCNany2spe"].sum()
        sum_ecn = cb.loc[:, "ECNany2spe"].sum()
        sum_ocs = cb.loc[:, "OCSany2spe"].sum()
        sum_ecs = cb.loc[:, "ECSany2spe"].sum()
        print(f"Total observed/expected nonsynonymous convergence at K = {arity}: {int(sum_ocn):,} / {int(sum_ecn):,}")
        print(f"Total observed/expected synonymous convergence at K = {arity}: {int(sum_ocs):,} / {int(sum_ecs):,}")
        for trait in trait_names:
            flag_zip = True
            out_name = (
                f"csubst_site_{trait}_K{arity}_minOCN{min_OCNany2spe}"
                f"_minomegaC{min_omegaCany2spe}_minOCNCoD{min_OCNCoD}"
                f"{csubst_nonsyn_recode_output_suffix(args.csubst_nonsyn_recode)}"
            )
            print(f"{datetime.datetime.now()}: Started processing: {out_name}", flush=True)
            file_trait_color = find_file_trait_color(trait)
            dir_out = os.path.realpath(os.path.join(args.dir_out, out_name))
            out_zip = os.path.join(cwd, dir_out + ".zip")
            quarantined_zip = None
            if os.path.exists(out_zip) and not csubst_site_archive_is_complete(out_zip):
                quarantined_zip = quarantine_unverified_csubst_site_archive(out_zip)
                print(
                    f"Rebuilding unverified CSUBST site ZIP; previous file retained at: {quarantined_zip}",
                    flush=True,
                )
            if (not os.path.exists(dir_out)) and (not os.path.exists(out_zip)):
                os.makedirs(dir_out)
            conditions = True
            conditions &= cb["OCNany2spe"] >= min_OCNany2spe
            conditions &= cb["omegaCany2spe"] >= min_omegaCany2spe
            conditions &= cb["OCNCoD"] >= min_OCNCoD
            if "is_fg" in cb.columns:
                conditions &= cb.loc[:, "is_fg"] == "Y"
            else:
                conditions &= cb.loc[:, "is_fg_" + trait] == "Y"
            print(
                f"Number of branch combinations for {trait} with OCNany2spe>={min_OCNany2spe}, "
                f"omegaCany2spe>={min_omegaCany2spe}, and OCNCoD>={min_OCNCoD} at K={arity}: {conditions.sum()} "
                f"(before removing branch combinations already analyzed in greater K)",
                flush=True,
            )
            if conditions.sum() == 0:
                outfile = os.path.join(dir_out, out_name + ".txt")
                print("No branch combinations to analyze.", flush=True)
                if csubst_site_archive_is_complete(out_zip):
                    print(f"Skipped arity = {arity}. Verified output ZIP already exists: {out_zip}", flush=True)
                elif os.path.exists(dir_out):
                    print(f"Writing a placeholder file to: {outfile}", flush=True)
                    with open(outfile, "w") as f:
                        f.write("No analyzable branch combinations.\n")
                    create_csubst_site_archive(dir_out)
                    shutil.rmtree(dir_out)
                    cleanup_quarantined_csubst_site_archives(out_zip)
                continue
            cb_passed = cb.loc[conditions, :].sort_values(by="orthogroup").reset_index(drop=True)
            # gg-cache-guard: audited - the outer csubst_site artifact contract removes this directory on rebuild.
            if csubst_site_archive_is_complete(out_zip):
                print(f"Skipped arity = {arity}. Verified output ZIP already exists: {out_zip}", flush=True)
                continue
            plot_scatter(
                cb=cb_passed,
                xcol="OCNany2spe",
                ycol="omegaCany2spe",
                outbase=os.path.join(dir_out, out_name + "_OCNany2spe-omegaC"),
                polygon_xmin=min_OCNany2spe,
                polygon_ymin=min_omegaCany2spe,
            )
            plot_scatter(
                cb=cb_passed,
                xcol="OCNany2spe",
                ycol="OCNCoD",
                outbase=os.path.join(dir_out, out_name + "_OCNany2spe-OCNCoD"),
                polygon_xmin=min_OCNany2spe,
                polygon_ymin=min_OCNCoD,
            )
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
            print("Starting analyzing individual branch combinations.", flush=True)
            branch_id_cols = cb_passed.columns[cb_passed.columns.str.startswith("branch_id_")].tolist()
            og_list = cb_passed["orthogroup"].tolist()
            branch_id_str_list = [
                ",".join(str(bid) for bid in branch_ids)
                for branch_ids in cb_passed.loc[:, branch_id_cols].itertuples(index=False, name=None)
            ]
            besthit_cols = ["besthit_0.05", "besthit_0.25", "besthit_0.5", "besthit_0.75", "besthit_0.95"]
            besthit_frame = cb_passed.reindex(columns=besthit_cols)
            observed_cols = [
                "OCNany2spe",
                "ECNany2spe",
                "OCSany2spe",
                "ECSany2spe",
                "omegaCany2spe",
                "omegaCany2any",
                "omegaCdif2spe",
                "OCNCoD",
            ]
            observed_frame = cb_passed.reindex(columns=observed_cols)
            annotation_text_list = []
            for og, branch_id_str, besthit_values_raw, observed_values_raw in zip(
                og_list,
                branch_id_str_list,
                besthit_frame.itertuples(index=False, name=None),
                observed_frame.itertuples(index=False, name=None),
                strict=True,
            ):
                besthit_values = [str(value) if pandas.notna(value) else "NA" for value in besthit_values_raw]
                observed_values = dict(zip(observed_cols, observed_values_raw, strict=True))
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
                        observed_values=observed_values,
                    )
                )
            batch_order = []
            branch_ids_by_og = {}
            annotations_by_og = {}
            for og, branch_id_str, annotation_text in zip(
                og_list,
                branch_id_str_list,
                annotation_text_list,
                strict=True,
            ):
                if og not in branch_ids_by_og:
                    batch_order.append(og)
                    branch_ids_by_og[og] = []
                    annotations_by_og[og] = []
                branch_ids_by_og[og].append(branch_id_str)
                annotations_by_og[og].append(annotation_text)
            with ProcessPoolExecutor(max_workers=args.ncpu) as executor:
                results = executor.map(
                    process_og_batch,
                    batch_order,
                    [branch_ids_by_og[og] for og in batch_order],
                    repeat(dir_out),
                    repeat(dir_og),
                    repeat(file_trait_color),
                    repeat(args.ncpu),
                    repeat(args.csubst_nonsyn_recode),
                    [annotations_by_og[og] for og in batch_order],
                )
                for batch_results in results:
                    for og, result in batch_results:
                        if isinstance(result, Exception):
                            print(f"Error in {og}: {result}")
                            flag_zip = False
                            processing_failures.append((og, result))
            cleanup_materialization_metadata(dir_out)
            os.chdir(cwd)
            if flag_zip:
                print(f"No error detected. Zipping {dir_out}", flush=True)
                create_csubst_site_archive(dir_out)
                shutil.rmtree(dir_out)
                cleanup_quarantined_csubst_site_archives(out_zip)
            else:
                print(f"Error detected. Skipping zipping {dir_out}", flush=True)
    conn.dispose()
    raise_on_processing_failures(processing_failures)
    print(f"{datetime.datetime.now()}: {__file__} completed!", flush=True)
