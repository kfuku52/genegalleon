#!/usr/bin/env python

import argparse
from pathlib import Path

import numpy
import pandas

REQUIRED_COLUMNS = (
    "qacc",
    "sacc",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qlen",
    "slen",
)

NUMERIC_COLUMNS = (
    "qstart",
    "qend",
    "sstart",
    "send",
    "qlen",
    "slen",
)

HIT_COLUMNS_TO_CONCAT = (
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "frames",
)

def _interval_union_coverage(lo, hi, length):
    if length <= 0:
        return 0.0
    if lo.size == 0:
        return 0.0
    if lo.size < 32:
        pairs = sorted(zip(lo.tolist(), hi.tolist()), key=lambda x: x[0])
        covered = 0
        cur_lo, cur_hi = pairs[0]
        for next_lo, next_hi in pairs[1:]:
            if next_lo <= cur_hi + 1:
                cur_hi = max(cur_hi, next_hi)
            else:
                covered += cur_hi - cur_lo + 1
                cur_lo = next_lo
                cur_hi = next_hi
        covered += cur_hi - cur_lo + 1
        return covered / float(length)

    order = numpy.argsort(lo, kind="mergesort")
    lo_sorted = lo[order]
    hi_sorted = hi[order]
    hi_cummax = numpy.maximum.accumulate(hi_sorted)

    breaks = lo_sorted[1:] > (hi_cummax[:-1] + 1)
    if not breaks.any():
        covered = int(hi_cummax[-1] - lo_sorted[0] + 1)
        return covered / float(length)

    break_idx = numpy.flatnonzero(breaks)
    seg_starts = numpy.concatenate(([0], break_idx + 1))
    seg_ends = numpy.concatenate((break_idx, [lo_sorted.size - 1]))
    covered = numpy.sum(hi_cummax[seg_ends] - lo_sorted[seg_starts] + 1, dtype=numpy.int64)
    return int(covered) / float(length)


def _concat_values(values):
    return ";".join(values.tolist())


def _prepare_concat_array(series):
    arr = series.to_numpy(copy=False)
    kind = arr.dtype.kind
    if kind in ("U", "S"):
        return arr
    if kind in ("i", "u", "b"):
        return arr.astype(str)
    if kind == "f":
        out = arr.astype(str)
        nan_mask = numpy.isnan(arr)
        if nan_mask.any():
            out[nan_mask] = ""
        return out

    out = numpy.empty(arr.shape[0], dtype=object)
    na_mask = pandas.isna(arr)
    out[na_mask] = ""
    if (~na_mask).any():
        out[~na_mask] = arr[~na_mask].astype(str)
    return out


def _format_float(values):
    return ";".join(f"{float(v):.15g}" for v in values)


def _prepare_blast_dataframe(df):
    missing_columns = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    out = df.copy()
    for column in NUMERIC_COLUMNS:
        out[column] = pandas.to_numeric(out[column], errors="raise")
    if "evalue" in out.columns:
        out["_evalue_num"] = pandas.to_numeric(out["evalue"], errors="coerce")
    if "bitscore" in out.columns:
        out["_bitscore_num"] = pandas.to_numeric(out["bitscore"], errors="coerce")
    return out


def _annotate_blast_coverage_serial(out):
    if out.empty:
        return out.assign(
            qhitcov=pandas.Series(dtype=float),
            shitcov=pandas.Series(dtype=float),
            qjointcov=pandas.Series(dtype=float),
            sjointcov=pandas.Series(dtype=float),
        )

    qstart_arr = out["qstart"].to_numpy(dtype=numpy.int64, copy=False)
    qend_arr = out["qend"].to_numpy(dtype=numpy.int64, copy=False)
    sstart_arr = out["sstart"].to_numpy(dtype=numpy.int64, copy=False)
    send_arr = out["send"].to_numpy(dtype=numpy.int64, copy=False)
    qlen_arr = out["qlen"].to_numpy(dtype=float, copy=False)
    slen_arr = out["slen"].to_numpy(dtype=float, copy=False)
    evalue_arr = out["_evalue_num"].to_numpy(dtype=float, copy=False) if "_evalue_num" in out.columns else None
    bitscore_arr = out["_bitscore_num"].to_numpy(dtype=float, copy=False) if "_bitscore_num" in out.columns else None
    q_lo_arr = numpy.minimum(qstart_arr, qend_arr)
    q_hi_arr = numpy.maximum(qstart_arr, qend_arr)
    s_lo_arr = numpy.minimum(sstart_arr, send_arr)
    s_hi_arr = numpy.maximum(sstart_arr, send_arr)
    qhitcov_arr = numpy.nan_to_num((q_hi_arr - q_lo_arr + 1) / qlen_arr, nan=0.0, posinf=0.0, neginf=0.0)
    shitcov_arr = numpy.nan_to_num((s_hi_arr - s_lo_arr + 1) / slen_arr, nan=0.0, posinf=0.0, neginf=0.0)

    hit_arrays = {}
    for col in HIT_COLUMNS_TO_CONCAT:
        if col in out.columns:
            hit_arrays[col] = _prepare_concat_array(out[col])

    agg_rows = []
    grouped_indices = out.groupby(["qacc", "sacc"], sort=False).indices
    for (qacc, sacc), idx in grouped_indices.items():
        idx = numpy.asarray(idx, dtype=numpy.int64)
        q_lo = q_lo_arr[idx]
        q_hi = q_hi_arr[idx]
        s_lo = s_lo_arr[idx]
        s_hi = s_hi_arr[idx]
        qlen_values = qlen_arr[idx]
        slen_values = slen_arr[idx]

        qlen = float(numpy.max(qlen_values)) if qlen_values.size else 0.0
        slen = float(numpy.max(slen_values)) if slen_values.size else 0.0

        qjointcov = _interval_union_coverage(q_lo, q_hi, qlen)
        sjointcov = _interval_union_coverage(s_lo, s_hi, slen)

        row = {
            "qacc": qacc,
            "sacc": sacc,
            "qlen": qlen,
            "slen": slen,
            "num_hits": int(idx.size),
            "qhitcov": _format_float(qhitcov_arr[idx]),
            "shitcov": _format_float(shitcov_arr[idx]),
            "qjointcov": qjointcov,
            "sjointcov": sjointcov,
        }
        if evalue_arr is not None:
            evalues = evalue_arr[idx]
            row["min_evalue"] = numpy.nan if numpy.isnan(evalues).all() else float(numpy.nanmin(evalues))
        if bitscore_arr is not None:
            bitscores = bitscore_arr[idx]
            row["max_bitscore"] = numpy.nan if numpy.isnan(bitscores).all() else float(numpy.nanmax(bitscores))

        for col, arr in hit_arrays.items():
            row[col] = _concat_values(arr[idx])

        agg_rows.append(row)

    aggregated = pandas.DataFrame(agg_rows)
    column_order = [
        "qacc",
        "sacc",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "frames",
        "qlen",
        "slen",
        "num_hits",
        "qhitcov",
        "shitcov",
        "qjointcov",
        "sjointcov",
        "min_evalue",
        "max_bitscore",
    ]
    available = [c for c in column_order if c in aggregated.columns]
    trailing = [c for c in aggregated.columns if c not in available]
    return aggregated.loc[:, available + trailing]


def _normalize_ncpu(ncpu):
    try:
        value = int(ncpu)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid ncpu value: {ncpu}")
    if value < 1:
        value = 1
    return value


def annotate_blast_coverage(df, ncpu=1):
    _normalize_ncpu(ncpu)  # Keep CLI/API validation backward-compatible.
    out = _prepare_blast_dataframe(df)
    return _annotate_blast_coverage_serial(out)


def _parse_outfmt_columns(value):
    if value is None:
        return None
    columns = [col for col in value.replace(",", " ").split() if col]
    if not columns:
        raise ValueError("--outfmt-columns was provided but no column names were parsed.")
    return columns


def read_blast_table(infile, outfmt_columns=None):
    if outfmt_columns is None:
        return pandas.read_csv(infile, sep="\t", header=0)
    return pandas.read_csv(infile, sep="\t", header=None, names=outfmt_columns)


def filter_by_frames(df, frame_filter=None):
    if frame_filter is None:
        return df
    if "frames" not in df.columns:
        raise ValueError("--frame-filter was set but 'frames' column was not found in input.")
    filtered = df.loc[df["frames"].astype(str) == frame_filter, :].copy()
    return filtered


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="infile", metavar="PATH", required=True, type=Path)
    parser.add_argument("--out", dest="outfile", metavar="PATH", required=True, type=Path)
    parser.add_argument("--ncpu", dest="ncpu", metavar="INT", default=1, type=int)
    # Backward-compatible alias.
    parser.add_argument("--threads", dest="ncpu", metavar="INT", type=int, help=argparse.SUPPRESS)
    parser.add_argument(
        "--outfmt-columns",
        dest="outfmt_columns",
        default=None,
        help="Whitespace-separated BLAST outfmt columns for headerless input.",
    )
    parser.add_argument(
        "--frame-filter",
        dest="frame_filter",
        default=None,
        help='Keep only rows with an exact match in the "frames" column (e.g., "0/1").',
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    print("Start", __file__)
    print("Setting parameters...\n")
    print("mode =", "batch", "\n")
    print("infile =", args.infile)
    print("outfile =", args.outfile, "\n")
    print("ncpu =", args.ncpu)
    if args.outfmt_columns is not None:
        print("outfmt_columns =", args.outfmt_columns)
    if args.frame_filter is not None:
        print("frame_filter =", args.frame_filter)

    outfmt_columns = _parse_outfmt_columns(args.outfmt_columns)
    df = read_blast_table(args.infile, outfmt_columns=outfmt_columns)
    df = filter_by_frames(df, frame_filter=args.frame_filter)
    out = annotate_blast_coverage(df, ncpu=args.ncpu)
    n_pairs = out[["qacc", "sacc"]].drop_duplicates().shape[0] if not out.empty else 0
    print(f"Finished checking all query:subject combinations ({n_pairs})")
    out.to_csv(args.outfile, sep="\t", index=False)
    print("End:", __file__)


if __name__ == "__main__":
    main()
