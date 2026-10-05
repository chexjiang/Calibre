#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def write_binary_matrix(path, header, mat, block_rows=5000):
    with gzip.open(path, "wb", compresslevel=5) as handle:
        handle.write(("\t".join(header) + "\n").encode())
        ncol = mat.shape[1]
        separators = np.full(2 * ncol, ord("\t"), dtype=np.uint8)
        separators[-1] = ord("\n")
        for start in range(0, mat.shape[0], block_rows):
            block = mat[start:start + block_rows]
            encoded = np.tile(separators, (block.shape[0], 1))
            encoded[:, 0::2] = block + ord("0")
            handle.write(encoded.tobytes())


def write_binary_vector(path, header, values, block_rows=250000):
    with gzip.open(path, "wb", compresslevel=5) as handle:
        handle.write((header + "\n").encode())
        for start in range(0, values.shape[0], block_rows):
            block = values[start:start + block_rows]
            encoded = np.empty(block.shape[0] * 2, dtype=np.uint8)
            encoded[0::2] = block + ord("0")
            encoded[1::2] = ord("\n")
            handle.write(encoded.tobytes())


def read_bed(path, chrom):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 2), dtype=np.int64)
    x = pd.read_csv(path, sep="\t", header=None, names=["chr", "start0", "end"])
    x = x.loc[x["chr"].astype(str).str.replace("chr", "", regex=False) == str(chrom), ["start0", "end"]]
    if x.empty:
        return np.empty((0, 2), dtype=np.int64)
    x = x.sort_values(["start0", "end"]).to_numpy(dtype=np.int64)
    merged = []
    for start0, end in x:
        if not merged or start0 > merged[-1][1]:
            merged.append([start0, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return np.asarray(merged, dtype=np.int64)


def annotate(bp, intervals):
    if intervals.shape[0] == 0:
        return np.zeros(bp.shape[0], dtype=np.uint8)
    idx = np.searchsorted(intervals[:, 0], bp - 1, side="right") - 1
    ok = idx >= 0
    out = np.zeros(bp.shape[0], dtype=np.uint8)
    out[ok] = ((bp[ok] - 1) < intervals[idx[ok], 1]).astype(np.uint8)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--bfile-prefix", required=True, help="Prefix ending before chromosome number")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--methods", default="epilink_nb,epilink_poisson,scmultimap,signac")
    ap.add_argument("--min-intervals", type=int, default=200)
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    methods = set(args.methods.split(","))
    manifest = manifest.loc[
        manifest["method"].isin(methods) &
        (manifest["n_hg19_intervals"] >= args.min_intervals)
    ].copy()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out / "sldsc_annotation_manifest.tsv", sep="\t", index=False)
    (out / "selection.txt").write_text(
        "methods=" + ",".join(sorted(methods)) + "\n" +
        f"min_hg19_intervals={args.min_intervals}\n" +
        f"n_annotations={len(manifest)}\n"
    )
    ids = manifest["annotation_id"].tolist()

    for chrom in range(1, 23):
        bim_path = Path(f"{args.bfile_prefix}{chrom}.bim")
        bim = pd.read_csv(bim_path, sep=r"\s+", header=None,
                          names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        bp = bim["BP"].to_numpy(dtype=np.int64)
        columns = []
        for ann_id in ids:
            bed = Path(args.annotation_dir) / f"{ann_id}.hg19.bed"
            columns.append(annotate(bp, read_bed(bed, chrom)))
        mat = np.column_stack(columns)
        joint = out / f"joint.{chrom}.annot.gz"
        write_binary_matrix(joint, ids, mat)
        for col, ann_id in enumerate(ids):
            ann_dir = out / "split" / ann_id
            ann_dir.mkdir(parents=True, exist_ok=True)
            write_binary_vector(ann_dir / f"{ann_id}.{chrom}.annot.gz", ann_id, mat[:, col])
        print(f"chr{chrom}: {len(bp):,} variants x {len(ids)} annotations", flush=True)
    (out / "annotations.complete").write_text("complete\n")


if __name__ == "__main__":
    main()
