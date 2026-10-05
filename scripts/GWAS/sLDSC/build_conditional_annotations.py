#!/usr/bin/env python3
"""Build joint binary annotations for linked peaks and tested-peak backgrounds."""

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
    ap.add_argument("--target-manifest", required=True)
    ap.add_argument("--background-manifest", required=True)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--background-dir", required=True)
    ap.add_argument("--bfile-prefix", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chunk-size", type=int, default=25)
    args = ap.parse_args()

    targets = pd.read_csv(args.target_manifest, sep="\t")
    backgrounds = pd.read_csv(args.background_manifest, sep="\t")
    target_ids = targets.loc[targets["n_hg19_intervals"] > 0, "annotation_id"].drop_duplicates().tolist()
    background_ids = backgrounds.loc[backgrounds["n_hg19_intervals"] > 0, "background_id"].drop_duplicates().tolist()

    rows = ([{"annotation_id": x, "annotation_type": "linked_peaks",
              "bed_path": str(Path(args.annotation_dir) / f"{x}.hg19.bed")} for x in target_ids] +
            [{"annotation_id": x, "annotation_type": "tested_peak_background",
              "bed_path": str(Path(args.background_dir) / f"{x}.hg19.bed")} for x in background_ids])
    ld_manifest = pd.DataFrame(rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ld_manifest.to_csv(out / "sldsc_ld_annotation_manifest.tsv", sep="\t", index=False)

    ids = ld_manifest["annotation_id"].tolist()
    paths = ld_manifest["bed_path"].tolist()
    chunk_root = out / "chunks"
    chunk_root.mkdir(parents=True, exist_ok=True)
    chunk_manifests = []
    for start in range(0, len(ids), args.chunk_size):
        chunk_ids = ids[start:start + args.chunk_size]
        chunk_paths = paths[start:start + args.chunk_size]
        chunk_name = f"chunk_{start // args.chunk_size + 1:03d}"
        chunk_dir = chunk_root / chunk_name
        chunk_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"annotation_id": chunk_ids}).to_csv(
            chunk_dir / "annotation_manifest.tsv", sep="\t", index=False)
        chunk_manifests.append({"chunk": chunk_name, "n_annotations": len(chunk_ids)})
        for chrom in range(1, 23):
            bim = pd.read_csv(Path(f"{args.bfile_prefix}{chrom}.bim"), sep=r"\s+", header=None,
                              names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
            bp = bim["BP"].to_numpy(dtype=np.int64)
            columns = [annotate(bp, read_bed(path, chrom)) for path in chunk_paths]
            mat = np.column_stack(columns)
            write_binary_matrix(chunk_dir / f"joint.{chrom}.annot.gz", chunk_ids, mat)
        print(f"{chunk_name}: {len(chunk_ids)} annotations", flush=True)
    pd.DataFrame(chunk_manifests).to_csv(out / "chunk_manifest.tsv", sep="\t", index=False)

    (out / "annotations.complete").write_text("complete\n")


if __name__ == "__main__":
    main()
