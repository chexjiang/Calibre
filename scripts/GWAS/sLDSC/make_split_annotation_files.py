#!/usr/bin/env python3
"""Materialize per-category .annot.gz files from chunked joint annotations.

LDSC needs these files when --overlap-annot is used.  The corresponding
.l2.ldscore.gz files are produced by split_chunk_ldscores.py.
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-root", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--bfile-prefix", required=True)
    ap.add_argument("--chunk-name", default=None)
    ap.add_argument("--chrom", type=int, default=None)
    args = ap.parse_args()

    chunk_root = Path(args.chunk_root)
    out = Path(args.out_dir)
    chunks = sorted(p for p in chunk_root.iterdir()
                    if p.is_dir() and p.name.startswith("chunk_"))
    if args.chunk_name is not None:
        chunks = [chunk_root / args.chunk_name]
    for chunk in chunks:
        manifest = pd.read_csv(chunk / "annotation_manifest.tsv", sep="\t")
        ids = manifest["annotation_id"].tolist()
        chroms = [args.chrom] if args.chrom is not None else range(1, 23)
        for chrom in chroms:
            src = chunk / f"joint.{chrom}.annot.gz"
            annot = pd.read_csv(src, sep=r"\s+")
            bim = pd.read_csv(
                f"{args.bfile_prefix}{chrom}.bim", sep=r"\s+", header=None,
                names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
            base = bim[["CHR", "SNP", "BP"]]
            for ann_id in ids:
                target_dir = out / "split" / ann_id
                target_dir.mkdir(parents=True, exist_ok=True)
                target = target_dir / f"{ann_id}.{chrom}.annot.gz"
                if target.exists() and target.stat().st_size > 0:
                    continue
                one = base.copy()
                one[ann_id] = annot[ann_id].astype("int8")
                one.to_csv(target, sep="\t", index=False, compression="gzip")
            print(f"annot {chunk.name} chr{chrom}", flush=True)


if __name__ == "__main__":
    main()
