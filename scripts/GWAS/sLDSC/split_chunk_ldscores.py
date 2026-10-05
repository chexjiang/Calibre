#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-root", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    chunk_root = Path(args.chunk_root)
    out = Path(args.out_dir)
    chunks = sorted(p for p in chunk_root.iterdir() if p.is_dir() and p.name.startswith("chunk_"))
    for chunk in chunks:
        manifest = pd.read_csv(chunk / "annotation_manifest.tsv", sep="\t")
        ids = manifest["annotation_id"].tolist()
        for chrom in range(1, 23):
            ld_path = chunk / f"joint.{chrom}.l2.ldscore.gz"
            ld = pd.read_csv(ld_path, sep=r"\s+")
            m = (chunk / f"joint.{chrom}.l2.M").read_text().strip().split()
            m5 = (chunk / f"joint.{chrom}.l2.M_5_50").read_text().strip().split()
            score_cols = [c for c in ld.columns if c not in ["CHR", "SNP", "BP"]]
            if len(score_cols) != len(ids):
                raise RuntimeError(f"{chunk.name} chr{chrom}: {len(score_cols)} columns vs {len(ids)} ids")
            for j, ann_id in enumerate(ids):
                ann_dir = out / "split" / ann_id
                ann_dir.mkdir(parents=True, exist_ok=True)
                one = ld[["CHR", "SNP", "BP", score_cols[j]]].copy()
                one.columns = ["CHR", "SNP", "BP", ann_id]
                one.to_csv(ann_dir / f"{ann_id}.{chrom}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
                (ann_dir / f"{ann_id}.{chrom}.l2.M").write_text(str(m[j]) + "\n")
                (ann_dir / f"{ann_id}.{chrom}.l2.M_5_50").write_text(str(m5[j]) + "\n")
            print(f"split {chunk.name} chr{chrom}")


if __name__ == "__main__":
    main()
