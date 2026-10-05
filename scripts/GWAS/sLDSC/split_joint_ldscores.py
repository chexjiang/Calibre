#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joint-prefix", required=True, help="Prefix ending before chromosome number")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    ids = manifest["annotation_id"].tolist()
    root = Path(args.out_dir)

    for chrom in range(1, 23):
        ld_path = Path(f"{args.joint_prefix}{chrom}.l2.ldscore.gz")
        ld = pd.read_csv(ld_path, sep=r"\s+")
        m = Path(f"{args.joint_prefix}{chrom}.l2.M").read_text().strip().split()
        m5 = Path(f"{args.joint_prefix}{chrom}.l2.M_5_50").read_text().strip().split()
        score_cols = [c for c in ld.columns if c not in ["CHR", "SNP", "BP"]]
        if len(score_cols) != len(ids):
            raise RuntimeError(f"chr{chrom}: {len(score_cols)} LD-score columns but {len(ids)} annotations")
        for j, ann_id in enumerate(ids):
            ann_dir = root / "split" / ann_id
            out_ld = ann_dir / f"{ann_id}.{chrom}.l2.ldscore.gz"
            one = ld[["CHR", "SNP", "BP", score_cols[j]]].copy()
            one.columns = ["CHR", "SNP", "BP", ann_id]
            one.to_csv(out_ld, sep="\t", index=False, compression="gzip")
            (ann_dir / f"{ann_id}.{chrom}.l2.M").write_text(str(m[j]) + "\n")
            (ann_dir / f"{ann_id}.{chrom}.l2.M_5_50").write_text(str(m5[j]) + "\n")
        print(f"split chr{chrom}")


if __name__ == "__main__":
    main()
