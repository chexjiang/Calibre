#!/usr/bin/env python3
"""Export the count layer of the processed Open Problems BMMC h5ad.

The h5ad already contains the final RNA and ATAC feature definitions.  This
script deliberately does not perform QC, normalization, or peak calling.
"""

from pathlib import Path
import argparse
import gzip

import anndata as ad
import pandas as pd
from scipy.io import mmwrite


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if "counts" not in adata.layers:
        raise RuntimeError("Expected layers['counts'] in the processed h5ad")

    feature_types = adata.var["feature_types"].astype(str)
    masks = {"GEX": feature_types.eq("GEX"), "ATAC": feature_types.eq("ATAC")}
    counts = adata.layers["counts"]

    # R/Matrix can read Matrix Market files directly from gzip connections.
    for assay, mask in masks.items():
        matrix = counts[:, mask.to_numpy()].T.tocsc()
        with gzip.open(args.outdir / f"{assay}_counts.mtx.gz", "wb") as handle:
            mmwrite(handle, matrix)

        var = adata.var.loc[mask, ["gene_id"]].copy()
        var.index.name = "feature"
        var.to_csv(args.outdir / f"{assay}_features.tsv", sep="\t")

    obs = adata.obs.copy()
    obs.index.name = "barcode"
    obs.to_csv(args.outdir / "metadata.tsv", sep="\t")

    print(f"cells={adata.n_obs} RNA={masks['GEX'].sum()} ATAC={masks['ATAC'].sum()}")


if __name__ == "__main__":
    main()
