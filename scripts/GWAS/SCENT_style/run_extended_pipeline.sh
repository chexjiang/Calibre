#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
NBOOT="${2:-1000}"
SEED="${3:-20261003}"

Rscript "$ROOT/scripts/GWAS/SCENT_style/prepare_extended_references.R" "$ROOT"
Rscript "$ROOT/scripts/GWAS/SCENT_style/build_extended_annotation_manifest.R" "$ROOT"
for source in FinnGen RA T1D IBD; do
  Rscript "$ROOT/scripts/GWAS/SCENT_style/run_extended_scent.R" "$ROOT" "$source" "$NBOOT" "$SEED"
  Rscript "$ROOT/scripts/GWAS/SCENT_style/plot_extended_scent.R" "$ROOT" "$source"
done

echo "SCENT-style FinnGen/RA/T1D/IBD benchmark complete."
