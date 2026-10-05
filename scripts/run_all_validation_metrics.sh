#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

for correction in F T; do
  for dataset in pbmc bmmc Dogma_CTRL; do
    has_results=false
    for pairs in results/real/${dataset}/*/epilink_poisson/all_pairs.tsv.gz; do
      [ -f "$pairs" ] || continue
      has_results=true
      break
    done
    if [ "$has_results" = false ]; then
      echo "[$(date '+%F %T')] skipping ${dataset}: no completed EpiLink-Poisson celltype results"
      continue
    fi
    echo "[$(date '+%F %T')] validating ${dataset} correction=${correction}"
    Rscript scripts/run_crispr_validation.R "$dataset" --correction="$correction"
  done
done
echo "[$(date '+%F %T')] all CRISPR validation metrics finished"
