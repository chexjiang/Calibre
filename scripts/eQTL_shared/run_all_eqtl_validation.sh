#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

Rscript scripts/eQTL_shared/prepare_reference.R
Rscript scripts/eQTL1_SCENT/run_eqtl1_scent.R
Rscript scripts/eQTL2_pgBoost/run_eqtl2_ranked.R
Rscript scripts/eQTL_shared/plot_eqtl_validation.R

