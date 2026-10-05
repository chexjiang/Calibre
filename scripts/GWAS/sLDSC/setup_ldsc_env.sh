#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
ENV_DIR="$ROOT/results/validation/GWAS/sLDSC/software/ldsc-env"

if [[ ! -d "$ENV_DIR" ]]; then
  mamba create --prefix "$ENV_DIR" -c conda-forge -c bioconda -y \
    python=2.7 numpy=1.16 pandas=0.24 scipy=1.2 bitarray=0.8 nose=1.3 pybedtools=0.7
elif ! "$ENV_DIR/bin/python" -c 'import numpy, pandas, scipy, bitarray' >/dev/null 2>&1; then
  mamba install --prefix "$ENV_DIR" -c conda-forge -c bioconda -y \
    numpy=1.16 pandas=0.24 scipy=1.2 bitarray=0.8 nose=1.3 pybedtools=0.7
fi

echo "LDSC environment: $ENV_DIR"
