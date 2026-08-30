#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
Rscript scripts/plot_figures4_5.R
Rscript scripts/plot_supplementary.R
python3 scripts/plot_figure6.py
