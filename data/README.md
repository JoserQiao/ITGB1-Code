# Data policy

Raw and processed datasets are not committed. `config/datasets.tsv` records each public source URL and expected local path. Running `scripts/run_all.ps1` downloads and preprocesses the inputs before analysis.

The committed `results/` directory contains small machine-readable estimates used in the manuscript; it does not contain patient-level expression matrices.
