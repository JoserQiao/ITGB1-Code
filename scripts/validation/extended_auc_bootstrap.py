"""Compute extended bootstrap confidence intervals for response AUCs."""

from __future__ import annotations

import pandas as pd

from scripts.core.cli import parse_project_root, project_root_parser
from scripts.core.data import load_analysis_features
from scripts.validation.robustness_analysis import bootstrap_auc


def main():
    parser = project_root_parser()
    parser.add_argument("--bootstrap", type=int, default=5000)
    root, args = parse_project_root(parser)
    rows = []
    endpoints = {
        "GSE91061": "ICI objective response",
        "GSE78220": "ICI objective response",
        "GSE176307": "ICI objective response",
        "GSE126044": "ICI objective response",
        "GSE190826": "Chemoradiotherapy pCR",
    }
    for index, (cohort, endpoint) in enumerate(endpoints.items()):
        frame = load_analysis_features(root, cohort)
        result = bootstrap_auc(frame, cohort, args.bootstrap, 20261000 + index)
        for row in result:
            row["endpoint"] = endpoint
        rows.extend(result)
    out = root / "results/robustness"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(
        out / "extended_bootstrap_auc_and_differences.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
