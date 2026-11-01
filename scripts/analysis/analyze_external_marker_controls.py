"""Compare ITGB1 with integrin and stromal-marker negative controls.

Each marker is residualised against the same full transcriptomic-composition
model used for ITGB1.  The marker itself is removed from every signature before
score construction to prevent self-inclusion.  Residual-marker associations
with prespecified Hallmark pathways are pooled across the three external
solid-tumour cohorts with modified Hartung-Knapp inference.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

from scripts.analysis.analyze_gdsc_lineage_adjusted import random_effects_fisher_z
from scripts.core.cli import parse_project_root, project_root_parser
from scripts.core.common import COMPOSITION_MODELS, bh_fdr, read_signatures, signature_score


COHORTS = ["GSE172356", "GSE113255", "GSE76427"]
MARKERS = ["ITGB1", "ITGAV", "ITGA5", "ITGB5", "SPARC", "FAP", "COL1A1"]
PATHWAYS = {
    "EMT": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "TGF_beta": "HALLMARK_TGF_BETA_SIGNALING",
    "Hypoxia": "HALLMARK_HYPOXIA",
}


def load_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            sets[fields[0]] = fields[2:]
    return sets


def residualise_marker(expression, signatures, marker):
    adjusted_signatures = {
        name: [gene for gene in genes if gene.upper() != marker.upper()]
        for name, genes in signatures.items()
    }
    frame = pd.DataFrame({marker: expression.loc[marker]})
    for name, genes in adjusted_signatures.items():
        frame[name] = signature_score(expression, genes)
    predictors = COMPOSITION_MODELS["FullComposition"]
    complete = frame[[marker, *predictors]].dropna()
    fit = sm.OLS(complete[marker], sm.add_constant(complete[predictors])).fit()
    frame[f"{marker}_residual"] = np.nan
    frame.loc[complete.index, f"{marker}_residual"] = fit.resid
    return frame, float(fit.rsquared)


def main() -> None:
    parser = project_root_parser(__doc__)
    parser.add_argument("--gene-set-root", default=None)
    root, args = parse_project_root(parser)
    gene_set_root = Path(args.gene_set_root).resolve() if args.gene_set_root else root
    signatures = read_signatures(root / "config/signatures.tsv")
    hallmark = load_gmt(gene_set_root / "data/raw/gene_sets/h.all.v2025.1.Hs.symbols.gmt")
    pathway_genes = {
        label: [gene for gene in hallmark[name] if gene.upper() not in {m.upper() for m in MARKERS}]
        for label, name in PATHWAYS.items()
    }
    rows = []
    for cohort in COHORTS:
        expression = pd.read_csv(
            root / f"data/processed/{cohort}/expression_log2.tsv.gz", sep="\t", index_col=0
        )
        for marker in MARKERS:
            if marker not in expression.index:
                continue
            frame, r_squared = residualise_marker(expression, signatures, marker)
            for pathway, genes in pathway_genes.items():
                score = signature_score(expression, genes)
                data = pd.concat([frame[f"{marker}_residual"], score.rename("pathway")], axis=1).dropna()
                rho, p_value = spearmanr(data.iloc[:, 0], data.pathway)
                rows.append({
                    "cohort": cohort, "marker": marker, "pathway": pathway,
                    "n": len(data), "composition_r_squared": r_squared,
                    "spearman_rho": rho, "p_value": p_value,
                })
    table = pd.DataFrame(rows)
    table["fdr"] = bh_fdr(table.p_value)
    out = root / "results/external_solid_tumor"
    table.to_csv(out / "external_marker_control_associations.tsv", sep="\t", index=False)

    meta_rows = []
    for (marker, pathway), group in table.groupby(["marker", "pathway"]):
        meta_rows.append({"marker": marker, "pathway": pathway, **random_effects_fisher_z(group)})
    meta = pd.DataFrame(meta_rows)
    meta["meta_fdr"] = bh_fdr(meta.meta_p_value)
    meta.to_csv(out / "external_marker_control_meta.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
