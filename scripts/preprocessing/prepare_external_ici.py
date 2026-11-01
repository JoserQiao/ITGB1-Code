"""Prepare independent immune-checkpoint inhibitor cohorts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from scripts.core.cli import parse_project_root
from scripts.core.data import align_expression_metadata
from scripts.preprocessing.prepare_geo import collapse_duplicate_genes, ensembl_symbol_map, parse_miniml


def save_prepared_cohort(
    root: Path,
    cohort: str,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    expression_filename: str,
) -> None:
    expression, metadata = align_expression_metadata(expression, metadata, "patient_id")
    output = root / "data" / "processed" / cohort
    output.mkdir(parents=True, exist_ok=True)
    expression.to_csv(output / expression_filename, sep="\t", compression="gzip")
    metadata.reset_index().to_csv(output / "metadata_v2.tsv", sep="\t", index=False)


def prepare_gse176307(root: Path) -> None:
    raw = root / "data/raw/GSE176307"
    expr = pd.read_csv(raw / "GSE176307_salmon_tpm_gene.matrix.tsv.gz", sep="\t", index_col=0)
    key = pd.read_csv(raw / "GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv.gz")
    mapping = dict(zip(key["Omniseq_RS_ID (RNAseq)"].astype(str), key["Sample ID"].astype(str)))
    keep = [c for c in expr.columns if c in mapping]
    expr = expr[keep]
    expr.columns = [mapping[c] for c in keep]
    expr = expr.T.groupby(level=0).mean().T
    expr = np.log2(collapse_duplicate_genes(expr) + 1)
    meta = parse_miniml(raw / "GSE176307_family.xml.tgz")
    meta["patient_id"] = meta.title.str.extract(r"(BACI\d+)", expand=False)
    meta = meta.drop_duplicates("patient_id")
    meta = meta[meta.io_response.str.fullmatch(r"CR|PR|SD|PD", case=False, na=False)].copy()
    meta["response_binary"] = meta.io_response.str.fullmatch(r"CR|PR", case=False).astype(int)
    save_prepared_cohort(root, "GSE176307", expr, meta, "expression_log2_tpm.tsv.gz")


def prepare_gse126044(root: Path) -> None:
    raw = root / "data/raw/GSE126044"
    counts = pd.read_csv(raw / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    counts = collapse_duplicate_genes(counts)
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    expr = np.log2(cpm + 1)
    meta = parse_miniml(raw / "GSE126044_family.xml.tgz")
    meta = meta[meta.treatment_status.str.contains("pre", case=False, na=False)].copy()
    meta["patient_id"] = meta.title.str.extract(r"(Dis_\d+)", expand=False)
    meta["response_binary"] = meta.patient_response.str.fullmatch("responder", case=False).astype(int)
    save_prepared_cohort(root, "GSE126044", expr, meta, "expression_log2_cpm.tsv.gz")


def prepare_gse135222(root: Path) -> None:
    raw = root / "data/raw/GSE135222"
    expr = pd.read_csv(raw / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0)
    mapping = ensembl_symbol_map(root)
    expr.index = [mapping.get(str(x).split(".")[0], str(x).split(".")[0]) for x in expr.index]
    expr = np.log2(collapse_duplicate_genes(expr) + 1)
    meta = parse_miniml(raw / "GSE135222_family.xml.tgz")
    meta["patient_id"] = "NSCLC" + meta.title.str.extract(r"NSCLC\s*(\d+)", expand=False)
    meta["pfs_time"] = pd.to_numeric(meta.pfs_time, errors="coerce")
    meta["pfs_event"] = pd.to_numeric(meta.progression_free_survival_pfs, errors="coerce")
    meta = meta.dropna(subset=["pfs_time", "pfs_event", "patient_id"])
    save_prepared_cohort(root, "GSE135222", expr, meta, "expression_log2.tsv.gz")


def main() -> None:
    root, _ = parse_project_root()
    prepare_gse176307(root)
    prepare_gse126044(root)
    prepare_gse135222(root)


if __name__ == "__main__":
    main()
