"""Quantify the cellular origin of ITGB1 in single-cell data."""

from __future__ import annotations

import gzip
import re
import tarfile

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix
from scipy.stats import wilcoxon

from scripts.core.cli import parse_project_root


MARKERS = {
    "Fibroblast_CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "FAP", "PDGFRA"],
    "Epithelial_malignant_candidate": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "SFN"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "CD247", "NKG7", "GNLY"],
    "Myeloid": ["LYZ", "FCER1G", "TYROBP", "LST1", "AIF1", "CTSS"],
    "B_Plasma": ["CD79A", "MS4A1", "CD37", "CD74", "MZB1", "JCHAIN"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG", "RAMP2"],
}


def read_member_gzip(tf, member_name, sep="\t", header=None):
    with gzip.GzipFile(fileobj=tf.extractfile(member_name)) as gh:
        return pd.read_csv(gh, sep=sep, header=header)


def analyze_sample(tf, prefix, sample_title):
    features_name = next(x.name for x in tf.getmembers() if x.name.startswith(prefix) and "features.tsv.gz" in x.name)
    matrix_name = next(x.name for x in tf.getmembers() if x.name.startswith(prefix) and "matrix.mtx.gz" in x.name)
    features = read_member_gzip(tf, features_name)
    genes = features.iloc[:, 1].astype(str).str.upper().to_numpy()
    with gzip.GzipFile(fileobj=tf.extractfile(matrix_name)) as gh:
        matrix = csc_matrix(mmread(gh))
    n_genes = np.asarray((matrix > 0).sum(axis=0)).ravel()
    total = np.asarray(matrix.sum(axis=0)).ravel()
    mt = np.flatnonzero(np.char.startswith(genes.astype(str), "MT-"))
    mt_counts = np.asarray(matrix[mt].sum(axis=0)).ravel() if len(mt) else np.zeros(matrix.shape[1])
    keep = (n_genes >= 300) & (n_genes <= 7000) & (total >= 1000) & ((mt_counts / np.maximum(total, 1)) <= 0.20)
    matrix = matrix[:, keep]
    total = total[keep]
    # log1p counts per 10,000 without densifying the full matrix.
    norm = matrix.multiply(1e4 / np.maximum(total, 1)).tocsc()
    norm.data = np.log1p(norm.data)
    scores = []
    for label, marker_genes in MARKERS.items():
        idx = np.flatnonzero(np.isin(genes, marker_genes))
        score = np.asarray(norm[idx].mean(axis=0)).ravel() if len(idx) else np.zeros(norm.shape[1])
        scores.append(score)
    score_matrix = np.vstack(scores)
    best = score_matrix.argmax(axis=0)
    best_score = score_matrix.max(axis=0)
    labels = np.array(list(MARKERS))[best]
    labels[best_score < 0.15] = "Unclassified"
    itgb_idx = np.flatnonzero(genes == "ITGB1")
    itgb = np.asarray(norm[itgb_idx[0], :].todense()).ravel() if len(itgb_idx) else np.zeros(norm.shape[1])
    patient = re.search(r"Pt\d+", sample_title).group(0)
    timepoint = "Pre" if " Pre" in sample_title else "Post"
    rows = []
    for label in list(MARKERS) + ["Unclassified"]:
        sel = labels == label
        if not sel.any():
            continue
        rows.append({"patient_id": patient, "timepoint": timepoint, "cell_type": label,
                     "n_cells": int(sel.sum()), "mean_ITGB1": float(itgb[sel].mean()),
                     "median_ITGB1": float(np.median(itgb[sel])),
                     "ITGB1_detected_fraction": float((itgb[sel] > 0).mean())})
    return rows, {"sample": sample_title, "cells_raw": int(len(keep)), "cells_qc": int(keep.sum())}


def main():
    root, _ = parse_project_root()
    tar_path = root / "data/raw/GSE195832/GSE195832_RAW.tar"
    meta = pd.read_csv(root / "data/processed/GSE195832/metadata_all_v2.tsv", sep="\t")
    sc = meta[meta.title.str.startswith("ScRNA", na=False)].set_index("gsm")
    rows, qc = [], []
    with tarfile.open(tar_path) as tf:
        for gsm, record in sc.iterrows():
            member = next(x.name for x in tf.getmembers() if x.name.startswith(gsm) and "features.tsv.gz" in x.name)
            prefix = member.split("features.tsv.gz")[0]
            sample_rows, sample_qc = analyze_sample(tf, prefix, record.title)
            rows.extend(sample_rows)
            qc.append(sample_qc)
    out = root / "results/single_cell"
    out.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "patient_celltype_ITGB1.tsv", sep="\t", index=False)
    pd.DataFrame(qc).to_csv(out / "qc_summary.tsv", sep="\t", index=False)
    # Primary source ranking uses pretreatment patient-level compartment summaries.
    pre = summary[summary.timepoint.eq("Pre")]
    ranking = pre.groupby("cell_type").agg(
        patients=("patient_id", "nunique"), cells=("n_cells", "sum"),
        mean_ITGB1=("mean_ITGB1", "mean"), detected_fraction=("ITGB1_detected_fraction", "mean")
    ).reset_index().sort_values("mean_ITGB1", ascending=False)
    ranking.to_csv(out / "pretreatment_celltype_ranking.tsv", sep="\t", index=False)
    # Paired patient-level pre/post summaries; n=4, so these are descriptive
    # treatment-associated shifts rather than response biomarkers.
    paired_rows = []
    for cell_type, tab in summary.groupby("cell_type"):
        wide = tab.pivot(index="patient_id", columns="timepoint", values="mean_ITGB1").dropna()
        if len(wide) < 2 or not {"Pre", "Post"}.issubset(wide.columns):
            continue
        diff = wide.Post - wide.Pre
        p = wilcoxon(wide.Post, wide.Pre).pvalue if np.any(diff != 0) else 1.0
        paired_rows.append({"cell_type": cell_type, "patients_paired": len(wide),
                            "mean_post_minus_pre": diff.mean(),
                            "median_post_minus_pre": diff.median(), "wilcoxon_p": p})
    pd.DataFrame(paired_rows).to_csv(out / "paired_pre_post_ITGB1.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
