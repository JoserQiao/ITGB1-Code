"""Acquire TCGA ITGB1 expression and GDC pathology composition variables."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
import requests
import xenaPython as xena

from scripts.core.cli import parse_project_root


XENA_HOST = "https://tcga.xenahubs.net"
GDC_CASES = "https://api.gdc.cancer.gov/cases"
PROJECTS = ["HNSC", "BLCA", "SKCM", "LUAD", "LUSC", "READ", "COAD", "PAAD", "STAD", "LIHC"]
PATHOLOGY_FIELDS = [
    "percent_tumor_nuclei", "percent_tumor_cells", "percent_stromal_cells",
    "percent_normal_cells", "percent_necrosis",
]


def fetch_gdc_project(project):
    filters = {"op": "=", "content": {"field": "project.project_id", "value": [f"TCGA-{project}"]}}
    fields = ["case_id", "submitter_id", "samples.submitter_id", "samples.sample_type"]
    fields.extend(f"samples.portions.slides.{field}" for field in PATHOLOGY_FIELDS)
    response = requests.get(GDC_CASES, params={
        "filters": json.dumps(filters), "fields": ",".join(fields),
        "expand": "samples.portions.slides", "size": 2000,
    }, timeout=300)
    response.raise_for_status()
    return response.json()


def flatten_pathology(payload, cancer):
    rows = []
    for case in payload["data"]["hits"]:
        for sample in case.get("samples", []):
            if sample.get("sample_type") != "Primary Tumor":
                continue
            slides = [slide for portion in sample.get("portions", [])
                      for slide in portion.get("slides", [])]
            row = {
                "cancer": cancer, "case_barcode": case.get("submitter_id"),
                "sample_barcode_gdc": sample.get("submitter_id"),
                "sample_barcode": str(sample.get("submitter_id"))[:15],
                "n_slides": len(slides),
            }
            for field in PATHOLOGY_FIELDS:
                values = [slide.get(field) for slide in slides if slide.get(field) is not None]
                row[field] = float(np.median(values)) if values else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def fetch_itgb1(cancer):
    dataset = f"TCGA.{cancer}.sampleMap/HiSeqV2_PANCAN"
    samples = xena.dataset_samples(XENA_HOST, dataset, None)
    result = xena.dataset_probe_values(XENA_HOST, dataset, samples, ["ITGB1"])
    if not result or not result[1] or len(result[1][0]) != len(samples):
        raise RuntimeError(f"Incomplete Xena query for {cancer}")
    return pd.DataFrame({
        "cancer": cancer, "sample_barcode": samples,
        "ITGB1": pd.to_numeric(result[1][0], errors="coerce"),
        "xena_dataset": dataset,
    })


def main():
    warnings.filterwarnings("ignore", category=FutureWarning, module="xenaPython")
    root, _ = parse_project_root()
    raw = root / "data/raw/TCGA_pathology_purity"
    processed = root / "data/processed/TCGA_pathology_purity"
    raw.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    tables = []
    for cancer in PROJECTS:
        payload = fetch_gdc_project(cancer)
        (raw / f"TCGA-{cancer}_gdc_cases.json").write_text(
            json.dumps(payload), encoding="utf-8")
        pathology = flatten_pathology(payload, cancer)
        expression = fetch_itgb1(cancer)
        merged = expression.merge(pathology, on=["cancer", "sample_barcode"], how="inner")
        merged.to_csv(processed / f"TCGA-{cancer}_ITGB1_pathology.tsv", sep="\t", index=False)
        tables.append(merged)
    pd.concat(tables, ignore_index=True).to_csv(
        processed / "TCGA_ITGB1_pathology_combined.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
