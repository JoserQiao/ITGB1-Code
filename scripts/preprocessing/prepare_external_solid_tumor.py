"""Prepare three prespecified external solid-tumor expression cohorts."""

from __future__ import annotations

import csv
import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.core.cli import parse_project_root


def _clean_expression(expression: pd.DataFrame) -> pd.DataFrame:
    expression.index = expression.index.astype(str).str.strip().str.upper()
    expression = expression[(expression.index != "") & expression.index.notna()]
    expression = expression.apply(pd.to_numeric, errors="coerce")
    expression = expression.groupby(level=0).median()
    expression = expression.dropna(how="all")
    if "ITGB1" not in expression.index:
        raise ValueError("ITGB1 is absent after expression preprocessing")
    return expression


def read_series_metadata(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key, values = row[0][1:], [value.strip('"') for value in row[1:]]
            if key == "Sample_characteristics_ch1":
                characteristics.append(values)
            elif key in {"Sample_geo_accession", "Sample_title"}:
                fields[key] = values
    metadata = pd.DataFrame({
        "sample_id": fields["Sample_geo_accession"],
        "title": fields["Sample_title"],
    })
    for values in characteristics:
        if len(values) != len(metadata):
            continue
        keys = {value.split(":", 1)[0].strip().lower() for value in values if ":" in value}
        if len(keys) != 1:
            continue
        key = next(iter(keys)).replace(" ", "_").replace("(", "").replace(")", "")
        metadata[key] = [value.split(":", 1)[1].strip() if ":" in value else "" for value in values]
    return metadata


def read_series_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
        table = pd.read_csv(handle, sep="\t", index_col=0, na_values=["null", "NA", ""])
    table = table.loc[~table.index.astype(str).str.startswith("!series_matrix_table_end")]
    table.columns = table.columns.astype(str).str.strip('"')
    table.index = table.index.astype(str).str.strip('"')
    return table


def read_gpl_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ID\tGene title\tGene symbol"):
                header = line
                break
        text = header + handle.read()
    annotation = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, low_memory=False)
    annotation = annotation[["ID", "Gene symbol"]].dropna()
    annotation = annotation[~annotation["Gene symbol"].str.contains("///", regex=False)]
    annotation["Gene symbol"] = annotation["Gene symbol"].str.strip().str.upper()
    return annotation.drop_duplicates("ID")


def prepare_gse172356(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(raw / "GSE172356_PDA_gene_expression_matrix.txt.gz", sep="\t", index_col=0)
    # The submitted matrix combines two gene-annotation universes (41 and 21
    # samples).  Restrict to genes quantified in every sample so missing
    # annotation is not mistaken for zero expression or propagated into ranks.
    expression = expression.dropna(how="any")
    expression = _clean_expression(np.log2(expression + 1))
    metadata = pd.DataFrame({
        "sample_id": expression.columns,
        "sample_type": "primary_pancreatic_tumor",
    })
    return expression, metadata


def prepare_gse113255(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = read_series_metadata(raw / "GSE113255_series_matrix.txt.gz")
    if "tissue" not in metadata:
        raise ValueError("GSE113255 tissue annotation was not parsed")
    metadata = metadata[metadata.tissue.str.lower().eq("gastric cancer")].copy()
    wanted = set(metadata.sample_id)
    columns = []
    with tarfile.open(raw / "GSE113255_RAW.tar") as archive:
        for member in archive.getmembers():
            sample_id = member.name.split("_", 1)[0]
            if sample_id not in wanted or not member.isfile():
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            with gzip.open(extracted, "rt", encoding="utf-8", errors="replace") as handle:
                table = pd.read_csv(handle, sep="\t", index_col=0)
            series = pd.to_numeric(table.iloc[:, 0], errors="coerce").rename(sample_id)
            columns.append(series)
    expression = _clean_expression(pd.concat(columns, axis=1))
    metadata = metadata[metadata.sample_id.isin(expression.columns)].copy()
    metadata["sample_type"] = "primary_gastric_tumor"
    return expression.loc[:, metadata.sample_id], metadata


def prepare_gse76427(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = read_series_metadata(raw / "GSE76427_series_matrix.txt.gz")
    metadata = metadata[metadata.title.str.startswith("Primary HCC tumor tissue")].copy()
    probe = read_series_matrix(raw / "GSE76427_series_matrix.txt.gz")
    annotation = read_gpl_annotation(raw / "GPL10558.annot.gz").set_index("ID")
    probe = probe.loc[probe.index.intersection(annotation.index), metadata.sample_id]
    probe["gene_symbol"] = annotation.loc[probe.index, "Gene symbol"]
    expression = probe.groupby("gene_symbol").median(numeric_only=True)
    # A small fraction of RSN-normalized background probes are negative.
    # They carry no interpretable abundance information; truncate them at zero
    # before the log transform rather than converting them to missing values.
    expression = _clean_expression(np.log2(expression.clip(lower=0) + 1))
    metadata["sample_type"] = "primary_hcc_tumor"
    return expression.loc[:, metadata.sample_id], metadata


def main() -> None:
    root, _ = parse_project_root()
    raw = root / "data/raw/external_solid_tumor"
    output_rows = []
    loaders = {
        "GSE172356": prepare_gse172356,
        "GSE113255": prepare_gse113255,
        "GSE76427": prepare_gse76427,
    }
    for cohort, loader in loaders.items():
        expression, metadata = loader(raw)
        if expression.columns.duplicated().any() or metadata.sample_id.duplicated().any():
            raise ValueError(f"Duplicate sample identifiers in {cohort}")
        out = root / f"data/processed/{cohort}"
        out.mkdir(parents=True, exist_ok=True)
        expression.to_csv(out / "expression_log2.tsv.gz", sep="\t", compression="gzip")
        metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)
        output_rows.append({
            "cohort": cohort, "n_samples": expression.shape[1],
            "n_genes": expression.shape[0], "itgb1_missing": int(expression.loc["ITGB1"].isna().sum()),
        })
    qc = root / "results/external_solid_tumor"
    qc.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output_rows).to_csv(qc / "preprocessing_qc.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
