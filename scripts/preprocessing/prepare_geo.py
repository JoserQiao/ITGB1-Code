"""Prepare GEO expression matrices and clinical metadata."""

from __future__ import annotations

import gzip
import re
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.core.cli import parse_project_root


def parse_miniml(path: Path) -> pd.DataFrame:
    with tarfile.open(path, "r:gz") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith(".xml"))
        root = ET.parse(tf.extractfile(member)).getroot()
    rows = []
    for sample in root.iter():
        if sample.tag.rsplit("}", 1)[-1] != "Sample":
            continue
        row = {"gsm": sample.attrib.get("iid", "")}
        for node in sample.iter():
            tag = node.tag.rsplit("}", 1)[-1]
            value = (node.text or "").strip()
            if not value:
                continue
            if tag in {"Title", "Source"} and tag.lower() not in row:
                row[tag.lower()] = value
            elif tag == "Characteristics":
                key = re.sub(r"[^a-z0-9]+", "_", node.attrib.get("tag", "characteristic").lower()).strip("_")
                row[key] = value
        rows.append(row)
    return pd.DataFrame(rows)


def collapse_duplicate_genes(frame: pd.DataFrame) -> pd.DataFrame:
    frame.index = frame.index.astype(str).str.upper().str.strip()
    frame = frame.loc[frame.index != ""]
    return frame.groupby(level=0).mean(numeric_only=True)


def map_entrez_to_symbol(frame: pd.DataFrame, root: Path) -> pd.DataFrame:
    annotation = pd.read_csv(root / "data/raw/annotation/Homo_sapiens.gene_info.gz", sep="\t", dtype=str)
    mapping = annotation.set_index("GeneID")["Symbol"].to_dict()
    frame.index = [mapping.get(str(x), str(x)) for x in frame.index]
    return collapse_duplicate_genes(frame)


def prepare_gse78220(root: Path, meta: pd.DataFrame) -> None:
    src = root / "data/raw/GSE78220/GSE78220_PatientFPKM.xlsx"
    expr = pd.read_excel(src, sheet_name="FPKM", index_col=0)
    expr.columns = [re.sub(r"\.(baseline|OnTx)$", "", str(x), flags=re.I) for x in expr.columns]
    expr = collapse_duplicate_genes(expr)
    # Pt27A/Pt27B are two baseline lesions from one patient; use patient-level mean.
    expr.columns = [re.sub(r"([A-Z])$", "", x) if re.match(r"^Pt\d+[A-Z]$", x) else x for x in expr.columns]
    expr = expr.T.groupby(level=0).mean().T
    response_col = next(c for c in meta.columns if "response" in c)
    pheno = meta[["patient_id", response_col]].drop_duplicates().rename(
        columns={response_col: "response_raw"})
    pheno["patient_id"] = pheno.patient_id.str.replace(r"[A-Z]$", "", regex=True)
    pheno = pheno.drop_duplicates("patient_id")
    pheno["response_binary"] = pheno.response_raw.str.contains("Complete|Partial", case=False, regex=True).astype(int)
    keep = [x for x in expr.columns if x in set(pheno.patient_id)]
    expr = np.log2(expr[keep] + 1)
    pheno = pheno.set_index("patient_id").loc[keep].reset_index()
    out = root / "data/processed/GSE78220"
    out.mkdir(parents=True, exist_ok=True)
    if not (out / "expression_log2_fpkm.tsv.gz").exists():
        expr.to_csv(out / "expression_log2_fpkm.tsv.gz", sep="\t", compression="gzip")
    pheno.to_csv(out / "metadata_v2.tsv", sep="\t", index=False)


def prepare_gse91061(root: Path, meta: pd.DataFrame) -> None:
    src = root / "data/raw/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    expr = pd.read_csv(src, index_col=0)
    expr = collapse_duplicate_genes(expr)
    if "ITGB1" not in expr.index:
        expr = map_entrez_to_symbol(expr, root)
    meta = meta.rename(columns={"visit_pre_or_on_treatment": "visit", "response": "response_raw"})
    meta["sample_name"] = meta.title.str.extract(r"(Pt\d+_(?:Pre|On)_[A-Za-z0-9-]+)", expand=False)
    meta["patient_id"] = meta.title.str.extract(r"(Pt\d+)", expand=False)
    pre = meta[meta.visit.str.contains("Pre", case=False, na=False)].copy()
    # RECIST labels must be explicit. Samples marked UNK are not assumed to be
    # non-responders, because that silently introduces outcome misclassification.
    pre = pre[pre.response_raw.str.contains(r"CR|PR|PD|SD|complete|partial|progressive|stable",
                                             case=False, regex=True, na=False)].copy()
    pre["response_binary"] = pre.response_raw.str.contains(r"CR|PR|complete|partial", case=False, regex=True).astype(int)
    # Match exact sample identifiers first, then GSM if GEO column labels are used.
    mapping = {}
    for c in expr.columns:
        hit = pre[(pre.sample_name == c) | (pre.gsm == c) | pre.title.str.contains(re.escape(str(c)), na=False)]
        if len(hit) == 1:
            mapping[c] = hit.iloc[0].patient_id
    if not mapping:
        raise ValueError("GSE91061 expression columns could not be matched to MINiML metadata")
    expr = expr[list(mapping)].copy()
    expr.columns = [mapping[x] for x in expr.columns]
    expr = expr.T.groupby(level=0).mean().T
    pheno = pre.drop_duplicates("patient_id").set_index("patient_id").loc[expr.columns].reset_index()
    if "index" in pheno.columns and "patient_id" not in pheno.columns:
        pheno = pheno.rename(columns={"index": "patient_id"})
    expr = np.log2(expr + 1)
    out = root / "data/processed/GSE91061"
    out.mkdir(parents=True, exist_ok=True)
    if not (out / "expression_log2_fpkm.tsv.gz").exists():
        expr.to_csv(out / "expression_log2_fpkm.tsv.gz", sep="\t", compression="gzip")
    pheno.to_csv(out / "metadata_v2.tsv", sep="\t", index=False)


def prepare_rectal_metadata(root: Path, meta: pd.DataFrame) -> None:
    response = next(c for c in meta.columns if "treatment_response" in c)
    meta = meta.rename(columns={response: "response_raw"})
    # The series includes 105 pretreatment biopsies and 12 post-CRT resections.
    # Only pretreatment samples are eligible for treatment-response prediction.
    time_col = next((c for c in meta.columns if c in {"time", "timepoint"}), None)
    if time_col is None:
        raise ValueError("GSE190826 pretreatment time field was not found")
    meta = meta[meta[time_col].str.contains("Pre", case=False, na=False)].copy()
    meta["response_binary"] = meta.response_raw.str.fullmatch("pCR", case=False).astype(int)
    out = root / "data/processed/GSE190826"
    out.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out / "metadata_v2.tsv", sep="\t", index=False)


def ensembl_symbol_map(root: Path) -> dict[str, str]:
    annotation = pd.read_csv(root / "data/raw/annotation/Homo_sapiens.gene_info.gz", sep="\t", dtype=str)
    pairs = {}
    for row in annotation.itertuples(index=False):
        dbx = getattr(row, "dbXrefs", "") or ""
        for token in dbx.split("|"):
            if token.startswith("Ensembl:ENSG"):
                pairs[token.split(":", 1)[1]] = row.Symbol
    return pairs


def prepare_rectal_counts(root: Path) -> None:
    tar_path = root / "data/raw/GSE190826/GSE190826_RAW.tar"
    if not tar_path.exists():
        return
    mapping = ensembl_symbol_map(root)
    series = []
    with tarfile.open(tar_path) as tf:
        for member in tf.getmembers():
            if not member.name.endswith(".counts.txt.gz"):
                continue
            with gzip.GzipFile(fileobj=tf.extractfile(member)) as gh:
                tab = pd.read_csv(gh, sep="\t", usecols=[0, 6])
            gsm = member.name.split("_", 1)[0]
            s = tab.set_index("Geneid").iloc[:, 0]
            s.name = gsm
            series.append(s)
    counts = pd.concat(series, axis=1).fillna(0)
    counts.index = [mapping.get(str(x).split(".")[0], str(x).split(".")[0]) for x in counts.index]
    counts = counts.groupby(level=0).sum()
    out = root / "data/processed/GSE190826"
    if not (out / "counts.tsv.gz").exists():
        counts.to_csv(out / "counts.tsv.gz", sep="\t", compression="gzip")


def main() -> None:
    root, _ = parse_project_root()
    for accession in ["GSE78220", "GSE91061", "GSE190826"]:
        meta = parse_miniml(root / f"data/raw/{accession}/{accession}_family.xml.tgz")
        out = root / f"data/processed/{accession}"
        out.mkdir(parents=True, exist_ok=True)
        meta.to_csv(out / "metadata_all_v2.tsv", sep="\t", index=False)
        if accession == "GSE78220":
            prepare_gse78220(root, meta)
        elif accession == "GSE91061":
            prepare_gse91061(root, meta)
        else:
            prepare_rectal_metadata(root, meta)
            prepare_rectal_counts(root)


if __name__ == "__main__":
    main()
