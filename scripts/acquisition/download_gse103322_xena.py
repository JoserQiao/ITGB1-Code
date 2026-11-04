"""Acquire annotated GSE103322 single-cell ITGB1 data from UCSC Xena."""

from __future__ import annotations

import warnings
from io import StringIO

import pandas as pd
import requests
import xenaPython as xena

from scripts.core.cli import parse_project_root


HOST = "https://ucscpublic.xenahubs.net"
EXPRESSION_DATASET = "HNSCC/GSE103322_HNSCC_exp.txt"
METADATA_URL = f"{HOST}/download/HNSCC/GSE103322_cell_metadata.txt"
SAMPLE_METADATA_URL = f"{HOST}/download/HNSCC/GSE103322_sample_metadata.txt"


def main():
    warnings.filterwarnings("ignore", category=FutureWarning, module="xenaPython")
    root, _ = parse_project_root()
    raw = root / "data/raw/GSE103322"
    processed = root / "data/processed/GSE103322"
    raw.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    response = requests.get(METADATA_URL, timeout=180)
    response.raise_for_status()
    metadata_path = raw / "GSE103322_cell_metadata.txt"
    metadata_path.write_bytes(response.content)
    metadata = pd.read_csv(StringIO(response.text), sep="\t")
    metadata = metadata.rename(columns={metadata.columns[0]: "cell_id"})

    sample_response = requests.get(SAMPLE_METADATA_URL, timeout=180)
    sample_response.raise_for_status()
    sample_path = raw / "GSE103322_sample_metadata.txt"
    sample_path.write_bytes(sample_response.content)
    sample_metadata = pd.read_csv(StringIO(sample_response.text), sep="\t").rename(
        columns={"cell": "cell_id"})
    metadata = metadata.merge(sample_metadata, on="cell_id", how="left")

    available = set(xena.dataset_samples(HOST, EXPRESSION_DATASET, None))
    cells = [cell for cell in metadata.cell_id.astype(str) if cell in available]
    values = []
    for start in range(0, len(cells), 1000):
        chunk = cells[start:start + 1000]
        result = xena.dataset_probe_values(HOST, EXPRESSION_DATASET, chunk, ["ITGB1"])
        if not result or not result[1] or len(result[1][0]) != len(chunk):
            raise RuntimeError("UCSC Xena returned an incomplete GSE103322 ITGB1 query")
        values.extend(result[1][0])
    expression = pd.DataFrame({"cell_id": cells, "ITGB1": values})
    metadata = metadata[metadata.cell_id.isin(cells)].copy()
    metadata.to_csv(processed / "metadata.tsv", sep="\t", index=False)
    expression.to_csv(processed / "itgb1_expression.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
