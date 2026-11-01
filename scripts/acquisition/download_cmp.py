"""Download Cell Model Passports pharmacogenomic data."""

from __future__ import annotations

import json
from pathlib import Path

import requests

from scripts.core.cli import parse_project_root


BASE = "https://api.cellmodelpassports.sanger.ac.uk"


def fetch(path: str, params: dict, dest: Path) -> dict:
    response = requests.get(BASE + path, params=params, timeout=300)
    response.raise_for_status()
    payload = response.json()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(payload), encoding="utf-8")
    return payload


def fetch_pages(path: str, params: dict, dest: Path, page_size: int = 10000,
                resume: bool = False) -> None:
    # JSON Lines is written incrementally so an interrupted 240k-row download
    # can resume without keeping the entire response set in memory.
    existing_rows = 0
    if resume and dest.exists():
        with dest.open(encoding="utf-8") as old:
            existing_rows = sum(1 for _ in old)
        page = existing_rows // page_size + 1
        mode = "a"
    else:
        page = 1
        mode = "w"
    with dest.open(mode, encoding="utf-8") as handle:
        while True:
            query = {**params, "page[size]": page_size, "page[number]": page}
            response = requests.get(BASE + path, params=query, timeout=300)
            response.raise_for_status()
            payload = response.json()
            for item in payload.get("data", []):
                handle.write(json.dumps(item) + "\n")
            if not payload.get("links", {}).get("next"):
                break
            page += 1


def main() -> None:
    root, _ = parse_project_root()
    raw = root / "data/raw/CellModelPassports"
    raw.mkdir(parents=True, exist_ok=True)
    fetch("/genes/SIDG13308/datasets/rnaseq", {"page[size]": 3000}, raw / "ITGB1_rnaseq.json")
    if not (raw / "drugs.jsonl").exists():
        fetch_pages("/drugs", {}, raw / "drugs.jsonl", page_size=1000)
    fetch_pages("/datasets/nlmes", {
        "fields[nlme_curves]": "drug_id,model_id,ln_ic50,auc,z_score",
    }, raw / "GDSC_fitted_curves.jsonl", page_size=10000, resume=True)


if __name__ == "__main__":
    main()
