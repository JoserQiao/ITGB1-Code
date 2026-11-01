"""Download Cell Model Passports pharmacogenomic data."""

from __future__ import annotations

import json
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from scripts.core.cli import parse_project_root


BASE = "https://api.cellmodelpassports.sanger.ac.uk"


def session_with_retries() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=5, connect=5, read=5, backoff_factor=1.0,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=("GET",))
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def fetch(path: str, params: dict, dest: Path) -> dict:
    response = session_with_retries().get(BASE + path, params=params, timeout=300)
    response.raise_for_status()
    payload = response.json()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(payload), encoding="utf-8")
    return payload


def fetch_collection(path: str, params: dict, dest: Path, page_size: int = 500) -> None:
    """Download a JSON:API collection while preserving included resources."""
    session = session_with_retries()
    data, included = [], {}
    page = 1
    while True:
        query = {**params, "page[size]": page_size, "page[number]": page}
        response = session.get(BASE + path, params=query, timeout=300)
        response.raise_for_status()
        payload = response.json()
        data.extend(payload.get("data", []))
        for item in payload.get("included", []):
            included[(item.get("type"), str(item.get("id")))] = item
        if not payload.get("links", {}).get("next"):
            break
        page += 1
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps({"data": data, "included": list(included.values())}), encoding="utf-8")


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
    # The public API is more reliable for this endpoint with small pages.
    # Model IDs are already present in each sample relationship, so only the
    # tissue resource needs to be included.
    fetch_collection("/samples", {"include": "tissue"}, raw / "samples.json", page_size=50)
    if not (raw / "drugs.jsonl").exists():
        fetch_pages("/drugs", {}, raw / "drugs.jsonl", page_size=1000)
    fetch_pages("/datasets/nlmes", {
        "fields[nlme_curves]": "drug_id,model_id,ln_ic50,auc,z_score",
    }, raw / "GDSC_fitted_curves.jsonl", page_size=10000, resume=True)


if __name__ == "__main__":
    main()
