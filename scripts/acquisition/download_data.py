"""Download public datasets listed in the project manifest."""

from __future__ import annotations

import csv
import hashlib
import json
import ssl
import subprocess
import urllib.request
from pathlib import Path

import certifi

from scripts.core.cli import parse_project_root, project_root_parser


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".part")
    headers = {}
    mode = "wb"
    existing = partial.stat().st_size if partial.exists() else 0
    if existing:
        headers["Range"] = f"bytes={partial.stat().st_size}-"
        mode = "ab"
    request = urllib.request.Request(url, headers=headers)
    context = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(request, timeout=120, context=context) as response, partial.open(mode) as out:
        # Some signed object-store redirects ignore Range. Restart instead of
        # appending a duplicate full object to a partial download.
        if existing and getattr(response, "status", 200) != 206:
            out.close()
            partial.unlink()
            return download(url, dest)
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            out.write(block)
    partial.replace(dest)


def main() -> None:
    parser = project_root_parser()
    parser.add_argument("--manifest", default="config/datasets.tsv")
    parser.add_argument("--only", nargs="*", default=[])
    root, args = parse_project_root(parser)
    records = list(csv.DictReader((root / args.manifest).open(encoding="utf-8"), delimiter="\t"))
    audit = []
    for rec in records:
        if args.only and rec["dataset"] not in args.only:
            continue
        # Paged Cell Model Passports endpoints are handled by download_cmp.py,
        # which preserves JSON:API pagination and included tissue resources.
        if (rec["role"].startswith("pharmacogenomic_")
                or rec["role"] in {"single_cell_external_validation",
                                   "pathology_composition_validation",
                                   "protein_validation"}):
            continue
        dest = root / rec["local_path"]
        if rec["source_url"].endswith(".git"):
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                subprocess.run(["git", "clone", "--depth", "1", rec["source_url"], str(dest)], check=True)
        elif not dest.exists():
            download(rec["source_url"], dest)
        audit.append({**rec, "bytes": dest.stat().st_size if dest.is_file() else None,
                      "sha256": sha256(dest) if dest.is_file() else None})
    log_dir = root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "download_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
