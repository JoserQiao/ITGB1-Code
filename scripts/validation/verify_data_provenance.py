"""Verify local artifacts against the versioned SHA256 provenance manifest."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from scripts.core.cli import parse_project_root, project_root_parser


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = project_root_parser("Verify SHA256 provenance records")
    parser.add_argument("--strict", action="store_true")
    root, args = parse_project_root(parser)

    manifest = pd.read_csv(root / "config/data_provenance.tsv", sep="\t", dtype=str)
    rows = []
    for record in manifest.to_dict("records"):
        artifact = root / record["artifact_path"]
        observed = sha256(artifact) if artifact.is_file() else ""
        status = "verified" if observed == record["sha256"] else (
            "missing" if not artifact.is_file() else "mismatch"
        )
        rows.append({
            "dataset": record["dataset"],
            "artifact_path": record["artifact_path"],
            "expected_sha256": record["sha256"],
            "observed_sha256": observed,
            "status": status,
        })

    report = pd.DataFrame(rows)
    output = root / "results/validation"
    output.mkdir(parents=True, exist_ok=True)
    report.to_csv(output / "data_provenance_checks.tsv", sep="\t", index=False)
    counts = report.status.value_counts().to_dict()
    print(f"Provenance verification: {counts}")
    if args.strict and report.status.ne("verified").any():
        failed = report.loc[report.status.ne("verified"), ["dataset", "status"]]
        raise RuntimeError("Provenance verification failed:\n" + failed.to_string(index=False))


if __name__ == "__main__":
    main()
