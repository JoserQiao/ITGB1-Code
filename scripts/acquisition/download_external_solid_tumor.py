"""Download the prespecified external solid-tumor validation datasets."""

from __future__ import annotations

import pandas as pd

from scripts.acquisition.download_data import download
from scripts.core.cli import parse_project_root


ROLES = {"external_pathway_validation", "external_pathway_metadata", "platform_annotation"}


def main() -> None:
    root, _ = parse_project_root()
    registry = pd.read_csv(root / "config/datasets.tsv", sep="\t")
    selected = registry[registry.role.isin(ROLES)]
    for row in selected.itertuples(index=False):
        destination = root / row.local_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and destination.stat().st_size > 0:
            continue
        print(f"Downloading {row.dataset} -> {destination}")
        download(row.source_url, destination)


if __name__ == "__main__":
    main()
