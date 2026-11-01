"""Download the public IMvigor210 DESeq2 object used in this study."""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

from scripts.core.cli import parse_project_root


URL = (
    "https://raw.githubusercontent.com/BioInfoCloud/"
    "IMvigor210CoreBiologies/master/data/cds.RData"
)


def main() -> None:
    root, _ = parse_project_root()
    output = root / "data/raw/IMvigor210/cds.RData"
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.stat().st_size > 1_000_000:
        print(f"Using existing {output}")
        return
    temporary = output.with_suffix(".RData.part")
    with urlopen(URL, timeout=120) as response, temporary.open("wb") as handle:
        while block := response.read(1024 * 1024):
            handle.write(block)
    if temporary.stat().st_size < 1_000_000:
        raise RuntimeError("Downloaded IMvigor210 object is unexpectedly small")
    temporary.replace(output)
    print(f"Downloaded {output} ({output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
