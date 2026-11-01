"""Shared command-line parsing helpers."""

from __future__ import annotations

import argparse
from pathlib import Path


def project_root_parser(description: str | None = None) -> argparse.ArgumentParser:
    """Create a parser with the project-root option used by every command."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--project-root", default=".")
    return parser


def parse_project_root(
    parser: argparse.ArgumentParser | None = None,
) -> tuple[Path, argparse.Namespace]:
    """Parse command-line arguments and return a resolved project root."""
    parser = parser or project_root_parser()
    args = parser.parse_args()
    return Path(args.project_root).resolve(), args
