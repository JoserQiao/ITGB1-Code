#!/usr/bin/env python3
"""Add reproducible A/B/C panel labels to supplied GEPIA2 survival figures.

Figures 7 and 8 in the manuscript are portal-export composites. Their source
tables are not part of this repository, so this utility performs only the
deterministic panel-labeling step on the supplied PNGs and never changes the
underlying plotted values.
"""

from argparse import ArgumentParser
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def add_labels(source: Path, output: Path, positions: tuple[tuple[int, int], ...]) -> None:
    image = Image.open(source).convert("RGB")
    draw = ImageDraw.Draw(image)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", max(32, image.width // 32))
    except OSError:
        font = ImageFont.load_default()

    cover_boxes = ((0, 0, 92, 135), (0, 465, 105, 560), (735, 465, 810, 560))
    for label, (x, y), cover in zip(("A", "B", "C"), positions, cover_boxes):
        # The supplied portal composites already contain lower-case tags. Paint
        # only their whitespace margin, leaving heatmap row labels untouched.
        draw.rectangle(cover, fill="white")
        draw.text((x, y), label, fill="black", font=font)

    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output, format="PNG", optimize=True)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--figure", choices=("7", "8"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    # The supplied 1470 x 1035 composites use one upper heatmap and two lower
    # summary panels; Figures 7 and 8 share this layout.
    add_labels(args.input, args.output, ((22, 28), (22, 470), (770, 470)))


if __name__ == "__main__":
    main()
