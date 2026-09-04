#!/usr/bin/env python3
"""Recolour the output-detail text in a rendered nf-metro SVG.

nf-metro has no per-station text-colour directive, so every station label
(tool name or output-detail box) renders in the same theme colour. This
patches the rendered SVG in place, giving every station whose id ends in
``_notes`` (our output-detail boxes) a distinct fill so they read as
reference annotations rather than pipeline steps.

Usage: colour_output_notes.py <svg-file> [colour]
Run after every `nf-metro render` of assembly_binning_route_map.mmd.
"""

import re
import sys

DEFAULT_COLOUR = "#94a3b8"  # matches the 'outputs' line colour in the .mmd


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(1)
    path = sys.argv[1]
    colour = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_COLOUR

    svg = open(path, encoding="utf-8").read()
    style = (
        f'\n<style>text.nf-metro-station-label[data-station-id$="_notes"] '
        f"{{ fill: {colour} !important; }}</style>\n"
    )
    if "</defs>" not in svg:
        raise SystemExit("expected a <defs>...</defs> block to anchor the style insert")
    svg = svg.replace("</defs>", "</defs>" + style, 1)
    open(path, "w", encoding="utf-8").write(svg)
    print(f"Coloured *_notes station labels {colour} in {path}")


if __name__ == "__main__":
    main()
