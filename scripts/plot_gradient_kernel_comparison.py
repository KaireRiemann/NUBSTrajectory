#!/usr/bin/env python3
"""Plot the NUBS gradient-kernel evolution against MINCO."""

import csv
import html
import math
import os
import sys
from collections import defaultdict


SERIES = [
    ("NUBS scalar Local-AD", "scalar", "#999999", "5 4"),
    ("NUBS fused Local-Jet", "fused", "#009E73", "7 4"),
    ("NUBS exact Local-Jet", "exact", "#0072B2", None),
    ("NUBS scalar reverse", "reverse", "#CC79A7", "9 4"),
    ("MINCO", "minco", "#D55E00", "2 3"),
]
STAGES = [("Direct internal gradient construction", "direct_gradient"),
          ("Full waypoint/time gradient propagation", "full_propagation")]


def usage():
    print("usage: plot_gradient_kernel_comparison.py <scalar.csv> <fused.csv> <exact.csv> <scalar_reverse.csv> <output.svg>",
          file=sys.stderr)


def read_rows(path):
    with open(path, newline="") as csv_file:
        return list(csv.DictReader(csv_file))


def log_map(value, lo, hi, output_lo, output_hi):
    return output_lo + ((math.log10(value) - math.log10(lo)) /
                        (math.log10(hi) - math.log10(lo)) *
                        (output_hi - output_lo))


def text(x, y, value, size=12, anchor="middle", weight="400", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
            f'font-size="{size}" font-weight="{weight}"{transform}>'
            f'{html.escape(value)}</text>')


def tick_label(value):
    return f"{value / 1000:g}k" if value >= 1000 else f"{value:g}"


def render(paths, output_path):
    scalar_rows, fused_rows, exact_rows, reverse_rows = (
        read_rows(path) for path in paths)
    raw_sets = {"scalar": scalar_rows, "fused": fused_rows,
                "exact": exact_rows, "reverse": reverse_rows}
    grouped = defaultdict(dict)
    all_times = []
    for kernel, rows in raw_sets.items():
        for row in rows:
            key = (int(row["order"]), int(row["piece_count"]))
            grouped[key][kernel] = row
            for _, stage in STAGES:
                all_times.append(float(row[f"nubs_{stage}_us"]))
            if kernel == "reverse":
                for _, stage in STAGES:
                    all_times.append(float(row[f"minco_{stage}_us"]))

    orders = sorted({key[0] for key in grouped})
    y_lo = 10 ** math.floor(math.log10(min(all_times)) - 0.15)
    y_hi = 10 ** math.ceil(math.log10(max(all_times)) + 0.15)
    y_ticks = [10 ** exponent for exponent in range(
        int(math.ceil(math.log10(y_lo))), int(math.floor(math.log10(y_hi))) + 1)]

    width, height = 1450, 700
    left_margin, right_margin = 90, 25
    top_margin, bottom_margin = 110, 65
    column_gap, row_gap = 42, 56
    panel_width = (width - left_margin - right_margin -
                   column_gap * (len(orders) - 1)) / len(orders)
    panel_height = (height - top_margin - bottom_margin - row_gap) / 2
    pieces = [2, 4, 8, 16, 32, 64]

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; fill: #17202a; }",
        ".axis { stroke: #17202a; stroke-width: 1.15; }",
        ".grid { stroke: #d6dde6; stroke-width: 0.75; }",
        ".series { fill: none; stroke-width: 2.35; stroke-linejoin: round; }",
        ".marker { stroke: white; stroke-width: 1.0; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        text(width / 2, 27, "Evolution of NUBS gradient kernels versus MINCO", 18, weight="700"),
    ]

    legend_x = 55
    for label, _, color, dash in SERIES:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(f'<line x1="{legend_x}" y1="56" x2="{legend_x + 33}" y2="56" stroke="{color}" stroke-width="2.8"{dash_attr}/>')
        parts.append(text(legend_x + 40, 60, label, 11, "start", "700"))
        legend_x += 270

    for row_index, (stage_title, stage_key) in enumerate(STAGES):
        top = top_margin + row_index * (panel_height + row_gap)
        bottom = top + panel_height
        for column_index, order in enumerate(orders):
            left = left_margin + column_index * (panel_width + column_gap)
            right = left + panel_width
            parts.append(f'<rect x="{left:.1f}" y="{top:.1f}" width="{panel_width:.1f}" height="{panel_height:.1f}" fill="#fbfcfe" stroke="#d6dde6"/>')
            if row_index == 0:
                parts.append(text(left + 7, top - 12, f"s = {order}", 14, "start", "700"))
            if column_index == 0:
                parts.append(text(left + 5, top + 18, stage_title, 11, "start", "700"))
            for tick in y_ticks:
                y = log_map(tick, y_lo, y_hi, bottom, top)
                parts.append(f'<line class="grid" x1="{left:.1f}" y1="{y:.1f}" x2="{right:.1f}" y2="{y:.1f}"/>')
                if column_index == 0:
                    parts.append(text(left - 9, y + 4, tick_label(tick), 10, "end"))
            for piece in pieces:
                x = log_map(piece, 2, 64, left, right)
                parts.append(f'<line class="grid" x1="{x:.1f}" y1="{top:.1f}" x2="{x:.1f}" y2="{bottom:.1f}"/>')
                if row_index == len(STAGES) - 1:
                    parts.append(text(x, bottom + 20, str(piece), 10))
            parts.append(f'<line class="axis" x1="{left:.1f}" y1="{bottom:.1f}" x2="{right:.1f}" y2="{bottom:.1f}"/>')
            parts.append(f'<line class="axis" x1="{left:.1f}" y1="{top:.1f}" x2="{left:.1f}" y2="{bottom:.1f}"/>')

            for label, kernel, color, dash in SERIES:
                points = []
                for piece in pieces:
                    row = grouped[(order, piece)]["reverse" if kernel == "minco" else kernel]
                    prefix = "minco" if kernel == "minco" else "nubs"
                    value = float(row[f"{prefix}_{stage_key}_us"])
                    points.append((log_map(piece, 2, 64, left, right),
                                   log_map(value, y_lo, y_hi, bottom, top)))
                dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
                polyline = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
                parts.append(f'<polyline class="series" points="{polyline}" stroke="{color}"{dash_attr}/>')
                for x, y in points:
                    parts.append(f'<circle class="marker" cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{color}"/>')

    parts.append(text(width / 2, height - 17, "number of pieces M", 13, weight="700"))
    parts.append(text(24, height / 2, "average time per call (us, log scale)", 13,
                      weight="700", rotate=-90))
    parts.append("</svg>")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as output:
        output.write("\n".join(parts) + "\n")


def main():
    if len(sys.argv) != 6:
        usage()
        return 2
    render(sys.argv[1:5], sys.argv[5])
    print(f"wrote SVG: {sys.argv[5]}")


if __name__ == "__main__":
    raise SystemExit(main())
