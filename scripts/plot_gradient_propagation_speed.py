#!/usr/bin/env python3
"""Render scalar-reverse NUBS and MINCO gradient timing panels."""

import csv
import html
import math
import os
import sys
from collections import defaultdict


SERIES = [
    ("NUBS derivative-control reverse", "nubs", "#0072B2", None),
    ("MINCO", "minco", "#D55E00", "8 5"),
]
STAGES = [
    ("Direct internal gradient construction", "direct_gradient"),
    ("Full waypoint/time gradient propagation", "full_propagation"),
]


def usage():
    print("usage: plot_gradient_propagation_speed.py <benchmark.csv> <output.svg>", file=sys.stderr)


def read_rows(path):
    rows = []
    with open(path, newline="") as csv_file:
        for raw in csv.DictReader(csv_file):
            rows.append({
                "order": int(raw["order"]),
                "piece_count": int(raw["piece_count"]),
                "nubs_direct_gradient_us": float(raw["nubs_direct_gradient_us"]),
                "minco_direct_gradient_us": float(raw["minco_direct_gradient_us"]),
                "nubs_full_propagation_us": float(raw["nubs_full_propagation_us"]),
                "minco_full_propagation_us": float(raw["minco_full_propagation_us"]),
            })
    if not rows:
        raise RuntimeError("CSV has no rows")
    return rows


def log_map(value, lo, hi, output_lo, output_hi):
    return output_lo + (math.log10(value) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * (output_hi - output_lo)


def svg_text(x, y, text, size=12, anchor="middle", weight="400", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}"{transform}>'
        f"{html.escape(text)}</text>"
    )


def tick_label(value):
    if value >= 1000:
        return f"{value / 1000:g}k"
    if value >= 1:
        return f"{value:g}"
    return f"{value:g}"


def render(rows, output_path):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["order"]].append(row)
    orders = sorted(grouped)
    all_times = [
        row[f"{method}_{stage}_us"]
        for row in rows
        for _, stage in STAGES
        for _, method, _, _ in SERIES
    ]
    y_lo = 10 ** math.floor(math.log10(min(all_times)) - 0.15)
    y_hi = 10 ** math.ceil(math.log10(max(all_times)) + 0.15)
    y_ticks = [10 ** e for e in range(int(math.ceil(math.log10(y_lo))), int(math.floor(math.log10(y_hi))) + 1)]

    width = 1230
    height = 680
    left_margin = 86
    right_margin = 28
    top_margin = 95
    bottom_margin = 65
    column_gap = 44
    row_gap = 58
    panel_width = (width - left_margin - right_margin - column_gap * (len(orders) - 1)) / len(orders)
    panel_height = (height - top_margin - bottom_margin - row_gap) / 2
    x_ticks = [2, 4, 8, 16, 32, 64]

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; fill: #17202a; }",
        ".axis { stroke: #17202a; stroke-width: 1.15; }",
        ".grid { stroke: #d6dde6; stroke-width: 0.75; }",
        ".series { fill: none; stroke-width: 2.4; stroke-linejoin: round; }",
        ".marker { stroke: white; stroke-width: 1.0; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        svg_text(width / 2, 27, "NUBS scalar-reverse versus MINCO gradient timing", 18, weight="700"),
    ]

    legend_x = width / 2 - 118
    for label, _, color, dash in SERIES:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(f'<line x1="{legend_x:.1f}" y1="{54:.1f}" x2="{legend_x + 34:.1f}" y2="54" stroke="{color}" stroke-width="2.8"{dash_attr}/>')
        parts.append(svg_text(legend_x + 43, 58, label, 12, "start", "700"))
        legend_x += 176

    for row_index, (stage_title, stage_key) in enumerate(STAGES):
        top = top_margin + row_index * (panel_height + row_gap)
        bottom = top + panel_height
        for column_index, order in enumerate(orders):
            left = left_margin + column_index * (panel_width + column_gap)
            right = left + panel_width
            panel_rows = sorted(grouped[order], key=lambda row: row["piece_count"])
            parts.append(f'<rect x="{left:.1f}" y="{top:.1f}" width="{panel_width:.1f}" height="{panel_height:.1f}" fill="#fbfcfe" stroke="#d6dde6"/>')
            if row_index == 0:
                parts.append(svg_text(left + 7, top - 12, f"s = {order}", 14, "start", "700"))
            if column_index == 0:
                parts.append(svg_text(left + 5, top + 18, stage_title, 11, "start", "700"))

            for tick in y_ticks:
                y = log_map(tick, y_lo, y_hi, bottom, top)
                parts.append(f'<line class="grid" x1="{left:.1f}" y1="{y:.1f}" x2="{right:.1f}" y2="{y:.1f}"/>')
                if column_index == 0:
                    parts.append(svg_text(left - 9, y + 4, tick_label(tick), 10, "end"))
            for tick in x_ticks:
                x = log_map(tick, 2, 64, left, right)
                parts.append(f'<line class="grid" x1="{x:.1f}" y1="{top:.1f}" x2="{x:.1f}" y2="{bottom:.1f}"/>')
                if row_index == len(STAGES) - 1:
                    parts.append(svg_text(x, bottom + 20, str(tick), 10))
            parts.append(f'<line class="axis" x1="{left:.1f}" y1="{bottom:.1f}" x2="{right:.1f}" y2="{bottom:.1f}"/>')
            parts.append(f'<line class="axis" x1="{left:.1f}" y1="{top:.1f}" x2="{left:.1f}" y2="{bottom:.1f}"/>')

            for _, method, color, dash in SERIES:
                key = f"{method}_{stage_key}_us"
                points = [
                    (log_map(row["piece_count"], 2, 64, left, right),
                     log_map(row[key], y_lo, y_hi, bottom, top))
                    for row in panel_rows
                ]
                dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
                parts.append('<polyline class="series" points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in points) + f'" stroke="{color}"{dash_attr}/>')
                for x, y in points:
                    parts.append(f'<circle class="marker" cx="{x:.1f}" cy="{y:.1f}" r="3.4" fill="{color}"/>')

    parts.append(svg_text(width / 2, height - 17, "number of pieces M", 13, weight="700"))
    parts.append(svg_text(23, height / 2, "average time per call (us, log scale)", 13, weight="700", rotate=-90))
    parts.append("</svg>")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as output:
        output.write("\n".join(parts) + "\n")


def main():
    if len(sys.argv) != 3:
        usage()
        return 2
    render(read_rows(sys.argv[1]), sys.argv[2])
    print(f"wrote SVG: {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
