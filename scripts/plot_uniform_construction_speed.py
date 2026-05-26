#!/usr/bin/env python3
"""Render a dependency-free SVG construction-speed plot from benchmark CSV."""

import csv
import html
import math
import os
import sys
from collections import defaultdict


SERIES = [
    ("UBS", "ubs_avg_us", "#0072B2"),
    ("NUBS", "nubs_avg_us", "#D55E00"),
    ("MINCO", "minco_avg_us", "#009E73"),
    ("large-scale", "large_scale_avg_us", "#CC79A7"),
]


def usage():
    print(
        "usage: plot_uniform_construction_speed.py "
        "<benchmark.csv> <output.svg>",
        file=sys.stderr,
    )


def read_rows(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            converted = dict(row)
            converted["order"] = int(row["order"])
            converted["piece_count"] = int(row["piece_count"])
            for _, key, _ in SERIES:
                converted[key] = float(row[key])
            rows.append(converted)
    if not rows:
        raise RuntimeError("CSV has no benchmark rows.")
    return rows


def log_map(value, lo, hi, out_lo, out_hi):
    a = math.log10(lo)
    b = math.log10(hi)
    t = (math.log10(value) - a) / (b - a)
    return out_lo + t * (out_hi - out_lo)


def log_ticks(lo, hi, preferred=None):
    if preferred is None:
        preferred = []
    ticks = [v for v in preferred if lo <= v <= hi]
    if ticks:
        return ticks

    result = []
    p_min = int(math.floor(math.log10(lo)))
    p_max = int(math.ceil(math.log10(hi)))
    for power in range(p_min, p_max + 1):
        for base in (1, 2, 5):
            value = base * (10**power)
            if lo <= value <= hi:
                result.append(value)
    return result


def format_tick(value):
    if value >= 1000 and abs(value % 1000) < 1.0e-9:
        return f"{int(value / 1000)}k"
    if value >= 1:
        return str(int(value))
    return f"{value:.1f}"


def svg_text(x, y, text, size=13, anchor="middle", weight="400", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}"'
        f' font-size="{size}" font-weight="{weight}"{transform}>'
        f"{html.escape(text)}</text>"
    )


def render_svg(rows, output_path):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["order"]].append(row)
    orders = sorted(grouped)

    width = 1280
    height = 660
    margin_left = 86
    margin_right = 36
    margin_top = 148
    margin_bottom = 82
    panel_gap = 64
    panel_width = (
        width - margin_left - margin_right - panel_gap * (len(orders) - 1)
    ) / len(orders)
    plot_height = height - margin_top - margin_bottom

    all_piece_counts = [row["piece_count"] for row in rows]
    x_min = min(all_piece_counts)
    x_max = max(all_piece_counts)
    all_values = []
    for row in rows:
        for _, key, _ in SERIES:
            all_values.append(row[key])
    y_min = 10 ** math.floor(math.log10(min(all_values) * 0.8))
    y_max = 10 ** math.ceil(math.log10(max(all_values) * 1.2))
    if y_min <= 0:
        y_min = min(all_values)

    x_ticks = log_ticks(
        x_min,
        x_max,
        [2, 4, 8, 16, 32, 64, 128, 256, 512, 1000, 2000, 5000],
    )
    y_ticks = log_ticks(y_min, y_max)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Inter, Segoe UI, Arial, sans-serif; fill: #1f2933; }",
        ".axis { stroke: #25313f; stroke-width: 1.3; }",
        ".grid { stroke: #d7dee8; stroke-width: 0.8; }",
        ".series { fill: none; stroke-width: 2.5; stroke-linejoin: round; }",
        ".marker { stroke: white; stroke-width: 1.2; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        svg_text(width / 2, 34, "Uniform-Time Trajectory Construction Speed", 22, weight="700"),
        svg_text(
            width / 2,
            58,
            "UBS is the reduced uniform-time B-spline path; NUBS uses the general non-uniform solver with uniform durations.",
            13,
        ),
    ]

    legend_width = 4 * 142
    legend_x = (width - legend_width) / 2
    legend_y = 86
    for label, _, color in SERIES:
        parts.append(
            f'<line x1="{legend_x:.1f}" y1="{legend_y:.1f}" '
            f'x2="{legend_x + 28:.1f}" y2="{legend_y:.1f}" '
            f'stroke="{color}" stroke-width="3"/>'
        )
        parts.append(
            f'<circle class="marker" cx="{legend_x + 14:.1f}" '
            f'cy="{legend_y:.1f}" r="4" fill="{color}"/>'
        )
        parts.append(svg_text(legend_x + 36, legend_y + 4, label, 13, "start"))
        legend_x += 142

    for panel_idx, order in enumerate(orders):
        panel_left = margin_left + panel_idx * (panel_width + panel_gap)
        panel_right = panel_left + panel_width
        panel_top = margin_top
        panel_bottom = panel_top + plot_height
        title = "s = 3  minimum jerk" if order == 3 else "s = 4  minimum snap"

        parts.append(
            f'<rect x="{panel_left:.1f}" y="{panel_top:.1f}" '
            f'width="{panel_width:.1f}" height="{plot_height:.1f}" '
            'fill="#fbfcfe" stroke="#d7dee8"/>'
        )
        parts.append(
            svg_text(
                (panel_left + panel_right) / 2,
                panel_top - 22,
                title,
                16,
                weight="700",
            )
        )

        for tick in y_ticks:
            y = log_map(tick, y_min, y_max, panel_bottom, panel_top)
            parts.append(
                f'<line class="grid" x1="{panel_left:.1f}" y1="{y:.1f}" '
                f'x2="{panel_right:.1f}" y2="{y:.1f}"/>'
            )
            if panel_idx == 0:
                parts.append(svg_text(panel_left - 10, y + 4, format_tick(tick), 12, "end"))

        for tick in x_ticks:
            x = log_map(tick, x_min, x_max, panel_left, panel_right)
            parts.append(
                f'<line class="grid" x1="{x:.1f}" y1="{panel_top:.1f}" '
                f'x2="{x:.1f}" y2="{panel_bottom:.1f}"/>'
            )
            parts.append(svg_text(x, panel_bottom + 22, format_tick(tick), 11))

        parts.append(
            f'<line class="axis" x1="{panel_left:.1f}" y1="{panel_bottom:.1f}" '
            f'x2="{panel_right:.1f}" y2="{panel_bottom:.1f}"/>'
        )
        parts.append(
            f'<line class="axis" x1="{panel_left:.1f}" y1="{panel_top:.1f}" '
            f'x2="{panel_left:.1f}" y2="{panel_bottom:.1f}"/>'
        )

        panel_rows = sorted(grouped[order], key=lambda r: r["piece_count"])
        for label, key, color in SERIES:
            points = []
            for row in panel_rows:
                x = log_map(row["piece_count"], x_min, x_max, panel_left, panel_right)
                y = log_map(row[key], y_min, y_max, panel_bottom, panel_top)
                points.append((x, y))
            point_attr = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
            parts.append(
                f'<polyline class="series" points="{point_attr}" '
                f'stroke="{color}"/>'
            )
            for x, y in points:
                parts.append(
                    f'<circle class="marker" cx="{x:.1f}" cy="{y:.1f}" '
                    f'r="3.7" fill="{color}"/>'
                )

    parts.append(svg_text(width / 2, height - 20, "piece count M (log scale)", 13))
    parts.append(svg_text(24, height / 2, "avg construction time, us (log scale)", 13, rotate=-90))
    parts.append("</svg>")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")


def main():
    if len(sys.argv) != 3:
        usage()
        return 2
    rows = read_rows(sys.argv[1])
    render_svg(rows, sys.argv[2])
    print(f"wrote SVG: {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
