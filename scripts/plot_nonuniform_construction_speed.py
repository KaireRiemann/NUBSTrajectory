#!/usr/bin/env python3
"""Render a NUBS-vs-MINCO non-uniform construction-speed plot from CSV."""

import csv
import html
import math
import os
import sys
from collections import defaultdict


SERIES = [
    ("NUBS", "nubs_avg_us", "#0072B2", None),
    ("MINCO", "minco_avg_us", "#D55E00", "8 5"),
]


def usage():
    print(
        "usage: plot_nonuniform_construction_speed.py "
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
            for _, key, _, _ in SERIES:
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
    preferred = preferred or []
    ticks = [v for v in preferred if lo <= v <= hi]
    if ticks:
        return ticks
    result = []
    for power in range(int(math.floor(math.log10(lo))),
                       int(math.ceil(math.log10(hi))) + 1):
        for base in (1, 2, 5):
            value = base * (10 ** power)
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
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}"{transform}>'
        f"{html.escape(text)}</text>"
    )


def render_svg(rows, output_path):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["order"]].append(row)
    orders = sorted(grouped)

    width = 980
    height = 440
    margin_left = 70
    margin_right = 28
    margin_top = 76
    margin_bottom = 68
    panel_gap = 58
    panel_width = (
        width - margin_left - margin_right - panel_gap * (len(orders) - 1)
    ) / len(orders)
    plot_height = height - margin_top - margin_bottom

    piece_counts = [row["piece_count"] for row in rows]
    x_min = min(piece_counts)
    x_max = max(piece_counts)
    values = [row[key] for row in rows for _, key, _, _ in SERIES]
    y_min = 10 ** math.floor(math.log10(min(values) * 0.8))
    y_max = 10 ** math.ceil(math.log10(max(values) * 1.2))

    x_ticks = log_ticks(x_min, x_max, [2, 4, 8, 16, 32, 64])
    y_ticks = log_ticks(y_min, y_max)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; fill: #17202a; }",
        ".axis { stroke: #17202a; stroke-width: 1.2; }",
        ".grid { stroke: #d6dde6; stroke-width: 0.75; }",
        ".series { fill: none; stroke-width: 2.5; stroke-linejoin: round; }",
        ".marker { stroke: white; stroke-width: 1.1; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
    ]

    legend_x = width / 2 - 118
    legend_y = 34
    for label, _, color, dash in SERIES:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(
            f'<line x1="{legend_x:.1f}" y1="{legend_y:.1f}" '
            f'x2="{legend_x + 38:.1f}" y2="{legend_y:.1f}" '
            f'stroke="{color}" stroke-width="2.8"{dash_attr}/>'
        )
        parts.append(svg_text(legend_x + 50, legend_y + 4, label, 13, "start", "700"))
        legend_x += 132

    for panel_idx, order in enumerate(orders):
        panel_left = margin_left + panel_idx * (panel_width + panel_gap)
        panel_right = panel_left + panel_width
        panel_top = margin_top
        panel_bottom = panel_top + plot_height
        title = "s = 3" if order == 3 else "s = 4"

        parts.append(
            f'<rect x="{panel_left:.1f}" y="{panel_top:.1f}" '
            f'width="{panel_width:.1f}" height="{plot_height:.1f}" '
            'fill="#fbfcfe" stroke="#d6dde6"/>'
        )
        parts.append(svg_text(panel_left + 8, panel_top - 15, title, 14, "start", "700"))

        for tick in y_ticks:
            y = log_map(tick, y_min, y_max, panel_bottom, panel_top)
            parts.append(
                f'<line class="grid" x1="{panel_left:.1f}" y1="{y:.1f}" '
                f'x2="{panel_right:.1f}" y2="{y:.1f}"/>'
            )
            if panel_idx == 0:
                parts.append(svg_text(panel_left - 9, y + 4, format_tick(tick), 11, "end"))

        for tick in x_ticks:
            x = log_map(tick, x_min, x_max, panel_left, panel_right)
            parts.append(
                f'<line class="grid" x1="{x:.1f}" y1="{panel_top:.1f}" '
                f'x2="{x:.1f}" y2="{panel_bottom:.1f}"/>'
            )
            parts.append(svg_text(x, panel_bottom + 20, format_tick(tick), 11))

        parts.append(
            f'<line class="axis" x1="{panel_left:.1f}" y1="{panel_bottom:.1f}" '
            f'x2="{panel_right:.1f}" y2="{panel_bottom:.1f}"/>'
        )
        parts.append(
            f'<line class="axis" x1="{panel_left:.1f}" y1="{panel_top:.1f}" '
            f'x2="{panel_left:.1f}" y2="{panel_bottom:.1f}"/>'
        )

        panel_rows = sorted(grouped[order], key=lambda r: r["piece_count"])
        for _, key, color, dash in SERIES:
            points = []
            for row in panel_rows:
                x = log_map(row["piece_count"], x_min, x_max, panel_left, panel_right)
                y = log_map(row[key], y_min, y_max, panel_bottom, panel_top)
                points.append((x, y))
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            parts.append(
                f'<polyline class="series" points="'
                + " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
                + f'" stroke="{color}"{dash_attr}/>'
            )
            for x, y in points:
                parts.append(
                    f'<circle class="marker" cx="{x:.1f}" cy="{y:.1f}" '
                    f'r="3.6" fill="{color}"/>'
                )

    parts.append(svg_text(width / 2, height - 18, "M", 13, weight="700"))
    parts.append(svg_text(20, height / 2, "construction time (us)", 13, rotate=-90, weight="700"))
    parts.append("</svg>")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")


def main():
    if len(sys.argv) != 3:
        usage()
        return 2
    render_svg(read_rows(sys.argv[1]), sys.argv[2])
    print(f"wrote SVG: {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
