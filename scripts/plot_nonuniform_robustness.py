#!/usr/bin/env python3
"""Render reviewer-facing non-uniform-duration robustness figures from CSV."""

import csv
import html
import math
import os
import sys
from collections import defaultdict


ERROR_SERIES = [
    ("position vs MINCO", "max_position_error_vs_minco", "#0072B2", None),
    ("relative energy vs MINCO", "energy_rel_error_vs_minco", "#D55E00", None),
    ("boundary residual", "max_boundary_residual", "#009E73", None),
    ("waypoint residual", "max_waypoint_residual", "#CC79A7", None),
    ("waypoint derivative vs MINCO", "max_waypoint_derivative_error_vs_minco", "#E69F00", "7 4"),
    ("continuity jump", "max_continuity_jump", "#333333", "2 3"),
]

CONDITION_SERIES = [
    ("full-system cond. number", "full_system_condition_number_2", "#D55E00", None),
    ("reduced uniform-system cond. number", "reduced_uniform_system_condition_number_2", "#0072B2", None),
    ("relative solve residual", "max_relative_linear_residual", "#333333", "7 4"),
    ("minimum pivot", "min_abs_pivot", "#009E73", "2 3"),
]


def usage():
    print("usage: plot_nonuniform_robustness.py <robustness.csv> <output-dir>", file=sys.stderr)


def read_rows(path):
    rows = []
    numeric_keys = {
        "duration_ratio",
        "max_position_error_vs_minco",
        "energy_rel_error_vs_minco",
        "max_boundary_residual",
        "max_waypoint_residual",
        "max_waypoint_derivative_error_vs_minco",
        "max_continuity_jump",
        "full_system_condition_number_2",
        "reduced_uniform_system_condition_number_2",
        "min_abs_pivot",
        "max_relative_linear_residual",
    }
    with open(path, newline="") as csv_file:
        for raw in csv.DictReader(csv_file):
            row = dict(raw)
            row["order"] = int(row["order"])
            for key in numeric_keys:
                row[key] = float(row[key])
            rows.append(row)
    if not rows:
        raise RuntimeError("CSV has no rows")
    return rows


def positive(value):
    return max(value, 1.0e-20)


def bounds(rows, series):
    values = [positive(row[key]) for row in rows for _, key, _, _ in series]
    lo = 10 ** math.floor(math.log10(min(values)) - 0.15)
    hi = 10 ** math.ceil(math.log10(max(values)) + 0.15)
    if lo == hi:
        lo /= 10.0
        hi *= 10.0
    return lo, hi


def log_map(value, lo, hi, out_lo, out_hi):
    numerator = math.log10(positive(value)) - math.log10(lo)
    denominator = math.log10(hi) - math.log10(lo)
    return out_lo + numerator / denominator * (out_hi - out_lo)


def x_map(value, left, right):
    return log_map(value, 1.0, 1000.0, left, right)


def tick_label(value):
    exponent = int(round(math.log10(value)))
    return f"1e{exponent}" if exponent else "1"


def svg_text(x, y, text, size=12, anchor="middle", weight="400", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}"{transform}>'
        f"{html.escape(text)}</text>"
    )


def y_ticks(lo, hi):
    start = int(math.ceil(math.log10(lo)))
    stop = int(math.floor(math.log10(hi)))
    return [10 ** exponent for exponent in range(start, stop + 1)]


def render(rows, series, output_path, title, y_label, note):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["order"]].append(row)
    orders = sorted(grouped)
    lo, hi = bounds(rows, series)

    width = 1260
    height = 510
    margin_left = 86
    margin_right = 26
    margin_top = 108
    margin_bottom = 90
    gap = 48
    panel_width = (width - margin_left - margin_right - gap * (len(orders) - 1)) / len(orders)
    panel_height = height - margin_top - margin_bottom
    bottom = margin_top + panel_height

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; fill: #17202a; }",
        ".axis { stroke: #17202a; stroke-width: 1.15; }",
        ".grid { stroke: #d6dde6; stroke-width: 0.8; }",
        ".series { fill: none; stroke-width: 2.35; stroke-linejoin: round; }",
        ".marker { stroke: white; stroke-width: 1.0; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        svg_text(width / 2, 27, title, 17, weight="700"),
    ]

    legend_x = 54
    legend_y = 54
    for label, _, color, dash in series:
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(
            f'<line x1="{legend_x:.1f}" y1="{legend_y:.1f}" x2="{legend_x + 28:.1f}" y2="{legend_y:.1f}" '
            f'stroke="{color}" stroke-width="2.6"{dash_attr}/>'
        )
        parts.append(svg_text(legend_x + 34, legend_y + 4, label, 11, "start", "600"))
        legend_x += 192

    for panel_index, order in enumerate(orders):
        left = margin_left + panel_index * (panel_width + gap)
        right = left + panel_width
        panel_rows = sorted(grouped[order], key=lambda row: row["duration_ratio"])
        parts.append(
            f'<rect x="{left:.1f}" y="{margin_top:.1f}" width="{panel_width:.1f}" height="{panel_height:.1f}" '
            'fill="#fbfcfe" stroke="#d6dde6"/>'
        )
        parts.append(svg_text(left + 7, margin_top - 14, f"minimum {'acceleration' if order == 2 else 'jerk' if order == 3 else 'snap'} (s = {order})", 13, "start", "700"))

        for tick in y_ticks(lo, hi):
            y = log_map(tick, lo, hi, bottom, margin_top)
            parts.append(f'<line class="grid" x1="{left:.1f}" y1="{y:.1f}" x2="{right:.1f}" y2="{y:.1f}"/>')
            if panel_index == 0:
                parts.append(svg_text(left - 9, y + 4, tick_label(tick), 10, "end"))

        for ratio in (1, 10, 100, 1000):
            x = x_map(ratio, left, right)
            parts.append(f'<line class="grid" x1="{x:.1f}" y1="{margin_top:.1f}" x2="{x:.1f}" y2="{bottom:.1f}"/>')
            parts.append(svg_text(x, bottom + 20, str(ratio), 11))

        parts.append(f'<line class="axis" x1="{left:.1f}" y1="{bottom:.1f}" x2="{right:.1f}" y2="{bottom:.1f}"/>')
        parts.append(f'<line class="axis" x1="{left:.1f}" y1="{margin_top:.1f}" x2="{left:.1f}" y2="{bottom:.1f}"/>')

        for _, key, color, dash in series:
            points = [
                (x_map(row["duration_ratio"], left, right),
                 log_map(row[key], lo, hi, bottom, margin_top))
                for row in panel_rows
            ]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            parts.append(
                '<polyline class="series" points="'
                + " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
                + f'" stroke="{color}"{dash_attr}/>'
            )
            for x, y in points:
                parts.append(f'<circle class="marker" cx="{x:.1f}" cy="{y:.1f}" r="3.4" fill="{color}"/>')

    parts.append(svg_text(width / 2, height - 42, "maximum-to-minimum duration ratio R", 13, weight="700"))
    parts.append(svg_text(23, margin_top + panel_height / 2, y_label, 13, weight="700", rotate=-90))
    parts.append(svg_text(width / 2, height - 13, note, 10, weight="400"))
    parts.append("</svg>")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as out:
        out.write("\n".join(parts) + "\n")


def main():
    if len(sys.argv) != 3:
        usage()
        return 2
    rows = read_rows(sys.argv[1])
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    errors_path = os.path.join(output_dir, "nonuniform_robustness_errors.svg")
    conditioning_path = os.path.join(output_dir, "nonuniform_robustness_conditioning.svg")
    render(rows, ERROR_SERIES, errors_path,
           "Accuracy under strongly non-uniform durations",
           "maximum absolute / relative error (log scale)",
           "Fixed total duration; maxima over deterministic cases; lower is better.")
    render(rows, CONDITION_SERIES, conditioning_path,
           "Conditioning and linear-solve diagnostics",
           "condition number / diagnostic magnitude (log scale)",
           "Reduced uniform-system conditioning is independent of R for a fixed piece count.")
    print(f"wrote SVG: {errors_path}")
    print(f"wrote SVG: {conditioning_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
