#!/usr/bin/env python3
"""Generate paper-style SVG figures for NUBSTrajectory."""

import argparse
import csv
import html
import math
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "docs" / "data"
DEFAULT_IMAGE_DIR = ROOT / "docs" / "images"


COLORS = {
    "ink": "#17202a",
    "muted": "#52616f",
    "grid": "#d6dde6",
    "panel": "#fbfcfe",
    "paper": "#ffffff",
    "blue": "#0072b2",
    "blue_light": "#eef6fb",
    "orange": "#d55e00",
    "orange_light": "#fff4eb",
    "green": "#009e73",
    "green_light": "#eaf7f2",
    "purple": "#8a63d2",
    "purple_light": "#f3effc",
    "red": "#b23b3b",
    "red_light": "#fffafa",
    "gray": "#697586",
    "gray_light": "#f4f6f8",
    "dark_grid": "#aeb8c5",
}


def esc(value):
    return html.escape(str(value), quote=True)


def text(x, y, value, size=14, anchor="middle", weight=400, fill=None, family=None):
    fill = fill or COLORS["ink"]
    family_attr = f' font-family="{family}"' if family else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
        f'letter-spacing="0"{family_attr}>{esc(value)}</text>'
    )


def multiline(x, y, lines, size=14, anchor="middle", weight=400, fill=None, line_height=20):
    fill = fill or COLORS["ink"]
    spans = []
    for idx, line in enumerate(lines):
        dy = 0 if idx == 0 else line_height
        spans.append(f'<tspan x="{x:.1f}" dy="{dy}">{esc(line)}</tspan>')
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
        f'letter-spacing="0">{"".join(spans)}</text>'
    )


def rect(x, y, w, h, fill, stroke, radius=8, sw=1.2, dash=None, opacity=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    opacity_attr = f' opacity="{opacity}"' if opacity is not None else ""
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'rx="{radius:.1f}" ry="{radius:.1f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"{dash_attr}{opacity_attr}/>'
    )


def line(x1, y1, x2, y2, color, sw=1.3, dash=None, marker=None, opacity=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    marker_attr = f' marker-end="url(#{marker})"' if marker else ""
    opacity_attr = f' opacity="{opacity}"' if opacity is not None else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{color}" stroke-width="{sw}" stroke-linecap="round"'
        f"{dash_attr}{marker_attr}{opacity_attr}/>"
    )


def path(d, color, sw=2.0, fill="none", dash=None, marker=None, opacity=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    marker_attr = f' marker-end="url(#{marker})"' if marker else ""
    opacity_attr = f' opacity="{opacity}"' if opacity is not None else ""
    return (
        f'<path d="{d}" fill="{fill}" stroke="{color}" stroke-width="{sw}" '
        f'stroke-linecap="round" stroke-linejoin="round"{dash_attr}'
        f"{marker_attr}{opacity_attr}/>"
    )


def polyline(points, color, sw=2.0, dash=None, opacity=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    opacity_attr = f' opacity="{opacity}"' if opacity is not None else ""
    point_attr = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return (
        f'<polyline points="{point_attr}" fill="none" stroke="{color}" '
        f'stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round"'
        f"{dash_attr}{opacity_attr}/>"
    )


def circle(cx, cy, r, fill, stroke="#ffffff", sw=1.0):
    return (
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"/>'
    )


def svg_open(width, height):
    return [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        "<defs>",
        '<marker id="arrow-gray" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">',
        f'<path d="M 0 0 L 10 4 L 0 8 z" fill="{COLORS["gray"]}"/>',
        "</marker>",
        '<marker id="arrow-green" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">',
        f'<path d="M 0 0 L 10 4 L 0 8 z" fill="{COLORS["green"]}"/>',
        "</marker>",
        '<marker id="arrow-blue" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">',
        f'<path d="M 0 0 L 10 4 L 0 8 z" fill="{COLORS["blue"]}"/>',
        "</marker>",
        "</defs>",
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{COLORS["paper"]}"/>',
    ]


def write_svg(parts, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(parts + ["</svg>"]) + "\n")


def read_csv(path):
    with path.open(newline="") as f:
        return [dict(row) for row in csv.DictReader(f)]


def fnum(row, key, default=0.0):
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def fmt_sci(value, digits=2):
    if value is None or value == "":
        return "-"
    value = float(value)
    if value == 0:
        return "0"
    return f"{value:.{digits}e}"


def nice_bounds(values, pad=0.08):
    lo = min(values)
    hi = max(values)
    if abs(hi - lo) < 1.0e-12:
        hi += 1.0
        lo -= 1.0
    delta = hi - lo
    return lo - pad * delta, hi + pad * delta


def ensure_equivalence_data(data_dir):
    needed = [
        data_dir / "fig2_equivalence_samples_s4.csv",
        data_dir / "fig2_equivalence_waypoints_s4.csv",
        data_dir / "fig2_equivalence_summary.csv",
    ]
    if all(path.exists() for path in needed):
        return
    exporter = ROOT / "bin" / "export_equivalence_figure_data"
    if exporter.exists():
        subprocess.run([str(exporter), str(data_dir)], cwd=ROOT, check=True)
        return
    missing = ", ".join(str(path) for path in needed if not path.exists())
    raise RuntimeError(
        "Missing Fig.2 data. Build and run export_equivalence_figure_data first: "
        f"{missing}"
    )


def draw_axes(parts, x, y, w, h, x_ticks, y_ticks, x_map, y_map, x_label, y_label):
    parts.append(rect(x, y, w, h, COLORS["panel"], COLORS["grid"], 4, 1.0))
    for tick in x_ticks:
        px = x_map(tick)
        parts.append(line(px, y, px, y + h, COLORS["grid"], 0.8))
        parts.append(text(px, y + h + 22, tick, 11, "middle", 400, COLORS["muted"]))
    for tick, label in y_ticks:
        py = y_map(tick)
        parts.append(line(x, py, x + w, py, COLORS["grid"], 0.8))
        parts.append(text(x - 9, py + 4, label, 11, "end", 400, COLORS["muted"]))
    parts.append(line(x, y + h, x + w, y + h, COLORS["ink"], 1.2))
    parts.append(line(x, y, x, y + h, COLORS["ink"], 1.2))
    parts.append(text(x + w / 2, y + h + 48, x_label, 13, "middle", 600, COLORS["ink"]))
    parts.append(text(x - 52, y + h / 2, y_label, 13, "middle", 600, COLORS["ink"]))


def render_fig2(data_dir, image_dir):
    ensure_equivalence_data(data_dir)
    samples = read_csv(data_dir / "fig2_equivalence_samples_s4.csv")
    waypoints = read_csv(data_dir / "fig2_equivalence_waypoints_s4.csv")
    summary = read_csv(data_dir / "fig2_equivalence_summary.csv")

    width, height = 1500, 700
    parts = svg_open(width, height)

    plot_x, plot_y, plot_w, plot_h = 78, 70, 610, 410
    xs = [fnum(row, "nubs_x") for row in samples] + [fnum(row, "minco_x") for row in samples]
    ys = [fnum(row, "nubs_y") for row in samples] + [fnum(row, "minco_y") for row in samples]
    x_lo, x_hi = nice_bounds(xs)
    y_lo, y_hi = nice_bounds(ys)
    ratio = (x_hi - x_lo) / (y_hi - y_lo)
    target_ratio = plot_w / plot_h
    if ratio > target_ratio:
        center = 0.5 * (y_lo + y_hi)
        half = 0.5 * (x_hi - x_lo) / target_ratio
        y_lo, y_hi = center - half, center + half
    else:
        center = 0.5 * (x_lo + x_hi)
        half = 0.5 * (y_hi - y_lo) * target_ratio
        x_lo, x_hi = center - half, center + half

    def mx(v):
        return plot_x + (v - x_lo) / (x_hi - x_lo) * plot_w

    def my(v):
        return plot_y + plot_h - (v - y_lo) / (y_hi - y_lo) * plot_h

    x_ticks = [round(x_lo + i * (x_hi - x_lo) / 4.0, 2) for i in range(5)]
    y_ticks = [
        (y_lo + i * (y_hi - y_lo) / 4.0, f"{y_lo + i * (y_hi - y_lo) / 4.0:.2f}")
        for i in range(5)
    ]
    draw_axes(parts, plot_x, plot_y, plot_w, plot_h, x_ticks, y_ticks, mx, my, "x", "y")
    parts.append(text(plot_x + 10, plot_y - 22, "(a)", 15, "start", 800))

    minco_points = [(mx(fnum(row, "minco_x")), my(fnum(row, "minco_y"))) for row in samples]
    nubs_points = [(mx(fnum(row, "nubs_x")), my(fnum(row, "nubs_y"))) for row in samples]
    waypoint_points = [(mx(fnum(row, "x")), my(fnum(row, "y"))) for row in waypoints]
    parts.append(polyline(minco_points, COLORS["orange"], 4.0, "10 7", 0.92))
    parts.append(polyline(nubs_points, COLORS["blue"], 2.4))
    for idx, (px, py) in enumerate(waypoint_points):
        parts.append(circle(px, py, 5.2, COLORS["green"], "#ffffff", 1.4))
        if idx in (0, len(waypoint_points) - 1):
            parts.append(text(px, py - 11, f"q{idx}", 11, "middle", 700, COLORS["green"]))

    legend_y = plot_y + 22
    legend_x = plot_x + plot_w - 230
    parts.append(rect(legend_x, legend_y - 18, 196, 86, "#ffffff", COLORS["grid"], 4, 0.9))
    parts.append(line(legend_x + 18, legend_y, legend_x + 70, legend_y, COLORS["orange"], 3.5, "9 6"))
    parts.append(text(legend_x + 82, legend_y + 4, "MINCO", 12, "start", 600, COLORS["ink"]))
    parts.append(line(legend_x + 18, legend_y + 28, legend_x + 70, legend_y + 28, COLORS["blue"], 2.4))
    parts.append(text(legend_x + 82, legend_y + 32, "NUBS", 12, "start", 600, COLORS["ink"]))
    parts.append(circle(legend_x + 44, legend_y + 56, 5.0, COLORS["green"]))
    parts.append(text(legend_x + 82, legend_y + 60, "q", 12, "start", 600, COLORS["ink"]))

    err_x, err_y, err_w, err_h = 820, 70, 590, 410
    total_t = fnum(samples[0], "total_time")
    t_values = [fnum(row, "t") for row in samples]
    err_values = [fnum(row, "err_d0") for row in samples]
    err_max = max(err_values)
    log_lo = 1.0e-17
    log_hi = 1.0e-13

    def mt(v):
        return err_x + v / total_t * err_w

    def me(v):
        v = max(v, log_lo)
        return err_y + err_h - (math.log10(v) - math.log10(log_lo)) / (math.log10(log_hi) - math.log10(log_lo)) * err_h

    y_ticks = [(10.0 ** k, f"1e{k}") for k in range(-17, -12)]
    x_ticks = [round(total_t * i / 4.0, 2) for i in range(5)]
    draw_axes(parts, err_x, err_y, err_w, err_h, x_ticks, y_ticks, mt, me, "t", "error")
    parts.append(text(err_x + 10, err_y - 22, "(b)", 15, "start", 800))
    err_points = [(mt(t), me(err)) for t, err in zip(t_values, err_values)]
    parts.append(polyline(err_points, COLORS["red"], 2.3))
    parts.append(rect(err_x + 32, err_y + 28, 185, 52, "#ffffff", COLORS["grid"], 4, 0.9))
    parts.append(line(err_x + 50, err_y + 54, err_x + 102, err_y + 54, COLORS["red"], 2.3))
    parts.append(text(err_x + 114, err_y + 58, "||p_N-p_M||", 12, "start", 600, COLORS["ink"]))
    parts.append(text(err_x + err_w - 18, err_y + 28, f"max = {fmt_sci(err_max, 2)}", 13, "end", 800, COLORS["red"]))

    table_x, table_y, table_w, table_h = 190, 560, 1120, 92
    parts.append(text(table_x - 24, table_y + 24, "(c)", 15, "end", 800))
    parts.append(rect(table_x, table_y, table_w, table_h, "#ffffff", COLORS["grid"], 4, 1.0))
    columns = [("s", 70), ("p", 130), ("v", 130), ("a", 130), ("j", 130), ("snap", 130), ("|dE|", 150), ("rel(E)", 150)]
    start_x = table_x + 24
    row_y0 = table_y + 27
    for row_idx, s in enumerate((3, 4)):
        rows = [row for row in summary if int(row["s"]) == s]
        agg = {}
        for key in [
            "max_position",
            "max_velocity",
            "max_acceleration",
            "max_jerk",
            "max_snap",
            "energy_abs_error",
            "energy_rel_error",
        ]:
            values = [fnum(row, key, None) for row in rows if row.get(key, "") != ""]
            agg[key] = max(values) if values else None
        values = [
            str(s),
            fmt_sci(agg["max_position"]),
            fmt_sci(agg["max_velocity"]),
            fmt_sci(agg["max_acceleration"]),
            fmt_sci(agg["max_jerk"]),
            fmt_sci(agg["max_snap"]),
            fmt_sci(agg["energy_abs_error"]),
            fmt_sci(agg["energy_rel_error"]),
        ]
        if row_idx == 0:
            x_cursor = start_x
            for label, w in columns:
                parts.append(text(x_cursor + w / 2, row_y0, label, 12, "middle", 700, COLORS["muted"]))
                x_cursor += w
            parts.append(line(table_x + 18, row_y0 + 12, table_x + table_w - 18, row_y0 + 12, COLORS["grid"], 1.0))
        y = row_y0 + 40 + 26 * row_idx
        x_cursor = start_x
        for value, (_, w) in zip(values, columns):
            parts.append(text(x_cursor + w / 2, y, value, 13, "middle", 700 if value == str(s) else 500))
            x_cursor += w

    write_svg(parts, image_dir / "fig2_nubs_minco_equivalence.svg")


def plot_line_panel(parts, x, y, w, h, title, s, m_max=64):
    y_max = 2 * s * m_max * 1.08

    def mx(v):
        return x + (v - 1) / (m_max - 1) * w

    def my(v):
        return y + h - v / y_max * h

    x_ticks = [1, 16, 32, 48, 64]
    y_step = 100 if s == 3 else 125
    y_ticks = [(v, str(v)) for v in range(0, int(y_max) + 1, y_step)]
    draw_axes(parts, x, y, w, h, x_ticks, y_ticks, mx, my, "M", "dim.")
    parts.append(text(x + 8, y - 20, title, 15, "start", 800))
    values = list(range(1, m_max + 1))
    series = [
        ("2Ms  MINCO / polynomial coeff.", lambda m: 2 * m * s, COLORS["orange"], None),
        ("M+2s-1  NUBS", lambda m: m + 2 * s - 1, COLORS["blue"], None),
        ("M-1  UBS online solve", lambda m: max(0, m - 1), COLORS["green"], "8 5"),
    ]
    for label, fn, color, dash in series:
        pts = [(mx(m), my(fn(m))) for m in values]
        parts.append(polyline(pts, color, 2.8, dash))


def render_fig3(image_dir):
    width, height = 1320, 520
    parts = svg_open(width, height)
    plot_line_panel(parts, 82, 84, 540, 340, "(a) s=3", 3)
    plot_line_panel(parts, 724, 84, 540, 340, "(b) s=4", 4)
    legend_x, legend_y = 345, 24
    parts.append(rect(legend_x, legend_y, 630, 34, "#ffffff", COLORS["grid"], 4, 0.9))
    entries = [
        ("2Ms", COLORS["orange"], None),
        ("M+2s-1", COLORS["blue"], None),
        ("M-1", COLORS["green"], "8 5"),
    ]
    x_cursor = legend_x + 30
    for label, color, dash in entries:
        parts.append(line(x_cursor, legend_y + 18, x_cursor + 55, legend_y + 18, color, 2.6, dash))
        parts.append(text(x_cursor + 68, legend_y + 22, label, 13, "start", 700, color))
        x_cursor += 195
    write_svg(parts, image_dir / "fig3_system_dimension_comparison.svg")


def render_matrix_grid(parts, x, y, n, cell, nonzeros, group_rows=None):
    size = n * cell
    parts.append(rect(x, y, size, size, "#ffffff", COLORS["grid"], 4, 1.0))
    for i in range(n + 1):
        color = COLORS["dark_grid"] if i % 4 == 0 else COLORS["grid"]
        sw = 0.75 if i % 4 == 0 else 0.45
        parts.append(line(x + i * cell, y, x + i * cell, y + size, color, sw))
        parts.append(line(x, y + i * cell, x + size, y + i * cell, color, sw))
    if group_rows:
        for start, count, color in group_rows:
            parts.append(rect(x - 18, y + start * cell, 10, count * cell, color, color, 2, 0.0))
    for r, c, color in nonzeros:
        parts.append(rect(x + c * cell + 2, y + r * cell + 2, cell - 4, cell - 4, color, color, 2, 0.0))


def render_fig4(image_dir):
    width, height = 1120, 620
    m, s = 13, 4
    n = m + 2 * s - 1
    cell = 23
    grid_x, grid_y = 190, 76
    parts = svg_open(width, height)

    nonzeros = []
    for r in range(n):
        if r < s:
            start = 0
            end = min(n - 1, s + r)
            color = COLORS["purple"]
        elif r < s + m - 1:
            i = r - s
            start = i
            end = min(n - 1, i + 2 * s - 1)
            color = COLORS["blue"]
        else:
            k = r - (s + m - 1)
            start = max(0, n - (2 * s - k))
            end = n - 1
            color = COLORS["green"]
        for c in range(start, end + 1):
            nonzeros.append((r, c, color))

    group_rows = [
        (0, s, COLORS["purple"]),
        (s, m - 1, COLORS["blue"]),
        (s + m - 1, s, COLORS["green"]),
    ]
    render_matrix_grid(parts, grid_x, grid_y, n, cell, nonzeros, group_rows)
    size = n * cell
    parts.append(text(grid_x + size / 2, grid_y - 22, "C_0 ... C_{N_c-1}", 14, "middle", 700))
    parts.append(text(grid_x - 42, grid_y + size / 2, "rows", 13, "middle", 700))
    parts.append(text(grid_x + size + 18, grid_y + 2.1 * cell, "head", 12, "start", 700, COLORS["purple"]))
    parts.append(text(grid_x + size + 18, grid_y + (s + (m - 1) / 2) * cell, "q", 12, "start", 700, COLORS["blue"]))
    parts.append(text(grid_x + size + 18, grid_y + (n - 1.4) * cell, "tail", 12, "start", 700, COLORS["green"]))
    parts.append(path(f"M {grid_x + 18:.1f} {grid_y + 22:.1f} L {grid_x + size - 20:.1f} {grid_y + size - 98:.1f}", COLORS["red"], 2.0, dash="8 6", opacity=0.85))
    parts.append(text(grid_x + size - 8, grid_y + size - 106, "band", 12, "end", 800, COLORS["red"]))

    legend_x, legend_y = 750, 126
    parts.append(rect(legend_x, legend_y, 265, 184, "#ffffff", COLORS["grid"], 4, 1.0))
    legend_entries = [
        ("boundary", COLORS["purple"]),
        ("waypoint", COLORS["blue"]),
        ("terminal", COLORS["green"]),
        ("bandwidth <= 2s", COLORS["red"]),
        ("N_c=M+2s-1", COLORS["ink"]),
    ]
    for idx, (label, color) in enumerate(legend_entries):
        yy = legend_y + 28 + idx * 30
        if idx < 3:
            parts.append(rect(legend_x + 22, yy - 12, 18, 18, color, color, 2, 0.0))
        elif idx == 3:
            parts.append(line(legend_x + 22, yy - 3, legend_x + 54, yy - 3, color, 2.0, "8 5"))
        else:
            parts.append(text(legend_x + 38, yy + 2, "A(T)", 12, "middle", 800, color))
        parts.append(text(legend_x + 68, yy + 2, label, 13, "start", 700, color))
    write_svg(parts, image_dir / "fig4_nubs_matrix_sparsity.svg")


def draw_small_banded_matrix(parts, x, y, n, cell, color):
    nonzeros = []
    bandwidth = 5
    for r in range(n):
        for c in range(max(0, r - 2), min(n, r + bandwidth - 2)):
            nonzeros.append((r, c, color))
    render_matrix_grid(parts, x, y, n, cell, nonzeros)


def block(parts, x, y, w, h, title, body, fill, stroke):
    parts.append(rect(x, y, w, h, fill, stroke, 8, 1.2))
    parts.append(text(x + w / 2, y + 28, title, 14, "middle", 800, stroke))
    if isinstance(body, list):
        parts.append(multiline(x + w / 2, y + 59, body, 15, "middle", 700, COLORS["ink"], 22))
    else:
        parts.append(text(x + w / 2, y + 62, body, 17, "middle", 700, COLORS["ink"]))


def render_fig5(image_dir):
    width, height = 1320, 520
    parts = svg_open(width, height)
    block(parts, 70, 70, 170, 70, "head", "x_0^{(0:s-1)}", COLORS["purple_light"], COLORS["purple"])
    block(parts, 70, 225, 170, 70, "q", "q_1...q_{M-1}", COLORS["green_light"], COLORS["green"])
    block(parts, 70, 380, 170, 70, "tail", "x_f^{(0:s-1)}", COLORS["purple_light"], COLORS["purple"])

    block(parts, 355, 70, 190, 70, "direct", "C_H", COLORS["gray_light"], COLORS["gray"])
    block(parts, 355, 380, 190, 70, "direct", "C_T", COLORS["gray_light"], COLORS["gray"])
    parts.append(line(240, 105, 355, 105, COLORS["gray"], 1.8, marker="arrow-gray"))
    parts.append(line(240, 415, 355, 415, COLORS["gray"], 1.8, marker="arrow-gray"))

    solve_x, solve_y = 360, 205
    parts.append(rect(solve_x, solve_y, 390, 95, COLORS["blue_light"], COLORS["blue"], 6, 1.2))
    parts.append(text(solve_x + 195, solve_y + 39, "A^u_{M,s} C_I = r", 22, "middle", 800, COLORS["ink"], "Cambria Math, STIX Two Math, Times New Roman, serif"))
    parts.append(text(solve_x + 195, solve_y + 70, "dim(C_I)=M-1", 13, "middle", 700, COLORS["blue"]))
    parts.append(line(240, 260, solve_x, solve_y + 48, COLORS["green"], 1.8, marker="arrow-green"))
    parts.append(line(450, 140, solve_x + 70, solve_y, COLORS["gray"], 1.3, marker="arrow-gray"))
    parts.append(line(450, 380, solve_x + 70, solve_y + 95, COLORS["gray"], 1.3, marker="arrow-gray"))

    mat_x, mat_y = 900, 103
    draw_small_banded_matrix(parts, mat_x, mat_y, 12, 20, COLORS["blue"])
    parts.append(text(mat_x + 120, mat_y - 22, "A^u_{M,s}", 18, "middle", 800, COLORS["blue"], "Cambria Math, STIX Two Math, Times New Roman, serif"))
    parts.append(line(solve_x + 390, solve_y + 48, mat_x - 32, mat_y + 120, COLORS["blue"], 1.8, marker="arrow-blue"))

    legend_x, legend_y = 865, 380
    parts.append(rect(legend_x, legend_y, 320, 58, "#ffffff", COLORS["grid"], 4, 1.0))
    parts.append(text(legend_x + 26, legend_y + 25, "cache:", 13, "start", 700, COLORS["ink"]))
    parts.append(text(legend_x + 94, legend_y + 25, "(M,s)", 13, "start", 800, COLORS["blue"]))
    parts.append(text(legend_x + 26, legend_y + 47, "independent of q, T, boundary values", 12, "start", 600, COLORS["muted"]))

    x0, y0 = 350, 462
    parts.append(rect(x0, y0, 150, 28, COLORS["purple_light"], COLORS["purple"], 4, 1.0))
    parts.append(text(x0 + 75, y0 + 19, "s", 12, "middle", 800, COLORS["purple"]))
    parts.append(rect(x0 + 150, y0, 300, 28, COLORS["blue_light"], COLORS["blue"], 4, 1.0))
    parts.append(text(x0 + 300, y0 + 19, "M-1", 12, "middle", 800, COLORS["blue"]))
    parts.append(rect(x0 + 450, y0, 150, 28, COLORS["purple_light"], COLORS["purple"], 4, 1.0))
    parts.append(text(x0 + 525, y0 + 19, "s", 12, "middle", 800, COLORS["purple"]))

    write_svg(parts, image_dir / "fig5_ubs_reduced_constant_matrix.svg")


def main():
    parser = argparse.ArgumentParser(description="Generate NUBSTrajectory paper figures.")
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR), help="CSV data directory")
    parser.add_argument("--image-dir", default=str(DEFAULT_IMAGE_DIR), help="SVG output directory")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    image_dir = Path(args.image_dir)
    render_fig2(data_dir, image_dir)
    render_fig3(image_dir)
    render_fig4(image_dir)
    render_fig5(image_dir)
    print(f"wrote SVGs to {image_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
