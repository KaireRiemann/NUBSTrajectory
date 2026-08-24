#!/usr/bin/env python3
"""Render a paper-style optimization-variable comparison figure as SVG."""

import argparse
import html
import os


DEFAULT_OUTPUT = "docs/images/optimization_variable_comparison.svg"


COLORS = {
    "ink": "#17202a",
    "muted": "#52616f",
    "grid": "#d6dde6",
    "panel": "#fbfcfe",
    "paper": "#ffffff",
    "orange": "#d55e00",
    "orange_light": "#fff4eb",
    "green": "#009e73",
    "green_light": "#eaf7f2",
    "blue": "#0072b2",
    "blue_light": "#eef6fb",
    "gray": "#697586",
    "gray_light": "#f4f6f8",
    "warn": "#b23b3b",
}


def esc(text):
    return html.escape(text, quote=True)


def text(
    x,
    y,
    value,
    size=16,
    anchor="middle",
    weight=400,
    fill=None,
    family=None,
    style="",
    letter_spacing=0,
):
    fill = fill or COLORS["ink"]
    family_attr = f' font-family="{family}"' if family else ""
    style_attr = f' style="{style}"' if style else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
        f'letter-spacing="{letter_spacing}"{family_attr}{style_attr}>'
        f"{esc(value)}</text>"
    )


def multiline(
    x,
    y,
    lines,
    size=16,
    anchor="middle",
    weight=400,
    fill=None,
    line_height=22,
    family=None,
):
    fill = fill or COLORS["ink"]
    family_attr = f' font-family="{family}"' if family else ""
    spans = []
    for idx, line in enumerate(lines):
        dy = 0 if idx == 0 else line_height
        spans.append(f'<tspan x="{x:.1f}" dy="{dy}">{esc(line)}</tspan>')
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}"'
        f"{family_attr}>{''.join(spans)}</text>"
    )


def rounded_rect(x, y, w, h, fill, stroke, radius=8, sw=1.4, dash=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'rx="{radius:.1f}" ry="{radius:.1f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"{dash_attr}/>'
    )


def line(x1, y1, x2, y2, color, sw=1.5, dash=None, marker=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    marker_attr = f' marker-end="url(#{marker})"' if marker else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{color}" stroke-width="{sw}" stroke-linecap="round"'
        f"{dash_attr}{marker_attr}/>"
    )


def path(d, color, sw=2.0, fill="none", dash=None, marker=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    marker_attr = f' marker-end="url(#{marker})"' if marker else ""
    return (
        f'<path d="{d}" fill="{fill}" stroke="{color}" stroke-width="{sw}" '
        f'stroke-linecap="round" stroke-linejoin="round"{dash_attr}{marker_attr}/>'
    )


def circle(cx, cy, r, fill, stroke="#ffffff", sw=1.2):
    return (
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"/>'
    )


def box_label(parts, x, y, w, h, title, body_lines, accent, fill, subtitle=None):
    parts.append(rounded_rect(x, y, w, h, fill, accent, 8, 1.6))
    parts.append(text(x + 18, y + 25, title.upper(), 11, "start", 700, accent))
    if subtitle:
        parts.append(text(x + w - 18, y + 25, subtitle, 11, "end", 600, accent))
    parts.append(multiline(x + w / 2, y + 54, body_lines, 20, "middle", 700, COLORS["ink"], 25))


def pill(parts, x, y, w, h, label, fill, stroke, text_color=None):
    parts.append(rounded_rect(x, y, w, h, fill, stroke, h / 2, 1.3))
    parts.append(text(x + w / 2, y + h / 2 + 5, label, 13, "middle", 700, text_color or stroke))


def arrow_down(parts, x, y1, y2, color=COLORS["gray"]):
    marker = "arrow-green" if color == COLORS["green"] else "arrow-gray"
    parts.append(line(x, y1, x, y2, color, 1.7, marker=marker))


def draw_traditional(parts, x, y, w, h):
    heading_y = y + 36
    parts.append(text(x + w / 2, heading_y, "Traditional B-spline planner", 24, "middle", 700))
    parts.append(text(x + w / 2, heading_y + 26, "control points -> trajectory", 14, "middle", 500, COLORS["muted"]))

    var_y = y + 79
    box_label(
        parts,
        x + 42,
        var_y,
        w - 84,
        86,
        "optimized variables",
        ["B-spline control points", "C"],
        COLORS["orange"],
        COLORS["orange_light"],
    )

    arrow_down(parts, x + w / 2, var_y + 88, var_y + 124)

    scene_x = x + 42
    scene_y = y + 213
    scene_w = w - 84
    scene_h = 176
    parts.append(rounded_rect(scene_x, scene_y, scene_w, scene_h, "#ffffff", COLORS["grid"], 8, 1.2))
    parts.append(text(scene_x + 18, scene_y + 27, "control polygon is optimized directly", 12, "start", 600, COLORS["muted"]))

    pts = [
        (scene_x + 40, scene_y + 128),
        (scene_x + 96, scene_y + 64),
        (scene_x + 160, scene_y + 102),
        (scene_x + 224, scene_y + 42),
        (scene_x + 286, scene_y + 98),
        (scene_x + 354, scene_y + 59),
        (scene_x + 420, scene_y + 118),
    ]
    poly = " ".join(f"{px:.1f},{py:.1f}" for px, py in pts)
    parts.append(
        f'<polyline points="{poly}" fill="none" stroke="{COLORS["orange"]}" '
        f'stroke-width="1.9" stroke-dasharray="7 5" stroke-linejoin="round"/>'
    )
    curve_d = (
        f"M {pts[0][0]:.1f} {pts[0][1]:.1f} "
        f"C {scene_x + 86:.1f} {scene_y + 95:.1f}, {scene_x + 125:.1f} {scene_y + 80:.1f}, "
        f"{scene_x + 168:.1f} {scene_y + 91:.1f} "
        f"S {scene_x + 253:.1f} {scene_y + 102:.1f}, {scene_x + 298:.1f} {scene_y + 82:.1f} "
        f"S {scene_x + 374:.1f} {scene_y + 77:.1f}, {scene_x + 424:.1f} {scene_y + 111:.1f}"
    )
    parts.append(path(curve_d, COLORS["blue"], 3.0))
    for px, py in pts:
        parts.append(circle(px, py, 6.0, COLORS["orange"]))
    parts.append(text(scene_x + 61, scene_y + 151, "C_i", 13, "middle", 700, COLORS["orange"]))
    parts.append(text(scene_x + scene_w - 25, scene_y + 111, "p(t)", 14, "end", 700, COLORS["blue"]))

    arrow_down(parts, x + w / 2, scene_y + scene_h + 12, scene_y + scene_h + 48)
    out_y = scene_y + scene_h + 57
    parts.append(rounded_rect(x + 42, out_y, w - 84, 74, COLORS["blue_light"], COLORS["blue"], 8, 1.5))
    parts.append(multiline(x + w / 2, out_y + 30, ["B-spline trajectory", "from free control points"], 17, "middle", 700, COLORS["ink"], 23))


def draw_minco(parts, x, y, w, h):
    heading_y = y + 36
    parts.append(text(x + w / 2, heading_y, "MINCO", 24, "middle", 700))
    parts.append(text(x + w / 2, heading_y + 26, "(q, T) -> polynomial trajectory", 14, "middle", 500, COLORS["muted"]))

    var_y = y + 79
    box_label(
        parts,
        x + 42,
        var_y,
        w - 84,
        86,
        "optimized variables",
        ["waypoints q", "segment times T"],
        COLORS["green"],
        COLORS["green_light"],
    )

    arrow_down(parts, x + w / 2, var_y + 88, var_y + 124, COLORS["green"])

    construct_y = y + 213
    parts.append(rounded_rect(x + 42, construct_y, w - 84, 104, "#ffffff", COLORS["grid"], 8, 1.2))
    parts.append(text(x + w / 2, construct_y + 32, "piecewise polynomial construction", 17, "middle", 700))
    parts.append(text(x + w / 2, construct_y + 60, "minimum-control-effort coefficients", 13, "middle", 500, COLORS["muted"]))
    pill(parts, x + 79, construct_y + 72, 94, 30, "q, T", COLORS["green_light"], COLORS["green"])
    parts.append(line(x + 184, construct_y + 87, x + 290, construct_y + 87, COLORS["gray"], 1.6, marker="arrow-gray"))
    pill(parts, x + 304, construct_y + 72, 142, 30, "polynomial", COLORS["gray_light"], COLORS["gray"])

    arrow_down(parts, x + w / 2, construct_y + 116, construct_y + 153, COLORS["gray"])

    scene_y = construct_y + 162
    scene_x = x + 42
    scene_w = w - 84
    scene_h = 88
    parts.append(rounded_rect(scene_x, scene_y, scene_w, scene_h, COLORS["blue_light"], COLORS["blue"], 8, 1.5))
    curve_d = (
        f"M {scene_x + 42:.1f} {scene_y + 58:.1f} "
        f"C {scene_x + 95:.1f} {scene_y + 18:.1f}, {scene_x + 132:.1f} {scene_y + 25:.1f}, {scene_x + 179:.1f} {scene_y + 47:.1f} "
        f"C {scene_x + 224:.1f} {scene_y + 68:.1f}, {scene_x + 262:.1f} {scene_y + 17:.1f}, {scene_x + 313:.1f} {scene_y + 39:.1f} "
        f"C {scene_x + 355:.1f} {scene_y + 57:.1f}, {scene_x + 381:.1f} {scene_y + 59:.1f}, {scene_x + 426:.1f} {scene_y + 27:.1f}"
    )
    parts.append(path(curve_d, COLORS["blue"], 3.0))
    for px, py in (
        (scene_x + 42, scene_y + 58),
        (scene_x + 179, scene_y + 47),
        (scene_x + 313, scene_y + 39),
        (scene_x + 426, scene_y + 27),
    ):
        parts.append(circle(px, py, 4.8, COLORS["green"]))
    parts.append(multiline(scene_x + scene_w / 2, scene_y + 106, ["polynomial minimum-control-effort", "trajectory"], 17, "middle", 700, COLORS["ink"], 23))


def draw_nubs(parts, x, y, w, h):
    heading_y = y + 36
    parts.append(text(x + w / 2, heading_y, "NUBSTrajectory", 24, "middle", 700))
    parts.append(text(x + w / 2, heading_y + 26, "(q, T) -> A(T) C = b(q) -> NUBS", 14, "middle", 500, COLORS["muted"]))

    var_y = y + 79
    box_label(
        parts,
        x + 42,
        var_y,
        w - 84,
        86,
        "optimized variables",
        ["waypoints q", "segment times T"],
        COLORS["green"],
        COLORS["green_light"],
        subtitle="same as MINCO",
    )

    arrow_down(parts, x + w / 2, var_y + 88, var_y + 124, COLORS["green"])

    solve_y = y + 213
    parts.append(rounded_rect(x + 42, solve_y, w - 84, 112, "#ffffff", COLORS["grid"], 8, 1.2))
    parts.append(text(x + w / 2, solve_y + 37, "A(T) C = b(q)", 25, "middle", 700, COLORS["ink"], family="Cambria Math, STIX Two Math, Times New Roman, serif"))
    parts.append(text(x + w / 2, solve_y + 66, "structured NUBS construction system", 13, "middle", 500, COLORS["muted"]))
    parts.append(rounded_rect(x + 146, solve_y + 80, w - 292, 28, COLORS["gray_light"], COLORS["gray"], 14, 1.1, dash="5 4"))
    parts.append(text(x + w / 2, solve_y + 99, "C is recovered, not optimized", 12, "middle", 700, COLORS["gray"]))

    arrow_down(parts, x + w / 2, solve_y + 124, solve_y + 160, COLORS["gray"])

    scene_y = solve_y + 169
    scene_x = x + 42
    scene_w = w - 84
    scene_h = 102
    parts.append(rounded_rect(scene_x, scene_y, scene_w, scene_h, COLORS["blue_light"], COLORS["blue"], 8, 1.5))

    net_pts = [
        (scene_x + 40, scene_y + 54),
        (scene_x + 98, scene_y + 20),
        (scene_x + 152, scene_y + 44),
        (scene_x + 218, scene_y + 17),
        (scene_x + 287, scene_y + 36),
        (scene_x + 352, scene_y + 24),
        (scene_x + 426, scene_y + 48),
    ]
    poly = " ".join(f"{px:.1f},{py:.1f}" for px, py in net_pts)
    parts.append(
        f'<polyline points="{poly}" fill="none" stroke="{COLORS["gray"]}" '
        f'stroke-width="1.5" stroke-dasharray="5 5" stroke-linejoin="round" opacity="0.75"/>'
    )
    for px, py in net_pts:
        parts.append(circle(px, py, 3.7, COLORS["gray"], "#ffffff", 1.0))

    curve_d = (
        f"M {scene_x + 42:.1f} {scene_y + 53:.1f} "
        f"C {scene_x + 94:.1f} {scene_y + 29:.1f}, {scene_x + 138:.1f} {scene_y + 32:.1f}, {scene_x + 183:.1f} {scene_y + 40:.1f} "
        f"S {scene_x + 268:.1f} {scene_y + 34:.1f}, {scene_x + 315:.1f} {scene_y + 34:.1f} "
        f"S {scene_x + 378:.1f} {scene_y + 34:.1f}, {scene_x + 427:.1f} {scene_y + 46:.1f}"
    )
    parts.append(path(curve_d, COLORS["blue"], 3.0))
    parts.append(text(scene_x + 66, scene_y + 16, "derived C", 12, "middle", 700, COLORS["gray"]))
    parts.append(multiline(scene_x + scene_w / 2, scene_y + 77, ["NUBS minimum-control-effort", "trajectory"], 16, "middle", 700, COLORS["ink"], 21))

    note_y = y + h - 46
    parts.append(rounded_rect(x + 58, note_y, w - 116, 36, "#fffafa", COLORS["warn"], 8, 1.2))
    parts.append(text(x + w / 2, note_y + 24, "optimizer variables are q and T, not C", 13, "middle", 700, COLORS["warn"]))


def draw_figure(output_path):
    width = 1320
    height = 470
    margin_x = 54
    panel_y = 70
    panel_h = 318
    gap = 34
    col_w = (width - 2 * margin_x - 2 * gap) / 3

    x1 = margin_x
    x2 = x1 + col_w + gap
    x3 = x2 + col_w + gap

    parts = [
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
        "</defs>",
        "<style>",
        "text { font-family: Inter, Helvetica Neue, Arial, sans-serif; }",
        ".panel-title { font-weight: 700; }",
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{COLORS["paper"]}"/>',
    ]

    legend_x = 350
    legend_y = 28
    parts.append(rounded_rect(legend_x, legend_y - 18, 620, 34, "#ffffff", COLORS["grid"], 4, 0.9))
    legend_entries = [
        ("optimized", COLORS["green"], None),
        ("internal", COLORS["gray"], None),
        ("control points", COLORS["orange"], "7 5"),
        ("trajectory", COLORS["blue"], None),
    ]
    x_cursor = legend_x + 28
    for label, color, dash in legend_entries:
        parts.append(line(x_cursor, legend_y, x_cursor + 46, legend_y, color, 2.5, dash))
        parts.append(text(x_cursor + 58, legend_y + 4, label, 12, "start", 700, color))
        x_cursor += 150

    bracket_y = 58
    bracket_x1 = x2 + 45
    bracket_x2 = x3 + col_w - 45
    parts.append(line(bracket_x1, bracket_y, bracket_x2, bracket_y, COLORS["green"], 2.0))
    parts.append(line(bracket_x1, bracket_y, bracket_x1, bracket_y + 10, COLORS["green"], 2.0))
    parts.append(line(bracket_x2, bracket_y, bracket_x2, bracket_y + 10, COLORS["green"], 2.0))
    parts.append(text((bracket_x1 + bracket_x2) / 2, bracket_y - 7, "same opt. vars", 12, "middle", 700, COLORS["green"]))

    for x in (x1, x2, x3):
        parts.append(rounded_rect(x, panel_y, col_w, panel_h, COLORS["panel"], COLORS["grid"], 6, 1.0))

    def mini_curve(px, py, color=COLORS["blue"], with_polygon=False):
        pts = [
            (px + 24, py + 72),
            (px + 72, py + 30),
            (px + 124, py + 53),
            (px + 184, py + 25),
            (px + 250, py + 62),
            (px + 314, py + 38),
        ]
        if with_polygon:
            poly = " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)
            parts.append(
                f'<polyline points="{poly}" fill="none" stroke="{COLORS["orange"]}" '
                f'stroke-width="1.4" stroke-dasharray="6 5" stroke-linejoin="round"/>'
            )
            for a, b in pts:
                parts.append(circle(a, b, 4.5, COLORS["orange"]))
        d = (
            f"M {pts[0][0]:.1f} {pts[0][1]:.1f} "
            f"C {px + 70:.1f} {py + 45:.1f}, {px + 113:.1f} {py + 42:.1f}, {px + 153:.1f} {py + 50:.1f} "
            f"S {px + 230:.1f} {py + 54:.1f}, {px + 315:.1f} {py + 45:.1f}"
        )
        parts.append(path(d, color, 2.7))

    panels = [
        (x1, "(a) B-spline", "C", COLORS["orange"], "C -> p(t)", True),
        (x2, "(b) MINCO", "(q,T)", COLORS["green"], "coeff.", False),
        (x3, "(c) NUBS", "(q,T)", COLORS["green"], "A(T)C=b(q)", False),
    ]
    for x, title, var_label, var_color, mid_label, polygon in panels:
        parts.append(text(x + 18, panel_y + 28, title, 15, "start", 800))
        var_y = panel_y + 56
        parts.append(rounded_rect(x + 54, var_y, col_w - 108, 54,
                                  COLORS["orange_light"] if var_color == COLORS["orange"] else COLORS["green_light"],
                                  var_color, 6, 1.2))
        parts.append(text(x + col_w / 2, var_y + 34, var_label, 21, "middle", 800, var_color,
                          family="Cambria Math, STIX Two Math, Times New Roman, serif"))
        arrow_down(parts, x + col_w / 2, var_y + 58, var_y + 88,
                   COLORS["green"] if var_color == COLORS["green"] else COLORS["gray"])

        mid_y = var_y + 94
        stroke = COLORS["gray"] if x != x1 else COLORS["orange"]
        fill = COLORS["gray_light"] if x != x1 else COLORS["orange_light"]
        parts.append(rounded_rect(x + 54, mid_y, col_w - 108, 58, fill, stroke, 6, 1.1,
                                  dash="5 4" if x == x3 else None))
        parts.append(text(x + col_w / 2, mid_y + 36, mid_label, 18, "middle", 800, stroke,
                          family="Cambria Math, STIX Two Math, Times New Roman, serif"))
        if x == x3:
            parts.append(text(x + col_w - 72, mid_y + 50, "C", 12, "middle", 800, COLORS["gray"],
                              family="Cambria Math, STIX Two Math, Times New Roman, serif"))
        arrow_down(parts, x + col_w / 2, mid_y + 62, mid_y + 92, COLORS["gray"])
        mini_curve(x + 32, mid_y + 94, with_polygon=polygon)
        parts.append(text(x + col_w / 2, panel_y + panel_h - 20, "p(t)", 16, "middle", 800, COLORS["blue"],
                          family="Cambria Math, STIX Two Math, Times New Roman, serif"))

    parts.append("</svg>")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Draw the NUBSTrajectory optimization-variable comparison figure."
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=DEFAULT_OUTPUT,
        help=f"SVG output path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()
    draw_figure(args.output)
    print(f"wrote SVG: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
