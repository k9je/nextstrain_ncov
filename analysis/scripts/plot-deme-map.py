#!/usr/bin/env python3
"""
Draw a publication-quality phylogeographic map from a Nextstrain/Auspice JSON.

Unlike Auspice, which scales deme circles continuously by tip count (so one
heavily-sequenced deme swamps the map), this draws:

  * one PIE per deme, sized from a 3-level BIN of its tip count, with slices
    coloured by Nextclade major lineage;
  * transmission lines whose widths come from a 3-level BIN of the number of
    inferred deme transitions, with deliberately exaggerated contrast.

Transmission lines are counted the same way Auspice draws them: a branch whose
parent and child carry different reconstructed `deme` values is one event.

Because dengue demes cluster tightly in Southeast Asia, overlapping pies are
pushed apart by a short force-directed relaxation and tied back to their true
centroid with a leader line and a dot. Labels are then placed by testing
candidate positions against the real rendered text extents.
"""
import argparse
import json
import sys
import textwrap
from collections import Counter, defaultdict
from itertools import combinations

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Wedge

# --- binning ---------------------------------------------------------------
# Pie size: tip-count bins. Radii (degrees) are deliberately NOT
# area-proportional -- that is the whole point of binning.
SIZE_BINS = [(1, 9, 1.7, "1-9"), (10, 99, 3.0, "10-99"), (100, 10**9, 4.6, "100+")]

# Line width: transition-count bins, with exaggerated contrast.
LINE_BINS = [
    (1, 2, 0.6, 0.26, False, "1-2"),
    (3, 9, 2.4, 0.58, True, "3-9"),
    (10, 10**9, 5.6, 0.88, True, "10+"),
]

LINE_COLOR = "#1f1f1f"
OTHER_COLOR = "#9e9e9e"   # merged small pie slices
# Focus deme highlighting: imports red, exports blue.
IN_COLOR = "#d62728"
OUT_COLOR = "#1f77b4"
FOCUS_MIN_WIDTH = 1.5   # keep highlighted lines visible even in the 1-2 bin
FOCUS_MIN_ALPHA = 0.92

PIE_PAD = 2.6          # extra degrees kept between pies (room for labels)
MAX_SHIFT = 30.0       # max degrees a pie may be displaced from its centroid
PUSH_GAIN = 0.35       # how hard overlapping pies repel
SPRING_GAIN = 0.0030   # how hard a pie is pulled back to its true centroid
LEADER_MIN = 0.9       # draw a leader line past this displacement

# Map extent and layout (inches). The axes box is sized to the map's own
# aspect ratio so cartopy does not letterbox it inside a taller box. The
# latitude window is derived from the data so no pie is clipped at an edge.
LAT_PAD = 6.0
LAT_MIN_BOUNDS = (-80.0, -42.0)   # (hard floor, latest possible top of range)
LAT_MAX_BOUNDS = (52.0, 84.0)
FIG_WIDTH = 16.5
AX_WIDTH_FRAC = 0.984
LEGEND_BAND = 2.0
TITLE_BAND = 0.5


def bin_for(value, bins):
    for b in bins:
        if b[0] <= value <= b[1]:
            return b
    return bins[-1]


def walk_tree(root, deme_key, lineage_key):
    """Return (tips-per-deme-per-lineage, directed transition counts)."""
    per_deme = defaultdict(Counter)
    transitions = Counter()

    def recurse(node, parent_deme):
        attrs = node.get("node_attrs", {})
        deme = attrs.get(deme_key, {}).get("value")
        children = node.get("children")
        if parent_deme and deme and parent_deme != deme:
            transitions[(parent_deme, deme)] += 1
        if not children:
            if deme:
                lineage = attrs.get(lineage_key, {}).get("value") or "unassigned"
                per_deme[deme][lineage] += 1
        for child in children or []:
            recurse(child, deme or parent_deme)

    recurse(root, None)
    return per_deme, transitions


def lineage_colors(meta, lineage_key, lineages):
    """Reuse the dataset's own colour scale so the figure matches Auspice."""
    colors = {}
    for coloring in meta.get("colorings", []):
        if coloring["key"] == lineage_key:
            for entry in coloring.get("scale", []) or []:
                colors[entry[0]] = entry[1]
    fallback = plt.get_cmap("tab20")
    for i, lineage in enumerate(sorted(lineages)):
        colors.setdefault(lineage, matplotlib.colors.to_hex(fallback(i % 20)))
    colors["unassigned"] = "#bdbdbd"
    return colors


def deme_coordinates(meta, deme_key):
    for resolution in meta.get("geo_resolutions", []):
        if resolution["key"] == deme_key:
            return {
                name: (float(v["longitude"]), float(v["latitude"]))
                for name, v in resolution["demes"].items()
            }
    raise SystemExit(f"No geo_resolution {deme_key!r} in this dataset")


def relax_positions(true_pos, radii, pad=PIE_PAD, max_shift=MAX_SHIFT,
                    iterations=1500):
    """Push overlapping pies apart, with a weak spring back to the centroid."""
    pos = {k: np.array(v, dtype=float) for k, v in true_pos.items()}
    names = list(pos)
    for _ in range(iterations):
        moved = False
        for a, b in combinations(names, 2):
            delta = pos[b] - pos[a]
            dist = float(np.hypot(*delta))
            want = radii[a] + radii[b] + pad
            if dist < 1e-6:
                delta = np.array([0.5, 0.5])
                dist = float(np.hypot(*delta))
            if dist < want:
                push = (want - dist) * PUSH_GAIN
                unit = delta / dist
                pos[a] -= unit * push
                pos[b] += unit * push
                moved = True
        # weak spring home, and hard cap on displacement
        for name in names:
            home = np.array(true_pos[name], dtype=float)
            pos[name] += (home - pos[name]) * SPRING_GAIN
            offset = pos[name] - home
            dist = float(np.hypot(*offset))
            if dist > max_shift:
                pos[name] = home + offset / dist * max_shift
        if not moved:
            break
    return pos


def measure_legend(fig_width, spec, font, title_font):
    """Render a legend on a throwaway canvas to get its size in inches.

    The legend band's depth has to be known before the real figure is built,
    and at large --font-scale the legends no longer fit on one row, so we
    measure them first and pack them into rows afterwards.
    """
    probe = plt.figure(figsize=(fig_width, 6))
    legend = probe.legend(
        handles=spec["handles"], title=spec["title"], loc="upper left",
        bbox_to_anchor=(0.0, 1.0), frameon=True, fontsize=font,
        title_fontsize=title_font, **spec.get("kw", {}),
    )
    probe.canvas.draw()
    bbox = legend.get_window_extent(probe.canvas.get_renderer())
    size = (bbox.width / probe.dpi, bbox.height / probe.dpi)
    plt.close(probe)
    return size


def pack_rows(specs, sizes, max_width, gap=0.35):
    """Greedily pack legends into rows no wider than max_width."""
    rows, current, used = [], [], 0.0
    for spec, (width, height) in zip(specs, sizes):
        advance = width if not current else gap + width
        if current and used + advance > max_width:
            rows.append(current)
            current, used = [], 0.0
            advance = width
        current.append((spec, width, height))
        used += advance
    if current:
        rows.append(current)
    return rows


def wrapped_offsets(x1, x2):
    """Longitude offsets needed to draw a line that may cross the map edge."""
    if x2 - x1 > 180:
        return [(0.0, -360.0), (360.0, 0.0)]
    if x2 - x1 < -180:
        return [(0.0, 360.0), (-360.0, 0.0)]
    return [(0.0, 0.0)]


def place_labels(ax, fig, entries, pie_display, fontsize):
    """Greedily place labels at the first candidate offset that does not collide.

    `entries` is a list of (name, text, data_xy, radius_deg) sorted so the
    most important labels are placed first. Collision testing uses the real
    rendered text extents in display coordinates.
    """
    renderer = fig.canvas.get_renderer()
    placed = []
    results = []

    # candidate directions, in order of preference
    directions = [
        (0, -1), (0, 1), (1, 0), (-1, 0),
        (0.85, -0.85), (-0.85, -0.85), (0.85, 0.85), (-0.85, 0.85),
    ]

    for name, text, (dx, dy), radius_deg in entries:
        cx, cy = ax.transData.transform((dx, dy))
        redge = abs(ax.transData.transform((dx + radius_deg, dy))[0] - cx)

        probe = ax.text(0, 0, text, fontsize=fontsize, ha="center", va="center")
        bbox = probe.get_window_extent(renderer)
        tw, th = bbox.width, bbox.height
        probe.remove()

        best, best_score = None, None
        for scale in (1.0, 1.5, 2.1, 2.9, 3.8):
            for ux, uy in directions:
                gap = redge + 8 + (th / 2 if uy else tw / 2) * 0.55
                px = cx + ux * gap * scale
                py = cy + uy * gap * scale
                rect = (px - tw / 2, py - th / 2, px + tw / 2, py + th / 2)
                score = 0.0
                for ox, oy, orad in pie_display:
                    nx = max(rect[0], min(ox, rect[2]))
                    ny = max(rect[1], min(oy, rect[3]))
                    if np.hypot(nx - ox, ny - oy) < orad:
                        score += orad
                for prect in placed:
                    ox1 = min(rect[2], prect[2]) - max(rect[0], prect[0])
                    oy1 = min(rect[3], prect[3]) - max(rect[1], prect[1])
                    if ox1 > 0 and oy1 > 0:
                        score += (ox1 * oy1) ** 0.5
                score += 0.35 * (scale - 1.0) * (tw + th) / 2
                if best_score is None or score < best_score:
                    best_score, best = score, (px, py, rect)
                if score == 0:
                    break
            if best_score == 0:
                break

        px, py, rect = best
        # Keep the label inside the axes, so a deme near the frame edge does
        # not get its label clipped off the figure.
        bounds = ax.get_window_extent()
        half_w, half_h = tw / 2, th / 2
        px = min(max(px, bounds.x0 + half_w + 2), bounds.x1 - half_w - 2)
        py = min(max(py, bounds.y0 + half_h + 2), bounds.y1 - half_h - 2)
        rect = (px - half_w, py - half_h, px + half_w, py + half_h)
        placed.append(rect)
        inv = ax.transData.inverted()
        results.append((name, text, inv.transform((px, py)), (dx, dy), redge))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auspice-json", required=True)
    parser.add_argument("--output", required=True, help="output path; .pdf also written")
    parser.add_argument("--title", default=None)
    parser.add_argument("--deme-key", default="deme")
    parser.add_argument("--lineage-key", default="major_lineage")
    parser.add_argument("--central-longitude", type=float, default=None,
                        help="default 150 globally, or the centre of --extent")
    parser.add_argument(
        "--extent", nargs=4, type=float, default=None,
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
        help="crop to this box instead of drawing the whole world",
    )
    parser.add_argument(
        "--pie-scale", type=float, default=1.0,
        help="multiply all pie radii (and the spacing/displacement that "
             "depend on them). Use well below 1 when cropping, since radii "
             "are in degrees and a cropped map spans far fewer of them.",
    )
    parser.add_argument("--line-scale", type=float, default=1.0,
                        help="multiply all transmission-line widths")
    parser.add_argument("--fig-width", type=float, default=FIG_WIDTH)
    parser.add_argument("--coastlines", default=None, choices=["110m", "50m", "10m"],
                        help="default 110m globally, 50m when cropped")
    parser.add_argument("--label-fontsize", type=float, default=6.6)
    parser.add_argument(
        "--line-bins", default=None,
        help="two comma-separated upper bounds for the three line-width bins, "
             "e.g. '9,49' gives 1-9 / 10-49 / 50+. The defaults (2,9) suit a "
             "sparse tree; a dense one saturates the top bin and every line "
             "ends up the same thick black, so rescale it to the data.")
    parser.add_argument(
        "--max-lines", type=int, default=0,
        help="draw only the N strongest transitions (0 = all). The single "
             "most effective declutter when a dense tree yields hundreds of "
             "weakly-supported lines.")
    parser.add_argument(
        "--focus-only", action="store_true",
        help="draw only transitions touching --focus-deme. Everything else "
             "is dropped, which isolates the focus deme's story completely.")
    parser.add_argument(
        "--min-slice-fraction", type=float, default=0.0,
        help="merge pie slices below this fraction into a single grey "
             "'other' wedge, so a pie shows structure instead of hairlines.")
    parser.add_argument(
        "--font-scale", type=float, default=1.0,
        help="multiply every font size (deme labels, legends, title, "
             "footnote) and the legend markers. The legend band is measured "
             "and re-packed to match, so legends stay on the page.",
    )
    parser.add_argument("--min-line-count", type=int, default=1)
    parser.add_argument(
        "--focus-deme", default="Singapore",
        help="highlight transitions INTO this deme in red and OUT of it in "
             "blue; pass an empty string to disable",
    )
    parser.add_argument("--no-labels", action="store_true")
    args = parser.parse_args()

    with open(args.auspice_json, encoding="utf-8") as fh:
        data = json.load(fh)
    meta = data["meta"]
    per_deme, transitions = walk_tree(data["tree"], args.deme_key, args.lineage_key)
    coords = deme_coordinates(meta, args.deme_key)
    per_deme = {d: c for d, c in per_deme.items() if d in coords}

    all_lineages = sorted({lin for c in per_deme.values() for lin in c})
    colors = lineage_colors(meta, args.lineage_key, all_lineages)

    # DENV-1 and DENV-2 have 24 Nextclade major lineages where DENV-4 has 7,
    # so the lineage legend's shape varies a lot between serotypes.
    lineage_ncol = max(3, min(5, -(-len(all_lineages) // 6)))

    if args.central_longitude is not None:
        central_longitude = args.central_longitude
    elif args.extent:
        central_longitude = (args.extent[0] + args.extent[1]) / 2.0
    else:
        central_longitude = 150.0
    projection = ccrs.PlateCarree(central_longitude=central_longitude)

    # Pie radii are in degrees, so a cropped map needs much smaller ones.
    scale = args.pie_scale

    # Lay the pies out before building the canvas, so the latitude window can
    # be sized to whatever the relaxation produced.
    totals = {d: sum(c.values()) for d, c in per_deme.items()}
    radii = {d: bin_for(n, SIZE_BINS)[2] * scale for d, n in totals.items()}
    true_pos = {
        d: projection.transform_point(*coords[d], ccrs.PlateCarree())
        for d in per_deme
    }
    pos = relax_positions(true_pos, radii, pad=PIE_PAD * scale,
                          max_shift=MAX_SHIFT * scale)
    leader_min = LEADER_MIN * scale

    # --- count the lines before drawing, so the legends can be built and
    #     measured before the canvas exists -------------------------------
    line_bins = LINE_BINS
    if args.line_bins:
        try:
            lo, hi = (int(x) for x in args.line_bins.split(","))
        except ValueError:
            sys.exit("--line-bins wants two integers, e.g. '9,49'")
        if not 0 < lo < hi:
            sys.exit("--line-bins needs 0 < first < second")
        line_bins = [
            (1, lo, *LINE_BINS[0][2:5], f"1-{lo}"),
            (lo + 1, hi, *LINE_BINS[1][2:5], f"{lo + 1}-{hi}"),
            (hi + 1, 10**9, *LINE_BINS[2][2:5], f"{hi + 1}+"),
        ]

    focus = args.focus_deme or None
    if focus and focus not in pos:
        print(f"  NOTE: focus deme {focus!r} is not on this tree; no highlighting")
        focus = None

    # Background lines first, highlighted ones on top.
    ordered = [
        (pair, count) for pair, count in sorted(
            transitions.items(),
            key=lambda kv: (focus in kv[0] if focus else False, kv[1]),
        )
        if count >= args.min_line_count
        and pair[0] in pos and pair[1] in pos
    ]
    if focus and args.focus_only:
        ordered = [(pair, n) for pair, n in ordered if focus in pair]
    if args.max_lines and len(ordered) > args.max_lines:
        # keep the strongest, but never drop a focus line
        strong = sorted(ordered, key=lambda kv: -kv[1])
        keep, extra = [], []
        for pair, n in strong:
            (keep if (focus and focus in pair) else extra).append((pair, n))
        keep += extra[:max(0, args.max_lines - len(keep))]
        kept_pairs = {pair for pair, _ in keep}
        ordered = [(pair, n) for pair, n in ordered if pair in kept_pairs]

    drawn, focus_counts = Counter(), Counter()
    for (source, target), count in ordered:
        drawn[bin_for(count, line_bins)[5]] += 1
        if focus and target == focus:
            focus_counts["in"] += count
        elif focus and source == focus:
            focus_counts["out"] += count

    # --- build the legends and work out how deep their band must be -------
    fs = args.font_scale
    leg_font, leg_title_font = 8.2 * fs, 9.2 * fs

    # Only legend the lineages that actually survive as their own slice
    # somewhere. With --min-slice-fraction the rest are merged into "other",
    # and listing all 56 SARS-CoV-2 clades would swamp the figure.
    legend_lineages = all_lineages
    if args.min_slice_fraction > 0:
        shown = set()
        for counts in per_deme.values():
            total_tips = sum(counts.values())
            for lineage, n in counts.items():
                if n / total_tips >= args.min_slice_fraction:
                    shown.add(lineage)
        legend_lineages = [l for l in all_lineages if l in shown]
    label_fontsize = args.label_fontsize * fs

    # A scaled-up title may wrap onto several lines; reserve room for them.
    title_size = 14 * fs
    title_wrap = max(20, int(args.fig_width * 72 / (0.52 * title_size)))
    title_text = textwrap.fill(args.title, title_wrap) if args.title else ""
    title_band = (0.18 + (title_text.count(chr(10)) + 1) * title_size / 72 * 1.3
                  + 6 * fs / 72) if title_text else 0.20

    legend_specs = [
        {
            "handles": [
                Line2D([], [], marker="o", linestyle="", markersize=7.5 * fs,
                       markerfacecolor=colors[l], markeredgecolor="white",
                       markeredgewidth=0.5 * fs, label=l)
                for l in legend_lineages
            ] + ([
                Line2D([], [], marker="o", linestyle="", markersize=7.5 * fs,
                       markerfacecolor=OTHER_COLOR, markeredgecolor="white",
                       markeredgewidth=0.5 * fs, label="other")
            ] if args.min_slice_fraction > 0 else []),
            "title": "Nextclade major lineage" + (
                "  (slices <%g%% merged)" % (args.min_slice_fraction * 100)
                if args.min_slice_fraction > 0 else ""),
            "kw": {"ncol": lineage_ncol, "borderpad": 0.8,
                   "handletextpad": 0.6, "columnspacing": 1.3},
        },
        {
            "handles": [
                Line2D([], [], marker="o", linestyle="",
                       markerfacecolor="#cfcfcf", markeredgecolor="#4a4a4a",
                       markeredgewidth=0.55 * fs,
                       markersize=(3.2 + 1.5 * r) * fs,
                       label=f"{lbl} tips")
                for (_lo, _hi, r, lbl) in SIZE_BINS
            ],
            "title": "Tips per deme (binned)",
            "kw": {"borderpad": 0.9, "labelspacing": 1.5, "handletextpad": 1.2},
        },
        {
            "handles": [
                Line2D([], [], color=LINE_COLOR, linewidth=w * args.line_scale,
                       alpha=a, label=f"{lbl} events  (n={drawn.get(lbl, 0)})")
                for (_lo, _hi, w, a, _arrow, lbl) in line_bins
            ],
            "title": "Inferred deme transitions (binned)",
            "kw": {"borderpad": 0.9, "labelspacing": 1.4, "handlelength": 3.6},
        },
    ]
    if focus:
        legend_specs.append({
            "handles": [
                Line2D([], [], color=IN_COLOR, linewidth=3.0 * args.line_scale,
                       alpha=FOCUS_MIN_ALPHA,
                       label=f"into {focus}  ({focus_counts['in']} events)"),
                Line2D([], [], color=OUT_COLOR, linewidth=3.0 * args.line_scale,
                       alpha=FOCUS_MIN_ALPHA,
                       label=f"out of {focus}  ({focus_counts['out']} events)"),
            ],
            "title": f"{focus} focus",
            "kw": {"borderpad": 0.9, "labelspacing": 1.4, "handlelength": 3.6},
        })

    legend_sizes = [
        measure_legend(args.fig_width, spec, leg_font, leg_title_font)
        for spec in legend_specs
    ]
    legend_rows = pack_rows(legend_specs, legend_sizes, args.fig_width * 0.96)
    row_heights = [max(h for _s, _w, h in row) for row in legend_rows]
    footnote_in = 2 * (7.4 * fs / 72.0) * 1.55 + 0.20
    legend_band = (sum(row_heights) + 0.30 * len(legend_rows)
                   + footnote_in + 0.18)

    if args.extent:
        lon_min, lon_max, lat_min, lat_max = args.extent
        lon_span = lon_max - lon_min
    else:
        lon_min = lon_max = None
        lon_span = 360.0
        tops = [pos[d][1] + radii[d] for d in pos] + [t[1] for t in true_pos.values()]
        bottoms = [pos[d][1] - radii[d] for d in pos] + [t[1] for t in true_pos.values()]
        lat_max = min(LAT_MAX_BOUNDS[1], max(LAT_MAX_BOUNDS[0], max(tops) + LAT_PAD))
        lat_min = max(LAT_MIN_BOUNDS[0], min(LAT_MIN_BOUNDS[1], min(bottoms) - LAT_PAD))

    # Size the canvas to the map's aspect ratio, plus bands for title/legends.
    ax_width_in = args.fig_width * AX_WIDTH_FRAC
    ax_height_in = ax_width_in / (lon_span / (lat_max - lat_min))
    fig_height = ax_height_in + legend_band + title_band
    legend_y = legend_band / fig_height
    fig = plt.figure(figsize=(args.fig_width, fig_height))
    ax = fig.add_axes(
        [(1 - AX_WIDTH_FRAC) / 2, legend_y, AX_WIDTH_FRAC,
         ax_height_in / fig_height],
        projection=projection,
    )
    resolution = args.coastlines or ("50m" if args.extent else "110m")
    ax.add_feature(cfeature.LAND.with_scale(resolution), facecolor="#f4f3ee",
                   edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor="#ffffff",
                   edgecolor="none", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale(resolution), linewidth=0.35,
                   edgecolor="#9a9a9a", zorder=1)
    ax.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.22,
                   edgecolor="#c9c9c9", zorder=1)
    if args.extent:
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    else:
        ax.set_global()
        ax.set_ylim(lat_min, lat_max)
    ax.spines["geo"].set_edgecolor("#999999")
    ax.spines["geo"].set_linewidth(0.6)

    # --- transmission lines (under the pies) -------------------------------
    # Settle the layout so transData is final before converting pie radii
    # (degrees) into the point units that shrinkA/shrinkB expect.
    fig.canvas.draw()

    def radius_points(deme):
        """The deme's drawn radius, including its white halo, in points."""
        cx, _cy = ax.transData.transform(pos[deme])
        ex, _ey = ax.transData.transform(
            (pos[deme][0] + radii[deme] * 1.1, pos[deme][1]))
        return abs(ex - cx) * 72.0 / fig.dpi

    def styled(source, target, width, alpha):
        """Colour, alpha, width and z-order for one transmission line."""
        if focus and target == focus:
            return IN_COLOR, max(alpha, FOCUS_MIN_ALPHA), max(width, FOCUS_MIN_WIDTH), 3.7
        if focus and source == focus:
            return OUT_COLOR, max(alpha, FOCUS_MIN_ALPHA), max(width, FOCUS_MIN_WIDTH), 3.6
        return LINE_COLOR, alpha, width, 2 + width / 10

    for (source, target), count in ordered:
        _lo, _hi, width, alpha, arrow, _label = bin_for(count, line_bins)
        color, alpha, width, zorder = styled(source, target, width, alpha)
        width *= args.line_scale
        shrink_a = radius_points(source) + 1.5
        shrink_b = radius_points(target) + 2.0 + (0.9 * width if arrow else 0.0)
        (x1, y1), (x2, y2) = pos[source], pos[target]
        for off1, off2 in wrapped_offsets(x1, x2):
            ax.add_patch(FancyArrowPatch(
                (x1 + off1, y1), (x2 + off2, y2),
                connectionstyle="arc3,rad=0.16",
                arrowstyle="-|>" if arrow else "-",
                mutation_scale=8 + 2.0 * width,
                linewidth=width, color=color, alpha=alpha,
                zorder=zorder, capstyle="round", joinstyle="round",
                # Stop short of the pies so the arrow tip stays visible
                # rather than being drawn underneath the target pie.
                shrinkA=shrink_a, shrinkB=shrink_b,
            ))

    # --- leader lines back to true centroids --------------------------------
    for deme in per_deme:
        (tx, ty), (px, py) = true_pos[deme], pos[deme]
        if np.hypot(px - tx, py - ty) > leader_min:
            ax.plot([tx, px], [ty, py], color="#7a7a7a", linewidth=0.5,
                    alpha=0.75, zorder=3, solid_capstyle="round")
            ax.plot([tx], [ty], marker="o", markersize=1.7, color="#5a5a5a",
                    zorder=3.5, linestyle="")

    # --- pies ---------------------------------------------------------------
    for deme in sorted(per_deme, key=lambda d: -totals[d]):
        counts, total, radius = per_deme[deme], totals[deme], radii[deme]
        x, y = pos[deme]
        ax.add_patch(Circle((x, y), radius * 1.1, facecolor="white",
                            edgecolor="none", alpha=0.8, zorder=4))
        slices = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        if args.min_slice_fraction > 0:
            big = [(k, v) for k, v in slices if v / total >= args.min_slice_fraction]
            small = sum(v for k, v in slices if v / total < args.min_slice_fraction)
            if small:
                big.append(("other", small))
            slices = big
        start = 90.0
        for lineage, n in slices:
            extent = 360.0 * n / total
            ax.add_patch(Wedge((x, y), radius, start - extent, start,
                               facecolor=colors.get(lineage, OTHER_COLOR),
                               edgecolor="white",
                               linewidth=0.4, zorder=5))
            start -= extent
        ax.add_patch(Circle((x, y), radius, facecolor="none",
                            edgecolor="#4a4a4a", linewidth=0.5, zorder=6))

    # --- labels -------------------------------------------------------------
    # When cropped, only label demes actually inside the frame; the rest are
    # off-map and only contribute lines running off the edge.
    x_lo, x_hi = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    visible = [
        d for d in per_deme
        if x_lo <= pos[d][0] <= x_hi and y_lo <= pos[d][1] <= y_hi
    ]

    if not args.no_labels:
        fig.canvas.draw()
        pie_display = []
        for deme in visible:
            cx, cy = ax.transData.transform(pos[deme])
            redge = abs(ax.transData.transform(
                (pos[deme][0] + radii[deme], pos[deme][1]))[0] - cx)
            pie_display.append((cx, cy, redge))

        entries = [
            (d, f"{d} ({totals[d]})", pos[d], radii[d])
            for d in sorted(visible, key=lambda d: -totals[d])
        ]
        for name, text, (lx, ly), (ax_, ay), redge in place_labels(
                ax, fig, entries, pie_display, fontsize=label_fontsize):
            ax.text(lx, ly, text, ha="center", va="center", fontsize=label_fontsize,
                    color="#141414", zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.1, foreground="white")])
            # short connector when the label sits well away from its pie
            dist = np.hypot(lx - ax_, ly - ay)
            if dist > radii[name] * 2.6:
                ux, uy = (lx - ax_) / dist, (ly - ay) / dist
                ax.plot([ax_ + ux * radii[name], lx - ux * 1.4],
                        [ay + uy * radii[name], ly - uy * 1.4],
                        color="#8a8a8a", linewidth=0.42, alpha=0.8, zorder=7)

    # --- legends, packed into the measured band below the map -------------
    x_margin = args.fig_width * 0.02
    y_cursor = legend_band - 0.10
    for row, row_height in zip(legend_rows, row_heights):
        x_cursor = x_margin
        for spec, width, _height in row:
            legend = fig.legend(
                handles=spec["handles"], title=spec["title"], loc="upper left",
                bbox_to_anchor=(x_cursor / args.fig_width,
                                y_cursor / fig_height),
                frameon=True, fontsize=leg_font, title_fontsize=leg_title_font,
                **spec.get("kw", {}),
            )
            legend.get_frame().set_edgecolor("#bbbbbb")
            legend.get_frame().set_linewidth(0.6 * fs)
            x_cursor += width + 0.35
        y_cursor -= row_height + 0.30

    fig.text(x_margin / args.fig_width, 0.10 / fig_height,
             "Pie area is binned, not proportional. Grey leader lines join a "
             "displaced pie to its true centroid (dot).\nArrowheads show "
             "inferred direction for bins of 3+ events.",
             fontsize=7.4 * fs, color="#555555", va="bottom")

    if title_text:
        ax.set_title(title_text, fontsize=title_size, pad=6 * fs)

    fig.savefig(args.output, dpi=330, facecolor="white")
    pdf = args.output.rsplit(".", 1)[0] + ".pdf"
    fig.savefig(pdf, facecolor="white")

    print(f"wrote {args.output} and {pdf}")
    print(f"  demes {len(per_deme)} | tips {sum(totals.values())} | "
          f"lineages {len(all_lineages)}")
    print(f"  lines drawn by bin: {dict(drawn)}  (total {sum(drawn.values())})")
    print(f"  demes by size bin: "
          f"{dict(Counter(bin_for(n, SIZE_BINS)[3] for n in totals.values()))}")
    if focus:
        print(f"  {focus}: {focus_counts['in']} events in (red), "
              f"{focus_counts['out']} events out (blue)")
    shifted = sum(1 for d in per_deme
                  if np.hypot(*(np.array(pos[d]) - np.array(true_pos[d]))) > leader_min)
    print(f"  pies displaced for legibility: {shifted}/{len(per_deme)}")


if __name__ == "__main__":
    main()
