"""Shared, statistics-preserving compact Granger rendering.

One row of six metrics, both directions in each panel. Plot inputs are the
existing group P values and significant-pair fractions, never re-estimated.
"""
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.lines import Line2D
import numpy as np

STYLE = {
    "layout": "1 row x 6 metrics; both directions in each panel",
    "figure_size_inches": [16.8, 3.0], "dpi": 300,
    "left_axis_limits": [0.0, 3.5], "right_axis_limits": [0.0, 0.85],
    "direction_x_offsets": {"1to2": -0.17, "2to1": 0.17},
    "direction_markers": {"1to2": "filled square", "2to1": "hollow square"},
    "line_width_pt": 0.8, "line_alpha": 0.4,
    "bar_width": 0.26, "bar_colour_strength": {"1to2": 0.28, "2to1": 0.09},
    "bar_hatching": "none", "main_figure_legend": False,
    "significance_line": "uncorrected P = 0.05",
    "roi_label_source": "config/project_config.json:roi_abbreviations",
}


def pale(colour, strength):
    rgb = np.asarray(colors.to_rgb(colour))
    return tuple(1 - strength * (1 - rgb))


def configure_font():
    available = {entry.name for entry in font_manager.fontManager.ttflist}
    plt.rcParams.update({
        "font.family": "Arial" if "Arial" in available else "DejaVu Sans",
        "font.size": 10.5, "axes.linewidth": 1.05,
        "axes.labelweight": "bold", "axes.titleweight": "bold",
        "xtick.labelsize": 9.5, "ytick.labelsize": 10,
        "savefig.facecolor": "white",
    })


def roi_label_lookup(config):
    """Canonical names only; normalize spaces/underscores, not abbreviations."""
    normalize = lambda value: " ".join(str(value).replace("_", " ").split()).casefold()
    labels = {}
    for entry in config["roi_abbreviations"]:
        key = normalize(entry["full_name"])
        if key in labels:
            raise ValueError("Duplicate normalized ROI key in shared project config")
        labels[key] = entry["abbreviation"]
    def resolve(name):
        key = normalize(name)
        if key not in labels:
            raise KeyError("ROI missing from shared abbreviation config: " + str(name))
        return labels[key]
    return resolve


def validate_plot_grid(plot_table, data_cuts, metrics):
    required = {"direction", "data_cut", "metric_label", "neglog10_p", "sig_proportion"}
    if not required.issubset(plot_table.columns):
        raise ValueError("Missing plot-table fields: " + str(required - set(plot_table.columns)))
    keys = ["metric_label", "direction", "data_cut"]
    if plot_table.duplicated(keys).any():
        raise ValueError("Duplicate direction / metric / segment plotting keys")
    actual = set(plot_table[keys].itertuples(index=False, name=None))
    expected = {(m["label"], d, c) for m in metrics
                for d in ("1to2", "2to1") for c in data_cuts}
    if actual != expected:
        raise ValueError("Incomplete direction / metric / segment plotting grid")
    for column, upper in (("neglog10_p", STYLE["left_axis_limits"][1]),
                          ("sig_proportion", STYLE["right_axis_limits"][1])):
        array = plot_table[column].to_numpy(dtype=float)
        finite = array[np.isfinite(array)]
        if np.isinf(array).any() or (finite < 0).any() or (finite > upper).any():
            raise ValueError("Plot range would truncate " + column)


def render_compact_figure(plot_table, data_cuts, metrics, metric_colors,
                          roi1_label, roi2_label, segment_labels, output_png,
                          dpi=300):
    validate_plot_grid(plot_table, data_cuts, metrics)
    if len(metrics) != 6 or len(segment_labels) != len(data_cuts):
        raise ValueError("Six metrics and matched segment labels are required")
    configure_font()
    fig, axes = plt.subplots(1, 6, figsize=STYLE["figure_size_inches"], sharey=True)
    # Same physical plotting height as the approved sample, minus its legend band.
    fig.subplots_adjust(left=0.053, right=0.957, bottom=0.20, top=0.82, wspace=0.17)
    x = np.arange(len(data_cuts), dtype=float)
    try:
        for panel, metric in enumerate(metrics):
            ax, colour = axes[panel], metric_colors[metric["label"]]
            ax2 = ax.twinx()
            for direction in ("1to2", "2to1"):
                sub = plot_table[(plot_table["metric_label"] == metric["label"])
                                 & (plot_table["direction"] == direction)].set_index("data_cut")
                y = sub.loc[data_cuts, "neglog10_p"].to_numpy(dtype=float)
                fractions = sub.loc[data_cuts, "sig_proportion"].to_numpy(dtype=float)
                positions = x + STYLE["direction_x_offsets"][direction]
                ax2.bar(positions, fractions, width=STYLE["bar_width"],
                        facecolor=pale(colour, STYLE["bar_colour_strength"][direction]),
                        edgecolor="none" if direction == "1to2" else pale(colour, 0.58),
                        linewidth=0.55, zorder=1)
                ax.plot(positions, y, color=colour, linewidth=STYLE["line_width_pt"],
                        alpha=STYLE["line_alpha"], solid_capstyle="round", zorder=4)
                ax.scatter(positions, y, marker="s", s=29,
                           facecolors=colour if direction == "1to2" else "white",
                           edgecolors="white" if direction == "1to2" else colour,
                           linewidths=0.45 if direction == "1to2" else 1.15, zorder=6)
            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)
            ax.axhline(-math.log10(0.05), color="#565656", linewidth=0.9,
                       linestyle=(0, (3.5, 2.5)), zorder=3)
            ax.set_ylim(*STYLE["left_axis_limits"])
            ax.set_yticks([0, 1, 2, 3])
            ax.set_xlim(-0.58, len(data_cuts) - 0.42)
            ax.set_xticks(x)
            ax.set_xticklabels(segment_labels)
            ax.set_title(metric.get("title", metric.get("legend", metric["label"])),
                         fontsize=12, color=colour, pad=10)
            ax.tick_params(axis="both", length=2.7, width=0.8, pad=3,
                           color="#393939", labelleft=panel == 0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color("#393939")
            ax2.set_ylim(*STYLE["right_axis_limits"])
            ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
            for side in ("top", "left", "bottom"):
                ax2.spines[side].set_visible(False)
            ax2.spines["right"].set_visible(panel == 5)
            ax2.spines["right"].set_color("#393939")
            ax2.tick_params(axis="y", length=2.7 if panel == 5 else 0, width=0.8,
                            labelright=panel == 5, pad=3)
            if panel == 0:
                ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=12, labelpad=7)
            if panel == 5:
                ax2.set_ylabel("Sig. pair proportion", fontsize=11.5, labelpad=8)
        fig.text(0.5, 0.965, roi1_label + " — " + roi2_label,
                 ha="center", va="center", fontsize=15, fontweight="bold")
        fig.text(0.5, 0.068, "Segment", ha="center", va="center",
                 fontsize=11.5, fontweight="bold")
        fig.canvas.draw()
        out = Path(output_png)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.07)
        if not out.is_file() or out.stat().st_size == 0:
            raise RuntimeError("PNG export failed: " + str(out))
        return {"png": str(out), "rows": len(plot_table), "style": STYLE}
    finally:
        plt.close(fig)


def render_direction_legend(roi1_label, output_png, dpi=300):
    """Two square markers only; the target is the second ROI in each title."""
    configure_font()
    handles = [Line2D([], [], linestyle="none", marker="s", markersize=6,
                      markerfacecolor="#424242", markeredgecolor="#424242"),
               Line2D([], [], linestyle="none", marker="s", markersize=6,
                      markerfacecolor="white", markeredgecolor="#424242", markeredgewidth=1.1)]
    labels = [roi1_label + " → Target ROI", "Target ROI → " + roi1_label]
    fig = plt.figure(figsize=(8, 0.55))
    try:
        legend = fig.legend(handles, labels, loc="center", ncol=2, frameon=False,
                            fontsize=11.5, handlelength=1.0, handletextpad=0.65,
                            columnspacing=2.7, borderaxespad=0)
        fig.canvas.draw()
        bbox = legend.get_window_extent(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        out = Path(output_png)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches=bbox.expanded(1.05, 1.25), pad_inches=0.025)
        return {"png": str(out), "labels": labels, "symbols": ["filled square", "hollow square"]}
    finally:
        plt.close(fig)
