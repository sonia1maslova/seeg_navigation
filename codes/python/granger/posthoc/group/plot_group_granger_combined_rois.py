"""Stack four existing Granger ROI results without recalculating statistics.

Read only the formal per-pair plotting tables. Preserve the approved compact
style, separate alignment branches/processes, and the existing shared legend.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from granger_compact_style import (
    STYLE, configure_font, pale, roi_label_lookup, validate_plot_grid,
)

SOURCE_ROI = "Temporal_Mid L_Anterior"
TARGET_ROIS = (
    "Hippocampus L_Anterior", "Hippocampus L_Posterior",
    "ParaHippocampal L", "OFCpost L",
)
BRANCHES = ("order_aligned", "order_aligned_and_terminal_aligned")
PROCESSES = ("corridor", "junction")
METRICS = [
    {"label": "T", "title": "Time-domain GC"},
    {"label": "D", "title": "0.5–4 Hz"},
    {"label": "TA", "title": "4–12 Hz"},
    {"label": "B", "title": "12–30 Hz"},
    {"label": "LG", "title": "30–70 Hz"},
    {"label": "HG", "title": "70–150 Hz"},
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_row(formal_root, branch, process, target):
    stem = (f"{process}_{SOURCE_ROI}_to_{target}_groupP_pairSigProp_byMetric"
            .replace(" ", "_"))
    base = formal_root / branch / process / "group_p_and_pair_sig_proportion"
    csv_path = base / "data" / f"{stem}_plot_data.csv"
    frame = pd.read_csv(csv_path)
    for column, expected in (("branch", branch), ("process_type", process),
                             ("brain_area1", SOURCE_ROI), ("brain_area2", target)):
        if frame.empty or not frame[column].eq(expected).all():
            raise ValueError(f"Incorrect {column} in {csv_path}")
    segments = frame[["data_cut", "segment_label"]].drop_duplicates()
    if segments["data_cut"].duplicated().any():
        raise ValueError(f"Ambiguous segment labels in {csv_path}")
    cuts = segments["data_cut"].tolist()
    labels = segments["segment_label"].tolist()
    validate_plot_grid(frame, cuts, METRICS)
    np.testing.assert_allclose(frame["neglog10_p"], -np.log10(frame["group_p"]),
                               rtol=1e-12, atol=1e-12, equal_nan=True)
    valid_pairs = frame["n_valid_pairs"].to_numpy(dtype=float)
    expected_fraction = np.divide(
        frame["n_sig_pairs"].to_numpy(dtype=float), valid_pairs,
        out=np.full_like(valid_pairs, np.nan), where=valid_pairs > 0,
    )
    np.testing.assert_allclose(frame["sig_proportion"], expected_fraction,
                               rtol=1e-12, atol=1e-12, equal_nan=True)
    return {"data": frame, "cuts": cuts, "labels": labels, "target": target,
            "csv": csv_path, "original_png": base / "fig" / f"{stem}.png"}


def draw_row(subfig, row, metric_colors, source_label, target_label):
    """Retain the original compact renderer's artists and physical row size."""
    axes = subfig.subplots(1, 6, sharey=True)
    subfig.subplots_adjust(left=0.053, right=0.957, bottom=0.20, top=0.82,
                           wspace=0.17)
    x = np.arange(len(row["cuts"]), dtype=float)
    for panel, metric in enumerate(METRICS):
        ax, colour = axes[panel], metric_colors[metric["label"]]
        ax2 = ax.twinx()
        for direction in ("1to2", "2to1"):
            sub = row["data"][(row["data"]["metric_label"] == metric["label"])
                              & (row["data"]["direction"] == direction)].set_index("data_cut")
            y = sub.loc[row["cuts"], "neglog10_p"].to_numpy(dtype=float)
            fractions = sub.loc[row["cuts"], "sig_proportion"].to_numpy(dtype=float)
            positions = x + STYLE["direction_x_offsets"][direction]
            bars = ax2.bar(
                positions, fractions, width=STYLE["bar_width"],
                facecolor=pale(colour, STYLE["bar_colour_strength"][direction]),
                edgecolor="none" if direction == "1to2" else pale(colour, 0.58),
                linewidth=0.55, zorder=1,
            )
            line, = ax.plot(positions, y, color=colour,
                            linewidth=STYLE["line_width_pt"],
                            alpha=STYLE["line_alpha"], solid_capstyle="round", zorder=4)
            points = ax.scatter(
                positions, y, marker="s", s=29,
                facecolors=colour if direction == "1to2" else "white",
                edgecolors="white" if direction == "1to2" else colour,
                linewidths=0.45 if direction == "1to2" else 1.15, zorder=6,
            )
            # Verify displayed values, including the two different y axes.
            np.testing.assert_array_equal(line.get_ydata(), y)
            offsets = points.get_offsets()
            finite = np.isfinite(y)
            np.testing.assert_array_equal(offsets[finite, 0], positions[finite])
            np.testing.assert_array_equal(offsets[finite, 1], y[finite])
            np.testing.assert_array_equal([bar.get_height() for bar in bars], fractions)
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.axhline(-math.log10(0.05), color="#565656", linewidth=0.9,
                   linestyle=(0, (3.5, 2.5)), zorder=3)
        ax.set_ylim(*STYLE["left_axis_limits"])
        ax.set_yticks([0, 1, 2, 3])
        ax.set_xlim(-0.58, len(row["cuts"]) - 0.42)
        ax.set_xticks(x)
        ax.set_xticklabels(row["labels"])
        ax.set_title(metric["title"], fontsize=12, color=colour, pad=10)
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
        ax2.tick_params(axis="y", length=2.7 if panel == 5 else 0,
                        width=0.8, labelright=panel == 5, pad=3)
        if panel == 0:
            ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=12, labelpad=7)
        if panel == 5:
            ax2.set_ylabel("Sig. pair proportion", fontsize=11.5, labelpad=8)
    subfig.text(0.5, 0.965, source_label + " — " + target_label,
                ha="center", va="center", fontsize=15, fontweight="bold")
    subfig.text(0.5, 0.068, "Segment", ha="center", va="center",
                fontsize=11.5, fontweight="bold")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.environ.get(
        "SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")))
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    formal = args.project_root / "figures/granger/formal_20261003"
    output_root = args.output_dir or formal / "combined_four_rois"
    config_path = args.project_root / "config/project_config.json"
    labels = roi_label_lookup(json.loads(config_path.read_text(encoding="utf-8")))
    existing_sources = sorted(
        path for branch in BRANCHES for path in (formal / branch).rglob("*")
        if path.is_file() and path.suffix in (".png", ".csv", ".json")
    ) + [config_path, Path(__file__).with_name("granger_compact_style.py")]
    protected = {str(path): digest(path) for path in existing_sources}
    output_paths = [output_root / branch /
                    f"{process}_aMTG_four_ROIs_groupP_pairSigProp_byMetric.png"
                    for branch in BRANCHES for process in PROCESSES]
    manifest_path = output_root / "combined_plot_manifest.json"
    for path in [*output_paths, manifest_path]:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite an existing combined output: {path}")
    configure_font()
    outputs = []
    width, height = STYLE["figure_size_inches"]
    for branch in BRANCHES:
        branch_manifest = formal / branch / "reference_plot_manifest.json"
        metric_colors = json.loads(branch_manifest.read_text(encoding="utf-8"))["metric_colors"]
        for process in PROCESSES:
            rows = [load_row(formal, branch, process, roi) for roi in TARGET_ROIS]
            if any(row["cuts"] != rows[0]["cuts"] or row["labels"] != rows[0]["labels"]
                   for row in rows[1:]):
                raise ValueError(f"Inconsistent segment order: {branch}, {process}")
            output = output_root / branch / f"{process}_aMTG_four_ROIs_groupP_pairSigProp_byMetric.png"
            fig = plt.figure(figsize=(width, height * len(rows)))
            try:
                subfigs = fig.subfigures(len(rows), 1, hspace=0.0)
                for subfig, row in zip(subfigs, rows):
                    draw_row(subfig, row, metric_colors, labels(SOURCE_ROI), labels(row["target"]))
                fig.canvas.draw()
                if fig.legends or any(ax.get_legend() for ax in fig.axes):
                    raise AssertionError("Unexpected embedded legend")
                output.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(output, dpi=STYLE["dpi"], bbox_inches="tight", pad_inches=0.07)
                outputs.append({"png": output.relative_to(output_root).as_posix(),
                                "branch": branch, "process": process, "roi_rows": list(TARGET_ROIS),
                                "data_sources": [row["csv"].relative_to(formal).as_posix() for row in rows],
                                "data_sha256": [digest(row["csv"]) for row in rows],
                                "rows_in_plot_tables": sum(len(row["data"]) for row in rows),
                                "displayed_values_exactly_equal_to_plot_tables": True,
                                "metric_colors": metric_colors, "sha256": digest(output)})
                print(f"Saved: {output}", flush=True)
            finally:
                plt.close(fig)
    for path, before_hash in protected.items():
        if digest(path) != before_hash:
            raise RuntimeError(f"Original source changed during rendering: {path}")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "layout": "4 ROI-pair rows x 6 metrics; both directions per panel",
              "target_roi_order": list(TARGET_ROIS), "style": STYLE,
              "figure_size_inches": [width, height * len(TARGET_ROIS)],
              "statistics_recomputed": False, "multiple_comparison_correction": "none",
              "original_sources_unchanged": True, "protected_source_sha256": protected,
              "shared_legend": "../direction_legend.png", "outputs": outputs,
              "code_sha256": digest(__file__)}
    manifest_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(outputs)} combined figures; {len(protected)} original files unchanged.")


if __name__ == "__main__":
    main()
