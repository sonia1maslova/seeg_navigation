#!/usr/bin/env python3
"""Render compact raw-PAC means; calculation lives in a separate script.

Plot styling and descriptive OLS trend/95% mean-CI calculations are preserved
from plot_lina_pac_by_segment_formal.py. The three trend lines are visual
summaries, not the LMM fixed effects. All actual positions remain on the axis;
missing positions are not filled or interpolated with synthetic PAC values.
Only PNG and their reproducible CSV tables are written (no PDF/SVG).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
from scipy.stats import linregress, t


COLORS = {"corridor": "#2A9D8F", "junction": "#7B61A8"}
MARKERS = {"corridor": "o", "junction": "s"}
BRANCHES = ("include_last", "exclude_terminal_corridor_and_junction")


def short_label(segment_id: str) -> str:
    number = re.search(r"\d+", segment_id).group()
    if "first" in segment_id:
        return f"C{number}-1st"
    if "second" in segment_id:
        return f"C{number}-2nd"
    return f"J{number}"


def add_fit(ax, rows: pd.DataFrame, color: str) -> None:
    x, y = rows.segment_position.to_numpy(float), rows.mean_raw_PAC.to_numpy(float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 2 or len(np.unique(x)) < 2:
        return
    fit = linregress(x, y)
    grid = np.linspace(float(x.min()) - 0.2, float(x.max()) + 0.2, 300)
    prediction = fit.intercept + fit.slope * grid
    ax.plot(grid, prediction, color=color, lw=1.65, ls="-", zorder=2)
    if len(x) > 2:
        residual = y - (fit.intercept + fit.slope * x)
        mse = float(np.sum(residual ** 2) / (len(x) - 2))
        sxx = float(np.sum((x - x.mean()) ** 2))
        se = np.sqrt(mse * (1 / len(x) + (grid - x.mean()) ** 2 / sxx))
        critical = float(t.ppf(0.975, len(x) - 2))
        ax.fill_between(grid, prediction - critical * se, prediction + critical * se,
                        color=color, alpha=0.10, linewidth=0, zorder=1)


def save_legend(output: Path) -> Path:
    handles = [
        Line2D([0], [0], color="#333333", lw=2.0, label="Overall trend"),
        Line2D([0], [0], color=COLORS["corridor"], marker="o", markerfacecolor="white",
               markeredgewidth=1.4, markersize=7, lw=1.4, label="Corridor"),
        Line2D([0], [0], color=COLORS["junction"], marker="s", markerfacecolor="white",
               markeredgewidth=1.4, markersize=7, lw=1.4, label="Junction"),
    ]
    fig = plt.figure(figsize=(5.8, 0.72), dpi=300)
    fig.legend(handles=handles, loc="center", ncol=3, frameon=False,
               fontsize=10.5, handlelength=1.7, columnspacing=1.2, handletextpad=0.55)
    output.mkdir(parents=True, exist_ok=True)
    path = output / "pac_by_segment_legend.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return path


def render(rows: pd.DataFrame, output: Path, roi_label: str, stem: str, branch_axis: pd.DataFrame) -> Path:
    rows = rows.sort_values("segment_position").copy()
    if rows.segment_position.duplicated().any():
        raise ValueError("Repeated plot positions")
    rows["short_label"] = rows.segment_id.map(short_label)
    fig, ax = plt.subplots(figsize=(4.55, 3.25), dpi=300)
    for kind in ("corridor", "junction"):
        subset = rows[rows.segment_type.eq(kind)]
        if subset.empty:
            continue
        ax.errorbar(
            subset.segment_position, subset.mean_raw_PAC, yerr=subset.sem_raw_PAC,
            fmt=MARKERS[kind], color=COLORS[kind], markerfacecolor="white",
            markeredgewidth=1.0, markersize=4.0, elinewidth=0.72, capsize=1.7, zorder=3,
        )
        add_fit(ax, subset, COLORS[kind])
    add_fit(ax, rows, "#333333")
    ax.set_xlabel("Segment", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("PAC", fontsize=10.5, fontweight="bold")
    ax.set_title(roi_label, fontsize=13.0, fontweight="bold", pad=5)
    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, -3))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.get_offset_text().set_fontsize(8.5)
    ax.margins(y=0.12)
    ax.set_xticks(branch_axis.segment_position)
    ax.set_xticklabels(branch_axis.segment_id.map(short_label), rotation=50, ha="right", fontsize=6.2)
    ax.set_xlim(float(branch_axis.segment_position.min()) - 0.5, float(branch_axis.segment_position.max()) + 0.5)
    ax.tick_params(axis="y", labelsize=8.2)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.35)
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"{stem}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    rows.to_csv(output / f"{stem}_values.csv", index=False)
    plt.close(fig)
    return path


def plot_all(project: Path, data_root: Path, output_root: Path) -> dict:
    with (project / "config/project_config.json").open(encoding="utf-8") as stream:
        settings = json.load(stream)
    labels = {item["full_name"]: item["abbreviation"] for item in settings["roi_abbreviations"]}
    files, counts = [], {}
    for branch in BRANCHES:
        df = pd.read_csv(data_root / branch / "group_mean_pac.csv")
        mapping = df[["segment_position", "segment_id"]].drop_duplicates()
        if mapping.segment_position.duplicated().any():
            raise ValueError("Conflicting chronological position identifiers")
        axis = mapping.sort_values("segment_position")
        count = 0
        for (roi, phase, amp), rows in df.groupby(["roi", "phase_band", "amp_band"], sort=True):
            if roi not in labels:
                raise ValueError(f"Missing standard abbreviation for {roi}")
            roi_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", roi).strip("_")
            stem = f"pac_by_segment_{roi_name}_{phase}_{amp}"
            path = render(rows, output_root / branch / roi_name, labels[roi], stem, axis)
            files.append(path.relative_to(output_root).as_posix())
            count += 1
        if count != 24:
            raise ValueError(f"Expected 4 ROIs x 6 combinations, got {count} in {branch}")
        counts[branch] = {"n_png": count, "n_axis_positions": len(axis), "first_position": int(axis.segment_position.min()), "last_position": int(axis.segment_position.max())}
    legend = save_legend(output_root / "legend")
    report = {"n_data_figures": len(files), "legend": legend.relative_to(output_root).as_posix(),
              "source_root": "results/pac/by_lina/pac_segment_effect/mean_trends/adjusted",
              "figure_size_inches": [4.55, 3.25], "dpi": 300, "branches": counts,
              "trends": "Descriptive OLS fit to subject-balanced group means; shaded 95% confidence interval of fitted mean",
              "error_bars": "Between-subject SEM; no FDR or per-segment test performed", "png_files": files}
    (output_root / "plot_manifest.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.getenv("SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")))
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    data = args.data_root or args.project_root / "results/pac/by_lina/pac_segment_effect/mean_trends/adjusted"
    output = args.output_root or args.project_root / "figures/pac/by_lina/pac_segment_effect/mean_trends/adjusted"
    report = plot_all(args.project_root, data, output)
    print(json.dumps({"n_png": report["n_data_figures"], "branches": report["branches"], "output_root": str(output)}, indent=2))


if __name__ == "__main__":
    main()
