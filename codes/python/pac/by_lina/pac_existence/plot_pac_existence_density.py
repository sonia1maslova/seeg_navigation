"""Preview pooled PAC/MI existence using saved surrogate density curves.

This is a visualization-only companion to plot_pac_existence_forest.py.
Each of the 24 cells uses a local raw-MI axis. KDE is a display aid only;
original empirical P and original BH q determine the significance marks.
No existing source result or forest figure is replaced.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

SCALE = 1000.0
BLUE = "#398EAD"
INK = "#242424"
FIGSIZE = (7.55, 5.15)
GRID_SIZE = 1024


def load_validated_data(project: Path, scope: str):
    """Reuse the established independent checks against Lina's originals."""
    companion_path = project / "codes/python/pac/by_lina/pac_existence/plot_pac_existence_forest.py"
    spec = importlib.util.spec_from_file_location("validated_mi_forest", companion_path)
    if spec is None or spec.loader is None:
        raise ImportError(str(companion_path))
    companion = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(companion)
    data, protected, labels = companion.load_data(project)
    protected[str(companion_path)] = companion.sha256(companion_path)
    selected = data.loc[data.scope.eq(scope)].copy()
    if len(selected) != 24 or selected.duplicated(companion.KEYS).any():
        raise ValueError("Expected 24 unique pooled results for the selected original scope")
    return companion, selected, protected, labels


def kde(grid: np.ndarray, samples: np.ndarray, bandwidth: float) -> np.ndarray:
    """Gaussian KDE, using Scott's n^(-1/5) bandwidth times sample SD."""
    if len(samples) < 2 or not np.isfinite(samples).all() or bandwidth <= 0:
        raise ValueError("Cannot estimate a finite surrogate density")
    u = (grid[:, None] - samples[None, :]) / bandwidth
    return np.exp(-0.5 * u * u).mean(axis=1) / (bandwidth * np.sqrt(2 * np.pi))


def prepare_curves(companion, data: pd.DataFrame, protected: dict):
    cache = {}
    curves = {}
    rows = []
    for _, row in data.iterrows():
        path = Path(row.source_null_file)
        if str(path) not in cache:
            cache[str(path)] = companion.read_csv(path, protected)
        table = cache[str(path)]
        subset = table.loc[table.trajectory_scope.eq(row.scope) & table.location1.eq(row.roi)]
        if len(subset) != 1000 or subset.surrogate_index.duplicated().any():
            raise ValueError("Missing or duplicate saved surrogate")
        if not subset.phase_band.eq(row.phase_band).all() or not subset.amp_band.eq(row.amp_band).all():
            raise ValueError("Surrogate frequency key mismatch")
        values = subset.group_surrogate_raw_pac.to_numpy(float)
        observed = float(row.group_observed_raw_pac)
        empirical_p = (1 + np.count_nonzero(values >= observed)) / (len(values) + 1)
        companion.assert_close(empirical_p, row.p_value, "density preview original empirical P")
        # Work in raw MI here. Scaling occurs only when rendering.
        bandwidth = float(np.std(values, ddof=1) * len(values) ** (-0.2))
        support_low = float(values.min() - 5 * bandwidth)
        support_high = float(values.max() + 5 * bandwidth)
        grid = np.linspace(support_low, support_high, GRID_SIZE)
        density = kde(grid, values, bandwidth)
        integral = float(np.trapz(density, grid))
        if not np.isclose(integral, 1.0, rtol=0, atol=1e-5):
            raise ValueError("KDE numerical integral does not equal one")
        lower = min(support_low, observed)
        upper = max(support_high, observed)
        span = upper - lower
        lower = max(0.0, lower - 0.09 * span)
        upper += 0.09 * span
        if not lower < observed < upper:
            raise ValueError("Observed MI would be outside the local plot range")
        key = (row.roi, row.phase_band, row.amp_band)
        curves[key] = (grid, density, lower, upper)
        record = row.to_dict()
        record.update({"density_method": "Gaussian KDE; Scott bandwidth; visual only",
                       "kde_bandwidth_raw_mi": bandwidth,
                       "density_integral": integral,
                       "x_min_raw_mi": lower, "x_max_raw_mi": upper,
                       "density_peak_per_raw_mi": float(density.max()),
                       "local_x_and_y_scales": True,
                       "x_statistic": "raw group-level Tort MI",
                       "density_determines_significance": False})
        rows.append(record)
    return pd.DataFrame(rows), curves


def format_ticks(ticks: np.ndarray) -> list[str]:
    gap = np.min(np.diff(ticks))
    decimals = max(0, int(np.ceil(-np.log10(gap))) + 1)
    decimals = min(decimals, 4)
    return [f"{v:.{decimals}f}" for v in ticks]


def render(companion, data, labels, curves, output: Path):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.text import Text
    from matplotlib.ticker import MaxNLocator

    fig, axes = plt.subplots(6, 4, figsize=FIGSIZE)
    fig.subplots_adjust(left=0.223, right=0.988, bottom=0.17, top=0.865,
                        wspace=0.32, hspace=0.58)
    for column, roi in enumerate(companion.ROIS):
        cell = data.loc[data.roi.eq(roi)].set_index(["phase_band", "amp_band"])
        for r, combination in enumerate(companion.COMBINATIONS):
            ax = axes[r, column]
            row = cell.loc[combination]
            grid, density, low, high = curves[(roi, *combination)]
            observed = float(row.group_observed_raw_pac) * SCALE
            x = grid * SCALE
            pdf = density / SCALE
            peak = float(pdf.max())
            ax.fill_between(x, 0, pdf, color=BLUE, alpha=0.25, linewidth=0)
            ax.plot(x, pdf, color=BLUE, lw=1.1)
            # Keep the observed marker equally legible even when its PDF is near zero.
            marker_top = peak * 1.04
            ax.vlines(observed, 0, marker_top, color=INK, linewidth=0.95, zorder=3)
            ax.plot(observed, marker_top, "v", color=INK, ms=3.0, zorder=4)
            symbol = str(row.plot_symbol)
            if symbol:
                ax.annotate(symbol, (observed, marker_top), xytext=(0, 3),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=8.6, fontweight="bold", color=INK,
                            annotation_clip=False)
            ax.set_xlim(low * SCALE, high * SCALE)
            ax.set_ylim(0, peak * 1.48)
            candidate = MaxNLocator(nbins=2, steps=[1, 2, 2.5, 5, 10]).tick_values(low * SCALE, high * SCALE)
            ticks = candidate[(candidate > low * SCALE) & (candidate < high * SCALE)]
            if len(ticks) < 2:
                ticks = np.linspace(low * SCALE, high * SCALE, 4)[1:3]
            ax.set_xticks(ticks)
            ax.set_xticklabels(format_ticks(ticks))
            ax.set_yticks([])
            ax.tick_params(axis="x", length=2.0, width=0.6, pad=1.8,
                           labelsize=7.1, color="#8E9498")
            for side in ("left", "right", "top"):
                ax.spines[side].set_visible(False)
            ax.spines["bottom"].set_color("#A9AEB2")
            ax.spines["bottom"].set_linewidth(0.65)
            if r == 0:
                ax.set_title(companion.compact_title(labels[roi]),
                             fontsize=10.2, fontweight="bold", pad=9, color=INK)
            if column == 0:
                ax.text(-0.12, 0.46, companion.ROW_LABELS[r], transform=ax.transAxes,
                        ha="right", va="center", fontsize=9.3, color=INK)
    fig.text(0.047, 0.517, "Phase / amplitude (Hz)", rotation=90,
             fontsize=10.4, fontweight="bold", ha="center", va="center")
    fig.text(0.604, 0.099, r"MI ($\times10^{-3}$)",
             fontsize=10.4, fontweight="bold", ha="center")
    legend = fig.legend(handles=[Patch(facecolor=BLUE, edgecolor=BLUE, alpha=0.35),
                                Line2D([], [], color=INK, marker="v", linewidth=0.95,
                                       markersize=3.0)],
                        labels=["Surrogate distribution", "Observed MI"],
                        loc="center", bbox_to_anchor=(0.605, 0.060), ncol=2,
                        frameon=False, fontsize=8.0, handlelength=1.4,
                        handletextpad=0.55, columnspacing=1.6)
    fig.text(0.604, 0.022, "Axis scales vary between panels", color="#626A6E",
             fontsize=7.0, ha="center")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for item in fig.findobj(Text):
        if not item.get_visible() or not item.get_text():
            continue
        box = item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            raise ValueError("Text outside figure: " + item.get_text())
    if len(axes.ravel()) != 24:
        raise ValueError("Incomplete density grid")
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-project-root", type=Path,
                        help="Preview destination; defaults to project-root")
    parser.add_argument("--scope", choices=["all_dynamic_split", "all_dynamic_split_with_terminal_C"],
                        default="all_dynamic_split")
    parser.add_argument("--replace-existing", action="store_true",
                        help="Replace only this script's density preview outputs")
    args = parser.parse_args()
    project = args.project_root.resolve()
    destination = (args.output_project_root or args.project_root).resolve()
    figure_root = destination / "figures/pac/by_lina/pac_existence"
    result_root = destination / "results/pac/by_lina/pac_existence"
    figure_path = figure_root / f"density_pac_existence_{args.scope}.png"
    table_path = result_root / f"pac_existence_density_plot_data_{args.scope}.csv"
    audit_path = figure_root / f"density_plot_audit_{args.scope}.json"
    notes_path = figure_root / f"DENSITY_PLOTTING_RULES_{args.scope}.txt"
    outputs = [figure_path, table_path, audit_path, notes_path]
    if not args.replace_existing and any(path.exists() for path in outputs):
        raise FileExistsError("Density output exists; explicit replacement required")
    companion, data, protected, labels = load_validated_data(project, args.scope)
    # Preserve the previous plot table, marker, figures and explanatory files too.
    for root in (project / "results/pac/by_lina/pac_existence",
                 project / "figures/pac/by_lina/pac_existence"):
        for path in root.iterdir():
            if path.is_file() and (path.name.startswith("forest") or path.name in
                                  ("pac_existence_forest_plot_data.csv", "_SUCCESS", "FDR_AND_PLOTTING_RULES.txt")):
                protected[str(path)] = companion.sha256(path)
    plot_data, curves = prepare_curves(companion, data, protected)
    figure_root.mkdir(parents=True, exist_ok=True)
    result_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(figure_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.titleweight": "bold",
                         "axes.labelweight": "bold", "axes.unicode_minus": False})
    render(companion, plot_data, labels, curves, figure_path)
    if any(companion.sha256(Path(p)) != digest for p, digest in protected.items()):
        raise RuntimeError("A protected original result or previous figure changed")
    plot_data.to_csv(table_path, index=False, encoding="utf-8-sig", float_format="%.17g")
    saved = pd.read_csv(table_path)
    for column in ("p_value", "q_value", "group_observed_raw_pac", "group_surrogate_raw_pac_mean",
                   "kde_bandwidth_raw_mi", "x_min_raw_mi", "x_max_raw_mi"):
        companion.assert_close(saved[column], plot_data[column], "density plot CSV " + column)
    if not saved.plot_symbol.fillna("").eq(plot_data.plot_symbol).all():
        raise ValueError("Density preview significance marks changed")
    notes = ("PAC/MI 本身显著性的密度图示例，仅复用原始保存结果，不重新计算统计检验。\n"
             f"分支：{args.scope}。all_dynamic_split仅排除终末corridor，保留终末junction；"
             "all_dynamic_split_with_terminal_C纳入终末corridor。\n"
             "蓝色：每组1000次组水平surrogate MI的Gaussian KDE；Scott带宽n^(-1/5)乘样本标准差。\n"
             "黑色竖线及顶部倒三角：真实组水平MI，不是减去零分布均值后的MI。\n"
             "每个小图分别使用局部横轴和局部密度纵轴范围；不将峰高或屏幕宽度用于跨组定量比较。\n"
             "横轴标注为MI×10^-3，仅显示缩放；保存的MI、带宽和横轴界限均为原始单位。\n"
             "显示全部surrogate样本的范围及至少5带宽尾部，所有真实值均包含在轴范围内；无断轴。\n"
             "KDE仅用于显示，不以KDE积分重新计算P。星号保持原表单侧经验P：*<.05、**<.01、***<.001。\n"
             "†保持源BH q<.05，每个原始分支aMTG6项、三个海马合并18项。\n"
             "各曲线为概率密度的平滑估计，不是参数或效应估计量的置信区间。\n"
             "原始结果、已有森林图及_SUCCESS未修改，未将此示例替换为唯一正式图。\n")
    notes_path.write_text(notes, encoding="utf-8")
    audit = {"scope": args.scope, "image": str(figure_path), "table": str(table_path),
             "figure_size_inches": FIGSIZE, "cells": len(plot_data), "surrogates_per_cell": 1000,
             "original_P_q_and_MI_preserved": True, "KDE_for_display_only": True,
             "bandwidth": "Scott: sample SD (ddof=1) times n^(-1/5)",
             "all_density_integrals_within_1e_minus_5_of_one": True,
             "local_axes": True, "all_observed_MI_inside_axes": True,
             "raw_significant": int(plot_data.raw_significant.sum()),
             "FDR_significant": int(plot_data.fdr_significant.sum()),
             "source_and_previous_output_hashes": protected, "protected_hashes_unchanged": True,
             "image_sha256": companion.sha256(figure_path),
             "table_sha256": companion.sha256(table_path), "all_text_inside_canvas": True,
             "no_new_surrogates_models_or_statistical_corrections": True,
             "ROI_labels": labels}
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"image": str(figure_path), "plot_table": str(table_path),
                      "cells": len(plot_data), "raw_significant": audit["raw_significant"],
                      "FDR_significant": audit["FDR_significant"],
                      "protected_files": len(protected)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
