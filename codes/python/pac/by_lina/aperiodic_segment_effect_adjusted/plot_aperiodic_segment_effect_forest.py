"""Compact coefficient plots of existing 1-150 Hz aperiodic exponent LMM effects.

Four ROIs on y, beta and horizontal 95% CIs on x. No model fitting or statistical changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROIS = (
    "Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
    "Hippocampus L_Posterior", "Hippocampus R_Anterior",
)
SCOPES = ("include_last", "exclude_terminal_corridor_and_junction")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value: str) -> str:
    return " ".join(value.replace("_", " ").split()).casefold()


def stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def bh_independent(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    adjusted = np.empty(len(p))
    adjusted[order] = ranked
    return adjusted


def plot_symbol(row: pd.Series) -> str:
    """Raw-P stars, plus a dagger only for an actually corrected hippocampal test."""
    dagger = "†" if row.roi in ROIS[1:] and float(row.p_fdr_bh) < 0.05 else ""
    return stars(float(row.p_value)) + dagger


def verify(data: pd.DataFrame, raw: pd.DataFrame) -> None:
    keys = ["scope", "roi"]
    expected = {(scope, roi) for scope in SCOPES for roi in ROIS}
    if len(data) != 8 or data.duplicated(keys).any():
        raise ValueError("Expected exactly eight unique segment effects")
    if set(data[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError("ROI/branch grid incomplete")
    if not data["term"].eq("segment_position").all() or not data["exponent_range"].eq("1-150 Hz").all():
        raise ValueError("Wrong term or spectral fitting range")
    merged = data.merge(raw, on=keys, suffixes=("", "_raw"), validate="one_to_one")
    for col in ["beta", "se", "t_value", "df", "ci_low", "ci_high", "p_value"]:
        if not np.allclose(merged[col], merged[col + "_raw"], rtol=1e-12, atol=1e-15):
            raise ValueError("Original fitted numbers changed: " + col)
    if not np.isfinite(data[["beta", "ci_low", "ci_high", "p_for_significance"]]).all().all():
        raise ValueError("Non-finite plotted values")
    if not (data["ci_low"].le(data["beta"]) & data["beta"].le(data["ci_high"])).all():
        raise ValueError("Estimate not within 95% CI")
    for scope in SCOPES:
        hipp = data.loc[data["scope"].eq(scope) & data["roi"].isin(ROIS[1:])]
        if len(hipp) != 3 or not hipp["correction_n_tests"].eq(3).all():
            raise ValueError("BH family must have three hippocampal tests")
        q = bh_independent(hipp["p_value"].to_numpy(float))
        if not np.allclose(q, hipp["p_fdr_bh"], rtol=1e-12, atol=1e-15):
            raise ValueError("BH verification failed")
        if not np.allclose(q, hipp["p_for_significance"], rtol=1e-12, atol=1e-15):
            raise ValueError("Hippocampus symbols must use corrected P")
        amt = data.loc[data["scope"].eq(scope) & data["roi"].eq(ROIS[0])]
        if len(amt) != 1 or not amt["p_fdr_bh"].isna().all():
            raise ValueError("aMTG must remain uncorrected")
        if not np.allclose(amt["p_for_significance"], amt["p_value"], rtol=1e-12, atol=1e-15):
            raise ValueError("aMTG symbols must use raw P")


def render(cell: pd.DataFrame, labels: dict[str, str], output: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    fig, ax = plt.subplots(figsize=(4.9, 2.55))
    fig.subplots_adjust(left=0.395, right=0.965, bottom=0.235, top=0.94)
    colour = "#398EAD"
    ax.axvline(0, color="#AEB5BA", lw=0.85, zorder=1)
    ax.grid(axis="x", color="#F0F1F2", lw=0.6, zorder=0)
    for y, roi in enumerate(ROIS):
        row = cell.loc[cell["roi"].eq(roi)].iloc[0]
        estimate, lo, hi = np.array([row.beta, row.ci_low, row.ci_high])
        ax.errorbar(estimate, y, xerr=np.array([[estimate - lo], [hi - estimate]]),
                    fmt="none", ecolor=colour, elinewidth=1.35, capsize=2.3,
                    capthick=1.05, alpha=0.56, zorder=2)
        ax.plot(estimate, y, "o", ms=5.4, mfc=colour, mec="white", mew=0.45, zorder=3)
        symbol = plot_symbol(row)
        if symbol:
            ax.annotate(symbol, (hi, y), xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=9.5, fontweight="bold", color="#242424")
    ax.set_xlim(-0.002, 0.0044)
    ax.set_ylim(3.45, -0.45)
    ax.set_yticks(range(4))
    ax.set_yticklabels([labels[roi] for roi in ROIS], rotation=0, ha="right",
                       fontsize=9, fontweight="bold")
    ax.set_xticks(np.arange(-0.001, 0.0041, 0.001))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, -3))
    formatter.set_useOffset(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlabel(r"$\beta$", fontsize=12, labelpad=4)
    ax.xaxis.get_offset_text().set_fontsize(9)
    ax.tick_params(axis="x", length=2.5, pad=3, width=0.8, color="#777777")
    ax.tick_params(axis="y", length=0, pad=7, width=0.8, color="#777777")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#777777")
        ax.spines[spine].set_linewidth(0.7)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in [ax.xaxis.label, ax.xaxis.get_offset_text(), *ax.texts,
                 *ax.get_xticklabels(), *ax.get_yticklabels()]:
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                raise ValueError("Text outside canvas: " + text.get_text())
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--replace-existing", action="store_true",
                        help="Re-render only this script's two explicitly named PNGs")
    args = parser.parse_args()
    project = args.project_root.resolve()
    result_root = project / "results/pac/by_lina/aperiodic_segment_effect_adjusted"
    output_root = args.output_root or project / "figures/pac/by_lina/aperiodic_segment_effect_adjusted"
    output_root.mkdir(parents=True, exist_ok=True)
    images = [output_root / scope / f"forest_aperiodic_segment_effect_{scope}.png" for scope in SCOPES]
    for path in images:
        if path.exists() and not args.replace_existing:
            raise FileExistsError("Will not overwrite existing plot: " + str(path))
        path.parent.mkdir(parents=True, exist_ok=True)
    if not (result_root / "_FDR_SUCCESS").exists():
        raise ValueError("Validated FDR post-processing is required before plotting")
    os.environ.setdefault("MPLCONFIGDIR", str(output_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {entry.name for entry in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.linewidth": 0.8, "axes.labelweight": "bold",
                         "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 9.5})
    paths = [result_root / "all_segment_effects_fdr_bh.csv", result_root / "all_segment_effects.csv",
             project / "config/project_config.json"]
    before = {str(path): sha256(path) for path in paths}
    data = pd.read_csv(paths[0])
    verify(data, pd.read_csv(paths[1]))
    cfg = json.loads(paths[2].read_text(encoding="utf-8"))
    lookup = {normalize(row["full_name"]): row["abbreviation"] for row in cfg["roi_abbreviations"]}
    labels = {roi: lookup[normalize(roi)] for roi in ROIS}
    symbols = []
    for scope, output in zip(SCOPES, images):
        cell = data.loc[data["scope"].eq(scope)].copy()
        render(cell, labels, output)
        for roi in ROIS:
            row = cell.loc[cell["roi"].eq(roi)].iloc[0]
            symbols.append({"scope": scope, "roi": roi, "roi_label": labels[roi],
                            "beta": float(row.beta), "ci_low": float(row.ci_low),
                            "ci_high": float(row.ci_high), "p_for_significance": float(row.p_for_significance),
                            "p_value_raw": float(row.p_value),
                            "p_fdr_bh": None if pd.isna(row.p_fdr_bh) else float(row.p_fdr_bh),
                            "raw_stars": stars(float(row.p_value)),
                            "fdr_dagger": "†" if row.roi in ROIS[1:] and float(row.p_fdr_bh) < 0.05 else "",
                            "significance_basis": "raw P stars; dagger denotes hippocampal BH q < 0.05",
                            "plot_symbol": plot_symbol(row)})
    if any(sha256(Path(path)) != expected for path, expected in before.items()):
        raise RuntimeError("Read-only plot source changed")
    audit = {"figure_files": [str(path) for path in images], "figure_size_inches": [4.9, 2.55],
             "original_sources_sha256": before, "input_hashes_unchanged": True,
             "BH_independently_verified": True, "model_fitting_in_plotter": False,
             "shared_xlim_raw_beta": [-0.002, 0.0044], "display_scale": "scientific axis offset 10^-3; beta unchanged",
             "x_axis": "beta and original 95% CI", "y_axis": "four ROIs",
             "ROI_order_top_to_bottom": list(ROIS), "ROI_label_rotation_degrees": 0,
             "symbol_offset_points": [4, 0], "symbols_right_of_upper_CI": True,
             "all_visible_text_checked_inside_canvas": True, "symbols": symbols}
    (output_root / "plot_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    notes = (
        "非周期性谱指数的 segment 序号效应：绘图规则\n"
        "仅呈现 1–150 Hz IRASA 指数，模型包含分类 condition 固定效应和被试随机截距。\n"
        "纵轴为四脑区，标准缩写水平标注；横轴为原模型beta，横线为原模型95%置信区间。\n"
        "横轴保留原始beta，科学计数法10^-3标注在轴右端；不改变模型估计值。\n"
        "标记在每条误差线右端附近：星号使用原始P，* P<0.05；** P<0.01；*** P<0.001。\n"
        "†表示三个海马共同BH-FDR（每分支3检验）后q<0.05；未通过则不加†。\n"
        "aMTG每分支仅一个检验，不校正、只标原始P星号，不加†；不加†不代表aMTG未通过校正。\n"
        "include_last 和排除末 corridor 两半及末 junction 的分支分别保存、分别校正。\n"
        "四脑区采用 config/project_config.json 的标准缩写；两分支共用横轴范围便于比较。\n"
        "PNG 不包含底部显著性说明，无 PDF/SVG；不重拟合、不修改结果数据。\n"
        "图中横向误差线是未作多重比较调整的95%区间。\n"
    )
    (output_root / "FDR_AND_PLOTTING_RULES.txt").write_text(notes, encoding="utf-8")
    print(json.dumps({"images": [str(path) for path in images], "star_policy": "raw P stars; hippocampus BH q dagger", "symbols": symbols}))


if __name__ == "__main__":
    main()
