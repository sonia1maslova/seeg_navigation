"""Compact Power LMM forest matrices from the corrected result table.

No model fitting. Verify R BH-FDR independently before displaying its symbols.
Keep original model estimates/intervals and shared project ROI abbreviations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

CORE_ROIS = (
    "Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
    "Hippocampus L_Posterior", "Hippocampus R_Anterior",
)
SUPPLEMENTARY_ROIS = (
    "Temporal_Mid L_Posterior", "Temporal_Sup L", "SupraMarginal L", "Precentral L",
)
BANDS = ("f0p5_4", "f4_12", "f12_30", "f30_70", "f70_150")
BAND_LABELS = ("0.5–4", "4–12", "12–30", "30–70", "70–150")
VARIANTS = ("include_last", "exclude_last")
PREDICTORS = {
    "Angle2orientation": ("Egocentric_bearing", "Egocentric bearing effect, β (dB/rad)", ("corridor",)),
    "Angle2Xaxis": ("Allocentric_bearing", "Allocentric bearing effect, β (dB/rad)", ("corridor", "junction", "all_segments")),
    "SegmentNumber": ("Segment_number", "Segment number effect, β (dB/segment)", ("corridor", "junction", "all_segments")),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value: str) -> str:
    return " ".join(value.replace("_", " ").split()).casefold()


def raw_stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def bh_independent(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    q = np.empty(len(p))
    q[order] = ranked
    return q


def verify_results(data: pd.DataFrame) -> None:
    keys = ["roi", "scope", "variant", "predictor", "band_name"]
    if data.duplicated(keys).any():
        raise ValueError("Duplicate model hypothesis keys")
    required = ["estimate", "conf_low_95", "conf_high_95", "p_value_raw", "p_fdr"]
    if not np.isfinite(data[required].to_numpy(float)).all():
        raise ValueError("Non-finite estimate, interval or p/q")
    if not data["fit_status"].eq("ok").all():
        raise ValueError("Not all fitted models passed the fitting stage")
    if not (data["conf_low_95"].le(data["estimate"]) & data["estimate"].le(data["conf_high_95"])).all():
        raise ValueError("Coefficient outside confidence interval")
    for _, cell in data.groupby(["predictor", "scope", "variant", "fdr_family"], sort=False):
        planned = cell["fdr_n_planned"].unique()
        if len(planned) != 1 or len(cell) != planned[0]:
            raise ValueError("Incomplete planned FDR family")
        expected = 15 if set(cell["roi"]) == set(CORE_ROIS[1:]) else 5
        if planned[0] != expected:
            raise ValueError("FDR family not 3 hippocampal regions x 5 or one other ROI x 5")
        q = bh_independent(cell["p_value_raw"].to_numpy(float))
        if not np.allclose(q, cell["p_fdr"], rtol=1e-12, atol=1e-14):
            raise ValueError("R and independent BH-FDR calculations disagree")
    expected_symbols = data["p_value_raw"].map(raw_stars) + np.where(data["p_fdr"] < 0.05, "†", "")
    if not data["plot_symbol"].fillna("").eq(expected_symbols).all():
        raise ValueError("Saved significance symbols disagree with p/q")


def select_cell(data: pd.DataFrame, rois: tuple[str, ...], predictor: str, scope: str, variant: str) -> pd.DataFrame:
    cell = data.loc[
        data["roi"].isin(rois) & data["predictor"].eq(predictor)
        & data["scope"].eq(scope) & data["variant"].eq(variant)
    ].copy()
    grid = {(roi, band) for roi in rois for band in BANDS}
    if set(cell[["roi", "band_name"]].itertuples(index=False, name=None)) != grid or len(cell) != 20:
        raise ValueError(f"Incomplete 4 ROI x 5 band grid: {predictor}/{scope}/{variant}")
    cell["roi_order"] = cell["roi"].map({roi: i for i, roi in enumerate(rois)})
    cell["band_order"] = cell["band_name"].map({band: i for i, band in enumerate(BANDS)})
    return cell.sort_values(["roi_order", "band_order"])


def shared_extent(data: pd.DataFrame) -> float:
    largest = float(np.abs(data[["conf_low_95", "conf_high_95"]].to_numpy()).max())
    if largest == 0:
        return 0.1
    padded = largest * 1.16
    scale = 10 ** math.floor(math.log10(padded))
    for multiplier in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if multiplier * scale >= padded:
            return multiplier * scale
    raise AssertionError("No axis extent")


def compact_title(label: str) -> str:
    if " Hippocampus " in label:
        return label.replace(" Hippocampus ", "\nHippocampus ")
    if label in ("Ant Temp Mid L", "Post Temp Mid L"):
        return label.replace(" Mid L", "\nMid L")
    return label


def render(cell: pd.DataFrame, rois: tuple[str, ...], labels: dict[str, str], xlabel: str, extent: float, output: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 4, figsize=(6.05, 2.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15, right=0.986, bottom=0.235, top=0.805, wspace=0.21)
    colour = "#398EAD"
    tick = extent / 2
    for ax, roi in zip(axes, rois):
        table = cell.loc[cell["roi"].eq(roi)].sort_values("band_order")
        ax.axvline(0, color="#BFC4C7", lw=0.85, zorder=0)
        for y, row in enumerate(table.itertuples(index=False)):
            ax.axhline(y, color="#F0F1F2", lw=0.6, zorder=0)
            error = np.array([[row.estimate - row.conf_low_95], [row.conf_high_95 - row.estimate]])
            ax.errorbar(row.estimate, y, xerr=error, fmt="none", ecolor=colour,
                        elinewidth=1.35, capsize=2.3, capthick=1.05, alpha=0.56, zorder=2)
            ax.plot(row.estimate, y, "o", ms=4.5, mfc=colour, mec="white", mew=0.4, zorder=3)
            if row.plot_symbol:
                ax.annotate(row.plot_symbol, (row.estimate, y), xytext=(0, 3.1),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=9, fontweight="bold", color="#242424")
        ax.set_title(compact_title(labels[roi]), fontsize=9.6, pad=8, color="#222222")
        ax.set_xlim(-extent, extent)
        ax.set_ylim(4.5, -0.56)
        ax.set_xticks([-tick, 0, tick])
        ax.set_xticklabels([f"−{tick:g}", "0", f"{tick:g}"])
        ax.set_yticks(np.arange(5))
        ax.tick_params(axis="x", length=2.5, pad=2.5, width=0.8, color="#777777")
        ax.tick_params(axis="y", length=0, pad=7)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#777777")
        ax.spines["bottom"].set_linewidth(0.7)
    axes[0].set_yticklabels(BAND_LABELS)
    axes[0].set_ylabel("Frequency (Hz)", fontsize=10, labelpad=6)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.text(0.56, 0.075, xlabel, ha="center", fontsize=10, fontweight="bold")
    # Check visible labels and symbols before saving every formal image.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        texts = [ax.title, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]
        for text in texts:
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                    raise ValueError("Text outside figure canvas: " + text.get_text())
    for text in fig.texts:
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            raise ValueError("Global label outside figure canvas")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.environ.get("SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")))
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output_root = args.output_root or project / "figures/power/encoding_lmm/forest_plots"
    output_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    available = {entry.name for entry in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.linewidth": 0.8, "axes.labelweight": "bold",
                         "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 9.5})
    model_path = project / "results/power/encoding_lmm/summaries/all_target_predictor_effects_fdr_bh.csv"
    config_path = project / "config/project_config.json"
    before = {str(path): sha256(path) for path in (model_path, config_path)}
    data = pd.read_csv(model_path).fillna({"plot_symbol": "", "raw_stars": "", "fdr_dagger": ""})
    verify_results(data)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    lookup = {normalize(row["full_name"]): row["abbreviation"] for row in config["roi_abbreviations"]}
    labels = {roi: lookup[normalize(roi)] for roi in CORE_ROIS + SUPPLEMENTARY_ROIS}
    plots = []
    for predictor, (slug, xlabel, scopes) in PREDICTORS.items():
        for scope in scopes:
            # The same predictor can use different indexing/ranges across
            # scopes. Share axes across ROIs and last-segment branches within
            # each scope, without flattening narrow CIs in unrelated scopes.
            extent = shared_extent(data.loc[data["roi"].isin(CORE_ROIS) & data["predictor"].eq(predictor) & data["scope"].eq(scope)])
            for variant in VARIANTS:
                cell = select_cell(data, CORE_ROIS, predictor, scope, variant)
                cell["roi_display"] = cell["roi"].map(labels)
                cell["figure_group"] = "core_four_rois"
                output = output_root / "core_four_rois" / slug / scope / variant / f"forest_{slug}_{scope}_{variant}.png"
                render(cell, CORE_ROIS, labels, xlabel, extent, output)
                cell.to_csv(output.with_name(output.stem + "_values.csv"), index=False, encoding="utf-8-sig")
                plots.append({"file": str(output), "group": "core_four_rois", "predictor": predictor,
                              "scope": scope, "variant": variant, "rows": len(cell),
                              "raw_significant": int(cell["raw_sig"].sum()), "fdr_significant": int(cell["fdr_sig"].sum()),
                              "axis_extent": extent})
    predictor = "Angle2orientation"
    slug, xlabel, _ = PREDICTORS[predictor]
    extent = shared_extent(data.loc[data["roi"].isin(SUPPLEMENTARY_ROIS) & data["predictor"].eq(predictor)])
    for variant in VARIANTS:
        cell = select_cell(data, SUPPLEMENTARY_ROIS, predictor, "corridor", variant)
        cell["roi_display"] = cell["roi"].map(labels)
        cell["figure_group"] = "supplementary_cortical_rois"
        output = output_root / "supplementary_cortical_rois" / slug / "corridor" / variant / f"forest_{slug}_corridor_{variant}.png"
        render(cell, SUPPLEMENTARY_ROIS, labels, xlabel, extent, output)
        cell.to_csv(output.with_name(output.stem + "_values.csv"), index=False, encoding="utf-8-sig")
        plots.append({"file": str(output), "group": "supplementary_cortical_rois", "predictor": predictor,
                      "scope": "corridor", "variant": variant, "rows": len(cell),
                      "raw_significant": int(cell["raw_sig"].sum()), "fdr_significant": int(cell["fdr_sig"].sum()),
                      "axis_extent": extent})
    pd.DataFrame(plots).to_csv(output_root / "figure_manifest.csv", index=False, encoding="utf-8-sig")
    audit = {"figures": len(plots), "all_plots_have_20_model_rows": all(p["rows"] == 20 for p in plots),
             "BH_independently_verified": True, "ROI_abbreviations": labels, "input_sha256": before,
             "original_start_to_position_bearing": True, "model_fitting_in_plotter": False,
             "figure_size_inches": [6.05, 2.6], "symbol_offset_points": 3.1,
             "axis_scale_shared_across_ROIs_and_variants_within_predictor_scope_and_figure_group": True,
             "all_visible_text_checked_inside_canvas": True}
    if any(sha256(Path(path)) != expected for path, expected in before.items()):
        raise RuntimeError("Read-only plot source changed")
    (output_root / "plot_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    notes = (
        "Power 线性混合模型：绘图与 FDR 校正规则（2026-10-10）\n"
        "每个预测特征、segment 范围和 include_last/exclude_last 分支分别进行 BH-FDR。\n"
        "三个海马子区的五个频段联合校正，共 15 个检验；其余每个脑区各自五频段校正，共 5 个检验。\n"
        "角度沿用起点指向当前位置的编码；STG 使用整区，不区分前后。\n"
        "点表示回归系数，横线表示原模型的 95% 置信区间；绘图不改变系数、区间或原始 P 值。\n"
        "星号表示原始 P：* <0.05、** <0.01、*** <0.001；† 表示 BH-FDR q <0.05。\n"
        "图内不打印底部显著性说明；现有 Segment number 热图已同步采用新校正规则。\n\n"
        "Power encoding compact forest plots\n"
        "Points: model coefficient; horizontal intervals: original model 95% confidence interval.\n"
        "Original bearing encoding: start-to-position. STG: whole Temporal_Sup L.\n"
        "Stars: raw p, * <.05, ** <.01, *** <.001. Dagger: BH q <.05.\n"
        "BH is separate per predictor x scope x last-segment branch.\n"
        "Three hippocampal subregions pooled over five bands: 15 tests. Every other ROI: 5 tests.\n"
        "No bottom significance footnote is printed in PNGs.\n"
        "Raw model estimates and p values are unchanged by FDR.\n"
        "Retained Segment-number heatmaps use the same updated correction families.\n"
    )
    (output_root / "FDR_AND_PLOTTING_RULES.txt").write_text(notes, encoding="utf-8")
    print(json.dumps({"output_root": str(output_root), "figures": len(plots), "BH_verified": True}))


if __name__ == "__main__":
    main()
