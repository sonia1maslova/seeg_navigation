"""Lina whole-trial PAC-behavior forest plot, with separate 6/18-test BH families.

Read-only plotting: no signal extraction, fitting, or source-results edits.
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

ROIS = ("Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
        "Hippocampus L_Posterior", "Hippocampus R_Anterior")
PHASES = ("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
AMPLITUDES = ("low_gamma_30_70Hz", "high_gamma_70_150Hz")
COMBINATIONS = tuple((phase, amp) for amp in AMPLITUDES for phase in PHASES)
ROW_LABELS = ("0.5–4 / 30–70", "4–12 / 30–70", "12–30 / 30–70",
              "0.5–4 / 70–150", "4–12 / 70–150", "12–30 / 70–150")
DISPLAY_DIVISOR = 10000.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value: str) -> str:
    return " ".join(value.replace("_", " ").split()).casefold()


def raw_stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(ranked, 1.0)
    return result


def verify(data: pd.DataFrame, source: pd.DataFrame) -> None:
    keys = ["roi", "phase_band", "amp_band"]
    expected = {(roi, phase, amp) for roi in ROIS for phase, amp in COMBINATIONS}
    if len(data) != 24 or data.duplicated(keys).any() or set(data[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError("Expected a complete four-ROI x six-combination grid")
    if not data["scope"].eq("whole_trial").all() or not data["analysis_family"].eq("04C1_pac_behavior").all():
        raise ValueError("Only whole-trial 04C1 behavior effects may be plotted")
    if not data["target_term"].eq("raw_PAC").all():
        raise ValueError("Wrong fixed-effect term")
    if not data["fit_status"].eq("ok").all() or not data["is_singular"].eq(0).all():
        raise ValueError("Model numerical checks failed")
    if not np.isfinite(data[["beta", "ci_low", "ci_high", "p_value", "q_value"]]).all().all():
        raise ValueError("Non-finite plotted statistic")
    if not (data["ci_low"].le(data["beta"]) & data["beta"].le(data["ci_high"])).all():
        raise ValueError("Estimate lies outside confidence interval")
    original = source.loc[source["analysis_family"].eq("04C1_pac_behavior")]
    if len(original) != 24:
        raise ValueError("Original behavior family must have 24 effects")
    joined = data.merge(original, on=keys, validate="one_to_one", suffixes=("", "_original"))
    for col in ("beta", "se", "t_value", "df", "ci_low", "ci_high", "p_value", "n_obs", "n_subjects"):
        if not np.allclose(joined[col], joined[col + "_original"], rtol=1e-12, atol=1e-15):
            raise ValueError("Original fitted number changed: " + col)
    if not np.allclose(joined["q_value_lina_24_tests"], joined["q_value_original"], rtol=1e-12, atol=1e-15):
        raise ValueError("Original correction provenance not retained")
    if data["fdr_family_id"].nunique() != 2:
        raise ValueError("Exactly two BH families required")
    for family, cell in data.groupby("fdr_family_id", sort=False):
        roi_set = set(cell["roi"])
        if roi_set == {ROIS[0]}:
            expected_n = 6
        elif roi_set == set(ROIS[1:]):
            expected_n = 18
        else:
            raise ValueError("Incorrect family ROI composition")
        if len(cell) != expected_n or not cell["m_planned"].eq(expected_n).all():
            raise ValueError("Incorrect six/eighteen-test family: " + family)
        if not np.allclose(bh(cell["p_value"].to_numpy(float)), cell["q_value"], rtol=1e-12, atol=1e-15):
            raise ValueError("Independent BH verification failed")
    expected_symbols = data["p_value"].map(raw_stars) + np.where(data["q_value"] < 0.05, "†", "")
    if not data["plot_symbol"].fillna("").eq(expected_symbols).all():
        raise ValueError("Symbols inconsistent with raw P or new BH q")


def compact_title(label: str) -> str:
    if " Hippocampus " in label:
        return label.replace(" Hippocampus ", "\nHippocampus ")
    if label == "Ant Temp Mid L":
        return label.replace(" Mid L", "\nMid L")
    return label


def axis_extent(data: pd.DataFrame) -> float:
    padded = float(np.abs(data[["ci_low", "ci_high"]].to_numpy()).max()) / DISPLAY_DIVISOR * 1.16
    if padded <= 0:
        return 1.0
    scale = 10 ** math.floor(math.log10(padded))
    for multiplier in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if multiplier * scale >= padded:
            return multiplier * scale
    raise AssertionError("No axis extent")


def render(data: pd.DataFrame, labels: dict[str, str], output: Path) -> float:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 4, figsize=(6.55, 3.15), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.245, right=0.988, bottom=0.255, top=0.825, wspace=0.21)
    colour = "#398EAD"
    extent = axis_extent(data)
    for ax, roi in zip(axes, ROIS):
        cell = data.loc[data["roi"].eq(roi)].set_index(["phase_band", "amp_band"])
        ax.axvline(0, color="#BFC4C7", lw=0.85, zorder=0)
        for y, combination in enumerate(COMBINATIONS):
            row = cell.loc[combination]
            estimate, low, high = np.array([row.beta, row.ci_low, row.ci_high]) / DISPLAY_DIVISOR
            ax.axhline(y, color="#F0F1F2", lw=0.6, zorder=0)
            ax.errorbar(estimate, y, xerr=np.array([[estimate - low], [high - estimate]]),
                        fmt="none", ecolor=colour, elinewidth=1.35, capsize=2.3,
                        capthick=1.05, alpha=0.56, zorder=2)
            ax.plot(estimate, y, "o", ms=4.5, mfc=colour, mec="white", mew=0.4, zorder=3)
            if row.plot_symbol:
                ax.annotate(row.plot_symbol, (estimate, y), xytext=(0, 3.1),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=9, fontweight="bold", color="#242424")
        ax.set_title(compact_title(labels[roi]), fontsize=9.6, pad=8, color="#222222")
        ax.set_xlim(-extent, extent)
        ax.set_ylim(5.5, -0.56)
        ax.set_xticks([-extent / 2, 0, extent / 2])
        ax.set_xticklabels([f"−{extent / 2:g}", "0", f"{extent / 2:g}"])
        ax.set_yticks(np.arange(6))
        ax.tick_params(axis="x", length=2.5, pad=2.5, width=0.8, color="#777777")
        ax.tick_params(axis="y", length=0, pad=7)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#777777")
        ax.spines["bottom"].set_linewidth(0.7)
    axes[0].set_yticklabels(ROW_LABELS)
    axes[0].set_ylabel("Phase / amplitude (Hz)", fontsize=10, labelpad=6)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.text(0.605, 0.13, r"PAC effect on pointing error, $\beta$ ($\times10^4$)",
             ha="center", fontsize=10, fontweight="bold")
    fig.text(0.605, 0.045, "* p < 0.05    ** p < 0.01    *** p < 0.001    † BH-FDR q < 0.05",
             ha="center", fontsize=8.2, color="#5D6266")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        for text in [ax.title, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]:
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                    raise ValueError("Text outside canvas: " + text.get_text())
    for text in fig.texts:
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            raise ValueError("Figure label outside canvas")
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return extent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    result_root = project / "results/pac/by_lina/pac_behavior"
    output_root = project / "figures/pac/by_lina/pac_behavior"
    output_root.mkdir(parents=True, exist_ok=True)
    output = output_root / "forest_pac_behavior.png"
    if output.exists() and not args.replace_existing:
        raise FileExistsError("Named forest plot exists; explicit replacement required")
    if not (result_root / "_SUCCESS").exists():
        raise ValueError("Completed six/eighteen-test correction required")
    os.environ.setdefault("MPLCONFIGDIR", str(output_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.linewidth": 0.8, "axes.labelweight": "bold",
                         "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 9.5})
    paths = [result_root / "pac_behavior_effects_fdr_bh.csv",
             project / "results/pac/by_lina/formal_raw_pac/04_amtg_followup_and_behavior/model_coefficients.csv",
             project / "config/project_config.json"]
    before = {str(path): sha256(path) for path in paths}
    data = pd.read_csv(paths[0]).fillna({"raw_stars": "", "fdr_dagger": "", "plot_symbol": ""})
    verify(data, pd.read_csv(paths[1]))
    cfg = json.loads(paths[2].read_text(encoding="utf-8"))
    lookup = {normalize(row["full_name"]): row["abbreviation"] for row in cfg["roi_abbreviations"]}
    labels = {roi: lookup[normalize(roi)] for roi in ROIS}
    if not data["roi_label"].eq(data["roi"].map(labels)).all():
        raise ValueError("ROI labels disagree with shared config")
    extent = render(data, labels, output)
    if any(sha256(Path(path)) != expected for path, expected in before.items()):
        raise RuntimeError("Read-only source changed")
    audit = {"figure_file": str(output), "figure_size_inches": [6.55, 3.15],
             "row_count": 24, "ROI_columns": list(ROIS), "row_combinations": list(COMBINATIONS),
             "ROI_labels": labels, "BH_independently_verified": True,
             "family_sizes": {"aMTG": 6, "hippocampus": 18},
             "original_sources_sha256": before, "input_hashes_unchanged": True,
             "model_fitting_in_plotter": False, "display_divisor": DISPLAY_DIVISOR,
             "shared_xlim_display_units": [-extent, extent], "all_CIs_fully_shown": True,
             "symbol_offset_points": [0, 3.1], "all_visible_text_checked_inside_canvas": True,
             "symbols": data[["roi", "phase_band", "amp_band", "p_value", "q_value", "plot_symbol"]].to_dict(orient="records")}
    (output_root / "forest_plot_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    notes = (
        "PAC–指向误差关联森林图\n"
        "复用李娜04C1整试次模型，不重算、不重拟合；非segment趋势或ROI交互检验。\n"
        "aMTG六组合单独BH-FDR；三个海马18组合共同BH-FDR，不按振幅组拆分校正。\n"
        "星号为原始P：*<0.05、**<0.01、***<0.001；†表示本次BH q<0.05。\n"
        "点为beta，横线为原模型未作多重比较调整的95%CI；四脑区横轴尺度相同。\n"
        "仅绘图时将beta和CI除以10^4，横轴单位×10^4；原始表未缩放、未修改。\n"
        "六行由上到下为30–70Hz振幅的三相位组合，再为70–150Hz振幅的三相位组合。\n"
        "负beta表示更强PAC与更小pointing_error相关；不表示因果效应。\n"
        "PNG采用项目统一脑区缩写。已有PAC行为热图未删除、未覆盖。\n"
    )
    (output_root / "FOREST_FDR_AND_PLOTTING_RULES.txt").write_text(notes, encoding="utf-8")
    print(json.dumps({"figure": str(output), "raw_significant": int(data["raw_significant"].sum()),
                      "fdr_significant": int(data["fdr_significant"].sum()), "BH_verified": True}))


if __name__ == "__main__":
    main()
