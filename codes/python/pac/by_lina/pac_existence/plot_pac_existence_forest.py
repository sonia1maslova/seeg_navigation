"""Plot Lina's saved pooled MI existence results without rerunning inference.

Blue points are observed MI minus the mean surrogate MI. Gray horizontal
segments show the centered surrogate distribution's central 95% range.
They are NOT confidence intervals for an estimated effect. Original P and
BH q values are preserved, with independent read-only validation.
"""
from __future__ import annotations

import argparse
import hashlib
import io
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
COMBINATIONS = tuple((p, a) for a in AMPLITUDES for p in PHASES)
ROW_LABELS = ("0.5–4 / 30–70", "4–12 / 30–70", "12–30 / 30–70",
              "0.5–4 / 70–150", "4–12 / 70–150", "12–30 / 70–150")
SCOPES = ("all_dynamic_split", "all_dynamic_split_with_terminal_C")
KEYS = ["scope", "roi", "phase_band", "amp_band"]
SCALE = 1000.0
N_SURROGATES = 1000


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
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def assert_close(actual, expected, label: str) -> None:
    if not np.allclose(actual, expected, rtol=1e-11, atol=1e-15):
        raise ValueError("Source validation failed: " + label)


def read_csv(path: Path, protected: dict[str, str]) -> pd.DataFrame:
    content = path.read_bytes()
    protected[str(path)] = hashlib.sha256(content).hexdigest()
    return pd.read_csv(io.BytesIO(content))


def load_data(project: Path) -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    protected: dict[str, str] = {}
    source_root = project / "results/pac/by_lina/formal_raw_pac"
    stages = ((source_root / "04_amtg_followup_and_behavior/04A_amtg_pac_existence", (ROIS[0],)),
              (source_root / "01_hpc_pac_existence", ROIS[1:]))
    rows = []
    for stage, stage_rois in stages:
        manifest_path = stage / "stage_manifest.json"
        protected[str(manifest_path)] = sha256(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["n_surrogates"] != N_SURROGATES or manifest["status"] != "completed":
            raise ValueError("Expected completed original 1000-surrogate analysis")
        source = read_csv(stage / "existence_results.csv", protected)
        source = source.loc[source.analysis_unit.eq("pooled_segment") & source.scope.isin(SCOPES)].copy()
        expected = {(s, r, p, a) for s in SCOPES for r in stage_rois for p, a in COMBINATIONS}
        if source.duplicated(KEYS).any() or set(source[KEYS].itertuples(index=False, name=None)) != expected:
            raise ValueError("Missing or duplicate pooled MI existence cell")
        if not source.existence_status.eq("ok").all():
            raise ValueError("Unusable original result")
        for scope, family in source.groupby("scope", sort=False):
            family_size = len(stage_rois) * len(COMBINATIONS)
            if family.fdr_family_id.nunique() != 1 or len(family) != family_size:
                raise ValueError("Incorrect source BH family")
            if not family.m_planned.eq(family_size).all() or not family.m_valid.eq(family_size).all():
                raise ValueError("Incorrect original family size")
            assert_close(family.q_value, bh(family.p_value.to_numpy(float)), "original BH q")
        for phase, amp in COMBINATIONS:
            combo = f"{phase}_amp_{amp}"
            combo_root = stage / "01_pac_features" / combo
            null = read_csv(combo_root / "formal_group_pooled_surrogate_distribution.csv", protected)
            subjects = read_csv(combo_root / "formal_subject_pooled_existence.csv", protected)
            if len(null) != len(stage_rois) * len(SCOPES) * N_SURROGATES:
                raise ValueError("Incomplete surrogate source")
            for scope in SCOPES:
                for roi in stage_rois:
                    selected = source.loc[source.scope.eq(scope) & source.roi.eq(roi)
                                          & source.phase_band.eq(phase) & source.amp_band.eq(amp)]
                    if len(selected) != 1:
                        raise ValueError("Ambiguous source result")
                    original = selected.iloc[0]
                    subset = null.loc[null.trajectory_scope.eq(scope) & null.location1.eq(roi)]
                    if len(subset) != N_SURROGATES or set(subset.surrogate_index) != set(range(1, N_SURROGATES + 1)):
                        raise ValueError("Missing or repeated surrogate index")
                    if not subset.phase_band.eq(phase).all() or not subset.amp_band.eq(amp).all():
                        raise ValueError("Surrogate combination mismatch")
                    values = subset.group_surrogate_raw_pac.to_numpy(float)
                    if not np.isfinite(values).all():
                        raise ValueError("Non-finite surrogate values")
                    observed = float(original.group_observed_raw_pac)
                    null_mean = float(np.mean(values))
                    null_sd = float(np.std(values, ddof=0))
                    p = (1 + np.count_nonzero(values >= observed)) / (len(values) + 1)
                    assert_close(p, original.p_value, "right-tail empirical P")
                    assert_close(null_mean, original.group_surrogate_raw_pac_mean, "surrogate mean")
                    assert_close(null_sd, original.group_surrogate_raw_pac_sd, "surrogate SD")
                    subject = subjects.loc[subjects.trajectory_scope.eq(scope) & subjects.location1.eq(roi)]
                    if subject.sub_id.duplicated().any() or len(subject) != int(original.n_subjects_usable):
                        raise ValueError("Subject composition mismatch")
                    assert_close(subject.pac_tort_mi.mean(), observed, "subject-balanced observed MI")
                    delta = observed - null_mean
                    assert_close(delta, original.group_raw_pac_delta_vs_surrogate_mean, "observed MI excess")
                    low, high = np.quantile(values, [0.025, 0.975], method="linear")
                    row = original.to_dict()
                    row.update({"null_q025_raw_mi": float(low), "null_q975_raw_mi": float(high),
                                "null_q95_raw_mi": float(np.quantile(values, 0.95, method="linear")),
                                "observed_delta_mi": delta,
                                "null_q025_centered": float(low - null_mean),
                                "null_q975_centered": float(high - null_mean),
                                "n_surrogates": len(values), "source_stage": str(stage),
                                "source_null_file": str(combo_root / "formal_group_pooled_surrogate_distribution.csv"),
                                "range_definition": "central_95_percent_of_centered_surrogate_distribution_not_estimator_CI",
                                "raw_stars": raw_stars(float(original.p_value)),
                                "fdr_dagger": "†" if original.q_value < 0.05 else "",
                                "raw_significant": bool(original.p_value < 0.05),
                                "fdr_significant": bool(original.q_value < 0.05)})
                    row["plot_symbol"] = row["raw_stars"] + row["fdr_dagger"]
                    rows.append(row)
    data = pd.DataFrame(rows)
    if len(data) != 48 or data.duplicated(KEYS).any():
        raise ValueError("Expected exactly 48 original pooled results")
    cfg_path = project / "config/project_config.json"
    content = cfg_path.read_bytes()
    protected[str(cfg_path)] = hashlib.sha256(content).hexdigest()
    config = json.loads(content.decode("utf-8"))
    lookup = {normalize(row["full_name"]): row["abbreviation"] for row in config["roi_abbreviations"]}
    labels = {roi: lookup[normalize(roi)] for roi in ROIS}
    data["roi_label"] = data.roi.map(labels)
    return data, protected, labels


def compact_title(label: str) -> str:
    if " Hippocampus " in label:
        return label.replace(" Hippocampus ", "\nHippocampus ")
    return label.replace(" Mid L", "\nMid L") if label == "Ant Temp Mid L" else label


def common_axis(data: pd.DataFrame) -> tuple[float, float, np.ndarray]:
    lower = min(0, float(data.null_q025_centered.min())) * SCALE
    upper = max(float(data.observed_delta_mi.max()), float(data.null_q975_centered.max())) * SCALE
    lower = min(lower, float(data.observed_delta_mi.min()) * SCALE)
    width = upper - lower
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Invalid plot range")
    target_step = width / 3
    exponent = 10 ** math.floor(math.log10(target_step))
    step = next(m * exponent for m in (1, 2, 2.5, 5, 10) if m * exponent >= target_step)
    low = math.floor((lower - width * 0.07) / step) * step
    high = math.ceil((upper + width * 0.12) / step) * step
    return low, high, np.arange(low, high + step * 0.1, step)


def render(data: pd.DataFrame, labels: dict[str, str], output: Path,
           axis: tuple[float, float, np.ndarray]) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    low, high, ticks = axis
    fig, axes = plt.subplots(1, 4, figsize=(6.55, 3.15), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.245, right=0.988, bottom=0.255, top=0.825, wspace=0.21)
    blue, gray = "#398EAD", "#929AA0"
    for ax, roi in zip(axes, ROIS):
        cell = data.loc[data.roi.eq(roi)].set_index(["phase_band", "amp_band"])
        ax.axvline(0, color="#BFC4C7", lw=0.8, zorder=0)
        for y, combination in enumerate(COMBINATIONS):
            row = cell.loc[combination]
            observed = float(row.observed_delta_mi) * SCALE
            null_low = float(row.null_q025_centered) * SCALE
            null_high = float(row.null_q975_centered) * SCALE
            ax.axhline(y, color="#F0F1F2", lw=0.6, zorder=0)
            # The null range is NOT centered on the observed point.
            ax.hlines(y, null_low, null_high, colors=gray, linewidth=1.35, alpha=0.8, zorder=2)
            ax.plot([null_low, null_high], [y, y], "|", color=gray, markersize=4.6,
                    markeredgewidth=1.05, alpha=0.8, zorder=2)
            ax.plot(observed, y, "o", ms=4.5, mfc=blue, mec="white", mew=0.4, zorder=3)
            if row.plot_symbol:
                ax.annotate(row.plot_symbol, (observed, y), xytext=(0, 3.1),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=9, fontweight="bold", color="#242424")
        ax.set_title(compact_title(labels[roi]), fontsize=9.6, pad=8, color="#222222")
        ax.set_xlim(low, high)
        ax.set_ylim(5.5, -0.56)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{x:g}".replace("-", "−") for x in ticks])
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
    fig.text(0.61, 0.135, r"MI − mean surrogate MI ($\times10^{-3}$)",
             ha="center", fontsize=10, fontweight="bold")
    legend = fig.legend(handles=[Line2D([], [], color=blue, marker="o", linestyle="none", markersize=4.5),
                                 Line2D([], [], color=gray, linewidth=1.35)],
                        labels=["Observed MI excess", "Central 95% surrogate range"],
                        loc="center", bbox_to_anchor=(0.61, 0.055), ncol=2,
                        frameon=False, fontsize=7.4, handlelength=1.4,
                        handletextpad=0.45, columnspacing=1.05)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_items = [*fig.texts, *legend.get_texts()]
    for ax in axes:
        text_items.extend([ax.title, ax.yaxis.label, *ax.texts,
                           *ax.get_xticklabels(), *ax.get_yticklabels()])
    for item in text_items:
        if item.get_visible() and item.get_text():
            box = item.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                raise ValueError("Text outside figure: " + item.get_text())
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    result_root = project / "results/pac/by_lina/pac_existence"
    figure_root = project / "figures/pac/by_lina/pac_existence"
    csv_path = result_root / "pac_existence_forest_plot_data.csv"
    images = {scope: figure_root / f"forest_pac_existence_{scope}.png" for scope in SCOPES}
    outputs = [csv_path, *images.values(), figure_root / "forest_plot_audit.json",
               figure_root / "FDR_AND_PLOTTING_RULES.txt", result_root / "_SUCCESS"]
    if not args.replace_existing and any(path.exists() for path in outputs):
        raise FileExistsError("Output already exists; explicit replacement required")
    data, protected, labels = load_data(project)
    result_root.mkdir(parents=True, exist_ok=True)
    figure_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(figure_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.linewidth": 0.8, "axes.labelweight": "bold",
                         "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 9.5})
    axis = common_axis(data)
    for scope, output in images.items():
        render(data.loc[data.scope.eq(scope)], labels, output, axis)
    if any(sha256(Path(path)) != value for path, value in protected.items()):
        raise RuntimeError("Read-only input changed during plotting")
    data.to_csv(csv_path, index=False, encoding="utf-8-sig", float_format="%.17g")
    saved = pd.read_csv(csv_path)
    for col in ("group_observed_raw_pac", "group_surrogate_raw_pac_mean", "observed_delta_mi",
                "p_value", "q_value", "null_q025_centered", "null_q975_centered"):
        assert_close(saved[col], data[col], "saved CSV round-trip " + col)
    if not saved.plot_symbol.fillna("").eq(data.plot_symbol).all():
        raise ValueError("Saved significance marks differ")
    counts = {scope: {"total": len(cell), "raw_significant": int(cell.raw_significant.sum()),
                      "fdr_significant": int(cell.fdr_significant.sum())}
              for scope, cell in data.groupby("scope", sort=False)}
    audit = {"figures": {scope: str(path) for scope, path in images.items()},
             "figure_size_inches": [6.55, 3.15], "row_count": len(data), "counts": counts,
             "ROI_labels": labels, "source_files_sha256": protected,
             "source_hashes_unchanged": True, "plot_data_file": str(csv_path),
             "plot_data_sha256": sha256(csv_path), "surrogates_per_result": N_SURROGATES,
             "source_P_and_BH_q_verified_and_preserved": True,
             "family_sizes_per_scope": {"aMTG": 6, "three_hippocampi": 18},
             "source_analysis_unit": "pooled_segment", "scale_for_display_only": SCALE,
             "observed_point": "observed MI minus the matching surrogate mean",
             "horizontal_line": "surrogate quantiles 0.025 and 0.975 minus the matching surrogate mean",
             "horizontal_line_is_estimator_confidence_interval": False,
             "single_sided_P_not_deduced_from_two_sided_visual_range": True,
             "no_new_surrogates_or_model_fits": True, "all_visible_text_inside_canvas": True,
             "same_axes_for_all_ROIs_and_both_scopes": [axis[0], axis[1]],
             "symbols": data[KEYS + ["p_value", "q_value", "plot_symbol"]].to_dict(orient="records")}
    (figure_root / "forest_plot_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    notes = ("PAC/MI 本身的组水平存在性检验\n"
             "只复用李娜原始 pooled_segment 结果和已保存的1000次组级surrogate，不重算MI、置换或统计检验。\n"
             "两条分支分别出图，四脑区共享同一横轴范围。\n"
             "all_dynamic_split：按原结果排除终末corridor两个子段，保留终末junction。\n"
             "all_dynamic_split_with_terminal_C：在以上范围内纳入终末corridor两个子段。\n"
             "未改用我们后续同时排除终末corridor和junction的新规则，也未重算这些原始MI结果。\n"
             "FDR保持源结果：每分支aMTG六组合单独BH；三个海马18组合共同BH，未再按振幅组拆分。\n"
             "星号来自原始单侧经验P：*<0.05，**<0.01，***<0.001；†为原始BH q<0.05。\n"
             "蓝点：真实组MI减去同组surrogate均值。灰线：surrogate的2.5%–97.5%分位减去该均值。\n"
             "灰线是中心化后零分布的中央95%参考范围，不是MI或MI增量的估计量95%置信区间。\n"
             "单侧P与双侧中央95%绘图区间不等价：显著性始终使用原表P/q，不从图上区间重判。\n"
             "MI与每次surrogate均先在触点内平均纳入的窗，再在被试内触点等权，最后被试等权。\n"
             "仅绘图时将MI乘1000，坐标单位×10^-3，结果CSV保留原始未缩放数值。\n"
             "六行：30–70Hz振幅的0.5–4、4–12、12–30Hz相位，再为70–150Hz振幅的三个相位。\n"
             "脑区简称从config/project_config.json读取。未删除或覆盖任何李娜源表和旧图。\n")
    (figure_root / "FDR_AND_PLOTTING_RULES.txt").write_text(notes, encoding="utf-8")
    (result_root / "_SUCCESS").write_text("Saved surrogate data and original P/q verified; two PNGs generated.\n", encoding="utf-8")
    print(json.dumps({"figures": audit["figures"], "counts": counts,
                      "source_files_protected": len(protected), "no_source_edits": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
