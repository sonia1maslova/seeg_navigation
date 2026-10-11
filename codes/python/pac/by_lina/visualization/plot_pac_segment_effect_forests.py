"""Read-only PAC segment-position forests for original and adjusted Lina models.

No signals are extracted and no models or statistical result files are changed.
Stored BH q-values are independently checked in their COMPLETE source families.
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
AMPS = ("low_gamma_30_70Hz", "high_gamma_70_150Hz")
PHASE_LABELS = ("0.5–4", "4–12", "12–30")
COLORS = ("#40AE80", "#6B4C9A")
DISPLAY_UNIT = 1e-5
ORIGINAL_SCOPES = {"include_last": "all_dynamic_split_with_terminal_C",
                   "exclude_last": "all_dynamic_split"}
ADJUSTED_SCOPES = ("include_last", "exclude_terminal_corridor_and_junction")
ANALYSES = ("main_pac", "power_control", "aperiodic_control_1_150Hz")
NUMERIC_COLUMNS = ("beta", "se", "t_value", "df", "ci_low", "ci_high",
                   "p_value", "q_value", "n_obs", "n_subjects", "n_segments")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value: str) -> str:
    return " ".join(value.replace("_", " ").split()).casefold()


def bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1.0)
    return result


def raw_stars(p: float) -> str:
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def validate_values(data: pd.DataFrame, context: str) -> None:
    if not data["fit_status"].eq("ok").all():
        raise ValueError(context + ": contains a failed fit")
    if not np.isfinite(data[list(NUMERIC_COLUMNS)].to_numpy(float)).all():
        raise ValueError(context + ": non-finite fitted number")
    if not (data["ci_low"].le(data["beta"]) & data["beta"].le(data["ci_high"])).all():
        raise ValueError(context + ": beta outside its stored CI")
    if not (data["se"].gt(0) & data["df"].gt(0)).all():
        raise ValueError(context + ": invalid SE or degrees of freedom")
    if not (data["p_value"].between(0, 1) & data["q_value"].between(0, 1)).all():
        raise ValueError(context + ": invalid p or q")
    if not data["predictor"].eq("segment_position").all() or not data["target_term"].eq("segment_position").all():
        raise ValueError(context + ": wrong fixed effect")


def verify_grid(data: pd.DataFrame, rois: tuple[str, ...], context: str,
                ranges: tuple[str, ...] | None = None) -> None:
    keys = ["roi", "phase_band", "amp_band"]
    expected = {(roi, phase, amp) for roi in rois for phase in PHASES for amp in AMPS}
    if ranges:
        keys.append("exponent_range")
        expected = {(*key, rng) for key in expected for rng in ranges}
    if len(data) != len(expected) or data.duplicated(keys).any():
        raise ValueError(context + ": incomplete or duplicate model grid")
    if set(data[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError(context + ": unexpected ROI/frequency/range combination")
    validate_values(data, context)


def verify_families(data: pd.DataFrame, expected_size: int, context: str,
                    audit: list[dict]) -> None:
    for family_id, family in data.groupby("fdr_family_id", sort=False):
        if len(family) != expected_size or not family["m_planned"].eq(expected_size).all() or not family["m_valid"].eq(expected_size).all():
            raise ValueError(context + ": incomplete FDR family " + family_id)
        if family["scope"].nunique() != 1:
            raise ValueError(context + ": scope mixed in FDR family")
        delta = float(np.max(np.abs(bh(family["p_value"].to_numpy(float)) - family["q_value"].to_numpy(float))))
        if delta > 1e-12:
            raise ValueError(context + ": stored q fails independent BH check")
        audit.append({"version_analysis": context, "fdr_family_id": family_id,
                      "scope": family["scope"].iloc[0], "tests": expected_size,
                      "max_abs_q_difference": delta,
                      "rois": list(family["roi"].unique()),
                      "ranges": list(family.get("exponent_range", pd.Series(dtype=str)).dropna().unique())})


def read_original(project: Path, protected: dict[str, str], audit: list[dict]) -> list[dict]:
    root = project / "results/pac/by_lina/formal_raw_pac"
    specs = (("02_hpc_path_accumulation", "02A_basic_segment", "single_roi", ROIS[1:], 18),
             ("04_amtg_followup_and_behavior", "04B_amtg_accumulation", "single_roi", ROIS[:1], 6),
             ("03a_hpc_power_control", "03A2_power_adjusted_PAC", "complete_case_adjusted", ROIS[1:], 18),
             ("03b_hpc_exponent_control", "03B2_exponent_adjusted_PAC", "complete_case_adjusted", ROIS[1:], 36))
    selected: dict[str, pd.DataFrame] = {}
    for folder, name, sample_mode, rois, family_size in specs:
        path = root / folder / "model_coefficients.csv"
        protected[str(path)] = sha256(path)
        source = pd.read_csv(path)
        data = source.loc[source["analysis_family"].eq(name) & source["sample_mode"].eq(sample_mode)
                          & source["scope"].isin(ORIGINAL_SCOPES.values())
                          & source["predictor"].eq("segment_position")
                          & source["target_term"].eq("segment_position")
                          & source["roi"].isin(rois)].copy()
        ranges = ("1-30 Hz", "1-150 Hz") if family_size == 36 else None
        for scope in ORIGINAL_SCOPES.values():
            verify_grid(data.loc[data["scope"].eq(scope)], rois, name + " / " + scope, ranges)
        if data["fdr_family_id"].nunique() != 2:
            raise ValueError("Expected two scope-specific original FDR families: " + name)
        verify_families(data, family_size, "lina_original / " + name, audit)
        data["source_file"] = str(path)
        selected[name] = data
    plots = []
    for shown_scope, stored_scope in ORIGINAL_SCOPES.items():
        main = pd.concat([selected["04B_amtg_accumulation"], selected["02A_basic_segment"]], ignore_index=True)
        datasets = {"main_pac": main,
                    "power_control": selected["03A2_power_adjusted_PAC"],
                    "aperiodic_control_1_150Hz": selected["03B2_exponent_adjusted_PAC"].loc[
                        selected["03B2_exponent_adjusted_PAC"]["exponent_range"].eq("1-150 Hz")]}
        for analysis, frame in datasets.items():
            data = frame.loc[frame["scope"].eq(stored_scope)].copy()
            verify_grid(data, ROIS if analysis == "main_pac" else ROIS[1:], analysis)
            plots.append({"version": "lina_original", "scope": shown_scope, "analysis": analysis, "data": data})
    return plots


def verify_adjusted_subtables(root: Path, data: pd.DataFrame, protected: dict[str, str]) -> None:
    keys = ["branch", "scope", "roi", "phase_band", "amp_band", "exponent_range"]
    for folder, branches in (("main_pac", ("hpc_main", "amtg_main")),
                             ("power_control", ("power_control",)),
                             ("aperiodic_control", ("aperiodic_control",))):
        path = root / folder / "model_coefficients.csv"
        protected[str(path)] = sha256(path)
        component = pd.read_csv(path)
        expected = data.loc[data["branch"].isin(branches)]
        if len(component) != len(expected):
            raise ValueError("Component model coefficients have different length: " + folder)
        joined = expected.merge(component, on=keys, validate="one_to_one", suffixes=("", "_component"))
        if len(joined) != len(expected):
            raise ValueError("Adjusted component keys do not match")
        for column in NUMERIC_COLUMNS:
            if not np.allclose(joined[column], joined[column + "_component"], rtol=1e-12, atol=1e-15):
                raise ValueError("Adjusted combined statistic differs from component: " + column)


def read_adjusted(project: Path, protected: dict[str, str], audit: list[dict]) -> list[dict]:
    plots = []
    for scope in ADJUSTED_SCOPES:
        root = project / "results/pac/by_lina/segment_effect_adjusted" / scope
        if not (root / "_SUCCESS").exists():
            raise ValueError("Incomplete adjusted model batch: " + scope)
        path = root / "all_segment_effects.csv"
        protected[str(path)] = sha256(path)
        data = pd.read_csv(path)
        if len(data) != 78 or not data["scope"].eq(scope).all() or not data["term"].eq("segment_position").all():
            raise ValueError("Expected 78 adjusted segment-position effects: " + scope)
        branches = {"hpc_main": (ROIS[1:], None), "amtg_main": (ROIS[:1], None),
                    "power_control": (ROIS[1:], None),
                    "aperiodic_control": (ROIS[1:], ("1-30 Hz", "1-150 Hz"))}
        if set(data["branch"]) != set(branches):
            raise ValueError("Unexpected adjusted model branch")
        for branch, (rois, ranges) in branches.items():
            frame = data.loc[data["branch"].eq(branch)]
            verify_grid(frame, rois, scope + " / " + branch, ranges)
        if data["fdr_family_id"].nunique() != 13:
            raise ValueError("Expected 13 adjusted six-test FDR families")
        for family_id, family in data.groupby("fdr_family_id", sort=False):
            if family["branch"].nunique() != 1 or family["roi"].nunique() != 1:
                raise ValueError("Adjusted six-test family mixed ROI or branch")
            if family["branch"].iloc[0] == "aperiodic_control" and family["exponent_range"].nunique() != 1:
                raise ValueError("Adjusted six-test family mixed exponent ranges")
            verify_grid(family, (family["roi"].iloc[0],), family_id)
        verify_families(data, 6, "adjusted / " + scope, audit)
        verify_adjusted_subtables(root, data, protected)
        data["source_file"] = str(path)
        datasets = {"main_pac": data.loc[data["branch"].isin(("hpc_main", "amtg_main"))],
                    "power_control": data.loc[data["branch"].eq("power_control")],
                    "aperiodic_control_1_150Hz": data.loc[data["branch"].eq("aperiodic_control")
                                                            & data["exponent_range"].eq("1-150 Hz")]}
        for analysis, frame in datasets.items():
            verify_grid(frame, ROIS if analysis == "main_pac" else ROIS[1:], analysis)
            plots.append({"version": "adjusted", "scope": scope, "analysis": analysis, "data": frame.copy()})
    return plots


def compact_title(label: str) -> str:
    if " Hippocampus " in label:
        return label.replace(" Hippocampus ", "\nHippocampus ")
    if label == "Ant Temp Mid L":
        return label.replace(" Mid L", "\nMid L")
    return label


def axis_extent(frames: list[pd.DataFrame]) -> float:
    largest = max(float(np.abs(x[["ci_low", "ci_high"]].to_numpy()).max()) for x in frames)
    padded = largest / DISPLAY_UNIT * 1.12
    scale = 10 ** math.floor(math.log10(padded))
    for multiplier in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if multiplier * scale >= padded:
            return multiplier * scale
    raise AssertionError("No extent")


def render(frame: pd.DataFrame, labels: dict[str, str], output: Path,
           extent: float) -> dict:
    import matplotlib.pyplot as plt

    rois = ROIS if len(frame) == 24 else ROIS[1:]
    width = 6.30 if len(rois) == 4 else 5.15
    fig, axes = plt.subplots(1, len(rois), figsize=(width, 2.42), sharex=True, sharey=True)
    left = .16 if len(rois) == 4 else .195
    fig.subplots_adjust(left=left, right=.99, bottom=.245, top=.79, wspace=.20)
    points: list[dict] = []
    for ax, roi in zip(axes, rois):
        data = frame.loc[frame["roi"].eq(roi)].set_index(["phase_band", "amp_band"])
        ax.axvline(0, lw=.8, color="#BDC3C7", zorder=0)
        for row_index, phase in enumerate(PHASES):
            ax.axhline(row_index, lw=.55, color="#EFF1F2", zorder=0)
            for amp_index, (amp, color) in enumerate(zip(AMPS, COLORS)):
                row = data.loc[(phase, amp)]
                y = row_index + (-.145 if amp_index == 0 else .145)
                estimate, low, high = np.array([row.beta, row.ci_low, row.ci_high]) / DISPLAY_UNIT
                alpha = 1.0 if row.p_value < .05 else .35
                ax.errorbar(estimate, y, xerr=np.array([[estimate - low], [high - estimate]]),
                            fmt="none", ecolor=color, elinewidth=1.25, capsize=2.1,
                            capthick=.95, alpha=alpha * .65, zorder=2)
                ax.plot(estimate, y, "o", ms=4.3, mfc=color, mec="white", mew=.35,
                        alpha=alpha, zorder=3)
                symbol = raw_stars(float(row.p_value)) + ("†" if row.q_value < .05 else "")
                if symbol:
                    # Put each marker's symbols on its outer side. Drawing the
                    # lower (purple) symbol above would obscure its green partner.
                    symbol_offset = 2.8 if amp_index == 0 else -2.8
                    ax.annotate(symbol, (estimate, y), xytext=(0, symbol_offset), textcoords="offset points",
                                ha="center", va="bottom" if amp_index == 0 else "top",
                                fontsize=8.8, fontweight="bold", color="#242424")
                points.append({"roi": roi, "phase_band": phase, "amp_band": amp,
                               "beta_raw": float(row.beta), "ci_low_raw": float(row.ci_low),
                               "ci_high_raw": float(row.ci_high), "beta_display": float(estimate),
                               "p_value": float(row.p_value), "q_value": float(row.q_value),
                               "symbol": symbol, "alpha": alpha, "y": y,
                               "symbol_offset_points": [0, 2.8 if amp_index == 0 else -2.8]})
        ax.set_title(compact_title(labels[roi]), fontsize=9.1, pad=7, fontweight="bold")
        ax.set_xlim(-extent, extent)
        ax.set_ylim(2.50, -.55)
        ax.set_xticks([-extent / 2, 0, extent / 2])
        ax.set_xticklabels([f"−{extent / 2:g}", "0", f"{extent / 2:g}"])
        ax.set_yticks(np.arange(3))
        ax.tick_params(axis="x", length=2.5, pad=2, width=.75, colors="#555555")
        ax.tick_params(axis="y", length=0, pad=7)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#777777")
        ax.spines["bottom"].set_linewidth(.7)
    axes[0].set_yticklabels(PHASE_LABELS)
    axes[0].set_ylabel("Phase frequency (Hz)", fontsize=9.5, labelpad=4)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.text((left + .99) / 2, .095, r"Segment number effect, $\beta$ ($\times10^{-5}$ MI/segment)",
             ha="center", fontsize=9.5, fontweight="bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        for item in [ax.title, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]:
            if item.get_visible() and item.get_text():
                box = item.get_window_extent(renderer)
                if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                    raise ValueError("Text outside canvas: " + item.get_text())
    for item in fig.texts:
        box = item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            raise ValueError("Axis label outside canvas")
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return {"size_inches": [width, 2.42], "shared_xlim_display_units": [-extent, extent],
            "display_unit": DISPLAY_UNIT, "points": points,
            "all_visible_text_inside_canvas": True, "all_stored_95_CIs_shown": True}


def render_legend(output: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(4.90, .45))
    ax.axis("off")
    handles = [Line2D([0], [0], color=color, marker="o", lw=1.2, markersize=4)
               for color in COLORS]
    fig.legend(handles, ["Low gamma (30–70 Hz)", "High gamma (70–150 Hz)"],
               loc="center", ncol=2, frameon=False, fontsize=9,
               handlelength=1.5, handletextpad=.5, columnspacing=1.6)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.getenv("SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")))
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    output_root = args.output_root or project / "figures/pac/by_lina/pac_segment_effect"
    output_root.mkdir(parents=True, exist_ok=True)
    protected: dict[str, str] = {}
    audit: list[dict] = []
    config_path = project / "config/project_config.json"
    protected[str(config_path)] = sha256(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    lookup = {normalize(x["full_name"]): x["abbreviation"] for x in config["roi_abbreviations"]}
    labels = {roi: lookup[normalize(roi)] for roi in ROIS}
    plots = read_original(project, protected, audit) + read_adjusted(project, protected, audit)
    if len(plots) != 12 or sum(len(x["data"]) for x in plots) != 240:
        raise ValueError("Expected 12 forests and 240 displayed model effects")
    os.environ.setdefault("MPLCONFIGDIR", str(output_root / ".mpl_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in available else "DejaVu Sans",
                         "font.size": 9, "axes.labelweight": "bold", "axes.titleweight": "bold",
                         "xtick.labelsize": 8.2, "ytick.labelsize": 9.3})
    extents = {analysis: axis_extent([x["data"] for x in plots if x["analysis"] == analysis])
               for analysis in ANALYSES}
    manifest = {"model_fitting": False, "source_statistics_modified": False,
                "shared_standard_ROI_labels": labels, "original_scope_mapping": ORIGINAL_SCOPES,
                "adjusted_scopes": ADJUSTED_SCOPES, "independent_BH_audit": audit, "figures": []}
    for item in plots:
        folder = output_root / "lmm_effects" / item["version"] / item["scope"] / item["analysis"] / "forest_plots"
        folder.mkdir(parents=True, exist_ok=True)
        output = folder / "forest_Segment_number_all_segments.png"
        if output.exists() and not args.replace_existing:
            raise FileExistsError("Forest exists; use --replace-existing explicitly: " + str(output))
        frame = item["data"].copy()
        frame["roi_display"] = frame["roi"].map(labels)
        frame["plot_symbol"] = frame["p_value"].map(raw_stars) + np.where(frame["q_value"] < .05, "†", "")
        info = render(frame, labels, output, extents[item["analysis"]])
        frame.to_csv(folder / "forest_plot_data.csv", index=False)
        info.update({"version": item["version"], "scope": item["scope"], "analysis": item["analysis"],
                     "figure_file": output.relative_to(output_root).as_posix(), "row_count": len(frame),
                     "raw_significant": int((frame["p_value"] < .05).sum()),
                     "FDR_significant": int((frame["q_value"] < .05).sum())})
        (folder / "forest_plot_manifest.json").write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
        manifest["figures"].append(info)
    legend_dir = output_root / "lmm_effects"
    legend_path = legend_dir / "forest_amplitude_legend.png"
    if legend_path.exists() and not args.replace_existing:
        raise FileExistsError("Common legend exists; explicit replacement required")
    render_legend(legend_path)
    notes = (
        "PAC segment-position forests: plotting and FDR provenance\n"
        "直接复用已拟合模型的 beta、95%CI、原始 p 与已保存 q；不重算信号、不重新拟合、不修改结果表。\n"
        "李娜原始主模型：aMTG 6 项、三个海马共同 18 项；功率控制 18 项；1/f 控制按两范围共同 36 项核验后只显示 1–150 Hz。\n"
        "调整版：每一分析范围、模型分支、脑区及 1/f 范围各自 6 项 BH-FDR；完整 78 项模型、13 个检验族均独立核验。\n"
        "两版的原 q 口径分别保留；不同范围、版本或展示方式不合并校正，也不重复进入 FDR。\n"
        "原始 exclude_last 只排除末 corridor 两子段；调整 exclude_terminal_corridor_and_junction 同时排除末 corridor 两子段和末 junction。\n"
        "主 PAC 4 脑区；控制模型仅 3 海马。脑区标签来自 config/project_config.json 的 roi_abbreviations。\n"
        "绿色与紫色分别为 30–70 和 70–150 Hz 振幅，三行分别为 0.5–4、4–12、12–30 Hz 相位。\n"
        "点为 beta，水平线为源模型未作多重比较调整的 95%CI；各对应分析跨版本、范围使用相同横轴。\n"
        "仅显示单位缩放为 10^-5 MI/segment，保存 CSV 中 beta 与 CI 保持原数值。\n"
        "*：原始 p<0.05；**：p<0.01；***：p<0.001；†：源表 BH-FDR q<0.05。\n"
        "原始 p<0.05 点完全不透明，否则半透明；星号/†贴近对应点，不在每张图下重复长说明。\n"
        "共同振幅图例单独保存为 forest_amplitude_legend.png。\n")
    (legend_dir / "FOREST_FDR_AND_PLOTTING_RULES.txt").write_text(notes, encoding="utf-8")
    if any(sha256(Path(path)) != digest for path, digest in protected.items()):
        raise RuntimeError("A protected input file changed during read-only plotting")
    manifest.update({"input_sha256": protected, "protected_input_hashes_unchanged": True,
                     "figure_count": 12, "common_legend_file": legend_path.relative_to(output_root).as_posix(),
                     "plotted_effect_count": 240})
    (legend_dir / "forest_generation_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"forest_count": 12, "common_legend_count": 1,
                      "effect_count": 240, "BH_families_verified": len(audit),
                      "protected_inputs": len(protected), "input_hashes_unchanged": True,
                      "output_root": str(output_root)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
