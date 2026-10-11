"""
plot_granger_group_p_and_pair_proportion.py
==========================================
For each ROI pair, plot group-level Granger permutation p-values and the
proportion of significant electrode pairs.
Figure design
-------------
One figure per ROI pair and process_type.
Each figure has one row of six metric panels: T, D, TA, B, LG, HG.
Both directions share each panel: filled square / hollow square.
The direction legend is saved separately; there is no bar hatching.
In each subplot:
  left y-axis  : -log10(group-level p_greater), thin translucent lines and square markers.
  right y-axis : significant electrode-pair proportion, drawn as transparent bars.
In each metric panel, the x-axis is segment:
  corridor: C1, C2, ..., Last
  junction    : J1, J2, ..., Last
Line connections are across segments within each metric panel.
Inputs
------
1) Group-level summary CSV produced by plot_granger_null_histograms...:
   columns expected:
     process_type, brain_area1, brain_area2, data_cut,
     metric_scope, band_name, direction, p_greater, n_permutations
2) Electrode-pair all-channel CSV:
   expected path from cfg.pair_csv_path(...), usually:
     <data>/granger/<process_type>/allpairs/
       <prefix>_all_channel_pairs.csv
   columns expected:
     brain_area1, brain_area2, data_cut, metric_scope, band_name,
     pval_shuffle_1to2, pval_shuffle_2to1
Author: generated for SEEG navigation Granger posthoc plotting.
"""
from __future__ import annotations
import json
import math
import re
import gc
import traceback
import os
import sys
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib
# Keep TkAgg for interactive use in PyCharm; change to "Agg" if running on server.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# =========================================================
# User settings
# =========================================================
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import granger_posthoc_config as cfg
from granger_layout import PROJECT_ROOT
from granger_compact_style import STYLE, render_compact_figure, render_direction_legend

# Explicit formal inputs bypass the old histogram-summary fallback search.
formal_group_root = os.environ.get("GRANGER_POSTHOC_RESULTS_ROOT")
formal_figure_root = os.environ.get("GRANGER_FORMAL_FIGURE_ROOT")
if bool(formal_group_root) != bool(formal_figure_root):
    raise ValueError("Formal plotting requires both results and figure root overrides")
formal_mode = bool(formal_group_root)
data_path = cfg.data_path
input_root_path = cfg.input_root_path
process_types = cfg.process_types
brain_pairs = cfg.brain_pairs
gc_params = cfg.gc_params
pair_agg_method = cfg.pair_agg_method
# Group-level p-value column.
group_p_col = "p_greater"
# Electrode-pair p-value threshold for the significant proportion.
pair_alpha = cfg.alpha
# Horizontal significance line on the left axis.
show_alpha_line = True
alpha_line_value = -np.log10(pair_alpha)
# Output options.
save_png = True
save_pdf = False
show_fig = False
close_after_save = True
dpi = 300
# Figure style.
fontsize = 11
figsize = tuple(STYLE["figure_size_inches"])
# Bar style for significant proportion.
bar_alpha = None  # Explicit light fill colours, not a single bar alpha.
bar_width = STYLE["bar_width"]
# Left-axis line style.
line_lw = STYLE["line_width_pt"]
marker_size = 4.5
# y-axis limits.
left_ylim_min = 0.0
left_ylim_default_max = STYLE["left_axis_limits"][1]
right_ylim = tuple(STYLE["right_axis_limits"])
# Time-domain color: low-saturation neutral tone, distinct from the five band colors.
time_color = "#7F8C8D"
# Metric order.
METRICS = [
    {"label": "T",  "metric_scope": "time_domain",   "band_name": "time_domain",       "legend": "Time",        "title": "Time-domain GC"},
    {"label": "D",  "metric_scope": "spectral_band", "band_name": "slow_0p5_4Hz",         "legend": "Delta",       "title": "0.5–4 Hz"},
    {"label": "TA", "metric_scope": "spectral_band", "band_name": "theta_alpha_4_12",  "legend": "Theta-alpha", "title": "4–12 Hz"},
    {"label": "B",  "metric_scope": "spectral_band", "band_name": "beta_12_30",        "legend": "Beta",        "title": "12–30 Hz"},
    {"label": "LG", "metric_scope": "spectral_band", "band_name": "low_gamma_30_70",   "legend": "Low gamma",   "title": "30–70 Hz"},
    {"label": "HG", "metric_scope": "spectral_band", "band_name": "high_gamma_70_150", "legend": "High gamma",  "title": "70–150 Hz"},
]
DIRECTIONS = [
    {"key": "1to2", "title": "ROI1 \u2192 ROI2", "pair_p_cols": ["pval_shuffle_1to2", "p_shuffle_1to2", "pval_perm_1to2"]},
    {"key": "2to1", "title": "ROI2 \u2192 ROI1", "pair_p_cols": ["pval_shuffle_2to1", "p_shuffle_2to1", "pval_perm_2to1"]},
]
# =========================================================
# Optional project plot style
# =========================================================
try:
    from functions.my_plot_setup import set_plot
except Exception:
    set_plot = None
# =========================================================
# Formatting helpers
# =========================================================
def make_safe_filename(s: Any) -> str:
    s = str(s)
    s = ''.join(ch if ord(ch) >= 32 else "_" for ch in s)
    s = re.sub(r'[\\/:*?"<>|]+', "_", s)
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")
def read_roi_abbreviations() -> Dict[str, str]:
    """Use the project-wide display configuration; no local name overrides."""
    project_map = PROJECT_ROOT / "config/project_config.json"
    rows = json.loads(project_map.read_text(encoding="utf-8"))["roi_abbreviations"]
    return {row["full_name"]: row["abbreviation"] for row in rows}
ROI_ABBREV = read_roi_abbreviations()
def short_roi_name(s: Any) -> str:
    normalize = lambda value: " ".join(str(value).replace("_", " ").split()).casefold()
    labels = {normalize(key): value for key, value in ROI_ABBREV.items()}
    key = normalize(s)
    if key not in labels:
        raise KeyError("ROI missing from shared abbreviation config: " + str(s))
    return labels[key]
def normalize_direction(x: Any) -> str:
    """
    Normalize common direction labels to either "1to2" or "2to1".
    """
    s = str(x).strip().lower()
    s = s.replace(" ", "").replace("_", "").replace("-", "")
    s = s.replace("->", "to")
    if s in {"1to2", "12", "dir1", "direction1", "gc1to2", "roi1to2"}:
        return "1to2"
    if s in {"2to1", "21", "dir2", "direction2", "gc2to1", "roi2to1"}:
        return "2to1"
    if "1to2" in s:
        return "1to2"
    if "2to1" in s:
        return "2to1"
    return s
def format_segment_label(process_type: str, data_cut: Any) -> str:
    dc = str(data_cut)
    m = re.match(r"^corr(\d+)$", dc)
    if m:
        return f"C{m.group(1)}"
    if dc in {"corr_end", "corrend", "end", "last"}:
        return "Last"
    m = re.match(r"^junction(\d+)$", dc)
    if m:
        return f"J{m.group(1)}"
    if dc in {"junction_end", "turnend"}:
        return "Last"
    return dc
def natural_cut_key(x: Any) -> Tuple[int, int, str]:
    s = str(x)
    if s in {"corr_end", "junction_end", "end", "last"}:
        return (999, 999, s)
    m = re.search(r"(\d+)", s)
    if m:
        return (0, int(m.group(1)), s)
    return (500, 500, s)
def selected_cuts_for_process(process_type: str, df: Optional[pd.DataFrame] = None) -> List[str]:
    """
    Prefer cfg.get_selected_data_cuts if available.
    Otherwise infer from the data.
    """
    if hasattr(cfg, "get_selected_data_cuts"):
        try:
            cuts = cfg.get_selected_data_cuts(process_type)
            if cuts is not None and len(cuts) > 0:
                return list(cuts)
        except Exception:
            pass
    if hasattr(cfg, "selected_data_cuts_by_process"):
        try:
            cuts = cfg.selected_data_cuts_by_process.get(process_type)
            if cuts is not None and len(cuts) > 0:
                return list(cuts)
        except Exception:
            pass
    if df is not None and "data_cut" in df.columns:
        cuts = [str(x) for x in df["data_cut"].dropna().unique()]
        return sorted(cuts, key=natural_cut_key)
    return []
def get_metric_colors() -> Dict[str, str]:
    """
    Use:
      T = previous heatmap red
      D-HG = first five colors from set_plot() if available,
             otherwise Matplotlib default color cycle.
    """
    if set_plot is not None:
        try:
            tmp_fig, tmp_axs = set_plot(1, 1, figsize=(1, 1), fontsize=10)
            tmp_ax = np.ravel(np.atleast_1d(tmp_axs))[0]
            band_colors = [tmp_ax._get_lines.get_next_color() for _ in range(5)]
            plt.close(tmp_fig)
        except Exception as exc:
            print(f"[Warning] Could not read colors from set_plot(): {exc}")
            band_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][:5]
    else:
        band_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][:5]
    colors = {"T": time_color}
    for label, color in zip(["D", "TA", "B", "LG", "HG"], band_colors):
        colors[label] = color
    print("Metric colors:")
    for m in METRICS:
        print(f"  {m['label']}: {colors[m['label']]}")
    return colors
METRIC_COLOR = get_metric_colors()
# =========================================================
# Input path helpers
# =========================================================
def null_hist_root(process_type: str) -> Path:
    """
    Use cfg.null_hist_dir when available.
    """
    if hasattr(cfg, "null_hist_dir"):
        try:
            return cfg.null_hist_dir(process_type, pair_agg_method)
        except Exception:
            pass
    # Fallback.
    if hasattr(cfg, "output_base_path") and hasattr(cfg, "analysis_slug"):
        return (
                cfg.output_base_path
                / cfg.analysis_slug()
                / process_type
                / f"null_histograms_pairagg_{pair_agg_method}"
        )
    return (
            data_path
            / "granger_posthoc_results"
            / process_type
            / f"null_histograms_pairagg_{pair_agg_method}"
    )
def group_summary_search_dirs(process_type: str) -> List[Path]:
    """
    Return group-summary directories in priority order.

    Important:
    - Prefer the CURRENT posthoc-result directory only.
    - Use legacy locations only as fallbacks when the current directory
      does not exist.

    This avoids mixing current order_aligned_and_terminal_aligned summaries with older order_aligned
    summaries in the same DataFrame.
    """
    current_dir = null_hist_root(process_type) / "group_level" / "data"
    if current_dir.exists():
        return [current_dir]

    fallback_dirs = [
        input_root_path
        / process_type
        / "allpairs"
        / f"null_histograms_pairagg_{pair_agg_method}"
        / "group_level"
        / "data"
    ]

    posthoc_base = getattr(cfg, "output_base_path", data_path / "granger_posthoc_results")
    fallback_dirs.extend(
        sorted(
            posthoc_base.glob(
                f"*/{process_type}/null_histograms_pairagg_{pair_agg_method}/group_level/data"
            )
        )
    )

    out = []
    seen = set()
    for d in fallback_dirs:
        key = str(d)
        if key not in seen:
            out.append(d)
            seen.add(key)

    return out

def load_group_summary(process_type: str) -> pd.DataFrame:
    """
    Load all group-level summary CSV files for the current process_type.
    """
    if formal_mode:
        folder = Path(formal_group_root) / process_type / "allpairs" / f"group_level_pairagg_{pair_agg_method}"
        files = [folder / (cfg.make_gc_prefix(process_type, pair["brain_area1"], pair["brain_area2"])
                          + f"_group_level_pairagg_{pair_agg_method}_summary.csv") for pair in brain_pairs]
        if not all(path.is_file() for path in files):
            raise FileNotFoundError(f"Missing formal S3 inputs: {[str(path) for path in files if not path.is_file()]}")
        frames = []
        for path, pair in zip(files, brain_pairs):
            frame = pd.read_csv(path)
            if frame.empty or not frame["prefix"].eq(path.name.removesuffix(f"_group_level_pairagg_{pair_agg_method}_summary.csv")).all():
                raise ValueError(f"Empty or mismatched S3 prefix: {path}")
            if (not frame["process_type"].eq(process_type).all()
                    or not frame["brain_area1"].eq(pair["brain_area1"]).all()
                    or not frame["brain_area2"].eq(pair["brain_area2"]).all()):
                raise ValueError(f"Mismatched S3 ROI/process identity: {path}")
            if not frame["pair_agg_method"].eq("mean").all() or not frame["subject_agg_method"].eq("mean").all():
                raise ValueError("Formal Granger statistics must retain mean/mean aggregation")
            frame = frame.rename(columns={"group_p_greater": "p_greater",
                                           "n_permutations_effective": "n_permutations"})
            frame["__source_summary_csv"] = str(path)
            frames.append(frame)
        frame = pd.concat(frames, ignore_index=True)
        keys = ["process_type", "brain_area1", "brain_area2", "data_cut", "metric_scope", "band_name", "direction"]
        if frame.duplicated(keys).any():
            raise ValueError("Duplicate formal S3 statistical keys")
        if not frame["p_greater"].between(0, 1, inclusive="right").all():
            raise ValueError("Invalid formal empirical P values")
        frame["direction_norm"] = frame["direction"].map(normalize_direction)
        return frame
    files: List[Path] = []
    for d in group_summary_search_dirs(process_type):
        if d.exists():
            files.extend(sorted(d.glob("*group*summary*.csv")))
            files.extend(sorted(d.glob("*group_hist_summary.csv")))
    # Deduplicate.
    files = list(dict.fromkeys(files))
    if len(files) == 0:
        print("\n[Warning] No group-level summary CSV files found.")
        print("Searched:")
        for d in group_summary_search_dirs(process_type):
            print(f"  {d}")
        return pd.DataFrame()
    print(f"\nGroup summary CSV files for {process_type}:")
    for f in files:
        print(f"  {f}")
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f)
            df["__source_summary_csv"] = str(f)
            dfs.append(df)
        except Exception as exc:
            print(f"[Warning] Failed to read {f}: {exc}")
    if not dfs:
        return pd.DataFrame()
    df = pd.concat(dfs, ignore_index=True)
    # Normalize directions once.
    if "direction" in df.columns:
        df["direction_norm"] = df["direction"].map(normalize_direction)
    else:
        df["direction_norm"] = ""
    return df
def pair_csv_path(process_type: str, brain_area1: str, brain_area2: str) -> Path:
    """
    Prefer cfg.pair_csv_path, which should know the exact filename convention.
    """
    if hasattr(cfg, "pair_csv_path"):
        try:
            return cfg.pair_csv_path(process_type, brain_area1, brain_area2)
        except Exception:
            pass
    raise RuntimeError(
        "cfg.pair_csv_path(...) is not available. Please add it to granger_posthoc_config.py "
        "or adapt pair_csv_path() in this script."
    )
# =========================================================
# Data extraction
# =========================================================
def filter_metric(df: pd.DataFrame, metric: Dict[str, str]) -> pd.Series:
    """
    Boolean mask for one metric definition.
    """
    if df.empty:
        return pd.Series([], dtype=bool)
    metric_scope = df.get("metric_scope", pd.Series([""] * len(df))).astype(str)
    band_name = df.get("band_name", pd.Series([""] * len(df))).astype(str)
    if metric["label"] == "T":
        return (metric_scope == "time_domain") | (band_name == "time_domain")
    return (metric_scope == metric["metric_scope"]) & (band_name == metric["band_name"])
def get_group_p_value(
        group_df: pd.DataFrame,
        process_type: str,
        brain_area1: str,
        brain_area2: str,
        data_cut: str,
        metric: Dict[str, str],
        direction_key: str,
) -> Tuple[float, float]:
    """
    Return p_greater and n_permutations for one condition.

    Memory note:
    Work directly on group_df; do not make a full DataFrame copy for every
    segment x metric x direction query.
    """
    if group_df.empty:
        return np.nan, np.nan

    df = group_df
    mask = pd.Series(True, index=df.index)

    if "process_type" in df.columns:
        mask &= df["process_type"].astype(str).eq(str(process_type))
    if "brain_area1" in df.columns:
        mask &= df["brain_area1"].astype(str).eq(str(brain_area1))
    if "brain_area2" in df.columns:
        mask &= df["brain_area2"].astype(str).eq(str(brain_area2))
    if "data_cut" in df.columns:
        mask &= df["data_cut"].astype(str).eq(str(data_cut))
    if "direction_norm" in df.columns:
        mask &= df["direction_norm"].eq(direction_key)

    mask &= filter_metric(df, metric)

    if not mask.any():
        return np.nan, np.nan

    if group_p_col not in df.columns:
        raise RuntimeError(
            f"Group summary does not contain p-value column '{group_p_col}'. "
            f"Available columns:\n{df.columns.tolist()}"
        )

    p_vals = pd.to_numeric(df.loc[mask, group_p_col], errors="coerce").dropna()
    if p_vals.empty:
        return np.nan, np.nan

    p = float(p_vals.iloc[0])

    if "n_permutations" in df.columns:
        nperm_vals = pd.to_numeric(
            df.loc[mask, "n_permutations"], errors="coerce"
        ).dropna()
        nperm = float(nperm_vals.iloc[0]) if not nperm_vals.empty else np.nan
    else:
        nperm = np.nan

    return p, nperm

def p_to_neglog10(p: float, nperm: float = np.nan) -> float:
    if not np.isfinite(p):
        return np.nan
    if p <= 0:
        if np.isfinite(nperm) and nperm > 0:
            p = 1.0 / (nperm + 1.0)
        else:
            p = 1e-12
    return -float(np.log10(p))
def get_pair_p_col(pair_df: pd.DataFrame, direction: Dict[str, Any]) -> str:
    for c in direction["pair_p_cols"]:
        if c in pair_df.columns:
            return c
    raise RuntimeError(
        f"Cannot find electrode-pair p-value column for direction {direction['key']}.\n"
        f"Tried: {direction['pair_p_cols']}\n"
        f"Available columns:\n{pair_df.columns.tolist()}"
    )
def get_significant_proportion(
        pair_df: pd.DataFrame,
        data_cut: str,
        metric: Dict[str, str],
        p_col: str,
) -> Tuple[float, int, int]:
    """
    Return significant proportion, n_significant, n_valid.

    Memory note:
    Do not copy the full all-channel DataFrame for every query.
    """
    if pair_df.empty:
        return np.nan, 0, 0

    df = pair_df

    mask = pd.Series(True, index=df.index)
    if "data_cut" in df.columns:
        mask &= df["data_cut"].astype(str).eq(str(data_cut))

    mask &= filter_metric(df, metric)

    if not mask.any() or p_col not in df.columns:
        return np.nan, 0, 0

    p = pd.to_numeric(df.loc[mask, p_col], errors="coerce")
    valid = p[np.isfinite(p)]

    n_valid = int(valid.size)
    if n_valid == 0:
        return np.nan, 0, 0

    n_sig = int((valid < pair_alpha).sum())
    prop = n_sig / n_valid

    return float(prop), n_sig, n_valid


def load_pair_csv_minimal(pair_csv: Path) -> pd.DataFrame:
    """
    Read only the columns that this plotting script actually needs.

    The all-channel Granger CSV can be very large. Reading every column is
    unnecessary and can substantially increase RAM use.
    """
    header = pd.read_csv(pair_csv, nrows=0)
    available = set(header.columns)

    required_base = ["data_cut", "metric_scope", "band_name"]
    candidate_p_cols = []
    for direction in DIRECTIONS:
        candidate_p_cols.extend(direction["pair_p_cols"])

    usecols = [c for c in required_base + candidate_p_cols if c in available]

    missing_base = [c for c in required_base if c not in available]
    if missing_base:
        raise RuntimeError(
            f"Pair CSV is missing required columns: {missing_base}\n"
            f"Available columns:\n{header.columns.tolist()}"
        )

    print(f"  Reading only {len(usecols)} needed columns from pair CSV")
    return pd.read_csv(pair_csv, usecols=usecols)

def collect_plot_data(
        group_df: pd.DataFrame,
        pair_df: pd.DataFrame,
        process_type: str,
        brain_area1: str,
        brain_area2: str,
        data_cuts: List[str],
        direction: Dict[str, Any],
) -> pd.DataFrame:
    """
    Build a long table used for plotting one ROI pair and one direction.
    """
    p_col = get_pair_p_col(pair_df, direction)
    rows = []
    for s_idx, data_cut in enumerate(data_cuts):
        for m_idx, metric in enumerate(METRICS):
            x = s_idx * len(METRICS) + m_idx
            p_group, nperm = get_group_p_value(
                group_df=group_df,
                process_type=process_type,
                brain_area1=brain_area1,
                brain_area2=brain_area2,
                data_cut=data_cut,
                metric=metric,
                direction_key=direction["key"],
            )
            prop, n_sig, n_valid = get_significant_proportion(
                pair_df=pair_df,
                data_cut=data_cut,
                metric=metric,
                p_col=p_col,
            )
            rows.append({
                "process_type": process_type,
                "brain_area1": brain_area1,
                "brain_area2": brain_area2,
                "direction": direction["key"],
                "data_cut": data_cut,
                "segment_label": format_segment_label(process_type, data_cut),
                "metric_label": metric["label"],
                "metric_legend": metric["legend"],
                "metric_scope": metric["metric_scope"],
                "band_name": metric["band_name"],
                "x": x,
                "group_p": p_group,
                "n_permutations": nperm,
                "neglog10_p": p_to_neglog10(p_group, nperm),
                "sig_proportion": prop,
                "n_sig_pairs": n_sig,
                "n_valid_pairs": n_valid,
                "pair_p_col": p_col,
            })
    return pd.DataFrame(rows)
# =========================================================
# Plotting
# =========================================================


def plot_one_roi_pair(
        process_type: str,
        brain_area1: str,
        brain_area2: str,
        group_df: pd.DataFrame,
        pair_df: pd.DataFrame,
        data_cuts: List[str],
        out_fig_dir: Path,
        out_data_dir: Path,
):
    prefix = make_safe_filename(
        f"{process_type}_{brain_area1}_to_{brain_area2}_groupP_pairSigProp_byMetric"
    )
    out_png = out_fig_dir / f"{prefix}.png"
    out_csv = out_data_dir / f"{prefix}_plot_data.csv"
    if formal_mode:
        for target in (out_csv, out_png):
            if target.exists():
                raise FileExistsError(f"Refusing to overwrite an existing formal plot: {target}")
    # Statistical extraction and aggregation remain identical to the original entry.
    all_plot_tables = [collect_plot_data(
        group_df=group_df, pair_df=pair_df, process_type=process_type,
        brain_area1=brain_area1, brain_area2=brain_area2, data_cuts=data_cuts,
        direction=direction,
    ) for direction in DIRECTIONS]
    plot_table = pd.concat(all_plot_tables, ignore_index=True)
    if formal_mode:
        plot_table["branch"] = cfg.analysis_kind
        plot_table["source_pair_csv"] = str(pair_csv_path(process_type, brain_area1, brain_area2))
        plot_table["source_group_csv"] = str(Path(formal_group_root) / process_type / "allpairs"
            / f"group_level_pairagg_{pair_agg_method}"
            / (cfg.make_gc_prefix(process_type, brain_area1, brain_area2)
               + f"_group_level_pairagg_{pair_agg_method}_summary.csv"))
    render_compact_figure(
        plot_table, data_cuts, METRICS, METRIC_COLOR,
        short_roi_name(brain_area1), short_roi_name(brain_area2),
        [format_segment_label(process_type, cut) for cut in data_cuts], out_png, dpi=dpi,
    )
    plot_table.to_csv(out_csv, index=False)
    print(f"  Saved compact PNG: {out_png}")
    print(f"  Saved plotting data: {out_csv}")


# =========================================================
# Main
# =========================================================
def main():
    formal_outputs = []
    for process_type in process_types:
        print("\n============================================================")
        print(f"Process type: {process_type}")
        print("============================================================")

        group_df = load_group_summary(process_type)
        if group_df.empty:
            print(f"[Skip] No group summary data for process_type={process_type}")
            continue

        out_root = (Path(formal_figure_root) / cfg.analysis_kind / process_type / "group_p_and_pair_sig_proportion"
                    if formal_mode else null_hist_root(process_type) / "group_p_and_pair_sig_proportion")
        out_fig_dir = out_root / "fig"
        out_data_dir = out_root / "data"
        out_fig_dir.mkdir(parents=True, exist_ok=True)
        out_data_dir.mkdir(parents=True, exist_ok=True)

        data_cuts = selected_cuts_for_process(process_type, group_df)
        if len(data_cuts) == 0:
            print(f"[Skip] Could not determine data cuts for {process_type}")
            continue

        print(f"Data cuts: {data_cuts}")
        print(f"Output dir: {out_root}")

        for pair in brain_pairs:
            brain_area1 = pair["brain_area1"]
            brain_area2 = pair["brain_area2"]

            print("\n------------------------------------------------------------")
            print(f"ROI pair: {brain_area1}  <->  {brain_area2}")

            pair_csv = pair_csv_path(process_type, brain_area1, brain_area2)
            print(f"Pair CSV: {pair_csv}")

            if not pair_csv.exists():
                if formal_mode:
                    raise FileNotFoundError(pair_csv)
                print(f"  [Skip] Pair CSV not found: {pair_csv}")
                continue

            pair_df = None

            try:
                pair_df = load_pair_csv_minimal(pair_csv)

                mem_mb = pair_df.memory_usage(deep=True).sum() / (1024 ** 2)
                print(
                    f"  Pair rows: {len(pair_df):,}; "
                    f"in-memory size used by plot: {mem_mb:.1f} MB"
                )

                plot_one_roi_pair(
                    process_type=process_type,
                    brain_area1=brain_area1,
                    brain_area2=brain_area2,
                    group_df=group_df,
                    pair_df=pair_df,
                    data_cuts=data_cuts,
                    out_fig_dir=out_fig_dir,
                    out_data_dir=out_data_dir,
                )
                if formal_mode:
                    prefix = make_safe_filename(f"{process_type}_{brain_area1}_to_{brain_area2}_groupP_pairSigProp_byMetric")
                    formal_outputs.append({"png": str(out_fig_dir / f"{prefix}.png"),
                                           "plot_data_csv": str(out_data_dir / f"{prefix}_plot_data.csv"),
                                           "group_source": str(Path(formal_group_root) / process_type / "allpairs"
                                               / f"group_level_pairagg_{pair_agg_method}"
                                               / (cfg.make_gc_prefix(process_type, brain_area1, brain_area2)
                                                  + f"_group_level_pairagg_{pair_agg_method}_summary.csv")),
                                           "pair_source": str(pair_csv)})

            except Exception as exc:
                if formal_mode:
                    raise
                print(f"  [Error] Failed to plot {brain_area1} -> {brain_area2}")
                print(f"          {type(exc).__name__}: {exc!r}")
                traceback.print_exc()

            finally:
                if pair_df is not None:
                    del pair_df
                plt.close("all")
                gc.collect()

        del group_df
        gc.collect()

    if formal_mode:
        if len(formal_outputs) != len(brain_pairs) * len(process_types):
            raise RuntimeError("Incomplete formal plot coverage")
        manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                    "branch": cfg.analysis_kind, "layout": STYLE["layout"],
                    "script": str(Path(__file__).resolve()),
                    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    "config_path": cfg._settings["config_path"],
                    "multiple_comparison_correction": "none", "group_p_column": "group_p_greater",
                    "pair_alpha": pair_alpha, "pair_proportion": "all finite pair P values pooled; n(P<0.05)/n(finite P)",
                    "group_aggregation": "equal pair means within subject, then equal subject means",
                    "png_count": len(formal_outputs), "bar_alpha": bar_alpha, "bar_width": bar_width,
                    "figure_size_inches": list(figsize), "dpi": dpi, "metric_colors": METRIC_COLOR,
                    "left_y_minimum_maximum": left_ylim_default_max,
                    "outputs": formal_outputs}
        common_sources = {short_roi_name(pair["brain_area1"]) for pair in brain_pairs}
        if len(common_sources) != 1:
            raise ValueError("A shared direction legend requires a common source ROI")
        legend_target = Path(formal_figure_root) / "direction_legend.png"
        render_direction_legend(next(iter(common_sources)), legend_target, dpi=dpi)
        helper_file = Path(__file__).with_name("granger_compact_style.py")
        roi_config_file = PROJECT_ROOT / "config/project_config.json"
        manifest.update(STYLE)
        manifest.update({"render_helper": str(helper_file),
                         "render_helper_sha256": hashlib.sha256(helper_file.read_bytes()).hexdigest(),
                         "roi_config_sha256": hashlib.sha256(roi_config_file.read_bytes()).hexdigest(),
                         "shared_direction_legend": str(legend_target), "legend_item_count": 2})
        target = Path(formal_figure_root) / cfg.analysis_kind / "reference_plot_manifest.json"
        with target.open("x", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    print("\nDone.")


if __name__ == "__main__":
    main()
