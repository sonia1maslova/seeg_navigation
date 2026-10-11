#!/usr/bin/env python3
"""
S01_plot_plv_group_significance.py
==================================

Publication-style PLV group-significance figure.

This is a structured rewrite of the historical Figure6A.py / Figure6B.py.
The visual layout intentionally remains very close to the original figure,
while data loading and ROI-label handling are generalized.

Input:
    results/plv/pairs/<PAIR>/<SEGMENTATION_MODE>/group_statistics.csv

Output:
    figures/plv/<SEGMENTATION_MODE>/group_significance/

Default statistic:
    -log10(p_empirical)

One figure is generated for corridor and one for junction, matching the old
Figure6A / Figure6B separation.

The plotted significance is fixed to the new empirical Monte-Carlo p-value.

Both order-aligned and terminal-aligned displays exclude right hippocampal ROIs
and use one separate horizontal legend shared by the corridor and junction
figures in each output directory. Analysis inputs are unchanged.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from functions.plv_plot_style import (
    FREQUENCY_BANDS,
    ROI1_ORDER_BY_PROCESS,
    ROI2_ORDER,
    create_publication_grid,
    publication_figsize,
    roi_label,
    set_publication_style,
)
from functions.plv_plot_palette import PLV_BAND_COLORS


PYTHON_CODE_ROOT = next(
    parent
    for parent in Path(__file__).resolve().parents
    if (parent / "project_config.py").exists()
)

if str(PYTHON_CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_CODE_ROOT))

from project_config import CFG


P_VALUE_THRESHOLD = 0.05
MIN_Y_MAX = 6.5

FONT_SIZE = 12
SAVE_DPI = 600

# Historical Figure-6 behavior for an absent ROI/segment/band cell was to put
# the line at zero (-log10(1) = 0). Keep this default for visual continuity.
MISSING_VALUE_AS_ZERO = True



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot publication-style PLV group-significance figures."
    )

    parser.add_argument(
        "--segmentation-mode",
        default="terminal_aligned",
        choices=["order_aligned", "terminal_aligned"],
        help="Formal navigation-segmentation scheme.",
    )

    parser.add_argument(
        "--process-type",
        default="both",
        choices=["corridor", "junction", "both"],
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=SAVE_DPI,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output override for a read-only plotting preview.",
    )

    return parser.parse_args()


def load_all_pair_stats(
    project_root: Path,
    segmentation_mode: str,
) -> pd.DataFrame:
    pair_root = project_root / "results" / "plv" / "pairs"

    files = sorted(
        pair_root.glob(
            f"*/{segmentation_mode}/group_statistics.csv"
        )
    )

    if not files:
        raise FileNotFoundError(
            "No group_statistics.csv files found under:\n"
            f"{pair_root}/*/{segmentation_mode}/"
        )

    frames = []

    for file in files:
        frame = pd.read_csv(file)
        frame["source_file"] = str(file)
        frames.append(frame)

    data = pd.concat(
        frames,
        ignore_index=True,
    )

    required = {
        "roi1",
        "roi2",
        "process_type",
        "segment_index",
        "segment_label",
        "band_low",
        "band_high",
    }

    missing = sorted(required - set(data.columns))

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
        )

    return data


def get_segment_table(
    process_data: pd.DataFrame,
) -> pd.DataFrame:
    segment_table = (
        process_data[
            ["segment_index", "segment_label"]
        ]
        .drop_duplicates()
        .sort_values("segment_index")
        .reset_index(drop=True)
    )

    return segment_table


def p_to_plot_value(
    p: float,
) -> float:
    if not np.isfinite(p):
        return 0.0 if MISSING_VALUE_AS_ZERO else np.nan

    if p <= 0:
        # The canonical +1 empirical p should never be zero. Keep a defensive
        # plotting clip so malformed external tables cannot create infinity.
        p = 1e-12

    return -np.log10(p)


def get_line_values(
    pair_data: pd.DataFrame,
    segment_indices: np.ndarray,
    low: float,
    high: float,
    p_column: str,
) -> np.ndarray:

    band_data = pair_data[
        (pair_data["band_low"] == low)
        & (pair_data["band_high"] == high)
    ]

    values = []

    for segment_index in segment_indices:

        row = band_data[
            band_data["segment_index"] == segment_index
        ]

        if row.empty:
            values.append(
                0.0 if MISSING_VALUE_AS_ZERO else np.nan
            )
            continue

        if len(row) > 1:
            raise ValueError(
                "Duplicate group-statistic rows found for "
                f"segment={segment_index}, band={low}-{high}."
            )

        p = float(row.iloc[0][p_column])
        values.append(p_to_plot_value(p))

    return np.asarray(values, dtype=float)


def panel_ymax(
    line_values: list[np.ndarray],
) -> float:
    finite = np.concatenate(
        [
            x[np.isfinite(x)]
            for x in line_values
            if len(x) > 0
        ]
    )

    if finite.size == 0:
        return MIN_Y_MAX

    panel_max = float(np.max(finite))

    if panel_max > MIN_Y_MAX:
        return float(int(panel_max) + 1)

    return MIN_Y_MAX


def displayed_navigation_rois(segmentation_mode: str) -> list[str]:
    """Use the same display filter for both modes; keep analysis ROIs intact."""
    return [roi for roi in ROI2_ORDER if not roi.startswith("Hippocampus R")]


def save_shared_legend(output_dir: Path, dpi: int) -> None:
    """One horizontal frequency/threshold key for corridor and junction figures."""
    set_publication_style(fontsize=FONT_SIZE)
    handles = [
        Line2D([], [], color=PLV_BAND_COLORS.get((low, high), "#4d4d4d"),
               marker="o", linewidth=1.45, markersize=4.5,
               label=f"{low:g}–{high:g} Hz")
        for low, high in FREQUENCY_BANDS
    ]
    handles.append(Line2D([], [], color="black", linestyle="dashed",
                          label="p = 0.05"))
    fig = plt.figure(figsize=(9.0, 0.42))
    fig.legend(handles=handles, loc="center", ncol=len(handles), frameon=False,
               fontsize=11, handlelength=2.0, handletextpad=0.55,
               columnspacing=1.35, borderaxespad=0)
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        path = output_dir / f"plv_group_significance_shared_legend.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.04,
                    transparent=True)
        print(f"Saved: {path}")
    plt.close(fig)


def plot_process(
    data: pd.DataFrame,
    process_type: str,
    segmentation_mode: str,
    p_column: str,
    output_dir: Path,
    dpi: int,
) -> None:

    process_data = data[
        data["process_type"] == process_type
    ].copy()

    if process_data.empty:
        print(
            f"[SKIP] No rows for process_type={process_type}"
        )
        return

    # Fixed publication grid, independent of which result files exist.
    roi1s = list(ROI1_ORDER_BY_PROCESS[process_type])
    roi2s = displayed_navigation_rois(segmentation_mode)

    nrows = len(roi1s)
    ncols = len(roi2s)

    if nrows == 0 or ncols == 0:
        print(f"[SKIP] Empty ROI grid for {process_type}")
        return

    segment_table = get_segment_table(process_data)

    segment_indices = (
        segment_table["segment_index"]
        .to_numpy()
    )

    segment_labels = (
        segment_table["segment_label"]
        .astype(str)
        .tolist()
    )

    x_positions = np.arange(
        len(segment_indices),
        dtype=float,
    )

    xlim = (
        -0.5,
        len(segment_indices) - 0.5,
    )

    xlabel = (
        "Corridor number"
        if process_type == "corridor"
        else "Junction number"
    )

    figsize = publication_figsize(
        nrows=nrows,
        ncols=ncols,
    )

    plt.close("all")

    fig, axes = create_publication_grid(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        fontsize=FONT_SIZE,
    )

    threshold_y = -np.log10(P_VALUE_THRESHOLD)

    for row_i, roi1 in enumerate(roi1s):

        for col_i, roi2 in enumerate(roi2s):

            ax = axes[row_i, col_i]

            pair_data = process_data[
                (process_data["roi1"] == roi1)
                & (process_data["roi2"] == roi2)
            ].copy()

            if pair_data.empty:
                # Preserve grid geometry but show no frame for unavailable pairs.
                ax.axis("off")
                continue

            line_values = []

            for low, high in FREQUENCY_BANDS:

                y = get_line_values(
                    pair_data=pair_data,
                    segment_indices=segment_indices,
                    low=low,
                    high=high,
                    p_column=p_column,
                )

                line_values.append(y)

                # Keep fractional band edges (for example 0.5-4 Hz) intact so
                # they resolve to the categorical color defined in the palette.
                band_key = (float(low), float(high))
                band_color = PLV_BAND_COLORS.get(band_key, "#4d4d4d")

                ax.plot(
                    x_positions,
                    y,
                    marker="o",
                    color=band_color,
                    linewidth=1.45,
                    markersize=4.5,
                    label=f"[{low}, {high}]",
                )

            # The historical code drew the same threshold repeatedly.
            # Draw it once here: visually identical, cleaner legend.
            ax.plot(
                x_positions,
                np.full(
                    len(x_positions),
                    threshold_y,
                    dtype=float,
                ),
                linestyle="dashed",
                color="black",
                label="p=0.05",
            )

            ax.set_xlim(xlim)
            ax.set_ylim(
                0,
                panel_ymax(line_values),
            )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(segment_labels)

            # Only the bottom row keeps x tick labels and x-axis title.
            if row_i == nrows - 1:
                ax.set_xlabel(
                    xlabel,
                    fontsize=FONT_SIZE + 2,
                    fontweight="bold",
                )
            else:
                ax.set_xticks([])

            # ROI2 labels are the column titles.
            if row_i == 0:
                ax.set_title(
                    roi_label(roi2),
                    fontsize=FONT_SIZE + 2,
                    pad=13,
                    fontweight="bold",
                )

            # ROI1 labels are vertically centered on the left margin.
            if col_i == 0:
                ax.text(
                    -0.34,
                    0.5,
                    roi_label(roi1),
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontweight="bold",
                    fontsize=FONT_SIZE + 1,
                    clip_on=False,
                )

                ax.set_ylabel(
                    "-log10(p)",
                    fontsize=FONT_SIZE + 2,
                    fontweight="bold",
                )

    # Each segmentation mode has a separately saved shared horizontal legend.
    fig.tight_layout(
        pad=1,
        rect=(0.0, 0.0, 0.99, 1),
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    stem = (
        f"{process_type}_"
        f"{segmentation_mode}_"
        f"plv_group_significance_"
        f"{p_column}"
    )

    png_file = output_dir / f"{stem}.png"
    pdf_file = output_dir / f"{stem}.pdf"

    fig.savefig(
        png_file,
        dpi=dpi,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")


def main() -> None:
    args = parse_args()

    project_root = CFG.root

    data = load_all_pair_stats(
        project_root=project_root,
        segmentation_mode=args.segmentation_mode,
    )

    if "p_empirical" not in data.columns:
        raise ValueError(
            "Required p column is absent: p_empirical"
        )

    output_dir = args.output_dir if args.output_dir is not None else (
        project_root
        / "figures"
        / "plv"
        / args.segmentation_mode
        / "group_significance"
    )

    process_types = (
        ["corridor", "junction"]
        if args.process_type == "both"
        else [args.process_type]
    )

    for process_type in process_types:
        plot_process(
            data=data,
            process_type=process_type,
            segmentation_mode=args.segmentation_mode,
            p_column="p_empirical",
            output_dir=output_dir,
            dpi=args.dpi,
        )

    save_shared_legend(output_dir=output_dir, dpi=args.dpi)


if __name__ == "__main__":
    main()
