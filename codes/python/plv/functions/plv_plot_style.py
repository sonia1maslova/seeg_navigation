"""
plv_plot_style.py
=================

Shared plotting configuration for PLV figures.

The visual language intentionally follows the original publication Figure 6:
- Arial
- bold axis labels / panel titles
- outward ticks
- no top/right spines
- ROI1 labels on the left side of the grid
- ROI2 labels as column titles
- frequency bands as line + circle markers
- dashed black p=0.05 reference line

The implementation is more robust than the historical `set_plot()` helper:
it always returns a 2-D axes array and centralizes ROI abbreviations.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PYTHON_CODE_ROOT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "project_config.py").exists()
)
if str(PYTHON_CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_CODE_ROOT))

from project_config import CFG, load_analysis_config

PLV_CONFIG = load_analysis_config("plv_config")


# ---------------------------------------------------------------------
# Shared project/analysis labels and ordering
# ---------------------------------------------------------------------
ROI_LABELS = dict(CFG.roi_abbreviations)
ROI1_ORDER_BY_PROCESS = {
    process_type: list(rois)
    for process_type, rois in PLV_CONFIG["cortical_rois"].items()
}
ROI2_ORDER = list(PLV_CONFIG["navigation_rois"])
FREQUENCY_BANDS = [
    tuple(float(v) for v in band)
    for band in PLV_CONFIG["frequency_bands_hz"]
]


def roi_label(roi: str) -> str:
    """Return publication-style ROI abbreviation."""
    return ROI_LABELS.get(str(roi), str(roi))


def order_rois(values, preferred_order):
    """Order known ROIs canonically, then append unknown ROIs alphabetically."""
    values = list(dict.fromkeys(str(x) for x in values))

    preferred = [roi for roi in preferred_order if roi in values]
    remaining = sorted(set(values) - set(preferred))

    return preferred + remaining


def publication_figsize(nrows: int, ncols: int) -> tuple[float, float]:
    """
    Preserve the historical Figure-6 size rule.

    Original:
        (5 + 4*ncols/2.54, (3 + 4*nrows)/2.54)
    """
    return (
        5 + 4 * ncols / 2.54,
        (3 + 4 * nrows) / 2.54,
    )


def set_publication_style(
    fontsize: float = 12,
    fontfamily: str = "Arial",
    axis_linewidth: float = 1.8,
) -> None:
    """Set the global matplotlib style used by the PLV publication figure."""

    plt.rcParams.update(
        {
            "font.family": fontfamily,
            "font.size": fontsize,
            "axes.labelsize": fontsize,
            "axes.titlesize": fontsize + 2,
            "xtick.labelsize": fontsize - 2,
            "ytick.labelsize": fontsize - 2,
            "legend.fontsize": fontsize - 1,
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "axes.linewidth": axis_linewidth,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(
    ax,
    axis_linewidth: float = 1.8,
) -> None:
    """Apply the original outward-tick / two-spine axis style."""

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=4,
        width=1.2,
        bottom=True,
        top=False,
        left=True,
        right=False,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)

    ax.spines["bottom"].set_linewidth(axis_linewidth)
    ax.spines["left"].set_linewidth(axis_linewidth)


def create_publication_grid(
    nrows: int,
    ncols: int,
    figsize=None,
    dpi: int = 100,
    fontsize: float = 12,
    fontfamily: str = "Arial",
    axis_linewidth: float = 1.8,
):
    """
    Create a Figure-6-style grid and always return a 2-D axes array.
    """

    if figsize is None:
        figsize = publication_figsize(nrows, ncols)

    set_publication_style(
        fontsize=fontsize,
        fontfamily=fontfamily,
        axis_linewidth=axis_linewidth,
    )

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        dpi=dpi,
        squeeze=False,
    )

    for ax in axes.ravel():
        style_axis(
            ax,
            axis_linewidth=axis_linewidth,
        )

    return fig, axes
