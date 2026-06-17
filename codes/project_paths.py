"""Project-wide paths for the shared plotting code."""

from __future__ import annotations

import re
from pathlib import Path

CODES_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CODES_DIR.parent

DATA_DIR = PROJECT_ROOT / "data"
FIGURE_DIR = PROJECT_ROOT / "figures"
LOG_DIR = PROJECT_ROOT / "logs"

BEHAVIOR_DIR = DATA_DIR / "behavior"
LMM_DIR = DATA_DIR / "LMM"
LASSO_DIR = DATA_DIR / "Lasso"
AAL3_DIR = DATA_DIR / "AAL3_data"
SPECTRUM_DIR = DATA_DIR / "spectrum"


def ensure_project_dirs() -> None:
    """Create shared output folders if they do not exist."""
    for directory in [DATA_DIR, FIGURE_DIR, LOG_DIR, BEHAVIOR_DIR, LMM_DIR, LASSO_DIR, AAL3_DIR, SPECTRUM_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def get_figure_group(script_file) -> str:
    """Return the output group folder name for a figure script."""
    stem = Path(script_file).stem
    main_match = re.match(r"Figure(\d+)", stem, flags=re.IGNORECASE)
    if main_match:
        return f"Figure{main_match.group(1)}"

    sup_match = re.match(r"Sup_fig(\d+)", stem, flags=re.IGNORECASE)
    if sup_match:
        return f"SupFigure{sup_match.group(1)}"

    return "Other"


def get_figure_output_dir(script_file) -> Path:
    """Return and create the output directory for one figure script."""
    ensure_project_dirs()
    out_dir = FIGURE_DIR / get_figure_group(script_file)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


ensure_project_dirs()
