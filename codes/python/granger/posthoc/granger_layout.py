"""Resolve the project's language-specific layout or a standalone source copy."""
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
SPLIT_LAYOUT = PYTHON_ROOT.name == 'granger' and PYTHON_ROOT.parent.name == 'python' and PYTHON_ROOT.parent.parent.name == 'codes'
PROJECT_ROOT = PYTHON_ROOT.parents[2] if SPLIT_LAYOUT else PYTHON_ROOT
MATLAB_ROOT = PROJECT_ROOT / 'codes/matlab/granger' if SPLIT_LAYOUT else PYTHON_ROOT / 'matlab'
CONFIG_PATH = PROJECT_ROOT / 'config/granger_config.json'
ARCHIVE_ROOT = PROJECT_ROOT / 'archive/granger_organization' if SPLIT_LAYOUT else PYTHON_ROOT / 'archive'


def source_path(relative):
    """Resolve paths recorded by the original organization manifest."""
    path = Path(relative)
    if path.parts[0] == 'matlab':
        return MATLAB_ROOT.joinpath(*path.parts[1:])
    if path.parts[0] == 'config':
        return PROJECT_ROOT / path
    if path.parts[0] == 'archive':
        return ARCHIVE_ROOT.joinpath(*path.parts[1:])
    return PYTHON_ROOT / path
