"""Shared project configuration loader for the SEEG navigation project.

The cross-language runtime project configuration is:
    config/project_config.json

The human-editable source that generates it is:
    config/project_config.py

Analysis-specific JSON files live beside it and can be read with
``load_analysis_config``.  Project-internal paths in JSON are relative to the
project root so the same code works on the server and on another machine after
setting ``SEEG_NAVIGATION_ROOT``.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

_DEFAULT_PROJECT_ROOT = Path(
    "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
)

PROJECT_ROOT = Path(
    os.environ.get("SEEG_NAVIGATION_ROOT", str(_DEFAULT_PROJECT_ROOT))
).expanduser().resolve()

CONFIG_DIR = PROJECT_ROOT / "config"
PROJECT_CONFIG_FILE = CONFIG_DIR / "project_config.json"


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found:\n{path}")
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


PROJECT_SETTINGS = _read_json(PROJECT_CONFIG_FILE)


class ModulePaths(SimpleNamespace):
    """Root plus convenience access to module-specific result/figure/log dirs."""

    def __init__(self, root: Path):
        super().__init__(root=root)

    def module(self, module: str) -> Path:
        return self.root / module

    def __getattr__(self, module: str) -> Path:
        if module.startswith("_"):
            raise AttributeError(module)
        return self.module(module)


def _build_path_namespace(section: dict[str, str]) -> SimpleNamespace:
    return SimpleNamespace(
        **{key: _resolve_path(value) for key, value in section.items()}
    )


@dataclass(frozen=True)
class ProjectConfig:
    root: Path
    config_dir: Path
    sampling_rate_hz: float
    dataset_modes: tuple[str, ...]
    process_types: tuple[str, ...]
    condition_order: tuple[str, ...]
    conditions: dict[str, Any]
    segment_windows_sec: dict[str, Any]
    data: SimpleNamespace
    files: SimpleNamespace
    results: ModulePaths
    figures: ModulePaths
    logs: ModulePaths
    roi_abbreviations: dict[str, str]
    ap_split_roi_definition: dict[str, Any]
    external: SimpleNamespace


_path_values = PROJECT_SETTINGS["paths"]
_data_namespace = SimpleNamespace(
    root=_resolve_path(_path_values["data"]),
    processed=_resolve_path(_path_values["processed"]),
    trial_timeseries=_resolve_path(_path_values["trial_timeseries"]),
    trial_timeseries_by_turn_number=_resolve_path(
        _path_values["trial_timeseries_by_turn_number"]
    ),
    derived=_resolve_path(_path_values["derived"]),
    derived_features=_resolve_path(_path_values["derived_features"]),
    behavior=_resolve_path(_path_values["behavior"]),
    encoding=_resolve_path(_path_values["encoding"]),
    metadata=_resolve_path(_path_values["metadata"]),
    feature_selection=_resolve_path(_path_values["feature_selection"]),
)

CFG = ProjectConfig(
    root=PROJECT_ROOT,
    config_dir=CONFIG_DIR,
    sampling_rate_hz=float(PROJECT_SETTINGS["sampling_rate_hz"]),
    dataset_modes=tuple(PROJECT_SETTINGS["dataset_modes"]),
    process_types=tuple(PROJECT_SETTINGS["process_types"]),
    condition_order=tuple(PROJECT_SETTINGS["condition_order"]),
    conditions=PROJECT_SETTINGS["conditions"],
    segment_windows_sec=PROJECT_SETTINGS["segment_windows_sec"],
    data=_data_namespace,
    files=_build_path_namespace(PROJECT_SETTINGS["files"]),
    results=ModulePaths(_resolve_path(_path_values["results"])),
    figures=ModulePaths(_resolve_path(_path_values["figures"])),
    logs=ModulePaths(_resolve_path(_path_values["logs"])),
    roi_abbreviations={
        item["full_name"]: item["abbreviation"]
        for item in PROJECT_SETTINGS["roi_abbreviations"]
    },
    ap_split_roi_definition=dict(PROJECT_SETTINGS["ap_split_roi_definition"]),
    external=_build_path_namespace(PROJECT_SETTINGS.get("external", {})),
)


def load_analysis_config(filename: str) -> dict[str, Any]:
    """Load an analysis-specific JSON from ``config/``.

    ``filename`` may be supplied with or without the ``.json`` suffix.
    """
    if not filename.endswith(".json"):
        filename = f"{filename}.json"
    return _read_json(CONFIG_DIR / filename)


def roi_abbreviation(roi_name: str) -> str:
    """Return the shared project abbreviation, falling back to the full name."""
    name = str(roi_name)
    return CFG.roi_abbreviations.get(name, name)
