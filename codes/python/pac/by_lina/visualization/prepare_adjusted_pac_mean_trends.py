#!/usr/bin/env python3
"""Prepare subject-balanced, descriptive PAC means using the adjusted LMM rows.

No PAC is recomputed and no model or significance test is fitted here. Each
input row is already the equal-contact ROI mean for one subject/trial/segment.
Trials are averaged within each subject/position, then subjects are averaged
equally; error bars are the between-subject SEM, as in Lina's original plots.

The include branch has up to 20 chronological positions. Excluding each
trial's terminal corridor (both halves) AND terminal junction leaves up to
17. The old 18-position display excluded only the terminal corridor; it must
not be used to truncate the adjusted include branch.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd


BRANCHES = ("include_last", "exclude_terminal_corridor_and_junction")
ROIS = (
    "Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
    "Hippocampus L_Posterior", "Hippocampus R_Anterior",
)
ROW_KEYS = ["sub_id", "trial", "roi", "phase_band", "amp_band", "segment_id"]
GROUP_KEYS = ["roi", "phase_band", "amp_band", "segment_position", "segment_id", "segment_type"]


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_same_rows(left: pd.DataFrame, right: pd.DataFrame, label: str) -> float:
    for frame in (left, right):
        if frame.duplicated(ROW_KEYS).any():
            raise ValueError(f"{label}: duplicate subject/trial/ROI/combination/segment keys")
    joined = left[ROW_KEYS + ["raw_PAC", "segment_position"]].merge(
        right[ROW_KEYS + ["raw_PAC", "segment_position"]], on=ROW_KEYS,
        how="outer", validate="one_to_one", indicator=True, suffixes=("_left", "_right"),
    )
    if not joined["_merge"].eq("both").all():
        raise ValueError(f"{label}: differing row keys: {joined['_merge'].value_counts().to_dict()}")
    if not np.array_equal(joined.segment_position_left, joined.segment_position_right):
        raise ValueError(f"{label}: segment positions were changed")
    np.testing.assert_allclose(joined.raw_PAC_left, joined.raw_PAC_right, rtol=1e-12, atol=1e-15)
    return float(np.max(np.abs(joined.raw_PAC_left - joined.raw_PAC_right)))


def terminal_mask(frame: pd.DataFrame) -> pd.Series:
    parsed = frame.segment_id.str.extract(r"^([CJ])(\d+)(?:_(first|second))?$")
    if parsed.isna().iloc[:, :2].any().any():
        raise ValueError("Unexpected segment identifier")
    number = parsed[1].astype(int)
    expected_position = np.where(
        parsed[0].eq("J"), 3 * number,
        3 * number - np.where(parsed[2].eq("first"), 2, 1),
    )
    if not np.array_equal(expected_position, frame.segment_position.to_numpy()):
        raise ValueError("Original chronological segment positions do not match segment identifiers")
    return (parsed[0].eq("C") & number.eq(frame.turns + 1)) | (parsed[0].eq("J") & number.eq(frame.turns))


def aggregate(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    subject = frame.groupby(GROUP_KEYS + ["sub_id"], as_index=False, sort=True).agg(
        mean_raw_PAC=("raw_PAC", "mean"), n_trials=("trial", "nunique"), n_trial_rows=("trial", "size"),
    )
    if not subject.n_trials.eq(subject.n_trial_rows).all():
        raise ValueError("A subject/position has duplicate trial rows")
    subject["subject_group_weight"] = 1 / subject.groupby(GROUP_KEYS).sub_id.transform("nunique")
    grouped = subject.groupby(GROUP_KEYS, as_index=False, sort=True).agg(
        mean_raw_PAC=("mean_raw_PAC", "mean"), sem_raw_PAC=("mean_raw_PAC", "sem"),
        n_subjects=("sub_id", "nunique"), n_trial_rows=("n_trial_rows", "sum"),
        min_trials_per_subject=("n_trials", "min"), max_trials_per_subject=("n_trials", "max"),
    )
    # Independent numerical checks, not another call to the groupby aggregator.
    checks = grouped.set_index(GROUP_KEYS)
    for key, rows in subject.groupby(GROUP_KEYS, sort=True):
        values = rows.mean_raw_PAC.to_numpy(float)
        observed = checks.loc[key]
        np.testing.assert_allclose(observed.mean_raw_PAC, np.sum(values) / len(values), atol=1e-15, rtol=1e-12)
        expected_sem = np.std(values, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else np.nan
        np.testing.assert_allclose(observed.sem_raw_PAC, expected_sem, atol=1e-15, rtol=1e-12, equal_nan=True)
        np.testing.assert_allclose(rows.subject_group_weight.sum(), 1, atol=1e-15)
    return subject, grouped


def check_original_means(root: Path, sources: pd.DataFrame) -> dict:
    """Reproduce all 12 existing mean/SEM tables without changing old files."""
    legacy = sources[sources.scope.eq("all_dynamic_split")].copy()
    _, summary = aggregate(legacy)
    max_mean, max_sem, n_files = 0.0, 0.0, 0
    original_root = root / "figures/pac/by_lina/pac_by_segment"
    # Future reruns after the reorganization locate the same preserved files.
    if not original_root.exists():
        original_root = root / "figures/pac/by_lina/pac_segment_effect/mean_trends/lina_original"
    for path in sorted(original_root.rglob("*_values.csv")):
        alias = path.parent.name
        if alias == "aMTG":
            roi = ROIS[0]
        elif alias == "Hipp_L_Ant":
            roi = ROIS[1]
        else:
            continue
        matching = summary[summary.roi.eq(roi)]
        combos = matching[["phase_band", "amp_band"]].drop_duplicates()
        combo = [(p, a) for p, a in combos.itertuples(index=False, name=None) if f"_{p}_{a}_values" in path.stem]
        if len(combo) != 1:
            raise ValueError(f"Cannot determine combination for {path}")
        p, a = combo[0]
        old = pd.read_csv(path)
        actual = matching[matching.phase_band.eq(p) & matching.amp_band.eq(a)]
        merged = old.merge(actual, on=["segment_position", "segment_id"], validate="one_to_one", suffixes=("_old", "_new"))
        if len(merged) != len(old) or len(merged) != len(actual):
            raise ValueError(f"Original graph position mismatch: {path}")
        for metric in ("mean_raw_PAC", "sem_raw_PAC"):
            np.testing.assert_allclose(merged[f"{metric}_old"], merged[f"{metric}_new"], atol=1e-15, rtol=1e-10)
        max_mean = max(max_mean, float(np.max(np.abs(merged.mean_raw_PAC_old - merged.mean_raw_PAC_new))))
        max_sem = max(max_sem, float(np.max(np.abs(merged.sem_raw_PAC_old - merged.sem_raw_PAC_new))))
        n_files += 1
    return {"files_verified": n_files, "max_absolute_mean_difference": max_mean,
            "max_absolute_sem_difference": max_sem,
            "aggregation": "equal trial mean within each subject/position, then equal subject mean and between-subject SEM"}


def prepare(root: Path, output: Path) -> dict:
    model_root = root / "results/pac/by_lina/segment_effect_adjusted"
    inputs = {branch: model_root / branch / "main_pac/analysis_dataset.csv" for branch in BRANCHES}
    frames = {branch: pd.read_csv(path) for branch, path in inputs.items()}
    original_paths = (
        root / "results/pac/by_lina/formal_raw_pac/02_hpc_path_accumulation/analysis_dataset.csv",
        root / "results/pac/by_lina/formal_raw_pac/04_amtg_followup_and_behavior/analysis_dataset_amtg_accumulation.csv",
    )
    sources = pd.concat([pd.read_csv(path) for path in original_paths], ignore_index=True)
    include = frames[BRANCHES[0]]
    source_include = sources[sources.scope.eq("all_dynamic_split_with_terminal_C")]
    source_diff = assert_same_rows(include, source_include, "include vs unmodified Lina source")
    derived_exclude = include.loc[~terminal_mask(include)].copy()
    filter_diff = assert_same_rows(derived_exclude, frames[BRANCHES[1]], "independent trialwise exclusion vs adjusted model rows")
    original_check = check_original_means(root, sources)
    with (root / "config/project_config.json").open(encoding="utf-8") as stream:
        cfg = json.load(stream)
    roi_labels = {item["full_name"]: item["abbreviation"] for item in cfg["roi_abbreviations"]}
    report = {
        "analysis": "Descriptive raw PAC means, not covariate-adjusted predictions or a new LMM",
        "aggregation": original_check["aggregation"], "error_bars": "SEM across subject means (ddof=1)",
        "preserve_segment_position": True, "old_display_position_count": 18,
        "include_vs_lina_source_max_absolute_PAC_difference": source_diff,
        "exclude_trialwise_vs_adjusted_model_max_absolute_PAC_difference": filter_diff,
        "original_figure_reproduction": original_check, "branches": {},
        "source_files": [{"path": str(path), "sha256": sha256(path)} for path in (*inputs.values(), *original_paths)],
    }
    for branch, frame in frames.items():
        if set(frame.roi) != set(ROIS) or frame.duplicated(ROW_KEYS).any():
            raise ValueError(f"Unexpected ROI coverage or duplicate rows: {branch}")
        if not np.isfinite(frame.raw_PAC).all() or (frame.raw_PAC < 0).any():
            raise ValueError(f"Non-finite/negative PAC in {branch}")
        for roi, rows in frame.groupby("roi"):
            if len(rows[["phase_band", "amp_band"]].drop_duplicates()) != 6:
                raise ValueError(f"Expected six combinations for {roi}")
        terminal_mask(frame)  # Validate chronological positions in each branch.
        if branch == BRANCHES[1] and terminal_mask(frame).any():
            raise ValueError("Excluded branch contains terminal corridor/junction rows")
        subject, group = aggregate(frame)
        subject.insert(0, "scope", branch)
        group.insert(0, "scope", branch)
        group["roi_label"] = group.roi.map(roi_labels)
        if group.roi_label.isna().any():
            raise ValueError("A required ROI is missing from the standard abbreviation configuration")
        destination = output / branch
        destination.mkdir(parents=True, exist_ok=True)
        subject.to_csv(destination / "subject_mean_pac.csv", index=False)
        group.to_csv(destination / "group_mean_pac.csv", index=False)
        coverage = frame.groupby(["roi", "phase_band", "amp_band", "turns"], as_index=False).agg(
            n_trial_rows=("raw_PAC", "size"), n_subjects=("sub_id", "nunique"),
            n_positions=("segment_position", "nunique"), max_position=("segment_position", "max"),
        )
        coverage.to_csv(destination / "coverage.csv", index=False)
        for (roi, phase, amp), rows in group.groupby(["roi", "phase_band", "amp_band"], sort=True):
            roi_dir = destination / safe_name(roi)
            roi_dir.mkdir(parents=True, exist_ok=True)
            rows.to_csv(roi_dir / f"pac_by_segment_{safe_name(roi)}_{phase}_{amp}_values.csv", index=False)
        report["branches"][branch] = {
            "input_trial_roi_rows": len(frame), "subject_mean_rows": len(subject), "group_mean_rows": len(group),
            "n_ROIs": frame.roi.nunique(), "n_combinations_per_ROI": 6,
            "n_positions": frame.segment_position.nunique(), "max_position": int(frame.segment_position.max()),
            "n_figures_expected": 24, "unique_trial_roi_segment_keys": True,
            "subject_equal_weights_verified": True, "SEM_independently_verified": True,
        }
    output.mkdir(parents=True, exist_ok=True)
    (output / "mean_trends_manifest.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.getenv("SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")))
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    output = args.output_root or args.project_root / "results/pac/by_lina/pac_segment_effect/mean_trends/adjusted"
    report = prepare(args.project_root, output)
    print(json.dumps({"output_root": str(output), "branches": report["branches"], "original_reproduction": report["original_figure_reproduction"]}, indent=2))


if __name__ == "__main__":
    main()
