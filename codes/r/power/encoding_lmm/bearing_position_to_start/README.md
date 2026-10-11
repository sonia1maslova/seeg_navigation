# Position-to-start bearing: Power sensitivity branch

This branch reverses both `Angle2Xaxis` and `Angle2orientation` by `(x + pi) %% (2*pi)` while retaining the raw-radian linear-predictor model. It does not convert the model to sine/cosine coding. Original data and results are preserved.

## Entry points

- `run_power_bearing_position_to_start.R [project_root]`
- `plot_power_bearing_position_to_start.R [project_root]`
- `bearing_branch_common.R`: shared helpers, not a standalone entry point.

The project root can instead be supplied with `SEEG_NAVIGATION_ROOT`. R packages are checked but never installed. Run in the existing `r_env` environment.

## Scientific scope

- Reuse existing ROI/contact-averaged `data/derived/power/segment_power/` values, not neural signal extraction.
- Import only explicitly allowed exact S01 preparation/S02 fitting function and constant assignments via R syntax-tree parsing. Never source those scripts' original main bodies.
- Same canonical five fitted ROIs, five frequency bands, raw angle units, raw power dB, condition fixed effects, subject random intercept, segment-type covariate for pooled scopes, ML, bobyqa and numerical checks.
- `Angle2orientation`: corridor only. `Angle2Xaxis`: corridor, junction and all_segments. Both include_last/exclude_last are retained. Terminal filtering uses the original S01 flags, with no new renumbering.
- 200 reversed-coding models plus 200 original-coding reproduction models. The reproduced originals are saved only inside this new branch, not over the original result.
- Formal angle heatmaps retain the established four-ROI display and FDR families: aMTG across 5 frequency bands; three hippocampal ROIs jointly across 15 tests. Families are separate for angle, scope and include/exclude branch. pMTG retains raw results, q is undefined, and no new family is silently invented.

## Paths

Deploy the three R files to `codes/r/power/encoding_lmm/bearing_position_to_start/`.

New encoding input:

- `data/encoding/corridor_subsegi_encoding_table_position_to_start.csv`
- `data/encoding/turn_encoding_table_position_to_start.csv`

Output roots:

- `data/derived/power/encoding_lmm/bearing_position_to_start/`
- `results/power/encoding_lmm/bearing_position_to_start/`
- `figures/power/encoding_lmm/bearing_position_to_start/`

Results contain `comparison_with_original.csv` and `comparison_with_reproduced_original.csv`, summaries and model QC. `qc/original_result_reproduction_check.csv` identifies numerical or sample-count differences from stored original results; these must not be misinterpreted as bearing-direction effects. CSV estimates, raw p, q and sample counts are retained. Figures contain 8 PNG heatmaps and their plotted values, plus an FDR readme; no PDF/SVG.

## Safeguards

Before new outputs are created, the runner checks the encoding reversal, original keys, unchanged power/condition/other covariates, and agreement of reconstructed original inputs with stored prepared inputs. Model QC independently checks identical old/new sample counts. MD5 records verify that audited original inputs and canonical scripts did not change during the run.

Existing supplemental output directories cause a stop: no automatic deletion, retry overwrite or cleanup. `_SUCCESS` means the requested model loop finished and input preservation checks passed, not that every model is non-singular or statistically significant. Inspect `model_qc.csv` and the reproduction audit before interpretation.
