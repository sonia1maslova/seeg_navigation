# Recent analysis and figure code snapshot (2026-10-11)

This update supplements the existing figure-reproduction repository. Earlier `codes/Figure*.py`, data-manifest files and other repository content are retained. It is not a complete replacement of the original repository or a release of participant-level data.

## Separate repositories

- Analysis and plotting source: https://github.com/sonia1maslova/seeg_navigation
- Formal figure gallery: https://github.com/sonia1maslova/seeg-navigation-figure-gallery
- Gallery website: https://sonia1maslova.github.io/seeg-navigation-figure-gallery/

The gallery contains rendered figures, not the analysis source or raw recordings.

## Scope of this source snapshot

- `codes/python/pac/by_lina/`: PAC–behaviour, MI-existence forest plots, the optional MI-null-density preview, and 1–150 Hz aperiodic-exponent coefficient plots.
- `codes/r/pac/by_lina/`: existing Lina-derived PAC segment/control models, separate reporting/FDR procedures, and the supplementary condition-adjusted aperiodic-exponent models. Lina's original input results remain read-only.
- `codes/python/power/encoding_lmm/` and `codes/r/power/encoding_lmm/`: Power encoding LMM and compact forest-plot reporting.
- `codes/python/granger/posthoc/`: current compact, two-direction group-P/significant-channel-pair-proportion plots and their configuration helpers.
- `codes/python/pac/test/`: separate calculation and plotting scripts for the simulated circular-vector-averaging equivalence check.

The source and dependency list is recorded in `docs/source_snapshot_20261011.json`. No trial-level tables, contact lists, individual identifiers, raw sEEG, fitted models, posterior samples, permutation arrays, private keys, tokens, or local work-in-progress folders are added by this update.

## Runtime and data requirements

Use the project's existing analysis environment. Python plotting requires NumPy, pandas and Matplotlib. R requirements are declared at the beginning of each script; they include lme4, dplyr, readr, jsonlite, digest and ggplot2 as applicable. These are requirements, not a newly validated environment lockfile.

Set `SEEG_NAVIGATION_ROOT` to the project root and make the required existing result tables available at the paths checked by each script. The shared project configuration retains the standard ROI abbreviations; external data/toolbox locations in the published configuration are portable relative placeholders and must be configured for an authorized local environment.

**Server-specific limitation:** `run_aperiodic_segment_effect_adjusted.R` additionally verifies the original Lina exponent table at its recorded NFS location against the saved snapshot. Archiving that script does not make it independently runnable without the protected inputs and source-location checks. No model was refitted just to publish this update.

Some archived scripts retain original server-path defaults for provenance. They are not credentials, and do not grant access to protected data. Do not use the public repository as a destination for private data.

For each Granger branch, choose the matching `config/granger_formal_20261003_<branch>.json`. The published `config/granger_config.json` defaults to the terminal-aligned formal configuration. Set the following environment variables before running the group plotting script:

```sh
export GRANGER_CONFIG="$SEEG_NAVIGATION_ROOT/config/granger_formal_20261003_<branch>.json"
export GRANGER_RESULTS_ROOT="$SEEG_NAVIGATION_ROOT/results/granger/formal_20261003"
export GRANGER_POSTHOC_RESULTS_ROOT="$SEEG_NAVIGATION_ROOT/results/granger/formal_20261003/posthoc/<branch>"
export GRANGER_FORMAL_FIGURE_ROOT="$SEEG_NAVIGATION_ROOT/figures/granger/reproduced_formal_20261003"
python codes/python/granger/posthoc/group/plot_group_granger_p_and_sig_pair_proportion.py
```

Replace `<branch>` with `order_aligned` or `order_aligned_and_terminal_aligned`. Use a fresh output root: the formal plotting entry point intentionally refuses to overwrite existing named outputs. Required group summaries and channel-pair CSVs are protected external inputs. The retained `gc_params.subject_summary_method=median` is a legacy auxiliary S1 summary setting; formal posthoc aggregation is `mean` at both subject and group levels, as recorded separately in `posthoc`.

## Current plotting/reporting conventions

- Aperiodic exponent: four ROI labels on the **y axis**, original beta and horizontal 95% confidence intervals on the **x axis**; scientific offset is `×10^-3`. The two terminal-inclusion branches remain separate. Raw-P stars are retained; the three hippocampal effects are BH-corrected within each branch, while the single aMTG effect is uncorrected.
- Power forest plots: three hippocampal subregions by five bands form a 15-test family; every other ROI has its own five-band family, separately by predictor, scope and terminal-inclusion branch.
- PAC–behaviour: aMTG has a six-combination family and the three hippocampal subregions a combined 18-combination family.
- Lina MI-existence plots reuse saved pooled-segment statistics and the 1,000 saved group surrogates. The original `all_dynamic_split` scope excludes terminal corridor only, **not** terminal junction. It must not be relabelled as the newer exclude-both-terminal-segments analysis.
- Adjusted PAC segment/control models use the reporting policy in their own README and script: six phase–amplitude combinations within each ROI, with branches/ranges kept separate.
- Granger: each metric panel contains both directions, filled and hollow square markers, thin translucent lines and unhatched paired proportion bars. A separate direction legend is retained. Both alignment branches and time-domain/five-band metrics are retained. No new FDR is applied in this update.

The MI-density figure is still a preview and is not part of this formal-gallery publication. A saved source script does not imply that its preview is a final publication figure.

Raw-P stars use `* P < 0.05`, `** P < 0.01`, `*** P < 0.001`; a dagger indicates the relevant BH-FDR result where used. Exact figure-specific rules remain in each output directory's existing TXT/audit files.

## Publication integrity

The code update and gallery update are separate commits in separate repositories. Gallery publication adds/replaces selected formal PNG files and merges their metadata into the existing index; it does not delete historical gallery content or replace the website application. The current analysis project and existing statistical results are not overwritten by publication.
