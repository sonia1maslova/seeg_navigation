# Independent bearing-direction sensitivity analysis. Reuse existing extracted
# power, exact canonical preparation functions, and exact canonical LMM fitting.
# Only angular predictors are fitted; no original result, table, or plot changes.

entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(entry)) stop("Run with Rscript.", call. = FALSE)
source(file.path(dirname(sub("^--file=", "", entry[1])), "bearing_branch_common.R"))
branch_require(c("dplyr", "readr", "tibble", "lme4", "lmerTest", "jsonlite"))
library(dplyr)
library(readr)
library(tibble)
library(lme4)
library(lmerTest)
cfg <- branch_config()
canonical_dir <- file.path(cfg$codes$r, "power", "encoding_lmm")
prep_file <- file.path(canonical_dir, "S01_prepare_power_encoding_data.R")
model_file <- file.path(canonical_dir, "S02_run_power_encoding_lmm.R")
prep <- branch_load_assignments(prep_file, c("safe_filename", "find_col", "as_numeric_safe",
  "normalize_half", "n_junctions_from_condition", "assert_unique_keys", "read_corridor_encoding",
  "read_junction_encoding", "power_file_for", "prepare_power_common", "prepare_corridor",
  "prepare_junction", "make_all_segments", "summarize_preparation"), "FORMAL_BANDS")
model <- branch_load_assignments(model_file, c("safe_filename", "make_condition_factor",
  "get_convergence_message", "collapse_warnings", "extract_random_effect_info",
  "extract_coef_table", "fit_one_model"), c("SELECTED_ROIS", "FORMAL_BANDS", "BAND_LABELS",
  "USE_REML", "MIN_SUBJECTS", "MIN_OBSERVATIONS", "SINGULAR_TOLERANCE", "OPTIMIZER",
  "MAXFUN", "RAW_SIGNIFICANCE_ALPHA"))
prep$power_root <- file.path(cfg$data$derived, "power", "segment_power")
model$condition_order <- as.character(cfg$condition_order)
if (!identical(model$FORMAL_BANDS, branch_bands) || !identical(model$USE_REML, FALSE) ||
    !identical(model$OPTIMIZER, "bobyqa") || model$MAXFUN != 200000L) {
  stop("Canonical model settings changed; re-audit before this sensitivity run.", call. = FALSE)
}
original_prepared_root <- file.path(cfg$data$derived, "power", "encoding_lmm")
prepared_root <- file.path(original_prepared_root, branch_name)
original_results_root <- file.path(cfg$results$power, "encoding_lmm")
results_root <- file.path(original_results_root, branch_name)
branch_assert_isolated(prepared_root, original_prepared_root)
branch_assert_isolated(results_root, original_results_root)
if (dir.exists(prepared_root) || dir.exists(results_root)) {
  stop("Supplemental output directory already exists. Review partial/completed outputs; no automatic overwrite or cleanup.", call. = FALSE)
}
original_corridor_path <- file.path(cfg$data$encoding, "corridor_subsegi_encoding_table.csv")
original_junction_path <- file.path(cfg$data$encoding, "turn_encoding_table.csv")
reversed_corridor_path <- file.path(cfg$data$encoding, "corridor_subsegi_encoding_table_position_to_start.csv")
reversed_junction_path <- file.path(cfg$data$encoding, "turn_encoding_table_position_to_start.csv")
original_summary_path <- file.path(original_results_root, "summaries", "all_target_predictor_effects.csv")
original_corridor <- prep$read_corridor_encoding(original_corridor_path)
original_junction <- prep$read_junction_encoding(original_junction_path)
reversed_corridor <- prep$read_corridor_encoding(reversed_corridor_path)
reversed_junction <- prep$read_junction_encoding(reversed_junction_path)
branch_check_reversal(original_corridor, reversed_corridor, c("trial", "segment_no", "part_no"),
                      c("Angle2Xaxis", "Angle2orientation"), "Corridor encoding")
branch_check_reversal(original_junction, reversed_junction, c("trial", "segment_no"),
                      "Angle2Xaxis", "Junction encoding")
if (!file.exists(original_summary_path)) stop("Original model summary absent; cannot compare.", call. = FALSE)
original <- read_csv(original_summary_path, show_col_types = FALSE, progress = FALSE) |>
  filter(predictor %in% c("Angle2Xaxis", "Angle2orientation"), roi %in% model$SELECTED_ROIS)
model_keys <- c("roi", "scope", "variant", "band_name", "predictor")
if (anyDuplicated(branch_key(original, model_keys))) stop("Original summary model keys duplicated.", call. = FALSE)
original <- branch_fdr(original)
tracked <- c(prep_file, model_file, original_corridor_path, original_junction_path,
             reversed_corridor_path, reversed_junction_path, original_summary_path,
             file.path(cfg$config_dir, "project_config.json"))
prepared <- list(); baseline_prepared <- list(); preparation_qc <- list()
for (roi in model$SELECTED_ROIS) {
  message("Auditing and preparing ", roi)
  old_c <- prep$prepare_corridor(roi, original_corridor)
  old_j <- prep$prepare_junction(roi, original_junction)
  new_c <- prep$prepare_corridor(roi, reversed_corridor)
  new_j <- prep$prepare_junction(roi, reversed_junction)
  old <- list(corridor = old_c, junction = old_j, all_segments = prep$make_all_segments(old_c, old_j))
  new <- list(corridor = new_c, junction = new_j, all_segments = prep$make_all_segments(new_c, new_j))
  tracked <- c(tracked, unique(c(old_c$source_power_file, old_j$source_power_file)))
  for (scope in names(branch_predictors)) {
    keys <- c("sub_id", "trial", "segment_type", "segment_no", "part_no", "band_name")
    old_path <- file.path(original_prepared_root, scope, paste0(prep$safe_filename(roi), ".csv"))
    if (!file.exists(old_path)) stop("Original prepared dataset absent: ", old_path, call. = FALSE)
    old_saved <- read_csv(old_path, show_col_types = FALSE, progress = FALSE)
    branch_same_data(old_saved, old[[scope]], keys,
      c("roi", "condition", "power_db", "Distance", "Angle2Xaxis", "Angle2orientation",
        "SegmentNumber", "is_last_segment", "band_low_hz", "band_high_hz"),
      paste("Original prepared data audit", roi, scope))
    branch_same_data(new[[scope]], old[[scope]], keys,
      setdiff(names(old[[scope]]), c("Angle2Xaxis", "Angle2orientation")),
      paste("Only bearing changed", roi, scope))
    prepared[[paste(roi, scope, sep = "__")]] <- new[[scope]]
    baseline_prepared[[paste(roi, scope, sep = "__")]] <- old[[scope]]
    tracked <- c(tracked, old_path)
  }
}
tracked <- unique(tracked)
before <- tools::md5sum(tracked)
dir.create(prepared_root, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
targets <- list(); fixed <- list(); qc <- list(); attempted <- 0L
baseline_targets <- list(); baseline_qc <- list()
for (roi in model$SELECTED_ROIS) {
  roi_fixed <- list()
  for (scope in names(branch_predictors)) {
    dat_scope <- prepared[[paste(roi, scope, sep = "__")]]
    baseline_scope <- baseline_prepared[[paste(roi, scope, sep = "__")]]
    out_path <- file.path(prepared_root, scope, paste0(prep$safe_filename(roi), ".csv"))
    branch_write_csv(dat_scope, out_path)
    preparation_qc[[length(preparation_qc) + 1L]] <- prep$summarize_preparation(dat_scope, scope, roi, out_path)
    for (variant in c("include_last", "exclude_last")) {
      dat_variant <- if (variant == "exclude_last") filter(dat_scope, !(is_last_segment %in% TRUE)) else dat_scope
      baseline_variant <- if (variant == "exclude_last") filter(baseline_scope, !(is_last_segment %in% TRUE)) else baseline_scope
      for (band in branch_bands) {
        dat_band <- dat_variant |> filter(band_name == .env$band)
        if (!nrow(dat_band) || !identical(unique(as.character(dat_band$band_name)), band)) {
          stop("Requested band missing or inconsistent: ", roi, " / ", scope, " / ", band, call. = FALSE)
        }
        for (predictor in branch_predictors[[scope]]) {
          baseline_fit <- model$fit_one_model(baseline_variant |> filter(band_name == .env$band),
                                              roi, scope, variant, band, predictor)
          fit <- model$fit_one_model(dat_band, roi, scope, variant, band, predictor)
          baseline_qc[[length(baseline_qc) + 1L]] <- baseline_fit$qc
          if (!is.null(baseline_fit$target)) baseline_targets[[length(baseline_targets) + 1L]] <- baseline_fit$target
          attempted <- attempted + 1L
          qc[[length(qc) + 1L]] <- fit$qc
          if (!is.null(fit$fixed)) roi_fixed[[length(roi_fixed) + 1L]] <- fit$fixed
          if (!is.null(fit$target)) targets[[length(targets) + 1L]] <- fit$target
          message("Model ", attempted, ": ", fit$qc$model_id, " [", fit$qc$fit_status, "]")
        }
      }
    }
  }
  fixed[[roi]] <- bind_rows(roi_fixed)
  branch_write_csv(fixed[[roi]], file.path(results_root, "by_roi", paste0(prep$safe_filename(roi), "_lmm_detailed.csv")))
}
model_qc <- bind_rows(qc)
baseline_model_qc <- bind_rows(baseline_qc)
expected_count <- length(model$SELECTED_ROIS) * length(branch_bands) * 2L * sum(lengths(branch_predictors))
if (attempted != expected_count || anyDuplicated(branch_key(model_qc, model_keys))) stop("Model coverage audit failed.", call. = FALSE)
branch_same_data(model_qc, baseline_model_qc, model_keys,
  c("n_observations", "n_subjects", "n_subject_trials", "n_conditions", "n_segment_types", "power_mean_db", "power_sd_db"),
  "Reproduced original vs reversed model sample audit")
targets <- branch_fdr(bind_rows(targets))
baseline_targets <- branch_fdr(bind_rows(baseline_targets))
comparison_metrics <- c("estimate", "std_error", "statistic", "p_value_raw", "p_fdr", "n_observations",
                        "n_subjects", "n_subject_trials", "singular_fit", "fit_status")
baseline_reproduction <- full_join(
  original |> select(all_of(c(model_keys, comparison_metrics))),
  baseline_targets |> select(all_of(c(model_keys, comparison_metrics))),
  by = model_keys, suffix = c("_stored", "_reproduced")) |>
  mutate(beta_difference = estimate_reproduced - estimate_stored,
    p_difference = p_value_raw_reproduced - p_value_raw_stored,
    sample_counts_identical = n_observations_stored == n_observations_reproduced &
      n_subjects_stored == n_subjects_reproduced & n_subject_trials_stored == n_subject_trials_reproduced,
    numerical_reproduction_matches = !is.na(sample_counts_identical) & sample_counts_identical &
      is.finite(beta_difference) & abs(beta_difference) <= 1e-7 * pmax(1, abs(estimate_stored)) &
      is.finite(p_difference) & abs(p_difference) <= 1e-7)
comparison <- full_join(
  original |> select(all_of(c(model_keys, comparison_metrics))),
  targets |> select(all_of(c(model_keys, comparison_metrics))),
  by = model_keys, suffix = c("_original", "_position_to_start")) |>
  mutate(beta_difference = estimate_position_to_start - estimate_original,
         p_difference = p_value_raw_position_to_start - p_value_raw_original,
         raw_significant_original = is.finite(p_value_raw_original) & p_value_raw_original < 0.05,
         raw_significant_position_to_start = is.finite(p_value_raw_position_to_start) & p_value_raw_position_to_start < 0.05,
         fdr_significant_original = is.finite(p_fdr_original) & p_fdr_original < 0.05,
         fdr_significant_position_to_start = is.finite(p_fdr_position_to_start) & p_fdr_position_to_start < 0.05,
         sample_counts_identical = n_observations_original == n_observations_position_to_start &
           n_subjects_original == n_subjects_position_to_start &
           n_subject_trials_original == n_subject_trials_position_to_start)
reproduced_vs_new <- full_join(
  baseline_targets |> select(all_of(c(model_keys, comparison_metrics))),
  targets |> select(all_of(c(model_keys, comparison_metrics))),
  by = model_keys, suffix = c("_original_reproduced", "_position_to_start"))
if (any(comparison$sample_counts_identical %in% FALSE)) stop("Old/new fitted model sample counts differ; review rather than interpret as coding effect.", call. = FALSE)
branch_write_csv(bind_rows(preparation_qc), file.path(prepared_root, "qc", "preparation_qc.csv"))
branch_write_csv(model_qc, file.path(results_root, "qc", "model_qc.csv"))
branch_write_csv(targets, file.path(results_root, "summaries", "all_target_predictor_effects.csv"))
branch_write_csv(baseline_targets, file.path(results_root, "qc", "original_coding_reproduced_effects.csv"))
branch_write_csv(baseline_model_qc, file.path(results_root, "qc", "original_coding_reproduced_model_qc.csv"))
branch_write_csv(baseline_reproduction, file.path(results_root, "qc", "original_result_reproduction_check.csv"))
branch_write_csv(reproduced_vs_new, file.path(results_root, "comparison_with_reproduced_original.csv"))
branch_write_csv(comparison, file.path(results_root, "comparison_with_original.csv"))
after <- tools::md5sum(tracked)
if (!identical(before, after)) stop("An audited input changed during the run; outputs not marked successful.", call. = FALSE)
branch_write_csv(tibble(path = names(before), md5_before = unname(before), md5_after = unname(after)),
                 file.path(results_root, "qc", "input_hashes.csv"))
branch_write_lines(c(
  "Bearing sensitivity: position-to-start rather than start-to-position.",
  "Angle2Xaxis and Angle2orientation replaced by (original + pi) modulo 2*pi, radians in [0, 2*pi).",
  "No neural power re-extraction; no change to model, covariates, filtering, raw units, or terminal-segment rules.",
  "Exact S01 preparation and S02 model functions imported via selected AST assignments; original main bodies never executed.",
  "FDR separately per predictor x scope x include/exclude branch: aMTG 5 bands; three hippocampal ROIs jointly 15 bands.",
  "pMTG retained among the canonical 5 fitted ROIs; raw p only, no newly invented FDR family and no formal angle heatmap row.",
  "Comparison contains original/new estimates, raw p, q, and sample counts; missing q is not a non-significant tested result.",
  "Original-coding models are independently reproduced inside this branch; QC compares them with stored original results, and a second comparison uses the reproduced baseline.",
  "All original inputs and outputs retained. See QC for numerical warnings and failed/singular fits."
), file.path(results_root, "ANALYSIS_README.txt"))
capture.output(sessionInfo(), file = file.path(results_root, "qc", "sessionInfo.txt"))
branch_write_lines(c(paste("models_attempted", attempted), paste("models_ok", nrow(targets)),
  paste("models_not_ok", sum(model_qc$fit_status != "ok")),
  paste("stored_baseline_reproduction_matches", sum(baseline_reproduction$numerical_reproduction_matches %in% TRUE)),
  paste("stored_baseline_reproduction_missing_or_mismatches", sum(!(baseline_reproduction$numerical_reproduction_matches %in% TRUE)))),
  file.path(results_root, "_SUCCESS"))
message("Supplemental Power completed: ", results_root)
