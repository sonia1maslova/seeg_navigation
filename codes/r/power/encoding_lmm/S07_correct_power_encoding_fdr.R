# BH-FDR for Power encoding LMMs; never refits or replaces raw model results.
# Each predictor x scope x variant is separate. Hippocampal subregions form
# one planned 3 x 5 family; every other ROI forms its own planned 5-band family.

required_packages <- c("jsonlite", "dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]), mustWork = TRUE)
source(file.path(dirname(dirname(dirname(script_file))), "project_config.R"))

root <- file.path(cfg$results$power, "encoding_lmm")
original <- file.path(root, "summaries", "all_target_predictor_effects.csv")
supplementary <- file.path(root, "supplementary_cortical_rois", "summaries", "all_target_predictor_effects.csv")
inputs <- c(original, supplementary)
if (!all(file.exists(inputs))) stop("Required model results missing: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
source_hashes <- tools::md5sum(inputs)

frames <- lapply(inputs, function(path) {
  dat <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  dat$source_result_file <- path
  dat
})
dat <- dplyr::bind_rows(frames)
bands <- c("f0p5_4", "f4_12", "f12_30", "f30_70", "f70_150")
hipp <- c("Hippocampus L_Anterior", "Hippocampus L_Posterior", "Hippocampus R_Anterior")
keys <- c("roi", "scope", "variant", "predictor", "band_name")
required <- c(keys, "estimate", "conf_low_95", "conf_high_95", "p_value_raw", "fit_status")
if (length(setdiff(required, names(dat)))) stop("Missing columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
if (any(duplicated(dat[keys]))) stop("Duplicate hypothesis keys")
if (!all(dat$band_name %in% bands)) stop("Unexpected frequency band")
if (!all(dat$fit_status == "ok")) stop("Some model fits failed; inspect model QC before correcting")
if (!all(is.finite(dat$p_value_raw) & dat$p_value_raw >= 0 & dat$p_value_raw <= 1)) stop("Invalid raw p values")

dat <- dat |>
  dplyr::mutate(
    fdr_family = ifelse(.data$roi %in% hipp, "hippocampal_subregions", paste0("roi:", .data$roi)),
    fdr_n_planned = ifelse(.data$roi %in% hipp, 15L, 5L)
  ) |>
  dplyr::group_by(.data$predictor, .data$scope, .data$variant, .data$fdr_family) |>
  dplyr::group_modify(function(.x, .y) {
    n_expected <- unique(.x$fdr_n_planned)
    expected_rois <- if (.y$fdr_family == "hippocampal_subregions") hipp else unique(.x$roi)
    grid <- expand.grid(roi = expected_rois, band_name = bands, stringsAsFactors = FALSE)
    actual_keys <- paste(.x$roi, .x$band_name, sep = "|")
    expected_keys <- paste(grid$roi, grid$band_name, sep = "|")
    if (length(n_expected) != 1L || nrow(.x) != n_expected || !setequal(actual_keys, expected_keys)) {
      stop("Incomplete FDR family: ", paste(unlist(.y), collapse = " / "))
    }
    .x$p_fdr <- stats::p.adjust(.x$p_value_raw, method = "BH", n = n_expected)
    .x$fdr_n_observed <- nrow(.x)
    .x
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    raw_sig = .data$p_value_raw < 0.05,
    fdr_sig = .data$p_fdr < 0.05,
    raw_stars = dplyr::case_when(
      .data$p_value_raw < 0.001 ~ "***",
      .data$p_value_raw < 0.01 ~ "**",
      .data$p_value_raw < 0.05 ~ "*",
      TRUE ~ ""
    ),
    fdr_dagger = ifelse(.data$fdr_sig, "\u2020", ""),
    plot_symbol = paste0(.data$raw_stars, .data$fdr_dagger)
  ) |>
  dplyr::arrange(.data$predictor, .data$scope, .data$variant, .data$roi, match(.data$band_name, bands))

family_audit <- dat |>
  dplyr::group_by(.data$predictor, .data$scope, .data$variant, .data$fdr_family) |>
  dplyr::summarise(
    n_planned = dplyr::first(.data$fdr_n_planned), n_observed = dplyr::n(),
    rois = paste(sort(unique(.data$roi)), collapse = ";"),
    raw_significant = sum(.data$raw_sig), fdr_significant = sum(.data$fdr_sig),
    .groups = "drop"
  )
out <- file.path(root, "summaries", "all_target_predictor_effects_fdr_bh.csv")
readr::write_csv(dat, out, na = "")
readr::write_csv(family_audit, file.path(root, "summaries", "power_encoding_fdr_families.csv"), na = "")
notes <- c(
  "Power encoding LMM BH-FDR rule, confirmed 2026-10-10",
  "Separate correction within each predictor x scope x include_last/exclude_last branch.",
  "Hippocampus L_Anterior, L_Posterior, R_Anterior: pooled 3 regions x 5 bands = 15 tests.",
  "Every non-hippocampal ROI: its own 5-band family; no pooling across cortical regions.",
  "All 5 planned bands must be present; family membership is independent of plot layout.",
  "Raw p values, model coefficients, and confidence intervals are unchanged.",
  "Stars denote raw p: * <0.05, ** <0.01, *** <0.001. Dagger denotes BH q <0.05.",
  "Bearing encoding: original start-to-position; not the reversed position-to-start branch.",
  "Whole STG is Temporal_Sup L; no anterior/posterior STG pooling of fitted results.",
  "New compact forest plots use this corrected table; retained Segment-number heatmaps are updated to the same families.",
  "Raw model source files:", inputs
)
writeLines(notes, file.path(root, "summaries", "FDR_CORRECTION_RULES.txt"), useBytes = TRUE)
if (!identical(unname(tools::md5sum(inputs)), unname(source_hashes))) stop("Raw source changed during correction")
jsonlite::write_json(list(
  confirmed_date = "2026-10-10", output = out, n_rows = nrow(dat), n_families = nrow(family_audit),
  input_md5 = as.list(source_hashes), source_files_unchanged = TRUE,
  raw_models_refit = FALSE, correction_method = "BH", hippocampal_tests = 15L,
  other_roi_tests = 5L
), file.path(root, "summaries", "power_encoding_fdr_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Saved corrected rows:", nrow(dat), "families:", nrow(family_audit), "\n", out, "\n")
