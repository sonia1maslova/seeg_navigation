#!/usr/bin/env Rscript
# Post-process Lina's existing whole-trial PAC-behavior effects; never refit.
suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
  library(digest)
})
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) args[[1]] else
  "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
source_root <- file.path(project_root,
  "results/pac/by_lina/formal_raw_pac/04_amtg_followup_and_behavior")
source_path <- file.path(source_root, "model_coefficients.csv")
output_root <- file.path(project_root, "results/pac/by_lina/pac_behavior")
outputs <- file.path(output_root, c("pac_behavior_effects_fdr_bh.csv",
  "FDR_CORRECTION_NOTES.txt", "fdr_analysis_manifest.json", "_SUCCESS"))
if (any(file.exists(outputs))) stop("Existing correction output will not be overwritten.")
rois <- c("Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
          "Hippocampus L_Posterior", "Hippocampus R_Anterior")
phase_bands <- c("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
amp_bands <- c("low_gamma_30_70Hz", "high_gamma_70_150Hz")
config_path <- file.path(project_root, "config/project_config.json")
protected <- c(source_path, file.path(source_root, c("fdr_manifest.csv",
                "analysis_dataset_whole_trial_behavior.csv")), config_path)
sha <- function(path) digest::digest(file = path, algo = "sha256")
stopifnot(all(file.exists(protected)))
before <- vapply(protected, sha, character(1))
all <- readr::read_csv(source_path, show_col_types = FALSE)
dat <- as.data.frame(all[all$analysis_family == "04C1_pac_behavior", ])
stopifnot(nrow(dat) == 24L, all(dat$scope == "whole_trial"),
          all(dat$target_term == "raw_PAC"), all(dat$fit_status == "ok"),
          all(dat$is_singular == 0),
          all(is.finite(dat$p_value)), all(dat$p_value >= 0 & dat$p_value <= 1),
          !anyDuplicated(paste(dat$roi, dat$phase_band, dat$amp_band)))
for (roi in rois) {
  cell <- dat[dat$roi == roi, ]
  planned <- expand.grid(phase_band = phase_bands, amp_band = amp_bands,
                         stringsAsFactors = FALSE)
  stopifnot(nrow(cell) == 6L,
    identical(sort(paste(cell$phase_band, cell$amp_band)),
              sort(paste(planned$phase_band, planned$amp_band))))
}
# Explicitly retain the original 24-test correction, under provenance columns.
names(dat)[names(dat) == "q_value"] <- "q_value_lina_24_tests"
names(dat)[names(dat) == "fdr_family_id"] <- "fdr_family_id_lina_24_tests"
names(dat)[names(dat) == "m_planned"] <- "m_planned_lina_24_tests"
names(dat)[names(dat) == "m_valid"] <- "m_valid_lina_24_tests"
dat$fdr_family_id <- ifelse(dat$roi == rois[[1]], "aMTG_6_phase_amplitude_combinations",
                          "three_hippocampi_18_phase_amplitude_combinations")
dat$m_planned <- ifelse(dat$roi == rois[[1]], 6L, 18L)
dat$m_valid <- dat$m_planned
dat$q_value <- NA_real_
for (family in unique(dat$fdr_family_id)) {
  idx <- which(dat$fdr_family_id == family)
  stopifnot(length(idx) == dat$m_planned[idx[[1]]])
  dat$q_value[idx] <- p.adjust(dat$p_value[idx], method = "BH", n = length(idx))
}
dat$correction_method <- "BH"
dat$raw_stars <- ifelse(dat$p_value < 0.001, "***", ifelse(dat$p_value < 0.01, "**",
                        ifelse(dat$p_value < 0.05, "*", "")))
dat$fdr_dagger <- ifelse(dat$q_value < 0.05, "\u2020", "")
dat$plot_symbol <- paste0(dat$raw_stars, dat$fdr_dagger)
dat$raw_significant <- dat$p_value < 0.05
dat$fdr_significant <- dat$q_value < 0.05
cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
roi_lookup <- setNames(vapply(cfg$roi_abbreviations, function(row) row$abbreviation, character(1)),
                       vapply(cfg$roi_abbreviations, function(row) row$full_name, character(1)))
dat$roi_label <- unname(roi_lookup[dat$roi])
stopifnot(!anyNA(dat$roi_label), all(is.finite(dat$q_value)),
          all(dat$q_value >= dat$p_value), all(dat$q_value <= 1))
dat <- dat[order(match(dat$roi, rois), match(dat$amp_band, amp_bands),
                 match(dat$phase_band, phase_bands)), ]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(dat, outputs[[1]], na = "")
saved <- readr::read_csv(outputs[[1]], show_col_types = FALSE)
saved$plot_symbol[is.na(saved$plot_symbol)] <- ""
stopifnot(isTRUE(all.equal(saved$p_value, dat$p_value, tolerance = 1e-14)),
          isTRUE(all.equal(saved$q_value, dat$q_value, tolerance = 1e-14)),
          identical(saved$plot_symbol, dat$plot_symbol))
after <- vapply(protected, sha, character(1))
stopifnot(identical(before, after))
notes <- c(
  "Lina-derived whole-trial PAC-behavior analysis: plotting correction policy",
  "Selected family: 04C1_pac_behavior; pointing_error ~ raw_PAC + (1|sub_id) + (1|junctions).",
  "Use the original estimates, standard errors, 95% confidence intervals, and raw P values; do not refit.",
  "BH-FDR family 1: aMTG, three phase bands x two amplitude bands, m=6.",
  "BH-FDR family 2: left anterior, left posterior, right anterior hippocampus, 3 ROIs x 6 combinations, m=18.",
  "Phase: 0.5-4, 4-12, 12-30 Hz. Amplitude: 30-70, 70-150 Hz.",
  "Include all six combinations in each ROI; do not correct amplitude panels separately.",
  "Do not include 04B segment trends or 04C2 ROI interaction effects in these two families.",
  "Whole-trial PAC has no include_last/exclude_last segment branches; do not create duplicate tests.",
  "Original q_value is retained as q_value_lina_24_tests; q_value contains the new requested family correction.",
  "Plot stars: raw P <.05/.01/.001; dagger: new BH q <.05. Original 95% intervals are not simultaneous intervals.",
  "Lina's original files remain unchanged. The new correction is a separately saved post-processing result."
)
writeLines(notes, outputs[[2]], useBytes = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
manifest <- list(
  analysis = "whole-trial PAC versus pointing error; original model retained",
  completed_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
  source_path = source_path, source_sha256 = before[[1]],
  script_sha256 = if (!is.na(script_path)) sha(script_path) else NA_character_,
  target_family = "04C1_pac_behavior", target_term = "raw_PAC", scope = "whole_trial",
  model_formula = "pointing_error ~ raw_PAC + (1|sub_id) + (1|junctions)",
  observation_unit = "subject x trial x ROI x phase-amplitude combination",
  PAC_definition = "whole-trial raw Tort MI, equal-averaged over within-ROI contacts in Lina source",
  fdr_method = "stats::p.adjust(method='BH')", family_sizes = c(aMTG = 6L, hippocampus = 18L),
  no_refitting = TRUE, raw_results_unchanged = TRUE, alpha = 0.05,
  protected_files = data.frame(path = protected, sha256_before = before,
    sha256_after = after, unchanged = before == after),
  effects = 24L, raw_significant_effects = sum(dat$raw_significant),
  fdr_significant_effects = sum(dat$fdr_significant),
  R_version = R.version.string
)
jsonlite::write_json(manifest, outputs[[3]], auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines("24 effects corrected; aMTG m=6; hippocampus m=18; original files unchanged", outputs[[4]])
print(dat[dat$fdr_significant, c("roi", "phase_band", "amp_band", "p_value", "q_value", "plot_symbol")],
      row.names = FALSE)
