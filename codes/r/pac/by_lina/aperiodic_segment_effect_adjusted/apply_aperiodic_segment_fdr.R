#!/usr/bin/env Rscript
# Post-process existing segment coefficients only. Never refit or alter raw results.
suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) args[[1]] else
  "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
root <- file.path(project_root, "results/pac/by_lina/aperiodic_segment_effect_adjusted")
input <- file.path(root, "all_segment_effects.csv")
outputs <- file.path(root, c("all_segment_effects_fdr_bh.csv",
                             "FDR_CORRECTION_NOTES.txt",
                             "fdr_analysis_manifest.json", "_FDR_SUCCESS"))
if (any(file.exists(outputs))) {
  stop("FDR output already exists; no existing output will be overwritten.")
}
stopifnot(file.exists(input), file.exists(file.path(root, "_SUCCESS")))
scopes <- c("include_last", "exclude_terminal_corridor_and_junction")
amtg <- "Temporal_Mid L_Anterior"
hipp <- c("Hippocampus L_Anterior", "Hippocampus L_Posterior",
          "Hippocampus R_Anterior")
rois <- c(amtg, hipp)
models <- unlist(lapply(scopes, function(scope) {
  list.files(file.path(root, scope, "models"), pattern = "\\.rds$", full.names = TRUE)
}), use.names = FALSE)
stopifnot(length(models) == 8L)
protected <- c(input, file.path(root, c("analysis_manifest.json", "_SUCCESS")), models)
sha <- function(path) digest::digest(file = path, algo = "sha256")
before <- vapply(protected, sha, character(1))
dat <- readr::read_csv(input, show_col_types = FALSE)
required <- c("scope", "roi", "exponent_range", "term", "p_value", "fdr_status")
stopifnot(all(required %in% names(dat)), nrow(dat) == 8L,
          identical(sort(unique(dat$scope)), sort(scopes)),
          all(dat$term == "segment_position"),
          all(dat$exponent_range == "1-150 Hz"),
          all(dat$fdr_status == "not_applied"),
          all(is.finite(dat$p_value)), all(dat$p_value >= 0 & dat$p_value <= 1),
          !anyDuplicated(paste(dat$scope, dat$roi)))
for (scope in scopes) {
  stopifnot(identical(sort(dat$roi[dat$scope == scope]), sort(rois)))
}

# Preserve all original coefficient and P-value columns. Add explicit policy columns.
dat$fdr_status <- ifelse(dat$roi == amtg, "not_required_single_aMTG_test",
                        "BH_applied_three_hippocampi")
dat$correction_method <- ifelse(dat$roi == amtg, "none", "BH")
dat$correction_family <- paste(dat$scope,
                               ifelse(dat$roi == amtg, "aMTG_single_test", "three_hippocampi"),
                               sep = "__")
dat$correction_n_tests <- ifelse(dat$roi == amtg, 1L, 3L)
dat$p_fdr_bh <- NA_real_
for (scope in scopes) {
  idx <- which(dat$scope == scope & dat$roi %in% hipp)
  stopifnot(length(idx) == 3L)
  dat$p_fdr_bh[idx] <- p.adjust(dat$p_value[idx], method = "BH", n = 3L)
}
dat$p_for_significance <- ifelse(dat$roi == amtg, dat$p_value, dat$p_fdr_bh)
dat$significance_basis <- ifelse(dat$roi == amtg, "raw_p", "BH_adjusted_p")
dat$significant_at_0p05 <- dat$p_for_significance < 0.05
stopifnot(all(is.na(dat$p_fdr_bh[dat$roi == amtg])),
          all(is.finite(dat$p_for_significance)),
          all(dat$p_fdr_bh[dat$roi %in% hipp] >= dat$p_value[dat$roi %in% hipp]),
          all(dat$p_fdr_bh[dat$roi %in% hipp] <= 1))

readr::write_csv(dat, outputs[[1]], na = "")
saved <- readr::read_csv(outputs[[1]], show_col_types = FALSE)
stopifnot(identical(saved$p_value, dat$p_value),
          isTRUE(all.equal(saved$p_for_significance, dat$p_for_significance,
                           tolerance = 1e-14)),
          isTRUE(all.equal(saved$p_fdr_bh, dat$p_fdr_bh, tolerance = 1e-14)),
          identical(saved$significant_at_0p05, dat$significant_at_0p05))
after <- vapply(protected, sha, character(1))
stopifnot(identical(before, after))

notes <- c(
  "Correction policy for 1-150 Hz IRASA exponent segment-position effects:",
  "Model: exponent_raw ~ segment_position + condition + (1|sub_id).",
  "aMTG: one target segment effect in each branch; no multiple-testing correction. Use raw P.",
  "Hippocampus: left anterior, left posterior, right anterior jointly form one BH-FDR family (m=3).",
  "Apply BH independently within include_last and exclude_terminal_corridor_and_junction.",
  "Do not combine the two branches. Do not include intercept or condition coefficients in these families.",
  "Only the selected 1-150 Hz exponent is included. No additional frequency-range tests are included.",
  "The exclude branch removes both halves of the last corridor and the last junction; original order is retained.",
  "p_value is the untouched original P; p_fdr_bh is blank for aMTG (not a corrected value).",
  "p_for_significance uses raw P for aMTG and BH-adjusted P for hippocampus; threshold < 0.05.",
  "The raw model tables, RDS files, fit manifest, and original _SUCCESS remain unchanged.",
  "The original fit manifest's no-FDR entry describes the fitting stage; this manifest records later correction.",
  "This step does not refit models, recompute IRASA, or generate figures."
)
writeLines(notes, outputs[[2]], useBytes = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
manifest <- list(
  analysis = "BH post-processing of existing aperiodic exponent segment effects",
  completed_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
  source_table = input, source_sha256 = before[[1]],
  script_sha256 = if (!is.na(script_path)) sha(script_path) else NA_character_,
  exponent_range = "1-150 Hz", target_term = "segment_position", alpha = 0.05,
  amtg_policy = "raw P; one target effect per branch; no correction",
  hippocampus_method = "stats::p.adjust(method='BH', n=3)",
  hippocampus_rois = hipp, scopes = scopes, families_kept_separate_by_scope = TRUE,
  hippocampus_family_size = 3L, hippocampus_families = 2L,
  no_refitting = TRUE, no_plotting = TRUE,
  protected_files = data.frame(path = protected, sha256_before = before,
                               sha256_after = after, unchanged = before == after),
  results = as.data.frame(dat[, c("scope", "roi", "p_value", "p_fdr_bh",
                                 "p_for_significance", "significant_at_0p05")]),
  R_version = R.version.string
)
jsonlite::write_json(manifest, outputs[[3]], pretty = TRUE, auto_unbox = TRUE,
                     digits = NA, na = "null")
writeLines("8 effects processed; two hippocampus BH families of 3; aMTG raw P; protected files unchanged",
           outputs[[4]], useBytes = TRUE)
print(as.data.frame(dat[, c("scope", "roi", "p_value", "p_fdr_bh",
                           "p_for_significance", "significant_at_0p05")]), row.names = FALSE)
