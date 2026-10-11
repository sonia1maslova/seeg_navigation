#!/usr/bin/env Rscript
# Lina-derived aperiodic exponent trends, not exponent-adjusted PAC.
# No FDR, signal extraction, IRASA recomputation, or plotting is performed.
suppressPackageStartupMessages({
  library(lme4)
  library(readr)
  library(dplyr)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) args[[1]] else "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
output_root <- file.path(project_root, "results/pac/by_lina/aperiodic_segment_effect_adjusted")
source_path <- if (length(args) >= 2L) args[[2]] else
  file.path(output_root, "inputs/lina_exponent_trial_level_dataset.csv")
source_original <- paste0(
  "/nfs/visu/workplace/na_li_/3dpointing/",
  "results_raw_pac_phase_set_0p5_4_1_4_4_12_12_30_change_name_roi_right_hpc_ap_full_rois/",
  "22_focused_hpc_amtg_aperiodic_control/00_joined_input/",
  "pipeline3_exponent_trial_level_dataset.csv")
scopes <- c("include_last", "exclude_terminal_corridor_and_junction")
rois <- c("Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
          "Hippocampus L_Posterior", "Hippocampus R_Anterior")
range_primary <- "1-150 Hz"
formula <- exponent_raw ~ segment_position + condition + (1 | sub_id)
formula_text <- "exponent_raw ~ segment_position + condition + (1|sub_id)"
config_path <- file.path(project_root, "config/project_config.json")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
started_utc <- format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
sha <- function(path) digest::digest(file = path, algo = "sha256")
stopifnot(file.exists(source_path), file.exists(source_original), file.exists(config_path))
if (file.exists(file.path(output_root, "_SUCCESS"))) {
  stop("A successful run already exists. Existing results will not be overwritten.")
}
for (scope in scopes) {
  if (dir.exists(file.path(output_root, scope)) &&
      length(list.files(file.path(output_root, scope), all.files = FALSE))) {
    stop("Nonempty result branch will not be overwritten: ", scope)
  }
}
dir.create(file.path(output_root, "audit"), recursive = TRUE, showWarnings = FALSE)
source_sha <- sha(source_path)
stopifnot(identical(source_sha, sha(source_original)))
cfg <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
condition_levels <- unlist(cfg$condition_order, use.names = FALSE)
stopifnot(identical(condition_levels, c("turns_num2", "turns_num4", "turns_num6")))
condition_mapping <- bind_rows(lapply(condition_levels, function(label) {
  data.frame(condition = label,
             trial = as.integer(unlist(cfg$conditions[[label]]$trial_ids)),
             turns = as.integer(cfg$conditions[[label]]$n_junctions))
}))
stopifnot(!anyDuplicated(condition_mapping$trial))
roi_labels <- setNames(
  vapply(cfg$roi_abbreviations, function(item) item$abbreviation, character(1)),
  vapply(cfg$roi_abbreviations, function(item) item$full_name, character(1)))
stopifnot(all(rois %in% names(roi_labels)))

# Protect both the Lina results and the previous PAC-as-outcome re-fits.
protected_paths <- c(source_original, source_path, config_path,
  file.path(project_root, "results/pac/by_lina/formal_raw_pac/03b_hpc_exponent_control",
            c("model_coefficients.csv", "analysis_dataset_exponent_tracking.csv")),
  file.path(project_root, "results/pac/by_lina/segment_effect_adjusted", scopes,
            "all_segment_effects.csv"))
protected_paths <- unique(protected_paths[file.exists(protected_paths)])
protected_before <- data.frame(path = protected_paths,
                               sha256_before = vapply(protected_paths, sha, character(1)))

source <- read_csv(source_path, show_col_types = FALSE, progress = FALSE,
                   col_types = cols(sub_id = col_character(), .default = col_guess()))
required <- c("sub_id", "trial", "location1", "turns", "segment_id", "segment_type",
              "segment_order", "exponent_range", "exponent_raw", "is_terminal")
if (!all(required %in% names(source))) stop("Missing source columns")
selected <- source |>
  filter(location1 %in% rois, exponent_range == range_primary) |>
  mutate(roi = location1, trial = as.integer(trial), turns = as.integer(turns),
         segment_order = as.integer(segment_order))
if (!setequal(unique(selected$roi), rois)) stop("Missing target ROI")
identity <- c("sub_id", "trial", "roi", "segment_id", "exponent_range")
if (anyDuplicated(as.data.frame(selected[identity]))) {
  stop("Duplicate independent exponent observations; PAC-combination copies are not allowed")
}
if (anyNA(selected[c("sub_id", "trial", "roi", "turns", "segment_id", "segment_type",
                     "segment_order")]) || any(!selected$turns %in% c(2L, 4L, 6L))) {
  stop("Missing or unknown observation identity")
}
trial_conditions <- selected |> distinct(sub_id, trial, turns) |> count(sub_id, trial)
if (any(trial_conditions$n != 1L)) stop("A trial has inconsistent conditions")
selected <- selected |>
  left_join(rename(condition_mapping, config_turns = turns), by = "trial")
if (anyNA(selected$condition) || any(selected$turns != selected$config_turns)) {
  stop("Input turns/trial do not match project condition definitions")
}

# Derive path positions from the named units, not the maximum observed row.
is_c <- grepl("^C[0-9]+_(first|second)$", selected$segment_id)
is_j <- grepl("^J[0-9]+$", selected$segment_id)
if (any(!(is_c | is_j))) stop("Unknown segment ID")
segment_index <- integer(nrow(selected))
segment_index[is_c] <- as.integer(sub("^C([0-9]+)_.*$", "\\1", selected$segment_id[is_c]))
segment_index[is_j] <- as.integer(sub("^J([0-9]+)$", "\\1", selected$segment_id[is_j]))
expected_order <- integer(nrow(selected))
expected_order[is_c] <- 3L * (segment_index[is_c] - 1L) + 1L +
  as.integer(grepl("_second$", selected$segment_id[is_c]))
expected_order[is_j] <- 3L * segment_index[is_j]
if (any(expected_order != selected$segment_order) ||
    any(is_c & selected$segment_type != "corridor") ||
    any(is_j & selected$segment_type != "junction") ||
    any(segment_index < 1L) || any(is_c & segment_index > selected$turns + 1L) ||
    any(is_j & segment_index > selected$turns)) stop("Navigation identity/order inconsistency")
selected <- selected |>
  mutate(segment_position = segment_order,
         is_terminal_corridor = is_c & segment_index == turns + 1L,
         is_terminal_junction = is_j & segment_index == turns,
         keep_after_terminal_exclusion = !(is_terminal_corridor | is_terminal_junction),
         valid_outcome = is.finite(exponent_raw))
if (anyNA(selected$is_terminal) ||
    any(as.logical(selected$is_terminal) != selected$is_terminal_corridor)) {
  stop("Original terminal flag disagrees with corridor-only flag")
}
finite <- selected |> filter(valid_outcome) |> arrange(roi, sub_id, trial, segment_order)
selection_qc <- selected |> group_by(roi, condition) |> summarise(
  rows_selected = n(), rows_nonfinite_outcome = sum(!valid_outcome),
  include_rows = sum(valid_outcome),
  terminal_corridor_rows_removed = sum(valid_outcome & is_terminal_corridor),
  terminal_junction_rows_removed = sum(valid_outcome & is_terminal_junction),
  adjusted_exclude_rows = sum(valid_outcome & keep_after_terminal_exclusion),
  n_subjects = n_distinct(sub_id), .groups = "drop")
write_csv(selection_qc, file.path(output_root, "audit/input_selection.csv"))
write_csv(selected |> select(all_of(identity), condition, turns, segment_position,
                            is_terminal_corridor, is_terminal_junction,
                            keep_after_terminal_exclusion, valid_outcome),
          file.path(output_root, "audit/row_selection.csv"))

fit_one <- function(dat, roi, scope, model_dir) {
  dat$sub_id <- factor(dat$sub_id)
  dat$condition <- factor(dat$condition, levels = condition_levels)
  contrasts(dat$condition) <- contr.treatment(condition_levels, base = 1)
  if (nlevels(droplevels(dat$condition)) != 3L || nlevels(dat$sub_id) < 2L) {
    stop("Model lacks required conditions/subjects")
  }
  X_expected <- model.matrix(~ segment_position + condition, data = dat)
  if (qr(X_expected)$rank != ncol(X_expected)) stop("Fixed effects not identifiable")
  messages <- character()
  fit <- withCallingHandlers(lmer(formula, data = dat, REML = TRUE,
    control = lmerControl(optimizer = "bobyqa",
                          optCtrl = list(maxfun = 200000, rhoend = 1e-9))),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  optimizer_code <- fit@optinfo$conv$opt
  convergence <- fit@optinfo$conv$lme4$messages
  if (any(optimizer_code != 0) || length(convergence)) {
    stop("Model convergence check failed: ", roi, " / ", scope, " / ",
         paste(c(optimizer_code, convergence), collapse = "; "))
  }
  X <- getME(fit, "X")
  if (ncol(X) != ncol(X_expected) || !is.null(attr(X, "col.dropped"))) stop("Fixed effects dropped")
  beta <- fixef(fit)
  se <- sqrt(diag(as.matrix(vcov(fit))))
  df <- nrow(dat) - qr(X)$rank
  stat <- beta / se
  p <- 2 * pt(-abs(stat), df = df)
  radius <- qt(0.975, df = df) * se
  if (any(!is.finite(c(beta, se, p, radius))) || any(se <= 0) || any(p < 0 | p > 1)) {
    stop("Invalid coefficient inference")
  }
  vc <- as.data.frame(VarCorr(fit))
  diagnostic <- data.frame(scope = scope, roi = roi,
    roi_label = unname(roi_labels[[roi]]), exponent_range = range_primary,
    model_formula = formula_text, n_obs = nrow(dat), n_subjects = nlevels(dat$sub_id),
    n_trials = nrow(unique(dat[c("sub_id", "trial")])),
    n_segments = length(unique(dat$segment_position)), n_conditions = 3L,
    fixed_effect_rank = qr(X)$rank, residual_df = df,
    condition_reference = condition_levels[[1]],
    subject_intercept_variance = vc$vcov[vc$grp == "sub_id"],
    residual_variance = sigma(fit)^2, is_singular = isSingular(fit, tol = 1e-4),
    optimizer_code = paste(optimizer_code, collapse = ";"),
    warnings = paste(unique(messages), collapse = "; "),
    log_likelihood = as.numeric(logLik(fit)), aic = AIC(fit), bic = BIC(fit),
    fit_status = "ok", estimator = "REML", fdr_status = "not_applied")
  fixed <- data.frame(term = names(beta), beta = unname(beta), se = unname(se),
    t_value = unname(stat), df = df, ci_low = unname(beta - radius),
    ci_high = unname(beta + radius), p_value = unname(p))
  fixed <- cbind(diagnostic[rep(1L, nrow(fixed)),
    c("scope", "roi", "roi_label", "exponent_range", "model_formula", "n_obs",
      "n_subjects", "condition_reference", "fdr_status")], fixed)
  slug <- gsub("[^A-Za-z0-9]+", "_", roi)
  saveRDS(fit, file.path(model_dir, paste0(slug, ".rds")))
  list(fixed = fixed, diagnostic = diagnostic)
}

all_target <- list(); all_fixed <- list(); all_diag <- list()
for (scope in scopes) {
  branch_dir <- file.path(output_root, scope)
  dir.create(file.path(branch_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  dat <- if (scope == "include_last") finite else
    finite |> filter(keep_after_terminal_exclusion)
  dat <- dat |> mutate(scope = scope)
  if (scope != "include_last" && (any(dat$is_terminal_corridor | dat$is_terminal_junction) ||
      any(dat$segment_position > 3L * dat$turns - 1L))) stop("Terminal segments survived")
  write_csv(dat, file.path(branch_dir, "analysis_dataset.csv"))
  coverage <- dat |> group_by(roi, condition) |> summarise(
    n_obs = n(), n_subjects = n_distinct(sub_id), n_trials = n_distinct(paste(sub_id, trial)),
    min_position = min(segment_position), max_position = max(segment_position),
    n_corridor = sum(segment_type == "corridor"), n_junction = sum(segment_type == "junction"),
    .groups = "drop")
  write_csv(coverage, file.path(branch_dir, "coverage_by_condition.csv"))
  fits <- lapply(rois, function(roi) {
    cat("Fitting ", scope, " / ", roi, "\n", sep = "")
    fit_one(dat |> filter(.data$roi == .env$roi), roi, scope, file.path(branch_dir, "models"))
  })
  fixed <- bind_rows(lapply(fits, `[[`, "fixed"))
  diagnostic <- bind_rows(lapply(fits, `[[`, "diagnostic"))
  target <- fixed |> filter(term == "segment_position")
  stopifnot(nrow(target) == 4L, nrow(fixed) == 16L)
  write_csv(target, file.path(branch_dir, "segment_effects.csv"))
  write_csv(fixed, file.path(branch_dir, "all_fixed_effects.csv"))
  write_csv(diagnostic, file.path(branch_dir, "model_diagnostics.csv"))
  all_target[[scope]] <- target; all_fixed[[scope]] <- fixed; all_diag[[scope]] <- diagnostic
}
target <- bind_rows(all_target); fixed <- bind_rows(all_fixed); diagnostics <- bind_rows(all_diag)
stopifnot(nrow(target) == 8L, nrow(fixed) == 32L, nrow(diagnostics) == 8L)
write_csv(target, file.path(output_root, "all_segment_effects.csv"))
write_csv(fixed, file.path(output_root, "all_fixed_effects.csv"))
write_csv(diagnostics, file.path(output_root, "model_diagnostics.csv"))
protected_after <- protected_before |>
  mutate(sha256_after = vapply(path, sha, character(1)), unchanged = sha256_before == sha256_after)
write_csv(protected_after, file.path(output_root, "audit/protected_files_sha256.csv"))
stopifnot(all(protected_after$unchanged))
manifest <- list(analysis = "aperiodic exponent segment trend adjusted for categorical condition",
  status = "completed", started_utc = started_utc,
  completed_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
  input_snapshot = source_path, source_original = source_original, source_sha256 = source_sha,
  script_sha256 = if (!is.na(script_path)) sha(script_path) else NULL,
  config_sha256 = sha(config_path), source_rows = nrow(source), selected_rows = nrow(selected),
  rows_with_finite_exponent = nrow(finite), rois = rois, exponent_range = range_primary,
  exponent_definition = "negative log10(power) versus log10(frequency) slope of the IRASA fractal spectrum",
  observation_unit = "subject x ROI x trial x 2.5-s segment; contact exponents already equal-averaged in source",
  outcome = "raw exponent; no normalization or rescaling", model_formula = formula_text,
  condition_levels = condition_levels, condition_reference = condition_levels[[1]],
  condition_encoding = "categorical treatment contrasts; verified against trial_ids and n_junctions in project_config.json",
  scopes = scopes, segment_position = "original navigation segment_order; corridor halves and junction interleaved; never renumber after exclusions",
  exclusion_rule = "per subject-trial, remove C(turns+1)_first, C(turns+1)_second, and J(turns)",
  estimator = "Gaussian subject-random-intercept linear mixed model; REML; lme4::lmer",
  optimizer = "bobyqa; maxfun=200000; rhoend=1e-9",
  inference = "two-sided coefficient Wald t test with residual df = n_obs - rank(X); t-based 95% confidence intervals",
  fdr = "not applied by user request; no q-values generated", figures = "not generated; deferred",
  no_signal_or_irasa_recomputation = TRUE, no_interaction_or_segment_type_covariate = TRUE,
  models_planned = 8L, models_completed = nrow(diagnostics),
  singular_models = sum(diagnostics$is_singular), models_with_warnings = sum(nzchar(diagnostics$warnings)),
  raw_p_lt_0p05_segment_effects = sum(target$p_value < 0.05),
  R_version = R.version.string,
  packages = setNames(lapply(c("lme4", "readr", "dplyr", "jsonlite", "digest"),
                             function(pkg) as.character(packageVersion(pkg))),
                      c("lme4", "readr", "dplyr", "jsonlite", "digest")))
write_json(manifest, file.path(output_root, "analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(c("All P values are uncorrected. No FDR/q values have been calculated.",
             "Multiple-comparison grouping will be decided separately at plotting/reporting stage.",
             "Positive segment beta means increasing exponent, hence a more negative spectral slope.",
             "See code README and analysis_manifest.json for model and input provenance."),
           file.path(output_root, "ANALYSIS_NOTES.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
writeLines("8 models completed; raw P values only; all source hashes unchanged", file.path(output_root, "_SUCCESS"))
print(target[, c("scope", "roi", "beta", "p_value", "n_obs", "n_subjects")], row.names = FALSE)
