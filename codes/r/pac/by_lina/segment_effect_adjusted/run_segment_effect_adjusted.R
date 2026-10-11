# Lina-derived PAC segment effects: exclude terminal C/J or include all C/J.
# Source files stay read-only. REML and residual-DF inference are reproduced
# against MATLAB; current BH families contain six combinations within each ROI.
required <- c("lme4", "dplyr", "readr", "jsonlite")
for (p in required) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else
  "/nfs/visu/workplace/yuhdua/projects/seeg_navigation", mustWork = TRUE)
mode <- if (length(args) >= 2L) args[[2]] else "full"
stopifnot(mode %in% c("pilot", "full", "recalculate_fdr"))
source_root <- file.path(root, "results/pac/by_lina/formal_raw_pac")
scope_new <- if (length(args) >= 3L) args[[3]] else "exclude_terminal_corridor_and_junction"
stopifnot(scope_new %in% c("exclude_terminal_corridor_and_junction", "include_last"))
source_scope <- if (scope_new == "include_last") "all_dynamic_split_with_terminal_C" else "all_dynamic_split"
output_root <- file.path(root, "results/pac/by_lina/segment_effect_adjusted", scope_new)
audit_root <- file.path(output_root, "audit")
if (mode != "recalculate_fdr" && file.exists(file.path(output_root, "_SUCCESS"))) {
  stop("Completed output already exists; refusing to overwrite")
}
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)
stages <- c(hpc_main = "02_hpc_path_accumulation",
            amtg_main = "04_amtg_followup_and_behavior",
            power_control = "03a_hpc_power_control",
            aperiodic_control = "03b_hpc_exponent_control")
datasets <- c(hpc_main = "analysis_dataset.csv",
              amtg_main = "analysis_dataset_amtg_accumulation.csv",
              power_control = "analysis_dataset_pac_power.csv",
              aperiodic_control = "analysis_dataset_pac_exponent.csv")
families <- c(hpc_main = "02A_basic_segment", amtg_main = "04B_amtg_accumulation",
              power_control = "03A2_power_adjusted_PAC",
              aperiodic_control = "03B2_exponent_adjusted_PAC")
planned <- c(hpc_main = 18L, amtg_main = 6L, power_control = 18L, aperiodic_control = 36L)
# The counts above describe models and the ORIGINAL MATLAB audit families only.
# Current reporting groups every ROI (and each exponent range) separately: m=6.
phase_levels <- c("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
amp_levels <- c("low_gamma_30_70Hz", "high_gamma_70_150Hz")
covariates <- list(hpc_main = character(), amtg_main = character(),
                  power_control = c("LF_power_db", "HF_power_db"),
                  aperiodic_control = "exponent_raw")
input_paths <- c(setNames(file.path(source_root, stages, datasets), paste0(names(stages), "_data")),
                 setNames(file.path(source_root, stages, "model_coefficients.csv"),
                          paste0(names(stages), "_coefficients")))
input_md5 <- tools::md5sum(input_paths)
write_csv <- function(x, path) readr::write_csv(x, path, na = "")
read_csv <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
  col_types = readr::cols(sub_id = readr::col_character()))
key_cols <- c("sub_id", "trial", "roi", "phase_band", "amp_band", "segment_id")
row_key <- function(x) do.call(paste, c(x[key_cols], sep = "|"))
fdr_policy <- list(method = "BH", family_size = 6L,
  grouping = c("scope", "analysis_branch", "roi", "exponent_range_if_applicable"),
  phase_bands_per_family = 3L, amplitude_bands_per_family = 2L,
  combine_rois = FALSE, combine_exponent_ranges = FALSE, combine_terminal_scopes = FALSE,
  target_term = "segment_position", expected_families_per_scope = 13L)
apply_roi6_fdr <- function(x) {
  x$exponent_range[is.na(x$exponent_range)] <- ""
  if (nrow(x) != 78L || any(x$scope != scope_new) || any(!is.finite(x$p_value)) ||
      any(x$p_value < 0 | x$p_value > 1)) stop("Invalid model grid for ROI-wise FDR")
  coverage <- x |> count(branch)
  if (!setequal(coverage$branch, names(planned)) ||
      any(coverage$n != unname(planned[coverage$branch]))) stop("Unexpected branch coverage")
  groups <- x |> group_by(branch, scope, roi, exponent_range) |>
    summarise(rows = n(), combinations = n_distinct(phase_band, amp_band), .groups = "drop")
  if (nrow(groups) != 13L || any(groups$rows != 6L | groups$combinations != 6L) ||
      any(!x$phase_band %in% phase_levels) || any(!x$amp_band %in% amp_levels)) {
    stop("Each FDR family must contain exactly all 3 phase x 2 amplitude combinations")
  }
  x |> group_by(branch, scope, roi, exponent_range) |> mutate(
    fdr_family_id = paste0("lina_roi6__", scope, "__", branch, "__",
      gsub("[^A-Za-z0-9]+", "_", roi), "__",
      if_else(exponent_range == "", "no_exponent_range", gsub("[^A-Za-z0-9]+", "_", exponent_range)),
      "__segment_position"),
    m_planned = 6L, m_valid = n(), q_value = p.adjust(p_value, "BH", n = 6L),
    raw_significant = p_value < 0.05, fdr_significant = q_value < 0.05) |> ungroup()
}
persist_statistics <- function(target) {
  summary <- target |> group_by(branch, exponent_range) |>
    summarise(models = n(), raw_significant = sum(raw_significant),
      fdr_significant = sum(fdr_significant), singular_models = sum(is_singular), .groups = "drop")
  family_summary <- target |> group_by(branch, scope, roi, exponent_range, fdr_family_id) |>
    summarise(models = n(), raw_significant = sum(raw_significant),
      fdr_significant = sum(fdr_significant), m_planned = first(m_planned), m_valid = first(m_valid),
      .groups = "drop")
  write_csv(summary, file.path(output_root, "analysis_summary.csv"))
  write_csv(family_summary, file.path(output_root, "fdr_family_summary.csv"))
  write_csv(target, file.path(output_root, "all_segment_effects.csv"))
  write_csv(target |> distinct(branch, scope, roi, exponent_range, fdr_family_id, m_planned, m_valid),
            file.path(output_root, "fdr_manifest.csv"))
  for (analysis in c("main_pac", "power_control", "aperiodic_control")) {
    branches <- if (analysis == "main_pac") c("hpc_main", "amtg_main") else analysis
    out <- file.path(output_root, analysis)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    write_csv(target |> filter(branch %in% branches), file.path(out, "model_coefficients.csv"))
  }
  print(summary)
}
make_comparison <- function(target, baseline) {
  target |> left_join(baseline |> select(branch, roi, phase_band, amp_band, exponent_range,
    original_beta, original_se, original_p_value, original_q_value, original_n_obs),
    by = c("branch", "roi", "phase_band", "amp_band", "exponent_range")) |>
    mutate(beta_change = beta - original_beta, p_change = p_value - original_p_value,
      q_change = q_value - original_q_value, removed_rows = original_n_obs - n_obs,
      q_change_includes_fdr_family_change = TRUE)
}

if (mode == "recalculate_fdr") {
  path <- file.path(output_root, "all_segment_effects.csv")
  if (!file.exists(file.path(output_root, "_SUCCESS")) || !file.exists(path)) {
    stop("Need a completed R result to recalculate FDR without refitting")
  }
  before <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (any(before$fit_status != "ok")) stop("Failed model cannot enter FDR update")
  before$exponent_range[is.na(before$exponent_range)] <- ""
  if (all(before$m_planned == 6L)) {
    updated <- apply_roi6_fdr(before)
    if (max(abs(updated$q_value - before$q_value)) > 1e-12) stop("Existing ROI6 q-values failed check")
    message("ROI-wise six-combination FDR already applied; no writes needed")
    quit(save = "no", status = 0)
  }
  previous_path <- file.path(audit_root, "segment_effects_before_roi6_fdr.csv")
  if (file.exists(previous_path)) stop("Prior FDR audit already exists; refusing to overwrite it")
  target <- before
  for (nm in c("q_value", "fdr_family_id", "m_planned", "m_valid", "fdr_significant")) {
    target[[paste0(nm, "_before_roi6")]] <- before[[nm]]
  }
  target <- apply_roi6_fdr(target)
  stable <- c("beta", "se", "t_value", "df", "p_value", "ci_low", "ci_high", "n_obs", "n_subjects", "n_segments")
  if (!identical(before[stable], as.data.frame(target)[stable])) {
    # readr tibbles and dplyr tibbles may differ in attributes; compare values.
    for (nm in stable) if (!identical(before[[nm]], target[[nm]])) stop("Raw statistic changed: ", nm)
  }
  baseline <- readr::read_csv(file.path(audit_root, "baseline_matlab_vs_r_full.csv"),
                             show_col_types = FALSE, progress = FALSE)
  baseline$exponent_range[is.na(baseline$exponent_range)] <- ""
  manifest <- jsonlite::read_json(file.path(output_root, "analysis_manifest.json"), simplifyVector = FALSE)
  manifest$previous_fdr_policy <- manifest$fdr
  manifest$fdr <- fdr_policy
  manifest$fdr_update_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$fdr_update_refit_models <- FALSE
  manifest$raw_statistics_unchanged <- TRUE
  write_csv(before, previous_path)
  file.copy(file.path(output_root, "analysis_manifest.json"),
            file.path(audit_root, "analysis_manifest_before_roi6_fdr.json"), overwrite = FALSE)
  success <- file.path(output_root, "_SUCCESS")
  if (!file.remove(success)) stop("Cannot clear prior success marker before authorized update")
  persist_statistics(target)
  write_csv(make_comparison(target, baseline), file.path(audit_root, "old_vs_adjusted_segment_effects.csv"))
  write_csv(target |> select(branch, scope, roi, phase_band, amp_band, exponent_range,
    p_value, q_value_before_roi6, q_value, m_planned_before_roi6, m_planned, fdr_significant),
    file.path(audit_root, "fdr_change_to_six_combinations.csv"))
  if (!identical(input_md5, tools::md5sum(input_paths))) stop("Original source files changed during FDR update")
  jsonlite::write_json(manifest, file.path(output_root, "analysis_manifest.json"),
                        pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines("78 R models unchanged; all 13 ROI/range-specific BH families passed, six tests each.", success)
  message("Completed ROI-wise six-test FDR update: ", output_root)
  quit(save = "no", status = 0)
}
prepared <- list()
original <- list()
selection_audit <- list()

for (branch in names(stages)) {
  path <- input_paths[[paste0(branch, "_data")]]
  dat <- read_csv(path)
  if (branch == "aperiodic_control") key_cols <- c(key_cols, "exponent_range")
  cols <- unique(c(key_cols, "turns", "segment_type", "segment_order", "segment_position",
                   "is_terminal", "scope", "raw_PAC", covariates[[branch]]))
  if (length(setdiff(cols, names(dat)))) stop("Missing columns: ", path)
  dat <- dat |> filter(phase_band %in% phase_levels, amp_band %in% amp_levels)
  if (anyNA(dat[cols]) || any(!dat$turns %in% c(2, 4, 6))) stop("Missing identity/unknown trial length")
  trial_lengths <- dat |> distinct(sub_id, trial, turns) |> count(sub_id, trial)
  if (any(trial_lengths$n != 1L)) stop("Inconsistent path length within subject-trial")
  if (any(dat$segment_position != dat$segment_order)) stop("Original position differs from navigation order")
  old <- dat |> filter(scope == "all_dynamic_split")
  included <- dat |> filter(scope == "all_dynamic_split_with_terminal_C")
  rows_old_exclude_C <- nrow(old)
  if (!nrow(old) || !nrow(included)) stop("Missing original terminal variant")
  if (anyDuplicated(row_key(old)) || anyDuplicated(row_key(included))) stop("Duplicate source observations")
  # Compare per row, not membership across trials with different path lengths.
  terminal_C <- included$segment_id == paste0("C", included$turns + 1L, "_first") |
                included$segment_id == paste0("C", included$turns + 1L, "_second")
  terminal_J <- included$segment_id == paste0("J", included$turns)
  keep <- !(terminal_C | terminal_J)
  excluded <- included[keep, ]
  old_kept <- old[old$segment_id != paste0("J", old$turns), ]
  if (!setequal(row_key(excluded), row_key(old_kept))) stop("New selection differs between source variants")
  if (any(excluded$segment_position > 3L * excluded$turns - 1L)) stop("Terminal rows survived selection")
  new <- if (scope_new == "include_last") included else excluded
  if (scope_new == "include_last") old <- included
  numeric_cols <- c("raw_PAC", "segment_position", covariates[[branch]])
  if (any(!is.finite(as.matrix(new[numeric_cols]))) ||
      any(!is.finite(as.matrix(old[numeric_cols])))) stop("Non-finite model input")
  new$source_scope <- new$scope
  new$scope <- scope_new
  new$source_dataset <- substring(path, nchar(root) + 2L)
  minimal <- unique(c(key_cols, "turns", "segment_type", "segment_order", "segment_position",
                       "raw_PAC", covariates[[branch]], "source_scope", "scope", "source_dataset"))
  prepared[[branch]] <- list(old = old, new = new[minimal])
  coeff <- readr::read_csv(input_paths[[paste0(branch, "_coefficients")]],
                           show_col_types = FALSE, progress = FALSE) |>
    filter(analysis_family == unname(families[[branch]]), scope == .env$source_scope,
           predictor == "segment_position", target_term == "segment_position",
           sample_mode == if (length(covariates[[branch]])) "complete_case_adjusted" else "single_roi")
  if (nrow(coeff) != planned[[branch]] || any(coeff$fit_status != "ok")) stop("Incomplete original model grid")
  if (any(coeff$m_planned != planned[[branch]]) || any(coeff$m_valid != planned[[branch]]) ||
      max(abs(coeff$q_value - p.adjust(coeff$p_value, "BH", n = planned[[branch]]))) > 1e-12) {
    stop("Original FDR family mismatch: ", branch)
  }
  original[[branch]] <- coeff
  selection_audit[[branch]] <- data.frame(
    branch = branch, scope = scope_new, rows_old_exclude_C = rows_old_exclude_C,
    rows_include_all = nrow(included), rows_new = nrow(new),
    terminal_corridor_rows_removed = if (scope_new == "include_last") 0L else sum(terminal_C),
    terminal_junction_rows_removed = if (scope_new == "include_last") 0L else sum(terminal_J),
    include_scope_preserved = scope_new == "include_last",
    excluded_selection_from_both_source_scopes_identical = TRUE)
}
write_csv(bind_rows(selection_audit), file.path(audit_root, "selection_audit.csv"))

# Internal outcome rescaling improves numerical optimization without changing
# the model, raw-unit coefficients, standard errors, t statistics or p values.
outcome_scale <- 1000
fit_one <- function(dat, branch, template, scope) {
  dat$subject <- factor(dat$sub_id)
  dat$y_scaled <- dat$raw_PAC * outcome_scale
  fixed_terms <- c("segment_position", covariates[[branch]])
  formula <- as.formula(paste("y_scaled ~", paste(fixed_terms, collapse = " + "), "+ (1|subject)"))
  warnings <- character()
  fit <- withCallingHandlers(lme4::lmer(formula, dat, REML = TRUE,
    control = lme4::lmerControl(optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000, rhoend = 1e-9))),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  optimizer_code <- fit@optinfo$conv$opt
  convergence_messages <- fit@optinfo$conv$lme4$messages
  if (any(optimizer_code != 0) || length(convergence_messages)) {
    stop("Model convergence check failed: ", paste(c(optimizer_code, convergence_messages), collapse = "; "))
  }
  beta <- lme4::fixef(fit) / outcome_scale
  se <- sqrt(diag(as.matrix(vcov(fit)))) / outcome_scale
  X <- lme4::getME(fit, "X")
  df <- nrow(dat) - qr(X)$rank
  t <- beta / se
  p <- 2 * pt(-abs(t), df = df)
  ci_radius <- qt(0.975, df = df) * se
  vc <- as.data.frame(lme4::VarCorr(fit))
  tau2 <- vc$vcov[vc$grp == "subject"] / outcome_scale^2
  sigma2 <- sigma(fit)^2 / outcome_scale^2
  diag <- data.frame(branch = branch, scope = scope, roi = template$roi,
    phase_band = template$phase_band, amp_band = template$amp_band,
    exponent_range = if ("exponent_range" %in% names(template)) template$exponent_range else "",
    n_obs = nrow(dat), n_subjects = nlevels(dat$subject),
    n_segments = length(unique(dat$segment_position)), fixed_effect_rank = qr(X)$rank,
    residual_df = df, residual_variance = sigma2, subject_intercept_variance = tau2,
    is_singular = lme4::isSingular(fit, tol = 1e-4), optimizer_code = optimizer_code,
    warnings = paste(unique(warnings), collapse = "; "),
    model_formula = paste("raw_PAC ~", paste(fixed_terms, collapse = " + "), "+ (1|sub_id)"))
  fixed <- data.frame(term = names(beta), beta = unname(beta), se = unname(se),
    t_value = unname(t), df = df, p_value = unname(p),
    ci_low = unname(beta - ci_radius), ci_high = unname(beta + ci_radius))
  fixed <- cbind(diag[rep(1L, nrow(fixed)), c("branch", "scope", "roi", "phase_band", "amp_band", "exponent_range")], fixed)
  target <- fixed[fixed$term == "segment_position", ]
  target <- cbind(target, diag[setdiff(names(diag), names(target))])
  list(target = target, fixed = fixed, diagnostics = diag)
}
subset_model <- function(dat, template, branch) {
  keep <- dat$roi == template$roi & dat$phase_band == template$phase_band & dat$amp_band == template$amp_band
  if (branch == "aperiodic_control") keep <- keep & dat$exponent_range == template$exponent_range
  dat[keep, ]
}

baseline <- list()
baseline_fits <- list()
new_fits <- list()
for (branch in names(stages)) {
  templates <- original[[branch]]
  if (mode == "pilot") templates <- templates[1L, ]
  for (i in seq_len(nrow(templates))) {
    key <- paste(branch, i, sep = "__")
    template <- templates[i, ]
    message("Baseline ", key, " ", template$roi, " ", template$phase_band, " ", template$amp_band)
    one <- fit_one(subset_model(prepared[[branch]]$old, template, branch), branch,
                   template, source_scope)
    baseline_fits[[key]] <- one
    row <- one$target
    for (nm in c("beta", "se", "t_value", "df", "ci_low", "ci_high", "p_value", "q_value", "n_obs", "n_subjects", "n_segments")) {
      row[[paste0("original_", nm)]] <- template[[nm]]
    }
    baseline[[key]] <- row
  }
}
baseline <- bind_rows(baseline) |>
  group_by(branch) |> mutate(q_value = p.adjust(p_value, "BH", n = planned[first(branch)])) |> ungroup()
baseline <- baseline |> mutate(
  beta_difference = beta - original_beta, se_difference = se - original_se,
  p_difference = p_value - original_p_value, q_difference = q_value - original_q_value,
  beta_matches = abs(beta_difference) <= 1e-10 + 2e-4 * abs(original_beta),
  se_matches = abs(se_difference) <= 1e-10 + 2e-4 * abs(original_se),
  p_matches = abs(p_difference) <= 1e-4,
  coverage_matches = n_obs == original_n_obs & n_subjects == original_n_subjects & n_segments == original_n_segments,
  df_matches = df == original_df,
  raw_significance_matches = (p_value < 0.05) == (original_p_value < 0.05),
  fdr_significance_matches = (q_value < 0.05) == (original_q_value < 0.05))
write_csv(baseline, file.path(audit_root, paste0("baseline_matlab_vs_r_", mode, ".csv")))
passed <- all(baseline$beta_matches & baseline$se_matches & baseline$p_matches & baseline$coverage_matches & baseline$df_matches & baseline$raw_significance_matches)
if (mode == "full") passed <- passed && all(baseline$fdr_significance_matches) && max(abs(baseline$q_difference)) <= 1e-4
if (!passed) stop("Original-rule reproduction failed; inspect baseline audit before applying new exclusions")
if (!identical(input_md5, tools::md5sum(input_paths))) stop("Source input changed during baseline audit")
message("PASS: MATLAB-to-R baseline verified for ", nrow(baseline), " models")
if (mode == "pilot") quit(save = "no", status = 0)

for (branch in names(stages)) {
  templates <- original[[branch]]
  for (i in seq_len(nrow(templates))) {
    key <- paste(branch, i, sep = "__")
    if (scope_new == "include_last") {
      # The baseline IS the newly fitted R include_last model. Reuse it, rather
      # than fitting the identical complete-case observations a second time.
      message("Reuse newly fitted R include_last model ", key)
      one <- baseline_fits[[key]]
      for (element in c("target", "fixed", "diagnostics")) one[[element]]$scope <- scope_new
      new_fits[[key]] <- one
    } else {
      message("Adjusted ", key)
      new_fits[[key]] <- fit_one(subset_model(prepared[[branch]]$new, templates[i, ], branch),
                                 branch, templates[i, ], scope_new)
    }
  }
}
target <- apply_roi6_fdr(bind_rows(lapply(new_fits, `[[`, "target")))
if (nrow(target) != 78L || any(target$m_valid != 6L | target$m_planned != 6L)) stop("Incomplete ROI-wise test families")
target <- target |> mutate(fit_status = "ok", predictor = "segment_position", target_term = "segment_position",
                           raw_significant = p_value < 0.05, fdr_significant = q_value < 0.05)
fixed <- bind_rows(lapply(new_fits, `[[`, "fixed"))
diagnostics <- bind_rows(lapply(new_fits, `[[`, "diagnostics"))
comparison <- make_comparison(target, baseline)
write_csv(comparison, file.path(audit_root, "old_vs_adjusted_segment_effects.csv"))
persist_statistics(target)
for (analysis in c("main_pac", "power_control", "aperiodic_control")) {
  branches <- if (analysis == "main_pac") c("hpc_main", "amtg_main") else analysis
  out <- file.path(output_root, analysis)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv(fixed |> filter(branch %in% branches), file.path(out, "all_fixed_effects.csv"))
  write_csv(diagnostics |> filter(branch %in% branches), file.path(out, "model_diagnostics.csv"))
  write_csv(bind_rows(lapply(prepared[branches], `[[`, "new")), file.path(out, "analysis_dataset.csv"))
}
capture.output(sessionInfo(), file = file.path(output_root, "session_info.txt"))
manifest <- list(analysis = "Lina-derived adjusted PAC segment-number effects",
  scope = scope_new, included = "corridor 2.5-s halves and junctions; conditions pooled",
  terminal_exclusion = if (scope_new == "include_last") "none; retain terminal corridor halves and terminal junction" else
    "per subject-trial, remove C(turns+1)_first, C(turns+1)_second and J(turns)",
  segment_numbering = "original navigation sequence retained, no valid-row reindexing",
  source_input_paths = as.list(input_paths), source_input_md5 = as.list(input_md5),
  estimator = "Gaussian linear mixed model, subject random intercept, REML, lme4::lmer",
  internal_outcome_scale = outcome_scale, reported_coefficient_scale = "raw Tort MI units",
  inference = "two-sided Wald t with residual DF n - rank(X), reproducing MATLAB fitlme default",
  fdr = fdr_policy,
  exponent_ranges_fitted = c("1-30 Hz", "1-150 Hz"), exponent_range_primary_plot = "1-150 Hz",
  matlab_baseline_models = nrow(baseline), adjusted_models = nrow(target),
  matlab_baseline_scope = source_scope,
  matlab_baseline_fdr_only_for_implementation_audit = as.list(planned),
  maximum_baseline_p_difference = max(abs(baseline$p_difference)),
  maximum_baseline_q_difference = max(abs(baseline$q_difference)),
  source_files_modified = FALSE, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
if (!identical(input_md5, tools::md5sum(input_paths))) stop("Source files changed during refit")
jsonlite::write_json(manifest, file.path(output_root, "analysis_manifest.json"),
                      pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("78 R models, MATLAB implementation audit, selection audit and 13 ROI/range-specific six-test BH families passed.",
           file.path(output_root, "_SUCCESS"))
message("Completed adjusted PAC segment analysis: ", output_root)
