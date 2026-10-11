# =============================================================================
# S02_run_power_encoding_lmm.R
#
# Formal Power encoding LMMs using existing ROI-averaged segment-level Power.
# Predictors are used in RAW ORIGINAL UNITS (no z-transformation).
#
# Primary models
# --------------
# corridor / junction:
#   power_db ~ predictor + condition + (1 | sub_id)
#
# all_segments:
#   power_db ~ predictor + segment_type + condition + (1 | sub_id)
#
# condition is modeled as a FIXED effect. If only one condition level is
# present in a particular model, the condition fixed effect is omitted.
#
# Representation predictors are entered in their ORIGINAL UNITS.
# No z-transformation is performed.
# Power (power_db) is also NOT standardized.
#
# Formal predictors
# -----------------
# corridor     : Distance, SegmentNumber, Angle2Xaxis, Angle2orientation
# junction     : Distance, SegmentNumber, Angle2Xaxis
# all_segments : Distance, SegmentNumber, Angle2Xaxis
#
# Variants
# --------
# include_last : retain every available segment
# exclude_last : remove the condition-specific final corridor and final junction
#
# Multiple comparisons
# --------------------
# NO FDR or other multiple-comparison correction is performed here.
# The significant summary uses the raw target-predictor p value < 0.05.
# Run with --supplementary-cortical to fit only corridor egocentric-bearing
# effects for whole left STG, SMG, and Precentral in a separate result folder.
# =============================================================================

required_packages <- c("dplyr", "readr", "tibble", "lme4", "lmerTest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(readr)
library(tibble)
library(lme4)
library(lmerTest)

# =============================================================================
# 1. Locate project_config.R
# =============================================================================

locate_this_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) if (is.null(x$ofile)) NA_character_ else as.character(x$ofile),
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)], mustWork = FALSE))
  }

  stop("Cannot determine current script path. Run with Rscript.", call. = FALSE)
}

script_file <- locate_this_script()
script_dir <- dirname(script_file)
r_code_root <- dirname(dirname(script_dir))
project_config_file <- file.path(r_code_root, "project_config.R")
if (!file.exists(project_config_file)) {
  stop("project_config.R not found:\n", project_config_file, call. = FALSE)
}
source(project_config_file)

# =============================================================================
# 2. Formal analysis definition
# =============================================================================

SELECTED_ROIS <- c(
  "Hippocampus L_Anterior",
  "Hippocampus L_Posterior",
  "Hippocampus R_Anterior",
  "Temporal_Mid L_Anterior",
 "Temporal_Mid L_Posterior"
)

FORMAL_BANDS <- c("f0p5_4", "f4_12", "f12_30", "f30_70", "f70_150")
BAND_LABELS <- c(
  f0p5_4 = "0.5-4 Hz",
  f4_12 = "4-12 Hz",
  f12_30 = "12-30 Hz",
  f30_70 = "30-70 Hz",
  f70_150 = "70-150 Hz"
)

PREDICTORS_BY_SCOPE <- list(
  corridor = c("Distance", "SegmentNumber", "Angle2Xaxis", "Angle2orientation"),
  junction = c("Distance", "SegmentNumber", "Angle2Xaxis"),
  all_segments = c("Distance", "SegmentNumber", "Angle2Xaxis")
)

VARIANTS <- c("include_last", "exclude_last")
SCOPES <- names(PREDICTORS_BY_SCOPE)

# Optional, narrowly scoped extension. The default five-ROI analysis is
# unchanged; this mode writes separately and never replaces its results.
SUPPLEMENTARY_CORTICAL <- "--supplementary-cortical" %in% commandArgs(trailingOnly = TRUE)
if (SUPPLEMENTARY_CORTICAL) {
  SELECTED_ROIS <- c("Temporal_Sup L", "SupraMarginal L", "Precentral L")
  PREDICTORS_BY_SCOPE <- list(corridor = "Angle2orientation")
  SCOPES <- names(PREDICTORS_BY_SCOPE)
}

# Model settings. No scientific thresholding is imposed beyond what is needed
# to make a mixed model mathematically meaningful.
USE_REML <- FALSE
MIN_SUBJECTS <- 2L
MIN_OBSERVATIONS <- 8L
SINGULAR_TOLERANCE <- 1e-4
OPTIMIZER <- "bobyqa"
MAXFUN <- 200000L
RAW_SIGNIFICANCE_ALPHA <- 0.05

prepared_root <- file.path(cfg$data$derived, "power", "encoding_lmm")
results_root <- file.path(cfg$results$power, "encoding_lmm")
if (SUPPLEMENTARY_CORTICAL) {
  results_root <- file.path(results_root, "supplementary_cortical_rois")
}
by_roi_root <- file.path(results_root, "by_roi")
summary_root <- file.path(results_root, "summaries")
qc_root <- file.path(results_root, "qc")

dir.create(by_roi_root, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

condition_order <- as.character(cfg$condition_order)

# =============================================================================
# 3. Helpers
# =============================================================================

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

make_condition_factor <- function(x) {
  x_chr <- as.character(x)
  present <- condition_order[condition_order %in% unique(x_chr)]
  extra <- setdiff(sort(unique(x_chr)), present)
  factor(x_chr, levels = c(present, extra))
}

get_convergence_message <- function(model) {
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) return(NA_character_)
  paste(msg, collapse = " | ")
}

collapse_warnings <- function(x) {
  x <- unique(x[nzchar(x)])
  if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
}

extract_random_effect_info <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  subject_var <- NA_real_
  residual_var <- NA_real_

  hit_sub <- which(vc$grp == "sub_id" & is.na(vc$var2))
  if (length(hit_sub) > 0) subject_var <- as.numeric(vc$vcov[hit_sub[1]])

  hit_res <- which(vc$grp == "Residual")
  if (length(hit_res) > 0) residual_var <- as.numeric(vc$vcov[hit_res[1]])

  c(subject_random_intercept_variance = subject_var, residual_variance = residual_var)
}

extract_coef_table <- function(model) {
  tab <- as.data.frame(summary(model)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  get_numeric_col <- function(candidates) {
    hit <- intersect(candidates, names(tab))
    if (length(hit) == 0) return(rep(NA_real_, nrow(tab)))
    suppressWarnings(as.numeric(tab[[hit[1]]]))
  }

  tibble::tibble(
    term = as.character(tab$term),
    estimate = get_numeric_col(c("Estimate")),
    std_error = get_numeric_col(c("Std. Error", "Std..Error")),
    df = get_numeric_col(c("df")),
    statistic = get_numeric_col(c("t value", "t.value")),
    p_value_raw = get_numeric_col(c("Pr(>|t|)", "Pr...t.."))
  ) |>
    dplyr::mutate(
      conf_low_95 = .data$estimate - 1.96 * .data$std_error,
      conf_high_95 = .data$estimate + 1.96 * .data$std_error
    )
}

fit_one_model <- function(dat, roi, scope, variant, band_name, predictor) {
  dat <- dat |>
    dplyr::mutate(
      sub_id = factor(as.character(.data$sub_id)),
      condition = make_condition_factor(.data$condition),
      segment_type = factor(as.character(.data$segment_type), levels = c("corridor", "junction")),
      power_db = suppressWarnings(as.numeric(.data$power_db)),
      predictor_value = suppressWarnings(as.numeric(.data[[predictor]]))
    ) |>
    dplyr::filter(
      is.finite(.data$power_db),
      is.finite(.data$predictor_value),
      !is.na(.data$sub_id),
      !is.na(.data$condition)
    ) |>
    droplevels()

  n_obs <- nrow(dat)
  n_subjects <- dplyr::n_distinct(dat$sub_id)
  n_trials <- nrow(dplyr::distinct(dat, .data$sub_id, .data$trial))
  n_conditions <- nlevels(dat$condition)
  n_segment_types <- nlevels(dat$segment_type)
  predictor_mean <- if (n_obs > 0) mean(dat$predictor_value) else NA_real_
  predictor_sd <- if (n_obs > 1) stats::sd(dat$predictor_value) else NA_real_
  predictor_min <- if (n_obs > 0) min(dat$predictor_value) else NA_real_
  predictor_max <- if (n_obs > 0) max(dat$predictor_value) else NA_real_
  power_mean_db <- if (n_obs > 0) mean(dat$power_db) else NA_real_
  power_sd_db <- if (n_obs > 1) stats::sd(dat$power_db) else NA_real_

  band_low <- if (n_obs > 0) unique(dat$band_low_hz)[1] else NA_real_
  band_high <- if (n_obs > 0) unique(dat$band_high_hz)[1] else NA_real_

  model_id <- paste(roi, scope, variant, band_name, predictor, sep = "__")

  base <- tibble::tibble(
    model_id = model_id,
    roi = roi,
    scope = scope,
    variant = variant,
    band_name = band_name,
    band_label = unname(BAND_LABELS[band_name]),
    band_low_hz = as.numeric(band_low),
    band_high_hz = as.numeric(band_high),
    predictor = predictor,
    predictor_scaling = "raw_original_units",
    n_observations = n_obs,
    n_subjects = n_subjects,
    n_subject_trials = n_trials,
    n_conditions = n_conditions,
    n_segment_types = n_segment_types,
    predictor_mean = predictor_mean,
    predictor_sd = predictor_sd,
    predictor_min = predictor_min,
    predictor_max = predictor_max,
    power_mean_db = power_mean_db,
    power_sd_db = power_sd_db
  )

  if (n_subjects < MIN_SUBJECTS || n_obs < MIN_OBSERVATIONS) {
    qc <- base |>
      dplyr::mutate(
        model_formula = NA_character_,
        singular_fit = NA,
        convergence_message = NA_character_,
        warning_message = NA_character_,
        fit_status = "insufficient_data"
      )
    return(list(fixed = NULL, target = NULL, qc = qc))
  }

  if (!is.finite(predictor_sd) || predictor_sd == 0) {
    qc <- base |>
      dplyr::mutate(
        model_formula = NA_character_,
        singular_fit = NA,
        convergence_message = NA_character_,
        warning_message = NA_character_,
        fit_status = "zero_predictor_variance"
      )
    return(list(fixed = NULL, target = NULL, qc = qc))
  }

  # Use the predictor directly in its original units.
  # The dependent variable power_db remains in dB and is NOT standardized.
  fixed_terms <- c("predictor_value")
  if (scope == "all_segments" && n_segment_types > 1) {
    fixed_terms <- c(fixed_terms, "segment_type")
  }
  if (n_conditions > 1) {
    fixed_terms <- c(fixed_terms, "condition")
  }

  random_terms <- c("(1 | sub_id)")

  formula_text <- paste0(
    "power_db ~ ",
    paste(fixed_terms, collapse = " + "),
    " + ",
    paste(random_terms, collapse = " + ")
  )
  formula_obj <- stats::as.formula(formula_text)

  warning_messages <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        formula_obj,
        data = dat,
        REML = USE_REML,
        control = lme4::lmerControl(
          optimizer = OPTIMIZER,
          optCtrl = list(maxfun = MAXFUN)
        )
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    qc <- base |>
      dplyr::mutate(
        model_formula = formula_text,
        singular_fit = NA,
        convergence_message = conditionMessage(fit),
        warning_message = collapse_warnings(warning_messages),
        fit_status = "fit_error"
      )
    return(list(fixed = NULL, target = NULL, qc = qc))
  }

  singular_fit <- lme4::isSingular(fit, tol = SINGULAR_TOLERANCE)
  convergence_message <- get_convergence_message(fit)
  warning_message <- collapse_warnings(warning_messages)
  random_info <- extract_random_effect_info(fit)

  coef_table <- extract_coef_table(fit)
  if (!"predictor_value" %in% coef_table$term) {
    qc <- base |>
      dplyr::mutate(
        model_formula = formula_text,
        singular_fit = singular_fit,
        convergence_message = convergence_message,
        warning_message = warning_message,
        fit_status = "target_term_missing"
      )
    return(list(fixed = NULL, target = NULL, qc = qc))
  }

  fixed <- base[rep(1, nrow(coef_table)), , drop = FALSE] |>
    dplyr::bind_cols(coef_table) |>
    dplyr::mutate(
      term_display = dplyr::case_when(
        .data$term == "predictor_value" ~ predictor,
        TRUE ~ .data$term
      ),
      is_target_predictor = .data$term == "predictor_value",
      model_formula = formula_text,
      reml = USE_REML,
      optimizer = OPTIMIZER,
      singular_fit = singular_fit,
      convergence_message = convergence_message,
      warning_message = warning_message,
      subject_random_intercept_variance = unname(random_info["subject_random_intercept_variance"]),
      residual_variance = unname(random_info["residual_variance"]),
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit),
      logLik = as.numeric(stats::logLik(fit)),
      fit_status = "ok"
    )

  target <- fixed |>
    dplyr::filter(.data$is_target_predictor) |>
    dplyr::transmute(
      dplyr::across(dplyr::everything())
    )

  qc <- base |>
    dplyr::mutate(
      model_formula = formula_text,
      singular_fit = singular_fit,
      convergence_message = convergence_message,
      warning_message = warning_message,
      fit_status = "ok"
    )

  list(fixed = fixed, target = target, qc = qc)
}

# =============================================================================
# 4. Run all models
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("S02 | Power encoding LMM\n")
cat(strrep("=", 78), "\n", sep = "")
cat("Dependent variable : power_db\n")
cat("Predictor scaling  : raw original units; NO z-transformation; power_db NOT standardized\n")
cat("Condition          : fixed effect\n")
cat("Subject            : random intercept\n")
cat("Estimation method  : ML (REML = FALSE)\n")
cat("Analysis mode      : ", if (SUPPLEMENTARY_CORTICAL) "supplementary cortical" else "original five ROI", "\n", sep = "")
if ("all_segments" %in% SCOPES) cat("All segments       : additionally controls segment_type\n")
cat("Bands              : 0.5-4, 4-12, 12-30, 30-70, 70-150 Hz\n")
cat("Variants           : include_last, exclude_last\n")
cat("Multiple correction: NONE\n")
cat("Significant summary: raw p < ", RAW_SIGNIFICANCE_ALPHA, "\n", sep = "")
cat(strrep("=", 78), "\n\n", sep = "")

all_target <- list()
all_qc <- list()
target_i <- 1L
qc_i <- 1L

for (roi in SELECTED_ROIS) {
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat("ROI: ", roi, "\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")

  roi_fixed <- list()
  roi_fixed_i <- 1L
  roi_model_count <- 0L

  for (scope in SCOPES) {
    input_file <- file.path(prepared_root, scope, paste0(safe_filename(roi), ".csv"))
    if (!file.exists(input_file)) {
      stop(
        "Prepared input not found:\n", input_file,
        "\nRun S01_prepare_power_encoding_data.R first.",
        call. = FALSE
      )
    }

    dat_scope <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)

    required <- c(
      "sub_id", "trial", "condition", "roi", "segment_type", "band_name",
      "band_low_hz", "band_high_hz", "power_db", "is_last_segment",
      PREDICTORS_BY_SCOPE[[scope]]
    )
    missing <- setdiff(required, names(dat_scope))
    if (length(missing) > 0) {
      stop(
        "Prepared file missing columns: ", paste(missing, collapse = ", "),
        "\n", input_file,
        call. = FALSE
      )
    }

    for (variant in VARIANTS) {
      dat_variant <- dat_scope
      if (variant == "exclude_last") {
        dat_variant <- dat_variant |>
          dplyr::filter(!(.data$is_last_segment %in% TRUE))
      }

      cat("  ", scope, " / ", variant, " ...\n", sep = "")

      for (band_name in FORMAL_BANDS) {
        dat_band <- dat_variant |>
          dplyr::filter(.data$band_name == .env$band_name)

        # Defensive QC: each model must contain exactly the requested band.
        observed_bands <- unique(as.character(dat_band$band_name))
        if (length(observed_bands) != 1L || observed_bands[1] != band_name) {
          stop(
            "Band filtering failed for ROI=", roi,
            ", scope=", scope,
            ", variant=", variant,
            ", requested band=", band_name,
            "; observed bands=", paste(observed_bands, collapse = ", "),
            call. = FALSE
          )
        }

        for (predictor in PREDICTORS_BY_SCOPE[[scope]]) {
          fit <- fit_one_model(
            dat = dat_band,
            roi = roi,
            scope = scope,
            variant = variant,
            band_name = band_name,
            predictor = predictor
          )

          roi_model_count <- roi_model_count + 1L

          if (!is.null(fit$fixed)) {
            roi_fixed[[roi_fixed_i]] <- fit$fixed
            roi_fixed_i <- roi_fixed_i + 1L
          }
          if (!is.null(fit$target)) {
            all_target[[target_i]] <- fit$target
            target_i <- target_i + 1L
          }
          all_qc[[qc_i]] <- fit$qc
          qc_i <- qc_i + 1L
        }
      }
    }
  }

  roi_detail <- if (length(roi_fixed) > 0) {
    dplyr::bind_rows(roi_fixed) |>
      dplyr::arrange(
        factor(.data$scope, levels = SCOPES),
        factor(.data$variant, levels = VARIANTS),
        factor(.data$band_name, levels = FORMAL_BANDS),
        .data$predictor,
        dplyr::desc(.data$is_target_predictor),
        .data$term
      )
  } else {
    tibble::tibble()
  }

  roi_out <- file.path(
    by_roi_root,
    paste0(safe_filename(roi), "_lmm_detailed.csv")
  )
  readr::write_csv(roi_detail, roi_out, na = "")

  n_sig_roi <- if (nrow(roi_detail) > 0) {
    sum(
      roi_detail$is_target_predictor %in% TRUE &
        is.finite(roi_detail$p_value_raw) &
        roi_detail$p_value_raw < RAW_SIGNIFICANCE_ALPHA
    )
  } else 0L

  cat(
    "  Saved: ", roi_out,
    "\n  Models attempted: ", roi_model_count,
    " | significant target effects (raw p<0.05): ", n_sig_roi,
    "\n",
    sep = ""
  )
}

# =============================================================================
# 5. Save global summaries
# =============================================================================

target_effects <- if (length(all_target) > 0) {
  dplyr::bind_rows(all_target) |>
    dplyr::arrange(
      factor(.data$roi, levels = SELECTED_ROIS),
      factor(.data$scope, levels = SCOPES),
      factor(.data$variant, levels = VARIANTS),
      factor(.data$band_name, levels = FORMAL_BANDS),
      .data$predictor
    )
} else {
  tibble::tibble()
}

model_qc <- dplyr::bind_rows(all_qc) |>
  dplyr::arrange(
    factor(.data$roi, levels = SELECTED_ROIS),
    factor(.data$scope, levels = SCOPES),
    factor(.data$variant, levels = VARIANTS),
    factor(.data$band_name, levels = FORMAL_BANDS),
    .data$predictor
  )

all_target_file <- file.path(summary_root, "all_target_predictor_effects.csv")
sig_file <- file.path(
  summary_root,
  "significant_target_effects_uncorrected_p_lt_0p05.csv"
)
sig_compact_file <- file.path(
  summary_root,
  "significant_target_effects_uncorrected_p_lt_0p05_compact.csv"
)
qc_file <- file.path(qc_root, "model_qc.csv")

readr::write_csv(target_effects, all_target_file, na = "")
readr::write_csv(model_qc, qc_file, na = "")

significant <- if (nrow(target_effects) > 0) {
  target_effects |>
    dplyr::filter(
      .data$fit_status == "ok",
      is.finite(.data$p_value_raw),
      .data$p_value_raw < RAW_SIGNIFICANCE_ALPHA
    ) |>
    dplyr::arrange(.data$p_value_raw)
} else {
  target_effects
}

readr::write_csv(significant, sig_file, na = "")

sig_compact <- if (nrow(significant) > 0) {
  significant |>
    dplyr::transmute(
      roi = .data$roi,
      scope = .data$scope,
      variant = .data$variant,
      band = .data$band_label,
      predictor = .data$predictor,
      beta = .data$estimate,
      SE = .data$std_error,
      df = .data$df,
      t = .data$statistic,
      p_raw = .data$p_value_raw,
      n = .data$n_observations,
      n_subjects = .data$n_subjects,
      singular_fit = .data$singular_fit
    )
} else {
  tibble::tibble(
    roi = character(), scope = character(), variant = character(),
    band = character(), predictor = character(), beta = numeric(),
    SE = numeric(), df = numeric(), t = numeric(), p_raw = numeric(),
    n = integer(), n_subjects = integer(), singular_fit = logical()
  )
}
readr::write_csv(sig_compact, sig_compact_file, na = "")

# =============================================================================
# 6. Print significant raw-p summary
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SIGNIFICANT TARGET EFFECT SUMMARY | uncorrected raw p < 0.05\n")
cat(strrep("=", 78), "\n", sep = "")

if (nrow(sig_compact) == 0) {
  cat("No target-predictor effects reached raw p < 0.05.\n")
} else {
  print(sig_compact, n = Inf, width = Inf)

  cat("\nCounts by ROI:\n")
  print(
    significant |>
      dplyr::count(.data$roi, sort = TRUE, name = "n_significant"),
    n = Inf,
    width = Inf
  )

  cat("\nCounts by scope / variant:\n")
  print(
    significant |>
      dplyr::count(.data$scope, .data$variant, sort = TRUE, name = "n_significant"),
    n = Inf,
    width = Inf
  )

  cat("\nCounts by predictor:\n")
  print(
    significant |>
      dplyr::count(.data$predictor, sort = TRUE, name = "n_significant"),
    n = Inf,
    width = Inf
  )

  n_sig_singular <- sum(significant$singular_fit %in% TRUE, na.rm = TRUE)
  cat("\nSignificant effects from singular fits: ", n_sig_singular, "\n", sep = "")
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("S02 completed\n")
cat("Detailed ROI files : ", by_roi_root, "\n", sep = "")
cat("All target effects : ", all_target_file, "\n", sep = "")
cat("Significant summary: ", sig_compact_file, "\n", sep = "")
cat("Model QC           : ", qc_file, "\n", sep = "")
cat("FDR                : NOT performed\n")
cat(strrep("=", 78), "\n", sep = "")
