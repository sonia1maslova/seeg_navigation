# Shared helpers for the independent position-to-start bearing sensitivity branch.
# Existing analysis scripts are PARSED, not sourced: only explicitly selected
# function/constant assignments are evaluated, so their original output code
# never runs. No packages are installed and no original artifact is written.

branch_name <- "bearing_position_to_start"
branch_rois_plot <- c("Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
                      "Hippocampus L_Posterior", "Hippocampus R_Anterior")
branch_bands <- c("f0p5_4", "f4_12", "f12_30", "f30_70", "f70_150")
branch_predictors <- list(corridor = c("Angle2Xaxis", "Angle2orientation"),
                         junction = "Angle2Xaxis", all_segments = "Angle2Xaxis")

branch_script_path <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Run this entry point with Rscript.", call. = FALSE)
  normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)
}

branch_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required packages unavailable: ", paste(missing, collapse = ", "), call. = FALSE)
}

branch_config <- function() {
  root <- Sys.getenv("SEEG_NAVIGATION_ROOT", unset = "/nfs/visu/workplace/yuhdua/projects/seeg_navigation")
  args <- commandArgs(TRUE)
  if (length(args) > 1L) stop("Usage: Rscript entry.R [project_root]", call. = FALSE)
  if (length(args)) root <- args[1]
  root <- normalizePath(root, mustWork = TRUE)
  Sys.setenv(SEEG_NAVIGATION_ROOT = root)
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "codes", "r", "project_config.R"), envir = env)
  env$cfg
}

branch_load_assignments <- function(path, functions, constants = character()) {
  ast <- parse(path, keep.source = FALSE)
  selected <- c(functions, constants)
  env <- new.env(parent = globalenv())
  found <- character()
  for (expr in ast) {
    if (!is.call(expr) || !identical(expr[[1]], as.name("<-")) || !is.symbol(expr[[2]])) next
    name <- as.character(expr[[2]])
    if (!name %in% selected) next
    if (name %in% found) stop("Duplicate canonical assignment: ", name, call. = FALSE)
    is_fun <- is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
    if (name %in% functions && !is_fun) stop("Expected a canonical function: ", name, call. = FALSE)
    if (name %in% constants && is_fun) stop("Expected a canonical constant: ", name, call. = FALSE)
    eval(expr, envir = env)
    found <- c(found, name)
  }
  if (length(setdiff(selected, found))) stop("Canonical definitions missing: ", paste(setdiff(selected, found), collapse = ", "), call. = FALSE)
  env
}

branch_fdr <- function(dat) {
  # Same established ANGLE heatmap families, separately for each predictor,
  # scope, and last-segment branch. pMTG has no established angle FDR family.
  dat$p_fdr <- NA_real_
  dat$fdr_family <- ifelse(dat$roi == "Temporal_Mid L_Anterior", "aMTG_5_bands",
                          ifelse(dat$roi %in% branch_rois_plot[-1], "hippocampus_3_ROIs_15_bands", "not_defined_original_angle_plot"))
  dat$fdr_planned_tests <- ifelse(dat$roi == "Temporal_Mid L_Anterior", 5L,
                                 ifelse(dat$roi %in% branch_rois_plot[-1], 15L, NA_integer_))
  families <- unique(dat[c("scope", "variant", "predictor", "fdr_family")])
  for (i in seq_len(nrow(families))) {
    f <- families[i, ]
    idx <- which(dat$scope == f$scope & dat$variant == f$variant &
                   dat$predictor == f$predictor & dat$fdr_family == f$fdr_family)
    n <- unique(dat$fdr_planned_tests[idx])
    if (length(n) == 1L && !is.na(n)) {
      if (length(idx) > n) stop("More results than planned FDR family tests.", call. = FALSE)
      dat$p_fdr[idx] <- stats::p.adjust(dat$p_value_raw[idx], method = "BH", n = n)
    }
  }
  dat
}

branch_key <- function(dat, keys) {
  if (any(!keys %in% names(dat))) stop("Audit key columns absent.", call. = FALSE)
  values <- lapply(dat[keys], function(x) ifelse(is.na(x), "<NA>", as.character(x)))
  do.call(paste, c(values, sep = "\037"))
}

branch_same_data <- function(actual, expected, keys, columns, label, tolerance = 1e-10) {
  a <- branch_key(actual, keys); b <- branch_key(expected, keys)
  if (anyDuplicated(a) || anyDuplicated(b) || length(a) != length(b) || !setequal(a, b)) {
    stop(label, ": row keys/counts differ.", call. = FALSE)
  }
  actual <- actual[match(b, a), , drop = FALSE]
  for (name in columns) {
    if (!name %in% names(actual) || !name %in% names(expected)) stop(label, ": missing ", name, call. = FALSE)
    x <- actual[[name]]; y <- expected[[name]]
    if (is.numeric(x) && is.numeric(y)) {
      ok <- identical(is.na(x), is.na(y)) && identical(is.finite(x), is.finite(y)) &&
        all(abs(x[is.finite(x)] - y[is.finite(y)]) <= tolerance)
    } else ok <- identical(as.character(x), as.character(y))
    if (!isTRUE(ok)) stop(label, ": values differ in ", name, call. = FALSE)
  }
  invisible(TRUE)
}

branch_check_reversal <- function(original, reversed, keys, angles, label) {
  branch_same_data(reversed, original, keys, setdiff(names(original), angles), label)
  rev_idx <- match(branch_key(original, keys), branch_key(reversed, keys))
  for (angle in angles) {
    x <- original[[angle]]; y <- reversed[[angle]][rev_idx]
    if (!identical(is.finite(x), is.finite(y))) stop(label, ": finite angle pattern changed.", call. = FALSE)
    k <- is.finite(x)
    if (any(x[k] < 0 | x[k] >= 2 * pi + 1e-10)) stop(label, ": original angle outside 0..2pi.", call. = FALSE)
    expected <- (x[k] + pi) %% (2 * pi)
    # Canonical numeric conversion rounds through as.character(). Values on
    # the pi boundary can land on opposite sides by a few floating-point ulps,
    # yielding 0 versus almost 2*pi. Compare circular distance, not linear gap.
    circular_error <- abs(atan2(sin(expected - y[k]), cos(expected - y[k])))
    if (any(circular_error > 1e-10) || any(y[k] < 0 | y[k] >= 2 * pi + 1e-10)) {
      stop(label, ": reversed angle formula mismatch.", call. = FALSE)
    }
  }
}

branch_write_csv <- function(dat, path) {
  if (file.exists(path)) stop("Refusing to overwrite branch output: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".bearing_", tmpdir = dirname(path), fileext = ".csv")
  readr::write_csv(dat, temporary, na = "")
  if (!file.rename(temporary, path)) stop("Cannot commit branch output: ", path, call. = FALSE)
}

branch_write_lines <- function(lines, path) {
  if (file.exists(path)) stop("Refusing to overwrite branch output: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_lines(lines, path)
}

branch_assert_isolated <- function(path, expected_parent) {
  resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent <- normalizePath(expected_parent, winslash = "/", mustWork = FALSE)
  expected <- paste0(parent, "/", branch_name)
  if (!identical(resolved, expected)) stop("Unexpected supplemental output root: ", resolved, call. = FALSE)
}
