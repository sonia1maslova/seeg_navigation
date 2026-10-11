# =============================================================================
# S01_prepare_power_encoding_data.R
#
# Prepare the existing ROI-averaged segment-level Power for encoding LMMs.
#
# IMPORTANT
# ---------
# - Power is NOT recomputed here.
# - Input Power comes from the new project pipeline:
#     data/derived/power/segment_power/
# - Only ROI/contact-averaged Power is used.
# - Corridor analysis uses corridor/half_2p5s.
# - Junction analysis uses junction/whole_2p5s.
# - all_segments pools the two 2.5-s observation types above.
# - Distance / Angle values are attached from the canonical encoding tables.
# - SegmentNumber is created directly from navigation order.
# - No predictor is z-scored.
# =============================================================================

required_packages <- c("dplyr", "readr", "tibble")
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
# codes/r/power/encoding_lmm -> codes/r
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
  "Hippocampus R_Anterior",
  "Hippocampus L_Anterior",
  "Hippocampus L_Posterior",
  "Temporal_Mid L_Anterior",
  "Temporal_Mid L_Posterior",
  "Temporal_Sup L",
  "ParaHippocampal L",
  "ParaHippocampal R",
  "Precentral L",
  "SupraMarginal L"
)

FORMAL_BANDS <- tibble::tribble(
  ~band_name, ~band_low_hz, ~band_high_hz,
  "f0p5_4",          0.5,           4,
  "f4_12",             4,          12,
  "f12_30",           12,          30,
  "f30_70",           30,          70,
  "f70_150",          70,         150
)

power_root <- file.path(cfg$data$derived, "power", "segment_power")
corridor_encoding_file <- file.path(
  cfg$data$encoding,
  "corridor_subsegi_encoding_table.csv"
)
junction_encoding_file <- file.path(
  cfg$data$encoding,
  "turn_encoding_table.csv"
)

output_root <- file.path(cfg$data$derived, "power", "encoding_lmm")
qc_root <- file.path(output_root, "qc")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 3. Helpers
# =============================================================================

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

find_col <- function(dat, candidates, required = TRUE) {
  nm <- names(dat)
  for (candidate in candidates) {
    hit <- which(tolower(nm) == tolower(candidate))
    if (length(hit) > 0) return(nm[hit[1]])
  }
  if (required) {
    stop(
      "Required column not found. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_half <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  txt <- tolower(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(txt))
  out[grepl("first", txt) | txt %in% c("f", "h1", "first_half")] <- 1
  out[grepl("second", txt) | txt %in% c("s", "h2", "second_half")] <- 2
  out
}

n_junctions_from_condition <- function(condition) {
  out <- rep(NA_integer_, length(condition))
  condition <- as.character(condition)
  out[condition == "turns_num2"] <- 2L
  out[condition == "turns_num4"] <- 4L
  out[condition == "turns_num6"] <- 6L
  out
}

assert_unique_keys <- function(dat, keys, label) {
  dup <- dat |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup) > 0) {
    preview <- utils::capture.output(print(utils::head(dup, 10)))
    stop(
      label, " contains duplicated keys: ", paste(keys, collapse = ", "), "\n",
      paste(preview, collapse = "\n"),
      call. = FALSE
    )
  }
}

read_corridor_encoding <- function(path) {
  if (!file.exists(path)) stop("Corridor encoding table not found:\n", path, call. = FALSE)
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

  trial_col <- find_col(raw, c("Trial", "trial"))
  seg_col <- find_col(raw, c("Corridor", "corridor", "corridor_number", "segment_no"))
  half_col <- find_col(raw, c("Half", "half", "part_no", "part_name"))
  dist_col <- find_col(raw, c("Distance", "distance"))
  anglex_col <- find_col(raw, c("Angle2Xaxis", "angle2xaxis"))
  ango_col <- find_col(raw, c("Angle2orientation", "angle2orientation"))

  out <- tibble::tibble(
    trial = as_numeric_safe(raw[[trial_col]]),
    segment_no = as_numeric_safe(raw[[seg_col]]),
    part_no = normalize_half(raw[[half_col]]),
    Distance = as_numeric_safe(raw[[dist_col]]),
    Angle2Xaxis = as_numeric_safe(raw[[anglex_col]]),
    Angle2orientation = as_numeric_safe(raw[[ango_col]])
  ) |>
    dplyr::filter(
      is.finite(.data$trial),
      is.finite(.data$segment_no),
      .data$part_no %in% c(1, 2)
    ) |>
    dplyr::mutate(
      trial = as.integer(.data$trial),
      segment_no = as.integer(.data$segment_no),
      part_no = as.integer(.data$part_no)
    )

  assert_unique_keys(out, c("trial", "segment_no", "part_no"), "Corridor encoding table")
  out
}

read_junction_encoding <- function(path) {
  if (!file.exists(path)) stop("Junction encoding table not found:\n", path, call. = FALSE)
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

  trial_col <- find_col(raw, c("Trial", "trial"))
  seg_col <- find_col(raw, c("Turn", "turn", "turn_number", "segment_no"))
  dist_col <- find_col(raw, c("Distance", "distance"))
  anglex_col <- find_col(raw, c("Angle2Xaxis", "angle2xaxis"))

  out <- tibble::tibble(
    trial = as_numeric_safe(raw[[trial_col]]),
    segment_no = as_numeric_safe(raw[[seg_col]]),
    Distance = as_numeric_safe(raw[[dist_col]]),
    Angle2Xaxis = as_numeric_safe(raw[[anglex_col]])
  ) |>
    dplyr::filter(
      is.finite(.data$trial),
      is.finite(.data$segment_no)
    ) |>
    dplyr::mutate(
      trial = as.integer(.data$trial),
      segment_no = as.integer(.data$segment_no)
    )

  assert_unique_keys(out, c("trial", "segment_no"), "Junction encoding table")
  out
}

power_file_for <- function(roi, process_type) {
  roi_tag <- safe_filename(roi)
  if (process_type == "corridor") {
    mode <- "half_2p5s"
  } else if (process_type == "junction") {
    mode <- "whole_2p5s"
  } else {
    stop("Unsupported process_type: ", process_type, call. = FALSE)
  }

  file.path(
    power_root,
    process_type,
    mode,
    roi_tag,
    "contact_averaged",
    paste0(
      "power_", process_type, "_", mode, "_", roi_tag,
      "_contact_averaged.csv"
    )
  )
}

prepare_power_common <- function(power_file, roi_expected, process_expected) {
  if (!file.exists(power_file)) {
    stop("Power file not found:\n", power_file, call. = FALSE)
  }

  dat <- readr::read_csv(power_file, show_col_types = FALSE, progress = FALSE)

  required <- c(
    "sub_id", "roi", "trial", "condition", "process_type", "segment_mode",
    "segment_no", "part_no", "part_name", "band_name", "band_low_hz",
    "band_high_hz", "power_db"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(
      "Power file is missing columns: ", paste(missing, collapse = ", "),
      "\n", power_file,
      call. = FALSE
    )
  }

  dat <- dat |>
    dplyr::mutate(
      sub_id = as.character(.data$sub_id),
      roi = as.character(.data$roi),
      trial = as.integer(as_numeric_safe(.data$trial)),
      condition = as.character(.data$condition),
      process_type = as.character(.data$process_type),
      segment_no = as.integer(as_numeric_safe(.data$segment_no)),
      part_no = as.integer(as_numeric_safe(.data$part_no)),
      band_name = as.character(.data$band_name),
      band_low_hz = as_numeric_safe(.data$band_low_hz),
      band_high_hz = as_numeric_safe(.data$band_high_hz),
      power_db = as_numeric_safe(.data$power_db)
    ) |>
    dplyr::filter(.data$band_name %in% FORMAL_BANDS$band_name) |>
    dplyr::left_join(FORMAL_BANDS, by = "band_name", suffix = c("", "_formal"))

  bad_band <- dat |>
    dplyr::filter(
      !is.finite(.data$band_low_hz) |
        !is.finite(.data$band_high_hz) |
        abs(.data$band_low_hz - .data$band_low_hz_formal) > 1e-9 |
        abs(.data$band_high_hz - .data$band_high_hz_formal) > 1e-9
    )
  if (nrow(bad_band) > 0) {
    stop("Band frequencies do not match the formal five-band definition in: ", power_file, call. = FALSE)
  }

  dat <- dat |>
    dplyr::select(-dplyr::all_of(c("band_low_hz_formal", "band_high_hz_formal")))

  wrong_roi <- unique(dat$roi[!is.na(dat$roi) & dat$roi != roi_expected])
  if (length(wrong_roi) > 0) {
    stop("Unexpected ROI label in ", power_file, ": ", paste(wrong_roi, collapse = ", "), call. = FALSE)
  }
  wrong_process <- unique(dat$process_type[!is.na(dat$process_type) & dat$process_type != process_expected])
  if (length(wrong_process) > 0) {
    stop("Unexpected process_type in ", power_file, ": ", paste(wrong_process, collapse = ", "), call. = FALSE)
  }

  key_cols <- c("sub_id", "trial", "segment_no", "part_no", "band_name")
  assert_unique_keys(dat, key_cols, paste0("Power file for ", roi_expected, " / ", process_expected))

  dat
}

prepare_corridor <- function(roi, encoding) {
  power_file <- power_file_for(roi, "corridor")
  dat <- prepare_power_common(power_file, roi, "corridor")

  dat <- dat |>
    dplyr::left_join(
      encoding,
      by = c("trial", "segment_no", "part_no")
    ) |>
    dplyr::mutate(
      segment_type = "corridor",
      SegmentNumber = 2L * (.data$segment_no - 1L) + .data$part_no,
      n_junctions = n_junctions_from_condition(.data$condition),
      is_last_segment = is.finite(.data$n_junctions) &
        .data$segment_no == (.data$n_junctions + 1L),
      source_power_file = power_file
    )

  dat
}

prepare_junction <- function(roi, encoding) {
  power_file <- power_file_for(roi, "junction")
  dat <- prepare_power_common(power_file, roi, "junction")

  dat <- dat |>
    dplyr::left_join(
      encoding,
      by = c("trial", "segment_no")
    ) |>
    dplyr::mutate(
      Angle2orientation = NA_real_,
      segment_type = "junction",
      SegmentNumber = .data$segment_no,
      n_junctions = n_junctions_from_condition(.data$condition),
      is_last_segment = is.finite(.data$n_junctions) &
        .data$segment_no == .data$n_junctions,
      source_power_file = power_file
    )

  dat
}

make_all_segments <- function(corridor, junction) {
  cdat <- corridor |>
    dplyr::mutate(
      SegmentNumber = 3L * (.data$segment_no - 1L) + .data$part_no
    )

  jdat <- junction |>
    dplyr::mutate(
      SegmentNumber = 3L * .data$segment_no
    )

  # bind_rows keeps union of columns and fills unavailable fields with NA.
  dplyr::bind_rows(cdat, jdat) |>
    dplyr::arrange(
      .data$sub_id,
      .data$trial,
      .data$band_name,
      .data$SegmentNumber,
      .data$segment_type
    )
}

save_prepared <- function(dat, scope, roi) {
  scope_dir <- file.path(output_root, scope)
  dir.create(scope_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(scope_dir, paste0(safe_filename(roi), ".csv"))
  readr::write_csv(dat, out_file, na = "")
  out_file
}

summarize_preparation <- function(dat, scope, roi, out_file) {
  tibble::tibble(
    roi = roi,
    scope = scope,
    output_file = out_file,
    n_rows = nrow(dat),
    n_subjects = dplyr::n_distinct(dat$sub_id),
    n_trials = dplyr::n_distinct(dat$trial),
    n_bands = dplyr::n_distinct(dat$band_name),
    n_finite_power = sum(is.finite(dat$power_db)),
    n_distance_missing = sum(!is.finite(dat$Distance)),
    n_angle_xaxis_missing = sum(!is.finite(dat$Angle2Xaxis)),
    n_angle_orientation_missing = if (scope == "corridor") {
      sum(!is.finite(dat$Angle2orientation))
    } else {
      NA_integer_
    },
    n_last_segment_rows = sum(dat$is_last_segment %in% TRUE, na.rm = TRUE)
  )
}

# =============================================================================
# 4. Read canonical encoding tables once
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("S01 | Prepare Power encoding LMM data\n")
cat(strrep("=", 78), "\n", sep = "")
cat("Power input root : ", power_root, "\n", sep = "")
cat("Prepared output  : ", output_root, "\n", sep = "")
cat("Corridor input   : half_2p5s / ROI average\n")
cat("Junction input   : whole_2p5s / ROI average\n")
cat("Bands            : 0.5-4, 4-12, 12-30, 30-70, 70-150 Hz\n")
cat("Predictor scaling: raw units; NO z-transformation\n")
cat(strrep("=", 78), "\n\n", sep = "")

corridor_encoding <- read_corridor_encoding(corridor_encoding_file)
junction_encoding <- read_junction_encoding(junction_encoding_file)

cat("Corridor encoding rows: ", nrow(corridor_encoding), "\n", sep = "")
cat("Junction encoding rows: ", nrow(junction_encoding), "\n\n", sep = "")

# =============================================================================
# 5. Prepare all ROIs
# =============================================================================

qc_rows <- list()
unmatched_rows <- list()
q_i <- 1L
u_i <- 1L

for (roi in SELECTED_ROIS) {
  cat("[ROI] ", roi, "\n", sep = "")

  corridor <- prepare_corridor(roi, corridor_encoding)
  junction <- prepare_junction(roi, junction_encoding)
  all_segments <- make_all_segments(corridor, junction)

  # Save unmatched encoding keys for audit. Power rows are retained; model-level
  # complete-case filtering happens in S02 separately for each predictor.
  unmatched_corridor <- corridor |>
    dplyr::filter(
      !is.finite(.data$Distance) |
        !is.finite(.data$Angle2Xaxis) |
        !is.finite(.data$Angle2orientation)
    ) |>
    dplyr::distinct(
      .data$roi, .data$trial, .data$condition, .data$segment_no,
      .data$part_no, .data$part_name
    ) |>
    dplyr::mutate(scope = "corridor", .before = 1)

  unmatched_junction <- junction |>
    dplyr::filter(
      !is.finite(.data$Distance) |
        !is.finite(.data$Angle2Xaxis)
    ) |>
    dplyr::distinct(
      .data$roi, .data$trial, .data$condition, .data$segment_no,
      .data$part_no, .data$part_name
    ) |>
    dplyr::mutate(scope = "junction", .before = 1)

  if (nrow(unmatched_corridor) > 0) {
    unmatched_rows[[u_i]] <- unmatched_corridor
    u_i <- u_i + 1L
  }
  if (nrow(unmatched_junction) > 0) {
    unmatched_rows[[u_i]] <- unmatched_junction
    u_i <- u_i + 1L
  }

  out_c <- save_prepared(corridor, "corridor", roi)
  out_j <- save_prepared(junction, "junction", roi)
  out_a <- save_prepared(all_segments, "all_segments", roi)

  qc_rows[[q_i]] <- summarize_preparation(corridor, "corridor", roi, out_c); q_i <- q_i + 1L
  qc_rows[[q_i]] <- summarize_preparation(junction, "junction", roi, out_j); q_i <- q_i + 1L
  qc_rows[[q_i]] <- summarize_preparation(all_segments, "all_segments", roi, out_a); q_i <- q_i + 1L

  cat(
    "  corridor rows=", nrow(corridor),
    " | junction rows=", nrow(junction),
    " | all_segments rows=", nrow(all_segments), "\n",
    sep = ""
  )
}

preparation_qc <- dplyr::bind_rows(qc_rows)
readr::write_csv(preparation_qc, file.path(qc_root, "preparation_qc.csv"), na = "")

if (length(unmatched_rows) > 0) {
  unmatched <- dplyr::bind_rows(unmatched_rows)
} else {
  unmatched <- tibble::tibble(
    scope = character(),
    roi = character(),
    trial = integer(),
    condition = character(),
    segment_no = integer(),
    part_no = integer(),
    part_name = character()
  )
}
readr::write_csv(unmatched, file.path(qc_root, "unmatched_encoding_keys.csv"), na = "")

readr::write_csv(FORMAL_BANDS, file.path(qc_root, "formal_bands.csv"))
readr::write_lines(SELECTED_ROIS, file.path(qc_root, "formal_rois.txt"))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("S01 completed\n")
cat("Prepared data: ", output_root, "\n", sep = "")
cat("QC file      : ", file.path(qc_root, "preparation_qc.csv"), "\n", sep = "")
cat("Unmatched    : ", nrow(unmatched), " unique encoding keys\n", sep = "")
cat(strrep("=", 78), "\n", sep = "")
