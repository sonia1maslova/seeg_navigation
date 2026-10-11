# Shared project configuration loader for the SEEG navigation project.
#
# Cross-language runtime project configuration:
#   config/project_config.json
# Human-editable source:
#   config/project_config.py
#
# Override project location without editing scripts:
#   Sys.setenv(SEEG_NAVIGATION_ROOT = "/path/to/seeg_navigation")

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required to read project_config.json.", call. = FALSE)
}

project_root <- Sys.getenv(
  "SEEG_NAVIGATION_ROOT",
  unset = "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
)
project_root <- normalizePath(project_root, mustWork = FALSE)
config_dir <- file.path(project_root, "config")
project_config_file <- file.path(config_dir, "project_config.json")

if (!file.exists(project_config_file)) {
  stop("Project config not found:\n", project_config_file, call. = FALSE)
}

.project_settings <- jsonlite::fromJSON(
  project_config_file,
  simplifyVector = FALSE
)

resolve_project_path <- function(path_value) {
  path_value <- as.character(path_value)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path_value)) {
    return(normalizePath(path_value, mustWork = FALSE))
  }
  normalizePath(file.path(project_root, path_value), mustWork = FALSE)
}

resolve_path_section <- function(section) {
  out <- lapply(section, resolve_project_path)
  names(out) <- names(section)
  out
}

path_cfg <- .project_settings$paths
file_cfg <- resolve_path_section(.project_settings$files)
external_cfg <- resolve_path_section(.project_settings$external)

cfg <- list(
  root = project_root,
  config_dir = config_dir,
  sampling_rate_hz = as.numeric(.project_settings$sampling_rate_hz),
  dataset_modes = unlist(.project_settings$dataset_modes, use.names = FALSE),
  process_types = unlist(.project_settings$process_types, use.names = FALSE),
  condition_order = unlist(.project_settings$condition_order, use.names = FALSE),
  conditions = .project_settings$conditions,
  segment_windows_sec = .project_settings$segment_windows_sec,
  ap_split_roi_definition = .project_settings$ap_split_roi_definition,
  roi_abbreviations = setNames(
    vapply(.project_settings$roi_abbreviations, function(x) x$abbreviation, character(1)),
    vapply(.project_settings$roi_abbreviations, function(x) x$full_name, character(1))
  ),
  data = list(
    root = resolve_project_path(path_cfg$data),
    processed = resolve_project_path(path_cfg$processed),
    trial_timeseries = resolve_project_path(path_cfg$trial_timeseries),
    trial_timeseries_by_turn_number = resolve_project_path(
      path_cfg$trial_timeseries_by_turn_number
    ),
    derived = resolve_project_path(path_cfg$derived),
    derived_features = resolve_project_path(path_cfg$derived_features),
    behavior = resolve_project_path(path_cfg$behavior),
    encoding = resolve_project_path(path_cfg$encoding),
    metadata = resolve_project_path(path_cfg$metadata),
    feature_selection = resolve_project_path(path_cfg$feature_selection)
  ),
  codes = list(
    root = resolve_project_path(path_cfg$codes),
    matlab = file.path(resolve_project_path(path_cfg$codes), "matlab"),
    python = file.path(resolve_project_path(path_cfg$codes), "python"),
    r = file.path(resolve_project_path(path_cfg$codes), "r")
  ),
  results = list(root = resolve_project_path(path_cfg$results)),
  figures = list(root = resolve_project_path(path_cfg$figures)),
  logs = list(root = resolve_project_path(path_cfg$logs)),
  files = file_cfg,
  external = external_cfg
)

modules <- c(
  "data_preparation", "behavior", "power", "feature_extraction",
  "lasso", "rsa", "plv", "pac", "granger", "amtg_hipp_similarity"
)
for (module_name in modules) {
  cfg$results[[module_name]] <- file.path(cfg$results$root, module_name)
  cfg$figures[[module_name]] <- file.path(cfg$figures$root, module_name)
  cfg$logs[[module_name]] <- file.path(cfg$logs$root, module_name)
}
cfg$modules <- modules

load_analysis_config <- function(filename) {
  if (!grepl("\\.json$", filename)) {
    filename <- paste0(filename, ".json")
  }
  path <- file.path(cfg$config_dir, filename)
  if (!file.exists(path)) {
    stop("Analysis config not found:\n", path, call. = FALSE)
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

roi_abbreviation <- function(roi_name) {
  key <- as.character(roi_name)
  value <- cfg$roi_abbreviations[key]
  ifelse(is.na(value), key, unname(value))
}
