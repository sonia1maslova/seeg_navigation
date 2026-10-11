#!/usr/bin/env Rscript

# ============================================================================
# S05_plot_corridor_bandpass_dfi_group.R
# ----------------------------------------------------------------------------
# Publication-style visualization for corridor band-pass waveform DFI.
#
# This script ONLY plots existing group-level results.
# It does not recompute DFI, permutation tests, or cluster inference.
#
# Final output policy:
#   - PNG only
#   - separate observed and point-wise significance figures
#   - per-figure display overrides can restrict observed frequency panels
#   - observed DFI and point-wise two-sided permutation P are saved separately
#   - no titles or x-axis labels inside result figures; major ticks are retained
#   - one shared standalone delay-axis reference
#   - no PDF
#   - no ROI subfolders
#   - no r_pretty directory
#
# Plot conventions:
#   - red  : aMTG -> Hipp
#   - blue : Hipp -> aMTG
#   - white-centered circles: ordinary delay points
#   - filled diamonds: point-wise two-sided P < 0.05 (uncorrected)
#   - translucent band: cluster-level P < 0.05, corrected across delays
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------
# Packages
# ----------------------------------------------------------------------------
required_packages <- c("ggplot2", "dplyr", "readr", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(ggplot2)
library(dplyr)
library(readr)

# ============================================================================
# USER SETTINGS -- normally only edit this section
# ============================================================================

# Permutation method options:
#   "matched_cross_trial"
#   "within_trial"
#   "within_subject_condition"
permutation_method <- "within_subject_condition"

pointwise_alpha <- 0.05
cluster_alpha <- 0.05

# ============================================================================
# Project configuration
# ============================================================================
project_root <- Sys.getenv("SEEG_NAVIGATION_ROOT")
if (project_root == "") {
  project_root <- "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
}

roi_pairs <- c(
  "aMTG_LAntHipp",
  "aMTG_LPostHipp",
  "aMTG_RAntHipp"
)

features <- c(
  "Distance",
  "Angle2Xaxis_raw",
  "Angle2orientation_raw",
  "SegmentNumber"
)

# Optional scoped redraw; empty selectors preserve the full original workflow.
selected_roi <- Sys.getenv("DFI_PLOT_ROI", "")
selected_feature <- Sys.getenv("DFI_PLOT_FEATURE", "")
plot_kind <- Sys.getenv("DFI_PLOT_KIND", "both")
if (!plot_kind %in% c("both", "observed", "pointwise")) {
  stop("DFI_PLOT_KIND must be both, observed, or pointwise.")
}
if (nzchar(selected_roi)) {
  if (!selected_roi %in% roi_pairs) stop("Unknown DFI_PLOT_ROI: ", selected_roi)
  roi_pairs <- selected_roi
}
if (nzchar(selected_feature)) {
  if (!selected_feature %in% features) stop("Unknown DFI_PLOT_FEATURE: ", selected_feature)
  features <- selected_feature
}
write_shared_axis <- !nzchar(selected_roi) && !nzchar(selected_feature) && plot_kind == "both"

# This requested observed figure shows only the three low-frequency bands.
# The original five-band data/inference checks and y-axis scaling are retained.
observed_display_band_overrides <- list(
  "within_subject_condition/Angle2orientation_raw/aMTG_LAntHipp" = 1:3
)

# Resolve all region labels from the shared project configuration.
project_config <- jsonlite::fromJSON(file.path(project_root, "config", "project_config.json"))
roi_abbreviations <- setNames(
  project_config$roi_abbreviations$abbreviation,
  project_config$roi_abbreviations$full_name
)
required_regions <- c(
  "Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
  "Hippocampus L_Posterior", "Hippocampus R_Anterior"
)
if (any(!required_regions %in% names(roi_abbreviations))) {
  stop("Shared configuration is missing one or more DFI region abbreviations.")
}

# Keep the same colors as the original plotting script.
color_forward <- "#EA415E"
color_reverse <- "#5A85E8"

# Approved compact split-figure settings.
fig_width <- 14.5
fig_height <- 4.4
fig_dpi <- 320

# ----------------------------------------------------------------------------
# Resolve result and figure roots from the permutation method.
# ----------------------------------------------------------------------------
if (permutation_method == "matched_cross_trial") {

  results_bandpass_root <- file.path(
    project_root, "results", "dfi", "bandpass"
  )

} else if (permutation_method == "within_trial") {

  results_bandpass_root <- file.path(
    project_root, "results", "dfi",
    "bandpass_within_trial_permutation"
  )

} else if (permutation_method == "within_subject_condition") {

  results_bandpass_root <- file.path(
    project_root, "results", "dfi",
    "bandpass_within_subject_condition_permutation"
  )

} else {

  stop(
    "Unknown permutation_method: ", permutation_method,
    "\nValid options:",
    "\n  matched_cross_trial",
    "\n  within_trial",
    "\n  within_subject_condition"
  )
}

# Final figure structure:
# figures/dfi/bandpass/corridor/group_level/<permutation_method>/<feature>/
figure_root <- file.path(
  project_root,
  "figures", "dfi", "bandpass", "corridor", "group_level",
  permutation_method
)
figure_root_override <- Sys.getenv("DFI_FIGURE_ROOT", "")
if (nzchar(figure_root_override)) figure_root <- figure_root_override

dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
direction_label <- function(direction, roi) {
  region_a <- unname(roi_abbreviations["Temporal_Mid L_Anterior"])
  partner_code <- sub("^aMTG_", "", roi)
  partner_name <- c(
    LAntHipp = "Hippocampus L_Anterior",
    LPostHipp = "Hippocampus L_Posterior",
    RAntHipp = "Hippocampus R_Anterior"
  )[[partner_code]]
  region_b <- unname(roi_abbreviations[partner_name])
  if (identical(direction, paste0("aMTG_to_", partner_code))) {
    return(paste0(region_a, " → ", region_b))
  }
  if (identical(direction, paste0(partner_code, "_to_aMTG"))) {
    return(paste0(region_b, " → ", region_a))
  }
  stop("Unrecognized direction for ", roi, ": ", direction)
}

superscript_integer <- function(x) {
  glyphs <- c(
    "-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³",
    "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹"
  )
  paste0(glyphs[strsplit(as.character(x), "", fixed = TRUE)[[1]]], collapse = "")
}

integer_breaks <- function(limits) {
  step <- max(1, ceiling(diff(range(limits)) / 6))
  seq(ceiling(min(limits) / step) * step,
      floor(max(limits) / step) * step, by = step)
}

make_symmetric_limits <- function(x, expansion = 1.08) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-1, 1))
  m <- max(abs(x), na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- 1
  c(-m * expansion, m * expansion)
}

base_theme <- function() {
  theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(
        size = 21, face = "bold", hjust = 0.5,
        margin = margin(b = 7)
      ),
      plot.subtitle = element_text(
        size = 13, color = "grey30", hjust = 0.5,
        margin = margin(b = 8)
      ),
      strip.background = element_blank(),
      strip.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.6, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      panel.grid.major.y = element_line(linewidth = 0.35, color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = grid::unit(0.8, "lines"),
      legend.position = "top",
      legend.justification = "center",
      legend.title = element_blank(),
      legend.text = element_text(size = 15.6),
      legend.key.width = grid::unit(1.4, "lines"),
      plot.margin = margin(10, 10, 8, 10)
    )
}

add_common_scales <- function(p, direction_colors) {
  p +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    scale_fill_manual(values = direction_colors, drop = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 1200, by = 200),
      expand = expansion(mult = c(0.01, 0.01))
    )
}

add_cluster_shading <- function(p, cluster_df) {
  if (is.null(cluster_df) || nrow(cluster_df) == 0) return(p)

  p + geom_rect(
    data = cluster_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = direction_display
    ),
    inherit.aes = FALSE,
    alpha = 0.11,
    color = NA,
    show.legend = FALSE
  )
}

add_pointwise_markers <- function(p, df, y_var) {
  sig_df <- df %>%
    filter(sig_raw_05 %in% TRUE, is.finite(.data[[y_var]]))

  if (nrow(sig_df) == 0) return(p)

  p + geom_point(
    data = sig_df,
    aes(
      x = delay_ms,
      y = .data[[y_var]],
      fill = direction_display
    ),
    shape = 23,
    size = 2.0,
    stroke = 0.55,
    color = "black",
    show.legend = FALSE
  )
}

# ----------------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------------
message("\n============================================================")
message("Corridor band-pass DFI plotting")
message("Permutation method : ", permutation_method)
message("Results root       : ", results_bandpass_root)
message("Figure root        : ", figure_root)
message("Output             : separate observed/significance PNGs and shared delay axis")
message("============================================================")

generated_files <- character()
expected_delay_grid <- as.numeric(seq(10, 1200, by = 10))

for (feature_name in features) {

  feature_out <- file.path(figure_root, feature_name)
  dir.create(feature_out, recursive = TRUE, showWarnings = FALSE)

  for (roi in roi_pairs) {

    message(
      "\n------------------------------------------------------------\n",
      "Feature: ", feature_name, " | ROI: ", roi
    )

    group_root <- file.path(
      results_bandpass_root,
      "corridor", roi,
      "significance", "group_level"
    )

    group_csv <- file.path(
      group_root,
      "group_level_significance.csv"
    )

    cluster_csv <- file.path(
      group_root,
      "clusters",
      "significant_clusters.csv"
    )

    if (!file.exists(group_csv)) {
      stop("Missing group CSV: ", group_csv)
    }

    df <- read_csv(group_csv, show_col_types = FALSE)

    required_cols <- c(
      "roi_pair", "band_index", "band_name",
      "band_low_hz", "band_high_hz",
      "feature", "direction", "delay_ms", "n_subjects",
      "DFI_group_mean", "null_mean",
      "p_two_sided", "sig_raw_05"
    )

    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(
        "Missing columns in ", group_csv, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }

    # ------------------------------------------------------------------------
    # Select current feature.
    # ------------------------------------------------------------------------
    df <- df %>%
      filter(feature == feature_name) %>%
      mutate(
        band_label_plot = sprintf(
          "%g-%g Hz", band_low_hz, band_high_hz
        ),
        direction_raw = direction,
        direction_display = vapply(
          direction,
          direction_label,
          character(1),
          roi = roi
        ),
        neglog10_p = -log10(
          pmax(p_two_sided, .Machine$double.xmin)
        ),
        sig_raw_05 = as.logical(sig_raw_05)
      )

    if (nrow(df) == 0) {
      stop(
        "No rows found for feature ", feature_name,
        " | ROI: ", roi
      )
    }

    for (delay_values in split(df$delay_ms, interaction(df$band_index, df$direction))) {
      if (!identical(sort(as.numeric(delay_values)), expected_delay_grid)) {
        stop("Incomplete or duplicated delay grid for ", feature_name, " | ", roi)
      }
    }
    if (length(unique(df$band_index)) != 5L || length(unique(df$direction)) != 2L) {
      stop("Expected five bands and both directions for ", feature_name, " | ", roi)
    }

    # Display scaling only: the original DFI and inference columns are untouched.
    observed_max <- max(abs(c(df$DFI_group_mean, df$null_mean)), na.rm = TRUE)
    if (!is.finite(observed_max) || observed_max <= 0) {
      stop("No finite nonzero DFI values for ", feature_name, " | ", roi)
    }
    observed_exponent <- floor(log10(observed_max))
    observed_scale <- 10^observed_exponent
    df$DFI_group_mean_plot <- df$DFI_group_mean / observed_scale
    df$null_mean_plot <- df$null_mean / observed_scale

    # Stable band order.
    band_levels <- df %>%
      distinct(
        band_index, band_low_hz,
        band_high_hz, band_label_plot
      ) %>%
      arrange(
        band_index, band_low_hz, band_high_hz
      ) %>%
      pull(band_label_plot)

    df <- df %>%
      mutate(
        band_label_plot = factor(
          band_label_plot,
          levels = band_levels
        )
      )

    # ------------------------------------------------------------------------
    # Direction colors.
    # ------------------------------------------------------------------------
    direction_levels <- unique(df$direction_display)

    direction_colors <- setNames(
      rep("#555555", length(direction_levels)),
      direction_levels
    )

    for (level in direction_levels) {
      raw_direction <- unique(df$direction_raw[
        df$direction_display == level
      ])[1]
      if (grepl("^aMTG_to_", raw_direction)) {
        direction_colors[level] <- color_forward
      } else {
        direction_colors[level] <- color_reverse
      }
    }

    if (
      length(direction_levels) >= 1 &&
      all(direction_colors == "#555555")
    ) {
      direction_colors[1] <- color_forward
      if (length(direction_levels) >= 2) {
        direction_colors[2] <- color_reverse
      }
    }

    # ------------------------------------------------------------------------
    # Cluster intervals.
    # ------------------------------------------------------------------------
    cluster_df <- NULL

    if (file.exists(cluster_csv)) {

      cluster_df <- read_csv(
        cluster_csv,
        show_col_types = FALSE
      )

      if (nrow(cluster_df) > 0) {

        cluster_df <- cluster_df %>%
          filter(
            feature == feature_name,
            as.logical(sig_cluster_05),
            p_cluster < cluster_alpha
          ) %>%
          mutate(
            direction_display = vapply(
              direction,
              direction_label,
              character(1),
              roi = roi
            )
          )

        if (nrow(cluster_df) > 0) {

          band_lookup <- df %>%
            distinct(
              band_name,
              band_label_plot
            )

          cluster_df <- cluster_df %>%
            left_join(
              band_lookup,
              by = "band_name"
            ) %>%
            mutate(
              # The delay grid is 10 ms, so extend by half a sample.
              xmin = start_delay_ms - 5,
              xmax = end_delay_ms + 5
            )
        }
      }

    } else {

      warning(
        "Cluster CSV not found. Figure will be made ",
        "without cluster shading: ",
        cluster_csv
      )
    }

    common_mapping <- aes(
      x = delay_ms,
      color = direction_display,
      group = direction_display
    )

    # ========================================================================
    # Top row: raw observed DFI
    # ========================================================================
    y_observed <- make_symmetric_limits(
      c(df$DFI_group_mean_plot, df$null_mean_plot)
    )

    observed_band_indices <- observed_display_band_overrides[[
      paste(permutation_method, feature_name, roi, sep = "/")
    ]]
    if (is.null(observed_band_indices)) observed_band_indices <- unique(df$band_index)
    if (!all(observed_band_indices %in% df$band_index)) {
      stop("Requested observed display bands are absent: ", feature_name, " | ", roi)
    }
    df_observed <- df %>%
      filter(band_index %in% observed_band_indices) %>%
      mutate(band_label_plot = droplevels(band_label_plot))
    cluster_df_observed <- cluster_df
    if (!is.null(cluster_df_observed) && nrow(cluster_df_observed) > 0) {
      cluster_df_observed <- cluster_df_observed %>%
        filter(band_name %in% unique(df_observed$band_name)) %>%
        mutate(band_label_plot = factor(
          as.character(band_label_plot), levels = levels(df_observed$band_label_plot)
        ))
    }
    observed_width <- fig_width * length(unique(df_observed$band_index)) / length(unique(df$band_index))
    # Leave sufficient room for the unchanged two-direction legend.
    if (length(unique(df_observed$band_index)) < length(unique(df$band_index))) {
      observed_width <- max(observed_width, 10)
    }
    message("Observed display bands: ", paste(levels(df_observed$band_label_plot), collapse = ", "))

    p_observed <- ggplot(df_observed, common_mapping) +
      geom_hline(
        yintercept = 0,
        linewidth = 0.65,
        color = "grey55"
      )

    p_observed <- add_cluster_shading(
      p_observed,
      cluster_df_observed
    )

    p_observed <- p_observed +
      geom_line(
        aes(y = null_mean_plot),
        linewidth = 0.75,
        linetype = "22",
        alpha = 0.40,
        show.legend = FALSE
      ) +
      geom_line(
        aes(y = DFI_group_mean_plot),
        linewidth = 1.05,
        lineend = "round"
      ) +
      geom_point(
        aes(y = DFI_group_mean_plot),
        shape = 21,
        size = 1.7,
        stroke = 0.7,
        fill = "white"
      )

    p_observed <- add_pointwise_markers(
      p_observed,
      df_observed,
      "DFI_group_mean_plot"
    )

    p_observed <- p_observed +
      facet_wrap(
        ~band_label_plot,
        nrow = 1,
        drop = FALSE
      ) +
      coord_cartesian(
        ylim = y_observed,
        clip = "off"
      ) +
      scale_y_continuous(
        breaks = integer_breaks,
        labels = function(x) format(x, trim = TRUE, scientific = FALSE)
      ) +
      labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Observed DFI"
      ) +
      base_theme() +
      theme(
        legend.position = "top"
      )

    p_observed <- add_common_scales(
      p_observed,
      direction_colors
    )

    scale_note <- data.frame(
      band_label_plot = factor(
        levels(df_observed$band_label_plot)[1],
        levels = levels(df_observed$band_label_plot)
      ),
      x = 0,
      y = Inf,
      label = paste0("×10", superscript_integer(observed_exponent))
    )
    p_observed <- p_observed +
      geom_text(
        data = scale_note,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = -0.35,
        size = 4.5,
        fontface = "bold"
      )

    # ========================================================================
    # Bottom row: -log10(point-wise P)
    # ========================================================================
    p_threshold <- -log10(pointwise_alpha)

    max_p <- max(
      df$neglog10_p,
      na.rm = TRUE
    )

    if (!is.finite(max_p)) {
      max_p <- p_threshold
    }

    y_p <- c(
      0,
      max(
        max_p * 1.08,
        p_threshold + 0.35
      )
    )

    p_pvalue <- ggplot(df, common_mapping) +
      geom_hline(
        yintercept = p_threshold,
        linetype = "dashed",
        linewidth = 0.75,
        color = "grey35"
      )

    p_pvalue <- add_cluster_shading(
      p_pvalue,
      cluster_df
    )

    p_pvalue <- p_pvalue +
      geom_line(
        aes(y = neglog10_p),
        linewidth = 1.05,
        lineend = "round"
      ) +
      geom_point(
        aes(y = neglog10_p),
        shape = 21,
        size = 1.7,
        stroke = 0.7,
        fill = "white"
      )

    p_pvalue <- add_pointwise_markers(
      p_pvalue,
      df,
      "neglog10_p"
    )

    p_pvalue <- p_pvalue +
      facet_wrap(
        ~band_label_plot,
        nrow = 1,
        drop = FALSE
      ) +
      coord_cartesian(
        ylim = y_p
      ) +
      scale_y_continuous(
        breaks = function(x) seq(0, ceiling(max(x)), by = 1),
        labels = function(x) format(round(x), trim = TRUE, scientific = FALSE)
      ) +
      labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = expression(-log[10](P))
      ) +
      base_theme() +
      theme(
        legend.position = "top"
      )

    p_pvalue <- add_common_scales(
      p_pvalue,
      direction_colors
    )

    # ------------------------------------------------------------------------
    # Save the two rows as separate final figures.
    # ------------------------------------------------------------------------
    observed_file <- file.path(
      feature_out,
      paste0(roi, "_observed.png")
    )
    pvalue_file <- file.path(
      feature_out,
      paste0(roi, "_pointwise_significance.png")
    )

    if (plot_kind %in% c("both", "observed")) {
      # Guard against excluded bands reappearing through a layer or facet level.
      observed_layout <- ggplot_build(p_observed)$layout$layout
      if (!identical(as.character(observed_layout$band_label_plot), levels(df_observed$band_label_plot))) {
        stop("Unexpected observed frequency panels: ", feature_name, " | ", roi)
      }
      ggsave(
      filename = observed_file,
      plot = p_observed,
      width = observed_width,
      height = fig_height,
      units = "in",
      dpi = fig_dpi,
      bg = "white"
      )
      message("Saved: ", observed_file)
      generated_files <- c(generated_files, observed_file)
    }
    if (plot_kind %in% c("both", "pointwise")) {
      ggsave(
      filename = pvalue_file,
      plot = p_pvalue,
      width = fig_width,
      height = fig_height,
      units = "in",
      dpi = fig_dpi,
      bg = "white"
      )
      message("Saved: ", pvalue_file)
      generated_files <- c(generated_files, pvalue_file)
    }
  }
}

expected_per_pair <- if (plot_kind == "both") 2L else 1L
if (length(generated_files) != length(features) * length(roi_pairs) * expected_per_pair) {
  stop("Not all requested figures were generated.")
}

# Standalone reference: the same delay scale as one frequency panel.
# Only short major ticks and labels every 200 ms are displayed.
if (write_shared_axis) {
axis_file <- file.path(figure_root, "delay_axis_reference.png")
png(axis_file, width = 3.4, height = 0.9, units = "in", res = fig_dpi,
    type = "cairo", bg = "white")
grid::grid.newpage()
axis_left <- 0.045
axis_right <- 0.93
axis_x <- function(ms) axis_left + (axis_right - axis_left) * ms / 1200
grid::grid.lines(
  x = grid::unit(c(axis_left, axis_right), "npc"),
  y = grid::unit(c(0.83, 0.83), "npc"),
  gp = grid::gpar(lwd = 1.2, col = "black")
)
delay_ticks <- seq(0, 1200, by = 200)
grid::grid.segments(
  x0 = axis_x(delay_ticks), x1 = axis_x(delay_ticks),
  y0 = 0.83, y1 = 0.76,
  gp = grid::gpar(lwd = 0.8, col = "black")
)
major_ticks <- seq(0, 1200, by = 200)
grid::grid.text(
  label = major_ticks, x = axis_x(major_ticks), y = 0.64,
  gp = grid::gpar(fontsize = 13, col = "black")
)
grid::grid.text(
  label = "Delay (ms)", x = (axis_left + axis_right) / 2, y = 0.39,
  gp = grid::gpar(fontsize = 15, fontface = "bold", col = "black")
)
dev.off()
message("Saved: ", axis_file)
generated_files <- c(generated_files, axis_file)
}
if (any(!file.exists(generated_files)) || any(file.info(generated_files)$size <= 0)) {
  stop("One or more generated figure files are missing or empty.")
}
message("Verified generated file count: ", length(generated_files))

message("\n============================================================")
message("R plotting completed.")
message("Permutation method : ", permutation_method)
message("Output root        : ", figure_root)
message("============================================================")
