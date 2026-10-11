# PAC segment-position effects, using the approved power-encoding heatmap style.
# Statistics come from formal_raw_pac. Stored q-values are verified, not replaced.
# Hippocampal effects form 18-test BH families; aMTG forms 6-test BH families.
# The two terminal-corridor variants remain separate statistical families.

required_packages <- c("dplyr", "readr", "ggplot2", "grid", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
library(dplyr)
library(readr)
library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(arguments)) arguments[[1]] else Sys.getenv(
  "SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
)
project_root <- normalizePath(project_root, mustWork = TRUE)
Sys.setenv(SEEG_NAVIGATION_ROOT = project_root)
source(file.path(project_root, "codes", "r", "project_config.R"))

result_root <- file.path(project_root, "results", "pac", "by_lina", "formal_raw_pac")
output_dir <- Sys.getenv("PAC_SEGMENT_HEATMAP_OUTPUT", unset = file.path(
  project_root, "figures", "pac", "by_lina", "pac_segment_effect", "heatmaps"
))
hpc_input <- file.path(result_root, "02_hpc_path_accumulation", "model_coefficients.csv")
amtg_input <- file.path(result_root, "04_amtg_followup_and_behavior", "model_coefficients.csv")
variants <- c(
  exclude_last = "all_dynamic_split",
  include_last = "all_dynamic_split_with_terminal_C"
)
selected_rois <- c(
  "Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
  "Hippocampus L_Posterior", "Hippocampus R_Anterior"
)
# Reuse the exact short labels in the user-selected power heatmap. The current
# generic configuration has longer labels, so do not alter it project-wide.
roi_display_map <- setNames(
  c("aMTG", "Hipp L Ant", "Hipp L Post", "Hipp R Ant"), selected_rois
)
phase_levels <- c("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
amp_levels <- c("low_gamma_30_70Hz", "high_gamma_70_150Hz")
combo_labels <- c(
  "0.5–4 → 30–70", "4–12 → 30–70", "12–30 → 30–70",
  "0.5–4 → 70–150", "4–12 → 70–150", "12–30 → 70–150"
)
alpha_raw <- 0.05
alpha_fdr <- 0.05
beta_display_unit <- 1e-5
q_tolerance <- 1e-12

read_coefficients <- function(path, analysis_family_name) {
  dat <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  required_columns <- c(
    "analysis_family", "scope", "sample_mode", "roi", "phase_band", "amp_band",
    "predictor", "target_term", "beta", "p_value", "q_value", "fit_status",
    "fdr_family_id", "m_planned", "m_valid", "model_formula"
  )
  absent <- setdiff(required_columns, names(dat))
  if (length(absent)) stop("Missing columns in ", path, ": ", paste(absent, collapse = ", "))
  dat |>
    filter(
      analysis_family == analysis_family_name,
      scope %in% unname(variants), sample_mode == "single_roi",
      predictor == "segment_position", target_term == "segment_position",
      roi %in% selected_rois, phase_band %in% phase_levels, amp_band %in% amp_levels
    )
}

statistics <- bind_rows(
  read_coefficients(hpc_input, "02A_basic_segment"),
  read_coefficients(amtg_input, "04B_amtg_accumulation")
)
if (nrow(statistics) != 48L) stop("Expected 48 segment effects: 24 per terminal variant.")
if (any(statistics$fit_status != "ok") ||
    any(!is.finite(statistics$beta)) || any(!is.finite(statistics$p_value)) ||
    any(!is.finite(statistics$q_value))) {
  stop("The requested main-effect grid contains failed fits or missing statistics.")
}
if (any(statistics$p_value < 0 | statistics$p_value > 1) ||
    any(statistics$q_value < 0 | statistics$q_value > 1)) stop("Invalid p or q values.")

# Independent audit: use each original family, never the combined display grid.
families <- split(statistics, statistics$fdr_family_id)
if (length(families) != 4L) stop("Expected two hippocampal and two aMTG BH families.")
fdr_audit <- lapply(families, function(one) {
  is_amtg <- all(one$roi == selected_rois[[1]])
  expected_n <- if (is_amtg) 6L else 18L
  expected_rois <- if (is_amtg) selected_rois[[1]] else selected_rois[-1]
  if (nrow(one) != expected_n || !setequal(unique(one$roi), expected_rois) ||
      length(unique(one$scope)) != 1L || any(one$m_planned != expected_n) ||
      any(one$m_valid != expected_n)) stop("Unexpected FDR family: ", one$fdr_family_id[[1]])
  difference <- max(abs(one$q_value - p.adjust(one$p_value, method = "BH", n = expected_n)))
  if (difference > q_tolerance) stop("Stored q-values fail the BH audit.")
  list(
    scope = one$scope[[1]], planned = expected_n, valid = nrow(one),
    max_abs_q_difference = difference,
    raw_p_lt_005 = sum(one$p_value < alpha_raw), q_lt_005 = sum(one$q_value < alpha_fdr)
  )
})
for (family in names(fdr_audit)) {
  one <- fdr_audit[[family]]
  message("Verified BH: ", family, "; n=", one$valid,
          "; max q difference=", format(one$max_abs_q_difference, scientific = TRUE))
}

duplicates <- statistics |> count(scope, roi, phase_band, amp_band) |> filter(n != 1L)
if (nrow(duplicates)) stop("Duplicated effect cells.")
statistics <- statistics |>
  mutate(
    roi_display = unname(roi_display_map[roi]),
    combo_index = match(phase_band, phase_levels) + 3L * (match(amp_band, amp_levels) - 1L),
    x_pos = 1.00 + 0.22 * (combo_index - 1L),
    y_pos = 2.08 - 0.22 * (match(roi, selected_rois) - 1L),
    beta_display = beta / beta_display_unit,
    raw_sig = p_value < alpha_raw, fdr_sig = q_value < alpha_fdr,
    sig_symbol = case_when(fdr_sig ~ "†", raw_sig ~ "*", TRUE ~ "")
  )
if (any(is.na(statistics$roi_display))) stop("Missing reference ROI abbreviation.")

# Use one symmetric scale for both variants; the source beta values stay raw.
fill_limit <- max(abs(statistics$beta_display))
legend_limit <- signif(fill_limit, digits = 2)
if (legend_limit < fill_limit) legend_limit <- signif(fill_limit * 1.05, digits = 2)
if (!is.finite(legend_limit) || legend_limit <= 0) legend_limit <- 1
heatmap_center_x <- mean(range(statistics$x_pos))
roi_label_df <- data.frame(
  label = unname(roi_display_map[selected_rois]), x = 0.87,
  y = 2.08 - 0.22 * (seq_along(selected_rois) - 1L)
)
combo_label_df <- data.frame(
  label = combo_labels, x = 1.00 + 0.22 * (seq_along(combo_labels) - 1L), y = 1.18
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
written_files <- character()

for (variant in names(variants)) {
  dat <- statistics |> filter(scope == unname(variants[[variant]]))
  if (nrow(dat) != 24L) stop("Incomplete effect grid for ", variant)
  plot <- ggplot(dat, aes(x = x_pos, y = y_pos, fill = beta_display)) +
    geom_tile(width = 0.22, height = 0.22, color = "white", linewidth = 0.85) +
    geom_text(aes(label = sig_symbol), size = 5.0, fontface = "bold", vjust = 0.5) +
    geom_text(data = roi_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 4.2) +
    geom_text(data = combo_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1,
              fontface = "bold", size = 3.4) +
    annotate("text", x = heatmap_center_x, y = 2.30,
             label = "Segment number", fontface = "bold", size = 6.2) +
    annotate("text", x = heatmap_center_x, y = 0.72,
             label = "PAC frequency combination (phase → amplitude, Hz)",
             fontface = "bold", size = 4.7) +
    annotate("text", x = heatmap_center_x, y = 0.53,
             label = "* = raw p < 0.05      † = BH-FDR q < 0.05",
             fontface = "bold", size = 3.4, color = "grey20") +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
      limits = c(-legend_limit, legend_limit), name = expression(beta / 10^{-5}),
      guide = guide_colorbar(direction = "vertical", title.position = "top",
                             title.hjust = 0.5, barwidth = grid::unit(0.55, "cm"),
                             barheight = grid::unit(5.1, "cm"), ticks = TRUE,
                             frame.colour = "grey55")
    ) +
    scale_x_continuous(breaks = NULL, limits = c(-0.35, 3.25),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = NULL, limits = c(0.32, 2.82),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = NULL, x = NULL, y = NULL, caption = NULL) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_blank(), plot.caption = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), legend.position = c(0.77, 0.56),
      legend.justification = c(0, 0.5), legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_blank(), legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
  png_file <- file.path(output_dir, paste0(
    "heatmap_Segment_number_all_segments_", variant, ".png"
  ))
  ggsave(png_file, plot, width = 8.0, height = 5.5, units = "in", dpi = 300, bg = "white")
  written_files <- c(written_files, png_file)
  message("Saved ", png_file)
}

manifest <- list(
  analysis = "PAC main segment-position effects", model = "raw_PAC ~ segment_position + (1|sub_id)",
  inputs = c(hpc_input, amtg_input), input_md5 = as.list(tools::md5sum(c(hpc_input, amtg_input))),
  reference_figure = file.path(project_root, "figures", "power", "encoding_lmm", "heatmaps",
                              "heatmap_Segment_number_all_segments_exclude_last.png"),
  roi_order = selected_rois, roi_labels = as.list(roi_display_map),
  scope_mapping = as.list(variants), beta_display_unit = beta_display_unit,
  color_scale_raw_beta = c(-legend_limit, legend_limit) * beta_display_unit,
  significance_rule = "dagger: stored BH q < 0.05; otherwise star: raw p < 0.05; otherwise blank",
  fdr_method = "Stored family-specific q-values; BH independently checked, never replaced",
  fdr_audit = fdr_audit, outputs = written_files
)
jsonlite::write_json(manifest, file.path(output_dir, "plot_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Done. Model coefficients and original FDR results were not modified.")
