# PAC segment-position effects after power or aperiodic-exponent adjustment.
# Reuse the approved power-encoding heatmap style; never refit models or
# replace the family-specific q-values in formal_raw_pac.

required_packages <- c("dplyr", "readr", "ggplot2", "grid", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}
library(dplyr)
library(readr)
library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(arguments)) arguments[[1]] else Sys.getenv(
  "SEEG_NAVIGATION_ROOT", "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
)
project_root <- normalizePath(project_root, mustWork = TRUE)
output_root <- Sys.getenv("PAC_CONTROL_HEATMAP_OUTPUT", unset = file.path(
  project_root, "figures", "pac", "by_lina", "pac_segment_effect", "lmm_effects", "lina_original"
))
input_paths <- c(
  power_control = file.path(project_root, "results", "pac", "by_lina", "formal_raw_pac",
                            "03a_hpc_power_control", "model_coefficients.csv"),
  aperiodic_control_1_150Hz = file.path(project_root, "results", "pac", "by_lina", "formal_raw_pac",
                                     "03b_hpc_exponent_control", "model_coefficients.csv")
)
input_md5_before <- tools::md5sum(input_paths)
variants <- c(
  exclude_last = "all_dynamic_split",
  include_last = "all_dynamic_split_with_terminal_C"
)
rois <- c("Hippocampus L_Anterior", "Hippocampus L_Posterior", "Hippocampus R_Anterior")
# Exact labels used in the user-selected power reference and main PAC heatmap.
# Do not change the project's longer generic labels for unrelated analyses.
roi_labels <- setNames(c("Hipp L Ant", "Hipp L Post", "Hipp R Ant"), rois)
phase_levels <- c("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
amp_levels <- c("low_gamma_30_70Hz", "high_gamma_70_150Hz")
combo_labels <- c(
  "0.5–4 → 30–70", "4–12 → 30–70", "12–30 → 30–70",
  "0.5–4 → 70–150", "4–12 → 70–150", "12–30 → 70–150"
)
beta_display_unit <- 1e-5
q_tolerance <- 1e-12
alpha <- 0.05

read_adjusted <- function(path, family_name, planned) {
  dat <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  required <- c("analysis_family", "scope", "sample_mode", "roi", "phase_band",
                "amp_band", "predictor", "target_term", "beta", "p_value", "q_value",
                "fit_status", "fdr_family_id", "m_planned", "m_valid", "model_formula")
  absent <- setdiff(required, names(dat))
  if (length(absent)) stop("Missing fields: ", paste(absent, collapse = ", "))
  dat <- dat |> filter(
    analysis_family == family_name, sample_mode == "complete_case_adjusted",
    predictor == "segment_position", target_term == "segment_position",
    scope %in% unname(variants), roi %in% rois,
    phase_band %in% phase_levels, amp_band %in% amp_levels
  )
  if (nrow(dat) != 2L * planned) stop("Unexpected adjusted effect count: ", path)
  if (any(dat$fit_status != "ok") || any(!is.finite(dat$beta)) ||
      any(!is.finite(dat$p_value)) || any(!is.finite(dat$q_value))) {
    stop("Requested effects contain failed fits or unavailable statistics.")
  }
  if (any(dat$p_value < 0 | dat$p_value > 1) ||
      any(dat$q_value < 0 | dat$q_value > 1)) stop("Invalid p or q values.")
  groups <- split(dat, dat$fdr_family_id)
  if (length(groups) != 2L) stop("Expected one FDR family per terminal variant.")
  audit <- lapply(groups, function(one) {
    if (nrow(one) != planned || length(unique(one$scope)) != 1L ||
        any(one$m_planned != planned) || any(one$m_valid != planned)) {
      stop("Unexpected original FDR family: ", one$fdr_family_id[[1]])
    }
    difference <- max(abs(one$q_value - p.adjust(one$p_value, method = "BH", n = planned)))
    if (difference > q_tolerance) stop("Stored q-values fail independent BH verification.")
    list(scope = one$scope[[1]], planned = planned, valid = nrow(one),
         max_abs_q_difference = difference,
         raw_p_lt_005 = sum(one$p_value < alpha), q_lt_005 = sum(one$q_value < alpha))
  })
  list(statistics = dat, audit = audit)
}

power <- read_adjusted(input_paths[["power_control"]], "03A2_power_adjusted_PAC", 18L)
aperiodic <- read_adjusted(input_paths[["aperiodic_control_1_150Hz"]],
                           "03B2_exponent_adjusted_PAC", 36L)
if (!"exponent_range" %in% names(aperiodic$statistics) ||
    !setequal(unique(aperiodic$statistics$exponent_range), c("1-30 Hz", "1-150 Hz"))) {
  stop("Cannot verify the full two-range exponent FDR family.")
}
statistics <- bind_rows(
  power$statistics |> mutate(control = "power_control"),
  aperiodic$statistics |> filter(exponent_range == "1-150 Hz") |>
    mutate(control = "aperiodic_control_1_150Hz")
)
duplicates <- statistics |> count(control, scope, roi, phase_band, amp_band) |> filter(n != 1L)
if (nrow(duplicates) || nrow(statistics) != 72L) stop("Incomplete or duplicated plotting grid.")
statistics <- statistics |> mutate(
  roi_display = unname(roi_labels[roi]),
  combo_index = match(phase_band, phase_levels) + 3L * (match(amp_band, amp_levels) - 1L),
  frequency_label = combo_labels[combo_index],
  x_pos = 1.00 + 0.22 * (combo_index - 1L),
  y_pos = 2.08 - 0.22 * (match(roi, rois) - 1L),
  beta_display = beta / beta_display_unit,
  raw_sig = p_value < alpha, fdr_sig = q_value < alpha,
  sig_symbol = case_when(fdr_sig ~ "†", raw_sig ~ "*", TRUE ~ "")
)

# Shared display scale across both controls and both terminal variants.
# Scaling changes only labels/colors, not the source regression coefficients.
fill_limit <- max(abs(statistics$beta_display))
legend_limit <- signif(fill_limit, digits = 2)
if (legend_limit < fill_limit) legend_limit <- signif(fill_limit * 1.05, digits = 2)
if (!is.finite(legend_limit) || legend_limit <= 0) legend_limit <- 1
heatmap_center_x <- mean(range(statistics$x_pos))
tile_bottom <- min(statistics$y_pos) - 0.11
frequency_label_gap <- 0.035
frequency_label_y <- tile_bottom - frequency_label_gap
roi_label_df <- data.frame(label = unname(roi_labels[rois]), x = 0.87,
                           y = 2.08 - 0.22 * (seq_along(rois) - 1L))
combo_label_df <- data.frame(label = combo_labels,
                             x = 1.00 + 0.22 * (seq_along(combo_labels) - 1L),
                             y = frequency_label_y)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
written_files <- character()
summaries <- list()

for (control in names(input_paths)) {
  for (variant in names(variants)) {
    output_dir <- file.path(output_root, variant, control, "heatmaps")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dat <- statistics |> filter(.data$control == .env$control,
                                scope == unname(variants[[variant]]))
    if (nrow(dat) != 18L) stop("Expected 3 ROI by 6 PAC combinations.")
    plot <- ggplot(dat, aes(x = x_pos, y = y_pos, fill = beta_display)) +
      geom_tile(width = 0.22, height = 0.22, color = "white", linewidth = 0.85) +
      geom_text(aes(label = sig_symbol), size = 5.0, fontface = "bold", vjust = 0.5) +
      geom_text(data = roi_label_df, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 4.2) +
      geom_text(data = combo_label_df, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1,
                fontface = "bold", size = 3.4) +
      annotate("text", x = heatmap_center_x, y = 2.30, label = "Segment number",
               fontface = "bold", size = 6.2) +
      annotate("text", x = heatmap_center_x, y = 0.90,
               label = "PAC frequency combination (phase → amplitude, Hz)",
               fontface = "bold", size = 4.7) +
      annotate("text", x = heatmap_center_x, y = 0.68,
               label = "* = raw p < 0.05      † = BH-FDR q < 0.05",
               fontface = "bold", size = 3.4, color = "grey20") +
      scale_fill_gradient2(
        low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
        limits = c(-legend_limit, legend_limit), name = expression(beta / 10^{-5}),
        guide = guide_colorbar(direction = "vertical", title.position = "top",
                               title.hjust = 0.5, barwidth = grid::unit(0.55, "cm"),
                               barheight = grid::unit(3.0, "cm"), ticks = TRUE,
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
        panel.grid = element_blank(), legend.position = c(0.77, 0.615),
        legend.justification = c(0, 0.5), legend.direction = "vertical",
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11, color = "black"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_blank(), legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
      )
    file_base <- paste0("heatmap_Segment_number_all_segments_", variant)
    png_file <- file.path(output_dir, paste0(file_base, ".png"))
    values_file <- file.path(output_dir, paste0(file_base, "_values.csv"))
    ggsave(png_file, plot, width = 8.0, height = 5.5, units = "in", dpi = 300, bg = "white")
    readr::write_csv(dat, values_file, na = "")
    written_files <- c(written_files, png_file, values_file)
    summaries[[paste(control, variant, sep = "__")]] <- list(
      effects = nrow(dat), raw_significant = sum(dat$raw_sig),
      fdr_significant = sum(dat$fdr_sig),
      singular_fits = sum(as.character(dat$is_singular) %in% c("1", "TRUE", "true")),
      stored_family_size = unique(dat$m_planned)
    )
    message("Saved ", png_file, "; raw significant=", sum(dat$raw_sig),
            "; FDR significant=", sum(dat$fdr_sig))
  }
}
if (!identical(input_md5_before, tools::md5sum(input_paths))) {
  stop("The source coefficient files changed during plotting.")
}
manifest <- list(
  analysis = "Power-adjusted and aperiodic-exponent-adjusted PAC segment effects",
  models = list(power_control = "raw_PAC ~ segment_position + LF_power_db + HF_power_db + (1|sub_id)",
                aperiodic_control_1_150Hz = "raw_PAC ~ segment_position + exponent_raw + (1|sub_id)"),
  inputs = as.list(input_paths), input_md5 = as.list(input_md5_before),
  reference = file.path(project_root, "figures", "power", "encoding_lmm", "heatmaps",
                       "heatmap_Segment_number_all_segments_exclude_last.png"),
  roi_order = rois, roi_labels = as.list(roi_labels), scope_mapping = as.list(variants),
  exponent_range_shown = "1-150 Hz", exponent_ranges_in_original_fdr = c("1-30 Hz", "1-150 Hz"),
  beta_display_unit = beta_display_unit,
  common_color_scale_raw_beta = c(-legend_limit, legend_limit) * beta_display_unit,
  frequency_label_gap = frequency_label_gap,
  significance_rule = "dagger: stored BH q < 0.05; otherwise star: raw p < 0.05; otherwise blank",
  original_fdr_preserved = TRUE,
  fdr_audit = list(power_control = power$audit, aperiodic_control = aperiodic$audit),
  summary = summaries, outputs = written_files
)
provenance_dir <- file.path(output_root, "provenance")
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(manifest, file.path(provenance_dir, "control_heatmap_plot_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Done. Source models, coefficients, and FDR results were not modified.")
