# Canonical Segment-number power-encoding heatmaps.
#
# The figure deliberately follows the established single-panel heatmap style:
# one beta map, one vertical beta color bar, and raw/FDR significance symbols
# in the tiles. A separate plot is written for each requested scope.

required_packages <- c("dplyr", "readr", "ggplot2", "grid")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)
library(readr)
library(ggplot2)

project_root <- "/nfs/visu/workplace/yuhdua/projects/seeg_navigation"
variant_to_plot <- "exclude_last"
scopes_to_plot <- c("corridor", "junction", "all_segments")
predictor_to_plot <- "SegmentNumber"
predictor_display_label <- "Segment number"
alpha_raw <- 0.05
alpha_fdr <- 0.05

band_levels <- c("f0p5_4", "f4_12", "f12_30", "f30_70", "f70_150")
band_labels <- c("0.5–4", "4–12", "12–30", "30–70", "70–150")

selected_rois <- c(
  "Temporal_Mid L_Anterior",
  "Temporal_Mid L_Posterior",
  "Hippocampus L_Anterior",
  "Hippocampus L_Posterior",
  "Hippocampus R_Anterior"
)

# For each predictor x scope x include/exclude branch, the three hippocampal
# subregions share a 15-test BH family. Every other ROI has its own 5-band family.
hipp_rois <- c(
  "Hippocampus L_Anterior",
  "Hippocampus L_Posterior",
  "Hippocampus R_Anterior"
)

# Keep the established display abbreviations used in the reference figures.
roi_display_map <- c(
  "Temporal_Mid L_Anterior" = "aMTG",
  "Temporal_Mid L_Posterior" = "pMTG",
  "Hippocampus L_Anterior" = "Hipp L Ant",
  "Hippocampus L_Posterior" = "Hipp L Post",
  "Hippocampus R_Anterior" = "Hipp R Ant"
)

tile_w <- 0.22
tile_h <- 0.22
band_x_map <- c("0.5–4" = 1.00, "4–12" = 1.22, "12–30" = 1.44,
                "30–70" = 1.66, "70–150" = 1.88)
roi_y_map <- c("aMTG" = 2.08, "pMTG" = 1.86, "Hipp L Ant" = 1.64,
               "Hipp L Post" = 1.42, "Hipp R Ant" = 1.20)
heatmap_center_x <- mean(range(unname(band_x_map)))

input_file <- file.path(project_root, "results", "power", "encoding_lmm",
                        "summaries", "all_target_predictor_effects.csv")
output_dir <- file.path(project_root, "figures", "power", "encoding_lmm", "heatmaps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input file not found:\n", input_file, call. = FALSE)
df <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)
required_cols <- c("roi", "scope", "variant", "band_name", "predictor",
                   "estimate", "p_value_raw", "fit_status")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

template <- expand.grid(
  roi = selected_rois,
  band_name = band_levels,
  stringsAsFactors = FALSE
) |>
  mutate(
    roi_display = unname(roi_display_map[roi]),
    band_label = unname(setNames(band_labels, band_levels)[band_name]),
    x_pos = unname(band_x_map[band_label]),
    y_pos = unname(roi_y_map[roi_display])
  )

make_scope_plot <- function(scope_name) {
  actual <- df |>
    filter(
      variant == variant_to_plot,
      scope == scope_name,
      predictor == predictor_to_plot,
      fit_status == "ok",
      roi %in% selected_rois,
      band_name %in% band_levels
    ) |>
    select(roi, band_name, estimate, p_value_raw)

  dup <- actual |>
    count(roi, band_name, name = "n") |>
    filter(n > 1L)
  if (nrow(dup) > 0) stop("Duplicated ROI x band rows in ", scope_name, call. = FALSE)

  dat <- template |>
    left_join(actual, by = c("roi", "band_name")) |>
    mutate(
      fdr_family = if_else(roi %in% hipp_rois, "hippocampus_3subregions", roi),
      fdr_planned_tests = if_else(roi %in% hipp_rois, 15L, 5L)
    ) |>
    group_by(fdr_family) |>
    mutate(
      p_fdr = p.adjust(p_value_raw, method = "BH", n = first(fdr_planned_tests))
    ) |>
    ungroup() |>
    mutate(
      raw_sig = is.finite(p_value_raw) & p_value_raw < alpha_raw,
      fdr_sig = is.finite(p_fdr) & p_fdr < alpha_fdr,
      sig_symbol = case_when(fdr_sig ~ "†", raw_sig ~ "*", TRUE ~ "")
    )

  finite_beta <- dat$estimate[is.finite(dat$estimate)]
  fill_limit <- if (length(finite_beta)) max(abs(finite_beta)) else 1
  if (!is.finite(fill_limit) || fill_limit == 0) fill_limit <- 1
  legend_limit <- signif(fill_limit, digits = 2)
  if (legend_limit < fill_limit) legend_limit <- signif(fill_limit * 1.05, digits = 2)
  if (!is.finite(legend_limit) || legend_limit <= 0) legend_limit <- fill_limit

  roi_label_df <- data.frame(
    label = names(roi_y_map), x = 0.87, y = unname(roi_y_map), stringsAsFactors = FALSE
  )
  band_label_df <- data.frame(
    label = band_labels, x = unname(band_x_map), y = 0.96, stringsAsFactors = FALSE
  )
  p <- ggplot(dat, aes(x = x_pos, y = y_pos, fill = estimate)) +
    geom_tile(width = tile_w, height = tile_h, color = "white", linewidth = 0.85) +
    geom_text(aes(label = sig_symbol), size = 5.0, fontface = "bold", vjust = 0.5) +
    geom_text(data = roi_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 4.2) +
    geom_text(data = band_label_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1,
              fontface = "bold", size = 3.4) +
    annotate("text", x = heatmap_center_x, y = 2.30,
             label = predictor_display_label, fontface = "bold", size = 6.2) +
    annotate("text", x = heatmap_center_x, y = 0.72,
             label = "Frequency band (Hz)", fontface = "bold", size = 4.7) +
    annotate("text", x = heatmap_center_x, y = 0.53,
             label = "* = raw p < 0.05      † = BH-FDR q < 0.05",
             fontface = "bold", size = 3.4, color = "grey20") +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
      limits = c(-legend_limit, legend_limit), name = expression(beta),
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
      panel.grid = element_blank(), legend.position = c(0.69, 0.56),
      legend.justification = c(0, 0.5), legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_blank(), legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )

  base_name <- paste0("heatmap_Segment_number_", scope_name, "_", variant_to_plot)
  png_file <- file.path(output_dir, paste0(base_name, ".png"))
  table_file <- file.path(output_dir, paste0(base_name, "_values.csv"))
  ggsave(png_file, p, width = 8.0, height = 5.5, units = "in", dpi = 300, bg = "white")
  write_csv(dat |>
              mutate(scope = scope_name, variant = variant_to_plot,
                     predictor = predictor_to_plot) |>
              select(roi, roi_display, scope, variant, predictor, band_name, band_label,
                     estimate, p_value_raw, p_fdr, fdr_family, fdr_planned_tests,
                     raw_sig, fdr_sig, sig_symbol),
            table_file, na = "")
  message("Saved ", png_file)
}

for (scope_name in scopes_to_plot) make_scope_plot(scope_name)
message("Done.")
