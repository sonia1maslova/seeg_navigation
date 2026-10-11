# PNG-only heatmaps: same compact S03 angle-heatmap geometry and significance
# notation. aMTG + 3 hippocampal ROIs only; pMTG raw effects remain in the tables.

entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(entry)) stop("Run with Rscript.", call. = FALSE)
source(file.path(dirname(sub("^--file=", "", entry[1])), "bearing_branch_common.R"))
branch_require(c("dplyr", "readr", "ggplot2", "jsonlite"))
library(dplyr)
library(readr)
library(ggplot2)
cfg <- branch_config()
results_root <- file.path(cfg$results$power, "encoding_lmm", branch_name)
output_dir <- file.path(cfg$figures$power, "encoding_lmm", branch_name)
branch_assert_isolated(output_dir, file.path(cfg$figures$power, "encoding_lmm"))
if (!file.exists(file.path(results_root, "_SUCCESS"))) stop("Supplemental model run not completed.", call. = FALSE)
if (dir.exists(output_dir)) stop("Supplemental figure directory exists; no automatic overwrite/cleanup.", call. = FALSE)
df <- read_csv(file.path(results_root, "summaries", "all_target_predictor_effects.csv"),
               show_col_types = FALSE, progress = FALSE)
band_labels <- c("0.5–4", "4–12", "12–30", "30–70", "70–150")
roi_display_map <- c("Temporal_Mid L_Anterior" = "aMTG",
  "Hippocampus L_Anterior" = "Hipp L Ant", "Hippocampus L_Posterior" = "Hipp L Post",
  "Hippocampus R_Anterior" = "Hipp R Ant")
tile_w <- 0.22; tile_h <- 0.22
band_x_map <- setNames(c(1.00, 1.22, 1.44, 1.66, 1.88), band_labels)
roi_y_map <- c("aMTG" = 1.98, "Hipp L Ant" = 1.65, "Hipp L Post" = 1.43, "Hipp R Ant" = 1.21)
template <- expand.grid(roi = branch_rois_plot, band_name = branch_bands, stringsAsFactors = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_plot <- function(predictor, scope, variant) {
  actual <- df |> filter(.data$predictor == .env$predictor, .data$scope == .env$scope,
    .data$variant == .env$variant, .data$roi %in% branch_rois_plot,
    .data$band_name %in% branch_bands, .data$fit_status == "ok")
  if (anyDuplicated(branch_key(actual, c("roi", "band_name")))) stop("Duplicate ROI/band result.", call. = FALSE)
  dat <- template |> left_join(actual, by = c("roi", "band_name")) |>
    mutate(scope = .env$scope, variant = .env$variant, predictor = .env$predictor,
      roi_display = unname(roi_display_map[roi]),
      band_label = factor(band_labels[match(band_name, branch_bands)], levels = band_labels),
      raw_sig = is.finite(p_value_raw) & p_value_raw < 0.05,
      fdr_sig = is.finite(p_fdr) & p_fdr < 0.05,
      sig_symbol = case_when(fdr_sig ~ "†", raw_sig ~ "*", TRUE ~ ""),
      y_pos = unname(roi_y_map[roi_display]), x_pos = unname(band_x_map[as.character(band_label)]))
  finite <- dat$estimate[is.finite(dat$estimate)]
  fill_limit <- if (length(finite)) max(abs(finite)) else 1
  if (fill_limit == 0) fill_limit <- 1
  legend_limit <- signif(fill_limit, 2)
  if (legend_limit < fill_limit) legend_limit <- signif(fill_limit * 1.05, 2)
  title <- if (predictor == "Angle2orientation") "Egocentric bearing" else "Allocentric bearing"
  roi_labels <- data.frame(label = names(roi_y_map), x = 0.87, y = unname(roi_y_map))
  band_text <- data.frame(label = band_labels, x = unname(band_x_map), y = 1.01)
  center_x <- mean(range(unname(band_x_map)))
  p <- ggplot(dat, aes(x = x_pos, y = y_pos, fill = estimate)) +
    geom_tile(width = tile_w, height = tile_h, color = "white", linewidth = 0.85) +
    geom_text(aes(label = sig_symbol), size = 5.0, fontface = "bold", vjust = 0.5) +
    annotate("segment", x = min(unname(band_x_map)) - tile_w / 2,
      xend = max(unname(band_x_map)) + tile_w / 2, y = 1.815, yend = 1.815,
      linewidth = 0.70, color = "grey42") +
    geom_text(data = roi_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE,
      hjust = 1, fontface = "bold", size = 4.2, color = "black") +
    geom_text(data = band_text, aes(x = x, y = y, label = label), inherit.aes = FALSE,
      angle = 45, hjust = 1, vjust = 1, fontface = "bold", size = 3.4, color = "black") +
    annotate("text", x = center_x, y = 2.20, label = title, fontface = "bold", size = 6.2) +
    annotate("text", x = center_x, y = 0.76, label = "Frequency band (Hz)", fontface = "bold", size = 4.7) +
    annotate("text", x = center_x, y = 0.57, label = "* = raw p < 0.05      † = BH-FDR q < 0.05",
      fontface = "bold", size = 3.4, color = "grey20") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
      limits = c(-legend_limit, legend_limit), name = expression(beta), na.value = "grey90",
      guide = guide_colorbar(direction = "vertical", title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(0.55, "cm"), barheight = grid::unit(5.1, "cm"), ticks = TRUE,
        frame.colour = "grey55")) +
    scale_x_continuous(breaks = NULL, limits = c(-0.35, 3.25), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = NULL, limits = c(0.39, 2.75), expand = expansion(mult = c(0, 0))) +
    labs(title = NULL, x = NULL, y = NULL, caption = NULL) +
    coord_fixed(ratio = 1, clip = "off") + theme_minimal(base_size = 14) +
    theme(plot.title = element_blank(), plot.caption = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), legend.position = c(0.69, 0.56),
      legend.justification = c(0, 0.5), legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11, color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_blank(), legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2))
  base <- paste0("heatmap_", gsub(" ", "_", title), "_", scope, "_", variant)
  png <- file.path(output_dir, "heatmaps", paste0(base, ".png"))
  if (file.exists(png)) stop("Refusing to overwrite: ", png, call. = FALSE)
  dir.create(dirname(png), recursive = TRUE, showWarnings = FALSE)
  ggsave(png, plot = p, width = 8.0, height = 5.5, units = "in", dpi = 300, bg = "white")
  branch_write_csv(dat |> select(roi, roi_display, scope, variant, predictor, band_name, band_label,
    estimate, p_value_raw, p_fdr, fdr_family, fdr_planned_tests, n_observations, n_subjects,
    raw_sig, fdr_sig, sig_symbol), file.path(output_dir, "heatmaps", paste0(base, "_values.csv")))
  message("Saved ", png)
}
for (scope in names(branch_predictors)) {
  for (predictor in branch_predictors[[scope]]) {
    for (variant in c("include_last", "exclude_last")) make_plot(predictor, scope, variant)
  }
}
branch_write_lines(c(
  "Position-to-start bearing sensitivity plots; original plots retained.",
  "PNG-only; exact established S03 angle-heatmap geometry, colors and symbols.",
  "FDR separately within each angle, scope and include/exclude branch: aMTG 5 band tests; 3 hippocampal ROIs jointly 15 tests.",
  "pMTG is fitted and compared in result tables, but has no established angle-plot FDR family and is not plotted.",
  "Gray cells denote missing/failed fits, not zero effects. No new scientific threshold or family is introduced.",
  "* denotes raw p<0.05 without BH significance; dagger denotes BH q<0.05."
), file.path(output_dir, "FDR_CORRECTION_README.txt"))
branch_write_lines("8 heatmaps written; models and comparison tables are in the corresponding results branch.",
                   file.path(output_dir, "_SUCCESS"))
