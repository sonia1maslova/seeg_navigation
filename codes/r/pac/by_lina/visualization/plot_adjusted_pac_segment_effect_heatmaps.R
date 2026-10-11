# Display adjusted PAC segment effects; never modify source statistics.
for (p in c("dplyr", "readr", "ggplot2", "jsonlite", "grid")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
}
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) stop("Usage: Rscript script.R PROJECT_ROOT [SCOPE] [replace]")
root <- normalizePath(if (length(args)) args[[1]] else
  "/nfs/visu/workplace/yuhdua/projects/seeg_navigation", mustWork = TRUE)
allowed_scopes <- c("exclude_terminal_corridor_and_junction", "include_last")
scope <- if (length(args) >= 2L) args[[2]] else allowed_scopes[[1]]
if (!scope %in% allowed_scopes) stop("Unsupported scope: ", scope)
replace <- length(args) >= 3L && identical(args[[3]], "replace")
if (length(args) >= 3L && !replace) stop("Third argument must be replace, or omitted")
result_root <- file.path(root, "results/pac/by_lina/segment_effect_adjusted", scope)
output_root <- file.path(root, "figures/pac/by_lina/segment_effect_adjusted", scope)
rois <- c("Temporal_Mid L_Anterior", "Hippocampus L_Anterior",
          "Hippocampus L_Posterior", "Hippocampus R_Anterior")
labels <- setNames(c("aMTG", "Hipp L Ant", "Hipp L Post", "Hipp R Ant"), rois)
phase_levels <- c("slow_0p5_4Hz", "theta_broad_4_12Hz", "beta_broad_12_30Hz")
amp_levels <- c("low_gamma_30_70Hz", "high_gamma_70_150Hz")
frequency_labels <- c("0.5–4 → 30–70", "4–12 → 30–70", "12–30 → 30–70",
                      "0.5–4 → 70–150", "4–12 → 70–150", "12–30 → 70–150")
expected_families <- c(
  paste("hpc_main", rois[-1L], "not_applicable", sep = "||"),
  paste("amtg_main", rois[[1]], "not_applicable", sep = "||"),
  paste("power_control", rois[-1L], "not_applicable", sep = "||"),
  unlist(lapply(c("1-30 Hz", "1-150 Hz"), function(rng)
    paste("aperiodic_control", rois[-1L], rng, sep = "||")), use.names = FALSE)
)

read_checked_statistics <- function(selected_scope) {
  selected_root <- file.path(root, "results/pac/by_lina/segment_effect_adjusted", selected_scope)
  if (!file.exists(file.path(selected_root, "_SUCCESS"))) {
    stop("Model analysis has not passed all checks for ", selected_scope)
  }
  path <- file.path(selected_root, "all_segment_effects.csv")
  md5 <- tools::md5sum(path)
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  required <- c("branch", "scope", "roi", "phase_band", "amp_band", "exponent_range",
                "term", "beta", "p_value", "q_value", "fit_status", "fdr_family_id",
                "m_planned", "m_valid", "predictor", "target_term")
  if (!all(required %in% names(x))) stop("Missing statistic columns in ", path)
  if (nrow(x) != 78L || anyNA(x[required[!required %in% "exponent_range"]]) ||
      any(x$fit_status != "ok") || any(x$scope != selected_scope) ||
      any(x$term != "segment_position") || any(x$target_term != "segment_position") ||
      any(x$predictor != "segment_position") || any(!is.finite(x$beta)) ||
      any(!is.finite(x$p_value)) || any(!is.finite(x$q_value)) ||
      any(x$p_value < 0 | x$p_value > 1 | x$q_value < 0 | x$q_value > 1)) {
    stop("Incomplete or invalid model results: ", selected_scope)
  }
  planned <- c(hpc_main = 18L, amtg_main = 6L, power_control = 18L, aperiodic_control = 36L)
  branch_counts <- table(factor(x$branch, levels = names(planned)))
  if (any(!x$branch %in% names(planned)) || any(as.integer(branch_counts) != planned)) {
    stop("Unexpected model branch grid: ", selected_scope)
  }
  if (any(x$branch == "aperiodic_control" &
          (is.na(x$exponent_range) | !x$exponent_range %in% c("1-30 Hz", "1-150 Hz"))) ||
      any(x$branch != "aperiodic_control" & !is.na(x$exponent_range) & nzchar(x$exponent_range))) {
    stop("Unexpected aperiodic ranges: ", selected_scope)
  }
  x <- x |> mutate(
    combo_index = match(phase_band, phase_levels) + 3L * (match(amp_band, amp_levels) - 1L),
    checked_family_key = paste(branch, roi,
      if_else(branch == "aperiodic_control", exponent_range, "not_applicable"), sep = "||")
  )
  if (anyNA(x$combo_index) || any(!nzchar(x$fdr_family_id)) ||
      !setequal(unique(x$checked_family_key), expected_families) ||
      length(unique(x$fdr_family_id)) != 13L) {
    stop("Expected 13 independent six-combination FDR families: ", selected_scope)
  }
  audit <- lapply(split(x, x$fdr_family_id), function(family) {
    if (nrow(family) != 6L || length(unique(family$checked_family_key)) != 1L ||
        !setequal(family$combo_index, seq_len(6L)) ||
        any(family$m_planned != 6L) || any(family$m_valid != 6L)) {
      stop("Invalid six-test FDR family: ", family$fdr_family_id[[1]])
    }
    delta <- max(abs(family$q_value - p.adjust(family$p_value, "BH", n = 6L)))
    if (delta > 1e-12) stop("FDR check failed: ", family$fdr_family_id[[1]])
    list(branch = family$branch[[1]], roi = family$roi[[1]],
         exponent_range = if (family$branch[[1]] == "aperiodic_control")
           family$exponent_range[[1]] else "not_applicable",
         tests = 6L, frequency_combinations_complete = TRUE, max_abs_q_difference = delta)
  })
  family_map <- x |> distinct(checked_family_key, fdr_family_id)
  if (nrow(family_map) != 13L) stop("FDR family identifiers are not one-to-one")
  list(statistics = x, fdr_audit = audit, path = path, md5 = md5)
}

make_plotting_data <- function(statistics) {
  shown <- bind_rows(
    statistics |> filter(branch %in% c("hpc_main", "amtg_main")) |> mutate(analysis = "main_pac"),
    statistics |> filter(branch == "power_control") |> mutate(analysis = "power_control"),
    statistics |> filter(branch == "aperiodic_control", exponent_range == "1-150 Hz") |>
      mutate(analysis = "aperiodic_control_1_150Hz")
  ) |> mutate(
    frequency_label = frequency_labels[combo_index], roi_display = unname(labels[roi]),
    beta_display = beta / 1e-5, raw_sig = p_value < 0.05, fdr_sig = q_value < 0.05,
    sig_symbol = case_when(fdr_sig ~ "†", raw_sig ~ "*", TRUE ~ "")
  )
  if (nrow(shown) != 60L || anyNA(shown$roi_display) || anyNA(shown$combo_index)) {
    stop("Invalid plotting grid")
  }
  if (any((shown |> count(analysis, roi, combo_index))$n != 1L)) stop("Duplicated heatmap cells")
  shown
}

checked_inputs <- setNames(list(read_checked_statistics(scope)), scope)
other_scope <- setdiff(allowed_scopes, scope)[[1]]
other_root <- file.path(root, "results/pac/by_lina/segment_effect_adjusted", other_scope)
if (file.exists(file.path(other_root, "_SUCCESS"))) {
  checked_inputs[[other_scope]] <- read_checked_statistics(other_scope)
}
shown <- make_plotting_data(checked_inputs[[scope]]$statistics)
all_color_data <- bind_rows(lapply(checked_inputs, function(x) make_plotting_data(x$statistics)))
max_beta_display <- max(abs(all_color_data$beta_display))
limit <- signif(max_beta_display, 2)
if (limit < max_beta_display) limit <- signif(max_beta_display * 1.05, 2)
if (!is.finite(limit) || limit <= 0) limit <- 1
if (!replace && dir.exists(output_root) &&
    length(list.files(output_root, recursive = TRUE, all.files = TRUE, no.. = TRUE))) {
  stop("Figure files exist; pass replace explicitly to overwrite this adjusted scope")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
success_path <- file.path(output_root, "_SUCCESS")
if (replace && file.exists(success_path) && unlink(success_path) != 0L) {
  stop("Cannot invalidate previous figure completion marker")
}

scope_description <- if (scope == "include_last") {
  "include_last：保留每个试次的末 corridor（两个 2.5 s 子段）及末 junction。"
} else {
  "exclude_terminal_corridor_and_junction：同时排除每个试次末 corridor 的两个 2.5 s 子段及末 junction，保留其余 segment 的原始导航序号。"
}
correction_text <- c(
  "PAC 随 segment 序号变化：FDR 校正与图中显著性说明", "",
  paste0("分析范围：", scope_description),
  "检验对象：线性混合效应模型中 segment_position 固定效应的双侧原始 p 值。",
  "校正方法：Benjamini–Hochberg（BH）；以 q < 0.05 为 FDR 显著性阈值。",
  "每个独立检验族固定为 m = 6：同一脑区内 3 个低频 phase 频段 × 2 个高频 amplitude 频段。",
  "phase 频段：0.5–4、4–12、12–30 Hz；amplitude 频段：30–70、70–150 Hz。",
  "检验族按分析类型、脑区及（适用时）1/f 指数估计范围分别建立；不同脑区、不同分析类型、不同估计范围及 include/exclude 范围均不合并校正。",
  "不按实验 condition 拆分，不将同一检验因不同展示视图重复计入 FDR。", "",
  "主 PAC：4 个脑区各自按 6 项校正，共 4 个检验族、24 项检验。",
  "功率控制：3 个海马脑区各自按 6 项校正，共 3 个检验族、18 项检验。",
  "非周期性谱成分控制：3 个海马脑区 × 2 个指数估计范围（1–30 Hz、1–150 Hz），各自按 6 项校正，共 6 个检验族、36 项检验。",
  "每个分析范围合计 13 个检验族、78 项检验；不是将全部 78 项一起校正。",
  "当前非周期性控制图仅呈现 1–150 Hz；该图每个脑区的校正仅包含同一 1–150 Hz 范围的 6 项，不与 1–30 Hz 合并。", "",
  "图中 †：BH-FDR q < 0.05；否则 *：原始 p < 0.05；否则无符号。† 与 * 不重复标注。",
  "颜色表示原始 segment_position 回归系数（展示单位 β / 10^-5），不表示 p 值或 q 值。",
  paste0("色阶：共同的对称零中心色阶；参与确定色阶的已完成范围为 ",
         paste(names(checked_inputs), collapse = "、"), "。"),
  "图为 PNG；CSV 为相应绘图数据。原始李娜结果和输入统计表未被修改。"
)
root_txt <- file.path(output_root, "FDR_CORRECTION_README.txt")
readr::write_lines(correction_text, root_txt)
summary <- list()
files <- root_txt
for (analysis in c("main_pac", "power_control", "aperiodic_control_1_150Hz")) {
  roi_order <- if (analysis == "main_pac") rois else rois[-1L]
  dat <- shown |> filter(.data$analysis == .env$analysis) |> mutate(
    x_pos = 1.00 + 0.22 * (combo_index - 1L),
    y_pos = 2.08 - 0.22 * (match(roi, roi_order) - 1L))
  if (nrow(dat) != 6L * length(roi_order)) stop("Missing cells for ", analysis)
  roi_text <- data.frame(x = 0.87, y = 2.08 - 0.22 * (seq_along(roi_order) - 1L),
                         label = unname(labels[roi_order]))
  frequency_text <- data.frame(x = 1.00 + 0.22 * (0:5),
                               y = min(dat$y_pos) - 0.11 - 0.035, label = frequency_labels)
  four_rows <- length(roi_order) == 4L
  plot <- ggplot(dat, aes(x = x_pos, y = y_pos, fill = beta_display)) +
    geom_tile(width = 0.22, height = 0.22, color = "white", linewidth = 0.85) +
    geom_text(aes(label = sig_symbol), size = 5.0, fontface = "bold") +
    geom_text(data = roi_text, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 4.2) +
    geom_text(data = frequency_text, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1, fontface = "bold", size = 3.4) +
    annotate("text", x = 1.55, y = 2.30, label = "Segment number", fontface = "bold", size = 6.2) +
    annotate("text", x = 1.55, y = if (four_rows) 0.72 else 0.90,
              label = "PAC frequency combination (phase → amplitude, Hz)", fontface = "bold", size = 4.7) +
    annotate("text", x = 1.55, y = if (four_rows) 0.53 else 0.68,
              label = "* = raw p < 0.05      † = BH-FDR q < 0.05", fontface = "bold", size = 3.4, color = "grey20") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
      limits = c(-limit, limit), name = expression(beta / 10^{-5}),
      guide = guide_colorbar(direction = "vertical", title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(0.55, "cm"), barheight = grid::unit(if (four_rows) 4.0 else 3.0, "cm"),
        ticks = TRUE, frame.colour = "grey55")) +
    scale_x_continuous(breaks = NULL, limits = c(-0.35, 3.25), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = NULL, limits = c(0.32, 2.82), expand = expansion(mult = c(0, 0))) +
    labs(title = NULL, x = NULL, y = NULL) + coord_fixed(ratio = 1, clip = "off") +
    theme_minimal(base_size = 14) + theme(
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "inside", legend.position.inside = c(0.77, if (four_rows) 0.58 else 0.615),
      legend.justification = c(0, 0.5), legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 13), legend.text = element_text(size = 11, color = "black"),
      legend.background = element_rect(fill = "white", color = NA), legend.box.background = element_blank(),
      legend.box.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 2, 2, 2))
  out <- file.path(output_root, analysis)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("heatmap_Segment_number_all_segments_", scope)
  png <- file.path(out, paste0(stem, ".png"))
  csv <- file.path(out, paste0(stem, "_values.csv"))
  ggsave(png, plot, width = 8, height = 5.5, units = "in", dpi = 300, bg = "white")
  readr::write_csv(dat, csv, na = "")
  analysis_description <- switch(analysis,
    main_pac = "本目录展示主 PAC 的 4 个脑区；每个脑区单独校正 6 项。",
    power_control = "本目录展示功率控制的 3 个海马脑区；每个脑区单独校正 6 项。",
    aperiodic_control_1_150Hz = "本目录展示 1–150 Hz 指数控制的 3 个海马脑区；每个脑区在本估计范围内单独校正 6 项，不与 1–30 Hz 合并。")
  txt <- file.path(out, "FDR_CORRECTION_README.txt")
  readr::write_lines(c(analysis_description, "", correction_text), txt)
  files <- c(files, png, csv, txt)
  summary[[analysis]] <- list(effects = nrow(dat), raw_significant = sum(dat$raw_sig),
                              fdr_significant = sum(dat$fdr_sig))
  message("Saved ", png)
}
for (checked in checked_inputs) {
  if (!identical(checked$md5, tools::md5sum(checked$path))) {
    stop("Source statistics changed during plotting: ", checked$path)
  }
}
jsonlite::write_json(list(input = checked_inputs[[scope]]$path,
  input_md5 = as.list(checked_inputs[[scope]]$md5), scope = scope,
  reference_figure = file.path(root, "figures/power/encoding_lmm/heatmaps/heatmap_Segment_number_all_segments_exclude_last.png"),
  roi_labels = as.list(labels), fdr_audit = checked_inputs[[scope]]$fdr_audit,
  fdr_method = "BH", fdr_family_size = 6L, fdr_family_count = 13L,
  fdr_family_dimensions = c("analysis branch", "ROI", "exponent range when applicable", "scope"),
  frequency_label_gap = 0.035,
  beta_display_unit = 1e-5, shared_raw_beta_limits = c(-limit, limit) * 1e-5,
  shared_color_scopes = names(checked_inputs),
  shared_color_sources = lapply(checked_inputs, function(x)
    list(input = x$path, input_md5 = as.list(x$md5))),
  significance = "dagger for BH q < .05; otherwise star for raw p < .05; otherwise blank",
  exponent_ranges_corrected_separately = c("1-30 Hz", "1-150 Hz"),
  exponent_range_shown = "1-150 Hz", replacement_explicitly_requested = replace,
  summary = summary, files = files), file.path(output_root, "plot_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("Three PNG heatmaps, plotting data and correction TXT files exported; 13 six-test BH families verified.",
           success_path)
message("Completed three adjusted heatmaps.")
