#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# plot_isec_results.R
#
# Read the TSV summaries produced by extract_isec_summary.sh and turn them into
# report-ready figures describing the SNP results.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "both")

# Consistent colors used across all plots.
col_froid  <- "#378ADD"
col_chaud  <- "#D85A30"
col_shared <- "#888780"
col_high     <- "#A32D2D"
col_moderate <- "#BA7517"
col_low      <- "#639922"
col_modifier <- "#888780"

impact_colors <- c(
  "HIGH" = col_high,
  "MODERATE" = col_moderate,
  "LOW" = col_low,
  "MODIFIER" = col_modifier
)

# -----------------------------------------------------------------------------
# Read all summary tables for one filter mode
# -----------------------------------------------------------------------------
read_summary <- function(filter_mode) {
  root <- sprintf("results/isec_%s/summary", filter_mode)
  if (!dir.exists(root)) {
    cat(sprintf("ERROR: %s not found\n", root))
    return(NULL)
  }

  list(
    counts = fread(file.path(root, "counts.tsv")),
    genes_chaud = if (file.exists(file.path(root, "genes_specific_chaud.tsv")))
      fread(file.path(root, "genes_specific_chaud.tsv")) else NULL,
    genes_froid = if (file.exists(file.path(root, "genes_specific_froid.tsv")))
      fread(file.path(root, "genes_specific_froid.tsv")) else NULL,
    genes_shared = if (file.exists(file.path(root, "genes_shared_froid_chaud.tsv")))
      fread(file.path(root, "genes_shared_froid_chaud.tsv")) else NULL,
    missense = fread(file.path(root, "missense_all.tsv")),
    individual = fread(file.path(root, "individual.tsv"))
  )
}

# -----------------------------------------------------------------------------
# FIGURE 1: Stacked barplot of specific/shared variants and impact categories
# -----------------------------------------------------------------------------
plot_variant_counts <- function(counts_dt, filter_label, outdir) {
  counts_dt[, category_label := fcase(
    category == "specific_froid", "Specific\ncold",
    category == "specific_chaud", "Specific\nhot",
    category == "shared_froid_chaud", "Shared"
  )]

  counts_dt[, impact := factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))]
  counts_dt[, category_label := factor(category_label,
    levels = c("Specific\ncold", "Shared", "Specific\nhot"))]

  totals <- counts_dt[, .(total = sum(count)), by = category_label]

  p <- ggplot(counts_dt, aes(x = category_label, y = count, fill = impact)) +
    geom_bar(stat = "identity", width = 0.65) +
    geom_text(
      data = totals,
      aes(x = category_label, y = total, label = total, fill = NULL),
      vjust = -0.5,
      size = 4.5,
      fontface = "bold",
      color = "grey30"
    ) +
    scale_fill_manual(values = impact_colors, name = "Impact") +
    labs(
      title = sprintf("Condition-specific variants — %s filter", filter_label),
      subtitle = "Thermal shock at P26: comparison between P25 (before) and P27 (after)",
      x = NULL,
      y = "Number of variants"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey50", size = 11),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12)
    ) +
    coord_cartesian(ylim = c(0, max(totals$total) * 1.15))

  ggsave(file.path(outdir, "fig1_variant_counts.png"), p, width = 8, height = 6, dpi = 300)
  cat("  -> fig1_variant_counts.png\n")
}

# -----------------------------------------------------------------------------
# FIGURE 2: Top affected ORFs by condition and impact type
# -----------------------------------------------------------------------------
plot_top_genes <- function(genes_chaud, genes_froid, filter_label, outdir) {
  parts <- list()
  if (!is.null(genes_chaud) && nrow(genes_chaud) > 0) {
    g <- genes_chaud[, .N, by = .(gene, impact)]
    g[, condition := "Hot"]
    parts[[length(parts) + 1]] <- g
  }
  if (!is.null(genes_froid) && nrow(genes_froid) > 0) {
    g <- genes_froid[, .N, by = .(gene, impact)]
    g[, condition := "Cold"]
    parts[[length(parts) + 1]] <- g
  }

  if (length(parts) == 0) return()
  combined <- rbindlist(parts)

  top_genes <- combined[, .(total = sum(N)), by = gene][order(-total)][1:min(15, .N)]
  combined <- combined[gene %in% top_genes$gene]
  combined[, gene := factor(gene, levels = rev(top_genes$gene))]
  combined[, impact := factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))]

  p <- ggplot(combined, aes(x = gene, y = N, fill = impact)) +
    geom_bar(stat = "identity") +
    facet_wrap(~condition, scales = "free_x") +
    scale_fill_manual(values = impact_colors, name = "Impact") +
    coord_flip() +
    labs(
      title = sprintf("Most affected ORFs — %s filter", filter_label),
      subtitle = "Number of variants per gene and impact class",
      x = NULL,
      y = "Number of variants"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey50", size = 11),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(outdir, "fig2_top_genes.png"), p, width = 12, height = 7, dpi = 300)
  cat("  -> fig2_top_genes.png\n")
}

# -----------------------------------------------------------------------------
# FIGURE 3: Genome-wide missense positions (lollipop-like plot)
# -----------------------------------------------------------------------------
plot_missense_genome <- function(missense_dt, filter_label, outdir, genome_length = 295146) {
  if (nrow(missense_dt) == 0) {
    cat("  [SKIP] No missense variants\n")
    return()
  }

  missense_dt[, category_label := fcase(
    category == "specific_froid", "Cold",
    category == "specific_chaud", "Hot",
    category == "shared_froid_chaud", "Shared"
  )]

  missense_dt[, pos := as.numeric(pos)]
  missense_dt[, category_label := factor(category_label, levels = c("Cold", "Shared", "Hot"))]

  condition_colors <- c("Cold" = col_froid, "Hot" = col_chaud, "Shared" = col_shared)

  p <- ggplot(missense_dt, aes(x = pos, y = category_label, color = category_label)) +
    geom_segment(aes(xend = pos, yend = as.numeric(category_label) - 0.3), linewidth = 0.4) +
    geom_point(size = 3) +
    geom_text(
      aes(label = gene),
      size = 2.5,
      hjust = -0.1,
      vjust = -1,
      color = "grey30",
      check_overlap = TRUE
    ) +
    scale_color_manual(values = condition_colors, name = "Condition") +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1000), "kb"),
      breaks = seq(0, genome_length, by = 50000),
      limits = c(0, genome_length)
    ) +
    labs(
      title = sprintf("Missense / high-impact variants across the CyHV-3 genome — %s filter", filter_label),
      subtitle = "Genomic positions of variants with moderate or high impact",
      x = "Genome position (kb)",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey50", size = 11),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold")
    )

  ggsave(file.path(outdir, "fig3_missense_genome.png"), p, width = 14, height = 5, dpi = 300)
  cat("  -> fig3_missense_genome.png\n")
}

# -----------------------------------------------------------------------------
# FIGURE 4: Per-sample comparison of private/shared variants
# -----------------------------------------------------------------------------
plot_individual <- function(indiv_dt, filter_label, outdir) {
  if (nrow(indiv_dt) == 0) return()

  indiv_dt <- indiv_dt[group != "exclu"]

  long <- melt(
    indiv_dt,
    id.vars = c("filter", "sample", "group"),
    measure.vars = c("private_P25", "private_P27", "shared"),
    variable.name = "type",
    value.name = "count"
  )

  long[, type_label := fcase(
    type == "private_P25", "Lost (private P25)",
    type == "private_P27", "Appeared (private P27)",
    type == "shared", "Shared"
  )]
  long[, type_label := factor(type_label,
    levels = c("Lost (private P25)", "Shared", "Appeared (private P27)"))]

  long[, sample_short := sub(".*-(\\d+)\\..*", "\\1", sample)]
  long[, sample_short := factor(sample_short, levels = as.character(sort(as.numeric(unique(sample_short)))))]

  type_colors <- c(
    "Lost (private P25)" = "#85B7EB",
    "Shared" = "#B4B2A9",
    "Appeared (private P27)" = "#F0997B"
  )

  p <- ggplot(long, aes(x = sample_short, y = count, fill = type_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(
      ~group,
      scales = "free_x",
      labeller = labeller(group = c("froid" = "Cold group", "chaud" = "Hot group"))
    ) +
    scale_fill_manual(values = type_colors, name = NULL) +
    labs(
      title = sprintf("P25 vs P27 comparison per sample — %s filter", filter_label),
      subtitle = "Variants lost, shared or gained after thermal shock at P26",
      x = "Sample",
      y = "Number of variants"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey50", size = 11),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(outdir, "fig4_individual_comparison.png"), p, width = 12, height = 6, dpi = 300)
  cat("  -> fig4_individual_comparison.png\n")
}

# -----------------------------------------------------------------------------
# FIGURE 5: Comparison between moderate and strict filters
# -----------------------------------------------------------------------------
plot_filter_comparison <- function(counts_strict, counts_moderate, outdir) {
  counts_strict[, filter_label := "Strict"]
  counts_moderate[, filter_label := "Moderate"]
  combined <- rbindlist(list(counts_strict, counts_moderate))

  combined[, category_label := fcase(
    category == "specific_froid", "Cold",
    category == "specific_chaud", "Hot",
    category == "shared_froid_chaud", "Shared"
  )]
  combined[, category_label := factor(category_label, levels = c("Cold", "Shared", "Hot"))]

  totals <- combined[, .(total = sum(count)), by = .(filter_label, category_label)]

  p <- ggplot(totals, aes(x = category_label, y = total, fill = filter_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(
      aes(label = total),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 4,
      fontface = "bold",
      color = "grey30"
    ) +
    scale_fill_manual(values = c("Strict" = "#E24B4A", "Moderate" = "#EF9F27"), name = "Filter") +
    labs(
      title = "Filter comparison: strict vs moderate",
      subtitle = "Total number of variants by category and filter level",
      x = NULL,
      y = "Number of variants"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey50", size = 11),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, max(totals$total) * 1.15))

  ggsave(file.path(outdir, "fig5_filter_comparison.png"), p, width = 8, height = 6, dpi = 300)
  cat("  -> fig5_filter_comparison.png\n")
}

# -----------------------------------------------------------------------------
# Main driver for one filter mode
# -----------------------------------------------------------------------------
generate_figures <- function(filter_mode) {
  data <- read_summary(filter_mode)
  if (is.null(data)) return(NULL)

  figdir <- sprintf("results/isec_%s/figures", filter_mode)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  label <- ifelse(filter_mode == "strict", "strict", "moderate")
  cat(sprintf("\n=== Generating figures [%s] ===\n", filter_mode))

  plot_variant_counts(data$counts, label, figdir)
  plot_top_genes(data$genes_chaud, data$genes_froid, label, figdir)
  plot_missense_genome(data$missense, label, figdir)
  plot_individual(data$individual, label, figdir)

  cat(sprintf("  Figures written to %s/\n", figdir))
  return(data)
}

# Run either one filter mode or both.
if (mode == "both") {
  data_s <- generate_figures("strict")
  data_m <- generate_figures("moderate")

  if (!is.null(data_s) && !is.null(data_m)) {
    figdir <- "results/figures_comparison"
    dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
    cat("\n=== Comparative figure ===\n")
    plot_filter_comparison(data_s$counts, data_m$counts, figdir)
    cat(sprintf("  -> %s/\n", figdir))
  }
} else {
  generate_figures(mode)
}

cat("\n[DONE]\n")
