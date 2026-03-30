#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# plot_volcano_aggregated_annotated.R
#
# Aggregate the per-sample volcano outputs and generate an annotated summary
# plot showing recurrent ORFs and the overall distribution of the signal.
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage:\n")
  cat("  Rscript plot_volcano_aggregated_annotated.R <volcano_root_dir> <outdir> [min_replicates] [filter_label]\n")
  quit(status = 1)
}

volcano_root <- args[1]
outdir <- args[2]
min_replicates <- ifelse(length(args) >= 3, as.integer(args[3]), 2)
filter_label <- ifelse(length(args) >= 4, args[4], "unspecified")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(volcano_root)) stop("Volcano root directory not found.")

infer_condition <- function(sample_name) {
  m <- regexec("P25-([0-9]+)", sample_name)
  reg <- regmatches(sample_name, m)[[1]]
  if (length(reg) < 2) return(NA_character_)
  idx <- suppressWarnings(as.integer(reg[2]))
  if (is.na(idx)) return(NA_character_)
  if (idx >= 1 && idx <= 5) return("cold")
  if (idx >= 6 && idx <= 10) return("hot")
  return(NA_character_)
}

table_files <- list.files(volcano_root, pattern = "volcano_table_improved\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(table_files) == 0) stop("No volcano_table_improved.tsv files found.")

all_dt <- rbindlist(lapply(table_files, function(f) {
  dt <- fread(f, sep = "\t")
  sample_name <- basename(dirname(f))
  dt[, sample := sample_name]
  dt[, condition := infer_condition(sample_name)]
  dt
}), use.names = TRUE, fill = TRUE)

all_dt <- all_dt[!is.na(condition)]
all_dt[, neglog10_padj := -log10(padj)]
all_dt[!is.finite(neglog10_padj), neglog10_padj := NA_real_]
plot_dt <- all_dt[!is.na(log2OR) & !is.na(neglog10_padj)]

orf_support <- plot_dt[
  region_type == "Genic" & is_sig == TRUE & !is.na(ORF_ID) & ORF_ID != "unknown" & ORF_ID != "intergenic",
  .(
    n_replicates = uniqueN(sample),
    n_points = .N,
    median_log2OR = median(log2OR, na.rm = TRUE),
    median_neglog10_padj = median(neglog10_padj, na.rm = TRUE)
  ),
  by = .(condition, ORF_ID)
]

fwrite(orf_support, file.path(outdir, paste0("aggregated_orf_replicate_support_", filter_label, ".tsv")), sep = "\t")
label_dt <- orf_support[n_replicates >= min_replicates]

for (cond in c("cold", "hot")) {
  sub <- plot_dt[condition == cond]
  lab_sub <- label_dt[condition == cond]
  if (nrow(sub) == 0) next

  p <- ggplot(sub, aes(x = log2OR, y = neglog10_padj)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.20, contour_var = "ndensity", show.legend = FALSE) +
    geom_point(aes(color = class), alpha = 0.45, size = 1.0) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
    theme_bw() +
    labs(
      title = paste("Aggregated annotated volcano -", cond),
      subtitle = paste("Density layer + recurrent ORFs (replicates >=", min_replicates, ") | filter:", filter_label),
      x = "log2(odds ratio) (P27 vs P25)",
      y = "-log10(padj)",
      color = "Class"
    )

  if (nrow(lab_sub) > 0) {
    p <- p +
      geom_point(data = lab_sub, aes(x = median_log2OR, y = median_neglog10_padj), inherit.aes = FALSE, color = "black", size = 2.8, shape = 21, fill = "white", stroke = 0.8)

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab_sub,
        aes(x = median_log2OR, y = median_neglog10_padj, label = ORF_ID),
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = 30,
        box.padding = 0.4,
        point.padding = 0.3,
        segment.color = "grey40"
      )
    } else {
      p <- p + geom_text(
        data = lab_sub,
        aes(x = median_log2OR, y = median_neglog10_padj, label = ORF_ID),
        inherit.aes = FALSE,
        size = 3,
        vjust = -0.6,
        check_overlap = TRUE
      )
    }
  }

  ggsave(file.path(outdir, paste0("volcano_", cond, "_annotated_density_", filter_label, ".png")), p, width = 10, height = 7, dpi = 300)
}

cat("[OK] Annotated aggregated volcano plots written to:", normalizePath(outdir), "\n")
