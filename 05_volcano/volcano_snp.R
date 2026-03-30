#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# volcano_snp.R
#
# Compare two depth-enriched SNP VCF files, test frequency shifts among shared
# variants, and export volcano plots plus detailed result tables.
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript volcano_snp.R <P25.vcf(.gz)> <P27.vcf(.gz)> <genome.gff3> <outdir> [minDP] [minALT] [minQUAL] [padj] [log2FC_threshold] [minVAF]\n")
  quit(status = 1)
}

p25_path <- args[1]
p27_path <- args[2]
gff3_path <- args[3]
outdir <- args[4]

# User-adjustable filtering thresholds.
minDP   <- ifelse(length(args) >= 5, as.numeric(args[5]), 50)
minALT  <- ifelse(length(args) >= 6, as.numeric(args[6]), 10)
minQUAL <- ifelse(length(args) >= 7, as.numeric(args[7]), 10)
padj_th <- ifelse(length(args) >= 8, as.numeric(args[8]), 0.05)
log2FC_threshold <- ifelse(length(args) >= 9, as.numeric(args[9]), 1)
minVAF <- ifelse(length(args) >= 10, as.numeric(args[10]), 0.10)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read a VCF body (ignore header lines)
# -----------------------------------------------------------------------------
read_vcf_body <- function(path) {
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  cmd <- if (is_gz) {
    sprintf("zcat %s | grep -v '^#'", shQuote(path))
  } else {
    sprintf("grep -v '^#' %s", shQuote(path))
  }
  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 10) stop("VCF: less than 10 columns, unexpected format.")
  setnames(dt, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"))
  dt[, POS := as.integer(POS)]
  dt[, QUAL := suppressWarnings(as.numeric(QUAL))]
  dt
}

# Extract one FORMAT field (such as DP or AD) from a sample column.
extract_format_field <- function(format_vec, sample_vec, field) {
  out <- rep(NA_character_, length(format_vec))
  for (i in seq_along(format_vec)) {
    fmt <- format_vec[i]
    smp <- sample_vec[i]
    if (is.na(fmt) || is.na(smp)) next
    keys <- strsplit(fmt, ":", fixed = TRUE)[[1]]
    vals <- strsplit(smp, ":", fixed = TRUE)[[1]]
    if (length(keys) != length(vals)) next
    idx <- match(field, keys)
    if (is.na(idx)) next
    out[i] <- vals[idx]
  }
  out
}

# Parse AD only for reliable biallelic cases (REF,ALT).
parse_biallelic_ad <- function(ad_str) {
  if (is.na(ad_str) || ad_str == "." || ad_str == "") return(list(ref=NA_real_, alt=NA_real_, ok=FALSE))
  nums <- suppressWarnings(as.numeric(strsplit(ad_str, ",", fixed = TRUE)[[1]]))
  if (length(nums) != 2 || any(is.na(nums))) return(list(ref=NA_real_, alt=NA_real_, ok=FALSE))
  list(ref=nums[1], alt=nums[2], ok=TRUE)
}

# Read a GFF3 file and keep a minimal set of fields for genomic annotation.
read_gff3 <- function(path) {
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  cmd <- if (is_gz) {
    sprintf("zcat %s | grep -v '^#' | grep -v '^$'", shQuote(path))
  } else {
    sprintf("grep -v '^#' %s | grep -v '^$'", shQuote(path))
  }

  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 9) stop("GFF3: less than 9 columns, unexpected format.")

  setnames(dt, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]

  # Extract a simple ID from the attribute column.
  dt[, ID := sapply(attributes, function(x) {
    if (grepl("ID=", x)) {
      id_part <- sub(".*ID=([^;]+).*", "\\1", x)
      return(id_part)
    }
    return(NA_character_)
  })]

  dt
}

# Build a standardized variant table from one VCF.
vcf_to_table <- function(path, label) {
  dt <- read_vcf_body(path)
  dt[, DP := suppressWarnings(as.numeric(extract_format_field(FORMAT, SAMPLE, "DP")))]
  dt[, AD := extract_format_field(FORMAT, SAMPLE, "AD")]

  dt[, N_ALT := lengths(strsplit(ALT, ",", fixed = TRUE))]
  ad_parsed <- lapply(dt$AD, parse_biallelic_ad)

  dt[, REF_COUNT := vapply(ad_parsed, function(x) x$ref, numeric(1))]
  dt[, ALT_COUNT := vapply(ad_parsed, function(x) x$alt, numeric(1))]
  dt[, BIALLELIC_OK := (N_ALT == 1) & vapply(ad_parsed, function(x) x$ok, logical(1))]

  dt[, DP_AD := REF_COUNT + ALT_COUNT]
  dt[, VAF := ifelse(!is.na(DP_AD) & DP_AD > 0, ALT_COUNT / DP_AD, NA_real_)]

  dt[, SAMPLE_LABEL := label]
  dt[, .(CHROM, POS, REF, ALT, QUAL, FILTER, DP, REF_COUNT, ALT_COUNT, DP_AD, VAF, BIALLELIC_OK, SAMPLE_LABEL)]
}

cat("Processing VCF files...\n")
p25 <- vcf_to_table(p25_path, "P25")
p27 <- vcf_to_table(p27_path, "P27")
cat(sprintf("P25: %d variants, P27: %d variants\n", nrow(p25), nrow(p27)))

cat("Reading the GFF3 file...\n")
gff_data <- read_gff3(gff3_path)
cat(sprintf("GFF3 loaded: %d features found\n", nrow(gff_data)))

# Merge matching variants between P25 and P27.
cat("Merging P25 and P27 data...\n")
m <- merge(
  p25, p27,
  by = c("CHROM","POS","REF","ALT"),
  suffixes = c("_P25","_P27"),
  all = FALSE
)
cat(sprintf("After merge: %d shared variants\n", nrow(m)))

# Safety check: verify expected columns are present.
cat("Checking required columns...\n")
required_cols <- c("REF_COUNT_P25", "ALT_COUNT_P25", "REF_COUNT_P27", "ALT_COUNT_P27",
                   "DP_AD_P25", "DP_AD_P27", "VAF_P25", "VAF_P27", "QUAL_P25", "QUAL_P27",
                   "BIALLELIC_OK_P25", "BIALLELIC_OK_P27")
missing_cols <- setdiff(required_cols, names(m))
if (length(missing_cols) > 0) {
  cat("ERROR: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(names(m), collapse = ", "), "\n")
  quit(status = 1)
}

# Keep only trustworthy biallelic variants.
m_orig <- nrow(m)
m <- m[BIALLELIC_OK_P25 == TRUE & BIALLELIC_OK_P27 == TRUE]
cat(sprintf("Biallelic filtering: %d -> %d variants\n", m_orig, nrow(m)))
if (nrow(m) == 0) {
  cat("ERROR: No reliable biallelic variants left\n")
  quit(status = 1)
}

# -----------------------------------------------------------------------------
# Depth normalization
# -----------------------------------------------------------------------------
# The idea is to reduce the effect of different global coverage levels between P25 and P27.
median_depth_P25 <- median(m$DP_AD_P25, na.rm = TRUE)
median_depth_P27 <- median(m$DP_AD_P27, na.rm = TRUE)
depth_ratio <- median_depth_P27 / median_depth_P25
cat(sprintf("Normalization: median depth P25=%.1f, P27=%.1f, ratio=%.3f\n",
            median_depth_P25, median_depth_P27, depth_ratio))

# Normalize P25 counts so that both samples are more comparable.
m[, ALT_COUNT_P25_norm := ALT_COUNT_P25 * depth_ratio]
m[, REF_COUNT_P25_norm := REF_COUNT_P25 * depth_ratio]

# -----------------------------------------------------------------------------
# Genomic annotation with GFF3 features
# -----------------------------------------------------------------------------
cat("Annotating variants with ORFs...\n")
orfs <- gff_data[feature %in% c("gene", "CDS", "mRNA", "ORF", "coding_sequence")]

m[, `:=`(
  ORF_ID = "intergenic",
  ORF_Name = "intergenic",
  ORF_product = "intergenic",
  region_type = "Intergenic"
)]

if (nrow(orfs) > 0) {
  annotated_count <- 0
  for (i in 1:nrow(m)) {
    chrom_i <- m$CHROM[i]
    pos_i <- m$POS[i]

    overlapping <- orfs[seqname == chrom_i & start <= pos_i & end >= pos_i]

    if (nrow(overlapping) > 0) {
      best_orf <- overlapping[1]
      m[i, `:=`(
        ORF_ID = ifelse(is.na(best_orf$ID), "unknown", best_orf$ID),
        ORF_Name = ifelse(is.na(best_orf$ID), "unknown", best_orf$ID),
        ORF_product = "coding_region",
        region_type = "Genic"
      )]
      annotated_count <- annotated_count + 1
    }

    if (i %% 1000 == 0) {
      cat(sprintf("  Annotation progress: %d/%d variants\n", i, nrow(m)))
    }
  }
  cat(sprintf("Annotation finished: %d genic variants, %d intergenic variants\n",
              annotated_count, nrow(m) - annotated_count))
} else {
  cat("No ORFs found in the GFF3 file\n")
}

# -----------------------------------------------------------------------------
# Statistics: Fisher test and log2 odds ratio
# -----------------------------------------------------------------------------
cat("Computing statistics...\n")
m[, `:=`(
  a = ALT_COUNT_P27,
  b = REF_COUNT_P27,
  c = ALT_COUNT_P25_norm,
  d = REF_COUNT_P25_norm
)]

m[, log2OR := log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))]

cat("Running Fisher tests...\n")
m[, pvalue := NA_real_]
for (i in 1:nrow(m)) {
  mat <- matrix(c(m$a[i], m$b[i], m$c[i], m$d[i]), nrow = 2, byrow = TRUE)
  tryCatch({
    test_result <- stats::fisher.test(mat)
    m[i, pvalue := test_result$p.value]
  }, error = function(e) {
    # Conservative fallback if the test fails.
    m[i, pvalue := 1.0]
  })

  if (i %% 1000 == 0) {
    cat(sprintf("  Fisher progress: %d/%d variants\n", i, nrow(m)))
  }
}

m[, padj := p.adjust(pvalue, method = "BH")]

# -----------------------------------------------------------------------------
# Filtering and classification
# -----------------------------------------------------------------------------
m[, PASS_PROOF := (
  DP_AD_P25 >= minDP & DP_AD_P27 >= minDP &
  ALT_COUNT_P25 >= minALT & ALT_COUNT_P27 >= minALT &
  QUAL_P25 >= minQUAL & QUAL_P27 >= minQUAL &
  VAF_P25 >= minVAF & VAF_P27 >= minVAF
)]

m[, is_sig := (padj <= padj_th)]

m[, class := fifelse(
  padj <= padj_th & log2OR > log2FC_threshold, "Enriched in P27 (sig)",
  fifelse(padj <= padj_th & log2OR < -log2FC_threshold, "Enriched in P25 (sig)", "Not significant")
)]

cat("Saving result tables...\n")
fwrite(m, file.path(outdir, "volcano_table_improved.tsv"), sep = "\t")

# Export only significant genic variants as a convenient summary.
significant_orfs <- m[
  is_sig == TRUE & region_type == "Genic",
  .(CHROM, POS, REF, ALT, ORF_ID, ORF_Name, ORF_product,
    log2OR, pvalue, padj, class, VAF_P25, VAF_P27,
    ALT_COUNT_P25, ALT_COUNT_P27, DP_AD_P25, DP_AD_P27)
]

if (nrow(significant_orfs) > 0) {
  significant_orfs <- significant_orfs[order(padj)]
  fwrite(significant_orfs, file.path(outdir, "significant_orfs.tsv"), sep = "\t")
  cat(sprintf("Significant ORFs exported: %s\n", file.path(outdir, "significant_orfs.tsv")))
} else {
  cat("No significant ORFs found\n")
}

# -----------------------------------------------------------------------------
# Volcano plot with ORF labels
# -----------------------------------------------------------------------------
cat("Generating volcano plot with ORF annotations...\n")

m[, should_label := (padj <= padj_th & abs(log2OR) >= log2FC_threshold & region_type == "Genic")]
m[, plot_label := fifelse(
  should_label,
  ifelse(ORF_ID != "unknown" & ORF_ID != "intergenic" & !is.na(ORF_ID), ORF_ID, ""),
  ""
)]

p_vol <- ggplot(m, aes(x = log2OR, y = -log10(padj))) +
  geom_point(aes(color = class, alpha = PASS_PROOF, shape = region_type), size = 1.3) +
  geom_hline(yintercept = -log10(padj_th), linetype = "dashed") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dotted") +
  scale_shape_manual(values = c("Genic" = 16, "Intergenic" = 17)) +
  theme_bw() +
  labs(
    color = "Class",
    shape = "Region",
    title = "Volcano P27 vs P25 (log2OR vs -log10(padj))",
    subtitle = sprintf("Labelled ORFs pass padj≤%.3f and |log2FC|≥%.1f", padj_th, log2FC_threshold),
    x = "log2(odds ratio) (P27 vs P25)",
    y = "-log10(padj)"
  )

variants_to_label <- m[should_label == TRUE & plot_label != ""]
if (nrow(variants_to_label) > 0) {
  cat(sprintf("Annotating %d ORFs above significance thresholds\n", nrow(variants_to_label)))

  tryCatch({
    library(ggrepel)
    p_vol <- p_vol +
      geom_text_repel(
        data = variants_to_label,
        aes(label = plot_label),
        size = 3,
        max.overlaps = 30,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.3,
        min.segment.length = 0.1,
        force = 2,
        fontface = "bold"
      )
  }, error = function(e) {
    cat("ggrepel not available, using plain geom_text instead\n")
    p_vol <<- p_vol +
      geom_text(
        data = variants_to_label,
        aes(label = plot_label),
        size = 3,
        nudge_y = 0.2,
        check_overlap = TRUE,
        fontface = "bold"
      )
  })
} else {
  cat("No ORFs passed the significance thresholds for labelling\n")
}

ggsave(file.path(outdir, "volcano_improved.png"), p_vol, width = 12, height = 8, dpi = 300)

# Detailed plot with more explicit labels.
if (nrow(variants_to_label) > 0) {
  variants_to_label[, detailed_label := sprintf("%s\n(FC=%.1f, p=%.2e)",
                                                ORF_ID,
                                                2^abs(log2OR),
                                                padj)]

  p_vol_detailed <- ggplot(m, aes(x = log2OR, y = -log10(padj))) +
    geom_point(aes(color = class, alpha = PASS_PROOF, shape = region_type), size = 1.3) +
    geom_hline(yintercept = -log10(padj_th), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dotted", color = "blue") +
    scale_shape_manual(values = c("Genic" = 16, "Intergenic" = 17)) +
    geom_point(data = variants_to_label, aes(x = log2OR, y = -log10(padj)),
               color = "red", size = 4, shape = 1, stroke = 2) +
    theme_bw() +
    labs(
      color = "Class",
      shape = "Region",
      title = "Volcano P27 vs P25 — detailed significant ORFs",
      subtitle = sprintf("Only ORFs passing padj≤%.3f and |log2FC|≥%.1f are labelled", padj_th, log2FC_threshold),
      x = "log2(odds ratio) (P27 vs P25)",
      y = "-log10(padj)"
    )

  tryCatch({
    library(ggrepel)
    p_vol_detailed <- p_vol_detailed +
      geom_text_repel(
        data = variants_to_label,
        aes(label = detailed_label),
        size = 2.5,
        max.overlaps = 20,
        box.padding = 0.8,
        point.padding = 0.5,
        segment.color = "darkred",
        segment.size = 0.4,
        fontface = "bold",
        color = "darkred"
      )
  }, error = function(e) {
    p_vol_detailed <- p_vol_detailed +
      geom_text(
        data = variants_to_label,
        aes(label = detailed_label),
        size = 2.5,
        check_overlap = TRUE,
        fontface = "bold",
        color = "darkred"
      )
  })

  ggsave(file.path(outdir, "volcano_significant_detailed.png"), p_vol_detailed, width = 14, height = 10, dpi = 300)
  cat("Detailed volcano plot generated\n")
}

# -----------------------------------------------------------------------------
# Final console summary
# -----------------------------------------------------------------------------
cat("\n=== ANALYSIS SUMMARY ===\n")
cat(sprintf("Variants analysed: %d\n", nrow(m)))
cat(sprintf("Genic variants: %d (%.1f%%)\n",
            sum(m$region_type == "Genic"),
            100 * sum(m$region_type == "Genic") / nrow(m)))
cat(sprintf("Intergenic variants: %d (%.1f%%)\n",
            sum(m$region_type == "Intergenic"),
            100 * sum(m$region_type == "Intergenic") / nrow(m)))
cat(sprintf("Significant variants: %d (%.1f%%)\n",
            sum(m$is_sig, na.rm = TRUE),
            100 * sum(m$is_sig, na.rm = TRUE) / nrow(m)))
cat(sprintf("Enriched in P27: %d, enriched in P25: %d\n",
            sum(m$class == "Enriched in P27 (sig)", na.rm = TRUE),
            sum(m$class == "Enriched in P25 (sig)", na.rm = TRUE)))

cat("\n=== GENERATED FILES ===\n")
cat("- volcano_improved.png : volcano plot with labelled significant ORFs\n")
if (file.exists(file.path(outdir, "volcano_significant_detailed.png"))) {
  cat("- volcano_significant_detailed.png : detailed volcano plot with fold-changes\n")
}
cat("- volcano_table_improved.tsv : full result table\n")
if (file.exists(file.path(outdir, "significant_orfs.tsv"))) {
  cat("- significant_orfs.tsv : significant ORFs only\n")
}

cat(sprintf("\n[OK] Analysis finished - results in: %s\n", normalizePath(outdir)))
