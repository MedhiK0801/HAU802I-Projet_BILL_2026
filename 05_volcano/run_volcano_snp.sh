#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_volcano_snp.sh
#
# Run the volcano-plot workflow on depth-enriched SNP VCF files.
# The script filters each P25/P27 pair, runs the per-pair volcano analysis, and
# finally generates the annotated aggregated volcano plots.
#
# Usage:
#   bash scripts/run_volcano_snp.sh strict
#   bash scripts/run_volcano_snp.sh moderate
#   bash scripts/run_volcano_snp.sh both
#   bash scripts/run_volcano_snp.sh custom "FILTER"
# -----------------------------------------------------------------------------

MODE="${1:-strict}"
CUSTOM_FILTER="${2:-}"

R_SCRIPT="scripts/volcano_snp.R"
GFF3_FILE="References/DQ657948.1.gff3"
DEPTH_DIR="results/depth"
AGG_SCRIPT=""
AGG_ANN_SCRIPT="scripts/plot_volcano_aggregated_annotated.R"

run_one_mode() {
  local mode="$1"
  local filter out_root min_dp min_alt min_qual min_vaf padj threshold

  case "$mode" in
    strict)
      out_root="results/volcano_strict"
      filter='QUAL>=20 & VAF>=0.15 & AD[0:1]>=50 & DP>=200'
      min_dp=200
      min_alt=50
      min_qual=20
      min_vaf=0.15
      padj=0.05
      threshold=1
      ;;
    moderate|modere)
      out_root="results/volcano_moderate"
      filter='QUAL>=15 & VAF>=0.10 & AD[0:1]>=20 & DP>=50'
      min_dp=50
      min_alt=20
      min_qual=15
      min_vaf=0.10
      padj=0.05
      threshold=1
      ;;
    custom)
      if [[ -z "$CUSTOM_FILTER" ]]; then
        echo "ERREUR: fournir un filtre custom en 2e argument." >&2
        exit 1
      fi
      out_root="results/volcano_custom"
      filter="$CUSTOM_FILTER"
      min_dp="${MIN_DP:-50}"
      min_alt="${MIN_ALT:-20}"
      min_qual="${MIN_QUAL:-15}"
      min_vaf="${MIN_VAF:-0.10}"
      padj="${PADJ:-0.05}"
      threshold="${THRESHOLD:-1}"
      ;;
    *)
      echo "Mode non reconnu: $mode" >&2
      exit 1
      ;;
  esac

  for f in "$R_SCRIPT" "$GFF3_FILE" "$AGG_ANN_SCRIPT"; do
    [[ -f "$f" ]] || { echo "ERREUR: fichier manquant: $f" >&2; exit 1; }
  done

  local p25_dir="$DEPTH_DIR/P25"
  local p27_dir="$DEPTH_DIR/P27"
  [[ -d "$p25_dir" ]] || { echo "ERREUR: $p25_dir introuvable" >&2; exit 1; }
  [[ -d "$p27_dir" ]] || { echo "ERREUR: $p27_dir introuvable" >&2; exit 1; }

  mkdir -p "$out_root/logs"
  shopt -s nullglob

  echo "==============================================="
  echo " Volcano mode      : $mode"
  echo " Output directory  : $out_root"
  echo " bcftools filter   : $filter"
  echo " R thresholds      : MIN_DP=$min_dp MIN_ALT=$min_alt MIN_QUAL=$min_qual PADJ=$padj THRESHOLD=$threshold MIN_VAF=$min_vaf"
  echo "==============================================="

  local p25_list=("$p25_dir"/*.snp.depth.vcf)
  if (( ${#p25_list[@]} == 0 )); then
    p25_list=("$p25_dir"/*.snp.vcf.gz)
  fi

  if (( ${#p25_list[@]} == 0 )); then
    echo "Aucun VCF trouvé dans $p25_dir" >&2
    exit 1
  fi

  for p25 in "${p25_list[@]}"; do
    local base p27_base p27 sample outdir logfile p25_filtered p27_filtered n_p25 n_p27

    base="$(basename "$p25")"
    p27_base="${base//P25/P27}"
    p27="$p27_dir/$p27_base"

    if [[ ! -f "$p27" ]]; then
      echo "[SKIP] Pas de paire P27 pour: $base" >&2
      continue
    fi

    sample="${base%.vcf.gz}"
    sample="${sample%.snp.depth.vcf}"

    outdir="$out_root/$sample"
    logfile="$out_root/logs/$sample.log"
    mkdir -p "$outdir"

    p25_filtered="$outdir/${sample}.P25.filtered.vcf"
    p27_filtered="$outdir/${sample}.P27.filtered.vcf"

    bcftools view -i "$filter" "$p25" -o "$p25_filtered"
    bcftools view -i "$filter" "$p27" -o "$p27_filtered"

    n_p25=$(grep -c -v '^#' "$p25_filtered" || true)
    n_p27=$(grep -c -v '^#' "$p27_filtered" || true)

    echo "[RUN] $sample"
    echo "  P25 filtered variants: $n_p25"
    echo "  P27 filtered variants: $n_p27"

    if (( n_p25 == 0 || n_p27 == 0 )); then
      echo "  [SKIP] Pas assez de variants après filtrage" >&2
      continue
    fi

    Rscript "$R_SCRIPT" \
      "$p25_filtered" "$p27_filtered" "$GFF3_FILE" "$outdir" \
      "$min_dp" "$min_alt" "$min_qual" "$padj" "$threshold" "$min_vaf" \
      2>&1 | tee "$logfile"
  done

  echo ""
  echo "[AGG] Annotated + density aggregated volcano plots"
  Rscript "$AGG_ANN_SCRIPT" "$out_root" "$out_root/aggregated_annotated" 2 "$mode"

  echo ""
  echo "[DONE] Results written to $out_root/"
}

case "$MODE" in
  both)
    run_one_mode strict
    echo ""
    run_one_mode moderate
    ;;
  strict|moderate|modere|custom)
    run_one_mode "$MODE"
    ;;
  *)
    echo "Usage: bash scripts/run_volcano_snp.sh {strict|moderate|both|custom \"FILTER\"}" >&2
    exit 1
    ;;
esac
