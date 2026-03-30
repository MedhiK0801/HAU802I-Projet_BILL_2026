#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# extract_isec_summary.sh
#
# Convert raw `bcftools isec` result directories into compact TSV summary files
# that can be plotted and interpreted more easily.
#
# Usage:
#   bash scripts/extract_isec_summary.sh strict
#   bash scripts/extract_isec_summary.sh moderate
#   bash scripts/extract_isec_summary.sh both
# -----------------------------------------------------------------------------

set -euo pipefail

# =============================================================================
# extract_isec_summary.sh — Extrait les données des résultats isec pour figures
#
# Usage (from BILL project root BILL/) :
#   bash scripts/extract_isec_summary.sh strict
#   bash scripts/extract_isec_summary.sh moderate
#   bash scripts/extract_isec_summary.sh both        # strict & moderate
#
# Produce TSV files in results/isec_<mode>/summary/ :
#   - counts.tsv         : variants by category and impact
#   - genes_chaud.tsv    : genes affected by heat shock
#   - genes_froid.tsv    : genes affected by cold shock
#   - missense_detail.tsv: details of missense variants
#   - individual.tsv     : results by singles pair
# =============================================================================

SNPEFF_DIR="$HOME/snpEff"

extract_mode() {
  local mode="$1"
  local ROOT="results/isec_${mode}"
  local ANN_DIR="$ROOT/annotated"
  local INDIV_DIR="$ROOT/individual"
  local SUMMARY="$ROOT/summary"

  if [[ ! -d "$ROOT" ]]; then
    echo "ERREUR: $ROOT non trouvé. Lancez d'abord run_isec_analysis.sh $mode" >&2
    return 1
  fi

  mkdir -p "$SUMMARY"
  echo "=== Extraction des données [$mode] ==="

  # -----------------------------------------------------------------
  # 1. Comptages par catégorie et impact
  # -----------------------------------------------------------------
  echo "filter	category	impact	count" > "$SUMMARY/counts.tsv"

  for group in specific_froid specific_chaud shared_froid_chaud; do
    vcf="$ANN_DIR/${group}.ann.vcf"
    if [[ ! -f "$vcf" ]]; then continue; fi

    # Compter par impact
    java -jar "$SNPEFF_DIR/SnpSift.jar" extractFields "$vcf" \
      "ANN[0].IMPACT" 2>/dev/null | tail -n+2 | sort | uniq -c | \
    while read -r count impact; do
      echo "${mode}	${group}	${impact}	${count}"
    done >> "$SUMMARY/counts.tsv"
  done
  echo "  -> counts.tsv"

  # -----------------------------------------------------------------
  # 2. Gènes touchés par groupe
  # -----------------------------------------------------------------
  for group in specific_chaud specific_froid shared_froid_chaud; do
    vcf="$ANN_DIR/${group}.ann.vcf"
    outfile="$SUMMARY/genes_${group}.tsv"

    if [[ ! -f "$vcf" ]]; then continue; fi

    echo "filter	gene	effect	impact	pos	ref	alt" > "$outfile"

    java -jar "$SNPEFF_DIR/SnpSift.jar" extractFields "$vcf" \
      CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" \
      2>/dev/null | tail -n+2 | \
    while IFS=$'\t' read -r chrom pos ref alt gene effect impact; do
      echo "${mode}	${gene}	${effect}	${impact}	${pos}	${ref}	${alt}"
    done >> "$outfile"

    echo "  -> genes_${group}.tsv"
  done

  # -----------------------------------------------------------------
  # 3. Détail des variants missense
  # -----------------------------------------------------------------
  echo "filter	category	pos	ref	alt	gene	effect	impact" > "$SUMMARY/missense_all.tsv"

  for group in specific_froid specific_chaud shared_froid_chaud; do
    vcf="$ANN_DIR/${group}.ann.vcf"
    if [[ ! -f "$vcf" ]]; then continue; fi

    java -jar "$SNPEFF_DIR/SnpSift.jar" extractFields "$vcf" \
      CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" \
      2>/dev/null | tail -n+2 | \
    while IFS=$'\t' read -r chrom pos ref alt gene effect impact hgvs; do
      if [[ "$impact" == "MODERATE" || "$impact" == "HIGH" ]]; then
        echo "${mode}	${group}	${pos}	${ref}	${alt}	${gene}	${effect}	${hgvs}"
      fi
    done >> "$SUMMARY/missense_all.tsv"
  done
  echo "  -> missense_all.tsv"

  # -----------------------------------------------------------------
  # 4. Résultats individuels (par paire P25-X vs P27-X)
  # -----------------------------------------------------------------
  echo "filter	sample	group	private_P25	private_P27	shared" > "$SUMMARY/individual.tsv"

  for indiv_dir in "$INDIV_DIR"/*/; do
    if [[ ! -d "$indiv_dir" ]]; then continue; fi
    base="$(basename "$indiv_dir")"
    num="${base##*-}"
    num="${num%%.*}"

    # Déterminer le groupe selon le mode
    local group_label
    case "$mode" in
      strict)
        case "$num" in 2|3|4) group_label="froid" ;; 7|8|9|10) group_label="chaud" ;; *) group_label="exclu" ;; esac
        ;;
      moderate|modere)
        case "$num" in 1|2|3|4) group_label="froid" ;; 7|8|9|10) group_label="chaud" ;; *) group_label="exclu" ;; esac
        ;;
    esac

    n_priv25=$(grep -c -v '^#' "$indiv_dir/private_P25.vcf" 2>/dev/null) || true
    n_priv27=$(grep -c -v '^#' "$indiv_dir/private_P27.vcf" 2>/dev/null) || true
    n_shared=$(grep -c -v '^#' "$indiv_dir/shared_P25.vcf" 2>/dev/null) || true

    echo "${mode}	${base}	${group_label}	${n_priv25}	${n_priv27}	${n_shared}"
  done >> "$SUMMARY/individual.tsv"
  echo "  -> individual.tsv"

  echo "  Résumés dans $SUMMARY/"
}

# -----------------------------------------------------------------
# Main
# -----------------------------------------------------------------
MODE="${1:-both}"

case "$MODE" in
  strict)   extract_mode "strict" ;;
  moderate) extract_mode "moderate" ;;
  both)
    extract_mode "strict"
    echo ""
    extract_mode "moderate"
    ;;
  *) echo "Usage: $0 {strict|moderate|both}" >&2; exit 1 ;;
esac

echo ""
echo "[DONE]"
