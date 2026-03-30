#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_plot_depth_all.sh
#
# Apply a strict SNP filter to every depth-enriched VCF file and generate the
# sample-level support plots used to inspect DP, AD, QUAL, and VAF.
#
# Usage:
#   bash scripts/run_plot_depth_all.sh [results/depth] [results/plots]
# -----------------------------------------------------------------------------

set -euo pipefail

# =============================================================================
# run_plot_depth_all.sh — Adapté pour l'arborescence results/depth/P*/
#
# Usage (depuis la racine du projet BILL/) :
#   bash scripts/run_plot_depth_all.sh [results/depth] [results/plots]
#
# Applique le filtre strict, puis génère les plots pour chaque échantillon.
# =============================================================================

DEPTH_DIR="${1:-results/depth}"
OUTPUT_DIR="${2:-results/plots}"
PY_SCRIPT="scripts/plot_variant_support.py"

# Filtre strict
FILTER="QUAL>=20 & VAF>=0.15 & AD[0:1]>=50 & DP>=200"

if [[ ! -d "$DEPTH_DIR" ]]; then
  echo "ERREUR: Répertoire d'entrée '$DEPTH_DIR' non trouvé" >&2
  exit 1
fi
if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "ERREUR: Script Python '$PY_SCRIPT' non trouvé" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob

echo "Filtre appliqué : $FILTER"
echo ""

for passage_dir in "$DEPTH_DIR"/P*/; do
  passage="$(basename "$passage_dir")"
  echo "=== Passage : $passage ==="

  vcf_files=("${passage_dir}"*.snp.depth.vcf)

  if (( ${#vcf_files[@]} == 0 )); then
    echo "  [SKIP] Aucun .snp.depth.vcf dans $passage_dir"
    continue
  fi

  # Créer un sous-dossier de sortie par passage
  passage_outdir="$OUTPUT_DIR/$passage"
  mkdir -p "$passage_outdir"

  for vcf in "${vcf_files[@]}"; do
    base="$(basename "$vcf" .snp.depth.vcf)"
    filtered_vcf="${passage_dir}${base}.snp.depth.filtered.vcf"

    # Appliquer le filtre strict
    bcftools view -i "$FILTER" "$vcf" -o "$filtered_vcf"

    n_before=$(grep -c -v '^#' "$vcf") || true
    n_after=$(grep -c -v '^#' "$filtered_vcf") || true
    echo "  $base : $n_before -> $n_after variants (filtre strict)"

    if (( n_after == 0 )); then
      echo "  [SKIP] Aucun variant après filtrage pour $base"
      rm -f "$filtered_vcf"
      continue
    fi

    python3 "$PY_SCRIPT" "$filtered_vcf" "$passage_outdir/"
  done
done

echo ""
echo "=== Tous les fichiers traités ==="
echo "Plots dans $OUTPUT_DIR/<passage>/"
