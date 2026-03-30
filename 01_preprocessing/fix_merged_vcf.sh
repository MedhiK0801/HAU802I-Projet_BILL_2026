#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# fix_merged_vcf.sh
#
# Detect VCF files that contain more than one sample column and collapse them
# into a single representative sample using `fix_merged_vcf.py`.
#
# Usage:
#   bash scripts/fix_merged_vcf.sh [inputs_all]
# -----------------------------------------------------------------------------

set -euo pipefail

# =============================================================================
# fix_merged_vcf.sh — Détecte et corrige les VCF multi-sample
#
# Usage (depuis la racine du projet BILL/) :
#   bash scripts/fix_merged_vcf.sh [inputs_all]
#
# Pour chaque .snp.vcf dans inputs_all/P*/, vérifie s'il est multi-sample.
# Si oui, fusionne en single-sample et remplace le fichier original
# (l'original est sauvegardé en .snp.vcf.bak).
# =============================================================================

INPUT_DIR="${1:-inputs_all}"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERREUR: Répertoire '$INPUT_DIR' non trouvé" >&2
  exit 1
fi

FIX_SCRIPT="scripts/fix_merged_vcf.py"
if [[ ! -f "$FIX_SCRIPT" ]]; then
  echo "ERREUR: Script '$FIX_SCRIPT' non trouvé" >&2
  exit 1
fi

shopt -s nullglob

fixed=0
skipped=0

for passage_dir in "$INPUT_DIR"/P*/; do
  passage="$(basename "$passage_dir")"

  for vcf in "${passage_dir}"*.snp.vcf; do
    n_samples=$(bcftools query -l "$vcf" | wc -l)

    if (( n_samples > 1 )); then
      base="$(basename "$vcf")"
      echo "[FIX] $passage/$base : $n_samples samples -> fusion en single-sample"

      # Sauvegarder l'original
      cp "$vcf" "${vcf}.bak"

      # Fixer
      python3 "$FIX_SCRIPT" "$vcf" "${vcf}.tmp"
      mv "${vcf}.tmp" "$vcf"

      fixed=$((fixed + 1))
    else
      skipped=$((skipped + 1))
    fi
  done
done

echo ""
echo "=== TERMINÉ ==="
echo "  Fichiers corrigés : $fixed"
echo "  Fichiers déjà OK  : $skipped"
if (( fixed > 0 )); then
  echo "  (les originaux sont sauvegardés en .snp.vcf.bak)"
fi
