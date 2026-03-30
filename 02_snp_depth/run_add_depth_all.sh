#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_add_depth_all.sh
#
# For every sample, run `bcftools mpileup` on the aligned BAM file and then use
# `add_variant_depth.py` to copy DP, AD, and VAF into the Medaka SNP VCF.
#
# Usage:
#   bash scripts/run_add_depth_all.sh [inputs_all] [References/KHV-U_trunc.fasta] [threads]
# -----------------------------------------------------------------------------

set -euo pipefail

# =============================================================================
# run_add_depth_all.sh — Adapté pour l'arborescence inputs_all/P*/
# 
# Usage (depuis la racine du projet BILL/) :
#   bash scripts/run_add_depth_all.sh [inputs_all] [References/KHV-U_trunc.fasta] [threads]
#
# Structure attendue :
#   BILL/
#   ├── inputs_all/
#   │   ├── P15/   (contient *.aligned.sorted.bam et *.snp.vcf)
#   │   ├── P25/
#   │   └── ...
#   ├── References/
#   │   └── KHV-U_trunc.fasta
#   ├── results/depth/
#   └── scripts/
# =============================================================================

INPUT_DIR="${1:-inputs_all}"
REF="${2:-References/KHV-U_trunc.fasta}"
THREADS="${3:-4}"
OUTDIR="results/depth"

mkdir -p "$OUTDIR"

# Vérifications
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERREUR: Répertoire d'entrée '$INPUT_DIR' non trouvé" >&2
  exit 1
fi
if [[ ! -f "$REF" ]]; then
  echo "ERREUR: Fichier de référence '$REF' non trouvé" >&2
  exit 1
fi

shopt -s nullglob

# Boucle sur chaque sous-dossier de passage (P15, P25, P27, ...)
for passage_dir in "$INPUT_DIR"/P*/; do
  passage="$(basename "$passage_dir")"
  echo ""
  echo "=============================="
  echo "=== Passage : $passage ==="
  echo "=============================="

  bams=("$passage_dir"*.aligned.sorted.bam)

  if (( ${#bams[@]} == 0 )); then
    echo "[SKIP] Aucun BAM '*.aligned.sorted.bam' trouvé dans $passage_dir" >&2
    continue
  fi

  # Créer un sous-dossier de résultats par passage
  passage_outdir="$OUTDIR/$passage"
  mkdir -p "$passage_outdir"

  for bam in "${bams[@]}"; do
    base="$(basename "$bam")"
    sample="${base%.aligned.sorted.bam}"
    snp_vcf="${passage_dir}${sample}.snp.vcf"

    # Fichiers de sortie
    mpileup_vcf="${passage_outdir}/${sample}.mpileup.vcf"
    depth_vcf="${passage_outdir}/${sample}.snp.depth.vcf"

    if [[ ! -f "$snp_vcf" ]]; then
      echo "[SKIP] Manque $snp_vcf" >&2
      continue
    fi

    echo "[${passage}/${sample}] bcftools mpileup"
    bcftools mpileup \
      --threads "$THREADS" \
      -I \
      -d 1000000 \
      -o "$mpileup_vcf" \
      -f "$REF" \
      -a AD,DP \
      -T "$snp_vcf" \
      "$bam"

    echo "[${passage}/${sample}] add_variant_depth.py -> $depth_vcf"
    python3 scripts/add_variant_depth.py "$snp_vcf" "$mpileup_vcf" "$depth_vcf"
  done
done

echo ""
echo "=== ANALYSE TERMINÉE ==="
echo "Les fichiers .depth.vcf sont dans $OUTDIR/<passage>/"
