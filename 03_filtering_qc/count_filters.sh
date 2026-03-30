#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# EDUCATIONAL COMMENTARY
# This shell script compares several variant-filtering strategies on all samples.
# It counts how many variants remain after each filter and prints a comparison
# table, while also saving the output into a log file.
#
# This helps to choose thresholds that remove noise without discarding too
# much biological signal.
#
# Usage (from BILL project root BILL/) :
#   bash scripts/count_filters.sh [results/depth]
# -----------------------------------------------------------------------------

set -euo pipefail

DEPTH_DIR="${1:-results/depth}"
LOGDIR="results/logs"
mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/count_filters_$(date +%Y%m%d_%H%M%S).txt"

if [[ ! -d "$DEPTH_DIR" ]]; then
  echo "ERREUR: Répertoire '$DEPTH_DIR' non trouvé" >&2
  exit 1
fi

# Tout ce qui suit est affiché à l'écran ET écrit dans le fichier log
exec > >(tee "$LOGFILE") 2>&1

echo "=== count_filters.sh — $(date) ==="
echo "Répertoire analysé : $DEPTH_DIR"
echo "Log sauvegardé dans : $LOGFILE"

shopt -s nullglob

# Filters definition
# Format : "NAME|BCFTOOLS_EXPRESSION"
FILTERS=(
  "Aucun filtre|."
  "Base (Q10,VAF0.1,AD10,DP30)|QUAL>=10 & VAF>=0.1 & AD[0:1]>=10 & DP>=30"
  "Modéré (Q15,VAF0.1,AD20,DP50)|QUAL>=15 & VAF>=0.1 & AD[0:1]>=20 & DP>=50"
  "Strict (Q20,VAF0.15,AD50,DP200)|QUAL>=20 & VAF>=0.15 & AD[0:1]>=50 & DP>=200"
)

# Fonction pour compter les variants
count_variants() {
  local vcf="$1"
  local filter="$2"
  local n
  if [[ "$filter" == "." ]]; then
    n=$(bcftools view "$vcf" 2>/dev/null | grep -c -v '^#') || true
  else
    n=$(bcftools view -i "$filter" "$vcf" 2>/dev/null | grep -c -v '^#') || true
  fi
  echo "${n:-0}"
}

# Header du tableau
printf "\n%-45s" "Échantillon"
for f in "${FILTERS[@]}"; do
  name="${f%%|*}"
  printf "%15s" "$name"
done
printf "%15s\n" "DP médian"

# Ligne de séparation
total_width=$((45 + 15 * (${#FILTERS[@]} + 1)))
printf '%*s\n' "$total_width" '' | tr ' ' '-'

for passage_dir in "$DEPTH_DIR"/P*/; do
  passage="$(basename "$passage_dir")"

  depth_files=("${passage_dir}"*.snp.depth.vcf)

  if (( ${#depth_files[@]} == 0 )); then
    continue
  fi

  # Normaliser d'abord si pas encore fait
  for file in "${depth_files[@]}"; do
    norm_file="${file%.vcf}.norm.vcf"
    if [[ ! -f "$norm_file" ]]; then
      bcftools norm -m-any "$file" -o "$norm_file" 2>/dev/null
    fi
  done

  for nfile in "${passage_dir}"*.norm.vcf; do
    sample_name="${passage}/$(basename "$nfile" .snp.depth.norm.vcf)"

    printf "%-45s" "$sample_name"

    for f in "${FILTERS[@]}"; do
      filter_expr="${f##*|}"
      count=$(count_variants "$nfile" "$filter_expr")
      printf "%15s" "$count"
    done

    # DP médian
    median_dp=$(bcftools query -f '[ %DP]\n' "$nfile" 2>/dev/null \
      | sort -n \
      | awk '{a[NR]=$1} END{if(NR%2==1) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')
    printf "%15s" "${median_dp:-N/A}"

    printf "\n"
  done
done

printf '%*s\n' "$total_width" '' | tr ' ' '-'

echo ""
echo "=== GUIDE D'INTERPRÉTATION ==="
echo ""
echo "  Base    : nettoyage minimal, garde tout ce qui est plausible"
echo "            (QUAL>=10, VAF>=0.1, AD_alt>=10, DP>=30)"
echo ""
echo "  Modéré  : bonne confiance, recommandé pour l'analyse principale"
echo "            (QUAL>=15, VAF>=0.1, AD_alt>=20, DP>=50)"
echo ""
echo "  Strict  : haute confiance, pour figures et conclusions"
echo "            (QUAL>=20, VAF>=0.15, AD_alt>=50, DP>=200)"
echo ""
echo "  Attention : si un échantillon perd >80% de ses variants entre"
echo "  'Base' et 'Modéré', sa profondeur est probablement trop faible."
echo "  Vérifiez le DP médian de cet échantillon."
echo ""
echo "=== Log sauvegardé dans : $LOGFILE ==="
