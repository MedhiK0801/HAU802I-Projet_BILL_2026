#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_isec_analysis.sh
#
# Run the isec-based SNP comparison workflow. The script applies a selected
# filter, prepares bgzipped/indexed VCF files, performs pairwise P25/P27
# comparisons, groups the results by thermal condition, and annotates them.
#
# Usage:
#   bash scripts/run_isec_analysis.sh strict
#   bash scripts/run_isec_analysis.sh moderate
#   bash scripts/run_isec_analysis.sh custom "QUAL>=15 & DP>=100" "1 5" results/isec_custom
# -----------------------------------------------------------------------------

set -euo pipefail

MODE="${1:-strict}"
DEPTH_DIR="results/depth"
SNPEFF_DIR="$HOME/snpEff"
GENOME_NAME="CyHV3_KHV_U"

# =========================================================================
# PRESETS
# =========================================================================
case "$MODE" in
  strict)
    FILTER="QUAL>=20 & VAF>=0.15 & AD[0:1]>=50 & DP>=200"
    EXCLUDE_SAMPLES="1 5 6"
    FROID_SAMPLES="2 3 4"
    CHAUD_SAMPLES="7 8 9 10"
    OUT_ROOT="results/isec_strict"
    ;;
  moderate|modere)
    FILTER="QUAL>=15 & VAF>=0.1 & AD[0:1]>=20 & DP>=50"
    # Avec le filtre modéré, P25-5 (DP=41) et P25-6 (DP=29) restent
    # problématiques (DP < 50), mais P25-1 (DP=103) et P27-1 (DP=103)
    # passent maintenant le seuil DP>=50
    EXCLUDE_SAMPLES="5 6"
    FROID_SAMPLES="1 2 3 4"
    CHAUD_SAMPLES="7 8 9 10"
    OUT_ROOT="results/isec_moderate"
    ;;
  custom)
    FILTER="${2:?Usage: $0 custom \"FILTER_EXPR\" \"EXCLUDE\" OUT_DIR}"
    EXCLUDE_SAMPLES="${3:-}"
    OUT_ROOT="${4:-results/isec_custom}"
    # Groupes par défaut, modifiables
    FROID_SAMPLES="${5:-2 3 4}"
    CHAUD_SAMPLES="${6:-7 8 9 10}"
    ;;
  *)
    echo "Usage: $0 {strict|moderate|custom}" >&2
    echo "  strict   : QUAL>=20, VAF>=0.15, AD>=50, DP>=200 (exclut 1,5,6)"
    echo "  moderate : QUAL>=15, VAF>=0.1, AD>=20, DP>=50 (exclut 5,6)"
    echo "  custom   : $0 custom \"FILTER\" \"EXCLUDE\" OUT_DIR" >&2
    exit 1
    ;;
esac

# =========================================================================
# Vérifications
# =========================================================================
for passage in P25 P27; do
  if [[ ! -d "$DEPTH_DIR/$passage" ]]; then
    echo "ERREUR: $DEPTH_DIR/$passage non trouvé" >&2
    exit 1
  fi
done

if [[ ! -f "$SNPEFF_DIR/snpEff.jar" ]]; then
  echo "ERREUR: snpEff non trouvé. Lancez d'abord : bash scripts/setup_snpeff.sh" >&2
  exit 1
fi

# Nettoyer les résultats précédents
rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT/logs"

echo "============================================="
echo "  Pipeline isec + snpEff [$MODE]"
echo "============================================="
echo "Filtre           : $FILTER"
echo "Exclus           : échantillons ${EXCLUDE_SAMPLES:-aucun}"
echo "Groupe froid     : échantillons $FROID_SAMPLES"
echo "Groupe chaud     : échantillons $CHAUD_SAMPLES"
echo "Sortie           : $OUT_ROOT"
echo ""

# =========================================================================
# Fonctions utilitaires
# =========================================================================
get_sample_num() {
  local name="$1"
  local num="${name##*-}"
  num="${num%%.*}"
  echo "$num"
}

is_excluded() {
  local num="$1"
  if [[ -z "$EXCLUDE_SAMPLES" ]]; then
    return 1
  fi
  for ex in $EXCLUDE_SAMPLES; do
    if [[ "$num" == "$ex" ]]; then
      return 0
    fi
  done
  return 1
}

prepare_vcf() {
  local input_vcf="$1"
  local output_dir="$2"
  local base
  base="$(basename "$input_vcf" .snp.depth.vcf)"
  local filtered="${output_dir}/${base}.filtered.vcf.gz"

  bcftools view -i "$FILTER" "$input_vcf" | bgzip -c > "$filtered"
  tabix -p vcf "$filtered"
  echo "$filtered"
}

# =========================================================================
# ÉTAPE 1 : Préparer les VCF filtrés et indexés
# =========================================================================
echo "=== ÉTAPE 1 : Filtrage + bgzip + tabix ==="

PREP_DIR="$OUT_ROOT/prepared"
mkdir -p "$PREP_DIR"/{P15,P25,P27,P30}

for passage in P15 P25 P27 P30; do
  passage_dir="$DEPTH_DIR/$passage"
  if [[ ! -d "$passage_dir" ]]; then
    echo "  [SKIP] $passage_dir non trouvé (contrôle optionnel)"
    continue
  fi

  echo "  Passage $passage :"
  shopt -s nullglob
  for vcf in "$passage_dir"/*.snp.depth.vcf; do
    base="$(basename "$vcf" .snp.depth.vcf)"
    sample_num=$(get_sample_num "$base")

    if [[ "$passage" == "P25" || "$passage" == "P27" ]] && is_excluded "$sample_num"; then
      echo "    [EXCLU] $(basename "$vcf") (échantillon $sample_num)"
      continue
    fi

    out=$(prepare_vcf "$vcf" "$PREP_DIR/$passage")
    n=$(bcftools view -H "$out" | wc -l)
    echo "    $(basename "$vcf") -> $n variants"
  done
done

# =========================================================================
# ÉTAPE 2 : isec individu par individu (P25-X vs P27-X)
# =========================================================================
echo ""
echo "=== ÉTAPE 2 : isec individu par individu (P25 vs P27) ==="

INDIV_DIR="$OUT_ROOT/individual"
mkdir -p "$INDIV_DIR"

shopt -s nullglob
for p25_vcf in "$PREP_DIR/P25"/*.filtered.vcf.gz; do
  base="$(basename "$p25_vcf" .filtered.vcf.gz)"
  sample_num=$(get_sample_num "$base")

  p27_base="${base//P25/P27}"
  p27_vcf="$PREP_DIR/P27/${p27_base}.filtered.vcf.gz"

  if [[ ! -f "$p27_vcf" ]]; then
    echo "  [SKIP] Pas de paire P27 pour $base"
    continue
  fi

  n_p25=$(bcftools view -H "$p25_vcf" | wc -l)
  n_p27=$(bcftools view -H "$p27_vcf" | wc -l)

  if (( n_p25 == 0 || n_p27 == 0 )); then
    echo "  [SKIP] $base : P25=$n_p25, P27=$n_p27 variants (insuffisant)"
    continue
  fi

  outdir="$INDIV_DIR/$base"
  mkdir -p "$outdir"

  echo "  $base vs $p27_base :"
  bcftools isec -p "$outdir" "$p25_vcf" "$p27_vcf" 2>/dev/null

  mv "$outdir/0000.vcf" "$outdir/private_P25.vcf"
  mv "$outdir/0001.vcf" "$outdir/private_P27.vcf"
  mv "$outdir/0002.vcf" "$outdir/shared_P25.vcf"
  mv "$outdir/0003.vcf" "$outdir/shared_P27.vcf"

  n_priv25=$(grep -c -v '^#' "$outdir/private_P25.vcf") || true
  n_priv27=$(grep -c -v '^#' "$outdir/private_P27.vcf") || true
  n_shared=$(grep -c -v '^#' "$outdir/shared_P25.vcf") || true

  echo "    Privé P25 (perdus): $n_priv25 | Privé P27 (apparus): $n_priv27 | Partagés: $n_shared"

  # Afficher le groupe
  local_is_froid=false
  for fn in $FROID_SAMPLES; do
    [[ "$sample_num" == "$fn" ]] && local_is_froid=true
  done
  if $local_is_froid; then
    echo "    -> Groupe FROID"
  else
    echo "    -> Groupe CHAUD"
  fi
done

# =========================================================================
# ÉTAPE 3 : Variants récurrents par groupe
# =========================================================================
echo ""
echo "=== ÉTAPE 3 : Variants récurrents par groupe ==="

GROUP_DIR="$OUT_ROOT/groups"
mkdir -p "$GROUP_DIR"/{froid,chaud}

combine_group_variants() {
  local group_name="$1"
  shift
  local sample_nums=("$@")
  local group_outdir="$GROUP_DIR/$group_name"
  local vcf_list=()

  for num in "${sample_nums[@]}"; do
    for indiv_dir in "$INDIV_DIR"/*; do
      if [[ -d "$indiv_dir" ]]; then
        dir_base="$(basename "$indiv_dir")"
        dir_num=$(get_sample_num "$dir_base")

        if [[ "$dir_num" == "$num" && -f "$indiv_dir/private_P27.vcf" ]]; then
          local n_var
          n_var=$(grep -c -v '^#' "$indiv_dir/private_P27.vcf") || true
          if (( n_var == 0 )); then
            echo "  [SKIP] Échantillon $num : 0 variants apparus"
            continue
          fi

          local tmp_gz="${indiv_dir}/private_P27.vcf.gz"
          bgzip -c "$indiv_dir/private_P27.vcf" > "$tmp_gz"
          tabix -p vcf "$tmp_gz"
          vcf_list+=("$tmp_gz")
        fi
      fi
    done
  done

  if (( ${#vcf_list[@]} == 0 )); then
    echo "  [WARN] Aucun fichier avec variants pour le groupe $group_name"
    return 1
  fi

  echo "  Groupe $group_name : ${#vcf_list[@]} échantillons"

  bcftools concat -a -D "${vcf_list[@]}" | bcftools sort | bgzip -c > "$group_outdir/all_appeared.vcf.gz"
  tabix -p vcf "$group_outdir/all_appeared.vcf.gz"

  if (( ${#vcf_list[@]} >= 2 )); then
    bcftools isec -n+2 -o "$group_outdir/recurrent_sites.txt" -O v "${vcf_list[@]}" 2>/dev/null || true
    bcftools merge "${vcf_list[@]}" 2>/dev/null | bgzip -c > "$group_outdir/merged_appeared.vcf.gz" || true
    if [[ -f "$group_outdir/merged_appeared.vcf.gz" ]]; then
      tabix -p vcf "$group_outdir/merged_appeared.vcf.gz"
    fi
  fi

  local n_total
  n_total=$(bcftools view -H "$group_outdir/all_appeared.vcf.gz" | wc -l)
  echo "    Total variants apparus (union) : $n_total"
}

combine_group_variants "froid" $FROID_SAMPLES || true
combine_group_variants "chaud" $CHAUD_SAMPLES || true

# =========================================================================
# ÉTAPE 4 : isec froid vs chaud
# =========================================================================
echo ""
echo "=== ÉTAPE 4 : isec froid vs chaud ==="

COMPARE_DIR="$OUT_ROOT/froid_vs_chaud"
mkdir -p "$COMPARE_DIR"

froid_vcf="$GROUP_DIR/froid/all_appeared.vcf.gz"
chaud_vcf="$GROUP_DIR/chaud/all_appeared.vcf.gz"

if [[ -f "$froid_vcf" && -f "$chaud_vcf" ]]; then
  bcftools isec -p "$COMPARE_DIR" "$froid_vcf" "$chaud_vcf" 2>/dev/null

  mv "$COMPARE_DIR/0000.vcf" "$COMPARE_DIR/specific_froid.vcf"
  mv "$COMPARE_DIR/0001.vcf" "$COMPARE_DIR/specific_chaud.vcf"
  mv "$COMPARE_DIR/0002.vcf" "$COMPARE_DIR/shared_froid.vcf"
  mv "$COMPARE_DIR/0003.vcf" "$COMPARE_DIR/shared_chaud.vcf"

  n_froid=$(grep -c -v '^#' "$COMPARE_DIR/specific_froid.vcf") || true
  n_chaud=$(grep -c -v '^#' "$COMPARE_DIR/specific_chaud.vcf") || true
  n_shared=$(grep -c -v '^#' "$COMPARE_DIR/shared_froid.vcf") || true

  echo "  Spécifiques froid : $n_froid"
  echo "  Spécifiques chaud : $n_chaud"
  echo "  Partagés froid/chaud : $n_shared"
else
  echo "  [SKIP] Fichiers de groupe manquants"
fi

# =========================================================================
# ÉTAPE 5 : Annotation snpEff
# =========================================================================
echo ""
echo "=== ÉTAPE 5 : Annotation snpEff ==="

ANN_DIR="$OUT_ROOT/annotated"
mkdir -p "$ANN_DIR"

annotate_vcf() {
  local input_vcf="$1"
  local label="$2"

  if [[ ! -f "$input_vcf" ]]; then
    echo "  [SKIP] $input_vcf non trouvé"
    return
  fi

  local n
  n=$(grep -c -v '^#' "$input_vcf") || true
  if (( n == 0 )); then
    echo "  [SKIP] $label : 0 variants"
    return
  fi

  local output_vcf="$ANN_DIR/${label}.ann.vcf"
  local stats_file="$ANN_DIR/${label}.stats.html"

  java -Xmx4g -jar "$SNPEFF_DIR/snpEff.jar" \
    -noLog \
    "$GENOME_NAME" \
    -stats "$stats_file" \
    "$input_vcf" > "$output_vcf" 2>"$OUT_ROOT/logs/snpeff_${label}.log"

  local n_ann
  n_ann=$(grep -c -v '^#' "$output_vcf") || true
  echo "  $label : $n_ann variants annotés -> $output_vcf"
}

annotate_vcf "$COMPARE_DIR/specific_froid.vcf" "specific_froid"
annotate_vcf "$COMPARE_DIR/specific_chaud.vcf" "specific_chaud"
annotate_vcf "$COMPARE_DIR/shared_froid.vcf" "shared_froid_chaud"

for ctrl_passage in P15 P30; do
  for ctrl_vcf in "$PREP_DIR/$ctrl_passage"/*.filtered.vcf.gz; do
    if [[ -f "$ctrl_vcf" ]]; then
      ctrl_base="$(basename "$ctrl_vcf" .filtered.vcf.gz)"
      annotate_vcf "$ctrl_vcf" "ctrl_${ctrl_base}"
    fi
  done
done

# =========================================================================
# RÉSUMÉ FINAL
# =========================================================================
echo ""
echo "============================================="
echo "=== RÉSUMÉ [$MODE] ==="
echo "============================================="
echo ""
echo "Filtre           : $FILTER"
echo "Échantillons exclus : ${EXCLUDE_SAMPLES:-aucun}"
echo "Groupe froid     : échantillons $FROID_SAMPLES"
echo "Groupe chaud     : échantillons $CHAUD_SAMPLES"
echo ""
if [[ -f "$COMPARE_DIR/specific_froid.vcf" ]]; then
  echo "Résultats :"
  echo "  Spécifiques froid : $(grep -c -v '^#' "$COMPARE_DIR/specific_froid.vcf" || true)"
  echo "  Spécifiques chaud : $(grep -c -v '^#' "$COMPARE_DIR/specific_chaud.vcf" || true)"
  echo "  Partagés          : $(grep -c -v '^#' "$COMPARE_DIR/shared_froid.vcf" || true)"
fi
echo ""
echo "Résultats dans : $OUT_ROOT/"
echo "[DONE]"
