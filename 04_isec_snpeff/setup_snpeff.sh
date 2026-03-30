#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# setup_snpeff.sh
#
# Install snpEff if necessary and build the custom CyHV-3 annotation database
# used by the project.
#
# Usage:
#   bash scripts/setup_snpeff.sh [References/KHV-U_trunc.fasta] [References/DQ657948.1.gff3]
# -----------------------------------------------------------------------------

set -euo pipefail

# =============================================================================
# setup_snpeff.sh — Installe snpEff et construit la base CyHV-3 (KHV-U)
#
# Usage (depuis la racine du projet BILL/) :
#   bash scripts/setup_snpeff.sh [References/KHV-U_trunc.fasta] [References/DQ657948.1.gff3]
#
# Prérequis : java (>= 11)
# =============================================================================

FASTA="${1:-References/KHV-U_trunc.fasta}"
GFF3="${2:-References/DQ657948.1.gff3}"
SNPEFF_VERSION="5.2a"
SNPEFF_DIR="$HOME/snpEff"
GENOME_NAME="CyHV3_KHV_U"

# Vérifier Java
if ! command -v java &>/dev/null; then
  echo "ERREUR: Java non trouvé. Installez Java >= 11 (ex: sudo apt install default-jre)" >&2
  exit 1
fi

echo "=== Installation de snpEff ==="

if [[ -f "$SNPEFF_DIR/snpEff.jar" ]]; then
  echo "snpEff déjà installé dans $SNPEFF_DIR"
else
  echo "Téléchargement de snpEff v${SNPEFF_VERSION}..."
  cd "$HOME"
  wget -q "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" -O snpEff_latest_core.zip
  unzip -qo snpEff_latest_core.zip
  rm -f snpEff_latest_core.zip
  echo "snpEff installé dans $SNPEFF_DIR"
fi

# Vérifier les fichiers de référence
cd - > /dev/null  # Retour au répertoire initial
if [[ ! -f "$FASTA" ]]; then
  echo "ERREUR: Fichier FASTA '$FASTA' non trouvé" >&2
  exit 1
fi
if [[ ! -f "$GFF3" ]]; then
  echo "ERREUR: Fichier GFF3 '$GFF3' non trouvé" >&2
  exit 1
fi

echo ""
echo "=== Construction de la base snpEff pour $GENOME_NAME ==="

# Ajouter le génome à snpEff.config s'il n'y est pas déjà
if ! grep -q "$GENOME_NAME" "$SNPEFF_DIR/snpEff.config"; then
  echo "" >> "$SNPEFF_DIR/snpEff.config"
  echo "# CyHV-3 KHV-U (DQ657948.1)" >> "$SNPEFF_DIR/snpEff.config"
  echo "${GENOME_NAME}.genome : Cyprinid herpesvirus 3 strain KHV-U" >> "$SNPEFF_DIR/snpEff.config"
  echo "Génome ajouté à snpEff.config"
else
  echo "Génome déjà présent dans snpEff.config"
fi

# Créer le répertoire de données
DATA_DIR="$SNPEFF_DIR/data/$GENOME_NAME"
mkdir -p "$DATA_DIR"

# Copier les fichiers
cp "$(realpath "$GFF3")" "$DATA_DIR/genes.gff"
cp "$(realpath "$FASTA")" "$DATA_DIR/sequences.fa"

echo "Fichiers copiés dans $DATA_DIR/"

# Construire la base
echo "Construction de la base de données..."
cd "$SNPEFF_DIR"
java -jar snpEff.jar build -gff3 -v "$GENOME_NAME" 2>&1 | tail -5

echo ""
echo "=== TERMINÉ ==="
echo "Base snpEff '$GENOME_NAME' construite avec succès."
echo ""
echo "Pour annoter un VCF :"
echo "  java -jar $SNPEFF_DIR/snpEff.jar $GENOME_NAME input.vcf > output.ann.vcf"
echo ""
echo "Nom du génome à utiliser dans les scripts : $GENOME_NAME"
