#!/usr/bin/env python3
"""
plot_variant_support.py

Read one depth-enriched SNP VCF file and generate exploratory plots showing how
variant support behaves in terms of DP, AD, QUAL, and VAF.
"""

import os
import sys
from typing import Sized
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps


def safe_len(obj: Sized | None) -> int:
    return len(obj) if obj is not None else 0


# Check if the inputs files are provided and has a .vcf extension
if len(sys.argv) < 3 or not sys.argv[1].lower().endswith(".vcf"):
    print("\nUsage: python plot_variant_support.py <input_file.vcf> <output_directory_figures>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = os.path.basename(input_file)
output_directory = sys.argv[2]

# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)
if not output_directory.endswith('/'):
    output_directory += '/'

# We use the pysam.VariantFile to manipulate vcf input files
vcf_in = ps.VariantFile(input_file)

pos = []
ad_total = []  # Support total du variant (REF + ALT)
ad_alt = []    # Support alternatif (ALT)
ad_ref = []    # Support référence (REF)
qual = []
vaf = []
dp_total = []  # DP pour comparaison

print("Lecture du fichier VCF - Focus sur le support des variants (AD)...")

for rec in vcf_in.fetch():
    # Skip if no samples
    if len(rec.samples) == 0:
        continue
    
    sample_name = list(rec.samples.keys())[0]
    sample = rec.samples[sample_name]
    
    # Get basic info
    pos.append(rec.pos)
    qual.append(rec.qual if rec.qual is not None else 0)
    
    # Get DP for comparison
    dp = sample.get("DP", 0)
    dp_total.append(dp if dp is not None else 0)
    
    # Focus on AD (Allelic Depth) - support direct du variant
    ad = sample.get("AD", None)
    if ad is not None and len(ad) >= 2:
        ref_count = ad[0] if ad[0] is not None else 0
        alt_count = ad[1] if ad[1] is not None else 0
        total_support = ref_count + alt_count
        
        ad_ref.append(ref_count)
        ad_alt.append(alt_count)
        ad_total.append(total_support)
        
        if total_support > 0:
            variant_frequency = alt_count / total_support
        else:
            variant_frequency = 0.0
    else:
        # Skip variants without AD information
        ad_ref.append(0)
        ad_alt.append(0)
        ad_total.append(0)
        variant_frequency = 0.0
    
    vaf.append(variant_frequency)

vcf_in.close()

# Convert to numpy arrays
pos = np.array(pos)
ad_total = np.array(ad_total)
ad_alt = np.array(ad_alt) 
ad_ref = np.array(ad_ref)
qual = np.array(qual)
vaf = np.array(vaf)
dp_total = np.array(dp_total)

# Calculer le ratio AD/DP immédiatement après conversion en arrays
ad_dp_ratio = np.where(dp_total > 0, ad_total / dp_total, 0)

print(f"Traitement terminé: {len(pos)} variants analysés")
print(f"Ratio AD/DP moyen: {np.mean(ad_dp_ratio):.3f}")
print(f"Support alternatif moyen: {np.mean(ad_alt):.1f}")
print(f"Support total moyen (AD): {np.mean(ad_total):.1f}")
print(f"Profondeur totale moyenne (DP): {np.mean(dp_total):.1f}")
print(f"VAF moyenne: {np.mean(vaf):.3f}")

# Set up matplotlib style
plt.style.use('seaborn-v0_8')
plt.rcParams['figure.dpi'] = 150

# *** 1. COMPARAISON DP vs AD ***
plt.figure(figsize=(12, 8))
plt.scatter(dp_total, ad_total, c=vaf, cmap='viridis', alpha=0.6, s=30)
plt.colorbar(label='VAF')
plt.xlabel("Profondeur totale (DP)")
plt.ylabel("Support du variant (AD total)")
plt.title("Support du variant (AD) vs Profondeur totale (DP)")
plt.plot([0, max(dp_total)], [0, max(dp_total)], 'r--', alpha=0.5, label='AD = DP')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.dp_vs_ad.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 2. SUPPORT ALTERNATIF vs QUALITÉ ***
plt.figure(figsize=(12, 8))
scatter = plt.scatter(qual, ad_alt, c=vaf, cmap='plasma', alpha=0.6, s=30)
plt.colorbar(scatter, label='VAF')
plt.xlabel("Qualité Phred")
plt.ylabel("Support alternatif (AD alt)")
plt.title("Support alternatif vs Qualité (coloré par VAF)")
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.qual_vs_ad_alt.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 3. SUPPORT TOTAL vs QUALITÉ ***
plt.figure(figsize=(12, 8))
scatter = plt.scatter(qual, ad_total, c=vaf, cmap='viridis', alpha=0.6, s=30)
plt.colorbar(scatter, label='VAF')
plt.xlabel("Qualité Phred")
plt.ylabel("Support total du variant (AD total)")
plt.title("Support total vs Qualité (coloré par VAF)")
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.qual_vs_ad_total.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 4. VAF vs SUPPORT ALTERNATIF ***
plt.figure(figsize=(12, 8))
scatter = plt.scatter(vaf, ad_alt, c=qual, cmap='coolwarm', alpha=0.6, s=30)
plt.colorbar(scatter, label='Qualité Phred')
plt.xlabel("VAF (Variant Allele Frequency)")
plt.ylabel("Support alternatif (AD alt)")
plt.title("Support alternatif vs VAF (coloré par Qualité)")
plt.axvline(x=0.05, color='red', linestyle='--', alpha=0.7, label='VAF = 5%')
plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='VAF = 50%')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.vaf_vs_ad_alt.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 5. HISTOGRAMMES DES SUPPORTS ***
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

# Histogramme support alternatif
ax1.hist(ad_alt, bins=50, edgecolor="black", alpha=0.7, color='skyblue')
ax1.set_xlabel("Support alternatif (AD alt)")
ax1.set_ylabel("Nombre de variants")
ax1.set_title("Distribution du support alternatif")
ax1.axvline(x=np.median(ad_alt), color='red', linestyle='--', 
           label=f'Médiane = {np.median(ad_alt):.0f}')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Histogramme support total
ax2.hist(ad_total, bins=50, edgecolor="black", alpha=0.7, color='lightgreen')
ax2.set_xlabel("Support total (AD total)")
ax2.set_ylabel("Nombre de variants")
ax2.set_title("Distribution du support total")
ax2.axvline(x=np.median(ad_total), color='red', linestyle='--', 
           label=f'Médiane = {np.median(ad_total):.0f}')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Histogramme VAF
ax3.hist(vaf, bins=50, edgecolor="black", alpha=0.7, color='lightcoral')
ax3.set_xlabel("VAF")
ax3.set_ylabel("Nombre de variants")
ax3.set_title("Distribution de la VAF")
ax3.axvline(x=0.05, color='red', linestyle='--', alpha=0.7, label='VAF = 5%')
ax3.axvline(x=np.median(vaf), color='blue', linestyle='--', 
           label=f'Médiane = {np.median(vaf):.3f}')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Histogramme qualité
ax4.hist(qual, bins=30, edgecolor="black", alpha=0.7, color='orange')
ax4.set_xlabel("Qualité Phred")
ax4.set_ylabel("Nombre de variants")
ax4.set_title("Distribution de la qualité")
ax4.axvline(x=np.median(qual), color='red', linestyle='--', 
           label=f'Médiane = {np.median(qual):.1f}')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_directory}{output_file}.distributions.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 6. GRAPHIQUES SUPPLÉMENTAIRES DU SCRIPT PRÉCÉDENT ADAPTÉS À AD ***

# Scatter plot VAF vs AD total (au lieu de DP)
plt.figure(figsize=(10, 8))
scatter = plt.scatter(vaf, ad_total, c=qual, cmap='plasma', alpha=0.6, s=30)
plt.colorbar(scatter, label='Qualité Phred')
plt.xlabel("VAF (Variant Allele Frequency)")
plt.ylabel("Support total du variant (AD total)")
plt.title("Support total vs VAF (coloré par Qualité)")
plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='VAF = 50%')
plt.axvline(x=0.05, color='orange', linestyle='--', alpha=0.7, label='VAF = 5%')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.vaf_vs_ad_total_qual.png", dpi=300, bbox_inches='tight')
plt.close()

# Scatter plot VAF vs Qualité coloré par AD total
plt.figure(figsize=(10, 8))
scatter = plt.scatter(vaf, qual, c=ad_total, cmap='coolwarm', alpha=0.6, s=30)
plt.colorbar(scatter, label='Support total (AD)')
plt.xlabel("VAF (Variant Allele Frequency)")
plt.ylabel("Qualité Phred")
plt.title("Qualité vs VAF (coloré par Support total)")
plt.axvline(x=0.5, color='black', linestyle='--', alpha=0.7, label='VAF = 50%')
plt.axvline(x=0.05, color='gray', linestyle='--', alpha=0.7, label='VAF = 5%')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.vaf_qual_ad.png", dpi=300, bbox_inches='tight')
plt.close()

# Graphique combiné 2x3 - Vue d'ensemble avec AD
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(20, 12))

# Support total vs Qualité coloré par VAF
scatter1 = ax1.scatter(qual, ad_total, c=vaf, cmap='viridis', alpha=0.6, s=20)
ax1.set_xlabel("Qualité Phred")
ax1.set_ylabel("Support total (AD)")
ax1.set_title("Support total vs Qualité (coloré par VAF)")
plt.colorbar(scatter1, ax=ax1, label='VAF')
ax1.grid(True, alpha=0.3)

# VAF vs Support total coloré par Qualité
scatter2 = ax2.scatter(vaf, ad_total, c=qual, cmap='plasma', alpha=0.6, s=20)
ax2.set_xlabel("VAF")
ax2.set_ylabel("Support total (AD)")
ax2.set_title("Support total vs VAF (coloré par Qualité)")
ax2.axvline(x=0.5, color='red', linestyle='--', alpha=0.7)
plt.colorbar(scatter2, ax=ax2, label='Qualité')
ax2.grid(True, alpha=0.3)

# Histogramme VAF
ax3.hist(vaf, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
ax3.set_xlabel("VAF")
ax3.set_ylabel("Nombre de variants")
ax3.set_title("Distribution de la VAF")
ax3.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='50%')
ax3.axvline(x=0.05, color='orange', linestyle='--', alpha=0.7, label='5%')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Support alternatif vs Qualité coloré par VAF
scatter4 = ax4.scatter(qual, ad_alt, c=vaf, cmap='viridis', alpha=0.6, s=20)
ax4.set_xlabel("Qualité Phred")
ax4.set_ylabel("Support alternatif (AD alt)")
ax4.set_title("Support alternatif vs Qualité (coloré par VAF)")
plt.colorbar(scatter4, ax=ax4, label='VAF')
ax4.grid(True, alpha=0.3)

# Comparaison DP vs AD total
ax5.scatter(dp_total, ad_total, c=vaf, cmap='coolwarm', alpha=0.6, s=20)
ax5.plot([0, max(dp_total)], [0, max(dp_total)], 'r--', alpha=0.5, label='AD = DP')
ax5.set_xlabel("Profondeur totale (DP)")
ax5.set_ylabel("Support total (AD)")
ax5.set_title("Comparaison DP vs AD total")
ax5.legend()
ax5.grid(True, alpha=0.3)

# Histogramme Support total
ax6.hist(ad_total, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
ax6.set_xlabel("Support total (AD)")
ax6.set_ylabel("Nombre de variants")
ax6.set_title("Distribution du Support total")
ax6.axvline(x=np.median(ad_total), color='red', linestyle='--', 
           label=f'Médiane = {np.median(ad_total):.0f}')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_directory}{output_file}.combined_ad_analysis.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 7. ANALYSE PAR CATÉGORIES DE VAF AVEC AD ***
# Définir des catégories de VAF pour l'analyse
vaf_categories = []
for v in vaf:
    if v < 0.1:
        vaf_categories.append("Très rare (< 10%)")
    elif v < 0.3:
        vaf_categories.append("Rare (10-30%)")
    elif v < 0.7:
        vaf_categories.append("Intermédiaire (30-70%)")
    else:
        vaf_categories.append("Majoritaire (> 70%)")

vaf_categories = np.array(vaf_categories)

# Box plot du support par catégorie VAF
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

unique_categories = ["Très rare (< 10%)", "Rare (10-30%)", "Intermédiaire (30-70%)", "Majoritaire (> 70%)"]

# Support total par catégorie VAF
ad_total_by_category = [ad_total[vaf_categories == cat] for cat in unique_categories]
box_plot1 = ax1.boxplot(ad_total_by_category, labels=unique_categories, patch_artist=True)
for patch, color in zip(box_plot1['boxes'], ['lightblue', 'lightgreen', 'lightyellow', 'lightcoral']):
    patch.set_facecolor(color)
ax1.set_ylabel("Support total (AD)")
ax1.set_xlabel("Catégorie de VAF")
ax1.set_title("Support total par catégorie de VAF")
ax1.tick_params(axis='x', rotation=45)
ax1.grid(True, alpha=0.3)

# Support alternatif par catégorie VAF
ad_alt_by_category = [ad_alt[vaf_categories == cat] for cat in unique_categories]
box_plot2 = ax2.boxplot(ad_alt_by_category, labels=unique_categories, patch_artist=True)
for patch, color in zip(box_plot2['boxes'], ['lightblue', 'lightgreen', 'lightyellow', 'lightcoral']):
    patch.set_facecolor(color)
ax2.set_ylabel("Support alternatif (AD alt)")
ax2.set_xlabel("Catégorie de VAF")
ax2.set_title("Support alternatif par catégorie de VAF")
ax2.tick_params(axis='x', rotation=45)
ax2.grid(True, alpha=0.3)

# Qualité par catégorie VAF
qual_by_category = [qual[vaf_categories == cat] for cat in unique_categories]
box_plot3 = ax3.boxplot(qual_by_category, labels=unique_categories, patch_artist=True)
for patch, color in zip(box_plot3['boxes'], ['lightblue', 'lightgreen', 'lightyellow', 'lightcoral']):
    patch.set_facecolor(color)
ax3.set_ylabel("Qualité Phred")
ax3.set_xlabel("Catégorie de VAF")
ax3.set_title("Qualité par catégorie de VAF")
ax3.tick_params(axis='x', rotation=45)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_directory}{output_file}.vaf_categories_ad_analysis.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 8. ANALYSE DE LA RELATION DP/AD ***
# Le ratio AD/DP a déjà été calculé plus haut

plt.figure(figsize=(12, 8))
scatter = plt.scatter(dp_total, ad_dp_ratio, c=vaf, cmap='viridis', alpha=0.6, s=30)
plt.colorbar(scatter, label='VAF')
plt.xlabel("Profondeur totale (DP)")
plt.ylabel("Ratio AD/DP")
plt.title("Efficacité du variant call (AD/DP) vs Profondeur")
plt.axhline(y=0.8, color='red', linestyle='--', alpha=0.7, label='Ratio 80% (bon)')
plt.axhline(y=0.5, color='orange', linestyle='--', alpha=0.7, label='Ratio 50% (moyen)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.dp_ad_efficiency.png", dpi=300, bbox_inches='tight')
plt.close()

# Distribution du ratio AD/DP
plt.figure(figsize=(10, 6))
plt.hist(ad_dp_ratio, bins=50, edgecolor="black", alpha=0.7, color='lightgreen')
plt.xlabel("Ratio AD/DP")
plt.ylabel("Nombre de variants")
plt.title("Distribution de l'efficacité du variant call (AD/DP)")
plt.axvline(x=np.median(ad_dp_ratio), color='red', linestyle='--', 
           label=f'Médiane = {np.median(ad_dp_ratio):.3f}')
plt.axvline(x=0.8, color='orange', linestyle='--', alpha=0.7, label='Efficacité 80%')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"{output_directory}{output_file}.ad_dp_ratio_dist.png", dpi=300, bbox_inches='tight')
plt.close()

# *** 9. ANALYSE PAR PERCENTILES POUR FILTRAGE ***
percentiles = [5, 10, 25, 50, 75, 90, 95]

print("\n" + "="*60)
print("ANALYSE STATISTIQUE POUR DÉFINIR DES FILTRES")
print("="*60)

print("\nSTATISTIQUES DESCRIPTIVES:")
print(f"Nombre total de variants: {len(pos):,}")

print(f"\nSupport alternatif (AD alt):")
for p in percentiles:
    value = np.percentile(ad_alt, p)
    count_above = np.sum(ad_alt >= value)
    percentage = count_above / len(ad_alt) * 100
    print(f"  P{p:2d}: {value:6.0f} - {count_above:6,} variants ({percentage:5.1f}%) ≥ cette valeur")

print(f"\nSupport total (AD total):")
for p in percentiles:
    value = np.percentile(ad_total, p)
    count_above = np.sum(ad_total >= value)
    percentage = count_above / len(ad_total) * 100
    print(f"  P{p:2d}: {value:6.0f} - {count_above:6,} variants ({percentage:5.1f}%) ≥ cette valeur")

print(f"\nQualité Phred:")
for p in percentiles:
    value = np.percentile(qual, p)
    count_above = np.sum(qual >= value)
    percentage = count_above / len(qual) * 100
    print(f"  P{p:2d}: {value:6.1f} - {count_above:6,} variants ({percentage:5.1f}%) ≥ cette valeur")

print(f"\nVAF:")
for p in percentiles:
    value = np.percentile(vaf, p)
    count_above = np.sum(vaf >= value)
    percentage = count_above / len(vaf) * 100
    print(f"  P{p:2d}: {value:6.3f} - {count_above:6,} variants ({percentage:5.1f}%) ≥ cette valeur")

print(f"\nRatio AD/DP (efficacité du variant call):")
for p in percentiles:
    value = np.percentile(ad_dp_ratio, p)
    count_above = np.sum(ad_dp_ratio >= value)
    percentage = count_above / len(ad_dp_ratio) * 100
    print(f"  P{p:2d}: {value:6.3f} - {count_above:6,} variants ({percentage:5.1f}%) ≥ cette valeur")

# Analyse de la qualité du séquençage
high_efficiency = np.sum(ad_dp_ratio >= 0.8)
medium_efficiency = np.sum((ad_dp_ratio >= 0.5) & (ad_dp_ratio < 0.8))
low_efficiency = np.sum(ad_dp_ratio < 0.5)

print(f"\nQUALITÉ DU SÉQUENÇAGE (basée sur ratio AD/DP):")
print(f"  Haute efficacité (≥0.8): {high_efficiency:6,} variants ({high_efficiency/len(ad_dp_ratio)*100:5.1f}%)")
print(f"  Efficacité moyenne (0.5-0.8): {medium_efficiency:6,} variants ({medium_efficiency/len(ad_dp_ratio)*100:5.1f}%)")
print(f"  Faible efficacité (<0.5): {low_efficiency:6,} variants ({low_efficiency/len(ad_dp_ratio)*100:5.1f}%)")

# *** 10. RECOMMANDATIONS DE FILTRAGE ***
print(f"\n" + "="*60)
print("RECOMMANDATIONS DE FILTRES (À analyser sur TOUS les échantillons)")
print("="*60)

print(f"\nFILTRES BASÉS SUR LE SUPPORT ALTERNATIF (AD alt):")
print(f"  - Très strict (top 5%): AD_alt ≥ {np.percentile(ad_alt, 95):.0f}")
print(f"  - Strict (top 10%): AD_alt ≥ {np.percentile(ad_alt, 90):.0f}")
print(f"  - Modéré (top 25%): AD_alt ≥ {np.percentile(ad_alt, 75):.0f}")
print(f"  - Permissif (médiane): AD_alt ≥ {np.percentile(ad_alt, 50):.0f}")

print(f"\nFILTRES BASÉS SUR LE SUPPORT TOTAL (AD total):")
print(f"  - Très strict: AD_total ≥ {np.percentile(ad_total, 95):.0f}")
print(f"  - Strict: AD_total ≥ {np.percentile(ad_total, 90):.0f}")
print(f"  - Modéré: AD_total ≥ {np.percentile(ad_total, 75):.0f}")
print(f"  - Permissif: AD_total ≥ {np.percentile(ad_total, 50):.0f}")

print(f"\nFILTRES COMBINÉS RECOMMANDÉS:")
print(f"  - Haute confiance: AD_alt ≥ {np.percentile(ad_alt, 90):.0f} ET QUAL ≥ {np.percentile(qual, 75):.0f} ET VAF ≥ 0.05")
print(f"  - Confiance modérée: AD_alt ≥ {np.percentile(ad_alt, 75):.0f} ET QUAL ≥ {np.percentile(qual, 50):.0f} ET VAF ≥ 0.02")
print(f"  - Détection sensible: AD_alt ≥ {np.percentile(ad_alt, 50):.0f} ET QUAL ≥ 20 ET VAF ≥ 0.01")

# *** 11. SAUVEGARDE DU RAPPORT ***
with open(f"{output_directory}{output_file}.support_analysis_report.txt", "w") as f:
    f.write("RAPPORT D'ANALYSE DU SUPPORT DES VARIANTS (AD)\n")
    f.write("="*50 + "\n\n")
    f.write(f"Fichier analysé: {input_file}\n")
    f.write(f"Nombre de variants: {len(pos):,}\n\n")
    
    f.write("STATISTIQUES DU SUPPORT ALTERNATIF (AD alt):\n")
    f.write("-"*40 + "\n")
    for p in percentiles:
        value = np.percentile(ad_alt, p)
        f.write(f"P{p:2d}: {value:6.0f}\n")
    
    f.write(f"\nMoyenne: {np.mean(ad_alt):.1f}\n")
    f.write(f"Médiane: {np.median(ad_alt):.0f}\n")
    f.write(f"Écart-type: {np.std(ad_alt):.1f}\n")
    
    f.write("\nRECOMMANDATIONS DE FILTRAGE:\n")
    f.write("-"*30 + "\n")
    f.write(f"Support alternatif strict: ≥ {np.percentile(ad_alt, 90):.0f}\n")
    f.write(f"Support alternatif modéré: ≥ {np.percentile(ad_alt, 75):.0f}\n")
    f.write(f"Support alternatif permissif: ≥ {np.percentile(ad_alt, 50):.0f}\n")

print(f"\n" + "="*60)
print("IMPORTANT: Analysez ces statistiques sur TOUS vos échantillons")
print("avant de fixer un filtre unique pour l'ensemble de votre étude!")
print("="*60)

print(f"\nFichiers générés dans: {output_directory}")
print("Fichiers créés:")
print("=== ANALYSES DE BASE ===")
print(f"  - {output_file}.dp_vs_ad.png (comparaison DP vs AD)")
print(f"  - {output_file}.qual_vs_ad_alt.png (qualité vs support alternatif)")
print(f"  - {output_file}.qual_vs_ad_total.png (qualité vs support total)")
print(f"  - {output_file}.vaf_vs_ad_alt.png (VAF vs support alternatif)")
print(f"  - {output_file}.distributions.png (histogrammes)")
print("=== ANALYSES SUPPLÉMENTAIRES (style précédent avec AD) ===")
print(f"  - {output_file}.vaf_vs_ad_total_qual.png (VAF vs support total coloré par qualité)")
print(f"  - {output_file}.vaf_qual_ad.png (VAF vs qualité coloré par support)")
print(f"  - {output_file}.combined_ad_analysis.png (vue d'ensemble 2x3)")
print(f"  - {output_file}.vaf_categories_ad_analysis.png (analyse par catégories VAF)")
print("=== ANALYSES DE QUALITÉ DU SÉQUENÇAGE ===")
print(f"  - {output_file}.dp_ad_efficiency.png (efficacité AD/DP vs profondeur)")
print(f"  - {output_file}.ad_dp_ratio_dist.png (distribution ratio AD/DP)")
print("=== RAPPORT ===")
print(f"  - {output_file}.support_analysis_report.txt (rapport détaillé)")

print(f"\n" + "="*80)
print("INTERPRÉTATION DU RATIO AD/DP:")
print("="*80)
print("• Ratio AD/DP élevé (>0.8) = Séquençage de bonne qualité")
print("• Ratio AD/DP moyen (0.5-0.8) = Qualité acceptable") 
print("• Ratio AD/DP faible (<0.5) = Beaucoup de lectures écartées")
print("→ Positions avec ratio faible = variants moins fiables")
print("="*80)
