#!/usr/bin/env python3
"""
analyze_snpeff_variants.py

Parse snpEff-annotated VCF files and extract the amino-acid changes that are
relevant for the candidate-variant and structural-analysis parts of the project.

Usage examples:
    python3 analyze_snpeff_variants.py --search-dir results/
    python3 analyze_snpeff_variants.py --file results/isec_moderate/annotated/specific_chaud.ann.vcf
"""

import argparse
import os
import re
import glob
from collections import defaultdict

def find_snpeff_files(search_dir):
    """Chercher récursivement les fichiers snpEff"""
    patterns = [
        '**/*snpeff*.vcf',
        '**/*annotated*.vcf',
        '**/isec*/**/snpeff*.vcf',
        '**/isec*/**/*annotated*.vcf'
    ]
    
    files = []
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(search_dir, pattern), recursive=True))
    
    # Éliminer les doublons et trier
    return sorted(set(files))

def parse_snpeff_vcf(vcf_file):
    """Parser un fichier VCF avec annotations snpEff"""
    variants = []
    
    print(f"\nAnalyse de : {vcf_file}")
    
    try:
        with open(vcf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                info = parts[7]
                
                # Chercher le champ ANN (annotations snpEff)
                ann_match = re.search(r'ANN=([^;]+)', info)
                if not ann_match:
                    continue
                
                # Parser les annotations
                for ann_text in ann_match.group(1).split(','):
                    ann_fields = ann_text.split('|')
                    if len(ann_fields) < 15:
                        continue
                    
                    annotation = ann_fields[1]
                    impact = ann_fields[2]
                    gene_name = ann_fields[3]
                    hgvs_p = ann_fields[10]
                    
                    # Nous intéresse : les mutations missense
                    if 'missense' in annotation.lower():
                        variant = {
                            'file': os.path.basename(vcf_file),
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'gene': gene_name,
                            'annotation': annotation,
                            'impact': impact,
                            'hgvs_p': hgvs_p,
                            'line_num': line_num
                        }
                        variants.append(variant)
                        
    except Exception as e:
        print(f"  Erreur lecture : {e}")
        return []
    
    print(f"  Trouvé : {len(variants)} mutations missense")
    return variants

def extract_mutation_info(hgvs_p):
    """Extraire les infos de mutation depuis HGVS.p"""
    if not hgvs_p or hgvs_p == '.':
        return None
    
    # Format attendu : p.Ala123Val ou p.A123V
    match = re.match(r'p\.([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})', hgvs_p)
    if not match:
        return None
    
    ref_aa = match.group(1)
    position = int(match.group(2))
    alt_aa = match.group(3)
    
    # Convertir en format 3-lettres si nécessaire
    aa_1to3 = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    
    if len(ref_aa) == 1:
        ref_aa = aa_1to3.get(ref_aa, ref_aa)
    if len(alt_aa) == 1:
        alt_aa = aa_1to3.get(alt_aa, alt_aa)
    
    return {
        'position': position,
        'ref_aa': ref_aa,
        'alt_aa': alt_aa,
        'mutation_string': f"{ref_aa}{position}{alt_aa}"
    }

def categorize_variants_by_condition(variants):
    """Catégoriser les variants par condition (froid/chaud/partagé)"""
    
    # Essayer d'inférer la condition depuis le nom de fichier
    by_condition = defaultdict(list)
    
    for variant in variants:
        filename = variant['file'].lower()
        
        # Patterns pour identifier la condition
        if 'froid' in filename or 'cold' in filename or 'specific_0' in filename:
            condition = 'froid'
        elif 'chaud' in filename or 'hot' in filename or 'specific_1' in filename:
            condition = 'chaud'
        elif 'shared' in filename or 'common' in filename or 'partage' in filename:
            condition = 'partage'
        else:
            # Essayer de deviner depuis le chemin
            if 'froid' in os.path.dirname(variant['file']):
                condition = 'froid'
            elif 'chaud' in os.path.dirname(variant['file']):
                condition = 'chaud'
            else:
                condition = 'inconnu'
        
        by_condition[condition].append(variant)
    
    return by_condition

def analyze_target_orfs(variants):
    """Analyser les ORFs qui nous intéressent spécifiquement"""
    
    target_orfs = ['CyHV3_ORF128', 'CyHV3_ORF154', 'CyHV3_ORF145', 'CyHV3_ORF89', 'CyHV3_ORF25', 'CyHV3_ORF52']
    
    results = {}
    position_697_candidates = []
    
    for variant in variants:
        gene = variant['gene']
        
        # Analyser nos ORFs cibles
        for target_orf in target_orfs:
            if target_orf in gene or gene in target_orf:
                if target_orf not in results:
                    results[target_orf] = []
                
                mut_info = extract_mutation_info(variant['hgvs_p'])
                if mut_info:
                    variant_with_mut = variant.copy()
                    variant_with_mut.update(mut_info)
                    results[target_orf].append(variant_with_mut)
        
        # Chercher spécifiquement position 697
        mut_info = extract_mutation_info(variant['hgvs_p'])
        if mut_info and mut_info['position'] == 697:
            variant_with_mut = variant.copy()
            variant_with_mut.update(mut_info)
            position_697_candidates.append(variant_with_mut)
    
    return results, position_697_candidates

def print_analysis_results(orf_results, position_697_candidates, by_condition):
    """Afficher les résultats d'analyse"""
    
    print("\n" + "="*80)
    print("ANALYSE DES VARIANTS SNPEFF")
    print("="*80)
    
    # Résumé par condition
    print(f"\n📊 RÉSUMÉ PAR CONDITION:")
    for condition, variants in by_condition.items():
        print(f"  {condition.upper()}: {len(variants)} variants")
    
    # Analyse des ORFs cibles
    print(f"\n🎯 ANALYSE DES ORFS CIBLES:")
    for orf, variants in orf_results.items():
        print(f"\n{orf}:")
        if variants:
            for var in variants:
                condition = "inconnu"
                filename = var['file'].lower()
                if 'froid' in filename or 'cold' in filename:
                    condition = "froid"
                elif 'chaud' in filename or 'hot' in filename:
                    condition = "chaud"
                
                print(f"  ✓ {var['mutation_string']} ({condition}) - {var['hgvs_p']}")
                print(f"    Fichier: {var['file']}")
        else:
            print(f"  ❌ Aucune mutation trouvée")
    
    # Position 697
    print(f"\n🔍 MUTATIONS POSITION 697:")
    if position_697_candidates:
        for candidate in position_697_candidates:
            print(f"  ⭐ {candidate['gene']}: {candidate['mutation_string']}")
            print(f"     {candidate['hgvs_p']} (fichier: {candidate['file']})")
    else:
        print(f"  ❌ Aucune mutation position 697 trouvée")
    
    # Recherche spécifique Trp->Ser
    print(f"\n🔍 MUTATIONS TRP -> SER:")
    trp_ser_mutations = []
    for variants in orf_results.values():
        for var in variants:
            if 'ref_aa' in var and 'alt_aa' in var:
                if var['ref_aa'] == 'Trp' and var['alt_aa'] == 'Ser':
                    trp_ser_mutations.append(var)
    
    if trp_ser_mutations:
        for mut in trp_ser_mutations:
            print(f"  ⭐ {mut['gene']}: {mut['mutation_string']}")
            print(f"     Position: {mut['position']} (fichier: {mut['file']})")
    else:
        print(f"  ❌ Aucune mutation Trp->Ser trouvée")

def generate_corrected_variants(orf_results, position_697_candidates):
    """Générer la liste de variants corrigée"""
    
    print(f"\n" + "="*80)
    print("VARIANTS CORRIGÉS POUR EXTRACTION")
    print("="*80)
    
    corrected_variants = []
    
    # ORFs à traiter
    target_orfs = {
        'CyHV3_ORF89': ['chaud', 'froid'],
        'CyHV3_ORF128': ['chaud'],
        'CyHV3_ORF154': ['chaud'],
        'CyHV3_ORF25': ['chaud'],
        'CyHV3_ORF52': ['chaud', 'froid']
    }
    
    for orf, expected_conditions in target_orfs.items():
        if orf in orf_results and orf_results[orf]:
            print(f"\n{orf}:")
            for var in orf_results[orf]:
                # Inférer condition
                condition = "inconnu"
                filename = var['file'].lower()
                if 'froid' in filename or 'cold' in filename or 'specific_0' in filename:
                    condition = "froid"
                elif 'chaud' in filename or 'hot' in filename or 'specific_1' in filename:
                    condition = "chaud"
                
                if condition in expected_conditions:
                    variant_tuple = (orf, var['position'], var['ref_aa'], var['alt_aa'], condition)
                    corrected_variants.append(variant_tuple)
                    print(f"  ✓ {var['mutation_string']} ({condition})")
                else:
                    print(f"  ⚠️ {var['mutation_string']} (condition {condition} inattendue)")
        else:
            print(f"\n{orf}: ❌ Pas trouvé dans snpEff")
    
    # Position 697 - traitement spécial
    if position_697_candidates:
        print(f"\n🔍 MUTATION POSITION 697:")
        best_candidate = position_697_candidates[0]  # Prendre le premier
        condition = "chaud"  # Par défaut, à ajuster selon le fichier
        
        filename = best_candidate['file'].lower()
        if 'froid' in filename:
            condition = "froid"
        
        variant_tuple = (best_candidate['gene'], 697, best_candidate['ref_aa'], best_candidate['alt_aa'], condition)
        corrected_variants.append(variant_tuple)
        print(f"  ✓ {best_candidate['gene']}: {best_candidate['mutation_string']} ({condition})")
    
    # Générer le code Python
    print(f"\n" + "="*80)
    print("CODE PYTHON CORRIGÉ:")
    print("="*80)
    
    print("VARIANTS = [")
    for variant in corrected_variants:
        orf, pos, ref_aa, alt_aa, condition = variant
        print(f'    ("{orf}", {pos}, "{ref_aa}", "{alt_aa}", "{condition}"),')
    print("]")
    
    return corrected_variants

def main():
    parser = argparse.ArgumentParser(description='Analyser les fichiers snpEff pour extraire les vraies mutations')
    parser.add_argument('--search-dir', default='.', help='Répertoire de recherche (défaut: répertoire courant)')
    parser.add_argument('--file', help='Fichier snpEff spécifique à analyser')
    parser.add_argument('--output', help='Fichier de sortie pour les variants corrigés')
    
    args = parser.parse_args()
    
    # Trouver les fichiers à analyser
    if args.file:
        if os.path.exists(args.file):
            snpeff_files = [args.file]
        else:
            print(f"Erreur: Fichier {args.file} introuvable")
            return
    else:
        print(f"Recherche des fichiers snpEff dans : {args.search_dir}")
        snpeff_files = find_snpeff_files(args.search_dir)
        
        if not snpeff_files:
            print("Aucun fichier snpEff trouvé!")
            print("Patterns recherchés:")
            print("  - **/*snpeff*.vcf")
            print("  - **/*annotated*.vcf") 
            print("  - **/isec*/**/snpeff*.vcf")
            print("\nVérifiez que vos fichiers snpEff sont bien présents.")
            return
    
    print(f"Fichiers trouvés: {len(snpeff_files)}")
    for f in snpeff_files:
        print(f"  - {f}")
    
    # Parser tous les fichiers
    all_variants = []
    for vcf_file in snpeff_files:
        variants = parse_snpeff_vcf(vcf_file)
        all_variants.extend(variants)
    
    if not all_variants:
        print("\nAucune mutation missense trouvée dans les fichiers snpEff!")
        return
    
    print(f"\nTotal mutations missense: {len(all_variants)}")
    
    # Analyser
    by_condition = categorize_variants_by_condition(all_variants)
    orf_results, position_697_candidates = analyze_target_orfs(all_variants)
    
    # Afficher les résultats
    print_analysis_results(orf_results, position_697_candidates, by_condition)
    
    # Générer les variants corrigés
    corrected_variants = generate_corrected_variants(orf_results, position_697_candidates)
    
    # Sauvegarder si demandé
    if args.output:
        with open(args.output, 'w') as f:
            f.write("# Variants corrigés d'après analyse snpEff\n")
            f.write("VARIANTS = [\n")
            for variant in corrected_variants:
                orf, pos, ref_aa, alt_aa, condition = variant
                f.write(f'    ("{orf}", {pos}, "{ref_aa}", "{alt_aa}", "{condition}"),\n')
            f.write("]\n")
        print(f"\nVariants sauvés dans: {args.output}")

if __name__ == "__main__":
    main()
