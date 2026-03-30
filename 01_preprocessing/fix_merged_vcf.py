#!/usr/bin/env python3
"""
fix_merged_vcf.py

Collapse a VCF file that contains multiple sample columns into a single-sample
VCF. This is useful when an upstream merge created one file containing several
sample-specific FORMAT columns that should be treated as one combined signal.

Usage:
    python3 fix_merged_vcf.py <input.vcf> [output.vcf]
"""

import sys
import pysam as ps


def safe_get(sample, field, default=None):
    """Récupère un champ FORMAT de manière sûre."""
    try:
        val = sample[field]
        if val is None:
            return default
        return val
    except KeyError:
        return default


def main():
    if len(sys.argv) < 2 or not sys.argv[1].lower().endswith(".vcf"):
        print("Usage: python fix_merged_vcf.py <input.vcf> [output.vcf]")
        sys.exit(1)

    input_file = sys.argv[1]
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        output_file = input_file.replace(".vcf", ".fixed.vcf")

    vcf_in = ps.VariantFile(input_file)
    sample_names = list(vcf_in.header.samples)
    n_samples = len(sample_names)

    if n_samples <= 1:
        print(f"[INFO] {input_file} : déjà single-sample ({n_samples} sample), rien à faire.")
        vcf_in.close()
        # Copier tel quel
        import shutil
        shutil.copy2(input_file, output_file)
        return

    print(f"[FIX] {input_file} : {n_samples} samples détectés ({', '.join(sample_names)})")
    print(f"      -> Fusion en un seul super-échantillon 'SAMPLE'")

    # Créer un nouveau header avec un seul sample
    new_header = ps.VariantHeader()

    # Copier les lignes de métadonnées (contigs, INFO, FORMAT, FILTER, etc.)
    for rec in vcf_in.header.records:
        new_header.add_record(rec)

    # Ajouter un seul sample
    new_header.add_sample("SAMPLE")

    vcf_out = ps.VariantFile(output_file, "w", header=new_header)

    n_variants = 0
    n_merged = 0

    for rec in vcf_in.fetch():
        n_variants += 1

        # Trouver le "meilleur" sample : celui avec le GQ le plus élevé
        best_idx = 0
        best_gq = -1

        for i, sname in enumerate(sample_names):
            gq = safe_get(rec.samples[sname], "GQ", 0)
            if isinstance(gq, tuple):
                gq = gq[0] if len(gq) > 0 else 0
            if gq is None:
                gq = 0
            if gq > best_gq:
                best_gq = gq
                best_idx = i

        best_sample = rec.samples[sample_names[best_idx]]

        # Vérifier si les samples divergent (au moins un a un GT différent)
        gts = set()
        for sname in sample_names:
            gt = safe_get(rec.samples[sname], "GT", None)
            if gt is not None:
                gts.add(gt)
        if len(gts) > 1:
            n_merged += 1

        # Créer le nouveau record
        new_rec = vcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual,
            filter=rec.filter,
            info=rec.info,
        )

        # Copier les champs FORMAT depuis le meilleur sample
        out_sample = new_rec.samples["SAMPLE"]
        for key in rec.format.keys():
            try:
                val = best_sample[key]
                if val is not None:
                    out_sample[key] = val
            except (KeyError, TypeError):
                pass

        vcf_out.write(new_rec)

    vcf_out.close()
    vcf_in.close()

    print(f"      {n_variants} variants traités, {n_merged} avaient des GT divergents entre samples")
    print(f"      -> {output_file}")


if __name__ == "__main__":
    main()
