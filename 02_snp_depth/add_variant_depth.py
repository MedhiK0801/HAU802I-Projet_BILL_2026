#!/usr/bin/env python3
"""
add_variant_depth.py

Merge a Medaka SNP VCF with a matching mpileup-derived VCF in order to copy the
read-support metrics needed for downstream filtering and plotting.

Added fields:
- DP  : total depth at the site
- AD  : REF and ALT allelic depth
- VAF : ALT allele frequency

Usage:
    python3 add_variant_depth.py <medaka.snp.vcf> <sample.mpileup.vcf> [output.depth.vcf]
"""

import sys
from typing import Sized

import pysam as ps


# A function to get the length of an object that is either a Sized object or None
def safe_len(obj: Sized | None) -> int:
    return len(obj) if obj is not None else 0


# Check if the inputs files are provided and has a .vcf extension
if (
    len(sys.argv) < 3
    or not sys.argv[1].lower().endswith(".vcf")
    or not sys.argv[2].lower().endswith(".vcf")
):
    print(
        "\nUsage: python depth_variant_adder.py <snp_input_file.vcf> <bcf_raw_input_file.vcf> <output_file>"
    )
    sys.exit(1)

# We add a depth in the name of the output file to differenciate it
snp_input_file = sys.argv[1]
bcf_input_file = sys.argv[2]
if len(sys.argv) >= 4:
    output_file = sys.argv[3]
else:
    output_file = snp_input_file.replace(".vcf", ".depth.vcf")

# We use the pysam.VariantFile to manipulate vcf input files
snp_vcf_in = ps.VariantFile(snp_input_file)
bcf_vcf_in = ps.VariantFile(bcf_input_file)

# We add the description for depth and allelic depth
snp_vcf_in.header.formats.add(
    "DP", number=1, type="Integer", description="Number of high-quality bases"
)

snp_vcf_in.header.formats.add(
    "AD", number="R", type="Integer", description="Allelic depths (high-quality bases)"
)

snp_vcf_in.header.formats.add(
    "VAF", number="A", type="Float", description="Variant allele frequency for each ALT allele"
)

# We open the output file in writing mode with the header being header of the first input file with the decription of DP and AD added
vcf_out = ps.VariantFile(output_file, "w", header=snp_vcf_in.header)

# We iterate through the each line (record) of both vcf (assumed to be both with the same position beacause the mpileup used the first input vcf as refernce but in case it's not the case, it throws a runtimeerror)
for rec_snp, rec_bcf in zip(snp_vcf_in.fetch(), bcf_vcf_in.fetch()):
    key_snp = (rec_snp.chrom, rec_snp.pos)
    key_bcf = (rec_bcf.chrom, rec_bcf.pos)
    if key_snp != key_bcf:
        raise RuntimeError(
            f"ERROR :Décalage entre deux lignes : snp file (",
            snp_input_file,
            ") at pos :",
            key_snp,
            " bcf file (",
            bcf_input_file,
            ") at pos :",
            key_bcf,
        )
        # Here, the file is supposed to have only one sample (not merged)
    if len(rec_snp.samples.keys()) > 1:
        raise RuntimeError("Too many samples in first input file(", snp_input_file, ")")
    if len(rec_bcf.samples.keys()) > 1:
        raise RuntimeError(
            "Too many samples in second input file(", bcf_input_file, ")"
        )
    if len(rec_snp.samples.keys()) == 0:
        raise RuntimeError("No samples in first input file(", snp_input_file, ")")
    if len(rec_bcf.samples.keys()) == 0:
        raise RuntimeError("No samples in second input file(", bcf_input_file, ")")

    # Get the depth (DP) and allelic depth (AD) from second file input (raw vcf from bcftool mpileup)
    bcf_sample_name = rec_bcf.samples.keys()[0]
    dp = int(rec_bcf.samples[bcf_sample_name]["DP"])
    ad = tuple(rec_bcf.samples[bcf_sample_name]["AD"])
    snp_nb_var = safe_len(rec_snp.alleles)
    bcf_nb_var = safe_len(rec_bcf.alleles)
    if snp_nb_var > bcf_nb_var:
        ad = tuple(list(ad) + [0] * (snp_nb_var - bcf_nb_var))
        print(
            "WARNING : alts in first input file (",
            snp_input_file,
            ") are greater than alts in second input file (",
            bcf_input_file,
            ") at position :",
            key_snp,
        )
    # Keep only the longer alleles coresponding to the one present in the dirst input file (snp.vcf from medaka)
    else:
        ad = tuple(ad[:snp_nb_var])
    
    # Calcul du VAF pour chaque ALT
    # ad[0] = REF, ad[1:] = ALT(s)
    if dp > 0:
        vaf = tuple(float(alt_depth) / dp for alt_depth in ad[1:])
    else:
        vaf = tuple(0.0 for _ in ad[1:])

    # Add Depth and allelice depth to the record of the .snp.vcf file to write it in the outputfile
    snp_sample_name = rec_snp.samples.keys()[0]
    rec_snp.samples[snp_sample_name]["DP"] = dp
    rec_snp.samples[snp_sample_name]["AD"] = ad
    rec_snp.samples[snp_sample_name]["VAF"] = vaf

    # Write the record(line) of the snp.vcf with the DP and AD added in the outputfile
    _ = vcf_out.write(rec_snp)

vcf_out.close()