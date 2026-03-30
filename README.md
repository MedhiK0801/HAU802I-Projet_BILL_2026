# CyHV-3 BILL pipeline : SNP variants analysis

## Purpose of the project

This pipeline was used to analyze **Cyprinid herpesvirus 3 (CyHV-3)** SNP variants from Nanopore sequencing data, with a focus on the comparison between viral populations before and after thermal stress.

The scripts in this repository were used to:
- merge and clean VCF-related files when needed
- enrich Medaka SNP VCF files with depth-aware metrics
- define filtering thresholds
- identify condition-specific and shared SNPs
- annotate variants with `snpEff`
- generate summary plots
- compare shared variants with volcano plots
- extract WT and mutant protein sequences for AlphaFold 3
- compare WT and mutant predicted structures in PyMOL

This README documents **only the scripts present in `scripts finaux.zip`**, using their **current filenames**, plus the PyMOL commands used for structural comparison.

---

# 1. Scripts included in this repository

The following scripts are included in the final script set:

- `files_merge.sh`
- `fix_merged_vcf.sh`
- `fix_merged_vcf.py`
- `run_add_depth_all.sh`
- `add_variant_depth.py`
- `run_plot_depth_all.sh`
- `plot_variant_support.py`
- `count_filters.sh`
- `plot_filter_impact.R`
- `run_isec_analysis.sh`
- `extract_isec_summary.sh`
- `setup_snpeff.sh`
- `analyze_snpeff_variants.py`
- `plot_isec_results.R`
- `run_volcano_snp.sh`
- `volcano_snp.R`
- `plot_volcano_aggregated_annotated.R`
- `extract_af3_sequences.py`

---

# 2. General order of the pipeline

The project follows the sequence below:

1. Merge / clean variant files if necessary
2. Add depth-related information to Medaka SNP VCFs
3. Explore support metrics and filtering thresholds
4. Identify condition-specific and shared variants
5. Annotate variants with `snpEff`
6. Summarize results and generate plots
7. Analyze shared variants with volcano plots
8. Extract WT and mutant protein sequences
9. Run AlphaFold 3 outside the repository
10. Compare WT and mutant structures in PyMOL

---

# 3. Step-by-step description of the workflow

## Step 1 - Merge and clean VCF files

### Goal
This step is used when VCF/BAM/BAI files from different sources must be merged and cleaned before downstream SNP analysis.

### Scripts used
- `files_merge.sh`
- `fix_merged_vcf.sh`
- `fix_merged_vcf.py`

---

### `files_merge.sh`

#### Role
This script compares / merges variant and alignment-related files from different sources and prepares unified files for downstream analysis.

#### Example command

```bash
bash scripts/files_merge.sh /path/to/base_directory
```

#### Optional mode

```bash
bash scripts/files_merge.sh /path/to/base_directory --merge-only
```

#### Expected result
Merged VCF-related files and summary tables are produced in the chosen base directory.

---

### `fix_merged_vcf.sh`

#### Role
This shell wrapper detects merged VCF files containing multiple sample columns and prepares them for correction.

#### Example command

```bash
bash scripts/fix_merged_vcf.sh inputs_all
```

#### Expected result
Multi-sample VCF files are converted into corrected single-sample VCFs, while the original files are backed up.

---

### `fix_merged_vcf.py`

#### Role
This Python script performs the actual VCF correction by collapsing multiple sample columns into one representative sample.

#### Example command

```bash
python3 scripts/fix_merged_vcf.py input.vcf output.fixed.vcf
```

#### Expected result
A cleaned VCF file in which duplicated sample columns are resolved.

---

## Step 2 - Add `DP`, `AD`, and `VAF` to Medaka SNP VCF files

### Goal
Medaka SNP VCF files identify variant positions, but they do not directly contain all the support metrics needed for robust downstream filtering and comparison.

This step adds:
- `DP` = total read depth at the site
- `AD` = allele depth
- `VAF` = variant allele frequency

### Scripts used
- `run_add_depth_all.sh`
- `add_variant_depth.py`

---

### `run_add_depth_all.sh`

#### Role
Batch launcher that processes all relevant sample directories and runs the depth-enrichment workflow.

#### Example command

```bash
bash scripts/run_add_depth_all.sh
```

#### Example with explicit arguments

```bash
bash scripts/run_add_depth_all.sh inputs_all References/KHV-U_trunc.fasta 4
```

#### Expected result
For each sample, a depth-enriched VCF is produced, typically named `*.snp.depth.vcf`.

---

### `add_variant_depth.py`

#### Role
This script merges:
- a Medaka SNP VCF
- a `bcftools mpileup`-derived VCF

It transfers `DP`, `AD`, and `VAF` information into the Medaka VCF.

#### Example command

```bash
python3 scripts/add_variant_depth.py \
  sample.snp.vcf \
  sample.mpileup.vcf \
  sample.snp.depth.vcf
```

#### Expected result
A new SNP VCF containing read-support information for each variant.

---

## Step 3 - Explore support metrics and filtering thresholds

### Goal
This step is used to inspect:
- variant support
- depth
- quality
- filtering impact

It helps choose which filtering level is appropriate for downstream analysis.

### Scripts used
- `run_plot_depth_all.sh`
- `plot_variant_support.py`
- `count_filters.sh`
- `plot_filter_impact.R`

---

### `run_plot_depth_all.sh`

#### Role
This wrapper applies the selected filter to depth-enriched VCF files and generates quality-control plots.

#### Example command

```bash
bash scripts/run_plot_depth_all.sh
```

#### Example with explicit directories

```bash
bash scripts/run_plot_depth_all.sh results/depth results/plots
```

#### Expected result
One set of support plots is generated per sample.

---

### `plot_variant_support.py`

#### Role
This Python script explores the relationship between:
- `DP`
- `AD`
- `VAF`
- `QUAL`

It is mainly used to visualize how well the variants are supported.

#### Example command

```bash
python3 scripts/plot_variant_support.py \
  results/depth/P25-1.trimed1000.snp.depth.vcf \
  results/plots/P25-1
```

#### Expected result
Plots illustrating the support and quality of variants in one sample.

---

### `count_filters.sh`

#### Role
This script compares different filtering levels on all samples and counts how many variants remain after each one.

#### Example command

```bash
bash scripts/count_filters.sh
```

#### Example with explicit directory

```bash
bash scripts/count_filters.sh results/depth
```

#### Expected result
A log file and summary table reporting retained variants under different filtering strategies.

---

### `plot_filter_impact.R`

#### Role
This script reads the output of `count_filters.sh` and generates summary plots showing the effect of filtering.

#### Example command

```bash
Rscript scripts/plot_filter_impact.R results/logs/count_filters_YYYYMMDD_HHMMSS.txt
```

#### Example with explicit output directory

```bash
Rscript scripts/plot_filter_impact.R \
  results/logs/count_filters_YYYYMMDD_HHMMSS.txt \
  results/filter_plots
```

#### Expected result
Figures showing:
- mean number of variants by passage and filter level
- retention percentages
- median depth summaries

---

## Step 4 - Identify condition-specific and shared variants

### Goal
This step identifies:
- variants specific to the hot condition
- variants specific to the cold condition
- variants shared between both conditions

This is the **presence/absence** branch of the SNP analysis.

### Scripts used
- `run_isec_analysis.sh`
- `extract_isec_summary.sh`

---

### `run_isec_analysis.sh`

#### Role
This script runs the complete `isec` workflow, using predefined or custom filtering modes.

#### Example command (strict)

```bash
bash scripts/run_isec_analysis.sh strict
```

#### Example command (moderate)

```bash
bash scripts/run_isec_analysis.sh moderate
```

#### Example command (custom filter)

```bash
bash scripts/run_isec_analysis.sh custom "QUAL>=15 & DP>=100" "1 5" results/isec_custom
```

#### Expected result
VCF subsets corresponding to:
- hot-specific variants
- cold-specific variants
- shared variants

---

### `extract_isec_summary.sh`

#### Role
This script converts raw `isec` outputs into summary TSV files that are easier to interpret and plot.

#### Example command (strict)

```bash
bash scripts/extract_isec_summary.sh strict
```

#### Example command (moderate)

```bash
bash scripts/extract_isec_summary.sh moderate
```

#### Example command (both)

```bash
bash scripts/extract_isec_summary.sh both
```

#### Expected result
Summary TSV files such as:
- counts by impact
- gene-level summaries
- missense detail tables
- per-sample summaries

---

## Step 5 - Annotate variants with `snpEff`

### Goal
This step assigns predicted biological consequences to variants:
- intergenic
- synonymous
- missense
- high impact
- affected ORF / gene

### Scripts used
- `setup_snpeff.sh`
- `analyze_snpeff_variants.py`

---

### `setup_snpeff.sh`

#### Role
This script installs or configures `snpEff` and builds the custom CyHV-3 / KHV-U annotation database.

#### Example command

```bash
bash scripts/setup_snpeff.sh
```

#### Example with explicit reference files

```bash
bash scripts/setup_snpeff.sh \
  References/KHV-U_trunc.fasta \
  References/DQ657948.1.gff3
```

#### Expected result
A working `snpEff` database for the CyHV-3 reference, ready for annotation.

---

### `analyze_snpeff_variants.py`

#### Role
This script parses `snpEff`-annotated VCF files and extracts the amino-acid changes of interest, especially missense variants.

#### Example command on one annotated file

```bash
python3 scripts/analyze_snpeff_variants.py \
  --file results/isec_moderate/annotated/specific_chaud.ann.vcf \
  --output results/isec_moderate/annotated/specific_chaud_variants.py
```

#### Example command on a directory

```bash
python3 scripts/analyze_snpeff_variants.py \
  --search-dir results/isec_moderate/annotated \
  --output results/isec_moderate/annotated/all_variants.py
```

#### Expected result
A cleaner mutation list extracted from the `ANN=` field.

---

## Step 6 - Generate summary figures for the isec/snpEff branch

### Goal
This step turns the summary tables into figures for the report.

### Script used
- `plot_isec_results.R`

---

### `plot_isec_results.R`

#### Role
This script generates the main summary plots from the TSV files produced by `extract_isec_summary.sh`.

#### Example command (strict)

```bash
Rscript scripts/plot_isec_results.R strict
```

#### Example command (moderate)

```bash
Rscript scripts/plot_isec_results.R moderate
```

#### Example command (both)

```bash
Rscript scripts/plot_isec_results.R both
```

#### Expected result
Figures such as:
- variant counts by category
- top affected ORFs
- distribution of missense / high-impact variants
- comparisons between filter levels

---

## Step 7 - Analyze shared variants with volcano plots

### Goal
This step focuses on **shared variants** and tests whether their allele balance changes between P25 and P27.

This is the **frequency-shift** branch of the SNP analysis.

### Scripts used
- `run_volcano_snp.sh`
- `volcano_snp.R`
- `plot_volcano_aggregated_annotated.R`

---

### `run_volcano_snp.sh`

#### Role
This shell wrapper runs the full volcano workflow:
- apply the selected filter
- generate per-sample volcano results
- generate aggregated annotated volcano plots

#### Example command (strict)

```bash
bash scripts/run_volcano_snp.sh strict
```

#### Example command (moderate)

```bash
bash scripts/run_volcano_snp.sh moderate
```

#### Example command (both)

```bash
bash scripts/run_volcano_snp.sh both
```

#### Example command (custom)

```bash
bash scripts/run_volcano_snp.sh custom "QUAL>=18 & VAF>=0.12 & AD[0:1]>=30 & DP>=100"
```

#### Expected result
A complete volcano-analysis output directory.

---

### `volcano_snp.R`

#### Role
This R script performs the actual pairwise comparison between two samples (typically P25 vs P27), runs Fisher’s exact test, computes odds ratios, and generates volcano plots.

#### Example command

```bash
Rscript scripts/volcano_snp.R \
  results/depth/P25-1.trimed1000.snp.depth.vcf.gz \
  results/depth/P27-1.trimed1000.snp.depth.vcf.gz \
  References/DQ657948.1.gff3 \
  results/volcano/P25-1_vs_P27-1 \
  1000 5 10 0.05
```

#### Expected result
Per-pair result tables and volcano plots.

---

### `plot_volcano_aggregated_annotated.R`

#### Role
This script combines all pairwise volcano results and produces an annotated aggregated volcano plot.

#### Example command

```bash
Rscript scripts/plot_volcano_aggregated_annotated.R \
  results/volcano_moderate \
  results/volcano_moderate/aggregated \
  2 \
  moderate
```

#### Expected result
Annotated aggregated volcano plots highlighting recurrent ORFs.

---

## Step 8 - Extract WT and mutant protein sequences

### Goal
After selecting candidate non-synonymous variants, this step generates:
- WT protein sequences
- mutant protein sequences
- AlphaFold 3 input files

### Script used
- `extract_af3_sequences.py`

---

### `extract_af3_sequences.py`

#### Role
This script reads:
- the reference genome FASTA
- the GFF3 annotation
- a TSV file listing candidate variants

and produces WT and mutant protein FASTA files.

#### Example candidate TSV

```tsv
Gene	Mutation	Condition
CyHV3_ORF149	p.Thr53Ala	hot
```

#### Example command

```bash
python3 scripts/extract_af3_sequences.py \
  --gff References/DQ657948.1.gff3 \
  --fasta References/KHV-U_trunc.fasta \
  --variants results/af3/orf149_candidate.tsv \
  --output results/af3/orf149
```

#### Expected result
Files such as:
- `*_WT.fasta`
- `*_MUT.fasta`
- `*_WT_config.json`
- `*_MUT_config.json`
- `mutations_summary.tsv`

---

# 4. PyMOL commands used for structural comparison

After AlphaFold 3 prediction, PyMOL was used to compare WT and mutant protein structures manually.

The goal was to determine whether a mutation was associated with a **credible local structural change**.

---

## Load WT and mutant structures

```pymol
load ORF149_WT.cif, wt
load ORF149_T53A_MUT.cif, mut
```

### Why this is used
- imports the structures
- creates short names for easier manipulation

---

## Simplify the display

```pymol
hide everything
show cartoon, wt
show cartoon, mut
color cyan, wt
color magenta, mut
```

### Why this is used
- removes visual clutter
- displays the backbone clearly
- distinguishes WT from mutant

---

## Align mutant onto WT

```pymol
align mut, wt
```

### Why this is used
- superposes both structures in the same coordinate frame
- returns an initial RMSD estimate

---

## Select the mutated residue

```pymol
select site_wt, wt and resi 53
select site_mut, mut and resi 53
show sticks, site_wt or site_mut
zoom site_wt or site_mut, 10
```

### Why this is used
- isolates the mutation site
- makes the side chains visible
- focuses the view on the biologically relevant position

---

## Define the local environment

```pymol
select local_wt, byres (wt within 8 of site_wt)
select local_mut, byres (mut within 8 of site_mut)
show sticks, local_wt or local_mut
```

### Why this is used
- defines the neighborhood around the mutation
- allows local structural comparison rather than whole-protein comparison only

---

## Measure the local RMSD

```pymol
rms_cur local_mut and name CA, local_wt and name CA
```

### Why this is used
- measures how much the local backbone differs
- focuses on alpha carbons for a robust backbone-based comparison
- was used to quantify local structural differences around the mutation

---

## Measure the global RMSD

```pymol
rms_cur mut and name CA, wt and name CA
```

### Why this is used
- measures the overall difference between WT and mutant after alignment
- must be interpreted carefully because flexible regions can inflate this value

---

## Read the local AlphaFold confidence score

```pymol
iterate wt and resi 53 and name CA, print("WT resid 53 pLDDT/B =", b)
iterate mut and resi 53 and name CA, print("MUT resid 53 pLDDT/B =", b)
```

### Why this is used
- AlphaFold confidence is often stored in the B-factor field
- this command prints the score at the mutated residue
- useful to judge whether the local comparison is reliable

---

## Compute the mean confidence in a local window

```pymol
python
vals=[]
cmd.iterate("wt and resi 40-65 and name CA", "vals.append(b)", space={"vals": vals})
print("WT mean pLDDT/B 40-65 =", sum(vals)/len(vals) if vals else "NA")
python end

python
vals=[]
cmd.iterate("mut and resi 40-65 and name CA", "vals.append(b)", space={"vals": vals})
print("MUT mean pLDDT/B 40-65 =", sum(vals)/len(vals) if vals else "NA")
python end
```

### Why this is used
- one residue alone can be misleading
- the local average gives a more stable confidence estimate for the whole region around the mutation

---

# 5. Scripts generating figures

The following scripts produce figures directly:

- `plot_variant_support.py`
- `plot_filter_impact.R`
- `plot_isec_results.R`
- `volcano_snp.R`
- `plot_volcano_aggregated_annotated.R`

These scripts cover:
- support and depth exploration
- filtering impact
- condition-specific SNP summaries
- volcano plots
- aggregated annotated volcano plots

---

# 6. Recommended practical order

A practical order of use is:

```text
1. files_merge.sh
2. fix_merged_vcf.sh
3. fix_merged_vcf.py
4. run_add_depth_all.sh
5. add_variant_depth.py
6. run_plot_depth_all.sh
7. plot_variant_support.py
8. count_filters.sh
9. plot_filter_impact.R
10. run_isec_analysis.sh
11. setup_snpeff.sh
12. extract_isec_summary.sh
13. analyze_snpeff_variants.py
14. plot_isec_results.R
15. run_volcano_snp.sh
16. volcano_snp.R
17. plot_volcano_aggregated_annotated.R
18. extract_af3_sequences.py
19. AlphaFold 3
20. PyMOL
```

---

# 7. Interpretation logic of the project

The analysis is based on two complementary SNP branches:

## Condition-specific / shared variant branch
Driven mainly by:
- `run_isec_analysis.sh`
- `extract_isec_summary.sh`
- `setup_snpeff.sh`
- `plot_isec_results.R`

This branch answers:
- which SNPs are hot-specific?
- which SNPs are cold-specific?
- which ORFs are recurrently affected?

## Shared-variant frequency branch
Driven mainly by:
- `run_volcano_snp.sh`
- `volcano_snp.R`
- `plot_volcano_aggregated_annotated.R`

This branch answers:
- among shared SNPs, which ones change in frequency between P25 and P27?
- which variants are enriched in P27 or lost after thermal stress?

## Structural branch
Driven mainly by:
- `extract_af3_sequences.py`
- AlphaFold 3
- PyMOL

This branch answers:
- can selected non-synonymous mutations be associated with plausible local structural effects?

---

# 8. Final note

This repository documents the final set of scripts used to move from:
- raw or merged VCF files
to:
- depth-aware SNP calls
- filtered condition-specific/shared variants
- annotated candidate mutations
- exploratory structural interpretation of selected protein variants

