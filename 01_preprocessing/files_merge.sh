#!/bin/bash
# -----------------------------------------------------------------------------
# files_merge.sh
#
# Compare Floris and Veli VCF/BAM/BAI files, then identify and merge
# complementary variant information when needed.
#
# Usage:
#   bash scripts/files_merge.sh <base_dir> [--merge-only]
# -----------------------------------------------------------------------------

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
PASSAGES=(15 30 50)
SAMPLES=($(seq 1 10))
BASE_DIR="${1:-.}"
MERGE_ONLY=false
[[ "${2:-}" == "--merge-only" ]] && MERGE_ONLY=true

LOG_FILE="${BASE_DIR}/vcf_compare_$(date +%Y%m%d_%H%M%S).log"
SUMMARY_FILE="${BASE_DIR}/veli_specific_counts.tsv"

# ── Logging ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG_FILE"; }
warn() { echo "[$(date '+%H:%M:%S')] WARNING: $*" | tee -a "$LOG_FILE" >&2; }

# ── Functions ─────────────────────────────────────────────────────────────────

# Find a VCF matching the given pattern (Floris or Veli) for a passage/sample
find_vcf() {
    local dir="$1" passage="$2" sample="$3" type="$4"
    local pattern

    if [[ "$type" == "Veli" ]]; then
        # Veli pattern: P{XX}-{S}_Veli*.snp.vcf (but NOT .vcf.gz)
        pattern="P${passage}-${sample}_Veli"
    else
        # Floris pattern: various forms like FlorisF1, f2FLoris, f2Floris03052022, etc.
        # Case-insensitive match on "floris"
        pattern="P${passage}-${sample}_"
    fi

    local found=""
    # Search both .snp.vcf and .snp.vcf.gz in one pass
    for f in "${dir}"/P${passage}-${sample}_*.snp.vcf "${dir}"/P${passage}-${sample}_*.snp.vcf.gz; do
        [[ -e "$f" ]] || continue
        # Skip .vcf.gz if the base .vcf was already matched
        # (avoid matching both file.vcf and file.vcf.gz)
        if [[ "$f" == *.vcf.gz && -n "$found" ]]; then
            continue
        fi
        # Skip files inside comparison_* output directories
        [[ "$f" == *comparison_* ]] && continue

        local basename_f
        basename_f=$(basename "$f")
        local lower_f
        lower_f=$(echo "$basename_f" | tr '[:upper:]' '[:lower:]')

        if [[ "$type" == "Veli" ]]; then
            # Must contain "veli" AND must NOT contain "floris"
            if [[ "$lower_f" == *veli* && "$lower_f" != *floris* ]]; then
                found="$f"
                break
            fi
        else
            # Must contain "floris" AND must NOT contain "veli"
            if [[ "$lower_f" == *floris* && "$lower_f" != *veli* ]]; then
                found="$f"
                break
            fi
        fi
    done

    echo "$found"
}

# Find a BAM matching the given pattern (Floris or Veli) for a passage/sample
find_bam() {
    local dir="$1" passage="$2" sample="$3" type="$4"

    local found=""
    for f in "${dir}"/P${passage}-${sample}_*.bam; do
        [[ -e "$f" ]] || continue
        # Skip .bam.bai index files
        [[ "$f" == *.bai ]] && continue
        # Skip files inside comparison_* output directories
        [[ "$f" == *comparison_* ]] && continue

        local basename_f
        basename_f=$(basename "$f")
        local lower_f
        lower_f=$(echo "$basename_f" | tr '[:upper:]' '[:lower:]')

        if [[ "$type" == "Veli" ]]; then
            if [[ "$lower_f" == *veli* && "$lower_f" != *floris* ]]; then
                found="$f"
                break
            fi
        else
            if [[ "$lower_f" == *floris* && "$lower_f" != *veli* ]]; then
                found="$f"
                break
            fi
        fi
    done

    echo "$found"
}

# Compress and index a VCF if needed. Returns path to .vcf.gz
compress_and_index() {
    local vcf="$1"

    if [[ "$vcf" == *.vcf.gz ]]; then
        # Already compressed — just make sure it's indexed
        if [[ ! -f "${vcf}.tbi" ]]; then
            log "  Indexing (tabix): $vcf" >&2
            bcftools index -t "$vcf"
        fi
        echo "$vcf"
        return
    fi

    local gz="${vcf}.gz"
    if [[ ! -f "$gz" ]]; then
        log "  Compressing (bgzip): $vcf" >&2
        bgzip -c "$vcf" > "$gz"
    fi
    if [[ ! -f "${gz}.tbi" ]]; then
        log "  Indexing (tabix): $gz" >&2
        bcftools index -t "$gz"
    fi
    echo "$gz"
}

# ── Main ──────────────────────────────────────────────────────────────────────

log "============================================================"
log "VCF Compare — Floris vs Veli"
log "Base directory : $BASE_DIR"
log "Passages       : ${PASSAGES[*]}"
log "Samples        : 1-10"
log "Mode           : $( $MERGE_ONLY && echo 'merge-only' || echo 'isec + merge' )"
log "============================================================"

# Header for summary TSV
echo -e "Passage\tSample\tFloris_VCF\tVeli_VCF\tVeli_specific_variants\tOutput_merged_VCF\tOutput_merged_BAM" > "$SUMMARY_FILE"

for passage in "${PASSAGES[@]}"; do
    # The VCFs may be directly in base_dir or in a P{XX} subdirectory — try both
    search_dir="$BASE_DIR"
    [[ -d "${BASE_DIR}/P${passage}" ]] && search_dir="${BASE_DIR}/P${passage}"

    for sample in "${SAMPLES[@]}"; do
        log "────────────────────────────────────────"
        log "Processing P${passage}-${sample}"

        floris_vcf=$(find_vcf "$search_dir" "$passage" "$sample" "Floris")
        veli_vcf=$(find_vcf "$search_dir" "$passage" "$sample" "Veli")

        # ── Check for pair existence ──────────────────────────────────────
        if [[ -z "$floris_vcf" && -z "$veli_vcf" ]]; then
            warn "P${passage}-${sample}: No Floris or Veli VCF found — skipping."
            continue
        fi
        if [[ -z "$floris_vcf" ]]; then
            warn "P${passage}-${sample}: No Floris VCF found (Veli exists: $(basename "$veli_vcf")) — skipping."
            continue
        fi
        if [[ -z "$veli_vcf" ]]; then
            warn "P${passage}-${sample}: No Veli VCF found (Floris exists: $(basename "$floris_vcf")) — skipping."
            continue
        fi

        log "  Floris: $(basename "$floris_vcf")"
        log "  Veli  : $(basename "$veli_vcf")"

        # ── Compress & index ──────────────────────────────────────────────
        floris_gz=$(compress_and_index "$floris_vcf")
        veli_gz=$(compress_and_index "$veli_vcf")

        # ── Output directory ──────────────────────────────────────────────
        out_dir="${search_dir}/comparison_P${passage}-${sample}"
        mkdir -p "$out_dir"

        veli_specific_count="N/A"
        merged_output="${out_dir}/P${passage}-${sample}_merged.vcf.gz"

        if $MERGE_ONLY; then
            # ── MERGE ONLY approach ───────────────────────────────────────
            log "  Running bcftools merge..."
            bcftools merge \
                --force-samples \
                -O z -o "$merged_output" \
                "$floris_gz" "$veli_gz"
            bcftools index -t "$merged_output"

            # Count variants present only in Veli (missing in Floris → "./." in Floris column)
            veli_specific_count=$(bcftools view "$merged_output" | \
                awk -F'\t' 'BEGIN{c=0} !/^#/{if($10 ~ /^\./) c++} END{print c}')
            log "  Merged output : $(basename "$merged_output")"
            log "  Veli-specific variants (approx.): $veli_specific_count"

        else
            # ── ISEC approach ─────────────────────────────────────────────
            isec_dir="${out_dir}/isec"
            rm -rf "$isec_dir"
            log "  Running bcftools isec..."
            bcftools isec \
                -p "$isec_dir" \
                -O z \
                "$floris_gz" "$veli_gz"
            # isec outputs:
            #   0000.vcf.gz = variants unique to Floris
            #   0001.vcf.gz = variants unique to Veli   ← this is what we want
            #   0002.vcf.gz = variants shared (Floris version)
            #   0003.vcf.gz = variants shared (Veli version)

            veli_only="${isec_dir}/0001.vcf.gz"
            if [[ -f "$veli_only" ]]; then
                veli_specific_count=$(bcftools view -H "$veli_only" | wc -l)
                log "  Veli-specific variants: $veli_specific_count"

                # ── Add Veli-specific variants to Floris ──────────────────
                if [[ "$veli_specific_count" -gt 0 ]]; then
                    # Index the Veli-only file
                    bcftools index -t "$veli_only"

                    log "  Merging Veli-specific variants into Floris..."
                    # Concat Floris + Veli-specific, then sort
                    bcftools concat \
                        -a -O z -o "${out_dir}/P${passage}-${sample}_floris_plus_veli_specific.vcf.gz" \
                        "$floris_gz" "$veli_only"
                    bcftools sort \
                        -O z -o "$merged_output" \
                        "${out_dir}/P${passage}-${sample}_floris_plus_veli_specific.vcf.gz"
                    bcftools index -t "$merged_output"
                    rm -f "${out_dir}/P${passage}-${sample}_floris_plus_veli_specific.vcf.gz"
                    log "  Output: $(basename "$merged_output")"
                else
                    log "  No Veli-specific variants to add."
                    cp "$floris_gz" "$merged_output"
                    bcftools index -t "$merged_output"
                fi
            else
                warn "  isec output 0001.vcf.gz not found."
                veli_specific_count=0
            fi
        fi

        # ── BAM merge ─────────────────────────────────────────────────────
        floris_bam=$(find_bam "$search_dir" "$passage" "$sample" "Floris")
        veli_bam=$(find_bam "$search_dir" "$passage" "$sample" "Veli")
        merged_bam=""

        if [[ -n "$floris_bam" && -n "$veli_bam" ]]; then
            merged_bam="${out_dir}/P${passage}-${sample}_merged.bam"
            log "  Merging BAMs..."
            log "    Floris BAM: $(basename "$floris_bam")"
            log "    Veli BAM  : $(basename "$veli_bam")"
            samtools merge -f "$merged_bam" "$floris_bam" "$veli_bam"
            samtools index "$merged_bam"
            log "  Merged BAM: $(basename "$merged_bam")"
        else
            [[ -z "$floris_bam" ]] && warn "P${passage}-${sample}: No Floris BAM found."
            [[ -z "$veli_bam" ]]   && warn "P${passage}-${sample}: No Veli BAM found."
            warn "P${passage}-${sample}: BAM merge skipped (missing file)."
        fi

        # ── Write summary line ────────────────────────────────────────────
        echo -e "P${passage}\t${sample}\t$(basename "$floris_vcf")\t$(basename "$veli_vcf")\t${veli_specific_count}\t$(basename "$merged_output")\t$(basename "${merged_bam:-N/A}")" \
            >> "$SUMMARY_FILE"
    done
done

log ""
log "============================================================"
log "Done. Summary written to: $SUMMARY_FILE"
log "Log written to: $LOG_FILE"
log "============================================================"

cat "$SUMMARY_FILE"
