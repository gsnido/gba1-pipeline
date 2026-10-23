#!/usr/bin/env bash
#
# LONG-NEXT.sh
#
# Implements the "tailored bioinformatic analysis" step of the LONG-NEXT
# method (Cuconato et al., Parkinsonism Relat Disord 2025) for GBA1
# sequencing analysis:
#
#   1) Sanitizes sample IDs (dashes/whitespace -> underscores).
#   2) Locates the paired FASTQ (R1/R2) files for each sample.
#   4) Aligns reads with bwa mem against the masked reference and produces
#      a sorted, indexed BAM (GATK4-compatible best-practices alignment).
#   5) Calls variants with GATK4 HaplotypeCaller.
#   6) Checks mean read depth over the GBA1 region against the paper's
#      minimum threshold of 500x, flagging samples that fall short as
#      possible allele-dropout (ADO) events requiring orthogonal
#      confirmation (e.g. MLPA).
#
# Usage:
#   ./LONG-NEXT.sh <input_id_list.txt> <fastq_dir> <reference.fa> <gba1_region_bed> <output_dir> [threads]
#
# Arguments:
#   input_id_list.txt   : one sample ID per line
#   fastq_dir           : directory containing *_R1_*/ *_R2_* FASTQ files
#   reference.fa        : unmasked hg38 reference FASTA (must be bwa- and
#                          samtools-indexed)
#   gba1_region_bed     : BED file with GBA1 gene coordinates (hg38), used
#                          for the 500x read-depth QC check
#   output_dir          : directory for masked reference, BAMs, VCFs, logs
#   threads (optional)  : number of threads for bwa/samtools/gatk (default 4)
#
# Requirements: bwa, samtools, bedtools, gatk (GATK4), bcftools must be in PATH.

set -euo pipefail

MIN_DEPTH=500      # LONG-NEXT minimum region read depth threshold (paper 2.3.3)
MIN_VARIANT_DP=30  # minimum per-variant depth for confident genotyping
MIN_VAF=0.15       # minimum variant allele frequency to trust a call
HET_VAF_LOW=0.30   # heterozygous VAF sanity window lower bound
HET_VAF_HIGH=0.70  # heterozygous VAF sanity window upper bound

usage() {
    echo "Usage: $0 <input_id_list.txt> <fastq_dir> <masked_reference.fa> <gba1_region_bed> <output_dir> [threads]" >&2
    exit 1
}

if [[ $# -lt 5 || $# -gt 6 ]]; then
    usage
fi

INPUT_FILE="$1"
FASTQ_DIR="$2"
MASKED_REFERENCE="$3"
GBA1_REGION_BED="$4"
OUT_DIR="$5"
THREADS="${6:-4}"

[[ -f "$INPUT_FILE" ]]      || { echo "Error: input file '$INPUT_FILE' not found." >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]]       || { echo "Error: fastq directory '$FASTQ_DIR' not found." >&2; exit 1; }
[[ -f "$MASKED_REFERENCE" ]]       || { echo "Error: GBA1LP-masked reference '$MASKED_REFERENCE' not found." >&2; exit 1; }
[[ -f "$GBA1_REGION_BED" ]] || { echo "Error: GBA1 region bed '$GBA1_REGION_BED' not found." >&2; exit 1; }

mkdir -p "$OUT_DIR"

sanitize_id() {
    local id="$1"
    id="${id//-/_}"
    id="$(echo "$id" | tr -s ' \t' '_')"
    echo "$id"
}


# run_bwa_alignment <sample_id> <r1_fastq> <r2_fastq> <masked_reference> <output_dir> <threads>
# Aligns R1/R2 against the GBA1LP-masked reference, sorts and indexes the BAM.
run_bwa_alignment() {
    local sample_id="$1"
    local r1="$2"
    local r2="$3"
    local reference="$4"
    local out_dir="$5"
    local threads="$6"

    local out_bam="${out_dir}/${sample_id}.sorted.bam"

    echo "Aligning sample '${sample_id}' against masked reference" >&2
    echo "  R1: $r1" >&2
    echo "  R2: $r2" >&2

    bwa mem -t "$threads" -R "@RG\tID:${sample_id}\tSM:${sample_id}" "$reference" "$r1" "$r2" \
        | samtools sort -@ "$threads" -o "$out_bam" -

    samtools index "$out_bam"

    echo "  -> BAM written to $out_bam" >&2
    echo "$out_bam"
}

# check_gba1_depth <sample_id> <bam> <gba1_region_bed> <min_depth>
# Computes mean read depth over the GBA1 region and flags samples below
# the LONG-NEXT minimum depth threshold (500x) as possible ADO events.
check_gba1_depth() {
    local sample_id="$1"
    local bam="$2"
    local gba1_region_bed="$3"
    local min_depth="$4"

    local mean_depth
    mean_depth="$(samtools depth -a -b "$gba1_region_bed" "$bam" \
        | awk '{sum+=$3; n++} END {if (n>0) printf "%.1f", sum/n; else print "0"}')"

    echo "Sample '${sample_id}': mean GBA1 depth = ${mean_depth}x" >&2
    echo "Sample '${sample_id}': mean GBA1 depth = ${mean_depth}x" >> ./log_GBA1_depth.txt

    if awk -v d="$mean_depth" -v m="$min_depth" 'BEGIN{exit !(d<m)}'; then
        echo "WARNING: sample '${sample_id}' mean depth (${mean_depth}x) is below the ${min_depth}x threshold." >&2
        echo "         Possible allele dropout (ADO)/complex recombination - consider orthogonal confirmation (e.g. MLPA)." >&2
    fi

    echo "$mean_depth"
}

# call_variants_gatk4 <sample_id> <bam> <reference> <output_dir>
# Calls variants following GATK4 best practices (HaplotypeCaller).
call_variants_gatk4() {
    local sample_id="$1"
    local bam="$2"
    local reference="$3"
    local gba1_region_bed="$4"
    local out_dir="$5"

    local out_vcf="${out_dir}/${sample_id}.g.vcf.gz"

    echo "Calling variants for sample '${sample_id}' with GATK4 HaplotypeCaller..." >&2

    gatk HaplotypeCaller \
        -R "$reference" \
        -I "$bam" \
        -O "$out_vcf" \
        -L "$gba1_region_bed" \
        -ip 100 \
        --emit-ref-confidence GVCF >&2

    echo "  -> gVCF written to $out_vcf" >&2
    echo "$out_vcf"
}

# filter_variants <sample_id> <gvcf> <reference> <output_dir>
# Genotypes the gVCF, applies GATK4 hard filters (SNPs/indels handled
# separately per GATK recommendations), then flags variants with
# suspiciously low VAF or low depth - the signature of pseudogene
# misalignment/allele-dropout described in the LONG-NEXT paper (e.g. the
# L444P false positive with 21-26% VAF from GBA1LP-derived reads, and the
# false negatives caused by low coverage/quality in GBA1 exons 10-11).
# Produces a final filtered, annotated VCF with FILTER column populated:
#   PASS            - passed all filters
#   GATK hard-filter labels (e.g. QD2, FS60, MQ40, ...) - failed GATK4
#                      hard filter thresholds
#   LOW_DEPTH        - per-variant depth below MIN_VARIANT_DP
#   LOW_VAF          - VAF below MIN_VAF (likely artifact/ADO)
#   HET_VAF_OUTLIER  - heterozygous call with VAF outside expected
#                      0.30-0.70 window (possible pseudogene contamination
#                      or unbalanced allele dropout)
filter_variants() {
    local sample_id="$1"
    local gvcf="$2"
    local reference="$3"
    local out_dir="$4"

    local genotyped_vcf="${out_dir}/${sample_id}.genotyped.vcf.gz"
    local snps_vcf="${out_dir}/${sample_id}.snps.vcf.gz"
    local indels_vcf="${out_dir}/${sample_id}.indels.vcf.gz"
    local snps_filtered="${out_dir}/${sample_id}.snps.hardfiltered.vcf.gz"
    local indels_filtered="${out_dir}/${sample_id}.indels.hardfiltered.vcf.gz"
    local merged_vcf="${out_dir}/${sample_id}.hardfiltered.vcf.gz"
    local final_vcf="${out_dir}/${sample_id}.final_filtered.vcf.gz"

    echo "Genotyping gVCF for sample '${sample_id}'..." >&2
    gatk GenotypeGVCFs -R "$reference" -V "$gvcf" -O "$genotyped_vcf" >&2

    echo "Applying GATK4 hard filters (SNPs and indels separately)..." >&2

    gatk SelectVariants -R "$reference" -V "$genotyped_vcf" --select-type-to-include SNP -O "$snps_vcf" >&2
    gatk SelectVariants -R "$reference" -V "$genotyped_vcf" --select-type-to-include INDEL -O "$indels_vcf" >&2

    gatk VariantFiltration -R "$reference" -V "$snps_vcf"         --filter-expression "QD < 2.0"              --filter-name "QD2"         --filter-expression "FS > 60.0"              --filter-name "FS60"         --filter-expression "MQ < 40.0"              --filter-name "MQ40"         --filter-expression "MQRankSum < -12.5"      --filter-name "MQRankSum-12.5"         --filter-expression "ReadPosRankSum < -8.0"  --filter-name "ReadPosRankSum-8"         --filter-expression "SOR > 3.0"              --filter-name "SOR3"         -O "$snps_filtered" >&2

    gatk VariantFiltration -R "$reference" -V "$indels_vcf"         --filter-expression "QD < 2.0"                --filter-name "QD2"         --filter-expression "FS > 200.0"              --filter-name "FS200"         --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20"         --filter-expression "SOR > 10.0"              --filter-name "SOR10"         -O "$indels_filtered" >&2

    echo "Merging hard-filtered SNPs and indels..." >&2
    gatk MergeVcfs -I "$snps_filtered" -I "$indels_filtered" -O "$merged_vcf" >&2

    echo "Flagging low-VAF / low-depth / heterozygous-VAF-outlier variants (ADO/pseudogene misalignment signature)..." >&2
    bcftools filter -e "FORMAT/DP < ${MIN_VARIANT_DP}" -s LOW_DEPTH -m+ "$merged_vcf" -Oz         | bcftools filter -e "(FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) < ${MIN_VAF}" -s LOW_VAF -m+ -Oz         | bcftools filter -e "GT[0]=\"het\" && ((FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) < ${HET_VAF_LOW} || (FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) > ${HET_VAF_HIGH})" -s HET_VAF_OUTLIER -m+ -Oz         -o "$final_vcf"

    bcftools index -t "$final_vcf" >&2

    echo "  -> Final filtered VCF written to $final_vcf" >&2
    echo "$final_vcf"
}


while IFS= read -r raw_id || [[ -n "$raw_id" ]]; do
    [[ -z "$raw_id" ]] && continue
    sanitized_id="$(sanitize_id "$raw_id")"

    matches=()
    while IFS= read -r -d '' f; do
        matches+=("$f")
    done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "${sanitized_id}_*" -print0)

    n_matches=${#matches[@]}

    if [[ "$n_matches" -ne 2 ]]; then
        echo "WARNING: expected exactly 2 matches for ID '${raw_id}' (sanitized: '${sanitized_id}'), found ${n_matches}. Skipping." >&2
        printf '%s\n' "${matches[@]}" >&2
        continue
    fi

    r1=""
    r2=""
    for f in "${matches[@]}"; do
        base="$(basename "$f")"
        if [[ "$base" == *_R1_* || "$base" == *_R1.* ]]; then
            r1="$f"
        elif [[ "$base" == *_R2_* || "$base" == *_R2.* ]]; then
            r2="$f"
        fi
    done

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: could not unambiguously assign R1/R2 for ID '${raw_id}' (sanitized: '${sanitized_id}'). Matches were:" >&2
        printf '%s\n' "${matches[@]}" >&2
        continue
    fi

    bam="$(run_bwa_alignment "$sanitized_id" "$r1" "$r2" "$MASKED_REFERENCE" "$OUT_DIR" "$THREADS")"
    check_gba1_depth "$sanitized_id" "$bam" "$GBA1_REGION_BED" "$MIN_DEPTH" >/dev/null
    gvcf="$(call_variants_gatk4 "$sanitized_id" "$bam" "$MASKED_REFERENCE" "$GBA1_REGION_BED" "$OUT_DIR")"
    filter_variants "$sanitized_id" "$gvcf" "$MASKED_REFERENCE" "$OUT_DIR" >/dev/null

done < "$INPUT_FILE"
