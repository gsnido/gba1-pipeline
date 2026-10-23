#!/usr/bin/env bash
#
# PREPARE_REF.sh
#
# Implements the masking of GBA1LP from the LONG-NEXT method (Cuconato et al.,
# Parkinsonism Relat Disord 2025) for GBA1 sequencing analysis:
#
#    - Masks the GBA1LP pseudogene sequence on the hg38 reference using
#      bedtools maskfasta, so reads align univocally to GBA1
#      (paper section 2.3.3).
#
# Usage:
#   ./PREPARE_REF.sh <reference.fa> <gba1lp_bed> <gba1_region_bed> [threads]
#
# Arguments:
#   reference.fa        : unmasked hg38 reference FASTA (must be bwa- and
#                          samtools-indexed, or indexes will be created)
#   gba1lp_bed          : BED file with GBA1LP pseudogene coordinates (hg38),
#                          used to mask the reference (paper 2.3.3)
#   gba1_region_bed     : BED file with GBA1 gene coordinates (hg38), used
#                          for the 500x read-depth QC check (paper 2.3.3)
#   threads (optional)  : number of threads for bwa/samtools/gatk (default 4)
#
# Requirements: bwa, samtools, bedtools, gatk (GATK4), bcftools must be in PATH.

set -euo pipefail


usage() {
    echo "Usage: $0 <reference.fa> <gba1lp_bed>" >&2
    exit 1
}

if [[ $# -lt 2 || $# -gt 2 ]]; then
    usage
fi

REFERENCE="$1"
GBA1LP_BED="$2"

[[ -f "$REFERENCE" ]]       || { echo "Error: reference '$REFERENCE' not found." >&2; exit 1; }
[[ -f "$GBA1LP_BED" ]]      || { echo "Error: GBA1LP bed '$GBA1LP_BED' not found." >&2; exit 1; }


# mask_reference <reference.fa> <gba1lp_bed> <output_dir>
# Masks the GBA1LP pseudogene sequence on the reference genome using
# bedtools maskfasta, so that reads cannot align to it (paper 2.3.3).
# Prints the path to the masked FASTA on stdout, indexed for bwa/samtools.
mask_reference() {
    local reference="$1"
    local gba1lp_bed="$2"

    local masked_fasta="${reference%.*}.GBA1LP_masked.fa"

    if [[ ! -f "$masked_fasta" ]]; then
        echo "1) Masking GBA1LP pseudogene sequence on reference (bedtools maskfasta)..." >&2
        bedtools maskfasta -fi "$reference" -bed "$gba1lp_bed" -fo "$masked_fasta" >&2

        echo "2) Indexing masked reference (bwa + samtools/faidx)..." >&2
        bwa index "$masked_fasta" >&2
        samtools faidx "$masked_fasta" >&2
        gatk CreateSequenceDictionary -R "$masked_fasta" >&2
    else
        echo "Masked reference already exists: $masked_fasta" >&2
    fi

    echo "$masked_fasta"
}


mask_reference "$REFERENCE" "$GBA1LP_BED"

