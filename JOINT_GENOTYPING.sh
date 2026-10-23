#!/bin/bash

usage() {
    echo "Usage: $0 <masked_reference.fa> <GVCF_list.list> <gba1_region_bed> <combined.g.vcf.gz> <joint.vcf.gz>" >&2
    exit 1
}

if [[ $# -lt 5 || $# -gt 5 ]]; then
    usage
fi

MASKED_REFERENCE="$1"
GVCF_LIST="$2"
GBA1_REGION_BED="$3"
OUT_FILE="$4"
JOINT_GT="$5"

find /data2/GBA1/results -type f -name "*.g.vcf.gz" > ${GVCF_LIST}

# 1) Combine GVCFs
gatk CombineGVCFs -R ${MASKED_REFERENCE} \
    -V ${GVCF_LIST} \
    -L ${GBA1_REGION_BED} \
    -ip 100 \
    -O ${OUT_FILE}

# 2) Joint genotype
gatk GenotypeGVCFs \
    -R ${MASKED_REFERENCE} \
    -V ${OUT_FILE} \
    -O ${JOINT_GT}

