#!/usr/bin/env bash
#
# LONG-NEXT-SV.sh
#
# Structural-variant (SV) and recombinant-allele screening for GBA1 long-range
# PCR amplicons sequenced with short reads, as a complement to LONG-NEXT.sh
# (Cuconato et al., Parkinsonism Relat Disord 2025). It re-uses the
# <sample>.sorted.bam files produced by LONG-NEXT.sh (aligned to the
# GBA1LP-masked reference) and does not realign reads.
#
# Per sample, four complementary lines of evidence are produced:
#
#   1) Delly (short-read mode) calls of DEL/DUP/INV/INS/BND restricted to the
#      amplicon, with extra FILTER labels for amplicon-specific artefacts.
#   2) A caller-independent scan of read alignments for mid-sized events that
#      sit between HaplotypeCaller indels and Delly paired-end calls:
#      CIGAR deletions/insertions >= MIN_EVENT_LEN and split-read (SA tag)
#      junctions. Known events such as c.1263_1317del (Rec(delta55)) are
#      labelled.
#   3) Gene-conversion / recombinant tracts: allele counts at GBA1/GBA1LP
#      paralogous sequence variants (PSVs). Because GBA1LP is masked,
#      pseudogene-derived sequence carried on a recombinant allele (e.g.
#      RecNciI, gene conversions, GBA1-GBA1LP fusions amplified by the
#      reverse primer that also anneals in the GBA1LP 3'-UTR) maps to GBA1
#      and appears as a run of consecutive PSVs carrying the GBA1LP allele.
#   4) A binned depth profile over the amplicon, normalised per sample and,
#      when >= MIN_COHORT samples are analysed, relative to the cohort median
#      profile (flags candidate deletions/duplications or allele dropout).
#
# Usage:
#   ./LONG-NEXT-SV.sh <input_id_list.txt> <bam_dir> <masked_reference.fa> <psv_table.tsv> <output_dir> [threads]
#
# Arguments:
#   input_id_list.txt : one sample ID per line (same list given to LONG-NEXT.sh)
#   bam_dir           : LONG-NEXT.sh output directory containing
#                       <sanitized_id>.sorted.bam (+ .bai)
#   masked_reference  : the GBA1LP-masked reference used by LONG-NEXT.sh
#                       (with .fai)
#   psv_table.tsv     : GBA1/GBA1LP PSV table
#                       (bedfiles/GBA1_GBA1LP_PSVs.hg38.tsv)
#   output_dir        : output directory (may be the LONG-NEXT.sh output dir)
#   threads           : optional, default 4
#
# Requirements: samtools, bcftools, delly (>= 1.2, uses `delly sr`) and awk.
#
# Limitations (read before interpreting results):
#   - Only events contained in the 6.5 kb amplicon and compatible with the
#     primer design can be seen. Deletions/reciprocal recombinants that remove
#     a primer-binding site are not amplified and appear as allele dropout
#     (apparent homozygosity), not as an SV. Copy number cannot be inferred
#     from PCR amplicons; use MLPA or long-read sequencing for confirmation.
#   - PCR favours shorter alleles, so allele fractions of deletion alleles may
#     be inflated and duplications deflated.
#   - Every flag here is a screening result requiring orthogonal confirmation.

set -euo pipefail

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
# LONG-NEXT amplicon (hg38), defined by the published primers:
#   GBA1_6.5kb_Forward TCCTAAAGTTGTCACCCATACATG (anneals at chr1:155241392-155241415, - strand)
#   GBA1_6.5kb_Reverse TAGTCACAGACAGCGTGTGAGC   (anneals at chr1:155234915-155234936, + strand)
AMP_CHROM="chr1"
AMP_START=155234915     # 1-based, inclusive
AMP_END=155241415       # 1-based, inclusive
EDGE_MARGIN=300         # breakpoints this close to amplicon ends are flagged (AMPLICON_EDGE)

MIN_MAPQ=20             # read MAPQ for the alignment scan, PSV counting and depth
MIN_BASEQ=20            # base quality for PSV counting
MIN_EVENT_LEN=30        # minimum CIGAR indel length reported by the alignment scan
MIN_EVENT_READS=10      # minimum supporting reads for an alignment-scan event
MIN_EVENT_FRAC=0.05     # minimum supporting reads / local depth for an event
JUNCTION_BIN=10         # split-read breakpoints are clustered to this resolution (bp)

MIN_SV_FRAC=0.05        # Delly: min (DV+RV)/(DR+DV+RR+RV) before LOW_SUPPORT flag

PSV_MIN_DP=50           # minimum depth for a PSV to be informative
PSV_MIN_FRAC=0.15       # GBA1LP-allele fraction at which a PSV counts as "pseudogene-like"
                        # (same value as MIN_VAF in LONG-NEXT.sh)
MIN_TRACT_PSVS=2        # minimum consecutive pseudogene-like PSVs to report a tract
TRACT_MAX_JUMP=0.25     # a PSV whose fraction differs from the running tract mean by more
                        # than this starts a new tract (keeps e.g. a homozygous common SNP
                        # from being merged with a heterozygous recombinant tract)
HET_FRAC_LOW=0.30       # tract mean fraction window consistent with heterozygosity
HET_FRAC_HIGH=0.70
HOM_FRAC=0.85           # tract mean fraction consistent with homozygosity / ADO

BIN_SIZE=100            # depth-profile bin size (bp)
MIN_COHORT=3            # minimum samples for the cohort-relative depth comparison
DEPTH_RATIO_LOW=0.60    # cohort-relative depth ratio flags
DEPTH_RATIO_HIGH=1.50
MIN_FLAG_BINS=3         # consecutive bins required for a depth flag

# GBA1 exons, NM_000157.4 (MANE Select), hg38, 1-based inclusive. GBA1 is on
# the minus strand, so exon 1 has the highest coordinates.
GBA1_EXONS="1:155241086-155241249 2:155240630-155240717 3:155239886-155240077 4:155239616-155239762 5:155238517-155238650 6:155238134-155238306 7:155237341-155237578 8:155236245-155236469 9:155235681-155235844 10:155235195-155235311 11:155234452-155235100"

# Known length-changing recombinant event: c.1263_1317del (55 bp, exon 9),
# left-normalised hg38 representation chr1:155235750-155235804 (+/- 5 bp).
DEL55_START=155235750
DEL55_LEN=55

usage() {
    echo "Usage: $0 <input_id_list.txt> <bam_dir> <masked_reference.fa> <psv_table.tsv> <output_dir> [threads]" >&2
    exit 1
}

if [[ $# -lt 5 || $# -gt 6 ]]; then
    usage
fi

INPUT_FILE="$1"
BAM_DIR="$2"
MASKED_REFERENCE="$3"
PSV_TABLE="$4"
OUT_DIR="$5"
THREADS="${6:-4}"

[[ -f "$INPUT_FILE" ]]              || { echo "Error: input file '$INPUT_FILE' not found." >&2; exit 1; }
[[ -d "$BAM_DIR" ]]                 || { echo "Error: BAM directory '$BAM_DIR' not found." >&2; exit 1; }
[[ -f "$MASKED_REFERENCE" ]]        || { echo "Error: masked reference '$MASKED_REFERENCE' not found." >&2; exit 1; }
[[ -f "${MASKED_REFERENCE}.fai" ]]  || { echo "Error: reference index '${MASKED_REFERENCE}.fai' not found." >&2; exit 1; }
[[ -f "$PSV_TABLE" ]]               || { echo "Error: PSV table '$PSV_TABLE' not found." >&2; exit 1; }

for tool in samtools bcftools delly awk; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: '$tool' not found in PATH." >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
SHARED_DIR="${OUT_DIR}/sv_shared"
mkdir -p "$SHARED_DIR"

AMP_REGION="${AMP_CHROM}:${AMP_START}-${AMP_END}"

sanitize_id() {
    local id="$1"
    id="${id//-/_}"
    id="$(echo "$id" | tr -s ' \t' '_')"
    echo "$id"
}

# exon_label <pos>  ->  "exon9" / "intron9" / "3'flank" / "5'flank"
# (awk-side copy of this logic is used in the tables; see EXON_AWK below)
EXON_AWK='
function exon_label(p,   i, n, a, b, e, s, idx) {
    n = split(EXONS, a, " ")
    for (i = 1; i <= n; i++) {
        split(a[i], b, /[:-]/); idx[i] = b[1]; s[i] = b[2]; e[i] = b[3]
        if (p >= s[i] && p <= e[i]) return "exon" b[1]
    }
    if (p > s[1]) return "5prime_flank"
    if (p < s[n]) return "3prime_flank"
    for (i = 1; i < n; i++) if (p < s[i] && p > e[i+1]) return "intron" idx[i]
    return "NA"
}
function median(arr, m,   i, j, t, w) {
    for (i = 1; i <= m; i++) w[i] = arr[i]
    for (i = 2; i <= m; i++) { t = w[i]; for (j = i - 1; j >= 1 && w[j] > t; j--) w[j + 1] = w[j]; w[j + 1] = t }
    if (m == 0) return 0
    return (m % 2) ? w[(m + 1) / 2] : (w[m / 2] + w[m / 2 + 1]) / 2
}'

# ----------------------------------------------------------------------------
# Shared resources (computed once)
# ----------------------------------------------------------------------------

# Delly exclusion file: every reference base outside the padded amplicon, so
# Delly only processes the target (and is fast on a full hg38 reference).
build_delly_exclude() {
    local excl="${SHARED_DIR}/delly_exclude.bed"
    local pad=1000
    awk -v c="$AMP_CHROM" -v s="$((AMP_START - pad))" -v e="$((AMP_END + pad))" 'BEGIN{OFS="\t"}
        $1 != c { print $1, 0, $2; next }
        { if (s > 1) print $1, 0, s - 1; if (e < $2) print $1, e, $2 }' \
        "${MASKED_REFERENCE}.fai" > "$excl"
    echo "$excl"
}

# Converts the PSV table into a bcftools regions file and checks that the
# GBA1 alleles agree with the reference (guards against assembly mismatch).
prepare_psv_sites() {
    local sites="${SHARED_DIR}/psv_sites.tsv"
    grep -v '^#' "$PSV_TABLE" | awk -F'\t' 'BEGIN{OFS="\t"} NF >= 4 {print $1, $2}' | sort -k1,1 -k2,2n -u > "$sites"

    local mismatches
    mismatches="$(grep -v '^#' "$PSV_TABLE" | while IFS=$'\t' read -r c p ref _rest; do
            b="$(samtools faidx "$MASKED_REFERENCE" "${c}:${p}-${p}" | tail -n +2 | tr 'a-z' 'A-Z')"
            [[ "$b" == "$ref" ]] || echo "${c}:${p}"
        done | wc -l)"
    if [[ "$mismatches" -gt 0 ]]; then
        echo "Error: ${mismatches} PSV(s) have a GBA1 allele that does not match the reference." >&2
        echo "       Is the reference hg38 with 'chr' contig names?" >&2
        exit 1
    fi
    echo "$sites"
}

# ----------------------------------------------------------------------------
# Per-sample analyses
# ----------------------------------------------------------------------------

# call_delly <sample_id> <bam> <exclude_bed> <out_dir>
# Runs Delly short-read discovery, keeps records with POS inside the amplicon
# and adds FILTER labels:
#   AMPLICON_EDGE - a breakpoint lies within EDGE_MARGIN bp of an amplicon end
#                   (read-end/primer artefacts are common there)
#   LOW_SUPPORT   - supporting fraction below MIN_SV_FRAC
#   (Delly's own LowQual label is kept)
call_delly() {
    local sample_id="$1"
    local bam="$2"
    local excl="$3"
    local out_dir="$4"

    local raw_bcf="${out_dir}/${sample_id}.delly.raw.bcf"
    local out_vcf="${out_dir}/${sample_id}.delly.sv.vcf.gz"

    echo "  [delly] calling SVs..." >&2
    delly sr -g "$MASKED_REFERENCE" -x "$excl" -q "$MIN_MAPQ" -h "$THREADS" \
        -o "$raw_bcf" "$bam" > "${out_dir}/${sample_id}.delly.log" 2>&1 || {
        echo "  [delly] WARNING: Delly failed (see ${out_dir}/${sample_id}.delly.log)." >&2
        return 0
    }

    # With no candidates Delly writes a header without sample columns.
    if [[ ! -s "$raw_bcf" ]] || [[ "$(bcftools view -H "$raw_bcf" | wc -l)" -eq 0 ]]; then
        echo "  [delly] no SV candidates." >&2
        rm -f "$raw_bcf" "${raw_bcf}.csi"
        return 0
    fi

    local lo_edge=$((AMP_START + EDGE_MARGIN)) hi_edge=$((AMP_END - EDGE_MARGIN))
    bcftools view -i "CHROM==\"${AMP_CHROM}\" && POS>=${AMP_START} && POS<=${AMP_END}" "$raw_bcf" -Ou \
        | bcftools filter -m+ -s AMPLICON_EDGE \
            -e "POS<${lo_edge} || POS>${hi_edge} || (INFO/END<${lo_edge} && INFO/SVTYPE!=\"BND\") || (INFO/END>${hi_edge} && INFO/SVTYPE!=\"BND\")" -Ou \
        | bcftools filter -m+ -s LOW_SUPPORT \
            -e "(FORMAT/DV[0]+FORMAT/RV[0]) / (FORMAT/DR[0]+FORMAT/DV[0]+FORMAT/RR[0]+FORMAT/RV[0]+0.0001) < ${MIN_SV_FRAC}" \
            -Oz -o "$out_vcf"
    bcftools index -t -f "$out_vcf"
    rm -f "$raw_bcf" "${raw_bcf}.csi"

    echo "  [delly] $(bcftools view -H "$out_vcf" | wc -l) record(s) in amplicon, $(bcftools view -H -f PASS "$out_vcf" | wc -l) PASS -> $out_vcf" >&2
}

# scan_alignment_events <sample_id> <bam> <out_dir>
# Caller-independent evidence table:
#   CIGAR_DEL / CIGAR_INS - indels >= MIN_EVENT_LEN inside read alignments
#   SPLIT_JUNCTION        - breakpoint pairs joined by a supplementary alignment
#                           (SA tag), with the implied orientation
# Columns: sample, type, chrom, start, end, length, support, depth, fraction,
#          region, label, flag
scan_alignment_events() {
    local sample_id="$1"
    local bam="$2"
    local out_dir="$3"

    local raw="${out_dir}/${sample_id}.events.raw.tsv"
    local out_tsv="${out_dir}/${sample_id}.alignment_events.tsv"

    # 1) extract events from primary alignments
    samtools view -q "$MIN_MAPQ" -F 0xF04 "$bam" "$AMP_REGION" \
        | awk -v minlen="$MIN_EVENT_LEN" -v jb="$JUNCTION_BIN" 'BEGIN{OFS="\t"}
        function reflen(cig,   n, l, op, t) {
            t = 0
            while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                l = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1)
                if (op ~ /[MDN=X]/) t += l
                cig = substr(cig, RLENGTH + 1)
            }
            return t
        }
        {
            pos = $4; cig = $6; r = pos
            # CIGAR indels
            c = cig
            while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                l = substr(c, 1, RLENGTH - 1) + 0; op = substr(c, RLENGTH, 1)
                if (op == "D" && l >= minlen) print "CIGAR_DEL", $3, r, r + l - 1, l
                if (op == "I" && l >= minlen) print "CIGAR_INS", $3, r - 1, r - 1, l
                if (op ~ /[MDN=X]/) r += l
                c = substr(c, RLENGTH + 1)
            }
            # split-read junctions: breakpoint on the clipped side of the primary
            sa = ""
            for (i = 12; i <= NF; i++) if ($i ~ /^SA:Z:/) { sa = substr($i, 6); break }
            if (sa == "") next
            split(sa, s1, ";"); split(s1[1], f, ",")
            if (f[1] != $3) next
            pend = pos + reflen(cig) - 1
            lead = (cig ~ /^[0-9]+[SH]/); trail = (cig ~ /[0-9]+[SH]$/)
            bp1 = (lead && !trail) ? pos : pend
            send = f[2] + reflen(f[4]) - 1
            slead = (f[4] ~ /^[0-9]+[SH]/); strail = (f[4] ~ /[0-9]+[SH]$/)
            bp2 = (slead && !strail) ? f[2] : send
            strand = (int($2 / 16) % 2 ? "-" : "+")
            ori = (strand == f[3]) ? "same_strand" : "inverted"
            a = (bp1 < bp2) ? bp1 : bp2; b = (bp1 < bp2) ? bp2 : bp1
            a = int(a / jb) * jb; b = int(b / jb) * jb
            if (b - a < minlen) next
            print "SPLIT_JUNCTION_" ori, $3, a, b, b - a
        }' | sort | uniq -c \
        | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $5, $6, $1}' > "$raw"

    # 2) attach local depth (at the event start), fraction, region and labels
    {
        printf "sample\ttype\tchrom\tstart\tend\tlength\tsupport\tdepth\tfraction\tregion\tlabel\tflag\n"
        while IFS=$'\t' read -r type chrom start end len support; do
            [[ "$support" -lt "$MIN_EVENT_READS" ]] && continue
            depth="$(samtools depth -a -J -Q "$MIN_MAPQ" -r "${chrom}:${start}-${start}" "$bam" | awk '{print $3+0}')"
            depth="${depth:-0}"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample_id" "$type" "$chrom" "$start" "$end" "$len" "$support" "$depth"
        done < "$raw" \
        | awk -F'\t' -v EXONS="$GBA1_EXONS" -v mf="$MIN_EVENT_FRAC" -v d55="$DEL55_START" -v l55="$DEL55_LEN" \
              -v as="$AMP_START" -v ae="$AMP_END" -v em="$EDGE_MARGIN" 'BEGIN{OFS="\t"}'"$EXON_AWK"'
        {
            frac = ($8 > 0) ? $7 / $8 : 0
            if (frac < mf) next
            reg = exon_label($4); if ($5 != $4) reg = reg "-" exon_label($5)
            label = "."
            if (($2 == "CIGAR_DEL" || $2 == "SPLIT_JUNCTION_same_strand") && $6 >= l55 - 10 && $6 <= l55 + 10 && $4 >= d55 - 10 && $4 <= d55 + 10)
                label = "c.1263_1317del_(Rec_delta55)"
            flag = "PASS"
            if ($4 < as + em || $5 > ae - em) flag = "AMPLICON_EDGE"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, frac, reg, label, flag
        }' | sort -t$'\t' -k4,4n
    } > "$out_tsv"
    rm -f "$raw"

    echo "  [scan] $(($(wc -l < "$out_tsv") - 1)) alignment event(s) -> $out_tsv" >&2
}

# analyse_psvs <sample_id> <bam> <psv_sites> <out_dir>
# Counts GBA1 and GBA1LP alleles at every PSV and segments runs of
# pseudogene-like PSVs into tracts.
#   <sample>.psv_counts.tsv : per-PSV counts and GBA1LP-allele fraction
#   <sample>.psv_tracts.tsv : tracts of >= MIN_TRACT_PSVS consecutive
#                             informative PSVs with fraction >= PSV_MIN_FRAC
# Tract classes:
#   HET_LIKE        mean fraction within HET_FRAC_LOW-HET_FRAC_HIGH
#   HOM_LIKE        mean fraction >= HOM_FRAC (homozygous or allele dropout)
#   LOW_FRACTION    mean fraction < HET_FRAC_LOW (pseudogene co-amplification,
#                   contamination or unbalanced amplification)
#   UNBALANCED      remaining values
# PSV status in psv_counts.tsv: GBA1_LIKE, GBA1LP_LIKE, LOW_DEPTH, or
# OTHER_ALLELE (enough reads, but mostly a third allele/indel at the site).
# "reaches_3prime_end"=yes means the tract includes the 3'-most informative
# PSV, as expected for a GBA1-GBA1LP fusion / RecNciI-type allele.
analyse_psvs() {
    local sample_id="$1"
    local bam="$2"
    local sites="$3"
    local out_dir="$4"

    local counts="${out_dir}/${sample_id}.psv_counts.tsv"
    local tracts="${out_dir}/${sample_id}.psv_tracts.tsv"

    {
        printf "sample\tchrom\tpos\tregion\tGBA1_allele\tGBA1LP_allele\tGBA1_count\tGBA1LP_count\tother_count\tdepth\tGBA1LP_fraction\tstatus\tsource\tlabel\n"
        bcftools mpileup -f "$MASKED_REFERENCE" -R "$sites" -a AD -d 1000000 \
            -q "$MIN_MAPQ" -Q "$MIN_BASEQ" --no-BAQ -I "$bam" 2>/dev/null \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' \
            | awk -F'\t' -v OFS='\t' -v sid="$sample_id" -v mindp="$PSV_MIN_DP" -v mf="$PSV_MIN_FRAC" \
                  -v EXONS="$GBA1_EXONS" "$EXON_AWK"'
              NR == FNR { if ($0 !~ /^#/) { lp[$1":"$2] = $4; src[$1":"$2] = $5; lab[$1":"$2] = $6 }; next }
              {
                  k = $1":"$2; if (!(k in lp)) next
                  n = split($4, alt, ","); split($5, ad, ",")
                  g = ad[1] + 0; p = 0; tot = 0
                  for (i = 1; i <= n + 1; i++) tot += ad[i]
                  for (i = 1; i <= n; i++) if (alt[i] == lp[k]) p = ad[i + 1] + 0
                  o = tot - g - p
                  f = (g + p > 0) ? p / (g + p) : 0
                  st = (g + p >= mindp) ? ((f >= mf) ? "GBA1LP_LIKE" : "GBA1_LIKE") : ((tot >= mindp) ? "OTHER_ALLELE" : "LOW_DEPTH")
                  printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%s\t%s\t%s\n", sid, $1, $2, exon_label($2), $3, lp[k], g, p, o, tot, f, st, src[k], lab[k]
              }' "$PSV_TABLE" -
    } > "$counts"

    # Segment tracts. PSVs are walked in genomic order; LOW_DEPTH and
    # OTHER_ALLELE PSVs neither extend nor break a tract.
    awk -F'\t' -v OFS='\t' -v sid="$sample_id" -v minn="$MIN_TRACT_PSVS" \
        -v hl="$HET_FRAC_LOW" -v hh="$HET_FRAC_HIGH" -v hom="$HOM_FRAC" -v jump="$TRACT_MAX_JUMP" -v EXONS="$GBA1_EXONS" "$EXON_AWK"'
        function flush() {
            if (n >= minn) {
                m = sum / n
                cls = (m >= hom) ? "HOM_LIKE" : ((m >= hl && m <= hh) ? "HET_LIKE" : ((m < hl) ? "LOW_FRACTION" : "UNBALANCED"))
                printf "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%s\t%s\t%s\n", sid, chr, s, e, exon_label(e), exon_label(s), n, m, mn, mx, cls, (s == first_inf ? "yes" : "no"), (labs == "" ? "." : labs)
            }
            n = 0; sum = 0; labs = ""
        }
        BEGIN { print "sample", "chrom", "start", "end", "region_5prime", "region_3prime", "n_psvs", "mean_GBA1LP_fraction", "min_fraction", "max_fraction", "class", "reaches_3prime_end", "markers" }
        NR == 1 { next }
        $12 ~ /GBA1/ && first_inf == "" { first_inf = $3 }
        $12 != "GBA1LP_LIKE" && $12 != "GBA1_LIKE" { next }
        $12 == "GBA1LP_LIKE" {
            if (n > 0 && ($11 - sum / n > jump || sum / n - $11 > jump)) flush()
            if (n == 0) { s = $3; mn = 1; mx = 0; chr = $2 }
            e = $3; n++; sum += $11
            if ($11 < mn) mn = $11; if ($11 > mx) mx = $11
            if ($14 != ".") labs = (labs == "" ? $14 : labs "," $14)
            next
        }
        { flush() }
        END { flush() }' "$counts" > "$tracts"

    local n_inf n_tr
    n_inf="$(awk -F'\t' 'NR > 1 && $12 ~ /GBA1/' "$counts" | wc -l)"
    n_tr="$(($(wc -l < "$tracts") - 1))"
    echo "  [psv] ${n_inf} informative PSVs, ${n_tr} pseudogene-like tract(s) -> $tracts" >&2
}

# depth_profile <sample_id> <bam> <bins_bed> <out_dir>
# Mean depth per BIN_SIZE bin, normalised to the sample's median bin depth.
depth_profile() {
    local sample_id="$1"
    local bam="$2"
    local bins="$3"
    local out_dir="$4"
    local out_tsv="${out_dir}/${sample_id}.depth_profile.tsv"

    samtools bedcov -Q "$MIN_MAPQ" "$bins" "$bam" \
        | awk -F'\t' -v OFS='\t' -v sid="$sample_id" -v EXONS="$GBA1_EXONS" "$EXON_AWK"'
            { c[NR] = $1; s[NR] = $2; e[NR] = $3; d[NR] = $4 / ($3 - $2) }
            END {
                med = median(d, NR)
                print "sample", "chrom", "start", "end", "mean_depth", "norm_depth"
                # start is reported 1-based for consistency with the other tables
                for (i = 1; i <= NR; i++) printf "%s\t%s\t%d\t%d\t%.1f\t%.3f\n", sid, c[i], s[i] + 1, e[i], d[i], (med > 0 ? d[i] / med : 0)
            }' > "$out_tsv"
    echo "$out_tsv"
}

# cohort_depth_flags <out_dir> <profile files...>
# Ratio of each sample's normalised depth to the cohort median per bin, and
# runs of >= MIN_FLAG_BINS bins outside DEPTH_RATIO_LOW/HIGH.
cohort_depth_flags() {
    local out_dir="$1"; shift
    local ratios="${out_dir}/SV_cohort_depth_ratios.tsv"
    local flags="${out_dir}/SV_cohort_depth_flags.tsv"

    awk -F'\t' -v OFS='\t' 'FNR == 1 { next } { print }' "$@" \
        | awk -F'\t' -v OFS='\t' -v EXONS="$GBA1_EXONS" "$EXON_AWK"'
            { key = $2 "\t" $3 "\t" $4; if (!(key in seen)) { seen[key] = 1; order[++nk] = key }
              val[key, ++cnt[key]] = $6; row[++nr] = $0; rk[nr] = key }
            END {
                for (k = 1; k <= nk; k++) {
                    key = order[k]; m = cnt[key]
                    for (i = 1; i <= m; i++) tmp[i] = val[key, i]
                    med[key] = median(tmp, m)
                }
                print "sample", "chrom", "start", "end", "mean_depth", "norm_depth", "cohort_median_norm_depth", "ratio"
                for (i = 1; i <= nr; i++) {
                    split(row[i], f, "\t"); md = med[rk[i]]
                    printf "%s\t%.3f\t%.3f\n", row[i], md, (md > 0 ? f[6] / md : 0)
                }
            }' > "$ratios"

    awk -F'\t' -v OFS='\t' -v lo="$DEPTH_RATIO_LOW" -v hi="$DEPTH_RATIO_HIGH" -v minb="$MIN_FLAG_BINS" -v EXONS="$GBA1_EXONS" "$EXON_AWK"'
        function flush() {
            if (n >= minb) printf "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%.3f\t%s\n", cs, cc, rs, re, exon_label(re), exon_label(rs), n, rsum / n, cur
            n = 0; rsum = 0; cur = ""
        }
        BEGIN { print "sample", "chrom", "start", "end", "region_5prime", "region_3prime", "n_bins", "mean_ratio", "flag" }
        NR == 1 { next }
        {
            st = ($8 < lo) ? "LOW_DEPTH_RATIO" : (($8 > hi) ? "HIGH_DEPTH_RATIO" : "")
            if ($1 != cs || st != cur || (n > 0 && $3 != re + 1)) flush()
            cs = $1; cc = $2
            if (st != "") { if (n == 0) { rs = $3; cur = st }; re = $4; n++; rsum += $8 }
        }
        END { flush() }' "$ratios" > "$flags"

    echo "Cohort depth ratios -> $ratios" >&2
    echo "Cohort depth flags  -> $flags ($(($(wc -l < "$flags") - 1)) region(s))" >&2
}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

echo "LONG-NEXT-SV: amplicon ${AMP_REGION}" >&2
EXCLUDE_BED="$(build_delly_exclude)"
PSV_SITES="$(prepare_psv_sites)"

BINS_BED="${SHARED_DIR}/amplicon_bins.bed"
awk -v c="$AMP_CHROM" -v s="$((AMP_START - 1))" -v e="$AMP_END" -v b="$BIN_SIZE" 'BEGIN{OFS="\t"
    for (x = s; x < e; x += b) print c, x, (x + b < e ? x + b : e) }' > "$BINS_BED"

SUMMARY="${OUT_DIR}/SV_summary.tsv"
printf "sample\tmean_amplicon_depth\tdelly_pass\tdelly_pass_types\talignment_events_pass\talignment_event_labels\tpsv_tracts\tpsv_tract_classes\tmarkers_in_tracts\treview\n" > "$SUMMARY"

PROFILES=()
SAMPLE_DIRS=()

while IFS= read -r raw_id || [[ -n "$raw_id" ]]; do
    [[ -z "$raw_id" ]] && continue
    sid="$(sanitize_id "$raw_id")"
    bam="${BAM_DIR}/${sid}.sorted.bam"

    if [[ ! -f "$bam" ]]; then
        echo "WARNING: BAM for '${raw_id}' not found at '${bam}'. Run LONG-NEXT.sh first. Skipping." >&2
        continue
    fi
    [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] || samtools index "$bam"

    echo "Sample '${sid}'" >&2
    sdir="${OUT_DIR}/${sid}_SV"
    mkdir -p "$sdir"

    call_delly "$sid" "$bam" "$EXCLUDE_BED" "$sdir"
    scan_alignment_events "$sid" "$bam" "$sdir"
    analyse_psvs "$sid" "$bam" "$PSV_SITES" "$sdir"
    profile="$(depth_profile "$sid" "$bam" "$BINS_BED" "$sdir")"
    PROFILES+=("$profile")
    SAMPLE_DIRS+=("$sdir")

    # --- per-sample summary line ---
    mean_depth="$(awk -F'\t' 'NR > 1 { s += $5; n++ } END { printf "%.1f", (n ? s / n : 0) }' "$profile")"
    dvcf="${sdir}/${sid}.delly.sv.vcf.gz"
    if [[ -f "$dvcf" ]]; then
        d_pass="$(bcftools view -H -f PASS "$dvcf" | wc -l)"
        d_types="$(bcftools query -i 'FILTER=="PASS"' -f '%INFO/SVTYPE:%POS-%INFO/END\n' "$dvcf" | paste -sd, -)"
    else
        d_pass=0; d_types=""
    fi
    ev="${sdir}/${sid}.alignment_events.tsv"
    e_pass="$(awk -F'\t' 'NR > 1 && $12 == "PASS"' "$ev" | wc -l)"
    e_labels="$(awk -F'\t' 'NR > 1 && $11 != "." { print $11 }' "$ev" | sort -u | paste -sd, -)"
    tr="${sdir}/${sid}.psv_tracts.tsv"
    t_n="$(($(wc -l < "$tr") - 1))"
    t_cls="$(awk -F'\t' 'NR > 1 { printf "%s%s:%s(%s-%s,%s)", (n++ ? "," : ""), $11, $8, $5, $6, ($12 == "yes" ? "3prime_end" : "internal") }' "$tr")"
    t_mk="$(awk -F'\t' 'NR > 1 && $13 != "." { print $13 }' "$tr" | tr ',' '\n' | sort -u | paste -sd, -)"
    review="no"
    if [[ "$d_pass" -gt 0 || "$e_pass" -gt 0 || "$t_n" -gt 0 ]]; then review="yes"; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sid" "$mean_depth" "$d_pass" "${d_types:-.}" "$e_pass" "${e_labels:-.}" \
        "$t_n" "${t_cls:-.}" "${t_mk:-.}" "$review" >> "$SUMMARY"
done < "$INPUT_FILE"

if [[ ${#PROFILES[@]} -ge $MIN_COHORT ]]; then
    cohort_depth_flags "$OUT_DIR" "${PROFILES[@]}"
else
    echo "Fewer than ${MIN_COHORT} samples: skipping cohort-relative depth comparison." >&2
fi

echo "Summary -> $SUMMARY" >&2
