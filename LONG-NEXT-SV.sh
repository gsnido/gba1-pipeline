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
# Per sample, five complementary lines of evidence are produced:
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
#   4) Read-backed phasing of heterozygous PSVs and the c.1263del55 site,
#      used to name recombinant alleles per haplotype (RecNciI, RecTL,
#      c.1263del55+RecTL, RecA456P, single GBAP1-like variants), following
#      the allele definitions used by Gauchian (Toffoli et al., Commun Biol
#      2022).
#   5) A binned depth profile over the amplicon, normalised per sample and,
#      when >= MIN_COHORT samples are analysed, relative to the cohort median
#      profile (flags candidate deletions/duplications or allele dropout).
#
# Usage:
#   ./LONG-NEXT-SV.sh <input_id_list.txt> <bam_dir> <masked_reference.fa> <psv_table.tsv> <output_dir> [threads] [reliable_sites]
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
#   reliable_sites    : optional file whose first two columns are chrom and
#                       position (e.g. gauchian/data/GBA_SNP_38.txt from your
#                       own copy of Gauchian). PSVs at these positions are
#                       marked in psv_counts.tsv and counted per tract. The
#                       file is not bundled because Gauchian is distributed
#                       under the PolyForm Strict licence.
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

PHASE_MIN_LINKS=10      # minimum fragments linking two heterozygous sites
PHASE_MIN_CIS=0.80      # link is CIS if >= this fraction of fragments are PP/GG,
                        # TRANS if <= 1 - this value, otherwise AMBIGUOUS

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

# Recombinant-allele markers (NM_000157.4, hg38), in the canonical order used
# to build allele names. Legacy names: D409H, L444P, A456P, V460V.
#   c.1263del55  c.1263_1317del (site = deletion start, see above)
#   D448H        c.1342G>C   chr1:155235727 C>G
#   L483P        c.1448T>C   chr1:155235252 A>G
#   A495P        c.1483G>C   chr1:155235217 C>G
#   V499=        c.1497G>C   chr1:155235203 C>G
RECOMBINANT_MARKERS="${DEL55_START}:c.1263del55 155235727:D448H 155235252:L483P 155235217:A495P 155235203:V499="

usage() {
    echo "Usage: $0 <input_id_list.txt> <bam_dir> <masked_reference.fa> <psv_table.tsv> <output_dir> [threads] [reliable_sites]" >&2
    exit 1
}

if [[ $# -lt 5 || $# -gt 7 ]]; then
    usage
fi

INPUT_FILE="$1"
BAM_DIR="$2"
MASKED_REFERENCE="$3"
PSV_TABLE="$4"
OUT_DIR="$5"
THREADS="${6:-4}"
RELIABLE_SITES="${7:-}"

[[ -f "$INPUT_FILE" ]]              || { echo "Error: input file '$INPUT_FILE' not found." >&2; exit 1; }
[[ -d "$BAM_DIR" ]]                 || { echo "Error: BAM directory '$BAM_DIR' not found." >&2; exit 1; }
[[ -f "$MASKED_REFERENCE" ]]        || { echo "Error: masked reference '$MASKED_REFERENCE' not found." >&2; exit 1; }
[[ -f "${MASKED_REFERENCE}.fai" ]]  || { echo "Error: reference index '${MASKED_REFERENCE}.fai' not found." >&2; exit 1; }
[[ -f "$PSV_TABLE" ]]               || { echo "Error: PSV table '$PSV_TABLE' not found." >&2; exit 1; }
[[ -z "$RELIABLE_SITES" || -f "$RELIABLE_SITES" ]] || { echo "Error: reliable-sites file '$RELIABLE_SITES' not found." >&2; exit 1; }

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
        printf "sample\tchrom\tpos\tregion\tGBA1_allele\tGBA1LP_allele\tGBA1_count\tGBA1LP_count\tother_count\tdepth\tGBA1LP_fraction\tstatus\tsource\tlabel\treliable_site\n"
        bcftools mpileup -f "$MASKED_REFERENCE" -R "$sites" -a AD -d 1000000 \
            -q "$MIN_MAPQ" -Q "$MIN_BASEQ" --no-BAQ -I "$bam" 2>/dev/null \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' \
            | awk -F'\t' -v OFS='\t' -v sid="$sample_id" -v mindp="$PSV_MIN_DP" -v mf="$PSV_MIN_FRAC" \
                  -v EXONS="$GBA1_EXONS" -v rel="$RELIABLE_SITES" "$EXON_AWK"'
              BEGIN { if (rel != "") while ((getline line < rel) > 0) { if (line ~ /^#/) continue; n = split(line, rf, /[ \t]+/); if (n >= 2) rs[rf[1] ":" rf[2]] = 1 } }
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
                  printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%s\t%s\t%s\t%s\n", sid, $1, $2, exon_label($2), $3, lp[k], g, p, o, tot, f, st, src[k], lab[k], (rel == "" ? "NA" : ((k in rs) ? "yes" : "no"))
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
                printf "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%s\t%s\t%s\t%s\n", sid, chr, s, e, exon_label(e), exon_label(s), n, m, mn, mx, cls, (s == first_inf ? "yes" : "no"), (labs == "" ? "." : labs), (relmode ? nrel : "NA")
            }
            n = 0; sum = 0; labs = ""; nrel = 0
        }
        BEGIN { print "sample", "chrom", "start", "end", "region_5prime", "region_3prime", "n_psvs", "mean_GBA1LP_fraction", "min_fraction", "max_fraction", "class", "reaches_3prime_end", "markers", "n_reliable_psvs" }
        NR == 1 { relmode = ($15 != "NA"); next }
        $12 ~ /GBA1/ && first_inf == "" { first_inf = $3 }
        $12 != "GBA1LP_LIKE" && $12 != "GBA1_LIKE" { next }
        $12 == "GBA1LP_LIKE" {
            if (n > 0 && ($11 - sum / n > jump || sum / n - $11 > jump)) flush()
            if (n == 0) { s = $3; mn = 1; mx = 0; chr = $2 }
            e = $3; n++; sum += $11
            if ($11 < mn) mn = $11; if ($11 > mx) mx = $11
            if ($14 != ".") labs = (labs == "" ? $14 : labs "," $14)
            if ($15 == "yes") nrel++
            next
        }
        { flush() }
        END { flush() }' "$counts" > "$tracts"

    local n_inf n_tr
    n_inf="$(awk -F'\t' 'NR > 1 && $12 ~ /GBA1/' "$counts" | wc -l)"
    n_tr="$(($(wc -l < "$tracts") - 1))"
    echo "  [psv] ${n_inf} informative PSVs, ${n_tr} pseudogene-like tract(s) -> $tracts" >&2
}

# phase_recombinant_markers <sample_id> <bam> <psv_counts> <out_dir>
# Read-backed phasing of heterozygous PSVs and the c.1263del55 site, and
# naming of recombinant alleles from the phased marker haplotypes.
#
#   1) For every primary read, the allele (GBA1 or GBA1LP) is recorded at each
#      PSV with base quality >= MIN_BASEQ. The c.1263del55 site is treated as
#      one extra biallelic site: GBA1LP-like if the read carries a CIGAR
#      deletion of ~55 bp there, GBA1-like if the read spans the 55 bp
#      without a deletion. Both mates of a pair are merged into one fragment
#      (sites where the mates disagree are dropped).
#   2) Heterozygous sites (GBA1LP fraction >= PSV_MIN_FRAC and < HOM_FRAC)
#      are linked to the next heterozygous site through fragments covering
#      both. A link is CIS if >= PHASE_MIN_CIS of informative fragments
#      carry the same parental type at both sites (GBA1LP-GBA1LP or
#      GBA1-GBA1), TRANS if <= 1 - PHASE_MIN_CIS, otherwise AMBIGUOUS. Links
#      with fewer than PHASE_MIN_LINKS fragments, or AMBIGUOUS links, break
#      the phase block.
#   3) Within each block two haplotypes (A and B) are built; homozygous
#      GBA1LP-like sites are added to both haplotypes.
#   4) Each haplotype is named from the marker set it carries:
#        c.1263del55 + D448H + L483P + A495P + V499=  -> c.1263del55+RecTL
#        D448H + L483P + A495P + V499=                -> RecTL
#        L483P + A495P + V499=                        -> RecNciI
#        A495P + V499=                                -> RecA456P
#        a single marker                              -> that variant
#        no marker                                    -> no_marker
#        any other combination                        -> complex(...)
#      If heterozygous markers fall in different phase blocks, the call is
#      reported as "unphased" with the list of markers.
#
# Outputs:
#   <sample>.phase_links.tsv          link counts between consecutive het sites
#   <sample>.phased_haplotypes.tsv    per-site haplotype A/B alleles and blocks
#   <sample>.recombinant_alleles.tsv  one line: named alleles and evidence
phase_recombinant_markers() {
    local sample_id="$1"
    local bam="$2"
    local counts="$3"
    local out_dir="$4"

    local sites="${out_dir}/${sample_id}.phase_sites.tmp"
    local frag="${out_dir}/${sample_id}.fragment_alleles.tmp"
    local links="${out_dir}/${sample_id}.phase_links.tsv"
    local haps="${out_dir}/${sample_id}.phased_haplotypes.tsv"
    local alleles="${out_dir}/${sample_id}.recombinant_alleles.tsv"

    # SNV PSVs (pos, GBA1 base, GBA1LP base)
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $3, $5, $6 }' "$counts" > "$sites"

    # 1) per-read alleles, merged per fragment: qname, pos, allele (G/P)
    samtools view -q "$MIN_MAPQ" -F 0xF04 "$bam" "$AMP_REGION" \
        | awk -F'\t' -v OFS='\t' -v sites="$sites" -v mbq="$MIN_BASEQ" -v d55s="$DEL55_START" -v d55l="$DEL55_LEN" '
            BEGIN {
                for (i = 33; i < 127; i++) qv[sprintf("%c", i)] = i - 33
                while ((getline line < sites) > 0) { split(line, f, "\t"); ga[f[1]] = f[2]; pa[f[1]] = f[3] }
            }
            {
                pos = $4; c = $6; seq = $10; qual = $11; r = pos; qi = 1; delfound = 0; bigdel = 0
                while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                    l = substr(c, 1, RLENGTH - 1) + 0; op = substr(c, RLENGTH, 1); c = substr(c, RLENGTH + 1)
                    if (op == "M" || op == "=" || op == "X") {
                        for (k = 0; k < l; k++) {
                            p = r + k
                            if (p in ga) {
                                if (qv[substr(qual, qi + k, 1)] >= mbq) {
                                    b = substr(seq, qi + k, 1)
                                    if (b == ga[p]) print $1, p, "G"
                                    else if (b == pa[p]) print $1, p, "P"
                                }
                            }
                        }
                        r += l; qi += l
                    } else if (op == "I" || op == "S") {
                        qi += l
                    } else if (op == "D" || op == "N") {
                        if (op == "D" && l >= d55l - 10 && l <= d55l + 10 && r >= d55s - 15 && r <= d55s + 15) delfound = 1
                        if (op == "D" && l >= 10 && r <= d55s + d55l + 5 && r + l - 1 >= d55s - 5) bigdel = 1
                        r += l
                    }
                }
                if (delfound) print $1, d55s, "P"
                else if (!bigdel && pos <= d55s - 5 && r - 1 >= d55s + d55l + 5) print $1, d55s, "G"
            }' \
        | sort -t$'\t' -k1,1 -k2,2n -k3,3 -u \
        | awk -F'\t' -v OFS='\t' '
            # drop sites where the two mates disagree
            { k = $1 "\t" $2 }
            k == pk { conflict = 1; next }
            { if (pk != "" && !conflict) print prow; prow = $0; pk = k; conflict = 0 }
            END { if (pk != "" && !conflict) print prow }' > "$frag"

    # 2-4) links, haplotypes and allele names
    awk -F'\t' -v OFS='\t' -v sid="$sample_id" -v counts="$counts" -v links="$links" -v haps="$haps" \
        -v minf="$PSV_MIN_FRAC" -v hom="$HOM_FRAC" -v mindp="$PSV_MIN_DP" \
        -v minlinks="$PHASE_MIN_LINKS" -v mincis="$PHASE_MIN_CIS" \
        -v d55s="$DEL55_START" -v markers="$RECOMBINANT_MARKERS" '
        function name_of(set,   n) {
            # set is a "+"-joined list in canonical marker order
            if (set == "") return "no_marker"
            if (set == "c.1263del55+D448H+L483P+A495P+V499=") return "c.1263del55+RecTL"
            if (set == "D448H+L483P+A495P+V499=") return "RecTL"
            if (set == "L483P+A495P+V499=") return "RecNciI"
            if (set == "A495P+V499=") return "RecA456P"
            if (set !~ /\+/) return (set ~ /^c\./ ? set : "p." set)
            return "complex(" set ")"
        }
        BEGIN {
            nm = split(markers, mk, " ")
            for (i = 1; i <= nm; i++) { split(mk[i], mt, ":"); mpos[i] = mt[1]; mlab[i] = mt[2]; ismarker[mt[1]] = mt[2] }
            # site fractions from the PSV counts table
            while ((getline line < counts) > 0) {
                split(line, f, "\t"); if (f[1] == "sample") continue
                if (f[12] == "GBA1_LIKE" || f[12] == "GBA1LP_LIKE") { frac[f[3]] = f[11] + 0; inf[f[3]] = 1 }
            }
        }
        # fragment alleles
        { fa[$1, $2] = $3; if (!($1 in seenq)) { seenq[$1] = 1; q[++nq] = $1 }
          if ($2 == d55s) { if ($3 == "P") d55p++; else d55g++ } }
        END {
            # c.1263del55 site fraction from spanning fragments
            if (d55p + d55g >= mindp) { frac[d55s] = d55p / (d55p + d55g); inf[d55s] = 1 }
            # ordered het and hom-GBA1LP sites
            nh = 0; nsites = 0
            for (s in inf) all[++nsites] = s + 0
            for (i = 2; i <= nsites; i++) { t = all[i]; for (j = i - 1; j >= 1 && all[j] > t; j--) all[j + 1] = all[j]; all[j + 1] = t }
            for (i = 1; i <= nsites; i++) {
                s = all[i]
                if (frac[s] >= minf && frac[s] < hom) het[++nh] = s
                else if (frac[s] >= hom) homp[s] = 1
            }
            # links between consecutive het sites
            print "sample", "site1", "site2", "PP", "GG", "PG", "GP", "cis_fraction", "relation" > links
            for (i = 1; i < nh; i++) {
                a = het[i]; b = het[i + 1]; pp = gg = pg = gp = 0
                for (k = 1; k <= nq; k++) {
                    x = fa[q[k], a]; y = fa[q[k], b]
                    if (x == "" || y == "") continue
                    if (x == "P" && y == "P") pp++; else if (x == "G" && y == "G") gg++
                    else if (x == "P" && y == "G") pg++; else gp++
                }
                tot = pp + gg + pg + gp; cf = (tot > 0) ? (pp + gg) / tot : 0
                rel[i] = (tot < minlinks) ? "NO_LINK" : ((cf >= mincis) ? "CIS" : ((cf <= 1 - mincis) ? "TRANS" : "AMBIGUOUS"))
                nlk[i] = tot
                printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%s\n", sid, a, b, pp, gg, pg, gp, cf, rel[i] > links
            }
            # haplotypes: first het site of each block carries GBA1LP on hap A
            print "sample", "pos", "marker", "GBA1LP_fraction", "zygosity", "block", "hapA", "hapB" > haps
            blk = 0
            for (i = 1; i <= nh; i++) {
                s = het[i]
                if (i == 1 || rel[i - 1] == "NO_LINK" || rel[i - 1] == "AMBIGUOUS") { blk++; ha[s] = "P" }
                else {
                    prev = het[i - 1]
                    ha[s] = (rel[i - 1] == "CIS") ? ha[prev] : (ha[prev] == "P" ? "G" : "P")
                }
                hb[s] = (ha[s] == "P") ? "G" : "P"; bl[s] = blk
            }
            for (i = 1; i <= nsites; i++) {
                s = all[i]
                if (s in bl) printf "%s\t%d\t%s\t%.3f\thet\t%d\t%s\t%s\n", sid, s, (s in ismarker ? ismarker[s] : "."), frac[s], bl[s], (ha[s] == "P" ? "GBA1LP" : "GBA1"), (hb[s] == "P" ? "GBA1LP" : "GBA1") > haps
                else if (s in homp) printf "%s\t%d\t%s\t%.3f\thom\t.\tGBA1LP\tGBA1LP\n", sid, s, (s in ismarker ? ismarker[s] : "."), frac[s] > haps
            }
            # allele naming from markers
            setA = setB = ""; hetblocks = ""; nhb = 0; missing = ""; unph = ""
            for (i = 1; i <= nm; i++) {
                s = mpos[i]; L = mlab[i]
                if (!(s in inf)) { missing = missing (missing == "" ? "" : ",") L; continue }
                if (s in homp) { setA = setA (setA == "" ? "" : "+") L; setB = setB (setB == "" ? "" : "+") L }
                else if (s in bl) {
                    if (!((bl[s]) in hbseen)) { hbseen[bl[s]] = 1; nhb++ }
                    unph = unph (unph == "" ? "" : "+") L
                    if (ha[s] == "P") setA = setA (setA == "" ? "" : "+") L; else setB = setB (setB == "" ? "" : "+") L
                }
            }
            # weakest link along the chain joining the heterozygous markers
            fm = lm = 0; minlk = -1
            for (i = 1; i <= nm; i++) { s = mpos[i]; if ((s in bl) && !(s in homp)) { if (fm == 0 || s < fm) fm = s; if (s > lm) lm = s } }
            for (i = 1; i < nh; i++) if (het[i] >= fm && het[i + 1] <= lm && (minlk < 0 || nlk[i] < minlk)) minlk = nlk[i]
            if (nhb > 1) { call = "unphased(" unph ")"; status = "UNPHASED" }
            else {
                n1 = name_of(setA); n2 = name_of(setB)
                if (n1 == "no_marker" && n2 != "no_marker") { tmp = n1; n1 = n2; n2 = tmp }
                call = n1 "/" n2
                status = (nhb == 1 && unph ~ /\+/) ? "PHASED" : "NOT_REQUIRED"
            }
            if (call == "no_marker/no_marker") call = "none"
            print "sample", "recombinant_alleles", "phasing", "min_link_fragments", "het_markers", "markers_not_assessed"
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", sid, call, status, (status == "PHASED" ? minlk : "."), (unph == "" ? "." : unph), (missing == "" ? "." : missing)
        }' "$frag" > "$alleles"

    rm -f "$sites" "$frag"
    echo "  [phase] $(awk -F'\t' 'NR == 2 { print $2 " (" $3 ")" }' "$alleles") -> $alleles" >&2
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
printf "sample\tmean_amplicon_depth\tdelly_pass\tdelly_pass_types\talignment_events_pass\talignment_event_labels\tpsv_tracts\tpsv_tract_classes\tmarkers_in_tracts\trecombinant_alleles\tphasing\treview\n" > "$SUMMARY"

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
    phase_recombinant_markers "$sid" "$bam" "${sdir}/${sid}.psv_counts.tsv" "$sdir"
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
    ra="${sdir}/${sid}.recombinant_alleles.tsv"
    r_call="$(awk -F'\t' 'NR == 2 { print $2 }' "$ra")"
    r_phase="$(awk -F'\t' 'NR == 2 { print $3 }' "$ra")"
    review="no"
    if [[ "$d_pass" -gt 0 || "$e_pass" -gt 0 || "$t_n" -gt 0 || "${r_call:-none}" != "none" ]]; then review="yes"; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sid" "$mean_depth" "$d_pass" "${d_types:-.}" "$e_pass" "${e_labels:-.}" \
        "$t_n" "${t_cls:-.}" "${t_mk:-.}" "${r_call:-.}" "${r_phase:-.}" "$review" >> "$SUMMARY"
done < "$INPUT_FILE"

if [[ ${#PROFILES[@]} -ge $MIN_COHORT ]]; then
    cohort_depth_flags "$OUT_DIR" "${PROFILES[@]}"
else
    echo "Fewer than ${MIN_COHORT} samples: skipping cohort-relative depth comparison." >&2
fi

echo "Summary -> $SUMMARY" >&2
