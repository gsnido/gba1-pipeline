# GBA1 variant calling pipeline

Bash pipeline for SNV and small-indel calling from paired-end sequencing of
PCR-amplified GBA1. It implements the tailored bioinformatic analysis of the
LONG-NEXT method (Cuconato et al., 2025). Structural-variant and
recombinant-allele screening is provided separately by `LONG-NEXT-SV.sh`
(see "Structural variants and recombinant alleles" below).


## Requirements

- Bash and standard Unix utilities, including awk and find.
- BWA, SAMtools, GATK4 (with a compatible Java runtime), and BCFtools.
- Delly (>= 1.2; `delly sr`) for `LONG-NEXT-SV.sh`.
- An hg38 reference FASTA with GBA1LP already masked, plus the BWA indexes, FASTA index, and GATK sequence dictionary. See section "Reference preparation" below.
- A GBA1 target-region BED file matching the reference assembly and contig names.

Additional software for building the masked reference:
- BEDTools
- Reference in fasta format

The reference can be obtained using Google cloud CLI:

```bash
gcloud storage cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy/v0/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta references/
```



## Inputs and usage

```bash
bash LONG-NEXT.sh \
  sample_ids.txt \
  /path/to/fastqs \
  /path/to/hg38.GBA1LP_masked.fa \
  /path/to/GBA1.hg38.bed \
  /path/to/results \
  8
```

Arguments, in order:

1. Sample ID list: one ID per line.
2. FASTQ directory.
3. Pre-masked reference FASTA.
4. GBA1 target BED.
5. Output directory.
6. Optional thread count (default: 4; passed to BWA and SAMtools sorting).

Sample IDs have dashes and whitespace replaced with underscores. The script
searches only the top level of the FASTQ directory and requires exactly two
files matching `<sanitized_id>_*`, identifiable by `_R1_`/`_R2_` or
`_R1.`/`_R2.` in their names. Samples with missing or ambiguous pairs are
skipped; __multi-lane inputs must be consolidated beforehand__.

Example command:
```bash
bash LONG-NEXT.sh examples/ids.txt examples/fastqfiles/ \
  references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.GBA1LP_masked.fa \
  bedfiles/GBA1.bed \
  examples/results \
  4
```

## Workflow

1. Align paired reads to the pre-masked reference with BWA-MEM; sort and index BAM files.
2. Calculate mean GBA1 depth, with an intended warning threshold of 500x.
3. Call variants using GATK HaplotypeCaller in gVCF mode over the target region with 100 bp interval padding.
4. Genotype each sample independently using GenotypeGVCFs.
5. Apply separate SNP/indel hard filters and merge the resulting records.
6. Add depth and allele-balance FILTER labels with BCFtools.

Additional FILTER labels:

| Label | Threshold in the script |
|---|---|
| `LOW_DEPTH` | FORMAT/DP < 30 |
| `LOW_VAF` | First alternate-allele fraction < 0.15 |
| `HET_VAF_OUTLIER` | Heterozygous call with first alternate-allele fraction outside 0.30–0.70 |

Allele fraction is calculated as AD[1] / (AD[0] + AD[1]); this is not a general
multiallelic calculation. Filters label records rather than removing them. The
depth and allele-balance flags are QC indicators, not confirmation of allele
dropout or recombination.

## Outputs

The output directory contains sorted/indexed BAMs, gVCFs, genotyped VCFs,
intermediate hard-filtered VCFs, and `<sample>.final_filtered.vcf.gz` with its
tabix index. The depth log is written to `log_GBA1_depth.txt` in the output
directory. Progress messages are written to stderr.

## Run an example

```bash
bash LONG-NEXT.sh examples/ids.txt examples/fastqfiles/ \
  /references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.GBA1LP_masked.fa \
  GBA1.bed \
  examples/results \
  4
```

## Structural variants and recombinant alleles

`LONG-NEXT-SV.sh` runs after `LONG-NEXT.sh` and reuses its
`<sample>.sorted.bam` files (no realignment):

```bash
bash LONG-NEXT-SV.sh \
  examples/ids.txt \
  examples/results \
  references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.GBA1LP_masked.fa \
  bedfiles/GBA1_GBA1LP_PSVs.hg38.tsv \
  examples/results \
  4 \
  [reliable_sites]
```

Arguments: sample ID list, LONG-NEXT.sh output directory (BAMs), masked
reference, PSV table, output directory, optional threads and an optional
`reliable_sites` file (see below).

All analyses are restricted to the LONG-NEXT amplicon
(chr1:155,234,915-155,241,415, hg38; defined by the published primers):

| Evidence | Method | Output (`<sample>_SV/`) |
|---|---|---|
| SV calls (DEL/DUP/INV/INS/BND) | Delly short-read mode; extra FILTER labels `AMPLICON_EDGE` (breakpoint < 300 bp from an amplicon end) and `LOW_SUPPORT` (supporting fraction < 0.05) | `<sample>.delly.sv.vcf.gz` |
| Mid-sized indels and junctions | CIGAR indels >= 30 bp and split-read (SA tag) junctions with >= 10 reads and >= 5% of local depth; `c.1263_1317del` (Rec(delta55)) is labelled | `<sample>.alignment_events.tsv` |
| Gene conversion / recombinant tracts | Allele counts at 134 GBA1/GBA1LP paralogous sequence variants (PSVs); runs of >= 2 consecutive PSVs with GBA1LP-allele fraction >= 0.15 are reported as tracts (`HET_LIKE`, `HOM_LIKE`, `LOW_FRACTION`, `UNBALANCED`) and whether they reach the 3' end (fusion/RecNciI-type) | `<sample>.psv_counts.tsv`, `<sample>.psv_tracts.tsv` |
| Phasing and named recombinant alleles | Read-backed phasing of heterozygous PSVs and the c.1263del55 site (mates merged into fragments); alleles named per haplotype from the marker set (table below) | `<sample>.phase_links.tsv`, `<sample>.phased_haplotypes.tsv`, `<sample>.recombinant_alleles.tsv` |
| Depth profile | 100 bp bins normalised to the sample median; with >= 3 samples, ratio to the cohort median and runs of >= 3 bins < 0.6 or > 1.5 | `<sample>.depth_profile.tsv`, `SV_cohort_depth_*.tsv` |

`SV_summary.tsv` has one line per sample with the counts of PASS
calls/events/tracts, the markers found inside tracts, the named
`recombinant_alleles` (e.g. `RecNciI/no_marker`), the `phasing` status and a
`review` column.

### Phasing and named alleles

Each primary read records the GBA1 or GBA1LP allele at every PSV
(base quality >= `MIN_BASEQ`). The c.1263del55 site counts as one extra
biallelic site: GBA1LP-like if the read carries a ~55 bp deletion there and
GBA1-like if it spans the region without one. Both mates are merged into one
fragment. Consecutive heterozygous sites are linked when >= `PHASE_MIN_LINKS`
(10) fragments cover both. A link is `CIS` if >= `PHASE_MIN_CIS` (0.80) of
fragments are GBA1LP-GBA1LP or GBA1-GBA1, `TRANS` if <= 0.20 and `AMBIGUOUS`
otherwise. `AMBIGUOUS` or missing links split phase blocks. Homozygous
GBA1LP-like sites are added to both haplotypes.

Haplotypes are named from the markers they carry, following the allele
definitions used by Gauchian ([Toffoli et al., Commun Biol 2022](https://www.nature.com/articles/s42003-022-03610-7)):

| Markers on the haplotype (HGVS, legacy) | Named allele |
|---|---|
| c.1263del55 + D448H + L483P + A495P + V499= | `c.1263del55+RecTL` |
| D448H (D409H) + L483P (L444P) + A495P (A456P) + V499= (V460V) | `RecTL` |
| L483P + A495P + V499= | `RecNciI` |
| A495P + V499= | `RecA456P` |
| a single marker | that variant (e.g. `p.L483P`, `c.1263del55`) |
| no marker | `no_marker` |
| any other combination | `complex(...)` |

Marker positions (hg38): c.1263del55 starts at chr1:155,235,750, D448H
155,235,727, L483P 155,235,252, A495P 155,235,217, V499= 155,235,203.

The `phasing` column is `PHASED` when two or more heterozygous markers fall in
one phase block, `NOT_REQUIRED` when the call does not depend on phasing (no
marker, a single heterozygous marker or only homozygous markers) and
`UNPHASED` when heterozygous markers fall in different blocks. In that case
the call is `unphased(<markers>)`. `recombinant_alleles.tsv` also reports the
weakest link along the phased chain (`min_link_fragments`) and markers without
coverage (`markers_not_assessed`). With 2x150 bp reads the c.1263del55 site
(about 500 bp from L483P) is phased through fragments spanning the insert, so
fragment counts across that link are lower.

Example calls on simulated data and the positive control:

| Sample | `recombinant_alleles` | `phasing` |
|---|---|---|
| Positive control | `RecNciI/no_marker` | `PHASED` |
| Heterozygous RecTL | `RecTL/no_marker` | `PHASED` |
| L483P in trans with RecA456P | `RecA456P/p.L483P` | `PHASED` (TRANS link) |
| Homozygous RecNciI | `RecNciI/RecNciI` | `NOT_REQUIRED` |

### Reliable sites (optional)

The optional 7th argument is a whitespace-separated file whose first two
columns are chromosome and hg38 position, for example Gauchian's
`gauchian/data/GBA_SNP_38.txt`. Matching PSVs are flagged in
`psv_counts.tsv` (`reliable_site`) and counted in `psv_tracts.tsv`
(`n_reliable_psvs`); calls are not changed. The Gauchian file is not bundled
because Gauchian is released under the PolyForm Strict 1.0.0 licence. Only
34 of its 83 hg38 sites lie inside the amplicon and 31 of those are already
in the bundled table. The other 3 are the c.1263del55 boundaries (handled by
the deletion site) and an indel PSV.

The PSV table (`bedfiles/GBA1_GBA1LP_PSVs.hg38.tsv`) was derived by tiling the
hg38 GBA1LP interval (`bedfiles/GBA1LP.bed`) onto the masked reference; the
derivation is described in the file header. Because the hg38 GBA1LP sequence
matches GBA1 at c.1483 and c.1497, the RecNciI markers c.1483G>C (A456P) and
c.1497G>C (V460V) were added explicitly as `known_marker`.

Limitations: phasing relies on fragments covering consecutive heterozygous
sites, so long PSV-free stretches can split blocks. Only events inside the amplicon and compatible with the primer
design are visible. Events removing a primer site (e.g. some reciprocal
recombinants or large deletions) are not amplified and appear as allele
dropout. PCR favours shorter alleles, so deletion-allele fractions may be
inflated. Copy number cannot be inferred from amplicons; confirm any flag
with MLPA or long-read sequencing.

## Docker image

It is possible to build a Docker image:

```bash
# Compile image (requires access to the internet)
docker build -t gba1-pipeline:dev .

# Test with commercial positive control sample
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD/examples/ids.txt,dst=/input/sample_ids.txt,readonly" \
  --mount "type=bind,src=$PWD/examples/fastqfiles,dst=/input/fastqfiles,readonly" \
  --mount "type=bind,src=$PWD/references,dst=/references,readonly" \
  --mount "type=bind,src=$PWD/examples/results,dst=/output" \
  --workdir /output \
  gba1-pipeline:dev \
  LONG-NEXT \
    /input/sample_ids.txt \
    /input/fastqfiles \
    /references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.GBA1LP_masked.fa \
    /opt/gba1/GBA1.bed \
    /output \
    4
```

To save the image into a tar:

```bash
docker save --output="image_name.tar" id_image
docker load --input image_name.tar
```


## Reference preparation

1. Download the hg38 reference FASTA from Google cloud storage:

```bash
# Get hg38 reference
gcloud storage cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy/v0/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta references/
```

2. Mask GBA1LP

```bash
# Mask pseudogene
bash PREPARE_REF.sh references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta bedfiles/GBA1LP.bed
```


## TO DO

1. Consider modifying the LONG-NEXT.sh script to analyse a single file to relax
   the requirement of a text file with IDs. Really annoying for Docker. Also,
   there is a script (CHECK_IDS.sh) that should probably be deleted since it's
   not used. Yet, also at this stage we should test wheter there are several
   pairs of files compatible with the ID and decide whether to consolidate
   fastq files from different lanes or exit with error.
2. Test the JOINT_GENOTYPING script.
3. Validate `LONG-NEXT-SV.sh` on samples with known recombinant alleles and
   deletions (thresholds were tuned on simulated data and one positive control).


