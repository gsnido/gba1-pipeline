# GBA1 variant calling pipeline

Bash pipeline for SNV and small-indel calling from paired-end sequencing of
PCR-amplified GBA1. It implements the tailored bioinformatic analysis of the
LONG-NEXT method (Cuconato et al., 2025). It does not perform
structural-variant calling.


## Requirements

- Bash and standard Unix utilities, including awk and find.
- BWA, SAMtools, GATK4 (with a compatible Java runtime), and BCFtools.
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
3. Create the structural variant caller.


