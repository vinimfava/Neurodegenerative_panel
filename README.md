# Neurodegenerative panel pipeline for leprosy Type-1 reaction

This repository contains a modular pipeline for processing targeted sequencing data using GATK, DeepVariant, XHMM, and Regenie.

The pipeline is structured in six stages:

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1️⃣ | [`1-QC_alignment_variant_calling.sh`](./1-QC_alignment_variant_calling.sh) | Raw FASTQ processing, BWA alignment, BQSR, GVCF calling (GATK + DeepVariant) |
| 2️⃣ | [`2-GATK_GVCF_DB.sh`](./2-GATK_GVCF_DB.sh) | GATK joint genotyping, filtering, and annotation |
| 3️⃣ | [`3-Deepvariant_GLNexus.sh`](./3-Deepvariant_GLNexus.sh) | DeepVariant + GLnexus joint calling and annotation |
| 4️⃣ | [`4-XHMM.sh`](./4-XHMM.sh) | CNV calling using XHMM from BAM alignments |
| 5️⃣ | [`5-Combine_VCFs.sh`](./5-Combine_VCFs.sh) | Merging overlapping SNVs from GATK & DeepVariant with CNVs |
| 6️⃣ | [`6-Analysis.sh`](./6-Analysis.sh) | Association testing using Regenie (single-variant & gene-wise) |

---

## Input 
Each script contains `WORKDIR` and `ASSEMBLY` placeholders and account was omited as the analysis was carried out in a private file system:

- Raw FASTQ files are not privided given conset to share patient's raw data.
- Reference genome: hg38

---

## ⚙️ Pipeline Stages

### 1. Quality Control & Alignment
`1-QC_alignment_variant_calling.sh`
- Trimming with Trim Galore
- Alignment using BWA
- BQSR and duplicate marking (Picard + GATK)
- Variant calling using GATK HaplotypeCaller and DeepVariant

### 2. GATK Genotyping & Annotation
`2-GATK_GVCF_DB.sh`
- Chromosome-level joint genotyping via GenomicsDBImport
- Variant filtration with VCFtools and GATK
- Functional annotation using ANNOVAR

### 3. DeepVariant Joint Calling
`3-Deepvariant_GLNexus.sh`
- Merges `.g.vcf.gz` using GLnexus
- Filters and annotates with GATK and ANNOVAR
- Final variant set includes protein-altering variants

### 4. CNV Calling with XHMM
`4-XHMM.sh`
- Depth normalization, PCA, Z-score matrix
- CNV discovery and genotyping
- Output is a compressed, indexed VCF of CNVs

### 5. Merge VCFs
`5-Combine_VCFs.sh`
- Uses `bcftools isec` to find shared variants from GATK & DeepVariant
- Concatenates them with CNVs from XHMM

### 6. Regenie-Based Analysis
`6-Analysis.sh`
- Converts VCF to PLINK2 `.pgen`
- Runs Regenie Step 1 and Step 2 for both:
  - Single-variant logistic regression
  - Gene-wise burden testing (SKAT-O / ACAT)

### 7. R script for meta-analysis
`7-Meta-analysis.R`
- Calulcate odds ratio (OR) and meta-analysis using betsa and standard errors from step 6

---

## 📌 Notes

- Most scripts are written for HPC use with SLURM.
- Some scripts use Dockerized tools (e.g. Regenie), others require module loading.

---

## 📬 Results

The pre-print of this study is availaible at:
https://www.medrxiv.org/content/10.1101/2025.06.05.25329059v1.article-metrics

---

