#!/bin/bash
#SBATCH --account=My_Account
#SBATCH --job-name=NeuroPanel
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load bcftools/1.17

# Step 1: Extract overlapping variants between GATK and DeepVariant
bcftools isec \
  -p $WORKDIR/annotation/isec_GATK_DeepVariant \
  $WORKDIR/annotation/Neuro_panel_annotated_final_GATK.vcf.gz \
  $WORKDIR/annotation/Neuro_panel_annotated_final_DeepVariant.vcf.gz

# Shared variants will be in 0003.vcf (both files)
bgzip -f $WORKDIR/annotation/isec_GATK_DeepVariant/0003.vcf
tabix -f $WORKDIR/annotation/isec_GATK_DeepVariant/0003.vcf.gz

# Step 2: Concatenate shared variants with CNVs from XHMM
bcftools concat -a -O z \
  $WORKDIR/annotation/isec_GATK_DeepVariant/0003.vcf.gz \
  $WORKDIR/annotation/Neuro_panel_XHMM_final.vcf.gz \
  -o $WORKDIR/annotation/Neuro_panel_final.vcf.gz
  
tabix -f $WORKDIR/annotation/Neuro_panel_final.vcf.gz