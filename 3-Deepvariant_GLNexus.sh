#!/bin/bash
#SBATCH --account=My_Account
#SBATCH --job-name=NeuroPanel
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

WORKDIR="/path/to/working_directory"
ASSEMBLY="/path/to/assembly_folder"
ANNOVAR="/path/to/annovar"

module load bcftools/1.17
module load vcftools/0.1.16
module load gatk/4.2.6.1
module load annovar/2020Jun08  
module load glnexus/1.4.1

# Step 1: Joint Genotyping with GLnexus
glnexus_cli --config DeepVariantWES \
  --bed $WORKDIR/Neuro_panel_gene_exons.bed \
  $WORKDIR/output/*.g.vcf.gz > $WORKDIR/Neuro_panel_joint_raw.bcf

# Convert BCF to compressed VCF
bcftools view $WORKDIR/Neuro_panel_joint_raw.bcf | bgzip -c > $WORKDIR/Neuro_panel_joint_raw.vcf.gz
tabix $WORKDIR/Neuro_panel_joint_raw.vcf.gz

# Step 2: Filtering with VCFtools
vcftools --gzvcf $WORKDIR/Neuro_panel_joint_raw.vcf.gz \
  --minQ 30 \
  --minDP 10 \
  --minGQ 10 \
  --min-meanDP 10 \
  --max-missing 0.2 \
  --recode --recode-INFO-all --stdout | bgzip > $WORKDIR/Neuro_panel_filtered.vcf.gz
tabix $WORKDIR/Neuro_panel_filtered.vcf.gz

# Step 3: GATK Variant Annotation
gatk VariantAnnotator \
  --reference $ASSEMBLY/hg38.fa \
  --variant $WORKDIR/Neuro_panel_filtered.vcf.gz \
  --annotation QualByDepth \
  --annotation StrandOddsRatio \
  --annotation MappingQuality \
  --annotation MappingQualityRankSumTest \
  --annotation FisherStrand \
  --annotation ReadPosRankSumTest \
  --annotation Coverage \
  --annotation GenotypeSummaries \
  --dbsnp $ASSEMBLY/Annotation/dbsnp_151.hg38.vcf.gz \
  --output $WORKDIR/Neuro_panel_annotated.vcf.gz

# Step 4: Functional Annotation with ANNOVAR
$ANNOVAR/table_annovar.pl \
  -vcfinput $WORKDIR/Neuro_panel_annotated.vcf.gz \
  -tempdir $WORKDIR/temp/ \
  $ASSEMBLY/Annotation/ \
  --buildver hg38 \
  -out $WORKDIR/annotation/Neuro_panel_annotated \
  -protocol refGene,avsnp150,dbnsfp35c,clinvar_20190305,revel,gnomad211_exome,gnomad30_genome \
  -operation g,f,f,f,f,f,f \
  -nastring .
bgzip $WORKDIR/annotation/Neuro_panel_annotated.hg38_multianno.vcf
tabix $WORKDIR/annotation/Neuro_panel_annotated.hg38_multianno.vcf.gz

# Step 5: Select Protein-Altering Variants
# select protein altering variants after the annotation
vcftools --gzvcf $WORKDIR/annotation/Neuro_panel_annotated.hg38_multianno.vcf.gz \
  --bed $WORKDIR/Protein_altering.bed \
  --recode --recode-INFO-all --stdout | bgzip > $WORKDIR/annotation/Neuro_panel_annotated_final_DeepVariant.vcf.gz
tabix -f $WORKDIR/annotation/Neuro_panel_annotated_final_DeepVariant.vcf.gz
