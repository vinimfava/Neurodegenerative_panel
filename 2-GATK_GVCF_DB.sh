#!/bin/bash
#SBATCH --account=My_Account
#SBATCH --job-name=DB
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-22%22

module load gatk/4.0.8.1
module load picard/2.18.9
module load annovar/2020Jun08  

# Define working directories (replace with actual paths)
WORKDIR='/path/to/working_directory'
ASSEMBLY='/path/to/assembly_folder'

# Get chromosome number from SLURM array index
chr=$SLURM_ARRAY_TASK_ID

# Step 1: Create a GenomicsDB workspace for this chromosome
gatk --java-options "-Xmx48G -Xms48g" GenomicsDBImport \
  -R $ASSEMBLY/hg38.fa \
  -L chr_${chr} \
  --batch-size 50 \
  --reader-threads 8 \
  --verbosity INFO \
  --sample-name-map $WORKDIR/Joined_genotyping/Samples_GVCFs.list \
  --genomicsdb-workspace-path $WORKDIR/Joined_genotyping/chr_${chr}

# Step 2: Genotype the variants using the GenomicsDB
gatk GenotypeGVCFs \
  -R $ASSEMBLY/hg38.fa \
  -V gendb://$WORKDIR/Joined_genotyping/chr_${chr} \
  -D $ASSEMBLY/Annotation/dbsnp_151.hg38.vcf.gz \
  -G StandardAnnotation -new-qual \
  -O $WORKDIR/Joined_genotyping/chr_${chr}_samples.vcf

# Step 3: Merge per chr VCFs
java -Xmx16G -jar $EBROOTPICARD/picard.jar MergeVcfs \
I=$WORKING/Joined_genotyping/chr1_samples.vcf \
I=$WORKING/Joined_genotyping/chr2_samples.vcf \
I=$WORKING/Joined_genotyping/chr3_samples.vcf \
I=$WORKING/Joined_genotyping/chr4_samples.vcf \
I=$WORKING/Joined_genotyping/chr5_samples.vcf \
I=$WORKING/Joined_genotyping/chr6_samples.vcf \
I=$WORKING/Joined_genotyping/chr7_samples.vcf \
I=$WORKING/Joined_genotyping/chr8_samples.vcf \
I=$WORKING/Joined_genotyping/chr9_samples.vcf \
I=$WORKING/Joined_genotyping/chr10_samples.vcf \
I=$WORKING/Joined_genotyping/chr11_samples.vcf \
I=$WORKING/Joined_genotyping/chr12_samples.vcf \
I=$WORKING/Joined_genotyping/chr13_samples.vcf \
I=$WORKING/Joined_genotyping/chr14_samples.vcf \
I=$WORKING/Joined_genotyping/chr15_samples.vcf \
I=$WORKING/Joined_genotyping/chr16_samples.vcf \
I=$WORKING/Joined_genotyping/chr17_samples.vcf \
I=$WORKING/Joined_genotyping/chr18_samples.vcf \
I=$WORKING/Joined_genotyping/chr19_samples.vcf \
I=$WORKING/Joined_genotyping/chr20_samples.vcf \
I=$WORKING/Joined_genotyping/chr21_samples.vcf \
I=$WORKING/Joined_genotyping/chr22_samples.vcf \
O=$WORKING/Joined_genotyping/GATK_joined_samples.vcf
bgzip $WORKING/Joined_genotyping/GATK_joined_samples.vcf
tabix $WORKING/Joined_genotyping/GATK_joined_samples.vcf.gz

# Step 4: Filtering
vcftools --gzvcf $WORKING/Joined_genotyping/GATK_joined_samples.vcf.gz \
--minQ 30 \
--minDP 10 \
--min-meanDP 10 \
--max-missing 0.2 \
--recode \
--recode-INFO-all \
--stdout | bgzip > $WORKING/Joined_genotyping/GATK_joined_samples_filter_1.vcf.gz

gatk VariantFiltration \
-R $ASSEMBLY/hg38.fa \
-V $WORKING/Joined_genotyping/GATK_joined_samples_filter_1.vcf.gz\
-O $WORKING/Joined_genotyping/GATK_joined_samples_final.vcf \
--filter-name "QD2.5" --filter-expression "QD < 2.5" \
--filter-name "SOR2" --filter-expression "SOR > 2.0" \
--filter-name "FS20" --filter-expression "FS > 20.0" \
--filter-name "MQ59" --filter-expression "MQ < 59.0" \
--filter-name "MQRankSum5" --filter-expression "MQRankSum < -5.0" \
--filter-name "ReadPosRankSum4" --filter-expression "ReadPosRankSum < -4.0"
bgzip $WORKING/Joined_genotyping/GATK_joined_samples_final.vcf
tabix $WORKING/Joined_genotyping/GATK_joined_samples_final.vcf.gz

# Step 5: Functional Annotation with ANNOVAR
$ANNOVAR/table_annovar.pl \
  -vcfinput $WORKING/Joined_genotyping/GATK_joined_samples_final.vcf.gz \
  -tempdir $WORKDIR/temp/ \
  $ASSEMBLY/Annotation/ \
  --buildver hg38 \
  -out $WORKDIR/annotation/Neuro_panel_annotated \
  -protocol refGene,avsnp150,dbnsfp35c,clinvar_20190305,revel,gnomad211_exome,gnomad30_genome \
  -operation g,f,f,f,f,f,f \
  -nastring .
bgzip $WORKDIR/annotation/GATK_joined_samples_final.hg38_multianno.vcf
tabix $WORKDIR/annotation/GATK_joined_samples_final.hg38_multianno.vcf.gz

# Step 6: Select Protein-Altering Variants
# select protein altering variants after the annotation
vcftools --gzvcf $WORKDIR/annotation/GATK_joined_samples_final.hg38_multianno.vcf.gz \
  --bed $WORKDIR/Protein_altering.bed \
  --recode --recode-INFO-all --stdout | bgzip > $WORKDIR/annotation/Neuro_panel_annotated_final_GATK.vcf.gz
tabix -f $WORKDIR/annotation/Neuro_panel_annotated_final_GATK.vcf.gz
