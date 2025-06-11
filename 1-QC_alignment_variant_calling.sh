#!/bin/bash
#SBATCH --account=My_Account
#SBATCH --job-name=NeuroPanel
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-940%940

# Load required modules
module load fastqc/0.11.5
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.18.9
module load gatk/4.0.8.1
module load singularity

# Define directory paths (update with your actual structure)
WORKDIR='/path/to/working_directory'       
SCRATCH='/path/to/scratch_directory'        
ASSEMBLY='/path/to/assembly_files'          
BIN_VERSION=1.0.0

# Get sample name from sample list (1-based index from SLURM array)
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $WORKDIR/sample.txt)

# Step 1: Adapter Trimming
/path/to/TrimGalore-0.4.5/trim_galore \
  --quality 20 \
  --nextera --gzip --paired \
  -o $SCRATCH \
  $WORKDIR/${sample}_R1.fastq.gz $WORKDIR/${sample}_R2.fastq.gz

# Step 2: Alignment with BWA
bwa mem -M -P \
  -t 4 \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:Neuro_Panel" \
  $ASSEMBLY/hg38.fa \
  $SCRATCH/${sample}_R1_val_1.fq.gz \
  $SCRATCH/${sample}_R2_val_2.fq.gz > $SCRATCH/${sample}.sam

# Convert and sort BAM
samtools view -@ 4 -b $SCRATCH/${sample}.sam | \
samtools sort -@ 4 -o $SCRATCH/${sample}_raw_sorted.bam

# Step 3: Mark Duplicates
java -Xmx16G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I=$SCRATCH/${sample}_raw_sorted.bam \
  O=$SCRATCH/${sample}_sorted_dups.bam \
  M=$WORKDIR/metrix/${sample}.marked_dup_metrics.txt \
  CREATE_INDEX=true

# Step 4: Base Recalibration
gatk BaseRecalibrator \
  -R $ASSEMBLY/hg38.fa \
  -I $SCRATCH/${sample}_sorted_dups.bam \
  --known-sites $ASSEMBLY/Annotation/dbsnp_151.hg38.vcf.gz \
  -O $WORKDIR/metrix/${sample}_recal_data.table

gatk ApplyBQSR \
  -R $ASSEMBLY/hg38.fa \
  -I $SCRATCH/${sample}_sorted_dups.bam \
  --bqsr-recal-file $WORKDIR/metrix/${sample}_recal_data.table \
  -O $WORKDIR/Alignment/${sample}_bqsr_sorted_marked_dups.bam

# Step 5: Variant Calling
gatk HaplotypeCaller \
  -R $ASSEMBLY/hg38.fa \
  -I $WORKDIR/Alignment/${sample}_bqsr_sorted_marked_dups.bam \
  --emit-ref-confidence GVCF \
  -D $ASSEMBLY/Annotation/dbsnp_151.hg38.vcf.gz \
  -L $WORKDIR/Targeted_regions.bed \
  -RF NotDuplicateReadFilter \
  -RF NotOpticalDuplicateReadFilter \
  -RF MappedReadFilter \
  -O $WORKDIR/GVCFs/${sample}_raw.g.vcf
  
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref=$ASSEMBLY/hg38.fa \
  --reads=$WORKDIR/Alignment/${sample}_bqsr_sorted_marked_dups.bam \
  --regions=$WORKDIR/Neuro_panel_gene_exons.bed \
  --output_vcf=$WORKDIR/output/${sample}_output.vcf.gz \
  --output_gvcf=$WORKDIR/output/${sample}_output.g.vcf.gz \
  --sample_name=${sample} \
  --intermediate_results_dir=$WORKDIR/temp/${sample} \
  --num_shards=4
