#!/bin/bash
#SBATCH --account=My_account
#SBATCH --job-name=CNV_XHMM
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

module load gatk/3.8
module load xhmm
module load plinkseq
module load tabix
module load bcftools

WORKDIR="/path/to/working_directory"
ASSEMBLY="/path/to/assembly_files/hg38.fa"    # Make sure to include full path to .fa
TARGETS="$WORKDIR/Targeted_probes_depth_hg38_prepared.list"
PARAMS="$WORKDIR/master/params.txt"
SEQDB="$WORKDIR/plinkseq/seqdb"
ANNOTATION_DIR="$WORKDIR/annotation"

# Step 1: Run DepthOfCoverage
for i in {1..22}; do
  java -Xmx48G -Xms48g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -I "$WORKDIR/batch${i}.bam.list" \
    -L "$TARGETS" \
    -R "$ASSEMBLY" \
    --downsampling_type BY_SAMPLE \
    --downsample_to_coverage 500 \
    --logging_level INFO \
    --minBaseQuality 0 \
    --minMappingQuality 20 \
    --start 1 --stop 101 --nBins 100 \
    --omitDepthOutputAtEachBase \
    --omitLocusTable \
    --includeRefNSites \
    --countType COUNT_FRAGMENTS \
    -o "$WORKDIR/batch${i}.DATA"
done

# Step 2: Merge all interval summaries
xhmm --mergeGATKdepths -o "$WORKDIR/Neuro.DATA.RD.txt" \
  $(for i in {1..22}; do echo "--GATKdepths $WORKDIR/batch${i}.DATA.sample_interval_summary"; done)

# Step 3: Compute GC content
java -Xmx48G -Xms48g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
  -T GCContentByInterval \
  -L "$TARGETS" \
  -R "$ASSEMBLY" \
  -o "$WORKDIR/Neuro.DATA.locus_GC.txt"

awk '$2 < 0.2 || $2 > 0.8 {print $1}' "$WORKDIR/Neuro.DATA.locus_GC.txt" > "$WORKDIR/Neuro.DATA.locus_GC_extreme.txt"

# Step 4: Compute target complexity with PLINK/Seq
plinkseq/interval_list_to_pseq_reg.sh "$TARGETS" > "$WORKDIR/Neuro.targets.reg"

plinkseq/pseq . loc-load \
  --locdb "$WORKDIR/Neuro.targets.LOCDB" \
  --file "$WORKDIR/Neuro.targets.reg" \
  --group targets \
  --out "$WORKDIR/Neuro.targets.LOCDB.loc-load"

plinkseq/pseq . loc-stats \
  --locdb "$WORKDIR/Neuro.targets.LOCDB" \
  --group targets \
  --seqdb "$SEQDB" | \
awk 'NR > 1 {print $0}' | \
sort -k1,1g | \
awk '{print $10}' | \
paste "$TARGETS" - | \
awk '{print $1"\t"$2}' > "$WORKDIR/Neuro.DATA.locus_complexity.txt"

awk '$2 > 0.25 {print $1}' "$WORKDIR/Neuro.DATA.locus_complexity.txt" > "$WORKDIR/Neuro.DATA.locus_complexity_low.txt"

# Step 5: Filter and mean-center data
xhmm --matrix -r "$WORKDIR/Neuro.DATA.RD.txt" --centerData --centerType target \
  -o "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt" \
  --outputExcludedTargets "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt.filtered_targets.txt" \
  --outputExcludedSamples "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt.filtered_samples.txt" \
  --excludeTargets "$WORKDIR/Neuro.DATA.locus_GC_extreme.txt" \
  --excludeTargets "$WORKDIR/Neuro.DATA.locus_complexity_low.txt" \
  --minMeanTargetRD 10 --maxMeanTargetRD 40 \
  --minMeanSampleRD 10 --maxMeanSampleRD 60 \
  --maxSdSampleRD 30

# Step 6: PCA
xhmm --PCA -r "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt" \
  --PCAfiles "$WORKDIR/Neuro.DATA.RD_PCA"

# Step 7: Normalize using PCA
xhmm --normalize -r "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt" \
  --PCAfiles "$WORKDIR/Neuro.DATA.RD_PCA" \
  --normalizeOutput "$WORKDIR/Neuro.DATA.PCA_normalized.txt" \
  --PCnormalizeMethod PVE_mean \
  --PVE_mean_factor 0.7

# Step 8: Z-score by sample
xhmm --matrix -r "$WORKDIR/Neuro.DATA.PCA_normalized.txt" \
  --centerData --centerType sample --zScoreData \
  -o "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt" \
  --outputExcludedTargets "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt" \
  --outputExcludedSamples "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt" \
  --maxSdTargetRD 30

# Step 9: Filter original RD matrix
xhmm --matrix -r "$WORKDIR/Neuro.DATA.RD.txt" \
  --excludeTargets "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt.filtered_targets.txt" \
  --excludeTargets "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt" \
  --excludeSamples "$WORKDIR/Neuro.DATA.filtered_centered.RD.txt.filtered_samples.txt" \
  --excludeSamples "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt" \
  -o "$WORKDIR/Neuro.DATA.same_filtered.RD.txt"

# Step 10: CNV Discovery
xhmm --discover -p "$PARAMS" \
  -r "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt" \
  -R "$WORKDIR/Neuro.DATA.same_filtered.RD.txt" \
  -c "$WORKDIR/Neuro.DATA.xcnv" \
  -a "$WORKDIR/Neuro.DATA.aux_xcnv" \
  -s "$WORKDIR/Neuro.DATA"

# Step 11: CNV Genotyping
xhmm --genotype -p "$PARAMS" \
  -r "$WORKDIR/Neuro.DATA.PCA_normalized.filtered.sample_zscores.RD.txt" \
  -R "$WORKDIR/Neuro.DATA.same_filtered.RD.txt" \
  -g "$WORKDIR/Neuro.DATA.xcnv" \
  -F "$ASSEMBLY" \
  -v "$WORKDIR/annotation/Neuro_panel_XHMM_final.vcf"

# Compress and index final VCF
bgzip -f "$WORKDIR/annotation/Neuro_panel_XHMM_final.vcf"
tabix -f "$WORKDIR/annotation/Neuro_panel_XHMM_final.vcf.gz"

