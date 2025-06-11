#!/bin/bash
#SBATCH --account=My_Account
#SBATCH --job-name=NeuroPanel
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load plink/2.00
WORKDIR="/path/to/working/dir"

# Step 1: Create PLINK2 binaries
plink2 \
  --vcf $WORKDIR/Neuro_panel_samples.vcf.gz \
  --max-alleles 2 \
  --make-pgen \
  --mac 1 \
  --sort-vars \
  --vcf-half-call missing \
  --pheno T1R_Neuro_samples.pheno \
  --pheno-name T1R \
  --1 \
  --update-sex T1R_Neuro_samples.pheno \
  --out Neuro_panel_samples

# Step 2: NULL model for T1R (single variant)
docker run --rm -v "$WORKDIR/Data:/data/input" -v "$WORKDIR/Results:/data/results" \
  ghcr.io/rgcgithub/regenie/regenie:v3.6.gz /usr/local/bin/regenie \
  --pgen /data/input/Neuro_panel_samples \
  --phenoFile /data/input/T1R_Neuro_samples.pheno \
  --covarFile /data/input/T1R_Neuro_samples.cov \
  --covarColList SEX \
  --step 1 \
  --phenoCol T1R \
  --bt \
  --bsize 5000 \
  --strict \
  --threads 8 \
  --out /data/results/STEP_1_T1R

# Step 3: Single variant analysis
docker run --rm -v "$WORKDIR/Data:/data/input" -v "$WORKDIR/Results:/data/results" \
  ghcr.io/rgcgithub/regenie/regenie:v3.6.gz /usr/local/bin/regenie \
  --pgen /data/input/Neuro_panel_samples \
  --phenoFile /data/input/T1R_Neuro_samples.pheno \
  --pred /data/results/STEP_1_T1R_pred.list \
  --phenoCol T1R \
  --step 2 \
  --bt \
  --af-cc \
  --bsize 5000 \
  --minMAC 2 \
  --test additive \
  --out /data/results/additive_single_variant

tr " " "\t" < "$WORKDIR/Results/additive_single_variant_T1R.regenie" > "$WORKDIR/Results/temp_file"
sed -i "s/\$/\tT1R_additive/" "$WORKDIR/Results/temp_file"
mv "$WORKDIR/Results/temp_file" "$WORKDIR/Results/additive_single_variant_T1R.regenie"

# Step 4: Gene-wise null model and test
PHENOS=("T1R" "T1R_PNAS" "T1R_New")

for PHENO in "${PHENOS[@]}"; do
  # Extract polymorphic variants
  plink2 \
    --pfile $WORKDIR/Data/Neuro_panel_samples \
    --keep $WORKDIR/Data/T1R_samples_${PHENO}.txt \
    --mac 2 \
    --max-alleles 2 \
    --write-snplist \
    --out $WORKDIR/Data/snps_pass_${PHENO}

  # Regenie Step 1 (null model)
  docker run --rm -v "$WORKDIR/Data:/data/input" -v "$WORKDIR/Results:/data/results" \
    ghcr.io/rgcgithub/regenie/regenie:v3.6.gz /usr/local/bin/regenie \
    --pgen /data/input/Neuro_panel_samples \
    --phenoFile /data/input/T1R_Neuro_samples.pheno \
    --extract /data/input/snps_pass_${PHENO}.snplist \
    --covarFile /data/input/T1R_Neuro_samples.cov \
    --phenoCol $PHENO \
    --bt \
    --bsize 5000 \
    --strict \
    --threads 8 \
    --out /data/results/STEP_1_${PHENO}

  # Gene-wise Step 2
  docker run --rm -v "$WORKDIR/Data:/data/input" -v "$WORKDIR/Results:/data/results" \
    ghcr.io/rgcgithub/regenie/regenie:v3.6.gz /usr/local/bin/regenie \
    --pgen /data/input/Neuro_panel_samples \
    --phenoFile /data/input/T1R_Neuro_samples_Borderline.pheno \
    --set-list /data/input/T1R_Neuro_samples.set \
    --anno-file /data/input/T1R_Neuro_samples.anno \
    --mask-def /data/input/T1R_Neuro_samples.mask \
    --pred /data/results/STEP_1_${PHENO}_pred.list \
    --step 2 \
    --phenoCol $PHENO \
    --bt \
    --write-samples \
    --bsize 5000 \
    --af-cc \
    --minMAC 1 \
    --rgc-gene-p \
    --write-mask-snplist \
    --vc-MACthr 1 \
    --vc-maxAAF 0.05 \
    --aaf-bins 0.01,0.05 \
    --test "additive" \
    --vc-tests skato-acat \
    --threads 4 \
    --out /data/results/Neuro_panel_samples_additive_Genewise_${PHENO}

  tr " " "\t" < "$WORKDIR/Results/Neuro_panel_samples_additive_Genewise_${PHENO}.regenie" > "$WORKDIR/Results/temp_file"
  sed -i "s/\$/\t${PHENO}/" "$WORKDIR/Results/temp_file"
  mv "$WORKDIR/Results/temp_file" "$WORKDIR/Results/Neuro_panel_samples_additive_Genewise_${PHENO}.regenie"
  mv "$WORKDIR/Results/Neuro_panel_samples_additive_Genewise_masks.snplist" "$WORKDIR/Results/Neuro_panel_samples_additive_Genewise_${PHENO}_masks.snplist"
done
