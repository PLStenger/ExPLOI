#!/bin/bash

# ========================================================================================================
# ExPLOI Project - 16S Metabarcoding Pipeline
# ========================================================================================================
#
# This script performs the complete metabarcoding analysis pipeline:
# 1. Quality Check (FastQC/MultiQC)
# 2. Adapter Trimming & Filtering (Trimmomatic)
# 3. QIIME2 Pipeline (Import, Denoise, Taxonomy, Rarefaction)
#
# CLUSTER CONFIGURATION:
# - Ensure Conda environments 'fastqc', 'multiqc', 'trimmomatic', and 'qiime2-2021.4' are available.
# - Adjust THREADS variable based on your cluster allocation.
#
# ========================================================================================================

# --- Configuration Variables ---
BASE_DIR="/nvme/bio/data_fungi/ExPLOI"
RAW_DATA_DIR="${BASE_DIR}/01_raw_data"
QC_DIR="${BASE_DIR}/02_quality_check"
CLEANED_DIR="${BASE_DIR}/03_cleaned_data"
POST_CLEAN_QC_DIR="${BASE_DIR}/04_post_clean_quality_check"
QIIME_DIR="${BASE_DIR}/05_QIIME2"
METADATA_FILE="${BASE_DIR}/metadata_ExPLOI.tsv"
MANIFEST_FILE="${BASE_DIR}/manifest_ExPLOI.txt"
ADAPTER_FILE="/nvme/bio/data_fungi/valormicro_nc/99_softwares/adapters/sequences.fasta"
THREADS=16 

# --- Create Directory Structure ---
mkdir -p "$QC_DIR" "$CLEANED_DIR" "$POST_CLEAN_QC_DIR" "$QIIME_DIR"
mkdir -p "$QIIME_DIR/core" "$QIIME_DIR/visual" "$QIIME_DIR/export"

# --- Create Metadata File ---

#17_Jourand	Negative_Control


echo "Creating Metadata File..."
cat <<EOF > "$METADATA_FILE"
sample-id	group
BL_PCR_Jourand	Negative_Control
T1_CO	Sample
T2_CO	Sample
T1_HO	Sample
T2_HO	Sample
T1_HE	Sample
T2_HE	Sample
T1_CIAN	Sample
T2_CIAN	Sample
T1_HIAN	Sample
MP1	Sample
MP2	Sample
MP3	Sample
MP4	Sample
MP5	Sample
MP6	Sample
MP7	Sample
EOF

# --- Create Manifest File ---
echo "Creating Manifest File..."
# IMPORTANT: Use TAB-separated format, not comma-separated
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_FILE"

#     ["17_Jourand_S110"]="17_Jourand"


# Mapping logic for renaming and manifest creation
declare -A SAMPLES
SAMPLES=(
    ["BL_PCR_Jourand_S151"]="BL_PCR_Jourand"
    ["BAR253117_S114"]="T1_CO"
    ["BAR253118_S115"]="T2_CO"
    ["BAR253119_S116"]="T1_HO"
    ["BAR253120_S117"]="T2_HO"
    ["BAR253121_S118"]="T1_HE"
    ["BAR253122_S119"]="T2_HE"
    ["BAR253123_S120"]="T1_CIAN"
    ["BAR253124_S121"]="T2_CIAN"
    ["BAR253125_S122"]="T1_HIAN"
    ["BAR253126_S123"]="MP1"
    ["BAR253127_S124"]="MP2"
    ["BAR253128_S125"]="MP3"
    ["BAR253129_S126"]="MP4"
    ["BAR253130_S127"]="MP5"
    ["BAR253131_S128"]="MP6"
    ["BAR253132_S129"]="MP7"
)

# --- 1. Quality Check (Raw Data) ---
echo "=========================================================="
echo "STEP 1: Quality Check on Raw Data"
echo "=========================================================="

eval "$(conda shell.bash hook)"
conda activate fastqc

echo "Running FastQC..."
fastqc -t $THREADS "$RAW_DATA_DIR"/*.fastq.gz -o "$QC_DIR"

conda deactivate
conda activate multiqc

echo "Running MultiQC..."
multiqc "$QC_DIR" -o "$QC_DIR"

conda deactivate

# --- 2. Trimmomatic Cleaning ---
echo "=========================================================="
echo "STEP 2: Adapter Trimming & Filtering"
echo "=========================================================="

conda activate trimmomatic

echo "Running Trimmomatic..."
# Loop through samples based on the associative array keys to process files
# Note: Using globbing to find files matching the pattern
for FILE_R1 in "$RAW_DATA_DIR"/*_R1_001.fastq.gz; do
    FILENAME=$(basename "$FILE_R1")
    BASE_NAME=${FILENAME%_L001_R1_001.fastq.gz} # Extract base name like BAR253117_S114
    
    # Check if this file is in our mapping list (optional safety)
    SAMPLE_ID=""
    for KEY in "${!SAMPLES[@]}"; do
        if [[ "$BASE_NAME" == *"$KEY"* ]]; then
            SAMPLE_ID="${SAMPLES[$KEY]}"
            break
        fi
    done

    if [ -z "$SAMPLE_ID" ]; then
        echo "Warning: No sample ID mapping found for $FILENAME. Skipping..."
        continue
    fi
    
    echo "Processing $SAMPLE_ID ($FILENAME)..."

    FILE_R2="${FILE_R1//_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    # Output filenames
    R1_PAIRED="$CLEANED_DIR/${SAMPLE_ID}_R1_paired.fastq.gz"
    R1_UNPAIRED="$CLEANED_DIR/${SAMPLE_ID}_R1_unpaired.fastq.gz"
    R2_PAIRED="$CLEANED_DIR/${SAMPLE_ID}_R2_paired.fastq.gz"
    R2_UNPAIRED="$CLEANED_DIR/${SAMPLE_ID}_R2_unpaired.fastq.gz"

    # Run Trimmomatic
    # Note: Ensure $ADAPTER_FILE exists. If generic Illumina adapters, standard Trimmomatic path usually works too.
    # If using custom adapter file, ensure it is created/downloaded.
    # Added SLIDINGWINDOW and MINLEN as requested.
    trimmomatic PE -threads $THREADS -phred33 \
        "$FILE_R1" "$FILE_R2" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100

    # Add to manifest (use TAB instead of comma)
    echo -e "${SAMPLE_ID}\t${R1_PAIRED}\t${R2_PAIRED}" >> "$MANIFEST_FILE"

done

conda deactivate

# --- 3. Post-Clean Quality Check ---
echo "=========================================================="
echo "STEP 3: Quality Check on Cleaned Data"
echo "=========================================================="

conda activate fastqc

echo "Running FastQC on paired cleaned reads..."
fastqc -t $THREADS "$CLEANED_DIR"/*_paired.fastq.gz -o "$POST_CLEAN_QC_DIR"

conda deactivate
conda activate multiqc

echo "Running MultiQC..."
multiqc "$POST_CLEAN_QC_DIR" -o "$POST_CLEAN_QC_DIR"

conda deactivate

# After the Trimmomatic loop, before QIIME2 import:
echo "Generated manifest:"
cat "$MANIFEST_FILE"
echo ""
echo "Number of samples in manifest: $(wc -l < "$MANIFEST_FILE")"


# --- 4. QIIME2 Pipeline ---
echo "=========================================================="
echo "STEP 4: QIIME2 Analysis"
echo "=========================================================="

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10
#conda activate qiime2-2021.4

#export PYTHONPATH="${PYTHONPATH}:/scratch_vol0/fungi/.local/lib/python3.9/site-packages/"
#echo $PYTHONPATH

# I'm doing this step in order to deal the no space left in cluster :
export TMPDIR='/scratch_vol0/fungi'
echo $TMPDIR

# Suppress Python warnings
export PYTHONWARNINGS="ignore"

# Fix TMPDIR for MAFFT
export TMPDIR='/scratch_vol0/fungi/tmp_ExPLOI'
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
echo "Using TMPDIR: $TMPDIR"

# 4.1 Import Data
echo "Importing data into QIIME2..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "Check visual/demux.qzv at https://view.qiime2.org to determine trim lengths."
# Assuming standard V3-V4 lengths or high quality, we set safe defaults (0 trimming). 
# ADJUST THESE VALUES if quality drops significantly at ends.

# 4.2 Denoising with DADA2
echo "Denoising with DADA2..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$QIIME_DIR/core/demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads $THREADS \
  --o-table "$QIIME_DIR/core/table.qza" \
  --o-representative-sequences "$QIIME_DIR/core/rep-seqs.qza" \
  --o-denoising-stats "$QIIME_DIR/core/dada2-stats.qza"

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/dada2-stats.qza" \
  --o-visualization "$QIIME_DIR/visual/dada2-stats.qzv"

# 4.2b Generate phylogenetic tree (required for diversity metrics)
echo "Building phylogenetic tree..."

# Use single thread to avoid memory issues on cluster
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME_DIR}/core/rep-seqs.qza" \
  --p-n-threads 1 \
  --o-alignment "${QIIME_DIR}/core/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME_DIR}/core/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME_DIR}/core/rooted-tree.qza" \
  --verbose

# If MAFFT still fails, try with a smaller dataset first (skip tree-based metrics)
if [ ! -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "WARNING: Phylogenetic tree creation failed. Skipping tree-based diversity metrics."
  echo "You can still analyze data without Faith PD and UniFrac distances."
fi

# 4.3 Filtering Contaminants (Negative Controls)
echo "Removing Contaminants based on Negative Controls..."

# Step 1: Extract ASV IDs present in negative controls
qiime feature-table filter-samples \
  --i-table "${QIIME_DIR}/core/table.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Negative_Control'" \
  --o-filtered-table "${QIIME_DIR}/core/neg-controls-table.qza"

# Step 2: Get list of ASV IDs from negative controls
qiime tools export \
  --input-path "${QIIME_DIR}/core/neg-controls-table.qza" \
  --output-path "${QIIME_DIR}/export/neg-controls"

biom convert \
  -i "${QIIME_DIR}/export/neg-controls/feature-table.biom" \
  -o "${QIIME_DIR}/export/neg-controls/feature-table.tsv" \
  --to-tsv

# Extract just the ASV IDs (first column, skip header)
tail -n +2 "${QIIME_DIR}/export/neg-controls/feature-table.tsv" | \
  cut -f1 > "${QIIME_DIR}/export/neg-controls/contamination_ids.txt"

# Step 3: Remove contaminating ASVs from main table
qiime feature-table filter-features \
  --i-table "${QIIME_DIR}/core/table.qza" \
  --m-metadata-file "${QIIME_DIR}/export/neg-controls/contamination_ids.txt" \
  --p-exclude-ids \
  --o-filtered-table "${QIIME_DIR}/core/table-decontam.qza"

# Step 4: Filter representative sequences to match
qiime feature-table filter-seqs \
  --i-data "${QIIME_DIR}/core/rep-seqs.qza" \
  --i-table "${QIIME_DIR}/core/table-decontam.qza" \
  --o-filtered-data "${QIIME_DIR}/core/rep-seqs-clean.qza"

# Step 5: Remove Negative Control samples from the table
qiime feature-table filter-samples \
  --i-table "${QIIME_DIR}/core/table-decontam.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Sample'" \
  --o-filtered-table "${QIIME_DIR}/core/table-final.qza"

echo "Decontamination complete. Final table: table-final.qza"


# 4.4 Rarefaction (Smart Auto-Depth)
echo "Running Rarefaction..."

# Export table to calculate depths
qiime tools export \
  --input-path "${QIIME_DIR}/core/table-final.qza" \
  --output-path "${QIIME_DIR}/export/table-final-temp"

biom convert \
  -i "${QIIME_DIR}/export/table-final-temp/feature-table.biom" \
  -o "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" \
  --to-tsv

MIN_DEPTH=1000

SAMPLING_DEPTH=$(tail -n +3 "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" | \
  awk -F'\t' '
    NR==1 {for(i=2; i<=NF; i++) header[i]=$i}
    {for(i=2; i<=NF; i++) sum[header[i]]+=$i}
    END {
      for(sample in sum) {
        if(sum[sample] > '$MIN_DEPTH') print sum[sample]
      }
    }' | \
  sort -n | head -1)

if [ -z "$SAMPLING_DEPTH" ] || [ "$SAMPLING_DEPTH" -lt "$MIN_DEPTH" ]; then
  echo "WARNING: No samples above $MIN_DEPTH reads. Using median depth instead."
  SAMPLING_DEPTH=$(tail -n +3 "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" | \
    awk -F'\t' '
      NR==1 {for(i=2; i<=NF; i++) header[i]=$i}
      {for(i=2; i<=NF; i++) sum[header[i]]+=$i}
      END {for(sample in sum) print sum[sample]}' | \
    sort -n | awk '{a[NR]=$0} END {print a[int(NR/2)]}')
fi

echo "Selected Sampling Depth: $SAMPLING_DEPTH"

# Check if phylogenetic tree exists
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "Using phylogenetic tree for diversity analysis..."
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME_DIR}/core/rooted-tree.qza" \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${QIIME_DIR}/core-metrics-results"
else
  echo "WARNING: No phylogenetic tree found. Using non-phylogenetic metrics only..."
  qiime diversity core-metrics \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${QIIME_DIR}/core-metrics-results"
fi


# 4.5 Taxonomy Classification (Optional - requires trained classifier)
 echo "Assigning Taxonomy..."
 CLASSIFIER_PATH="/scratch_vol0/fungi/dugong_microbiome/05_QIIME2/silva-138.2-ssu-nr99-341f-805r-classifier.qza"
 qiime feature-classifier classify-sklearn \
   --i-classifier "$CLASSIFIER_PATH" \
   --i-reads "$QIIME_DIR/core/rep-seqs-clean.qza" \
   --o-classification "$QIIME_DIR/core/taxonomy.qza"

 qiime metadata tabulate \
   --m-input-file "$QIIME_DIR/core/taxonomy.qza" \
   --o-visualization "$QIIME_DIR/visual/taxonomy.qzv"
 
 # =======================================================================
 # STEP 5: EXPORTS FOR DOWNSTREAM ANALYSES
 # =======================================================================

 EXPORT_DIR="${QIIME_DIR}/export"
 mkdir -p "${EXPORT_DIR}"

 echo "Exporting core objects..."

 # 5.1 Final feature table (ASV abundance per sample) ---------------------

 # Export BIOM
 qiime tools export \
   --input-path "${QIIME_DIR}/core/table-final.qza" \
   --output-path "${EXPORT_DIR}/feature_table"

 # Convert BIOM to TSV
 biom convert \
   -i "${EXPORT_DIR}/feature_table/feature-table.biom" \
   -o "${EXPORT_DIR}/feature_table/feature-table.tsv" \
   --to-tsv

# Rename OTU to ASV for clarity in exported files
sed -i 's/#OTU ID/#ASV_ID/g' "${EXPORT_DIR}/feature_table/feature-table.tsv" 2>/dev/null || true

 # 5.2 Representative sequences (ASV sequences in fasta) -----------------

 qiime tools export \
   --input-path "${QIIME_DIR}/core/rep-seqs-clean.qza" \
   --output-path "${EXPORT_DIR}/rep_seqs"

 # 5.3 Taxonomy table -----------------------------------------------------

 # (Uncomment taxonomy steps in the main script if not already done)
 qiime tools export \
   --input-path "${QIIME_DIR}/core/taxonomy.qza" \
   --output-path "${EXPORT_DIR}/taxonomy"

 # Convert taxonomy.qza into a clean TSV with headers:
 # ASV_ID  Kingdom Phylum Class Order Family Genus Species
 qiime metadata tabulate \
   --m-input-file "${QIIME_DIR}/core/taxonomy.qza" \
   --o-visualization "${QIIME_DIR}/visual/taxonomy-table.qzv"

 # 5.5 DADA2 stats: reads before/after filtering per sample --------------

 qiime tools export \
   --input-path "${QIIME_DIR}/core/dada2-stats.qza" \
   --output-path "${EXPORT_DIR}/dada2_stats"

 # Convert to TSV
 qiime metadata tabulate \
   --m-input-file "${QIIME_DIR}/core/dada2-stats.qza" \
   --o-visualization "${QIIME_DIR}/visual/dada2-stats.qzv"

 # 5.6 Diversity indices from core-metrics -------------------------------

 CORE_METRICS_DIR="${QIIME_DIR}/core-metrics-results"
 mkdir -p "${EXPORT_DIR}/diversity"

 # Observed ASVs
 qiime tools export \
   --input-path "${CORE_METRICS_DIR}/observed_features_vector.qza" \
   --output-path "${EXPORT_DIR}/diversity/observed_ASVs"

 # Shannon
 qiime tools export \
   --input-path "${CORE_METRICS_DIR}/shannon_vector.qza" \
   --output-path "${EXPORT_DIR}/diversity/shannon"

 # Evenness
 qiime tools export \
   --input-path "${CORE_METRICS_DIR}/evenness_vector.qza" \
   --output-path "${EXPORT_DIR}/diversity/evenness"

 # Faith PD (if tree was computed)
 qiime tools export \
   --input-path "${CORE_METRICS_DIR}/faith_pd_vector.qza" \
   --output-path "${EXPORT_DIR}/diversity/faith_pd"
 
# 5.4 Merge ASV abundance + taxonomy into one table ---------------------

python3 << 'EOF'
import pandas as pd
import os

export_dir = os.getenv("QIIME_DIR") + "/export"

# Load abundance table
tab = pd.read_csv(
    os.path.join(export_dir, "feature_table", "feature-table.tsv"),
    sep="\t", comment="#", index_col=0
)

# Load taxonomy
tax = pd.read_csv(
    os.path.join(export_dir, "taxonomy", "taxonomy.tsv"),
    sep="\t", index_col=0
)

# Merge
merged = tax.join(tab, how="inner")
merged.to_csv(os.path.join(export_dir, "ASV_abundance_taxonomy.tsv"), sep="\t")
print(f"Merged table saved: {os.path.join(export_dir, 'ASV_abundance_taxonomy.tsv')}")
EOF

# 5.7 Sample depths (total reads per sample after all filters) ----------

python3 << 'EOF'
import pandas as pd
import os

export_dir = os.getenv("QIIME_DIR") + "/export"

tab = pd.read_csv(
    os.path.join(export_dir, "feature_table", "feature-table.tsv"),
    sep="\t", comment="#", index_col=0
)

depths = tab.sum(axis=0)
depths.to_csv(os.path.join(export_dir, "sample_read_depths_final.tsv"), sep="\t", header=["reads"])
print(f"Sample depths saved: {os.path.join(export_dir, 'sample_read_depths_final.tsv')}")
EOF


# =======================================================================
# STEP 6: COMPREHENSIVE DIVERSITY INDICES CALCULATION
# =======================================================================

echo "=========================================================="
echo "STEP 6: Calculating All Diversity Indices"
echo "=========================================================="

DIVERSITY_DIR="${QIIME_DIR}/diversity_indices"
mkdir -p "${DIVERSITY_DIR}"

# 6.1 Alpha diversity indices (beyond core-metrics) ----------------------

echo "Calculating additional alpha diversity indices..."

# Simpson index
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric simpson \
  --o-alpha-diversity "${DIVERSITY_DIR}/simpson_vector.qza"

# Inverse Simpson (Simpson's Diversity Index)
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric simpson_e \
  --o-alpha-diversity "${DIVERSITY_DIR}/simpson_evenness_vector.qza"

# Chao1 richness estimator
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric chao1 \
  --o-alpha-diversity "${DIVERSITY_DIR}/chao1_vector.qza"

# ACE richness estimator
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric ace \
  --o-alpha-diversity "${DIVERSITY_DIR}/ace_vector.qza"

# Good's coverage
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric goods_coverage \
  --o-alpha-diversity "${DIVERSITY_DIR}/goods_coverage_vector.qza"

# Fisher's alpha
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric fisher_alpha \
  --o-alpha-diversity "${DIVERSITY_DIR}/fisher_alpha_vector.qza"

# Berger-Parker dominance
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric berger_parker_d \
  --o-alpha-diversity "${DIVERSITY_DIR}/berger_parker_vector.qza"

# Dominance (Gini-Simpson)
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric gini_index \
  --o-alpha-diversity "${DIVERSITY_DIR}/gini_index_vector.qza"

# Brillouin's diversity index
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric brillouin_d \
  --o-alpha-diversity "${DIVERSITY_DIR}/brillouin_vector.qza"

# Strong's dominance
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric strong \
  --o-alpha-diversity "${DIVERSITY_DIR}/strong_vector.qza"

# McIntosh diversity
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric mcintosh_d \
  --o-alpha-diversity "${DIVERSITY_DIR}/mcintosh_d_vector.qza"

# McIntosh evenness
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric mcintosh_e \
  --o-alpha-diversity "${DIVERSITY_DIR}/mcintosh_e_vector.qza"

# Margalef richness index
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric margalef \
  --o-alpha-diversity "${DIVERSITY_DIR}/margalef_vector.qza"

# Menhinick richness index
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric menhinick \
  --o-alpha-diversity "${DIVERSITY_DIR}/menhinick_vector.qza"

echo "Alpha diversity indices calculated."

# 6.2 Export all indices to TSV -------------------------------------------

echo "Exporting all diversity indices to TSV format..."

mkdir -p "${EXPORT_DIR}/diversity_all"

# Function to export and rename columns for clarity
export_diversity_metric() {
    local metric_name=$1
    local qza_file=$2
    local output_name=$3
    
    if [ -f "$qza_file" ]; then
        qiime tools export \
          --input-path "$qza_file" \
          --output-path "${EXPORT_DIR}/diversity_all/${output_name}_temp"
        
        # Rename column header from generic to specific metric name
        sed "1s/.*/sample-id\t${output_name}/" \
          "${EXPORT_DIR}/diversity_all/${output_name}_temp/alpha-diversity.tsv" > \
          "${EXPORT_DIR}/diversity_all/${output_name}.tsv"
        
        rm -rf "${EXPORT_DIR}/diversity_all/${output_name}_temp"
    else
        echo "Warning: ${qza_file} not found, skipping ${metric_name}..."
    fi
}

# Export core metrics (from rarefaction)
export_diversity_metric "Observed ASVs" \
  "${QIIME_DIR}/core-metrics-results/observed_features_vector.qza" \
  "observed_asvs"

export_diversity_metric "Shannon" \
  "${QIIME_DIR}/core-metrics-results/shannon_vector.qza" \
  "shannon"

export_diversity_metric "Pielou Evenness" \
  "${QIIME_DIR}/core-metrics-results/evenness_vector.qza" \
  "pielou_evenness"

export_diversity_metric "Faith PD" \
  "${QIIME_DIR}/core-metrics-results/faith_pd_vector.qza" \
  "faith_pd"

# Export additional metrics
export_diversity_metric "Simpson" \
  "${DIVERSITY_DIR}/simpson_vector.qza" \
  "simpson"

export_diversity_metric "Simpson Evenness" \
  "${DIVERSITY_DIR}/simpson_evenness_vector.qza" \
  "simpson_evenness"

export_diversity_metric "Chao1" \
  "${DIVERSITY_DIR}/chao1_vector.qza" \
  "chao1"

export_diversity_metric "ACE" \
  "${DIVERSITY_DIR}/ace_vector.qza" \
  "ace"

export_diversity_metric "Goods Coverage" \
  "${DIVERSITY_DIR}/goods_coverage_vector.qza" \
  "goods_coverage"

export_diversity_metric "Fisher Alpha" \
  "${DIVERSITY_DIR}/fisher_alpha_vector.qza" \
  "fisher_alpha"

export_diversity_metric "Berger Parker" \
  "${DIVERSITY_DIR}/berger_parker_vector.qza" \
  "berger_parker"

export_diversity_metric "Gini Index" \
  "${DIVERSITY_DIR}/gini_index_vector.qza" \
  "gini_index"

export_diversity_metric "Brillouin" \
  "${DIVERSITY_DIR}/brillouin_vector.qza" \
  "brillouin"

export_diversity_metric "Strong" \
  "${DIVERSITY_DIR}/strong_vector.qza" \
  "strong"

export_diversity_metric "McIntosh D" \
  "${DIVERSITY_DIR}/mcintosh_d_vector.qza" \
  "mcintosh_d"

export_diversity_metric "McIntosh E" \
  "${DIVERSITY_DIR}/mcintosh_e_vector.qza" \
  "mcintosh_e"

export_diversity_metric "Margalef" \
  "${DIVERSITY_DIR}/margalef_vector.qza" \
  "margalef"

export_diversity_metric "Menhinick" \
  "${DIVERSITY_DIR}/menhinick_vector.qza" \
  "menhinick"

# 6.3 Merge all indices into a single comprehensive table -----------------

echo "Merging all diversity indices into a single table..."

python3 << 'EOFPYTHON'
import pandas as pd
import os
from pathlib import Path

diversity_dir = Path(os.getenv("QIIME_DIR")) / "export" / "diversity_all"

# List all TSV files
tsv_files = sorted(diversity_dir.glob("*.tsv"))

if not tsv_files:
    print("ERROR: No diversity TSV files found!")
    exit(1)

# Read first file as base
df_merged = pd.read_csv(tsv_files[0], sep="\t", index_col=0)

# Merge all other files
for tsv_file in tsv_files[1:]:
    try:
        df_temp = pd.read_csv(tsv_file, sep="\t", index_col=0)
        df_merged = df_merged.join(df_temp, how="outer")
    except Exception as e:
        print(f"Warning: Could not merge {tsv_file.name}: {e}")

# Sort by sample name
df_merged = df_merged.sort_index()

# Save comprehensive table
output_file = diversity_dir.parent / "diversity_indices_all.tsv"
df_merged.to_csv(output_file, sep="\t")

print(f"\n✓ Comprehensive diversity table saved:")
print(f"  {output_file}")
print(f"\n✓ Number of samples: {len(df_merged)}")
print(f"✓ Number of indices: {len(df_merged.columns)}")
print(f"\nColumns included:")
for col in df_merged.columns:
    print(f"  - {col}")

EOFPYTHON

echo ""
echo "=========================================================="
echo "All diversity indices calculated and exported!"
echo "Main output: ${EXPORT_DIR}/diversity_indices_all.tsv"
echo "=========================================================="



echo "Pipeline Completed Successfully!"
echo "Outputs are located in $QIIME_DIR"
