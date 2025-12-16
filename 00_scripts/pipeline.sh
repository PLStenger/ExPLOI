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
echo "Creating Metadata File..."
cat <<EOF > "$METADATA_FILE"
sample-id	group
17_Jourand	Negative_Control
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
echo "sample-id,forward-absolute-filepath,reverse-absolute-filepath" > "$MANIFEST_FILE"

# Mapping logic for renaming and manifest creation
declare -A SAMPLES
SAMPLES=(
    ["17_Jourand_S110"]="17_Jourand"
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

    # Add to manifest
    echo "${SAMPLE_ID},${R1_PAIRED},${R2_PAIRED}" >> "$MANIFEST_FILE"

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


# --- 4. QIIME2 Pipeline ---
echo "=========================================================="
echo "STEP 4: QIIME2 Analysis"
echo "=========================================================="

# Activate QIIME2 environment
# Ensure proper environment name (e.g., qiime2-2021.4 or newer)
source activate qiime2-2021.4 || conda activate qiime2-2021.4

# Set TMPDIR to handle potential space issues
export TMPDIR="$BASE_DIR/tmp"
mkdir -p "$TMPDIR"

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

# 4.3 Filtering Contaminants (Negative Controls)
echo "Removing Contaminants based on Negative Controls..."

# Identify ASVs present in Negative Controls
qiime feature-table filter-samples \
  --i-table "$QIIME_DIR/core/table.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Negative_Control'" \
  --o-filtered-table "$QIIME_DIR/core/neg-controls-table.qza"

# Extract sequences found in Negative Controls
qiime feature-table filter-seqs \
  --i-data "$QIIME_DIR/core/rep-seqs.qza" \
  --i-table "$QIIME_DIR/core/neg-controls-table.qza" \
  --o-filtered-data "$QIIME_DIR/core/neg-control-seqs.qza"

# Remove these exact sequences from the main table (Strict decontamination)
# Alternatively, you could use 'exclude-seqs' with alignment if you suspect variations
qiime feature-table filter-features \
  --i-table "$QIIME_DIR/core/table.qza" \
  --i-data "$QIIME_DIR/core/neg-control-seqs.qza" \
  --p-exclude-ids \
  --o-filtered-table "$QIIME_DIR/core/table-clean.qza"

# Filter representative sequences to match cleaned table
qiime feature-table filter-seqs \
  --i-data "$QIIME_DIR/core/rep-seqs.qza" \
  --i-table "$QIIME_DIR/core/table-clean.qza" \
  --o-filtered-data "$QIIME_DIR/core/rep-seqs-clean.qza"

# Remove Negative Control Samples from the table now that they've served their purpose
qiime feature-table filter-samples \
  --i-table "$QIIME_DIR/core/table-clean.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Sample'" \
  --o-filtered-table "$QIIME_DIR/core/table-final.qza"

# 4.4 Rarefaction (Smart Auto-Depth)
echo "Running Rarefaction..."

# Export table to get sample depths
qiime tools export \
  --input-path "$QIIME_DIR/core/table-final.qza" \
  --output-path "$QIIME_DIR/export/table-final"

biom convert -i "$QIIME_DIR/export/table-final/feature-table.biom" \
             -o "$QIIME_DIR/export/table-final/feature-table.tsv" \
             --to-tsv

# Calculate smart sampling depth using python
# Logic: Find the median depth, or exclude only extremely low outliers (e.g., < 1000 reads)
# Here we set a minimum threshold. If a sample has fewer reads, it will be dropped.
# We aim to keep maximum samples.
MIN_DEPTH=1000 # Minimum acceptable reads per sample

# Get depth of the sample with the LOWEST count above the threshold to maximize sample retention
# This is a robust way to automate depth selection without losing samples unnecessarily
SAMPLING_DEPTH=$(grep -v "^#" "$QIIME_DIR/export/table-final/feature-table.tsv" | \
    awk -F'\t' '{for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) if(sum[i] > '$MIN_DEPTH') print sum[i]}' | \
    sort -n | head -1)

echo "Selected Sampling Depth: $SAMPLING_DEPTH"

if [ -z "$SAMPLING_DEPTH" ]; then
    echo "Error: Could not calculate sampling depth. Check if samples have enough reads."
    exit 1
fi

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$QIIME_DIR/core/rooted-tree.qza" \
  --i-table "$QIIME_DIR/core/table-final.qza" \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file "$METADATA_FILE" \
  --output-dir "$QIIME_DIR/core-metrics-results"

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

 # Join feature table and taxonomy on ASV IDs
 # Creates: ASV, taxonomy, confidence, Sample1, Sample2, ...
 python - << 'EOF'
 import pandas as pd
 import os

 export_dir = "${EXPORT_DIR}"

 # Load abundance table
 tab = pd.read_csv(
     os.path.join(export_dir, "feature_table", "feature-table.tsv"),
     sep="\t", comment="#", index_col=0
 )

 # Load taxonomy
 tax = pd.read_csv(
     os.path.join(export_dir, "taxonomy", "taxonomy.tsv"),
     sep="\t", comment="#"
 )
 tax = tax.rename(columns={"Feature ID": "ASV_ID", "Taxon": "taxonomy", "Confidence": "confidence"})
 tax = tax.set_index("ASV_ID")

 # Merge
 merged = tax.join(tab, how="inner")
 merged.to_csv(os.path.join(export_dir, "ASV_abundance_taxonomy.tsv"), sep="\t")
 EOF

 # 5.7 Sample depths (total reads per sample after all filters) ----------
 

 python - << 'EOF'
 import pandas as pd
 import os

 export_dir = "${EXPORT_DIR}"

 tab = pd.read_csv(
     os.path.join(export_dir, "feature_table", "feature-table.tsv"),
     sep="\t", comment="#", index_col=0
 )

 depths = tab.sum(axis=0)
 depths.to_csv(os.path.join(export_dir, "sample_read_depths_final.tsv"), sep="\t", header=["reads"])
 EOF

 echo "All export files are in: ${EXPORT_DIR}"


echo "Pipeline Completed Successfully!"
echo "Outputs are located in $QIIME_DIR"
