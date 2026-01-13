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
# - Ensure Conda environments 'fastqc', 'multiqc', 'trimmomatic', and 'qiime2-amplicon-2024.10' are available.
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
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_FILE"

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

# =======================================================================
# STEP 1: Quality Check (Raw Data)
# =======================================================================
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

# =======================================================================
# STEP 2: Adapter Trimming & Filtering
# =======================================================================
echo "=========================================================="
echo "STEP 2: Adapter Trimming & Filtering"
echo "=========================================================="

conda activate trimmomatic

echo "Running Trimmomatic..."
for FILE_R1 in "$RAW_DATA_DIR"/*_R1_001.fastq.gz; do
    FILENAME=$(basename "$FILE_R1")
    BASE_NAME=${FILENAME%_L001_R1_001.fastq.gz}
    
    # Check if this file is in our mapping list
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
    trimmomatic PE -threads $THREADS -phred33 \
        "$FILE_R1" "$FILE_R2" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100

    # Add to manifest
    echo -e "${SAMPLE_ID}\t${R1_PAIRED}\t${R2_PAIRED}" >> "$MANIFEST_FILE"
done

conda deactivate

# =======================================================================
# STEP 3: Post-Clean Quality Check
# =======================================================================
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

# Verify manifest
echo ""
echo "Generated manifest:"
cat "$MANIFEST_FILE"
echo ""
echo "Number of samples in manifest: $(wc -l < "$MANIFEST_FILE")"

# =======================================================================
# STEP 4: QIIME2 Analysis
# =======================================================================
echo "=========================================================="
echo "STEP 4: QIIME2 Analysis"
echo "=========================================================="

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

# Suppress Python warnings
export PYTHONWARNINGS="ignore"

# Use project directory for temp files
export TMPDIR="${BASE_DIR}/tmp"
mkdir -p "$TMPDIR"
echo "Using TMPDIR: $TMPDIR"

# 4.1 Import Data --------------------------------------------------------
echo "Importing data into QIIME2..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "✓ Data imported. Check visual/demux.qzv at https://view.qiime2.org"

# 4.2 Denoising with DADA2 ----------------------------------------------
echo ""
echo "Denoising with DADA2 (generating ASVs)..."
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

echo "✓ DADA2 denoising complete (ASVs generated)"

# 4.3 Filtering Contaminants (Negative Controls) ------------------------
echo ""
echo "Removing Contaminants based on Negative Controls..."

# Step 1: Extract ASV IDs present in negative controls
qiime feature-table filter-samples \
  --i-table "${QIIME_DIR}/core/table.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Negative_Control'" \
  --o-filtered-table "${QIIME_DIR}/core/neg-controls-table.qza"

# Step 1b: Check if negative control has any sequences
qiime feature-table summarize \
  --i-table "${QIIME_DIR}/core/neg-controls-table.qza" \
  --o-visualization "${QIIME_DIR}/visual/neg-controls-summary.qzv"

# Step 2: Get list of ASV IDs from negative controls
qiime tools export \
  --input-path "${QIIME_DIR}/core/neg-controls-table.qza" \
  --output-path "${QIIME_DIR}/export/neg-controls"

biom convert \
  -i "${QIIME_DIR}/export/neg-controls/feature-table.biom" \
  -o "${QIIME_DIR}/export/neg-controls/feature-table.tsv" \
  --to-tsv

# Extract ASV IDs (skip header and first comment line)
tail -n +3 "${QIIME_DIR}/export/neg-controls/feature-table.tsv" | \
  cut -f1 > "${QIIME_DIR}/export/neg-controls/contamination_ids.txt"

# Step 2b: Check if we found contaminant ASVs
NUM_CONTAMINANTS=$(wc -l < "${QIIME_DIR}/export/neg-controls/contamination_ids.txt")
echo "Found ${NUM_CONTAMINANTS} contaminant ASVs in negative controls"

if [ "$NUM_CONTAMINANTS" -eq 0 ]; then
  echo "WARNING: No contaminants found in negative control!"
  echo "Skipping decontamination step."
  
  # Just copy the original table
  cp "${QIIME_DIR}/core/table.qza" "${QIIME_DIR}/core/table-decontam.qza"
  cp "${QIIME_DIR}/core/rep-seqs.qza" "${QIIME_DIR}/core/rep-seqs-clean.qza"
else
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
fi

# Step 5: Remove Negative Control samples from the table
qiime feature-table filter-samples \
  --i-table "${QIIME_DIR}/core/table-decontam.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Sample'" \
  --o-filtered-table "${QIIME_DIR}/core/table-final.qza"

echo "✓ Decontamination complete. Final table: table-final.qza"

# 4.4 Build Phylogenetic Tree --------------------------------------------
echo ""
echo "Building phylogenetic tree..."

# Create temp directory in writable location
MAFFT_TMPDIR="${BASE_DIR}/tmp_mafft"
mkdir -p "$MAFFT_TMPDIR"
export TMPDIR="$MAFFT_TMPDIR"
export TMP="$MAFFT_TMPDIR"

# Try align-to-tree-mafft-fasttree first (faster, standard)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME_DIR}/core/rep-seqs-clean.qza" \
  --p-n-threads 1 \
  --o-alignment "${QIIME_DIR}/core/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME_DIR}/core/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME_DIR}/core/rooted-tree.qza" 2>&1

TREE_STATUS=$?

# If MAFFT fails, try alternative with parttree (for large datasets)
if [ $TREE_STATUS -ne 0 ] || [ ! -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "WARNING: Standard MAFFT alignment failed. Trying with --p-parttree..."
  
  # Clean up failed files
  rm -f "${QIIME_DIR}/core/aligned-rep-seqs.qza"
  rm -f "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza"
  rm -f "${QIIME_DIR}/core/unrooted-tree.qza"
  rm -f "${QIIME_DIR}/core/rooted-tree.qza"
  
  # Retry with parttree
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${QIIME_DIR}/core/rep-seqs-clean.qza" \
    --p-n-threads 1 \
    --p-parttree \
    --o-alignment "${QIIME_DIR}/core/aligned-rep-seqs.qza" \
    --o-masked-alignment "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza" \
    --o-tree "${QIIME_DIR}/core/unrooted-tree.qza" \
    --o-rooted-tree "${QIIME_DIR}/core/rooted-tree.qza" 2>&1
  
  TREE_STATUS=$?
fi

# Final check and summary
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "✓ Phylogenetic tree created successfully!"
  qiime tools peek "${QIIME_DIR}/core/rooted-tree.qza"
else
  echo "✗ Phylogenetic tree creation FAILED"
  echo "  Continuing without Faith PD (tree-based metrics will be skipped)"
fi

# Cleanup temp directory
rm -rf "$MAFFT_TMPDIR"

# 4.5 Rarefaction & Diversity Metrics ------------------------------------
echo ""
echo "Running Rarefaction Analysis..."

# Remove old core-metrics results if they exist
if [ -d "${QIIME_DIR}/core-metrics-results" ]; then
  rm -rf "${QIIME_DIR}/core-metrics-results"
fi

# Export table to calculate depths
qiime tools export \
  --input-path "${QIIME_DIR}/core/table-final.qza" \
  --output-path "${QIIME_DIR}/export/table-final-temp"

biom convert \
  -i "${QIIME_DIR}/export/table-final-temp/feature-table.biom" \
  -o "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" \
  --to-tsv

MIN_DEPTH=1000

# Calculate sampling depth (10th percentile to be conservative)
SAMPLING_DEPTH=$(tail -n +3 "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" | \
  awk -F'\t' '
    NR==1 {for(i=2; i<=NF; i++) header[i]=$i}
    {for(i=2; i<=NF; i++) sum[header[i]]+=$i}
    END {
      count = 0
      for(sample in sum) {
        if(sum[sample] > '$MIN_DEPTH') {
          depths[count++] = sum[sample]
        }
      }
      if(count > 0) {
        asort(depths)
        idx = int(count * 0.1)
        if(idx < 1) idx = 1
        print depths[idx]
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

# Create output directory
mkdir -p "${QIIME_DIR}/core-metrics-results"

# Check if phylogenetic tree exists
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "Using phylogenetic tree for diversity analysis (includes Faith PD)..."
  
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME_DIR}/core/rooted-tree.qza" \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${QIIME_DIR}/core-metrics-results" \
    --verbose
    
  METRICS_STATUS=$?
  
  if [ $METRICS_STATUS -ne 0 ]; then
    echo "WARNING: Phylogenetic metrics failed. Falling back to non-phylogenetic metrics..."
    rm -rf "${QIIME_DIR}/core-metrics-results"
    mkdir -p "${QIIME_DIR}/core-metrics-results"
    
    qiime diversity core-metrics \
      --i-table "${QIIME_DIR}/core/table-final.qza" \
      --p-sampling-depth "$SAMPLING_DEPTH" \
      --m-metadata-file "$METADATA_FILE" \
      --output-dir "${QIIME_DIR}/core-metrics-results"
  fi
else
  echo "WARNING: No phylogenetic tree found. Using non-phylogenetic metrics only..."
  
  qiime diversity core-metrics \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${QIIME_DIR}/core-metrics-results"
fi

echo "✓ Diversity metrics calculated!"

# 4.6 Generate Rarefaction Curves ----------------------------------------
echo ""
echo "Generating rarefaction curves..."

# Calculate max depth (for rarefaction curve x-axis)
MAX_DEPTH=$(tail -n +3 "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" | \
  awk -F'\t' '
    NR==1 {for(i=2; i<=NF; i++) header[i]=$i}
    {for(i=2; i<=NF; i++) sum[header[i]]+=$i}
    END {
      max = 0
      for(sample in sum) {
        if(sum[sample] > max) max = sum[sample]
      }
      print max
    }')

echo "Max depth: ${MAX_DEPTH}"

# Rarefaction curves for observed ASVs
qiime diversity alpha-rarefaction \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-min-depth 1 \
  --p-max-depth "$MAX_DEPTH" \
  --p-steps 20 \
  --m-metadata-file "$METADATA_FILE" \
  --o-visualization "${QIIME_DIR}/visual/rarefaction-curves.qzv"

echo "✓ Rarefaction curves generated: rarefaction-curves.qzv"

# Rarefaction curves with Faith PD (if tree exists)
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "Generating Faith PD rarefaction curves..."
  
  qiime diversity alpha-rarefaction \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --i-phylogeny "${QIIME_DIR}/core/rooted-tree.qza" \
    --p-min-depth 1 \
    --p-max-depth "$MAX_DEPTH" \
    --p-steps 20 \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv"
  
  echo "✓ Faith PD rarefaction curves generated: rarefaction-curves-phylogenetic.qzv"
fi

echo ""
echo "========== RAREFACTION SUMMARY =========="
echo "Sampling depth used: $SAMPLING_DEPTH reads"
echo "Max depth in dataset: $MAX_DEPTH reads"
echo "Visualizations:"
echo "  - rarefaction-curves.qzv (Observed ASVs + Shannon)"
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "  - rarefaction-curves-phylogenetic.qzv (Faith PD)"
fi
echo "=========================================="

# 4.7 Taxonomy Classification --------------------------------------------
echo ""
echo "Assigning Taxonomy..."
CLASSIFIER_PATH="/scratch_vol0/fungi/dugong_microbiome/05_QIIME2/silva-138.2-ssu-nr99-341f-805r-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_PATH" \
  --i-reads "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --o-classification "$QIIME_DIR/core/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/taxonomy.qza" \
  --o-visualization "$QIIME_DIR/visual/taxonomy.qzv"

echo "✓ Taxonomy assigned"

# =======================================================================
# STEP 5: EXPORTS FOR DOWNSTREAM ANALYSES
# =======================================================================
echo ""
echo "=========================================================="
echo "STEP 5: Exporting Results"
echo "=========================================================="

EXPORT_DIR="${QIIME_DIR}/export"
mkdir -p "${EXPORT_DIR}"

# 5.1 Final feature table (ASV abundance per sample) --------------------
echo "Exporting feature table..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/table-final.qza" \
  --output-path "${EXPORT_DIR}/feature_table"

biom convert \
  -i "${EXPORT_DIR}/feature_table/feature-table.biom" \
  -o "${EXPORT_DIR}/feature_table/feature-table.tsv" \
  --to-tsv

# Rename OTU to ASV for clarity
sed -i 's/#OTU ID/#ASV_ID/g' "${EXPORT_DIR}/feature_table/feature-table.tsv" 2>/dev/null || true

echo "✓ Feature table exported"

# 5.2 Representative sequences (ASV sequences in fasta) ------------------
echo "Exporting representative sequences..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/rep-seqs-clean.qza" \
  --output-path "${EXPORT_DIR}/rep_seqs"

echo "✓ Representative sequences exported"

# 5.3 Taxonomy table -----------------------------------------------------
echo "Exporting taxonomy..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/taxonomy.qza" \
  --output-path "${EXPORT_DIR}/taxonomy"

qiime metadata tabulate \
  --m-input-file "${QIIME_DIR}/core/taxonomy.qza" \
  --o-visualization "${QIIME_DIR}/visual/taxonomy-table.qzv"

echo "✓ Taxonomy exported"

# 5.4 DADA2 stats --------------------------------------------------------
echo "Exporting DADA2 stats..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/dada2-stats.qza" \
  --output-path "${EXPORT_DIR}/dada2_stats"

qiime metadata tabulate \
  --m-input-file "${QIIME_DIR}/core/dada2-stats.qza" \
  --o-visualization "${QIIME_DIR}/visual/dada2-stats.qzv"

echo "✓ DADA2 stats exported"

# 5.5 Diversity indices from core-metrics --------------------------------
echo "Exporting diversity indices..."

CORE_METRICS_DIR="${QIIME_DIR}/core-metrics-results"
mkdir -p "${EXPORT_DIR}/diversity"

# Observed ASVs
if [ -f "${CORE_METRICS_DIR}/observed_features_vector.qza" ]; then
  qiime tools export \
    --input-path "${CORE_METRICS_DIR}/observed_features_vector.qza" \
    --output-path "${EXPORT_DIR}/diversity/observed_ASVs"
fi

# Shannon
if [ -f "${CORE_METRICS_DIR}/shannon_vector.qza" ]; then
  qiime tools export \
    --input-path "${CORE_METRICS_DIR}/shannon_vector.qza" \
    --output-path "${EXPORT_DIR}/diversity/shannon"
fi

# Evenness
if [ -f "${CORE_METRICS_DIR}/evenness_vector.qza" ]; then
  qiime tools export \
    --input-path "${CORE_METRICS_DIR}/evenness_vector.qza" \
    --output-path "${EXPORT_DIR}/diversity/evenness"
fi

# Faith PD (if tree was computed)
if [ -f "${CORE_METRICS_DIR}/faith_pd_vector.qza" ]; then
  qiime tools export \
    --input-path "${CORE_METRICS_DIR}/faith_pd_vector.qza" \
    --output-path "${EXPORT_DIR}/diversity/faith_pd"
  echo "✓ Faith PD exported"
fi

echo "✓ Diversity indices exported"

# 5.6 Merge ASV abundance + taxonomy -------------------------------------
echo "Merging ASV abundance and taxonomy..."

python3 << 'EOFPYTHON'
import pandas as pd
import os
import sys

qiime_dir = os.getenv("QIIME_DIR")
if not qiime_dir:
    print("ERROR: QIIME_DIR environment variable not set!")
    sys.exit(1)

export_dir = os.path.join(qiime_dir, "export")

try:
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
    output_file = os.path.join(export_dir, "ASV_abundance_taxonomy.tsv")
    merged.to_csv(output_file, sep="\t")
    print(f"✓ Merged table saved: {output_file}")
except Exception as e:
    print(f"ERROR in merging tables: {e}")
    sys.exit(1)
EOFPYTHON

# 5.7 Sample depths ------------------------------------------------------
echo "Calculating sample depths..."

python3 << 'EOFPYTHON'
import pandas as pd
import os
import sys

qiime_dir = os.getenv("QIIME_DIR")
if not qiime_dir:
    print("ERROR: QIIME_DIR environment variable not set!")
    sys.exit(1)

export_dir = os.path.join(qiime_dir, "export")

try:
    tab = pd.read_csv(
        os.path.join(export_dir, "feature_table", "feature-table.tsv"),
        sep="\t", comment="#", index_col=0
    )

    depths = tab.sum(axis=0)
    output_file = os.path.join(export_dir, "sample_read_depths_final.tsv")
    depths.to_csv(output_file, sep="\t", header=["reads"])
    print(f"✓ Sample depths saved: {output_file}")
except Exception as e:
    print(f"ERROR in calculating depths: {e}")
    sys.exit(1)
EOFPYTHON

# =======================================================================
# STEP 6: COMPREHENSIVE DIVERSITY INDICES CALCULATION
# =======================================================================
echo ""
echo "=========================================================="
echo "STEP 6: Calculating All Diversity Indices"
echo "=========================================================="

DIVERSITY_DIR="${QIIME_DIR}/diversity_indices"
mkdir -p "${DIVERSITY_DIR}"

echo "Calculating additional alpha diversity indices..."

# Simpson index
qiime diversity alpha \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-metric simpson \
  --o-alpha-diversity "${DIVERSITY_DIR}/simpson_vector.qza"

# Simpson evenness
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

# Gini index
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

echo "✓ All alpha diversity indices calculated"

# Export all indices to TSV ----------------------------------------------
echo ""
echo "Exporting all diversity indices to TSV format..."

mkdir -p "${EXPORT_DIR}/diversity_all"

# Function to export and rename columns
export_diversity_metric() {
    local metric_name=$1
    local qza_file=$2
    local output_name=$3
    
    if [ -f "$qza_file" ]; then
        qiime tools export \
          --input-path "$qza_file" \
          --output-path "${EXPORT_DIR}/diversity_all/${output_name}_temp"
        
        # Rename column header
        sed "1s/.*/sample-id\t${output_name}/" \
          "${EXPORT_DIR}/diversity_all/${output_name}_temp/alpha-diversity.tsv" > \
          "${EXPORT_DIR}/diversity_all/${output_name}.tsv"
        
        rm -rf "${EXPORT_DIR}/diversity_all/${output_name}_temp"
    else
        echo "Warning: ${qza_file} not found, skipping ${metric_name}..."
    fi
}

# Export core metrics
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

# Merge all indices into comprehensive table -----------------------------
echo ""
echo "Merging all diversity indices into a single table..."

python3 << 'EOFPYTHON'
import pandas as pd
import os
import sys
from pathlib import Path

qiime_dir = os.getenv("QIIME_DIR")
if not qiime_dir:
    print("ERROR: QIIME_DIR environment variable not set!")
    sys.exit(1)

diversity_dir = Path(qiime_dir) / "export" / "diversity_all"

# List all TSV files
tsv_files = sorted(diversity_dir.glob("*.tsv"))

if not tsv_files:
    print("ERROR: No diversity TSV files found!")
    sys.exit(1)

print(f"Found {len(tsv_files)} diversity files to merge")

try:
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

except Exception as e:
    print(f"ERROR in merging diversity tables: {e}")
    sys.exit(1)

EOFPYTHON

# =======================================================================
# PIPELINE COMPLETION
# =======================================================================
echo ""
echo "=========================================================="
echo "Pipeline Completed Successfully!"
echo "=========================================================="
echo ""
echo "Main outputs:"
echo "  - Feature table: ${EXPORT_DIR}/feature_table/feature-table.tsv"
echo "  - ASV sequences: ${EXPORT_DIR}/rep_seqs/dna-sequences.fasta"
echo "  - Taxonomy: ${EXPORT_DIR}/taxonomy/taxonomy.tsv"
echo "  - ASV + Taxonomy merged: ${EXPORT_DIR}/ASV_abundance_taxonomy.tsv"
echo "  - All diversity indices: ${EXPORT_DIR}/diversity_indices_all.tsv"
echo "  - Rarefaction curves: ${QIIME_DIR}/visual/rarefaction-curves.qzv"
if [ -f "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv" ]; then
  echo "  - Faith PD rarefaction: ${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv"
fi
echo ""
echo "All outputs are in: $QIIME_DIR"
echo "=========================================================="
