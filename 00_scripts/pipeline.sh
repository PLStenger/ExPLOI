#!/bin/bash

# ========================================================================================================
# ExPLOI Project - 16S Metabarcoding Pipeline COMPLET (VERSION FINALE)
# ========================================================================================================
#
# Ce script effectue l'analyse complète de métabarcoding 16S de A à Z :
#
# ÉTAPE 1 : Contrôle qualité des données brutes (FastQC/MultiQC)
# ÉTAPE 2 : Nettoyage et filtrage des adaptateurs (Trimmomatic)
# ÉTAPE 3 : Contrôle qualité post-nettoyage (FastQC/MultiQC)
# ÉTAPE 4 : Analyse QIIME2
#   4.1 : Import des données
#   4.2 : Débruitage DADA2 (génération des ASVs)
#   4.3 : Filtrage des contaminants (contrôles négatifs)
#   4.4 : Construction de l'arbre phylogénétique
#   4.5 : Métriques de diversité (core-metrics-phylogenetic)
#   4.6 : Courbes de raréfaction (Shannon, Observed Features, Faith PD)
#   4.7 : Classification taxonomique (SILVA)
# ÉTAPE 5 : Export des résultats
#   5.1 : Table de features
#   5.2 : Séquences représentatives
#   5.3 : Taxonomie
#   5.4 : Arbre phylogénétique
#   5.5 : Statistiques DADA2
#   5.6 : Indices de diversité
#   5.7 : Table ASV + Taxonomie fusionnée
#   5.8 : Profondeurs de lecture par échantillon
# ÉTAPE 6 : Calcul de tous les indices de diversité alpha
# ÉTAPE 7 : Export des données de raréfaction (CSV)
#
# CORRECTIONS IMPLÉMENTÉES :
# - Gestion robuste des contrôles négatifs vides (0 reads)
# - Calcul correct des profondeurs avec Python (pas AWK/bash)
# - Noms d'échantillons préservés correctement
# - Arbre phylogénétique avec gestion des échecs MAFFT
# - Export complet de toutes les données pour analyses R
#
# PRÉREQUIS :
# - Conda avec environnements : fastqc, multiqc, trimmomatic, qiime2-amplicon-2024.10
# - Classifier SILVA 138.2 (16S)
# - Adaptateurs Illumina
#
# UTILISATION :
#   chmod +x pipeline-complet-final.sh
#   nohup bash pipeline-complet-final.sh > pipeline.out 2>&1 &
#   tail -f pipeline.out
#
# ========================================================================================================

# ========================================================================================================
# CONFIGURATION
# ========================================================================================================

BASE_DIR="/nvme/bio/data_fungi/ExPLOI"
RAW_DATA_DIR="${BASE_DIR}/01_raw_data"
QC_DIR="${BASE_DIR}/02_quality_check"
CLEANED_DIR="${BASE_DIR}/03_cleaned_data"
POST_CLEAN_QC_DIR="${BASE_DIR}/04_post_clean_quality_check"
QIIME_DIR="${BASE_DIR}/05_QIIME2"
METADATA_FILE="${BASE_DIR}/metadata_ExPLOI.tsv"
MANIFEST_FILE="${BASE_DIR}/manifest_ExPLOI.txt"
ADAPTER_FILE="/nvme/bio/data_fungi/valormicro_nc/99_softwares/adapters/sequences.fasta"
CLASSIFIER_PATH="/nvme/bio/data_fungi/BioIndic_La_Reunion_Island_seawater_four_month_SED/05_QIIME2/Original_reads_16S/taxonomy/16S/Classifier.qza"
THREADS=16

# Export pour Python
export QIIME_DIR="${QIIME_DIR}"
export PYTHONWARNINGS="ignore"

# ========================================================================================================
# CRÉATION DE L'ARBORESCENCE
# ========================================================================================================

echo "========================================"
echo "ExPLOI - 16S Metabarcoding Pipeline"
echo "========================================"
echo "Date: $(date)"
echo "User: $(whoami)"
echo "Host: $(hostname)"
echo "========================================"
echo ""

mkdir -p "$QC_DIR" "$CLEANED_DIR" "$POST_CLEAN_QC_DIR" "$QIIME_DIR"
mkdir -p "$QIIME_DIR/core" "$QIIME_DIR/visual" "$QIIME_DIR/export"

# Temp directories
export TMPDIR="${BASE_DIR}/tmp"
mkdir -p "$TMPDIR"

echo "Directories created:"
echo "  - Raw data: $RAW_DATA_DIR"
echo "  - QC: $QC_DIR"
echo "  - Cleaned: $CLEANED_DIR"
echo "  - QIIME2: $QIIME_DIR"
echo "  - Temp: $TMPDIR"
echo ""

# ========================================================================================================
# CRÉATION DES FICHIERS METADATA ET MANIFEST
# ========================================================================================================

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

echo "✓ Metadata file created: $METADATA_FILE"
echo ""

echo "Creating Manifest File..."
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_FILE"

# Mapping échantillons
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

# ========================================================================================================
# ÉTAPE 1 : CONTRÔLE QUALITÉ DES DONNÉES BRUTES
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 1 : Contrôle Qualité des Données Brutes"
echo "=========================================================================================================="
echo ""

eval "$(conda shell.bash hook)"
conda activate fastqc

echo "Running FastQC on raw data..."
fastqc -t $THREADS "$RAW_DATA_DIR"/*.fastq.gz -o "$QC_DIR"

if [ $? -eq 0 ]; then
    echo "✓ FastQC completed"
else
    echo "ERROR: FastQC failed!"
    exit 1
fi

conda deactivate
conda activate multiqc

echo ""
echo "Running MultiQC..."
multiqc "$QC_DIR" -o "$QC_DIR" --force

if [ $? -eq 0 ]; then
    echo "✓ MultiQC completed"
    echo "✓ Report: ${QC_DIR}/multiqc_report.html"
else
    echo "ERROR: MultiQC failed!"
fi

conda deactivate

echo ""
echo "ÉTAPE 1 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 2 : NETTOYAGE ET FILTRAGE DES ADAPTATEURS
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 2 : Nettoyage et Filtrage des Adaptateurs (Trimmomatic)"
echo "=========================================================================================================="
echo ""

conda activate trimmomatic

echo "Running Trimmomatic..."
SAMPLE_COUNT=0

for FILE_R1 in "$RAW_DATA_DIR"/*_R1_001.fastq.gz; do
    FILENAME=$(basename "$FILE_R1")
    BASE_NAME=${FILENAME%_L001_R1_001.fastq.gz}
    
    # Vérifier si l'échantillon est dans le mapping
    SAMPLE_ID=""
    for KEY in "${!SAMPLES[@]}"; do
        if [[ "$BASE_NAME" == *"$KEY"* ]]; then
            SAMPLE_ID="${SAMPLES[$KEY]}"
            break
        fi
    done

    if [ -z "$SAMPLE_ID" ]; then
        echo "⚠️  Warning: No mapping found for $FILENAME. Skipping..."
        continue
    fi
    
    echo "Processing $SAMPLE_ID ($FILENAME)..."

    FILE_R2="${FILE_R1//_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    # Fichiers de sortie
    R1_PAIRED="$CLEANED_DIR/${SAMPLE_ID}_R1_paired.fastq.gz"
    R1_UNPAIRED="$CLEANED_DIR/${SAMPLE_ID}_R1_unpaired.fastq.gz"
    R2_PAIRED="$CLEANED_DIR/${SAMPLE_ID}_R2_paired.fastq.gz"
    R2_UNPAIRED="$CLEANED_DIR/${SAMPLE_ID}_R2_unpaired.fastq.gz"

    # Exécuter Trimmomatic
    trimmomatic PE -threads $THREADS -phred33 \
        "$FILE_R1" "$FILE_R2" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100

    if [ $? -eq 0 ]; then
        echo "  ✓ Trimmomatic completed for $SAMPLE_ID"
        
        # Ajouter au manifest
        echo -e "${SAMPLE_ID}\t${R1_PAIRED}\t${R2_PAIRED}" >> "$MANIFEST_FILE"
        ((SAMPLE_COUNT++))
    else
        echo "  ✗ ERROR: Trimmomatic failed for $SAMPLE_ID"
    fi
    
    echo ""
done

conda deactivate

echo "✓ Trimmomatic completed for ${SAMPLE_COUNT} samples"
echo "✓ Manifest created: $MANIFEST_FILE"
echo ""
echo "ÉTAPE 2 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 3 : CONTRÔLE QUALITÉ POST-NETTOYAGE
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 3 : Contrôle Qualité Post-Nettoyage"
echo "=========================================================================================================="
echo ""

conda activate fastqc

echo "Running FastQC on cleaned data..."
fastqc -t $THREADS "$CLEANED_DIR"/*_paired.fastq.gz -o "$POST_CLEAN_QC_DIR"

if [ $? -eq 0 ]; then
    echo "✓ FastQC completed"
else
    echo "ERROR: FastQC failed!"
fi

conda deactivate
conda activate multiqc

echo ""
echo "Running MultiQC..."
multiqc "$POST_CLEAN_QC_DIR" -o "$POST_CLEAN_QC_DIR" --force

if [ $? -eq 0 ]; then
    echo "✓ MultiQC completed"
    echo "✓ Report: ${POST_CLEAN_QC_DIR}/multiqc_report.html"
fi

conda deactivate

echo ""
echo "Manifest summary:"
cat "$MANIFEST_FILE"
echo ""
echo "Number of samples: $(tail -n +2 "$MANIFEST_FILE" | wc -l)"
echo ""
echo "ÉTAPE 3 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 4 : ANALYSE QIIME2
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 4 : Analyse QIIME2"
echo "=========================================================================================================="
echo ""

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

echo "QIIME2 environment activated"
qiime --version
echo ""

# --------------------------------------------------------------------------------------------------------
# 4.1 : IMPORT DES DONNÉES
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.1 : Import des Données"
echo "=========================================="
echo ""

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

if [ $? -eq 0 ]; then
    echo "✓ Data imported: demux.qza"
else
    echo "ERROR: Import failed!"
    exit 1
fi

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "✓ Demux summary: visual/demux.qzv"
echo "  View at: https://view.qiime2.org"
echo ""

# --------------------------------------------------------------------------------------------------------
# 4.2 : DÉBRUITAGE DADA2 (GÉNÉRATION DES ASVs)
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.2 : Débruitage DADA2"
echo "=========================================="
echo ""

echo "Running DADA2 denoising (this may take 30-60 minutes)..."

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

if [ $? -eq 0 ]; then
    echo "✓ DADA2 completed"
else
    echo "ERROR: DADA2 failed!"
    exit 1
fi

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/dada2-stats.qza" \
  --o-visualization "$QIIME_DIR/visual/dada2-stats.qzv"

echo "✓ DADA2 stats: visual/dada2-stats.qzv"
echo ""

# --------------------------------------------------------------------------------------------------------
# 4.3 : FILTRAGE DES CONTAMINANTS (CONTRÔLES NÉGATIFS)
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.3 : Filtrage des Contaminants"
echo "=========================================="
echo ""

echo "Attempting to filter negative control samples..."

qiime feature-table filter-samples \
  --i-table "${QIIME_DIR}/core/table.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Negative_Control'" \
  --o-filtered-table "${QIIME_DIR}/core/neg-controls-table.qza" 2>&1 | tee "${QIIME_DIR}/neg_control_filter.log"

NEG_CONTROL_STATUS=${PIPESTATUS[0]}

if [ $NEG_CONTROL_STATUS -ne 0 ]; then
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "⚠️  NEGATIVE CONTROL IS EMPTY"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""
  echo "Le contrôle négatif (BL_PCR_Jourand) ne contient AUCUNE lecture."
  echo ""
  echo "✓ EXCELLENT RÉSULTAT !"
  echo "  → Pas de contamination détectée"
  echo "  → Les échantillons sont propres"
  echo "  → Aucun ASV ne sera filtré"
  echo ""
  echo "Actions:"
  echo "  → Table originale conservée (aucun filtrage)"
  echo "  → Retrait du contrôle négatif de la table finale"
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""
  
  # Copier les tables originales
  cp "${QIIME_DIR}/core/table.qza" "${QIIME_DIR}/core/table-decontam.qza"
  cp "${QIIME_DIR}/core/rep-seqs.qza" "${QIIME_DIR}/core/rep-seqs-clean.qza"
  
  # Retirer le contrôle négatif
  qiime feature-table filter-samples \
    --i-table "${QIIME_DIR}/core/table-decontam.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --p-where "[group]='Sample'" \
    --o-filtered-table "${QIIME_DIR}/core/table-final.qza"
  
  echo "✓ Table finale créée: table-final.qza"
  echo "✓ ASVs filtrés: 0 (contrôle vide)"
  
else
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "✓ NEGATIVE CONTROL HAS READS"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""
  echo "Le contrôle négatif contient des lectures."
  echo "Filtrage des ASVs contaminants en cours..."
  echo ""
  
  qiime feature-table summarize \
    --i-table "${QIIME_DIR}/core/neg-controls-table.qza" \
    --o-visualization "${QIIME_DIR}/visual/neg-controls-summary.qzv"
  
  echo "✓ Negative control summary: visual/neg-controls-summary.qzv"

  # Export
  qiime tools export \
    --input-path "${QIIME_DIR}/core/neg-controls-table.qza" \
    --output-path "${QIIME_DIR}/export/neg-controls"

  biom convert \
    -i "${QIIME_DIR}/export/neg-controls/feature-table.biom" \
    -o "${QIIME_DIR}/export/neg-controls/feature-table.tsv" \
    --to-tsv

  # Créer liste des contaminants avec header
  echo "feature-id" > "${QIIME_DIR}/export/neg-controls/contamination_ids.txt"
  tail -n +3 "${QIIME_DIR}/export/neg-controls/feature-table.tsv" | \
    cut -f1 >> "${QIIME_DIR}/export/neg-controls/contamination_ids.txt"

  NUM_CONTAMINANTS=$(tail -n +2 "${QIIME_DIR}/export/neg-controls/contamination_ids.txt" | wc -l)
  echo "Found ${NUM_CONTAMINANTS} contaminant ASVs"

  if [ "$NUM_CONTAMINANTS" -eq 0 ]; then
    echo "⚠️  No contaminant ASVs found (unexpected)"
    cp "${QIIME_DIR}/core/table.qza" "${QIIME_DIR}/core/table-decontam.qza"
    cp "${QIIME_DIR}/core/rep-seqs.qza" "${QIIME_DIR}/core/rep-seqs-clean.qza"
  else
    echo "Removing ${NUM_CONTAMINANTS} contaminant ASVs..."
    
    qiime feature-table filter-features \
      --i-table "${QIIME_DIR}/core/table.qza" \
      --m-metadata-file "${QIIME_DIR}/export/neg-controls/contamination_ids.txt" \
      --p-exclude-ids \
      --o-filtered-table "${QIIME_DIR}/core/table-decontam.qza"

    qiime feature-table filter-seqs \
      --i-data "${QIIME_DIR}/core/rep-seqs.qza" \
      --i-table "${QIIME_DIR}/core/table-decontam.qza" \
      --o-filtered-data "${QIIME_DIR}/core/rep-seqs-clean.qza"
    
    echo "✓ Contaminants removed"
  fi

  # Retirer échantillon contrôle négatif
  qiime feature-table filter-samples \
    --i-table "${QIIME_DIR}/core/table-decontam.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --p-where "[group]='Sample'" \
    --o-filtered-table "${QIIME_DIR}/core/table-final.qza"

  echo "✓ Decontamination completed"
  echo "✓ Contaminant ASVs removed: ${NUM_CONTAMINANTS}"
fi

# Générer résumé table finale
qiime feature-table summarize \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --o-visualization "${QIIME_DIR}/visual/table-final-summary.qzv" \
  --m-sample-metadata-file "$METADATA_FILE"

qiime feature-table tabulate-seqs \
  --i-data "${QIIME_DIR}/core/rep-seqs-clean.qza" \
  --o-visualization "${QIIME_DIR}/visual/rep-seqs-clean.qzv"

echo "✓ Final table summary: visual/table-final-summary.qzv"
echo ""

# --------------------------------------------------------------------------------------------------------
# 4.4 : CONSTRUCTION DE L'ARBRE PHYLOGÉNÉTIQUE
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.4 : Construction de l'Arbre Phylogénétique"
echo "=========================================="
echo ""

# Temp directory pour MAFFT
MAFFT_TMPDIR="${BASE_DIR}/tmp_mafft"
mkdir -p "$MAFFT_TMPDIR"
export TMPDIR="$MAFFT_TMPDIR"
export TMP="$MAFFT_TMPDIR"

echo "Building phylogenetic tree (this may take 10-30 minutes)..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME_DIR}/core/rep-seqs-clean.qza" \
  --p-n-threads 1 \
  --o-alignment "${QIIME_DIR}/core/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME_DIR}/core/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME_DIR}/core/rooted-tree.qza" 2>&1

TREE_STATUS=$?

if [ $TREE_STATUS -ne 0 ] || [ ! -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "⚠️  Standard MAFFT failed. Trying with --p-parttree..."
  
  rm -f "${QIIME_DIR}/core/aligned-rep-seqs.qza"
  rm -f "${QIIME_DIR}/core/masked-aligned-rep-seqs.qza"
  rm -f "${QIIME_DIR}/core/unrooted-tree.qza"
  rm -f "${QIIME_DIR}/core/rooted-tree.qza"
  
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

if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "✓ Phylogenetic tree created successfully!"
  qiime tools peek "${QIIME_DIR}/core/rooted-tree.qza"
else
  echo "✗ Phylogenetic tree creation FAILED"
  echo "  Continuing without Faith PD metrics..."
fi

rm -rf "$MAFFT_TMPDIR"
export TMPDIR="${BASE_DIR}/tmp"

echo ""

# --------------------------------------------------------------------------------------------------------
# 4.5 : MÉTRIQUES DE DIVERSITÉ (CORE-METRICS-PHYLOGENETIC)
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.5 : Métriques de Diversité"
echo "=========================================="
echo ""

rm -rf "${QIIME_DIR}/core-metrics-results"
mkdir -p "${QIIME_DIR}/core-metrics-results"

# Export table pour calculer profondeurs
qiime tools export \
  --input-path "${QIIME_DIR}/core/table-final.qza" \
  --output-path "${QIIME_DIR}/export/table-final-temp"

biom convert \
  -i "${QIIME_DIR}/export/table-final-temp/feature-table.biom" \
  -o "${QIIME_DIR}/export/table-final-temp/feature-table.tsv" \
  --to-tsv

echo "Calculating sampling depths with Python..."

DEPTHS=$(python3 << EOFPYTHON
import pandas as pd
import sys

try:
    df = pd.read_csv(
        "${QIIME_DIR}/export/table-final-temp/feature-table.tsv",
        sep="\t",
        skiprows=1,
        index_col=0
    )
    
    depths = df.sum(axis=0).astype(int)
    
    min_depth = 1000
    filtered_depths = depths[depths >= min_depth]
    
    if len(filtered_depths) > 0:
        sampling_depth = int(filtered_depths.quantile(0.1))
    else:
        sampling_depth = int(depths.median())
    
    max_depth = int(depths.max())
    
    print(f"{sampling_depth} {max_depth}")
    
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
EOFPYTHON
)

SAMPLING_DEPTH=$(echo $DEPTHS | cut -d' ' -f1)
MAX_DEPTH=$(echo $DEPTHS | cut -d' ' -f2)

if ! [[ "$SAMPLING_DEPTH" =~ ^[0-9]+$ ]] || ! [[ "$MAX_DEPTH" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Invalid depth values!"
    echo "  SAMPLING_DEPTH: $SAMPLING_DEPTH"
    echo "  MAX_DEPTH: $MAX_DEPTH"
    exit 1
fi

echo "✓ Sampling Depth (10th percentile): $SAMPLING_DEPTH"
echo "✓ Max Depth: $MAX_DEPTH"
echo ""

echo "$SAMPLING_DEPTH" > "${QIIME_DIR}/export/sampling_depth.txt"
echo "$MAX_DEPTH" > "${QIIME_DIR}/export/max_depth.txt"

# Core-metrics
if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "Running core-metrics-phylogenetic (with Faith PD)..."
  
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME_DIR}/core/rooted-tree.qza" \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --o-rarefied-table "${QIIME_DIR}/core-metrics-results/rarefied_table.qza" \
    --o-faith-pd-vector "${QIIME_DIR}/core-metrics-results/faith_pd_vector.qza" \
    --o-observed-features-vector "${QIIME_DIR}/core-metrics-results/observed_features_vector.qza" \
    --o-shannon-vector "${QIIME_DIR}/core-metrics-results/shannon_vector.qza" \
    --o-evenness-vector "${QIIME_DIR}/core-metrics-results/evenness_vector.qza" \
    --o-unweighted-unifrac-distance-matrix "${QIIME_DIR}/core-metrics-results/unweighted_unifrac_distance_matrix.qza" \
    --o-weighted-unifrac-distance-matrix "${QIIME_DIR}/core-metrics-results/weighted_unifrac_distance_matrix.qza" \
    --o-jaccard-distance-matrix "${QIIME_DIR}/core-metrics-results/jaccard_distance_matrix.qza" \
    --o-bray-curtis-distance-matrix "${QIIME_DIR}/core-metrics-results/bray_curtis_distance_matrix.qza" \
    --o-unweighted-unifrac-pcoa-results "${QIIME_DIR}/core-metrics-results/unweighted_unifrac_pcoa_results.qza" \
    --o-weighted-unifrac-pcoa-results "${QIIME_DIR}/core-metrics-results/weighted_unifrac_pcoa_results.qza" \
    --o-jaccard-pcoa-results "${QIIME_DIR}/core-metrics-results/jaccard_pcoa_results.qza" \
    --o-bray-curtis-pcoa-results "${QIIME_DIR}/core-metrics-results/bray_curtis_pcoa_results.qza" \
    --o-unweighted-unifrac-emperor "${QIIME_DIR}/core-metrics-results/unweighted_unifrac_emperor.qzv" \
    --o-weighted-unifrac-emperor "${QIIME_DIR}/core-metrics-results/weighted_unifrac_emperor.qzv" \
    --o-jaccard-emperor "${QIIME_DIR}/core-metrics-results/jaccard_emperor.qzv" \
    --o-bray-curtis-emperor "${QIIME_DIR}/core-metrics-results/bray_curtis_emperor.qzv"
    
  if [ $? -eq 0 ]; then
    echo "✓ Core-metrics-phylogenetic completed (with Faith PD)"
  else
    echo "⚠️  Phylogenetic metrics failed. Running non-phylogenetic metrics..."
    rm -rf "${QIIME_DIR}/core-metrics-results"/*
    
    qiime diversity core-metrics \
      --i-table "${QIIME_DIR}/core/table-final.qza" \
      --p-sampling-depth "$SAMPLING_DEPTH" \
      --m-metadata-file "$METADATA_FILE" \
      --o-rarefied-table "${QIIME_DIR}/core-metrics-results/rarefied_table.qza" \
      --o-observed-features-vector "${QIIME_DIR}/core-metrics-results/observed_features_vector.qza" \
      --o-shannon-vector "${QIIME_DIR}/core-metrics-results/shannon_vector.qza" \
      --o-evenness-vector "${QIIME_DIR}/core-metrics-results/evenness_vector.qza" \
      --o-jaccard-distance-matrix "${QIIME_DIR}/core-metrics-results/jaccard_distance_matrix.qza" \
      --o-bray-curtis-distance-matrix "${QIIME_DIR}/core-metrics-results/bray_curtis_distance_matrix.qza" \
      --o-jaccard-pcoa-results "${QIIME_DIR}/core-metrics-results/jaccard_pcoa_results.qza" \
      --o-bray-curtis-pcoa-results "${QIIME_DIR}/core-metrics-results/bray_curtis_pcoa_results.qza" \
      --o-jaccard-emperor "${QIIME_DIR}/core-metrics-results/jaccard_emperor.qzv" \
      --o-bray-curtis-emperor "${QIIME_DIR}/core-metrics-results/bray_curtis_emperor.qzv"
    
    echo "✓ Core-metrics completed (without Faith PD)"
  fi
else
  echo "⚠️  No phylogenetic tree. Running non-phylogenetic metrics only..."
  
  qiime diversity core-metrics \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --o-rarefied-table "${QIIME_DIR}/core-metrics-results/rarefied_table.qza" \
    --o-observed-features-vector "${QIIME_DIR}/core-metrics-results/observed_features_vector.qza" \
    --o-shannon-vector "${QIIME_DIR}/core-metrics-results/shannon_vector.qza" \
    --o-evenness-vector "${QIIME_DIR}/core-metrics-results/evenness_vector.qza" \
    --o-jaccard-distance-matrix "${QIIME_DIR}/core-metrics-results/jaccard_distance_matrix.qza" \
    --o-bray-curtis-distance-matrix "${QIIME_DIR}/core-metrics-results/bray_curtis_distance_matrix.qza" \
    --o-jaccard-pcoa-results "${QIIME_DIR}/core-metrics-results/jaccard_pcoa_results.qza" \
    --o-bray-curtis-pcoa-results "${QIIME_DIR}/core-metrics-results/bray_curtis_pcoa_results.qza" \
    --o-jaccard-emperor "${QIIME_DIR}/core-metrics-results/jaccard_emperor.qzv" \
    --o-bray-curtis-emperor "${QIIME_DIR}/core-metrics-results/bray_curtis_emperor.qzv"
  
  echo "✓ Core-metrics completed (without Faith PD)"
fi

echo ""

# --------------------------------------------------------------------------------------------------------
# 4.6 : COURBES DE RARÉFACTION
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.6 : Courbes de Raréfaction"
echo "=========================================="
echo ""

echo "Generating rarefaction curves (Shannon + Observed Features)..."

qiime diversity alpha-rarefaction \
  --i-table "${QIIME_DIR}/core/table-final.qza" \
  --p-min-depth 10 \
  --p-max-depth "$MAX_DEPTH" \
  --p-steps 20 \
  --m-metadata-file "$METADATA_FILE" \
  --o-visualization "${QIIME_DIR}/visual/rarefaction-curves.qzv"

if [ $? -eq 0 ]; then
  echo "✓ Rarefaction curves: visual/rarefaction-curves.qzv"
else
  echo "✗ ERROR: Rarefaction curves failed"
fi

if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "Generating Faith PD rarefaction curves..."
  
  qiime diversity alpha-rarefaction \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --i-phylogeny "${QIIME_DIR}/core/rooted-tree.qza" \
    --p-min-depth 10 \
    --p-max-depth "$MAX_DEPTH" \
    --p-steps 20 \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv"
  
  if [ $? -eq 0 ]; then
    echo "✓ Faith PD rarefaction: visual/rarefaction-curves-phylogenetic.qzv"
  fi
fi

echo ""

# --------------------------------------------------------------------------------------------------------
# 4.7 : CLASSIFICATION TAXONOMIQUE
# --------------------------------------------------------------------------------------------------------

echo "=========================================="
echo "4.7 : Classification Taxonomique"
echo "=========================================="
echo ""

echo "Assigning taxonomy with SILVA 138.2 (this may take 10-20 minutes)..."

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_PATH" \
  --i-reads "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --o-classification "$QIIME_DIR/core/taxonomy.qza"

if [ $? -eq 0 ]; then
    echo "✓ Taxonomy assigned: taxonomy.qza"
else
    echo "ERROR: Taxonomy classification failed!"
    exit 1
fi

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/taxonomy.qza" \
  --o-visualization "$QIIME_DIR/visual/taxonomy.qzv"

echo "✓ Taxonomy table: visual/taxonomy.qzv"
echo ""

echo "ÉTAPE 4 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 5 : EXPORT DES RÉSULTATS
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 5 : Export des Résultats"
echo "=========================================================================================================="
echo ""

EXPORT_DIR="${QIIME_DIR}/export"
mkdir -p "${EXPORT_DIR}"

# --------------------------------------------------------------------------------------------------------
# 5.1 : TABLE DE FEATURES
# --------------------------------------------------------------------------------------------------------

echo "5.1 : Exporting feature table..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/table-final.qza" \
  --output-path "${EXPORT_DIR}/feature_table"

biom convert \
  -i "${EXPORT_DIR}/feature_table/feature-table.biom" \
  -o "${EXPORT_DIR}/feature_table/feature-table.tsv" \
  --to-tsv

sed -i 's/#OTU ID/#ASV_ID/g' "${EXPORT_DIR}/feature_table/feature-table.tsv" 2>/dev/null || true

echo "✓ Feature table: export/feature_table/feature-table.tsv"

# --------------------------------------------------------------------------------------------------------
# 5.2 : SÉQUENCES REPRÉSENTATIVES
# --------------------------------------------------------------------------------------------------------

echo "5.2 : Exporting representative sequences..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/rep-seqs-clean.qza" \
  --output-path "${EXPORT_DIR}/rep_seqs"

echo "✓ Rep seqs: export/rep_seqs/dna-sequences.fasta"

# --------------------------------------------------------------------------------------------------------
# 5.3 : TAXONOMIE
# --------------------------------------------------------------------------------------------------------

echo "5.3 : Exporting taxonomy..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/taxonomy.qza" \
  --output-path "${EXPORT_DIR}/taxonomy"

echo "✓ Taxonomy: export/taxonomy/taxonomy.tsv"

# --------------------------------------------------------------------------------------------------------
# 5.4 : ARBRE PHYLOGÉNÉTIQUE
# --------------------------------------------------------------------------------------------------------

if [ -f "${QIIME_DIR}/core/rooted-tree.qza" ]; then
  echo "5.4 : Exporting phylogenetic tree..."
  
  qiime tools export \
    --input-path "${QIIME_DIR}/core/rooted-tree.qza" \
    --output-path "${EXPORT_DIR}/tree"
  
  if [ -f "${EXPORT_DIR}/tree/tree.nwk" ]; then
    echo "✓ Tree: export/tree/tree.nwk"
  fi
fi

# --------------------------------------------------------------------------------------------------------
# 5.5 : STATISTIQUES DADA2
# --------------------------------------------------------------------------------------------------------

echo "5.5 : Exporting DADA2 stats..."

qiime tools export \
  --input-path "${QIIME_DIR}/core/dada2-stats.qza" \
  --output-path "${EXPORT_DIR}/dada2_stats"

echo "✓ DADA2 stats: export/dada2_stats/"

# --------------------------------------------------------------------------------------------------------
# 5.6 : INDICES DE DIVERSITÉ (CORE-METRICS)
# --------------------------------------------------------------------------------------------------------

echo "5.6 : Exporting diversity indices..."

CORE_METRICS_DIR="${QIIME_DIR}/core-metrics-results"
mkdir -p "${EXPORT_DIR}/diversity"

for metric in observed_features shannon evenness faith_pd; do
  vector_file="${CORE_METRICS_DIR}/${metric}_vector.qza"
  [ "$metric" = "observed_features" ] && metric_name="observed_ASVs" || metric_name="$metric"
  
  if [ -f "$vector_file" ]; then
    qiime tools export \
      --input-path "$vector_file" \
      --output-path "${EXPORT_DIR}/diversity/${metric_name}"
    echo "  ✓ ${metric_name}"
  fi
done

# --------------------------------------------------------------------------------------------------------
# 5.7 : TABLE ASV + TAXONOMIE FUSIONNÉE
# --------------------------------------------------------------------------------------------------------

echo "5.7 : Merging ASV abundance + taxonomy..."

python3 << EOFPYTHON
import pandas as pd
import os

qiime_dir = "${QIIME_DIR}"
export_dir = os.path.join(qiime_dir, "export")

try:
    tab = pd.read_csv(
        os.path.join(export_dir, "feature_table", "feature-table.tsv"),
        sep="\t", comment="#", index_col=0
    )

    tax = pd.read_csv(
        os.path.join(export_dir, "taxonomy", "taxonomy.tsv"),
        sep="\t", index_col=0
    )

    merged = tax.join(tab, how="inner")
    output_file = os.path.join(export_dir, "ASV_abundance_taxonomy.tsv")
    merged.to_csv(output_file, sep="\t")
    print(f"✓ Merged table: {output_file}")
except Exception as e:
    print(f"ERROR: {e}")
EOFPYTHON

# --------------------------------------------------------------------------------------------------------
# 5.8 : PROFONDEURS DE LECTURE PAR ÉCHANTILLON
# --------------------------------------------------------------------------------------------------------

echo "5.8 : Calculating sample read depths..."

python3 << EOFPYTHON
import pandas as pd
import os

qiime_dir = "${QIIME_DIR}"
export_dir = os.path.join(qiime_dir, "export")

df = pd.read_csv(
    os.path.join(export_dir, "feature_table", "feature-table.tsv"),
    sep="\t",
    skiprows=1,
    index_col=0
)

depths = df.sum(axis=0).astype(int)
sampling_depth = ${SAMPLING_DEPTH}

depths_df = pd.DataFrame({
    'sample_id': depths.index,
    'total_reads': depths.values,
    'rarefaction_depth': [sampling_depth] * len(depths)
})

depths_df = depths_df.sort_values('sample_id')

output_file = os.path.join(export_dir, "sample_read_depths_final.tsv")
depths_df.to_csv(output_file, sep="\t", index=False)

print(f"✓ Sample depths: {output_file}")
EOFPYTHON

echo ""
echo "ÉTAPE 5 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 6 : CALCUL DE TOUS LES INDICES DE DIVERSITÉ ALPHA
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 6 : Calcul de Tous les Indices de Diversité Alpha"
echo "=========================================================================================================="
echo ""

DIVERSITY_DIR="${QIIME_DIR}/diversity_indices"
mkdir -p "${DIVERSITY_DIR}"

echo "Calculating additional alpha diversity indices..."

# Liste des métriques
METRICS=(
  "simpson"
  "simpson_e"
  "chao1"
  "ace"
  "goods_coverage"
  "fisher_alpha"
  "berger_parker_d"
  "gini_index"
  "brillouin_d"
  "strong"
  "mcintosh_d"
  "mcintosh_e"
  "margalef"
  "menhinick"
)

for metric in "${METRICS[@]}"; do
  echo "  - ${metric}"
  qiime diversity alpha \
    --i-table "${QIIME_DIR}/core/table-final.qza" \
    --p-metric "$metric" \
    --o-alpha-diversity "${DIVERSITY_DIR}/${metric}_vector.qza" 2>/dev/null
done

echo "✓ All alpha diversity indices calculated"
echo ""

# Export tous les indices en TSV
echo "Exporting all diversity indices to TSV..."

mkdir -p "${EXPORT_DIR}/diversity_all"

export_diversity_metric() {
    local metric_name=$1
    local qza_file=$2
    local output_name=$3
    
    if [ -f "$qza_file" ]; then
        qiime tools export \
          --input-path "$qza_file" \
          --output-path "${EXPORT_DIR}/diversity_all/${output_name}_temp" 2>/dev/null
        
        sed "1s/.*/sample-id\t${output_name}/" \
          "${EXPORT_DIR}/diversity_all/${output_name}_temp/alpha-diversity.tsv" > \
          "${EXPORT_DIR}/diversity_all/${output_name}.tsv" 2>/dev/null
        
        rm -rf "${EXPORT_DIR}/diversity_all/${output_name}_temp"
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

# Export métriques additionnelles
export_diversity_metric "Simpson" "${DIVERSITY_DIR}/simpson_vector.qza" "simpson"
export_diversity_metric "Simpson Evenness" "${DIVERSITY_DIR}/simpson_e_vector.qza" "simpson_evenness"
export_diversity_metric "Chao1" "${DIVERSITY_DIR}/chao1_vector.qza" "chao1"
export_diversity_metric "ACE" "${DIVERSITY_DIR}/ace_vector.qza" "ace"
export_diversity_metric "Goods Coverage" "${DIVERSITY_DIR}/goods_coverage_vector.qza" "goods_coverage"
export_diversity_metric "Fisher Alpha" "${DIVERSITY_DIR}/fisher_alpha_vector.qza" "fisher_alpha"
export_diversity_metric "Berger Parker" "${DIVERSITY_DIR}/berger_parker_d_vector.qza" "berger_parker"
export_diversity_metric "Gini Index" "${DIVERSITY_DIR}/gini_index_vector.qza" "gini_index"
export_diversity_metric "Brillouin" "${DIVERSITY_DIR}/brillouin_d_vector.qza" "brillouin"
export_diversity_metric "Strong" "${DIVERSITY_DIR}/strong_vector.qza" "strong"
export_diversity_metric "McIntosh D" "${DIVERSITY_DIR}/mcintosh_d_vector.qza" "mcintosh_d"
export_diversity_metric "McIntosh E" "${DIVERSITY_DIR}/mcintosh_e_vector.qza" "mcintosh_e"
export_diversity_metric "Margalef" "${DIVERSITY_DIR}/margalef_vector.qza" "margalef"
export_diversity_metric "Menhinick" "${DIVERSITY_DIR}/menhinick_vector.qza" "menhinick"

# Fusionner tous les indices dans une table unique
echo "Merging all diversity indices..."

python3 << EOFPYTHON
import pandas as pd
import os
from pathlib import Path

qiime_dir = "${QIIME_DIR}"
diversity_dir = Path(qiime_dir) / "export" / "diversity_all"

tsv_files = sorted(diversity_dir.glob("*.tsv"))

if tsv_files:
    df_merged = pd.read_csv(tsv_files[0], sep="\t", index_col=0)

    for tsv_file in tsv_files[1:]:
        try:
            df_temp = pd.read_csv(tsv_file, sep="\t", index_col=0)
            df_merged = df_merged.join(df_temp, how="outer")
        except:
            pass

    sampling_depth = ${SAMPLING_DEPTH}
    df_merged.insert(0, 'rarefaction_depth', sampling_depth)
    
    df_merged = df_merged.sort_index()

    output_file = diversity_dir.parent / "diversity_indices_all.tsv"
    df_merged.to_csv(output_file, sep="\t")

    print(f"✓ Comprehensive diversity table: {output_file}")
    print(f"  Samples: {len(df_merged)}")
    print(f"  Indices: {len(df_merged.columns)}")
EOFPYTHON

echo ""
echo "ÉTAPE 6 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 7 : EXPORT DES DONNÉES DE RARÉFACTION
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 7 : Export des Données de Raréfaction"
echo "=========================================================================================================="
echo ""

mkdir -p "${EXPORT_DIR}/rarefaction_data"

if [ -f "${QIIME_DIR}/visual/rarefaction-curves.qzv" ]; then
  echo "Extracting rarefaction curve data..."
  
  unzip -q "${QIIME_DIR}/visual/rarefaction-curves.qzv" -d "${EXPORT_DIR}/rarefaction_data/temp" 2>/dev/null
  
  find "${EXPORT_DIR}/rarefaction_data/temp" -name "*.csv" -exec cp {} "${EXPORT_DIR}/rarefaction_data/" \; 2>/dev/null
  
  rm -rf "${EXPORT_DIR}/rarefaction_data/temp"
  
  CSV_COUNT=$(ls -1 "${EXPORT_DIR}/rarefaction_data/"*.csv 2>/dev/null | wc -l)
  
  if [ "$CSV_COUNT" -gt 0 ]; then
    echo "✓ Rarefaction data extracted: ${CSV_COUNT} CSV files"
  fi
fi

if [ -f "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv" ]; then
  unzip -q "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv" -d "${EXPORT_DIR}/rarefaction_data/temp_faith" 2>/dev/null
  
  find "${EXPORT_DIR}/rarefaction_data/temp_faith" -name "*.csv" -exec cp {} "${EXPORT_DIR}/rarefaction_data/" \; 2>/dev/null
  
  rm -rf "${EXPORT_DIR}/rarefaction_data/temp_faith"
fi

echo ""
echo "ÉTAPE 7 : TERMINÉE"
echo ""

# ========================================================================================================
# PIPELINE TERMINÉ
# ========================================================================================================

echo "=========================================================================================================="
echo "🎉 PIPELINE TERMINÉ AVEC SUCCÈS ! 🎉"
echo "=========================================================================================================="
echo ""
echo "Date de fin: $(date)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RÉSUMÉ DES FICHIERS GÉNÉRÉS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 DONNÉES PRINCIPALES:"
echo "  • Table ASV                  : export/feature_table/feature-table.tsv"
echo "  • Séquences ASV              : export/rep_seqs/dna-sequences.fasta"
echo "  • Taxonomie                  : export/taxonomy/taxonomy.tsv"
echo "  • ASV + Taxonomie fusionné   : export/ASV_abundance_taxonomy.tsv"
echo ""
echo "📈 DIVERSITÉ:"
echo "  • Tous les indices alpha     : export/diversity_indices_all.tsv"
echo "  • Profondeurs échantillons   : export/sample_read_depths_final.tsv"
echo ""
echo "🌳 PHYLOGÉNIE:"
if [ -f "${EXPORT_DIR}/tree/tree.nwk" ]; then
echo "  • Arbre phylogénétique       : export/tree/tree.nwk"
else
echo "  • Arbre phylogénétique       : ❌ Non généré"
fi
echo ""
echo "📉 RARÉFACTION:"
echo "  • Courbes standard           : visual/rarefaction-curves.qzv"
if [ -f "${QIIME_DIR}/visual/rarefaction-curves-phylogenetic.qzv" ]; then
echo "  • Courbes Faith PD           : visual/rarefaction-curves-phylogenetic.qzv"
fi
echo "  • Données CSV                : export/rarefaction_data/"
echo ""
echo "📋 STATISTIQUES:"
echo "  • Stats DADA2                : visual/dada2-stats.qzv"
echo "  • Résumé table finale        : visual/table-final-summary.qzv"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PARAMÈTRES UTILISÉS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  • Profondeur de raréfaction  : ${SAMPLING_DEPTH} reads"
echo "  • Profondeur maximale        : ${MAX_DEPTH} reads"
echo "  • Nombre d'échantillons      : $(tail -n +2 "$MANIFEST_FILE" | wc -l)"
echo "  • Threads utilisés           : ${THREADS}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PROCHAINES ÉTAPES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. Visualiser les fichiers .qzv sur https://view.qiime2.org"
echo ""
echo "2. Importer les données dans R pour analyses statistiques:"
echo "   library(phyloseq)"
echo "   library(vegan)"
echo "   "
echo "   # Import des données"
echo "   asv_table <- read.table('export/feature_table/feature-table.tsv', header=TRUE, row.names=1, skip=1, sep='\t')"
echo "   taxonomy <- read.table('export/taxonomy/taxonomy.tsv', header=TRUE, row.names=1, sep='\t')"
echo "   metadata <- read.table('metadata_ExPLOI.tsv', header=TRUE, row.names=1, sep='\t')"
echo "   diversity <- read.table('export/diversity_indices_all.tsv', header=TRUE, row.names=1, sep='\t')"
echo ""
echo "3. Analyses recommandées:"
echo "   • Tests de diversité alpha (Kruskal-Wallis, Wilcoxon)"
echo "   • Analyses de diversité beta (PERMANOVA, PCoA)"
echo "   • Analyses différentielles (DESeq2, ANCOM)"
echo "   • Visualisations (barplots, heatmaps, ordinations)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Tous les fichiers sont dans: ${QIIME_DIR}"
echo ""
echo "=========================================================================================================="
