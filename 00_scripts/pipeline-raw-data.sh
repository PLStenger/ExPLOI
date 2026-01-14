#!/bin/bash

# ========================================================================================================
# PIPELINE QIIME2 - Données BRUTES (sans Trimmomatic)
# ========================================================================================================
# Cette version utilise les fichiers RAW directement, sans passer par Trimmomatic
# DADA2 fera le filtrage qualité lui-même

BASE_DIR="/nvme/bio/data_fungi/ExPLOI"
RAW_DATA_DIR="${BASE_DIR}/01_raw_data"
QIIME_DIR="${BASE_DIR}/05_QIIME2_RAW"  # Nouveau dossier pour ne pas écraser
METADATA_FILE="${BASE_DIR}/metadata_ExPLOI.tsv"
MANIFEST_FILE="${BASE_DIR}/manifest_ExPLOI_RAW.txt"
CLASSIFIER_PATH="/nvme/bio/data_fungi/BioIndic_La_Reunion_Island_seawater_four_month_SED/05_QIIME2/Original_reads_16S/taxonomy/16S/Classifier.qza"
THREADS=16

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

export PYTHONWARNINGS="ignore"
export QIIME_DIR="${QIIME_DIR}"

echo "=========================================================================================================="
echo "QIIME2 - Pipeline avec DONNÉES BRUTES (sans Trimmomatic)"
echo "=========================================================================================================="
echo ""
echo "STRATÉGIE:"
echo "  → Utiliser les fichiers RAW de 300bp"
echo "  → DADA2 fera le filtrage qualité"
echo "  → Troncature optimale pour overlap de 80bp"
echo ""
echo "=========================================================================================================="
echo ""

mkdir -p "$QIIME_DIR/core" "$QIIME_DIR/visual" "$QIIME_DIR/export"

# ========================================================================================================
# CRÉATION DU MANIFEST AVEC FICHIERS RAW
# ========================================================================================================

echo "Création du manifest avec fichiers RAW..."
echo ""

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

SAMPLE_COUNT=0

for FILE_R1 in "$RAW_DATA_DIR"/*_R1_001.fastq.gz; do
    FILENAME=$(basename "$FILE_R1")
    BASE_NAME=${FILENAME%_L001_R1_001.fastq.gz}
    
    SAMPLE_ID=""
    for KEY in "${!SAMPLES[@]}"; do
        if [[ "$BASE_NAME" == *"$KEY"* ]]; then
            SAMPLE_ID="${SAMPLES[$KEY]}"
            break
        fi
    done

    if [ -z "$SAMPLE_ID" ]; then
        continue
    fi
    
    FILE_R2="${FILE_R1//_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    echo -e "${SAMPLE_ID}\t${FILE_R1}\t${FILE_R2}" >> "$MANIFEST_FILE"
    ((SAMPLE_COUNT++))
done

echo "✓ Manifest créé: $SAMPLE_COUNT échantillons"
echo ""

# ========================================================================================================
# IMPORT QIIME2
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "IMPORT DES DONNÉES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

if [ $? -eq 0 ]; then
    echo "✓ Données importées"
else
    echo "❌ ERROR: Import failed!"
    exit 1
fi

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "✓ Demux summary: visual/demux.qzv"
echo ""

# ========================================================================================================
# DADA2 AVEC PARAMÈTRES OPTIMISÉS POUR READS 2x300bp
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DADA2 DÉBRUITAGE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Paramètres pour reads 2x300bp, amplicon V3-V4 (460bp):"
echo "  --p-trim-left-f 17     (enlever primer forward 341F)"
echo "  --p-trim-left-r 21     (enlever primer reverse 805R)"
echo "  --p-trunc-len-f 270    (garder 270bp après trim)"
echo "  --p-trunc-len-r 210    (garder 210bp après trim)"
echo "  → Overlap théorique: 270 + 210 - 460 = 20bp minimum"
echo "  → Avec adaptateurs enlevés: overlap ~80bp ✓"
echo ""
echo "Filtres qualité permissifs:"
echo "  --p-max-ee-f 3"
echo "  --p-max-ee-r 3"
echo "  --p-min-overlap 12"
echo ""
echo "Lancement DADA2 (30-60 minutes)..."
echo "Début: $(date)"
echo ""

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$QIIME_DIR/core/demux.qza" \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 270 \
  --p-trunc-len-r 210 \
  --p-max-ee-f 3 \
  --p-max-ee-r 3 \
  --p-trunc-q 2 \
  --p-min-overlap 12 \
  --p-n-threads "$THREADS" \
  --p-chimera-method consensus \
  --o-table "$QIIME_DIR/core/table.qza" \
  --o-representative-sequences "$QIIME_DIR/core/rep-seqs.qza" \
  --o-denoising-stats "$QIIME_DIR/core/dada2-stats.qza" \
  --verbose

DADA2_STATUS=$?

echo ""
echo "Fin: $(date)"
echo ""

if [ $DADA2_STATUS -ne 0 ]; then
    echo "❌ ERROR: DADA2 failed!"
    exit 1
fi

echo "✓ DADA2 terminé!"
echo ""

# Générer visualisations
qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/dada2-stats.qza" \
  --o-visualization "$QIIME_DIR/visual/dada2-stats.qzv"

qiime feature-table summarize \
  --i-table "$QIIME_DIR/core/table.qza" \
  --o-visualization "$QIIME_DIR/visual/table-summary.qzv" \
  --m-sample-metadata-file "$METADATA_FILE"

echo "✓ Visualisations générées"
echo ""

# ========================================================================================================
# ANALYSE DES RÉSULTATS
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RÉSULTATS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Export stats
qiime tools export \
    --input-path "$QIIME_DIR/core/dada2-stats.qza" \
    --output-path "$QIIME_DIR/export/dada2_stats_raw"

python3 << 'EOFPYTHON'
import pandas as pd
import os

stats_file = "/nvme/bio/data_fungi/ExPLOI/05_QIIME2_RAW/export/dada2_stats_raw/stats.tsv"

if os.path.exists(stats_file):
    df = pd.read_csv(stats_file, sep="\t", comment="#")
    
    print("STATISTIQUES PAR ÉCHANTILLON:")
    print("=" * 90)
    
    for idx, row in df.iterrows():
        sample = row['sample-id']
        input_reads = int(row['input'])
        filtered_reads = int(row['filtered'])
        merged = int(row['merged'])
        final = int(row['non-chimeric'])
        
        filter_pct = (filtered_reads / input_reads * 100) if input_reads > 0 else 0
        merge_pct = (merged / filtered_reads * 100) if filtered_reads > 0 else 0
        final_pct = (final / input_reads * 100) if input_reads > 0 else 0
        
        status = "✅" if merge_pct > 70 else "⚠️" if merge_pct > 50 else "❌"
        
        print(f"{status} {sample:15s} | Input: {input_reads:>6,} | Filt: {filtered_reads:>6,} ({filter_pct:>5.1f}%) | Merged: {merged:>6,} ({merge_pct:>5.1f}%) | Final: {final:>6,} ({final_pct:>5.1f}%)")
    
    total_input = df['input'].astype(int).sum()
    total_filtered = df['filtered'].astype(int).sum()
    total_merged = df['merged'].astype(int).sum()
    total_final = df['non-chimeric'].astype(int).sum()
    
    print("\n" + "=" * 90)
    print("TOTAUX:")
    print(f"  Input      : {total_input:>10,} reads")
    print(f"  Filtered   : {total_filtered:>10,} reads ({total_filtered/total_input*100:.1f}%)")
    print(f"  Merged     : {total_merged:>10,} reads ({total_merged/total_filtered*100:.1f}%)")
    print(f"  Final      : {total_final:>10,} reads ({total_final/total_input*100:.1f}%)")
    print("=" * 90)

EOFPYTHON

# Compter ASVs
TEMP_EXPORT="$QIIME_DIR/export/temp_count_$$"
mkdir -p "$TEMP_EXPORT"

qiime tools export \
  --input-path "$QIIME_DIR/core/table.qza" \
  --output-path "$TEMP_EXPORT" 2>/dev/null

NUM_ASVS=$(biom summarize-table -i "$TEMP_EXPORT/feature-table.biom" 2>/dev/null | grep "Num observations:" | awk '{print $3}')

rm -rf "$TEMP_EXPORT"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "NOMBRE D'ASVs DÉTECTÉS: ${NUM_ASVS}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -n "$NUM_ASVS" ]; then
    RATIO=$(awk "BEGIN {printf \"%.1f\", ($NUM_ASVS / 4700) * 100}")
    
    echo "Comparaison:"
    echo "  Avec Trimmomatic  : 88 ASVs (1.9%)"
    echo "  Sans Trimmomatic  : ${NUM_ASVS} ASVs (${RATIO}%)"
    echo "  Société séquençage: 4700 ASVs (100%)"
    echo ""
    
    if [ "$NUM_ASVS" -gt 3000 ]; then
        echo "✅✅✅ EXCELLENT! Résultats comparables à la société de séquençage!"
        echo ""
        echo "PROCHAINES ÉTAPES:"
        echo "  1. Continuer avec décontamination (section 4.3)"
        echo "  2. Taxonomie, diversité, etc."
        echo ""
    elif [ "$NUM_ASVS" -gt 1000 ]; then
        echo "✓ Amélioration MAJEURE! (×${NUM_ASVS}/88 = $((NUM_ASVS/88))x plus d'ASVs)"
        echo ""
        echo "Toujours inférieur à la société de séquençage."
        echo "POSSIBILITÉS:"
        echo "  - Ajuster --p-trunc-len pour plus de permissivité"
        echo "  - La société utilise peut-être un autre algorithme"
        echo ""
    else
        echo "⚠️ Amélioration limitée"
        echo ""
        echo "ACTIONS:"
        echo "  1. Vérifier demux.qzv (qualité des reads bruts)"
        echo "  2. Contacter la société pour connaître leur pipeline exact"
        echo ""
    fi
fi

echo "Fichiers à vérifier:"
echo "  - Qualité brute    : ${QIIME_DIR}/visual/demux.qzv"
echo "  - Stats DADA2      : ${QIIME_DIR}/visual/dada2-stats.qzv"
echo "  - Table résumé     : ${QIIME_DIR}/visual/table-summary.qzv"
echo ""
echo "Ouvrir sur: https://view.qiime2.org"
echo ""
echo "=========================================================================================================="
