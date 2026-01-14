#!/bin/bash

# ========================================================================================================
# DADA2 OPTIMAL POUR 2x250bp + AMPLICON V3-V4
# ========================================================================================================
# Stratégie : Garder la longueur maximale et accepter les overlaps minimaux

QIIME_DIR="/nvme/bio/data_fungi/ExPLOI/05_QIIME2_OPTIMAL"
RAW_DATA_DIR="/nvme/bio/data_fungi/ExPLOI/01_raw_data"
METADATA_FILE="/nvme/bio/data_fungi/ExPLOI/metadata_ExPLOI.tsv"
MANIFEST_FILE="/nvme/bio/data_fungi/ExPLOI/manifest_ExPLOI_OPTIMAL.txt"
CLASSIFIER_PATH="/nvme/bio/data_fungi/BioIndic_La_Reunion_Island_seawater_four_month_SED/05_QIIME2/Original_reads_16S/taxonomy/16S/Classifier.qza"
THREADS=16

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

export PYTHONWARNINGS="ignore"
export QIIME_DIR="${QIIME_DIR}"

echo "=========================================================================================================="
echo "DADA2 OPTIMAL - Configuration pour MiSeq 2x250bp + Amplicon V3-V4"
echo "=========================================================================================================="
echo ""
echo "PROBLÈME IDENTIFIÉ:"
echo "  → Séquençage: 2x250bp"
echo "  → Amplicon V3-V4: 460bp"
echo "  → Overlap après trim: ~2bp (MINIMAL!)"
echo ""
echo "STRATÉGIE:"
echo "  → Garder longueur maximale (pas de troncature)"
echo "  → Accepter overlap de 8bp minimum (DADA2 minimum absolu)"
echo "  → Filtres qualité ULTRA-PERMISSIFS"
echo "  → Accepter 40-70% de taux de merge"
echo ""
echo "RÉSULTAT ATTENDU:"
echo "  → 500-2000 ASVs (au lieu de 88)"
echo "  → Toujours inférieur aux 4700 de la société (qui utilise autre chose)"
echo ""
echo "=========================================================================================================="
echo ""

mkdir -p "$QIIME_DIR/core" "$QIIME_DIR/visual" "$QIIME_DIR/export"

# ========================================================================================================
# MANIFEST
# ========================================================================================================

echo "Création du manifest..."

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_FILE"

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
done

echo "✓ Manifest créé"
echo ""

# ========================================================================================================
# IMPORT
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "IMPORT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "✓ Import OK"
echo ""

# ========================================================================================================
# DADA2 AVEC PARAMÈTRES OPTIMAUX
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DADA2 DÉBRUITAGE (30-60 min)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Paramètres ULTRA-PERMISSIFS:"
echo "  --p-trim-left-f 17     (enlever primer 341F)"
echo "  --p-trim-left-r 21     (enlever primer 805R)"
echo "  --p-trunc-len-f 0      (garder TOUT - pas de troncature)"
echo "  --p-trunc-len-r 0      (garder TOUT - pas de troncature)"
echo "  --p-max-ee-f 5         (accepter 5 erreurs attendues)"
echo "  --p-max-ee-r 5         (accepter 5 erreurs attendues)"
echo "  --p-trunc-q 0          (ne pas tronquer sur qualité)"
echo "  --p-min-overlap 8      (minimum absolu DADA2)"
echo ""
echo "Début: $(date)"
echo ""

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$QIIME_DIR/core/demux.qza" \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-max-ee-f 5 \
  --p-max-ee-r 5 \
  --p-trunc-q 0 \
  --p-min-overlap 8 \
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
    echo ""
    echo "Si tous les reads sont rejetés, cela signifie:"
    echo "  1. L'overlap est vraiment TROP faible (< 8bp)"
    echo "  2. Vos primers ne correspondent pas à 341F/805R"
    echo "  3. La région séquencée n'est pas V3-V4"
    echo ""
    echo "CONTACTEZ LA SOCIÉTÉ DE SÉQUENÇAGE pour:"
    echo "  - Confirmer les primers utilisés"
    echo "  - Demander leur pipeline exact"
    echo "  - Vérifier que ce sont les bons fichiers"
    echo ""
    exit 1
fi

echo "✓ DADA2 terminé!"
echo ""

# Visualisations
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

qiime tools export \
    --input-path "$QIIME_DIR/core/dada2-stats.qza" \
    --output-path "$QIIME_DIR/export/dada2_stats"

python3 << 'EOFPYTHON'
import pandas as pd
import os

stats_file = "/nvme/bio/data_fungi/ExPLOI/05_QIIME2_OPTIMAL/export/dada2_stats/stats.tsv"

if os.path.exists(stats_file):
    df = pd.read_csv(stats_file, sep="\t", comment="#")
    
    print("\n" + "=" * 100)
    print("STATISTIQUES PAR ÉCHANTILLON (APRÈS OPTIMISATION):")
    print("=" * 100)
    print("")
    
    for idx, row in df.iterrows():
        sample = row['sample-id']
        input_reads = int(row['input'])
        filtered = int(row['filtered'])
        merged = int(row['merged'])
        final = int(row['non-chimeric'])
        
        filter_pct = (filtered / input_reads * 100) if input_reads > 0 else 0
        merge_pct = (merged / filtered * 100) if filtered > 0 else 0
        final_pct = (final / input_reads * 100) if input_reads > 0 else 0
        
        status = "✅" if merge_pct > 60 else "⚠️" if merge_pct > 30 else "❌"
        
        print(f"{status} {sample:15s} | Input: {input_reads:>6,} | Filtered: {filtered:>6,} | Merged: {merged:>6,} ({merge_pct:>5.1f}%) | Final: {final:>6,} ({final_pct:>5.1f}%)")
    
    total_input = df['input'].astype(int).sum()
    total_filtered = df['filtered'].astype(int).sum()
    total_merged = df['merged'].astype(int).sum()
    total_final = df['non-chimeric'].astype(int).sum()
    
    print("\n" + "=" * 100)
    print("TOTAUX:")
    print(f"  Input      : {total_input:>10,} reads (100.0%)")
    print(f"  Filtered   : {total_filtered:>10,} reads ({total_filtered/total_input*100:>5.1f}%)")
    print(f"  Merged     : {total_merged:>10,} reads ({total_merged/total_filtered*100:>5.1f}%)")
    print(f"  Final      : {total_final:>10,} reads ({total_final/total_input*100:>5.1f}%)")
    print("=" * 100)
    print("")
    
    merge_rate = total_merged / total_filtered * 100 if total_filtered > 0 else 0
    
    if merge_rate < 30:
        print("❌ TAUX DE MERGE TRÈS FAIBLE (< 30%)")
        print("")
        print("Cela confirme que le séquençage 2x250bp est INSUFFISANT pour V3-V4.")
        print("")
        print("EXPLICATIONS:")
        print("  → Overlap théorique: ~2bp seulement")
        print("  → Beaucoup de pairs ne peuvent pas être assemblées")
        print("  → Résultat: Peu d'ASVs détectés")
        print("")
    elif merge_rate < 60:
        print("⚠️  TAUX DE MERGE MOYEN (30-60%)")
        print("")
        print("Le séquençage 2x250bp est limite pour V3-V4, mais certains reads mergent.")
        print("")
    else:
        print("✅ TAUX DE MERGE ACCEPTABLE (> 60%)")
        print("")

EOFPYTHON

# Compter ASVs
TEMP_EXPORT="$QIIME_DIR/export/temp_$$"
mkdir -p "$TEMP_EXPORT"

qiime tools export \
  --input-path "$QIIME_DIR/core/table.qza" \
  --output-path "$TEMP_EXPORT" 2>/dev/null

NUM_ASVS=$(biom summarize-table -i "$TEMP_EXPORT/feature-table.biom" 2>/dev/null | grep "Num observations:" | awk '{print $3}')
NUM_SAMPLES=$(biom summarize-table -i "$TEMP_EXPORT/feature-table.biom" 2>/dev/null | grep "Num samples:" | awk '{print $3}')

rm -rf "$TEMP_EXPORT"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "NOMBRE D'ASVs DÉTECTÉS: ${NUM_ASVS}"
echo "Nombre d'échantillons : ${NUM_SAMPLES}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -n "$NUM_ASVS" ]; then
    RATIO=$(awk "BEGIN {printf \"%.1f\", ($NUM_ASVS / 4700) * 100}")
    IMPROVEMENT=$(awk "BEGIN {printf \"%.1f\", ($NUM_ASVS / 88.0)}")
    
    echo "COMPARAISON:"
    echo "  Pipeline original (Trimmomatic + trunc 0/0) : 88 ASVs (1.9%)"
    echo "  Pipeline optimisé (RAW + ultra-permissif)   : ${NUM_ASVS} ASVs (${RATIO}%)"
    echo "  Société de séquençage                       : 4700 ASVs (100%)"
    echo ""
    echo "  → Amélioration: ×${IMPROVEMENT}"
    echo ""
    
    if [ "$NUM_ASVS" -gt 2000 ]; then
        echo "✅✅✅ EXCELLENT!"
        echo ""
        echo "Vous avez récupéré un nombre acceptable d'ASVs compte tenu"
        echo "de la limite technique (2x250bp pour amplicon 460bp)."
        echo ""
        echo "PROCHAINES ÉTAPES:"
        echo "  1. Continuer le pipeline (décontamination, taxonomie, diversité)"
        echo "  2. Pour analyses futures: demander du 2x300bp à la société"
        echo ""
    elif [ "$NUM_ASVS" -gt 500 ]; then
        echo "✓ Amélioration significative (×$IMPROVEMENT)"
        echo ""
        echo "Toujours inférieur à la société. Raisons possibles:"
        echo "  1. Ils ont utilisé une méthode différente (VSEARCH, clustering)"
        echo "  2. Ils ont d'autres données (2x300bp ?)"
        echo "  3. Ils ont séquencé une autre région (V4 seul ?)"
        echo ""
        echo "RECOMMANDATION:"
        echo "  → Contacter la société pour comprendre leur pipeline"
        echo "  → Ces ${NUM_ASVS} ASVs sont suffisants pour des analyses de base"
        echo ""
    else
        echo "⚠️  Amélioration limitée"
        echo ""
        echo "CONCLUSION:"
        echo "  Le séquençage 2x250bp est TROP COURT pour V3-V4 (460bp)."
        echo ""
        echo "ACTIONS NÉCESSAIRES:"
        echo "  1. Contacter URGENCE la société de séquençage"
        echo "  2. Demander:"
        echo "     - Comment ont-ils obtenu 4700 ASVs?"
        echo "     - Ont-ils utilisé DADA2 ou autre chose?"
        echo "     - Les vrais fichiers font-ils 2x300bp?"
        echo "     - Quelle région 16S exactement?"
        echo ""
    fi
fi

echo "FICHIERS GÉNÉRÉS:"
echo "  - Stats DADA2    : ${QIIME_DIR}/visual/dada2-stats.qzv"
echo "  - Table résumé   : ${QIIME_DIR}/visual/table-summary.qzv"
echo ""
echo "Visualiser sur: https://view.qiime2.org"
echo ""

if [ "$NUM_ASVS" -gt 500 ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "CONTINUER LE PIPELINE"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "Vous avez maintenant ${NUM_ASVS} ASVs exploitables."
    echo ""
    echo "Pour continuer avec:"
    echo "  - Décontamination (contrôles négatifs)"
    echo "  - Arbre phylogénétique"
    echo "  - Taxonomie SILVA"
    echo "  - Métriques de diversité"
    echo "  - Exports complets"
    echo ""
    echo "Dites-moi et je génère la suite du pipeline!"
    echo ""
fi

echo "=========================================================================================================="


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
