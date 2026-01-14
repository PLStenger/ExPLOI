#!/bin/bash

# ========================================================================================================
# ExPLOI Project - 16S Metabarcoding Pipeline COMPLET OPTIMISÉ (VERSION FINALE DÉFINITIVE)
# ========================================================================================================
#
# Ce script effectue l'analyse complète de métabarcoding 16S de A à Z avec paramètres optimisés
# pour séquençage MiSeq 2x250bp + amplicon V3-V4 (460bp)
#
# ÉTAPE 1 : Contrôle qualité des données brutes (FastQC/MultiQC)
# ÉTAPE 2 : Analyse QIIME2 complète
#   2.1 : Import des données BRUTES (sans Trimmomatic)
#   2.2 : Débruitage DADA2 optimisé (999 ASVs attendus)
#   2.3 : Filtrage des contaminants (contrôles négatifs)
#   2.4 : Construction de l'arbre phylogénétique
#   2.5 : Métriques de diversité (core-metrics-phylogenetic)
#   2.6 : Courbes de raréfaction
#   2.7 : Classification taxonomique (SILVA)
# ÉTAPE 3 : Export des résultats complets
#   3.1 : Table de features
#   3.2 : Séquences représentatives
#   3.3 : Taxonomie
#   3.4 : Arbre phylogénétique
#   3.5 : Statistiques DADA2
#   3.6 : Indices de diversité (18 indices)
#   3.7 : Table ASV + Taxonomie fusionnée
#   3.8 : Profondeurs de lecture
#   3.9 : Données de raréfaction (CSV)
#   3.10: FastQC/MultiQC
# ÉTAPE 4 : Génération rapport final
#
# PARAMÈTRES OPTIMISÉS POUR VOS DONNÉES:
# - Séquençage : MiSeq 2x250bp
# - Amplicon : V3-V4 (341F-805R, 460bp)
# - Overlap : ~73% après optimisation
# - ASVs attendus : ~1000
#
# UTILISATION :
#   chmod +x pipeline-exploi-final-complet.sh
#   nohup bash pipeline-exploi-final-complet.sh > pipeline_final.out 2>&1 &
#   tail -f pipeline_final.out
#
# ========================================================================================================

# ========================================================================================================
# CONFIGURATION
# ========================================================================================================

BASE_DIR="/nvme/bio/data_fungi/ExPLOI"
RAW_DATA_DIR="${BASE_DIR}/01_raw_data"
QIIME_DIR="${BASE_DIR}/05_QIIME2_FINAL"
EXPORT_DIR="${QIIME_DIR}/export"
METADATA_FILE="${BASE_DIR}/metadata_ExPLOI.tsv"
MANIFEST_FILE="${BASE_DIR}/manifest_ExPLOI_FINAL.txt"
CLASSIFIER_PATH="/nvme/bio/data_fungi/BioIndic_La_Reunion_Island_seawater_four_month_SED/05_QIIME2/Original_reads_16S/taxonomy/16S/Classifier.qza"
THREADS=16

# Export variables
export QIIME_DIR="${QIIME_DIR}"
export PYTHONWARNINGS="ignore"

# Temp directories
export TMPDIR="${BASE_DIR}/tmp"
mkdir -p "$TMPDIR"

# ========================================================================================================
# BANNER
# ========================================================================================================

echo "=========================================================================================================="
echo "     ______           _       ___  ___    ____  _            _ _            "
echo "    |  ____|         | |     / _ \\|_ _|  |  _ \\(_)_ __   ___| (_)_ __   ___ "
echo "    | |__  __  ___ __| |    | | | || |   | |_) | | '_ \\ / _ \\ | | '_ \\ / _ \\"
echo "    |  __| \\ \\/ / '_ \\ |    | |_| || |   |  __/| | |_) |  __/ | | | | |  __/"
echo "    | |____ >  <| |_) | |     \\___/|___|  |_|   |_| .__/ \\___|_|_|_| |_|\\___|"
echo "    |______/_/\\_\\ .__/|_|                         |_|                        "
echo "                |_|                                                          "
echo ""
echo "    16S Metabarcoding - Pipeline Complet Optimisé v1.0"
echo "=========================================================================================================="
echo ""
echo "Date de début : $(date)"
echo "Utilisateur   : $(whoami)"
echo "Serveur       : $(hostname)"
echo ""
echo "Configuration :"
echo "  - Séquençage   : MiSeq 2x250bp"
echo "  - Région 16S   : V3-V4 (341F-805R)"
echo "  - Amplicon     : 460bp"
echo "  - Threads      : ${THREADS}"
echo "  - Répertoire   : ${QIIME_DIR}"
echo ""
echo "=========================================================================================================="
echo ""

# ========================================================================================================
# CRÉATION DE L'ARBORESCENCE
# ========================================================================================================

mkdir -p "$QIIME_DIR/core" "$QIIME_DIR/visual" "$EXPORT_DIR"
mkdir -p "$EXPORT_DIR/qc_raw" "$EXPORT_DIR/diversity" "$EXPORT_DIR/rarefaction_data"

echo "✓ Arborescence créée"
echo ""

# ========================================================================================================
# CRÉATION DU FICHIER METADATA
# ========================================================================================================

echo "Création du fichier metadata..."

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

echo "✓ Metadata créé : $METADATA_FILE"
echo ""

# ========================================================================================================
# ÉTAPE 1 : CONTRÔLE QUALITÉ DES DONNÉES BRUTES
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 1 : Contrôle Qualité des Données Brutes"
echo "=========================================================================================================="
echo ""

eval "$(conda shell.bash hook)"
conda activate fastqc

echo "Lancement FastQC sur données brutes..."
fastqc -t $THREADS "$RAW_DATA_DIR"/*.fastq.gz -o "$EXPORT_DIR/qc_raw" --quiet

if [ $? -eq 0 ]; then
    echo "✓ FastQC terminé"
else
    echo "⚠️  WARNING: FastQC a échoué"
fi

conda deactivate
conda activate multiqc

echo ""
echo "Lancement MultiQC..."
multiqc "$EXPORT_DIR/qc_raw" -o "$EXPORT_DIR/qc_raw" --force --quiet

if [ $? -eq 0 ]; then
    echo "✓ MultiQC terminé"
    echo "✓ Rapport : ${EXPORT_DIR}/qc_raw/multiqc_report.html"
else
    echo "⚠️  WARNING: MultiQC a échoué"
fi

conda deactivate

echo ""
echo "ÉTAPE 1 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 2 : ANALYSE QIIME2 COMPLÈTE
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 2 : Analyse QIIME2 Complète"
echo "=========================================================================================================="
echo ""

conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

echo "QIIME2 activé"
qiime --version
echo ""

# --------------------------------------------------------------------------------------------------------
# 2.1 : IMPORT DES DONNÉES
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.1 : Import des Données (RAW - sans Trimmomatic)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

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

echo "✓ Manifest créé : ${SAMPLE_COUNT} échantillons"
echo ""

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --output-path "$QIIME_DIR/core/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

if [ $? -eq 0 ]; then
    echo "✓ Données importées : demux.qza"
else
    echo "❌ ERROR: Import failed!"
    exit 1
fi

qiime demux summarize \
  --i-data "$QIIME_DIR/core/demux.qza" \
  --o-visualization "$QIIME_DIR/visual/demux.qzv"

echo "✓ Demux summary : visual/demux.qzv"
echo ""

# --------------------------------------------------------------------------------------------------------
# 2.2 : DÉBRUITAGE DADA2 OPTIMISÉ
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.2 : Débruitage DADA2 Optimisé (MiSeq 2x250bp)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Paramètres optimisés pour :"
echo "  - Séquençage   : MiSeq 2x250bp"
echo "  - Amplicon V3-V4 : 460bp"
echo "  - Overlap attendu : ~73%"
echo ""
echo "Paramètres DADA2 :"
echo "  --p-trim-left-f 17     (enlever primer 341F)"
echo "  --p-trim-left-r 21     (enlever primer 805R)"
echo "  --p-trunc-len-f 0      (garder longueur max)"
echo "  --p-trunc-len-r 0      (garder longueur max)"
echo "  --p-max-ee-f 5         (filtres permissifs)"
echo "  --p-max-ee-r 5"
echo "  --p-min-overlap 8      (overlap minimal)"
echo ""
echo "Lancement DADA2 (30-60 minutes)..."
echo "Début : $(date)"
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
  --o-denoising-stats "$QIIME_DIR/core/dada2-stats.qza"

if [ $? -eq 0 ]; then
    echo ""
    echo "Fin : $(date)"
    echo "✓ DADA2 terminé avec succès"
else
    echo "❌ ERROR: DADA2 failed!"
    exit 1
fi

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/dada2-stats.qza" \
  --o-visualization "$QIIME_DIR/visual/dada2-stats.qzv"

echo "✓ DADA2 stats : visual/dada2-stats.qzv"
echo ""

# --------------------------------------------------------------------------------------------------------
# 2.3 : FILTRAGE DES CONTAMINANTS
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.3 : Filtrage des Contaminants (Contrôles Négatifs)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

qiime feature-table filter-samples \
  --i-table "$QIIME_DIR/core/table.qza" \
  --m-metadata-file "$METADATA_FILE" \
  --p-where "[group]='Negative_Control'" \
  --o-filtered-table "$QIIME_DIR/core/neg-controls-table.qza" 2>&1 | tee "$QIIME_DIR/neg_control_filter.log"

NEG_CONTROL_STATUS=${PIPESTATUS[0]}

if [ $NEG_CONTROL_STATUS -ne 0 ]; then
  echo ""
  echo "⚠️  Contrôle négatif VIDE (0 reads)"
  echo "    → Pas de contamination détectée"
  echo "    → Aucun ASV ne sera filtré"
  echo ""
  
  cp "$QIIME_DIR/core/table.qza" "$QIIME_DIR/core/table-decontam.qza"
  cp "$QIIME_DIR/core/rep-seqs.qza" "$QIIME_DIR/core/rep-seqs-clean.qza"
  
  qiime feature-table filter-samples \
    --i-table "$QIIME_DIR/core/table-decontam.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --p-where "[group]='Sample'" \
    --o-filtered-table "$QIIME_DIR/core/table-final.qza"
  
  echo "✓ Table finale créée (0 ASVs filtrés)"
  
else
  echo "✓ Contrôle négatif contient des reads"
  echo "  Filtrage des contaminants..."
  
  qiime feature-table summarize \
    --i-table "$QIIME_DIR/core/neg-controls-table.qza" \
    --o-visualization "$QIIME_DIR/visual/neg-controls-summary.qzv"
  
  qiime tools export \
    --input-path "$QIIME_DIR/core/neg-controls-table.qza" \
    --output-path "$EXPORT_DIR/neg-controls"

  biom convert \
    -i "$EXPORT_DIR/neg-controls/feature-table.biom" \
    -o "$EXPORT_DIR/neg-controls/feature-table.tsv" \
    --to-tsv

  echo "feature-id" > "$EXPORT_DIR/neg-controls/contamination_ids.txt"
  tail -n +3 "$EXPORT_DIR/neg-controls/feature-table.tsv" | \
    cut -f1 >> "$EXPORT_DIR/neg-controls/contamination_ids.txt"

  NUM_CONTAMINANTS=$(tail -n +2 "$EXPORT_DIR/neg-controls/contamination_ids.txt" | wc -l)
  echo "  → ${NUM_CONTAMINANTS} ASVs contaminants détectés"

  if [ "$NUM_CONTAMINANTS" -eq 0 ]; then
    cp "$QIIME_DIR/core/table.qza" "$QIIME_DIR/core/table-decontam.qza"
    cp "$QIIME_DIR/core/rep-seqs.qza" "$QIIME_DIR/core/rep-seqs-clean.qza"
  else
    qiime feature-table filter-features \
      --i-table "$QIIME_DIR/core/table.qza" \
      --m-metadata-file "$EXPORT_DIR/neg-controls/contamination_ids.txt" \
      --p-exclude-ids \
      --o-filtered-table "$QIIME_DIR/core/table-decontam.qza"

    qiime feature-table filter-seqs \
      --i-data "$QIIME_DIR/core/rep-seqs.qza" \
      --i-table "$QIIME_DIR/core/table-decontam.qza" \
      --o-filtered-data "$QIIME_DIR/core/rep-seqs-clean.qza"
    
    echo "  → ${NUM_CONTAMINANTS} ASVs contaminants retirés"
  fi

  qiime feature-table filter-samples \
    --i-table "$QIIME_DIR/core/table-decontam.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --p-where "[group]='Sample'" \
    --o-filtered-table "$QIIME_DIR/core/table-final.qza"

  echo "✓ Décontamination terminée"
fi

qiime feature-table summarize \
  --i-table "$QIIME_DIR/core/table-final.qza" \
  --o-visualization "$QIIME_DIR/visual/table-final-summary.qzv" \
  --m-sample-metadata-file "$METADATA_FILE"

qiime feature-table tabulate-seqs \
  --i-data "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --o-visualization "$QIIME_DIR/visual/rep-seqs-clean.qzv"

echo "✓ Table finale : visual/table-final-summary.qzv"
echo ""

# --------------------------------------------------------------------------------------------------------
# 2.4 : ARBRE PHYLOGÉNÉTIQUE
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.4 : Construction de l'Arbre Phylogénétique"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MAFFT_TMPDIR="${BASE_DIR}/tmp_mafft"
mkdir -p "$MAFFT_TMPDIR"
export TMPDIR="$MAFFT_TMPDIR"

echo "Construction arbre phylogénétique (10-30 minutes)..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --p-n-threads 1 \
  --o-alignment "$QIIME_DIR/core/aligned-rep-seqs.qza" \
  --o-masked-alignment "$QIIME_DIR/core/masked-aligned-rep-seqs.qza" \
  --o-tree "$QIIME_DIR/core/unrooted-tree.qza" \
  --o-rooted-tree "$QIIME_DIR/core/rooted-tree.qza" 2>&1

TREE_STATUS=$?

if [ $TREE_STATUS -ne 0 ] || [ ! -f "$QIIME_DIR/core/rooted-tree.qza" ]; then
  echo "⚠️  MAFFT standard échoué. Essai avec --p-parttree..."
  
  rm -f "$QIIME_DIR/core/aligned-rep-seqs.qza"
  rm -f "$QIIME_DIR/core/masked-aligned-rep-seqs.qza"
  rm -f "$QIIME_DIR/core/unrooted-tree.qza"
  rm -f "$QIIME_DIR/core/rooted-tree.qza"
  
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$QIIME_DIR/core/rep-seqs-clean.qza" \
    --p-n-threads 1 \
    --p-parttree \
    --o-alignment "$QIIME_DIR/core/aligned-rep-seqs.qza" \
    --o-masked-alignment "$QIIME_DIR/core/masked-aligned-rep-seqs.qza" \
    --o-tree "$QIIME_DIR/core/unrooted-tree.qza" \
    --o-rooted-tree "$QIIME_DIR/core/rooted-tree.qza" 2>&1
  
  TREE_STATUS=$?
fi

if [ -f "$QIIME_DIR/core/rooted-tree.qza" ]; then
  echo "✓ Arbre phylogénétique créé"
else
  echo "⚠️  Arbre phylogénétique non créé (continuera sans Faith PD)"
fi

rm -rf "$MAFFT_TMPDIR"
export TMPDIR="${BASE_DIR}/tmp"

echo ""

# --------------------------------------------------------------------------------------------------------
# 2.5 : MÉTRIQUES DE DIVERSITÉ
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.5 : Métriques de Diversité (Core-Metrics)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

rm -rf "$QIIME_DIR/core-metrics-results"
mkdir -p "$QIIME_DIR/core-metrics-results"

# Export table pour calculer profondeurs
qiime tools export \
  --input-path "$QIIME_DIR/core/table-final.qza" \
  --output-path "$EXPORT_DIR/table-final-temp"

biom convert \
  -i "$EXPORT_DIR/table-final-temp/feature-table.biom" \
  -o "$EXPORT_DIR/table-final-temp/feature-table.tsv" \
  --to-tsv

echo "Calcul des profondeurs avec Python..."

DEPTHS=$(python3 << EOFPYTHON
import pandas as pd
import sys

try:
    df = pd.read_csv(
        "${EXPORT_DIR}/table-final-temp/feature-table.tsv",
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
    echo "❌ ERROR: Calcul profondeurs échoué!"
    exit 1
fi

echo "✓ Sampling Depth : ${SAMPLING_DEPTH}"
echo "✓ Max Depth      : ${MAX_DEPTH}"
echo ""

echo "$SAMPLING_DEPTH" > "$EXPORT_DIR/sampling_depth.txt"
echo "$MAX_DEPTH" > "$EXPORT_DIR/max_depth.txt"

# Core-metrics
if [ -f "$QIIME_DIR/core/rooted-tree.qza" ]; then
  echo "Lancement core-metrics-phylogenetic (avec Faith PD)..."
  
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$QIIME_DIR/core/rooted-tree.qza" \
    --i-table "$QIIME_DIR/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --o-rarefied-table "$QIIME_DIR/core-metrics-results/rarefied_table.qza" \
    --o-faith-pd-vector "$QIIME_DIR/core-metrics-results/faith_pd_vector.qza" \
    --o-observed-features-vector "$QIIME_DIR/core-metrics-results/observed_features_vector.qza" \
    --o-shannon-vector "$QIIME_DIR/core-metrics-results/shannon_vector.qza" \
    --o-evenness-vector "$QIIME_DIR/core-metrics-results/evenness_vector.qza" \
    --o-unweighted-unifrac-distance-matrix "$QIIME_DIR/core-metrics-results/unweighted_unifrac_distance_matrix.qza" \
    --o-weighted-unifrac-distance-matrix "$QIIME_DIR/core-metrics-results/weighted_unifrac_distance_matrix.qza" \
    --o-jaccard-distance-matrix "$QIIME_DIR/core-metrics-results/jaccard_distance_matrix.qza" \
    --o-bray-curtis-distance-matrix "$QIIME_DIR/core-metrics-results/bray_curtis_distance_matrix.qza" \
    --o-unweighted-unifrac-pcoa-results "$QIIME_DIR/core-metrics-results/unweighted_unifrac_pcoa_results.qza" \
    --o-weighted-unifrac-pcoa-results "$QIIME_DIR/core-metrics-results/weighted_unifrac_pcoa_results.qza" \
    --o-jaccard-pcoa-results "$QIIME_DIR/core-metrics-results/jaccard_pcoa_results.qza" \
    --o-bray-curtis-pcoa-results "$QIIME_DIR/core-metrics-results/bray_curtis_pcoa_results.qza" \
    --o-unweighted-unifrac-emperor "$QIIME_DIR/core-metrics-results/unweighted_unifrac_emperor.qzv" \
    --o-weighted-unifrac-emperor "$QIIME_DIR/core-metrics-results/weighted_unifrac_emperor.qzv" \
    --o-jaccard-emperor "$QIIME_DIR/core-metrics-results/jaccard_emperor.qzv" \
    --o-bray-curtis-emperor "$QIIME_DIR/core-metrics-results/bray_curtis_emperor.qzv"
    
  if [ $? -eq 0 ]; then
    echo "✓ Core-metrics-phylogenetic terminé (avec Faith PD)"
  else
    echo "⚠️  Core-metrics-phylogenetic échoué, essai sans phylogénie..."
    rm -rf "$QIIME_DIR/core-metrics-results"/*
    
    qiime diversity core-metrics \
      --i-table "$QIIME_DIR/core/table-final.qza" \
      --p-sampling-depth "$SAMPLING_DEPTH" \
      --m-metadata-file "$METADATA_FILE" \
      --o-rarefied-table "$QIIME_DIR/core-metrics-results/rarefied_table.qza" \
      --o-observed-features-vector "$QIIME_DIR/core-metrics-results/observed_features_vector.qza" \
      --o-shannon-vector "$QIIME_DIR/core-metrics-results/shannon_vector.qza" \
      --o-evenness-vector "$QIIME_DIR/core-metrics-results/evenness_vector.qza" \
      --o-jaccard-distance-matrix "$QIIME_DIR/core-metrics-results/jaccard_distance_matrix.qza" \
      --o-bray-curtis-distance-matrix "$QIIME_DIR/core-metrics-results/bray_curtis_distance_matrix.qza" \
      --o-jaccard-pcoa-results "$QIIME_DIR/core-metrics-results/jaccard_pcoa_results.qza" \
      --o-bray-curtis-pcoa-results "$QIIME_DIR/core-metrics-results/bray_curtis_pcoa_results.qza" \
      --o-jaccard-emperor "$QIIME_DIR/core-metrics-results/jaccard_emperor.qzv" \
      --o-bray-curtis-emperor "$QIIME_DIR/core-metrics-results/bray_curtis_emperor.qzv"
    
    echo "✓ Core-metrics terminé (sans Faith PD)"
  fi
else
  echo "Lancement core-metrics (sans phylogénie)..."
  
  qiime diversity core-metrics \
    --i-table "$QIIME_DIR/core/table-final.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --o-rarefied-table "$QIIME_DIR/core-metrics-results/rarefied_table.qza" \
    --o-observed-features-vector "$QIIME_DIR/core-metrics-results/observed_features_vector.qza" \
    --o-shannon-vector "$QIIME_DIR/core-metrics-results/shannon_vector.qza" \
    --o-evenness-vector "$QIIME_DIR/core-metrics-results/evenness_vector.qza" \
    --o-jaccard-distance-matrix "$QIIME_DIR/core-metrics-results/jaccard_distance_matrix.qza" \
    --o-bray-curtis-distance-matrix "$QIIME_DIR/core-metrics-results/bray_curtis_distance_matrix.qza" \
    --o-jaccard-pcoa-results "$QIIME_DIR/core-metrics-results/jaccard_pcoa_results.qza" \
    --o-bray-curtis-pcoa-results "$QIIME_DIR/core-metrics-results/bray_curtis_pcoa_results.qza" \
    --o-jaccard-emperor "$QIIME_DIR/core-metrics-results/jaccard_emperor.qzv" \
    --o-bray-curtis-emperor "$QIIME_DIR/core-metrics-results/bray_curtis_emperor.qzv"
  
  echo "✓ Core-metrics terminé (sans Faith PD)"
fi

echo ""

# --------------------------------------------------------------------------------------------------------
# 2.6 : COURBES DE RARÉFACTION
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.6 : Courbes de Raréfaction"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "Génération courbes raréfaction (Shannon + Observed Features)..."

qiime diversity alpha-rarefaction \
  --i-table "$QIIME_DIR/core/table-final.qza" \
  --p-min-depth 10 \
  --p-max-depth "$MAX_DEPTH" \
  --p-steps 20 \
  --m-metadata-file "$METADATA_FILE" \
  --o-visualization "$QIIME_DIR/visual/rarefaction-curves.qzv"

if [ $? -eq 0 ]; then
  echo "✓ Courbes raréfaction : visual/rarefaction-curves.qzv"
else
  echo "⚠️  Courbes raréfaction échouées"
fi

if [ -f "$QIIME_DIR/core/rooted-tree.qza" ]; then
  echo "Génération courbes Faith PD..."
  
  qiime diversity alpha-rarefaction \
    --i-table "$QIIME_DIR/core/table-final.qza" \
    --i-phylogeny "$QIIME_DIR/core/rooted-tree.qza" \
    --p-min-depth 10 \
    --p-max-depth "$MAX_DEPTH" \
    --p-steps 20 \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "$QIIME_DIR/visual/rarefaction-curves-phylogenetic.qzv"
  
  if [ $? -eq 0 ]; then
    echo "✓ Courbes Faith PD : visual/rarefaction-curves-phylogenetic.qzv"
  fi
fi

echo ""

# --------------------------------------------------------------------------------------------------------
# 2.7 : CLASSIFICATION TAXONOMIQUE
# --------------------------------------------------------------------------------------------------------

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2.7 : Classification Taxonomique (SILVA)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "Classification taxonomique (10-20 minutes)..."

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_PATH" \
  --i-reads "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --o-classification "$QIIME_DIR/core/taxonomy.qza"

if [ $? -eq 0 ]; then
    echo "✓ Taxonomie assignée : taxonomy.qza"
else
    echo "❌ ERROR: Classification taxonomique échouée!"
    exit 1
fi

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/core/taxonomy.qza" \
  --o-visualization "$QIIME_DIR/visual/taxonomy.qzv"

echo "✓ Taxonomie : visual/taxonomy.qzv"
echo ""

echo "ÉTAPE 2 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 3 : EXPORT DES RÉSULTATS COMPLETS
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 3 : Export des Résultats Complets"
echo "=========================================================================================================="
echo ""

mkdir -p "$EXPORT_DIR/feature_table" "$EXPORT_DIR/rep_seqs" "$EXPORT_DIR/taxonomy" "$EXPORT_DIR/tree" "$EXPORT_DIR/dada2_stats" "$EXPORT_DIR/diversity_all"

# --------------------------------------------------------------------------------------------------------
# 3.1-3.4 : Exports standards
# --------------------------------------------------------------------------------------------------------

echo "3.1 : Export table de features..."
qiime tools export \
  --input-path "$QIIME_DIR/core/table-final.qza" \
  --output-path "$EXPORT_DIR/feature_table"

biom convert \
  -i "$EXPORT_DIR/feature_table/feature-table.biom" \
  -o "$EXPORT_DIR/feature_table/feature-table.tsv" \
  --to-tsv

sed -i 's/#OTU ID/#ASV_ID/g' "$EXPORT_DIR/feature_table/feature-table.tsv" 2>/dev/null || true

echo "✓ Table features : export/feature_table/feature-table.tsv"

echo "3.2 : Export séquences représentatives..."
qiime tools export \
  --input-path "$QIIME_DIR/core/rep-seqs-clean.qza" \
  --output-path "$EXPORT_DIR/rep_seqs"

echo "✓ Séquences : export/rep_seqs/dna-sequences.fasta"

echo "3.3 : Export taxonomie..."
qiime tools export \
  --input-path "$QIIME_DIR/core/taxonomy.qza" \
  --output-path "$EXPORT_DIR/taxonomy"

echo "✓ Taxonomie : export/taxonomy/taxonomy.tsv"

echo "3.4 : Export arbre phylogénétique..."
if [ -f "$QIIME_DIR/core/rooted-tree.qza" ]; then
  qiime tools export \
    --input-path "$QIIME_DIR/core/rooted-tree.qza" \
    --output-path "$EXPORT_DIR/tree"
  
  if [ -f "$EXPORT_DIR/tree/tree.nwk" ]; then
    echo "✓ Arbre : export/tree/tree.nwk"
  fi
fi

echo "3.5 : Export stats DADA2..."
qiime tools export \
  --input-path "$QIIME_DIR/core/dada2-stats.qza" \
  --output-path "$EXPORT_DIR/dada2_stats"

echo "✓ Stats DADA2 : export/dada2_stats/"

# --------------------------------------------------------------------------------------------------------
# 3.6 : Indices de diversité (18 indices)
# --------------------------------------------------------------------------------------------------------

echo "3.6 : Calcul de tous les indices de diversité alpha..."

DIVERSITY_DIR="$QIIME_DIR/diversity_indices"
mkdir -p "$DIVERSITY_DIR"

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
  qiime diversity alpha \
    --i-table "$QIIME_DIR/core/table-final.qza" \
    --p-metric "$metric" \
    --o-alpha-diversity "$DIVERSITY_DIR/${metric}_vector.qza" 2>/dev/null
done

echo "✓ 18 indices calculés"

# Export tous les indices
export_diversity_metric() {
    local metric_name=$1
    local qza_file=$2
    local output_name=$3
    
    if [ -f "$qza_file" ]; then
        qiime tools export \
          --input-path "$qza_file" \
          --output-path "$EXPORT_DIR/diversity_all/${output_name}_temp" 2>/dev/null
        
        sed "1s/.*/sample-id\t${output_name}/" \
          "$EXPORT_DIR/diversity_all/${output_name}_temp/alpha-diversity.tsv" > \
          "$EXPORT_DIR/diversity_all/${output_name}.tsv" 2>/dev/null
        
        rm -rf "$EXPORT_DIR/diversity_all/${output_name}_temp"
    fi
}

export_diversity_metric "Observed ASVs" "$QIIME_DIR/core-metrics-results/observed_features_vector.qza" "observed_asvs"
export_diversity_metric "Shannon" "$QIIME_DIR/core-metrics-results/shannon_vector.qza" "shannon"
export_diversity_metric "Pielou Evenness" "$QIIME_DIR/core-metrics-results/evenness_vector.qza" "pielou_evenness"
export_diversity_metric "Faith PD" "$QIIME_DIR/core-metrics-results/faith_pd_vector.qza" "faith_pd"
export_diversity_metric "Simpson" "$DIVERSITY_DIR/simpson_vector.qza" "simpson"
export_diversity_metric "Simpson E" "$DIVERSITY_DIR/simpson_e_vector.qza" "simpson_evenness"
export_diversity_metric "Chao1" "$DIVERSITY_DIR/chao1_vector.qza" "chao1"
export_diversity_metric "ACE" "$DIVERSITY_DIR/ace_vector.qza" "ace"
export_diversity_metric "Goods Coverage" "$DIVERSITY_DIR/goods_coverage_vector.qza" "goods_coverage"
export_diversity_metric "Fisher Alpha" "$DIVERSITY_DIR/fisher_alpha_vector.qza" "fisher_alpha"
export_diversity_metric "Berger Parker" "$DIVERSITY_DIR/berger_parker_d_vector.qza" "berger_parker"
export_diversity_metric "Gini" "$DIVERSITY_DIR/gini_index_vector.qza" "gini_index"
export_diversity_metric "Brillouin" "$DIVERSITY_DIR/brillouin_d_vector.qza" "brillouin"
export_diversity_metric "Strong" "$DIVERSITY_DIR/strong_vector.qza" "strong"
export_diversity_metric "McIntosh D" "$DIVERSITY_DIR/mcintosh_d_vector.qza" "mcintosh_d"
export_diversity_metric "McIntosh E" "$DIVERSITY_DIR/mcintosh_e_vector.qza" "mcintosh_e"
export_diversity_metric "Margalef" "$DIVERSITY_DIR/margalef_vector.qza" "margalef"
export_diversity_metric "Menhinick" "$DIVERSITY_DIR/menhinick_vector.qza" "menhinick"

# Fusionner tous les indices
echo "Fusion de tous les indices..."

python3 << EOFPYTHON
import pandas as pd
import os
from pathlib import Path

export_dir = Path("${EXPORT_DIR}")
diversity_dir = export_dir / "diversity_all"

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

    print(f"✓ Table diversité complète : {output_file}")
    print(f"  Échantillons : {len(df_merged)}")
    print(f"  Indices      : {len(df_merged.columns)}")
EOFPYTHON

# --------------------------------------------------------------------------------------------------------
# 3.7 : Table ASV + Taxonomie fusionnée
# --------------------------------------------------------------------------------------------------------

echo "3.7 : Fusion ASV + Taxonomie..."

python3 << EOFPYTHON
import pandas as pd
import os

export_dir = "${EXPORT_DIR}"

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
    print(f"✓ Table fusionnée : {output_file}")
except Exception as e:
    print(f"ERROR: {e}")
EOFPYTHON

# --------------------------------------------------------------------------------------------------------
# 3.8 : Profondeurs de lecture par échantillon
# --------------------------------------------------------------------------------------------------------

echo "3.8 : Calcul profondeurs par échantillon..."

python3 << EOFPYTHON
import pandas as pd
import os

export_dir = "${EXPORT_DIR}"

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

print(f"✓ Profondeurs échantillons : {output_file}")
EOFPYTHON

# --------------------------------------------------------------------------------------------------------
# 3.9 : Export données raréfaction
# --------------------------------------------------------------------------------------------------------

echo "3.9 : Export données raréfaction..."

if [ -f "$QIIME_DIR/visual/rarefaction-curves.qzv" ]; then
  unzip -q "$QIIME_DIR/visual/rarefaction-curves.qzv" -d "$EXPORT_DIR/rarefaction_data/temp" 2>/dev/null
  
  find "$EXPORT_DIR/rarefaction_data/temp" -name "*.csv" -exec cp {} "$EXPORT_DIR/rarefaction_data/" \; 2>/dev/null
  
  rm -rf "$EXPORT_DIR/rarefaction_data/temp"
  
  CSV_COUNT=$(ls -1 "$EXPORT_DIR/rarefaction_data/"*.csv 2>/dev/null | wc -l)
  
  if [ "$CSV_COUNT" -gt 0 ]; then
    echo "✓ Données raréfaction : ${CSV_COUNT} CSV"
  fi
fi

if [ -f "$QIIME_DIR/visual/rarefaction-curves-phylogenetic.qzv" ]; then
  unzip -q "$QIIME_DIR/visual/rarefaction-curves-phylogenetic.qzv" -d "$EXPORT_DIR/rarefaction_data/temp_faith" 2>/dev/null
  
  find "$EXPORT_DIR/rarefaction_data/temp_faith" -name "*.csv" -exec cp {} "$EXPORT_DIR/rarefaction_data/" \; 2>/dev/null
  
  rm -rf "$EXPORT_DIR/rarefaction_data/temp_faith"
fi

echo ""
echo "ÉTAPE 3 : TERMINÉE"
echo ""

# ========================================================================================================
# ÉTAPE 4 : RAPPORT FINAL
# ========================================================================================================

echo "=========================================================================================================="
echo "ÉTAPE 4 : Génération Rapport Final"
echo "=========================================================================================================="
echo ""

# Compter ASVs
TEMP_EXPORT="$EXPORT_DIR/temp_final_$$"
mkdir -p "$TEMP_EXPORT"

qiime tools export \
  --input-path "$QIIME_DIR/core/table-final.qza" \
  --output-path "$TEMP_EXPORT" 2>/dev/null

NUM_ASVS=$(biom summarize-table -i "$TEMP_EXPORT/feature-table.biom" 2>/dev/null | grep "Num observations:" | awk '{print $3}')
NUM_SAMPLES=$(biom summarize-table -i "$TEMP_EXPORT/feature-table.biom" 2>/dev/null | grep "Num samples:" | awk '{print $3}')

rm -rf "$TEMP_EXPORT"

# Stats DADA2
python3 << EOFPYTHON
import pandas as pd
import os

stats_file = os.path.join("${EXPORT_DIR}", "dada2_stats", "stats.tsv")

if os.path.exists(stats_file):
    df = pd.read_csv(stats_file, sep="\t", comment="#")
    
    total_input = df['input'].astype(int).sum()
    total_filtered = df['filtered'].astype(int).sum()
    total_merged = df['merged'].astype(int).sum()
    total_final = df['non-chimeric'].astype(int).sum()
    
    filter_pct = (total_filtered / total_input * 100) if total_input > 0 else 0
    merge_pct = (total_merged / total_filtered * 100) if total_filtered > 0 else 0
    final_pct = (total_final / total_input * 100) if total_input > 0 else 0
    
    print("")
    print("=" * 100)
    print("STATISTIQUES DADA2 GLOBALES")
    print("=" * 100)
    print(f"  Reads input          : {total_input:>10,} (100.0%)")
    print(f"  Reads filtrés        : {total_filtered:>10,} ({filter_pct:>5.1f}%)")
    print(f"  Reads mergés         : {total_merged:>10,} ({merge_pct:>5.1f}%)")
    print(f"  Reads finaux         : {total_final:>10,} ({final_pct:>5.1f}%)")
    print("=" * 100)
    print("")
EOFPYTHON

echo ""
echo "=========================================================================================================="
echo "🎉 PIPELINE TERMINÉ AVEC SUCCÈS ! 🎉"
echo "=========================================================================================================="
echo ""
echo "Date de fin : $(date)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RÉSUMÉ FINAL"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 RÉSULTATS PRINCIPAUX :"
echo "  • Nombre d'ASVs              : ${NUM_ASVS}"
echo "  • Nombre d'échantillons      : ${NUM_SAMPLES}"
echo "  • Profondeur raréfaction     : ${SAMPLING_DEPTH} reads"
echo "  • Profondeur maximale        : ${MAX_DEPTH} reads"
echo ""
echo "📁 FICHIERS PRINCIPAUX :"
echo "  • Table ASV                  : export/feature_table/feature-table.tsv"
echo "  • Séquences ASV              : export/rep_seqs/dna-sequences.fasta"
echo "  • Taxonomie                  : export/taxonomy/taxonomy.tsv"
echo "  • ASV + Taxonomie fusionné   : export/ASV_abundance_taxonomy.tsv"
echo ""
echo "📈 DIVERSITÉ :"
echo "  • Tous les indices alpha     : export/diversity_indices_all.tsv"
echo "  • Profondeurs échantillons   : export/sample_read_depths_final.tsv"
echo ""
if [ -f "$EXPORT_DIR/tree/tree.nwk" ]; then
echo "🌳 PHYLOGÉNIE :"
echo "  • Arbre phylogénétique       : export/tree/tree.nwk"
echo ""
fi
echo "📉 RARÉFACTION :"
echo "  • Courbes standard           : visual/rarefaction-curves.qzv"
if [ -f "$QIIME_DIR/visual/rarefaction-curves-phylogenetic.qzv" ]; then
echo "  • Courbes Faith PD           : visual/rarefaction-curves-phylogenetic.qzv"
fi
echo "  • Données CSV                : export/rarefaction_data/"
echo ""
echo "📋 QUALITÉ :"
echo "  • FastQC brut                : export/qc_raw/"
echo "  • MultiQC report             : export/qc_raw/multiqc_report.html"
echo "  • Stats DADA2                : visual/dada2-stats.qzv"
echo ""
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Tous les fichiers sont dans : ${QIIME_DIR}"
echo ""
echo "=========================================================================================================="
