#!/bin/bash

# ========================================================================================================
# Script de NETTOYAGE COMPLET - À exécuter AVANT le pipeline
# ========================================================================================================
# Ce script supprime tous les fichiers générés pour repartir sur une base propre
# Il résout le problème de l'échantillon fantôme "17_Jourand"

BASE_DIR="/nvme/bio/data_fungi/ExPLOI"

echo "========================================"
echo "NETTOYAGE COMPLET DES RÉSULTATS QIIME2"
echo "========================================"
echo ""
echo "⚠️  ATTENTION: Cette opération va supprimer:"
echo "  - Tous les fichiers QIIME2"
echo "  - Les manifests et métadonnées générés"
echo "  - Les fichiers temporaires"
echo ""
read -p "Voulez-vous continuer? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Annulé par l'utilisateur."
    exit 0
fi

echo ""
echo "Suppression en cours..."

# 1. Supprimer TOUS les fichiers QIIME2
if [ -d "${BASE_DIR}/05_QIIME2" ]; then
    echo "✓ Suppression de 05_QIIME2/..."
    rm -rf "${BASE_DIR}/05_QIIME2"
fi
mkdir -p "${BASE_DIR}/05_QIIME2"

# 2. Supprimer les fichiers temporaires
if [ -d "${BASE_DIR}/tmp" ]; then
    echo "✓ Suppression de tmp/..."
    rm -rf "${BASE_DIR}/tmp"
fi

if [ -d "${BASE_DIR}/tmp_mafft" ]; then
    echo "✓ Suppression de tmp_mafft/..."
    rm -rf "${BASE_DIR}/tmp_mafft"
fi

# 3. Supprimer l'ancien manifest et metadata
if [ -f "${BASE_DIR}/manifest_ExPLOI.txt" ]; then
    echo "✓ Suppression de manifest_ExPLOI.txt..."
    rm -f "${BASE_DIR}/manifest_ExPLOI.txt"
fi

if [ -f "${BASE_DIR}/metadata_ExPLOI.tsv" ]; then
    echo "✓ Suppression de metadata_ExPLOI.tsv..."
    rm -f "${BASE_DIR}/metadata_ExPLOI.tsv"
fi

# 4. Vérifier les fichiers nettoyés (03_cleaned_data)
echo ""
echo "Vérification des fichiers nettoyés disponibles:"
echo "-----------------------------------------------"

if [ -d "${BASE_DIR}/03_cleaned_data" ]; then
    CLEANED_FILES=$(ls -1 "${BASE_DIR}/03_cleaned_data"/*_paired.fastq.gz 2>/dev/null | wc -l)
    
    if [ "$CLEANED_FILES" -eq 0 ]; then
        echo "⚠️  ATTENTION: Aucun fichier nettoyé trouvé dans 03_cleaned_data/"
        echo "   Vous devrez peut-être réexécuter les étapes 1-3 du pipeline."
    else
        echo "✓ Trouvé ${CLEANED_FILES} fichiers paired dans 03_cleaned_data/"
        echo ""
        echo "Échantillons détectés:"
        ls -1 "${BASE_DIR}/03_cleaned_data"/*_R1_paired.fastq.gz | \
            xargs -I {} basename {} _R1_paired.fastq.gz | \
            sort | \
            awk '{print "  - " $0}'
    fi
else
    echo "⚠️  ATTENTION: Le dossier 03_cleaned_data/ n'existe pas!"
fi

echo ""
echo "========================================"
echo "✓ NETTOYAGE TERMINÉ"
echo "========================================"
echo ""
echo "Vous pouvez maintenant relancer le pipeline complet:"
echo "  bash pipeline-final-fixed.sh"
echo ""
echo "OU relancer uniquement à partir de QIIME2 (étape 4) si les"
echo "fichiers nettoyés (03_cleaned_data/) sont déjà prêts."
echo "========================================"
