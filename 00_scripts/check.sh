#!/bin/bash

# ========================================================================================================
# DÉTECTION AUTOMATIQUE DE LA LONGUEUR DES READS
# ========================================================================================================

CLEANED_DIR="/nvme/bio/data_fungi/ExPLOI/03_cleaned_data"

echo "=========================================================================================================="
echo "DÉTECTION DE LA LONGUEUR DES READS"
echo "=========================================================================================================="
echo ""

if [ ! -d "$CLEANED_DIR" ]; then
    echo "❌ ERROR: Répertoire introuvable: $CLEANED_DIR"
    exit 1
fi

echo "Analyse des reads nettoyés (après Trimmomatic)..."
echo ""

# Prendre le premier fichier R1 et R2
R1_FILE=$(find "$CLEANED_DIR" -name "*_R1_paired.fastq.gz" | head -1)
R2_FILE=$(find "$CLEANED_DIR" -name "*_R2_paired.fastq.gz" | head -1)

if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ]; then
    echo "❌ ERROR: Aucun fichier paired trouvé dans $CLEANED_DIR"
    exit 1
fi

echo "Fichiers analysés:"
echo "  R1: $(basename $R1_FILE)"
echo "  R2: $(basename $R2_FILE)"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ANALYSE DE LA LONGUEUR DES READS (1000 premiers reads)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Analyser R1
echo "R1 (Forward):"
zcat "$R1_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | tail -20

R1_MIN=$(zcat "$R1_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | head -1)
R1_MAX=$(zcat "$R1_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | tail -1)
R1_MEDIAN=$(zcat "$R1_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | awk '{arr[NR]=$1} END {print arr[int(NR/2)]}')

echo ""
echo "R2 (Reverse):"
zcat "$R2_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | tail -20

R2_MIN=$(zcat "$R2_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | head -1)
R2_MAX=$(zcat "$R2_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | tail -1)
R2_MEDIAN=$(zcat "$R2_FILE" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | awk '{arr[NR]=$1} END {print arr[int(NR/2)]}')

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RÉSUMÉ"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "R1 (Forward):"
echo "  Min     : ${R1_MIN} bp"
echo "  Max     : ${R1_MAX} bp"
echo "  Médiane : ${R1_MEDIAN} bp"
echo ""
echo "R2 (Reverse):"
echo "  Min     : ${R2_MIN} bp"
echo "  Max     : ${R2_MAX} bp"
echo "  Médiane : ${R2_MEDIAN} bp"
echo ""

# Calculer paramètres optimaux
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PARAMÈTRES DADA2 RECOMMANDÉS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Hypothèse : Amplicon V3-V4 = 460bp
AMPLICON_LENGTH=460

# Calculer troncature pour avoir overlap de 60-80bp
# Formule: overlap = trunc_f + trunc_r - amplicon
# On veut overlap = 70bp
# Donc: trunc_f + trunc_r = amplicon + 70 = 530

# Stratégie : garder au maximum sans dépasser la longueur médiane
if [ "$R1_MEDIAN" -lt 250 ]; then
    TRUNC_F=$((R1_MEDIAN - 10))
else
    TRUNC_F=240
fi

if [ "$R2_MEDIAN" -lt 250 ]; then
    TRUNC_R=$((R2_MEDIAN - 10))
else
    TRUNC_R=240
fi

OVERLAP=$((TRUNC_F + TRUNC_R - AMPLICON_LENGTH))

echo "Hypothèse: Amplicon V3-V4 = ${AMPLICON_LENGTH}bp"
echo ""

if [ "$R1_MEDIAN" -lt 150 ] || [ "$R2_MEDIAN" -lt 150 ]; then
    echo "⚠️  ALERTE CRITIQUE: Reads TRÈS COURTS après Trimmomatic!"
    echo ""
    echo "Vos reads sont < 150bp → IMPOSSIBLE de couvrir un amplicon de 460bp"
    echo ""
    echo "CAUSES POSSIBLES:"
    echo "  1. Trimmomatic trop agressif (MINLEN trop élevé?)"
    echo "  2. Données de mauvaise qualité"
    echo "  3. Mauvais primers utilisés"
    echo ""
    echo "SOLUTION:"
    echo "  → Refaire Trimmomatic avec MINLEN:50 au lieu de MINLEN:100"
    echo "  → Contacter la société de séquençage"
    echo ""
elif [ "$R1_MEDIAN" -lt 250 ] || [ "$R2_MEDIAN" -lt 250 ]; then
    echo "⚠️  ATTENTION: Reads courts après Trimmomatic"
    echo ""
    echo "Vos reads médians:"
    echo "  R1: ${R1_MEDIAN}bp"
    echo "  R2: ${R2_MEDIAN}bp"
    echo ""
    echo "Paramètres DADA2 recommandés:"
    echo "  --p-trunc-len-f $((R1_MEDIAN - 20))"
    echo "  --p-trunc-len-r $((R2_MEDIAN - 20))"
    echo ""
    if [ "$OVERLAP" -lt 40 ]; then
        echo "⚠️  Overlap calculé: ${OVERLAP}bp (TRÈS FAIBLE!)"
        echo ""
        echo "PROBLÈME:"
        echo "  → Overlap insuffisant pour merger correctement"
        echo "  → Vous allez perdre beaucoup de reads"
        echo ""
        echo "SOLUTIONS:"
        echo "  1. Refaire Trimmomatic avec MINLEN plus bas"
        echo "  2. Utiliser --p-trunc-len 0 et accepter le faible taux de merge"
        echo "  3. Vérifier que les bons primers ont été utilisés"
        echo ""
    fi
else
    echo "✓ Reads de longueur acceptable"
    echo ""
    echo "Paramètres DADA2 OPTIMAUX:"
    echo ""
    echo "  qiime dada2 denoise-paired \\"
    echo "    --p-trunc-len-f $TRUNC_F \\"
    echo "    --p-trunc-len-r $TRUNC_R \\"
    echo "    --p-max-ee-f 3 \\"
    echo "    --p-max-ee-r 3 \\"
    echo "    --p-trunc-q 2 \\"
    echo "    --p-min-overlap 12"
    echo ""
    echo "  → Overlap théorique: ${OVERLAP}bp"
    echo ""
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DIAGNOSTIC TERMINÉ"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "PROCHAINES ACTIONS:"
echo ""
echo "1. Si reads < 250bp:"
echo "   → Refaire Trimmomatic avec MINLEN:50"
echo "   → Ou utiliser les données BRUTES (pas Trimmomatic)"
echo ""
echo "2. Si reads 250-300bp:"
echo "   → Relancer DADA2 avec les paramètres ci-dessus"
echo ""
echo "3. Si toujours problème:"
echo "   → Envoyer les screenshots de demux.qzv"
echo "   → Vérifier avec la société de séquençage"
echo ""
echo "=========================================================================================================="
