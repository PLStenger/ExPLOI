#!/bin/bash

# ========================================================================================================
# VÉRIFICATION LONGUEUR READS BRUTS
# ========================================================================================================

RAW_DATA_DIR="/nvme/bio/data_fungi/ExPLOI/01_raw_data"

echo "=========================================================================================================="
echo "VÉRIFICATION : Longueur des reads BRUTS (avant tout traitement)"
echo "=========================================================================================================="
echo ""

# Prendre le premier fichier R1 et R2
R1_FILE=$(find "$RAW_DATA_DIR" -name "*_R1_001.fastq.gz" | head -1)
R2_FILE=$(find "$RAW_DATA_DIR" -name "*_R2_001.fastq.gz" | head -1)

if [ -z "$R1_FILE" ]; then
    echo "❌ ERROR: Aucun fichier trouvé dans $RAW_DATA_DIR"
    exit 1
fi

echo "Fichier analysé: $(basename $R1_FILE)"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "LONGUEUR DES READS BRUTS (500 premiers reads)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "R1 (Forward) - Distribution des longueurs:"
zcat "$R1_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | tail -10

R1_MIN=$(zcat "$R1_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | head -1)
R1_MAX=$(zcat "$R1_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | tail -1)
R1_MODE=$(zcat "$R1_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')

echo ""
echo "R2 (Reverse) - Distribution des longueurs:"
zcat "$R2_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | tail -10

R2_MIN=$(zcat "$R2_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | head -1)
R2_MAX=$(zcat "$R2_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | tail -1)
R2_MODE=$(zcat "$R2_FILE" | head -2000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RÉSUMÉ"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "R1 (Forward):"
echo "  Min        : ${R1_MIN} bp"
echo "  Max        : ${R1_MAX} bp"
echo "  Mode       : ${R1_MODE} bp (longueur la plus fréquente)"
echo ""
echo "R2 (Reverse):"
echo "  Min        : ${R2_MIN} bp"
echo "  Max        : ${R2_MAX} bp"
echo "  Mode       : ${R2_MODE} bp (longueur la plus fréquente)"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DIAGNOSTIC"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ "$R1_MODE" -eq 251 ] || [ "$R1_MODE" -eq 301 ]; then
    echo "✓ SÉQUENÇAGE STANDARD DÉTECTÉ"
    echo ""
    if [ "$R1_MODE" -eq 251 ]; then
        echo "  Plateforme: Illumina MiSeq 2x250bp"
        echo ""
        echo "  ⚠️  PROBLÈME CRITIQUE:"
        echo "  → Amplicon V3-V4 = 460bp"
        echo "  → Reads 2x250bp = 500bp total"
        echo "  → Overlap théorique: 250 + 250 - 460 = 40bp"
        echo ""
        echo "  PARAMÈTRES DADA2 RECOMMANDÉS:"
        echo "    --p-trim-left-f 17"
        echo "    --p-trim-left-r 21"
        echo "    --p-trunc-len-f 0  (garder tout)"
        echo "    --p-trunc-len-r 0  (garder tout)"
        echo "    --p-min-overlap 12"
        echo ""
        echo "  → Après trim: (250-17) + (250-21) - 460 = 2bp overlap"
        echo "  → TRÈS LIMITE ! Taux de merge attendu: 40-60%"
        echo ""
    else
        echo "  Plateforme: Illumina MiSeq 2x300bp"
        echo ""
        echo "  ✓ PARFAIT pour amplicon V3-V4"
        echo "  → Overlap théorique: 300 + 300 - 460 = 140bp"
        echo ""
        echo "  PARAMÈTRES DADA2 RECOMMANDÉS:"
        echo "    --p-trim-left-f 17"
        echo "    --p-trim-left-r 21"
        echo "    --p-trunc-len-f 0"
        echo "    --p-trunc-len-r 0"
        echo ""
    fi
else
    echo "⚠️  LONGUEUR INHABITUELLE DÉTECTÉE: ${R1_MODE}bp"
    echo ""
    
    if [ "$R1_MODE" -lt 200 ]; then
        echo "  ❌ ALERTE: Reads TRÈS COURTS (< 200bp)"
        echo ""
        echo "  CAUSES POSSIBLES:"
        echo "    1. Mauvaise région amplifiée"
        echo "    2. Problème lors du séquençage"
        echo "    3. Mauvais fichiers fournis"
        echo ""
        echo "  IMPOSSIBLE de couvrir un amplicon V3-V4 de 460bp avec 2x${R1_MODE}bp"
        echo ""
        echo "  ACTIONS:"
        echo "    → Contacter URGENCE la société de séquençage"
        echo "    → Vérifier que ce sont les bons fichiers"
        echo "    → Demander quelle région 16S a été séquencée"
        echo ""
    else
        echo "  PARAMÈTRES DADA2 À TESTER:"
        echo "    --p-trim-left-f 17"
        echo "    --p-trim-left-r 21"
        echo "    --p-trunc-len-f $((R1_MODE - 20))"
        echo "    --p-trunc-len-r $((R2_MODE - 20))"
        echo "    --p-min-overlap 12"
        echo ""
        
        OVERLAP=$((R1_MODE - 17 + R2_MODE - 21 - 460))
        echo "  → Overlap estimé: ${OVERLAP}bp"
        
        if [ "$OVERLAP" -lt 20 ]; then
            echo "  → ⚠️  TRÈS FAIBLE! Beaucoup de pertes attendues"
        fi
    fi
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PROCHAINES ÉTAPES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ "$R1_MODE" -lt 200 ]; then
    echo "❌ STOP - Reads trop courts"
    echo "   → Contacter la société de séquençage"
    echo ""
elif [ "$R1_MODE" -eq 251 ]; then
    echo "1. Relancer DADA2 avec --p-trunc-len 0 (garder longueur maximale)"
    echo "2. Accepter un taux de merge faible (40-70%)"
    echo "3. Espérer ~500-2000 ASVs au lieu de 4700"
    echo ""
    echo "OU"
    echo ""
    echo "→ Demander à la société comment ils ont obtenu 4700 ASVs"
    echo "  (Peut-être ont-ils utilisé une méthode différente?)"
    echo ""
elif [ "$R1_MODE" -eq 301 ]; then
    echo "✓ Longueur OK pour V3-V4"
    echo ""
    echo "1. Relancer DADA2 avec:"
    echo "     --p-trim-left-f 17"
    echo "     --p-trim-left-r 21"
    echo "     --p-trunc-len-f 0"
    echo "     --p-trunc-len-r 0"
    echo ""
    echo "2. Devrait donner 3000-4500 ASVs"
    echo ""
fi

echo "=========================================================================================================="
