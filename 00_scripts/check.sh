#!/bin/bash

# ========================================================================================================
# DADA2 EMERGENCY FIX - Paramètres ultra-permissifs pour problème de MERGE
# ========================================================================================================

QIIME_DIR="/nvme/bio/data_fungi/ExPLOI/05_QIIME2"
METADATA_FILE="/nvme/bio/data_fungi/ExPLOI/metadata_ExPLOI.tsv"
THREADS=16

eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

export PYTHONWARNINGS="ignore"

echo "=========================================================================================================="
echo "DADA2 EMERGENCY FIX - Ultra-permissif pour MERGE"
echo "=========================================================================================================="
echo ""
echo "PROBLÈME IDENTIFIÉ:"
echo "  → Le merge des reads R1/R2 échoue (0% à 87% merged)"
echo "  → Cause probable: Overlap insuffisant ou qualité en fin de read"
echo ""
echo "SOLUTION:"
echo "  → Tronquer davantage pour FORCER un overlap"
echo "  → Assouplir tous les paramètres de qualité"
echo "  → Réduire min-overlap à 8bp (minimum absolu)"
echo ""
echo "=========================================================================================================="
echo ""

# Sauvegarde
BACKUP_DIR="${QIIME_DIR}/core_backup_emergency_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"

for file in table.qza rep-seqs.qza dada2-stats.qza; do
    if [ -f "${QIIME_DIR}/core/${file}" ]; then
        cp "${QIIME_DIR}/core/${file}" "$BACKUP_DIR/"
    fi
done

echo "✓ Sauvegarde: $BACKUP_DIR"
echo ""

# ========================================================================================================
# TEST 1 : Troncature agressive pour forcer overlap
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 1 : Troncature agressive (forcer overlap de ~100bp)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Hypothèse: Amplicon V3-V4 = ~460bp, séquençage 2x300bp"
echo ""
echo "Paramètres:"
echo "  --p-trunc-len-f 240  (garder 240bp pour R1)"
echo "  --p-trunc-len-r 200  (garder 200bp pour R2)"
echo "  → Overlap théorique: 240 + 200 - 460 = -20bp... NON!"
echo ""
echo "CORRECTION: Si amplicon = 460bp"
echo "  --p-trunc-len-f 280  (garder 280bp)"
echo "  --p-trunc-len-r 240  (garder 240bp)"
echo "  → Overlap théorique: 280 + 240 - 460 = 60bp ✓"
echo ""
echo "Paramètres de qualité ULTRA-PERMISSIFS:"
echo "  --p-max-ee-f 5       (accepter 5 erreurs attendues)"
echo "  --p-max-ee-r 5"
echo "  --p-trunc-q 0        (ne pas tronquer sur qualité)"
echo "  --p-min-overlap 8    (minimum absolu pour merge)"
echo ""

read -p "Lancer DADA2 avec ces paramètres ? (y/n): " CONFIRM

if [[ ! "$CONFIRM" =~ ^[yY]$ ]]; then
    echo "Abandon."
    exit 0
fi

echo ""
echo "Lancement DADA2 (30-60 minutes)..."
echo "Début: $(date)"
echo ""

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${QIIME_DIR}/core/demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 5 \
  --p-max-ee-r 5 \
  --p-trunc-q 0 \
  --p-min-overlap 8 \
  --p-n-threads "$THREADS" \
  --p-chimera-method consensus \
  --o-table "${QIIME_DIR}/core/table.qza" \
  --o-representative-sequences "${QIIME_DIR}/core/rep-seqs.qza" \
  --o-denoising-stats "${QIIME_DIR}/core/dada2-stats.qza" \
  --verbose

DADA2_STATUS=$?

echo ""
echo "Fin: $(date)"
echo ""

if [ $DADA2_STATUS -eq 0 ]; then
    echo "✓ DADA2 terminé!"
    echo ""
    
    # Export stats pour analyse
    qiime tools export \
        --input-path "${QIIME_DIR}/core/dada2-stats.qza" \
        --output-path "${QIIME_DIR}/export/dada2_stats_test1"
    
    echo "Analyse des résultats..."
    
    python3 << 'EOFPYTHON'
import pandas as pd
import os

stats_file = "/nvme/bio/data_fungi/ExPLOI/05_QIIME2/export/dada2_stats_test1/stats.tsv"

if os.path.exists(stats_file):
    df = pd.read_csv(stats_file, sep="\t", comment="#")
    
    print("\n" + "=" * 80)
    print("RÉSULTATS PAR ÉCHANTILLON:")
    print("=" * 80)
    
    for idx, row in df.iterrows():
        sample = row['sample-id']
        input_reads = int(row['input'])
        merged = int(row['merged'])
        final = int(row['non-chimeric'])
        
        merge_pct = (merged / input_reads * 100) if input_reads > 0 else 0
        final_pct = (final / input_reads * 100) if input_reads > 0 else 0
        
        status = "✓" if merge_pct > 70 else "⚠️" if merge_pct > 30 else "❌"
        
        print(f"{status} {sample:20s} | Input: {input_reads:>6,} | Merged: {merged:>6,} ({merge_pct:>5.1f}%) | Final: {final:>6,} ({final_pct:>5.1f}%)")
    
    total_input = df['input'].astype(int).sum()
    total_merged = df['merged'].astype(int).sum()
    total_final = df['non-chimeric'].astype(int).sum()
    
    print("\n" + "=" * 80)
    print("TOTAUX:")
    print(f"  Input      : {total_input:>10,} reads")
    print(f"  Merged     : {total_merged:>10,} reads ({total_merged/total_input*100:.1f}%)")
    print(f"  Final      : {total_final:>10,} reads ({total_final/total_input*100:.1f}%)")
    print("=" * 80)

EOFPYTHON
    
    # Compter ASVs
    TEMP_EXPORT="${QIIME_DIR}/export/temp_count_$$"
    mkdir -p "$TEMP_EXPORT"
    
    qiime tools export \
      --input-path "${QIIME_DIR}/core/table.qza" \
      --output-path "$TEMP_EXPORT" 2>/dev/null
    
    NUM_ASVS=$(biom summarize-table -i "${TEMP_EXPORT}/feature-table.biom" 2>/dev/null | grep "Num observations:" | awk '{print $3}')
    
    rm -rf "$TEMP_EXPORT"
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "NOMBRE D'ASVs DÉTECTÉS: ${NUM_ASVS}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    if [ -n "$NUM_ASVS" ]; then
        RATIO=$(awk "BEGIN {printf \"%.1f\", ($NUM_ASVS / 4700) * 100}")
        
        echo "Comparaison:"
        echo "  Avant (trunc 0/0)  : 88 ASVs (1.9%)"
        echo "  Après (trunc 280/240): ${NUM_ASVS} ASVs (${RATIO}%)"
        echo ""
        
        if [ "$NUM_ASVS" -gt 1000 ]; then
            echo "✅ EXCELLENT! Amélioration majeure!"
        elif [ "$NUM_ASVS" -gt 300 ]; then
            echo "✓ Amélioration significative (mais encore insuffisant)"
            echo ""
            echo "PROCHAINE ACTION: Essayer troncature moins agressive"
            echo "  --p-trunc-len-f 290"
            echo "  --p-trunc-len-r 250"
        else
            echo "⚠️  Amélioration limitée"
            echo ""
            echo "PROBLÈME POSSIBLE:"
            echo "  → L'amplicon est peut-être TROP LONG pour 2x300bp"
            echo "  → Vérifier les primers utilisés"
            echo "  → Contacter la société de séquençage"
        fi
    fi
    
    echo ""
    echo "Vérifier les stats détaillées:"
    echo "  ${QIIME_DIR}/export/dada2_stats_test1/stats.tsv"
    echo ""
    
else
    echo "❌ ERROR: DADA2 a échoué!"
    exit 1
fi

echo "=========================================================================================================="
