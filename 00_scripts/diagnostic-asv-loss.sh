#!/bin/bash

# ========================================================================================================
# DIAGNOSTIC COMPLET - Perte d'ASVs dans le pipeline
# ========================================================================================================
# Ce script analyse chaque étape pour identifier où les ASVs sont perdus

QIIME_DIR="/nvme/bio/data_fungi/ExPLOI/05_QIIME2"
EXPORT_DIR="${QIIME_DIR}/export"

# Activer QIIME2
eval "$(conda shell.bash hook)"
conda activate /scratch_vol0/fungi/envs/qiime2-amplicon-2024.10

export PYTHONWARNINGS="ignore"

echo "=========================================================================================================="
echo "DIAGNOSTIC : Analyse de la perte d'ASVs"
echo "=========================================================================================================="
echo ""
echo "Date: $(date)"
echo ""

# ========================================================================================================
# ÉTAPE 1 : STATISTIQUES DADA2
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 1 : Statistiques DADA2 (Débruitage)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -f "${QIIME_DIR}/core/dada2-stats.qza" ]; then
    # Export stats DADA2
    rm -rf "${EXPORT_DIR}/dada2_stats_diagnostic"
    qiime tools export \
        --input-path "${QIIME_DIR}/core/dada2-stats.qza" \
        --output-path "${EXPORT_DIR}/dada2_stats_diagnostic"
    
    echo "Analyse des stats DADA2..."
    echo ""
    
    python3 << 'EOFPYTHON'
import pandas as pd
import os

export_dir = "/nvme/bio/data_fungi/ExPLOI/05_QIIME2/export"
stats_file = os.path.join(export_dir, "dada2_stats_diagnostic", "stats.tsv")

if os.path.exists(stats_file):
    df = pd.read_csv(stats_file, sep="\t")
    
    print("Colonnes disponibles:")
    print(df.columns.tolist())
    print("")
    
    # Calculer statistiques par étape
    cols = df.columns.tolist()
    
    # Identifier colonnes de comptage
    count_cols = [col for col in cols if col not in ['sample-id', 'sample_id']]
    
    print("STATISTIQUES PAR ÉTAPE DADA2:")
    print("=" * 80)
    
    for col in count_cols:
        if df[col].dtype in ['int64', 'float64']:
            total = df[col].sum()
            mean = df[col].mean()
            min_val = df[col].min()
            max_val = df[col].max()
            print(f"\n{col}:")
            print(f"  Total    : {total:>12,.0f} reads")
            print(f"  Moyenne  : {mean:>12,.1f} reads/sample")
            print(f"  Min      : {min_val:>12,.0f} reads")
            print(f"  Max      : {max_val:>12,.0f} reads")
    
    # Calculer pertes par étape
    print("\n" + "=" * 80)
    print("PERTES PAR ÉTAPE:")
    print("=" * 80)
    
    if 'input' in cols and 'filtered' in cols:
        input_total = df['input'].sum()
        filtered_total = df['filtered'].sum()
        loss = input_total - filtered_total
        pct = (loss / input_total * 100) if input_total > 0 else 0
        print(f"\n1. Filtrage qualité:")
        print(f"   Input      : {input_total:>12,.0f} reads")
        print(f"   Filtré     : {filtered_total:>12,.0f} reads")
        print(f"   Perte      : {loss:>12,.0f} reads ({pct:.1f}%)")
    
    if 'filtered' in cols and 'denoised' in cols:
        filtered_total = df['filtered'].sum()
        denoised_total = df['denoised'].sum()
        loss = filtered_total - denoised_total
        pct = (loss / filtered_total * 100) if filtered_total > 0 else 0
        print(f"\n2. Débruitage:")
        print(f"   Filtré     : {filtered_total:>12,.0f} reads")
        print(f"   Débruité   : {denoised_total:>12,.0f} reads")
        print(f"   Perte      : {loss:>12,.0f} reads ({pct:.1f}%)")
    
    if 'denoised' in cols and 'merged' in cols:
        denoised_total = df['denoised'].sum()
        merged_total = df['merged'].sum()
        loss = denoised_total - merged_total
        pct = (loss / denoised_total * 100) if denoised_total > 0 else 0
        print(f"\n3. Merge des pairs:")
        print(f"   Débruité   : {denoised_total:>12,.0f} reads")
        print(f"   Mergé      : {merged_total:>12,.0f} reads")
        print(f"   Perte      : {loss:>12,.0f} reads ({pct:.1f}%)")
    
    if 'merged' in cols and 'non-chimeric' in cols:
        merged_total = df['merged'].sum()
        nonchimeric_total = df['non-chimeric'].sum()
        loss = merged_total - nonchimeric_total
        pct = (loss / merged_total * 100) if merged_total > 0 else 0
        print(f"\n4. Détection chimères:")
        print(f"   Mergé      : {merged_total:>12,.0f} reads")
        print(f"   Non-chimère: {nonchimeric_total:>12,.0f} reads")
        print(f"   Perte      : {loss:>12,.0f} reads ({pct:.1f}%)")
    
    print("\n" + "=" * 80)
    
    # Échantillons avec peu de reads
    if 'non-chimeric' in cols:
        low_samples = df[df['non-chimeric'] < 1000].sort_values('non-chimeric')
        if len(low_samples) > 0:
            print("\n⚠️  ÉCHANTILLONS AVEC < 1000 READS FINAUX:")
            print("-" * 80)
            for idx, row in low_samples.iterrows():
                sample_id = row.get('sample-id', row.get('sample_id', 'unknown'))
                final_reads = row['non-chimeric']
                input_reads = row.get('input', 0)
                retention = (final_reads / input_reads * 100) if input_reads > 0 else 0
                print(f"  {sample_id:20s} : {final_reads:>8,.0f} reads ({retention:.1f}% rétention)")
    
    print("\n")

else:
    print("ERROR: Fichier stats.tsv introuvable!")

EOFPYTHON

else
    echo "❌ ERROR: dada2-stats.qza introuvable!"
fi

echo ""

# ========================================================================================================
# ÉTAPE 2 : NOMBRE D'ASVs PAR TABLE
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 2 : Nombre d'ASVs à chaque étape"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

count_features() {
    local qza_file=$1
    local label=$2
    
    if [ -f "$qza_file" ]; then
        # Export temporaire
        local temp_dir="${EXPORT_DIR}/temp_count_$$"
        mkdir -p "$temp_dir"
        
        qiime tools export \
            --input-path "$qza_file" \
            --output-path "$temp_dir" 2>/dev/null
        
        if [ -f "${temp_dir}/feature-table.biom" ]; then
            # Compter features
            local num_features=$(biom summarize-table -i "${temp_dir}/feature-table.biom" 2>/dev/null | grep "Num observations:" | awk '{print $3}')
            local num_samples=$(biom summarize-table -i "${temp_dir}/feature-table.biom" 2>/dev/null | grep "Num samples:" | awk '{print $3}')
            
            echo "${label}:"
            echo "  ASVs       : ${num_features}"
            echo "  Échantillons: ${num_samples}"
            echo ""
        fi
        
        rm -rf "$temp_dir"
    else
        echo "${label}: ❌ Fichier introuvable"
        echo ""
    fi
}

count_features "${QIIME_DIR}/core/table.qza" "1. Table DADA2 (après débruitage)"
count_features "${QIIME_DIR}/core/table-decontam.qza" "2. Table après décontamination"
count_features "${QIIME_DIR}/core/table-final.qza" "3. Table finale (sans contrôle négatif)"

# ========================================================================================================
# ÉTAPE 3 : ANALYSE DE LA QUALITÉ DES READS (DEMUX)
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 3 : Qualité des reads (recommandations troncature)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -f "${QIIME_DIR}/visual/demux.qzv" ]; then
    echo "Pour analyser la qualité des reads:"
    echo "  1. Télécharger: ${QIIME_DIR}/visual/demux.qzv"
    echo "  2. Ouvrir sur: https://view.qiime2.org"
    echo "  3. Regarder la section 'Interactive Quality Plot'"
    echo ""
    echo "QUESTIONS À VÉRIFIER:"
    echo "  • À quelle position le score qualité tombe sous Q30 pour R1 ?"
    echo "  • À quelle position le score qualité tombe sous Q30 pour R2 ?"
    echo "  • Quelle est la longueur de vos reads ?"
    echo ""
    echo "RECOMMANDATIONS GÉNÉRALES (16S V3-V4):"
    echo "  • Si Illumina MiSeq 2x300bp:"
    echo "    --p-trunc-len-f 280  (garder 280bp pour R1)"
    echo "    --p-trunc-len-r 220  (garder 220bp pour R2)"
    echo ""
    echo "  • Si qualité reste bonne jusqu'à la fin:"
    echo "    --p-trunc-len-f 0"
    echo "    --p-trunc-len-r 0"
    echo ""
else
    echo "❌ demux.qzv introuvable!"
fi

echo ""

# ========================================================================================================
# ÉTAPE 4 : COMPARAISON AVEC SOCIÉTÉ DE SÉQUENÇAGE
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 4 : Comparaison avec société de séquençage"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -f "${EXPORT_DIR}/feature_table/feature-table.tsv" ]; then
    NUM_ASVS_PIPELINE=$(tail -n +3 "${EXPORT_DIR}/feature_table/feature-table.tsv" | wc -l)
    
    echo "Résultats:"
    echo "  • Votre pipeline       : ${NUM_ASVS_PIPELINE} ASVs"
    echo "  • Société séquençage   : 4700 ASVs (déclaré)"
    echo "  • Différence           : $(( 4700 - NUM_ASVS_PIPELINE )) ASVs manquants"
    echo "  • Ratio                : $(awk "BEGIN {printf \"%.1f\", ($NUM_ASVS_PIPELINE / 4700) * 100}")%"
    echo ""
    
    if [ "$NUM_ASVS_PIPELINE" -lt 500 ]; then
        echo "⚠️  ALERTE CRITIQUE: Vous avez perdu >90% des ASVs!"
        echo ""
        echo "CAUSES PROBABLES (par ordre de probabilité):"
        echo ""
        echo "1. ❌ TRONCATURE INCORRECTE (--p-trunc-len)"
        echo "   → La qualité des reads chute et DADA2 rejette tout"
        echo "   → Solution: Ajuster --p-trunc-len-f et --p-trunc-len-r"
        echo ""
        echo "2. ❌ MERGE IMPOSSIBLE"
        echo "   → Les reads R1 et R2 ne se chevauchent pas assez"
        echo "   → Solution: Ajuster --p-trunc-len ou --p-min-overlap"
        echo ""
        echo "3. ❌ DÉCONTAMINATION AGRESSIVE"
        echo "   → Le contrôle négatif contient beaucoup d'ASVs communs"
        echo "   → Solution: Vérifier neg-controls-summary.qzv"
        echo ""
        echo "4. ❌ PARAMÈTRES DADA2 TROP STRICTS"
        echo "   → max-ee trop bas, min-overlap trop haut"
        echo "   → Solution: Ajuster les paramètres (voir ci-dessous)"
        echo ""
    fi
fi

echo ""

# ========================================================================================================
# ÉTAPE 5 : RECOMMANDATIONS
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 5 : Recommandations pour corriger le pipeline"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "PARAMÈTRES DADA2 À TESTER:"
echo ""
echo "A) AJUSTER LA TRONCATURE (selon demux.qzv):"
echo "   qiime dada2 denoise-paired \\"
echo "     --p-trunc-len-f 280 \\"
echo "     --p-trunc-len-r 220 \\"
echo ""
echo "B) ASSOUPLIR LES FILTRES DE QUALITÉ:"
echo "   qiime dada2 denoise-paired \\"
echo "     --p-max-ee-f 3 \\"
echo "     --p-max-ee-r 3 \\"
echo "     --p-trunc-q 2 \\"
echo ""
echo "C) AJUSTER LE MERGE:"
echo "   qiime dada2 denoise-paired \\"
echo "     --p-min-overlap 12 \\"
echo ""
echo "D) DÉSACTIVER LA DÉTECTION DE CHIMÈRES (test):"
echo "   qiime dada2 denoise-paired \\"
echo "     --p-chimera-method none \\"
echo ""

echo ""

# ========================================================================================================
# ÉTAPE 6 : FICHIERS À VÉRIFIER
# ========================================================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "ÉTAPE 6 : Fichiers à vérifier manuellement"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "1. Qualité des reads bruts:"
echo "   ${QIIME_DIR}/visual/demux.qzv"
echo ""
echo "2. Stats DADA2 détaillées:"
echo "   ${QIIME_DIR}/visual/dada2-stats.qzv"
echo ""
echo "3. Résumé table finale:"
echo "   ${QIIME_DIR}/visual/table-final-summary.qzv"
echo ""
echo "4. Contrôle négatif (si existe):"
echo "   ${QIIME_DIR}/visual/neg-controls-summary.qzv"
echo ""

echo "Tous ces fichiers peuvent être ouverts sur: https://view.qiime2.org"
echo ""

echo "=========================================================================================================="
echo "DIAGNOSTIC TERMINÉ"
echo "=========================================================================================================="
echo ""
echo "PROCHAINES ACTIONS:"
echo "  1. Vérifier demux.qzv pour déterminer les bonnes valeurs de troncature"
echo "  2. Vérifier dada2-stats.qzv pour voir où les reads sont perdus"
echo "  3. Relancer DADA2 avec les paramètres ajustés"
echo ""
echo "=========================================================================================================="
