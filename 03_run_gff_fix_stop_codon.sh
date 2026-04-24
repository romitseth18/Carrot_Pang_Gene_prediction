#!/usr/bin/env bash
set -euo pipefail

# ============================================
# CONFIG
# ============================================
GFF_LIST="gff.txt"
FASTA_DIR="/home/wolftech/rseth2/data/00_final_GFFs_fasta"
SCRIPT="gff_rebuild_phase_correct.py"

OUTDIR="rebuild_results"
LOGDIR="$OUTDIR/logs"

mkdir -p "$OUTDIR" "$LOGDIR"

# ============================================
# LOOP THROUGH PAIRS
# ============================================
tail -n +2 "$GFF_LIST" | while read -r gff1 gff2; do

    prefix=$(basename "$gff2" .gff)
    fasta="$FASTA_DIR/${prefix}.fasta"

    echo "======================================="
    echo "[INFO] Processing: $prefix"
    echo "======================================="

    # ---- sanity checks ----
    if [[ ! -f "$gff1" ]]; then
        echo "[WARN] Missing gff1: $gff1 → skipping"
        continue
    fi

    if [[ ! -f "$gff2" ]]; then
        echo "[WARN] Missing gff2: $gff2 → skipping"
        continue
    fi

    if [[ ! -f "$fasta" ]]; then
        echo "[WARN] Missing FASTA: $fasta → skipping"
        continue
    fi

    # ---- step 1: ORF-safe rebuild ----
    raw_out="$OUTDIR/${prefix}.raw.gff"

    echo "[STEP1] Rebuilding annotation..."
    python3 "$SCRIPT" "$gff1" "$gff2" "$fasta" "$raw_out" \
        > "$LOGDIR/${prefix}.log" 2>&1

    # ---- quick sanity check ----
    gene_count=$(grep -c $'\tgene\t' "$raw_out" || true)

    if [[ "$gene_count" -lt 1000 ]]; then
        echo "[ERROR] Too few genes ($gene_count). Skipping $prefix"
        continue
    fi

    # ---- step 2: AGAT normalize ----
    agat_out="$OUTDIR/${prefix}.agat.gff"

    echo "[STEP2] AGAT cleanup..."
    agat_convert_sp_gxf2gxf.pl \
        -g "$raw_out" \
        -o "$agat_out" \
        >> "$LOGDIR/${prefix}.log" 2>&1

    # ---- step 3: GT tidy/sort ----
    final_out="$OUTDIR/${prefix}.final.gff"

    echo "[STEP3] Sorting + tidying..."
	# Proper sort (genomic order)
	gt gff3 -sort -tidy -retainids "$agat_out" > "$final_out"


    echo "[DONE] $prefix completed"

done

echo "======================================="
echo "ALL DONE"
echo "Results in: $OUTDIR"
