#!/bin/bash
set -euo pipefail
shopt -s nullglob

INPUT_DIR="./"
OUTPUT_DIR="$INPUT_DIR/01_rename_cds_pep_gff"
FASTA_DIR="/home/wolftech/rseth2/data/02_EC/01_EC_GFF_PROTEINS/fasta/"  #Change this to fasta file location
mkdir -p "$OUTPUT_DIR"

QC_FILE="$OUTPUT_DIR/QC_summary.csv"

# Initialize CSV header
echo "Prefix,Genes_Count,Genes_AvgLen,CDS_Count,CDS_AvgLen,Proteins_Count,Proteins_AvgLen" > "$QC_FILE"

for gff in "$INPUT_DIR"/*.gff; do
    base=$(basename "$gff" .gff)
    echo "Processing $gff ..."

    # Normalize prefix (e.g., PG4_H1)
    prefix=$(echo "$base" | sed -E 's/(PG[0-9]+)(_H[0-9]+).*/\1\2/; s/(B[0-9]+)/\1/')

    final_gff="$OUTPUT_DIR/${prefix}.gff"

    # Locate genome FASTA
    genome_name=$(echo "$base" | sed -E 's/(_cat|_ss.*)$//')
    fasta="$FASTA_DIR/EC_${genome_name}.fasta"

    if [[ ! -f "$fasta" ]]; then
        echo "⚠️ Genome FASTA $fasta not found, skipping $base"
        continue
    fi

    # Step 0: Pre-clean GFF
    echo "Pre-cleaning GFF..."
    agat_convert_sp_gxf2gxf.pl --gff "$gff" -o "$OUTPUT_DIR/${prefix}.clean.gff"

    # Step 1: Filter masked genes + sequential renaming
    echo "Filtering masked genes and renaming IDs..."
    python3 /home/wolftech/rseth2/rst/Scripts/Pangenome/sequential_gene_extraction/filter_and_rename_gff.py "$OUTPUT_DIR/${prefix}.clean.gff" "$prefix" "$fasta" "$final_gff"

    python3 /home/wolftech/rseth2/rst/Scripts/Pangenome/sequential_gene_extraction/validate_fasta_gff.py "$final_gff" "$fasta"

    rm -f "$OUTPUT_DIR/${prefix}.clean.gff"

    # Rewrap FASTA and reindex
    echo "Re-wrapping FASTA..."
    seqtk seq -l 60 "$fasta" > "${fasta%.fasta}.wrapped.fasta"
    mv "${fasta%.fasta}.wrapped.fasta" "$fasta"
    samtools faidx "$fasta"

    # Step 2: Extract sequences
    echo "Extracting gene sequences..."
    agat_sp_extract_sequences.pl --gff "$final_gff" --fasta "$fasta" --type gene --out "$OUTPUT_DIR/${prefix}_genes.fna"

    echo "Extracting CDS sequences..."
    agat_sp_extract_sequences.pl --gff "$final_gff" --fasta "$fasta" --type CDS --out "$OUTPUT_DIR/${prefix}_cds.fna"

    echo "Extracting protein sequences..."
    agat_sp_extract_sequences.pl --gff "$final_gff" --fasta "$fasta" --type CDS --protein --out "$OUTPUT_DIR/${prefix}_proteins.faa"

    # Step 3: QC metrics
    echo "Running QC checks..."
    genes_count=$(grep -c "^>" "$OUTPUT_DIR/${prefix}_genes.fna")
    genes_avg=$(awk 'NR%2==0 {sum+=length($0); n++} END {if(n>0) print int(sum/n); else print 0}' "$OUTPUT_DIR/${prefix}_genes.fna")

    cds_count=$(grep -c "^>" "$OUTPUT_DIR/${prefix}_cds.fna")
    cds_avg=$(awk 'NR%2==0 {sum+=length($0); n++} END {if(n>0) print int(sum/n); else print 0}' "$OUTPUT_DIR/${prefix}_cds.fna")

    prot_count=$(grep -c "^>" "$OUTPUT_DIR/${prefix}_proteins.faa")
    prot_avg=$(awk 'NR%2==0 {sum+=length($0); n++} END {if(n>0) print int(sum/n); else print 0}' "$OUTPUT_DIR/${prefix}_proteins.faa")

    echo "  - genes: $genes_count seqs, avg length = $genes_avg"
    echo "  - cds: $cds_count seqs, avg length = $cds_avg"
    echo "  - proteins: $prot_count seqs, avg length = $prot_avg"

    # Append results to CSV
    echo "$prefix,$genes_count,$genes_avg,$cds_count,$cds_avg,$prot_count,$prot_avg" >> "$QC_FILE"

    echo "✅ Finished $prefix"
done

echo "📊 QC summary saved to $QC_FILE"

