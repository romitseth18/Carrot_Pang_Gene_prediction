#!/bin/bash
set -euo pipefail

# ========== USER CONFIGURATION ==========
THREADS=24
BRAKER_THREADS=24

# Paths
GENOME_DIR="./"                          # Folder with *.softmasked.fasta
FASTQ_DIR="/home/wolftech/rseth2/pang/RNA-Seq/02_Filtered_reads"                            # PG2_R1.fastq.gz, PG2_R2.fastq.gz
PROTEINS="/home/wolftech/rseth2/data/blastdb/swissprot/sp.fa"          # Add the path to swissprot file
BRAKER_SCRIPT="/home/wolftech/rseth2/tools/BRAKER/scripts/braker.pl"    # Modify the path to braker.pl
TSEBRA_SCRIPT="/home/wolftech/rseth2/tools/BRAKER/TSEBRA/bin/tsebra.py" # Modify the path to tsebra.py
REF_GFF="/home/wolftech/rseth2/tools/chromeister/bin/DCv3.gff3"    # Add the path to reference gff      
REF_GENOME="/home/wolftech/rseth2/tools/chromeister/bin/DCv3.fa"   # Add the path to your reference genome
tsebra_cfg="/home/wolftech/rseth2/tools/BRAKER/TSEBRA/config/tsebra.cfg"  # Path to tsebra.cfg

# Optional paths for BRAKER
export GENEMARK_PATH="/home/wolftech/rseth2/tools/BRAKER/gmes_linux_64"   # Modify the ath to genemark
export PROTHINT_PATH="/home/wolftech/rseth2/tools/BRAKER/gmes_linux_64/ProtHint/bin" # Modify the ath to ProtHint

# Increase open file limit (for STAR temporary BAM sorting)
ulimit -n 65536


# ========== LOOP THROUGH ALL GENOMES ==========
for GENOME in *.softmasked.fasta; do
    BASE=$(basename "$GENOME" .softmasked.fasta)   # e.g., PG4_H1
    PREFIX=${BASE%%_*}                              # e.g., PG4

    echo "========== PROCESSING: $BASE (RNA-seq: ${PREFIX}_R1/2) =========="

    # Locate FASTQ files for this prefix
    R1="${FASTQ_DIR}/${PREFIX}_R1.flt.fastq.gz"
    R2="${FASTQ_DIR}/${PREFIX}_R2.flt.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "❌ RNA-seq files not found for $PREFIX — skipping $BASE"
        continue
    fi

    # Output directories
    STAR_INDEX="${BASE}_STARindex"
    STAR_OUT="${BASE}_STARout"
    BRAKER_ET_DIR="braker_${BASE}_ET"
    BRAKER_EP_DIR="braker_${BASE}_EP"
    TSEBRA_OUT="${BASE}_tsebra.gtf"
    TSEBRA_GFF="${BASE}_tsebra.gff"
    TSEBRA_mRNA="${BASE}_tsebra_mRNA.gff"
    GEMOMA_OUT="${BASE}_GeMoMa"

    # STAR: Index genome
    mkdir -p "$STAR_INDEX"
    if [[ ! -f "${STAR_INDEX}/Genome" ]]; then
        echo "Building STAR genome index for $BASE..."
    STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX" \
         --genomeFastaFiles "$GENOME" \
         --runThreadN "$THREADS" \
         --genomeSAindexNbases 13
    else
        echo "STAR index for $BASE already exists. Skipping."
    fi

    # STAR: Align RNA-seq reads
    mkdir -p "$STAR_OUT"
    if [[ ! -f "${STAR_OUT}/Aligned.sortedByCoord.out.bam" ]]; then
    echo "Running STAR alignment for $BASE..."  
    STAR --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${STAR_OUT}/" \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterType BySJout \
         --outSAMstrandField intronMotif \
         --outFilterMultimapNmax 1 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --sjdbOverhang 100
        
     else
         echo "STAR alignment BAM already exists for $BASE. Skipping STAR."
     fi
    
    # Step 3: Index BAM
    if [[ ! -f "${STAR_OUT}/Aligned.sortedByCoord.out.bam.bai" ]]; then
         echo "Indexing BAM for $BASE..."
    samtools index -@ "$THREADS" "${STAR_OUT}/Aligned.sortedByCoord.out.bam"
    else
         echo "BAM index for $BASE already exists. Skipping indexing."
    fi

    # BRAKER: ET mode
    perl "$BRAKER_SCRIPT" \
        --genome="$GENOME" \
        --bam="${STAR_OUT}/Aligned.sortedByCoord.out.bam" \
        --species="${BASE}_IC" \
        --threads="$BRAKER_THREADS" \
        --gff3 \
        --softmasking \
        --workingdir="$BRAKER_ET_DIR" \
        --useexisting

    # BRAKER: EP mode
    perl "$BRAKER_SCRIPT" \
        --genome="$GENOME" \
        --prot_seq="$PROTEINS" \
        --species="$BASE" \
        --threads="$BRAKER_THREADS" \
        --gff3 \
        --softmasking \
        --workingdir="$BRAKER_EP_DIR" \
        --useexisting

    # TSEBRA: Merge predictions
    python3 "$TSEBRA_SCRIPT" \
        --gtf "${BRAKER_EP_DIR}/braker.gtf" \
        --gtf "${BRAKER_ET_DIR}/braker.gtf" \
        -e "${BRAKER_EP_DIR}/hintsfile.gff" \
        -e "${BRAKER_ET_DIR}/hintsfile.gff" \
        -c "$tsebra_cfg" \
        -o "$TSEBRA_OUT"

    echo "✅ Done: $BASE → TSEBRA output: $TSEBRA_OUT"
    echo
    
    ##Covert gtf to gff
    
    agat_convert_sp_gxf2gxf.pl --gff "$TSEBRA_OUT" -o "$TSEBRA_GFF"

    ##Rename transcripts to mRNA
    sed 's/transcript/mRNA/g' "$TSEBRA_GFF" > "$TSEBRA_mRNA"
    rm -f "$TSEBRA_GFF"

     # Liftoff
    liftoff "$GENOME" "$REF_GENOME" \
        -g "$REF_GFF" \
        -o "${BASE}_liftoff.gff" \
        -p "$THREADS" \
         2> "${BASE}_liftoff.log"

    echo "✅ Completed: $BASE (TSEBRA + GeMoMa + Liftoff outputs ready)"
 
    # Organize outputs into folder
    OUTDIR="${BASE}"
    mkdir -p "$OUTDIR"

    mv "$STAR_INDEX" "$OUTDIR/"
    mv "$STAR_OUT" "$OUTDIR/"
    mv "$BRAKER_ET_DIR" "$OUTDIR/"
    mv "$BRAKER_EP_DIR" "$OUTDIR/"
    mv "$TSEBRA_OUT" "$OUTDIR/"
    mv "$TSEBRA_mRNA" "$OUTDIR/"
    mv "${BASE}_liftoff.gff" "$OUTDIR/"
    mv "${BASE}_liftoff.log" "$OUTDIR/"
    
    echo "📁 All outputs for $BASE moved to folder: $OUTDIR"
 
 #agat merge all gff and remove duplicates   
   agat_sp_merge_annotations.pl \
      --gff "${OUTDIR}/$TSEBRA_mRNA" \
      --gff "${OUTDIR}/${BASE}_liftoff.gff" \
      --out "${BASE}.gff"
      
#add start and stop codon in genomes
    agat_sp_add_start_and_stop.pl --gff "${BASE}.gff" --fasta "$GENOME" --out "${BASE}_ss.gff"

 #Step 1: Filter GFFs for Chr1–9
    awk -v pre="$BASE" '$1 ~ "^"pre"_Chr[0-9]" || $0 ~ /^#/ { print }' "${BASE}_ss.gff" > "${BASE}_ss_chr.gff" 
  
  # Step 3: Run AGAT statistics
     agat_sp_statistics.pl --gff "${BASE}_ss.gff" -f "$GENOME" -o "${BASE}_ss.txt" -p
     agat_sp_statistics.pl --gff "${BASE}_ss_chr.gff" -f "$GENOME" -o "${BASE}_ss_chr.txt" -p

  # Step 4: Convert stats to TSV
        awk '{val=$NF; $NF=""; sub(/[ \t]+$/, ""); print $0 "\t" val}' "${BASE}_ss.txt" > "${BASE}_ss.tsv"
        awk '{val=$NF; $NF=""; sub(/[ \t]+$/, ""); print $0 "\t" val}' "${BASE}_ss_chr.txt" > "${BASE}_ss_chr.tsv"
  
  # Step 5: Clean up
        rm -f "${BASE}_ss.txt"
        rm -f "${BASE}_ss_chr.txt"

#Move the output file to output dir
    mv "${BASE}.gff" "$OUTDIR/"
    mv "${BASE}_ss.gff" "$OUTDIR/"
    mv "${BASE}_ss_chr.gff" "$OUTDIR/"
    mv "${BASE}_ss.tsv" "$OUTDIR/"
    mv "${BASE}_ss_chr.tsv" "$OUTDIR/"
 
   echo "✅ Annotation stats for $BASE"     
   echo  
  
done
