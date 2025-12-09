**Carrot Pangenome Annotation Pipeline**

**Step 1: 01_STAR_braker_lift_ss.sh**
This repository contains a automated workflow for structural gene annotation across multiple carrot genomes (haplotype-resolved assemblies). The pipeline integrates RNA-seq–guided prediction (BRAKER ET), protein-homology–guided prediction (BRAKER EP), evidence reconciliation with TSEBRA, reference-guided projection with Liftoff, and final merging, polishing, and diagnostics via AGAT.
The script is intended for large-scale pangenomic annotation where each genome is processed identically and independently.

**Overview**
For every *.softmasked.fasta genome in the working directory, the pipeline performs:

**RNA-seq alignment (STAR)**
Creates a STAR index → aligns paired-end reads → produces sorted, indexed BAM.

**Gene prediction (BRAKER)**
ET (RNA-seq guided)
EP (SwissProt-guided protein homology)
Model reconciliation (TSEBRA)
Combines ET and EP predictions into a final structural annotation (GTF → GFF3).
Reference annotation lift-over (Liftoff)
Projects DCv3 gene models onto each genome assembly.
Merging and refinement (AGAT)
Merge structural models (TSEBRA + Liftoff)
Add start/stop codons
Restrict to primary chromosomes (Chr1–9)
Generate QC statistics (TSV)

**Organization**
All outputs for each genome are placed in their own folder.

**Key Features**
Fully resumable at multiple steps (STAR index, BAM, BRAKER dirs).
Detects missing RNA-seq files and gracefully skips incomplete samples.
Produces polished .gff and .tsv QC outputs for downstream pangenome integration.
Automatic model standardization (mRNA naming, merging, filtering).
Scalable to tens or hundreds of assemblies.

**Inputs**
Genome assemblies
Softmasked FASTA files in the current directory:
PG4_H1.softmasked.fasta
PG4_H2.softmasked.fasta
…

**RNA-seq reads**
Filtered paired-end reads named as:
PG4_R1.flt.fastq.gz
PG4_R2.flt.fastq.gz

**Reference materials**
SwissProt protein set
Carrot reference genome (DCv3) + GFF
TSEBRA configuration file

**Outputs**
For each genome PG4_H1, a directory PG4_H1/ is created containing:

**Category	Files**
STAR	STAR index, sorted BAM, BAM index
BRAKER ET/EP	braker.gtf, hints, logs
TSEBRA	PG4_H1_tsebra.gtf, PG4_H1_tsebra_mRNA.gff
Liftoff	PG4_H1_liftoff.gff, logs
Merged models	PG4_H1.gff
Final refined models	PG4_H1_ss.gff, PG4_H1_ss_chr.gff
Stats	PG4_H1_ss.tsv, PG4_H1_ss_chr.tsv

Clean, standardized, chromosome-restricted GFFs are ready for downstream pangenome graph integration, orthogroup assignment, or expression analyses.

**Software Requirements**
The script expects the following software in PATH (or paths provided):
STAR
BRAKER ≥2.1.6
Genemark-ES / ProtHint
TSEBRA
Samtools
Liftoff
AGAT (all agat_sp_*.pl tools)
Perl, Python3, zlib, gzip

**Recommended system:**
Linux cluster / workstation
≥ 12 threads
ulimit ≥ 65536 (required by STAR)


**Step  2: 02_rename_gff_ID_gene_extraction.sh :: GFF ID Renaming & Gene/CDS/Protein Extraction**
This step performs post-processing and normalization of genome annotations generated in Step 1 (BRAKER ET, BRAKER EP, TSEBRA, Liftoff, AGAT merge).
The script standardizes GFF identifiers, validates compatibility with the genome FASTA, extracts gene/CDS/protein sequences, and produces a QC report for downstream pangenome integration.

**Overview**
02_rename_gff_ID_gene_extraction.sh automates the following tasks:

_**Prefix normalization**_
Ensures consistent identifiers across assemblies (e.g., PG4_H1, PG6_H2, B123).

**_GFF cleaning_**
Converts and normalizes annotation structure using AGAT.

**_Masked-gene filtering + sequential ID renaming_**
Uses Python scripts to:
Remove masked/invalid models

_**Assign strict sequential gene and transcript IDs**_
Validate feature boundaries against the genome FASTA
FASTA sanitation
Rewraps line lengths to 60 bp
Indexes the FASTA with samtools faidx

_**Sequence extraction**_
Produces:
<prefix>_genes.fna
<prefix>_cds.fna
<prefix>_proteins.faa

_**QC metrics summary**_
Computes:
Gene / CDS / protein counts

_**Average sequence lengths**_
Outputs: QC_summary.csv

_**Input Requirements**_
1. GFF files
Place all final annotations from Step 1 in the working directory:
*.gff

2. Genome FASTA files
Expected path structure (editable inside script):
/home/wolftech/rseth2/data/02_EC/01_EC_GFF_PROTEINS/fasta/EC_<GenomeName>.fasta  #Change the path according to data
The <GenomeName> is automatically parsed from the GFF filename.

3. Python support scripts
Required:
filter_and_rename_gff.py
validate_fasta_gff.py
Both must reside in your working dir. 


Output
All processed annotations and extracted sequences are written to:
01_rename_cds_pep_gff/

## Workflow Diagram
The following diagram summarizes the complete two-step pangenome annotation workflow, including ID normalization, gene extraction, and QC:

flowchart TD

    %% === INPUTS ===
    A[Softmasked Genomes] 
    B[RNA-seq Reads R1/R2]
    C[SwissProt Proteins]
    D[DCv3 Genome + GFF3]
    E[TSEBRA Config]

    %% === STEP 1 ANNOTATION ===
    A --> STAR1[STAR: Build Index]
    B --> STAR2[STAR: Align Reads → BAM]

    STAR1 --> BRAKER_ET[BRAKER ET - RNA-seq guided]
    STAR2 --> BRAKER_ET

    C --> BRAKER_EP[BRAKER EP - protein guided]

    BRAKER_ET --> TSEBRA[TSEBRA Model Selection]
    BRAKER_EP --> TSEBRA

    D --> LIFTOFF[Liftoff: Reference Projection]

    TSEBRA --> MERGE[AGAT: Merge Models]
    LIFTOFF --> MERGE

    MERGE --> POLISH[AGAT: Add start/stop codons]
    POLISH --> CHR[Filter Chr1–9]

    CHR --> OUT1[Final Annotation (*.gff)]
    CHR --> OUT2[Statistics (*.tsv)]

    %% === STEP 2 RENAME + EXTRACTION ===
    OUT1 --> CLEAN[Convert + Normalize GFF]

    CLEAN --> RENAME[Python: Rename Sequential Gene IDs]
    RENAME --> VALIDATE[Python: Validate GFF vs FASTA]

    RENAME --> EXTRACT_GENES[Extract Gene Sequences]
    RENAME --> EXTRACT_CDS[Extract CDS Sequences]
    RENAME --> EXTRACT_PROT[Extract Protein Sequences]

    EXTRACT_GENES --> QC[QC Metrics]
    EXTRACT_CDS --> QC
    EXTRACT_PROT --> QC

    QC --> OUT3[QC Summary CSV]
    EXTRACT_GENES --> OUT4[Genes FNA]
    EXTRACT_CDS --> OUT5[CDS FNA]
    EXTRACT_PROT --> OUT6[Proteins FAA]




