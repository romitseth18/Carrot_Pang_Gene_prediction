# Carrot_PanGenome
Genome Sequencing and Pangenome Analysis 

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


**Step  2: 02_rename_gff_ID_gene_extraction.sh**

## Workflow Diagram
The following diagram summarizes the complete two-step pangenome annotation workflow, including ID normalization, gene extraction, and QC:

```mermaid
[flowchart TD

    %% STYLE
    classDef step fill:#f4f4f4,stroke:#333,stroke-width:1px,border-radius:4px;
    classDef tool fill:#e8f0ff,stroke:#3366cc,stroke-width:1px,border-radius:4px;
    classDef out fill:#eaffea,stroke:#228b22,stroke-width:1px,border-radius:4px;

    %% INPUTS
    A[Softmasked Genomes (*.softmasked.fasta)]:::step
    B[RNA-seq R1/R2 reads]:::step
    C[SwissProt proteins]:::step
    D[DCv3 genome + GFF3]:::step
    E[TSEBRA config]:::step

    %% STEP 1: ANNOTATION PIPELINE
    A --> STAR1[Build STAR Index]:::tool
    B --> STAR2[STAR Alignment → sorted BAM]:::tool
    STAR1 --> BRAKER_ET[BRAKER ET (RNA-seq guided)]:::tool
    STAR2 --> BRAKER_ET

    C --> BRAKER_EP[BRAKER EP (protein homology)]:::tool

    BRAKER_ET --> TSEBRA[TSEBRA Model Selection]:::tool
    BRAKER_EP --> TSEBRA

    D --> LIFTOFF[Liftoff: Reference-guided Projection]:::tool

    TSEBRA --> MERGE[AGAT Merge (TSEBRA + Liftoff)]:::tool
    LIFTOFF --> MERGE

    MERGE --> POLISH[AGAT Fixes (codons, consistency)]:::tool
    POLISH --> CHR[Chromosome filtering (Chr1–9)]:::tool

    CHR --> OUT1[Final annotation → *.gff]:::out
    CHR --> OUT2[Stats → *.tsv]:::out

    %% STEP 2: RENAMING + SEQUENCE EXTRACTION
    OUT1 --> CLEAN[Normalize + convert to standard GFF]:::tool
    CLEAN --> RENAME[Python: Filter masked genes + Rename IDs sequentially]:::tool

    RENAME --> VALIDATE[Python: Validate GFF–FASTA consistency]:::tool

    RENAME --> EXTRACT_GENES[AGAT Extract gene sequences]:::tool
    RENAME --> EXTRACT_CDS[AGAT Extract CDS sequences]:::tool
    RENAME --> EXTRACT_PROT[AGAT Extract Protein sequences]:::tool

    EXTRACT_GENES --> QC1[Compute QC metrics\n(Counts + Avg length)]:::tool
    EXTRACT_CDS --> QC1
    EXTRACT_PROT --> QC1

    QC1 --> OUT3[QC Summary CSV]:::out
    EXTRACT_GENES --> OUT4[*.fna (genes)]:::out
    EXTRACT_CDS --> OUT5[*.cds.fna]:::out
    EXTRACT_PROT --> OUT6[*.faa]:::out


