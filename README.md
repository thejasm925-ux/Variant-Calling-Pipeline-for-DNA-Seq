Variant Calling Pipeline for DNA-Seq
An automated Bash-based suite for prokaryotic and targeted eukaryotic variant calling. This repository provides multiple workflows to transition from raw sequencing data (SRA or local FASTQ) to annotated VCF files.

## Features
End-to-End Automation: Handles data retrieval, quality control, mapping, and variant detection.

Reference Flexibility: Automatically fetches the latest RefSeq genomes via NCBI Entrez Direct.

Versatile Applications: Optimized for bacterial genomics, yet compatible with Human Exome Sequencing (WES) and Prokaryotic Transcriptomics.

High-Fidelity Filtering: Implements strict quality filtering (Q≥30) using VCFtools.

Functional Annotation: Integrated with SnpEff to predict biological impacts.

## Scripts Included
## Research Applications
Bacterial Genomics: Rapidly detecting SNPs in genes associated with antibiotic resistance.

Human Exome Analysis: Capable of processing targeted human coding regions (Exomes) using BWA-MEM alignment.

Prokaryotic Transcriptomics: Identifying expressed variants and RNA-editing sites in bacterial systems.

## System Requirements & Configuration
To ensure the pipelines run efficiently, the following environment is recommended:

1. Hardware Recommendations
Processor: Minimum 4 cores (8+ recommended for human data).

RAM: 16GB minimum; 32GB+ required for indexing large eukaryotic genomes.

Storage: 50GB–100GB+ free space for intermediate BAM and SRA files.

2. Software Dependencies
The following tools must be in your $PATH:

Data Retrieval: SRA Toolkit, Entrez Direct.

Read Processing: fastp, Trimmomatic, FastQC.

Alignment & Discovery: BWA, Samtools, BCFtools, VCFtools.

Annotation: SnpEff (Java 8+).

## Quick Start
Clone the repository:

Run the full automated pipeline:

Example: ./automated_full_variant_pipeline.sh SRR12345678 "Escherichia coli K12".

## License
This project is licensed under the MIT License.
