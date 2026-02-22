#!/bin/bash

# Exit on error
set -e

# ========== USAGE ==========
# ./variant_pipeline.sh SRA_ID reference.fasta snpeff_genome_id output_prefix
# Example:
# ./variant_pipeline.sh SRR1234567 seq.fasta NZ_CP139575.1 ecoli

# ========== INPUT ==========
SRA_ID=$1
REF=$2
SNPEFF_GENOME=$3
OUT_PREFIX=$4

# ========== DIRECTORIES ==========
RAW_DIR="raw_data"
TRIM_DIR="trimmed_data"
QC_DIR="fastqc_reports"
RES_DIR="results"

mkdir -p $RAW_DIR $TRIM_DIR $QC_DIR $RES_DIR
cd $RES_DIR

# ========== STEP 1: Prefetch SRA ==========
echo "🔄 Prefetching SRA data..."
prefetch $SRA_ID --output-directory ../$RAW_DIR

# ========== STEP 2: Extract FASTQ ==========
echo "📦 Extracting FASTQ files..."
fasterq-dump ../$RAW_DIR/$SRA_ID -O ../$RAW_DIR --split-files

# ========== STEP 3: Run FastQC ==========
echo "🔬 Running FastQC..."
fastqc ../$RAW_DIR/${SRA_ID}_1.fastq -o ../$QC_DIR
fastqc ../$RAW_DIR/${SRA_ID}_2.fastq -o ../$QC_DIR

# ========== STEP 4: fastp ==========
echo "✨ Running fastp (trimming and filtering)..."
fastp -i ../$RAW_DIR/${SRA_ID}_1.fastq -I ../$RAW_DIR/${SRA_ID}_2.fastq \
      -o ../$TRIM_DIR/${SRA_ID}_1.trimmed.fq -O ../$TRIM_DIR/${SRA_ID}_2.trimmed.fq \
      -h ../$QC_DIR/${SRA_ID}_fastp.html -j ../$QC_DIR/${SRA_ID}_fastp.json

# ========== STEP 5: BWA Index ==========
echo "📚 Indexing reference genome..."
bwa index ../$REF

# ========== STEP 6: BWA MEM ==========
echo "🔗 Aligning reads..."
bwa mem ../$REF ../$TRIM_DIR/${SRA_ID}_1.trimmed.fq ../$TRIM_DIR/${SRA_ID}_2.trimmed.fq > ${OUT_PREFIX}.sam

# ========== STEP 7: Convert SAM → BAM ==========
echo "🔄 Converting SAM to BAM..."
samtools view -Sb ${OUT_PREFIX}.sam > ${OUT_PREFIX}.bam

# ========== STEP 8: Sort BAM ==========
echo "📏 Sorting BAM..."
samtools sort ${OUT_PREFIX}.bam -o ${OUT_PREFIX}_sort.bam

# ========== STEP 9: Index BAM ==========
echo "🗂️ Indexing BAM..."
samtools index ${OUT_PREFIX}_sort.bam

# ========== STEP 10: Variant Calling ==========
echo "🔍 Calling variants with bcftools..."
bcftools mpileup -f ../$REF ${OUT_PREFIX}_sort.bam | \
    bcftools call -mv -Oz -o ${OUT_PREFIX}.vcf.gz

# ========== STEP 11: Index VCF ==========
echo "📦 Indexing VCF..."
bcftools index ${OUT_PREFIX}.vcf.gz

# ========== STEP 12: VCF Filtering ==========
echo "⚗️ Filtering variants (min Q = 30)..."
vcftools --gzvcf ${OUT_PREFIX}.vcf.gz \
         --recode --recode-INFO-all --out ${OUT_PREFIX}_filtered \
         2> ${OUT_PREFIX}_vcftools.log

# ========== STEP 13: Annotation with SnpEff ==========
echo "🧬 Annotating with SnpEff..."
java -Xmx4g -jar ../snpEff/snpEff.jar $SNPEFF_GENOME \
     ${OUT_PREFIX}_filtered.recode.vcf > ${OUT_PREFIX}_ann.vcf

echo "✅ Pipeline complete. Results saved in ./results/"

