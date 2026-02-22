#!/bin/bash

# ========== User Input ==========
REF="seq.fasta"
READ1="file1"
READ2="file2"
SNPEFF_GENOME="NZ_CP045110.1"

# ========== Filenames ==========
SAM="out.sam"
BAM="out.bam"
SORT_BAM="out_sort.bam"
PILEUP="out.pileup"
VCF_GZ="out.vcf.gz"
VCF_INDEX="${VCF_GZ}.csi"
FILTERED_VCF="out_filtered.recode.vcf"
SNPEFF_OUT="out_ann.vcf"
VCFTOOLS_LOG="vcftools_error.log"

# ========== Pipeline Steps ==========

echo "📌 [1/9] Indexing reference..."
bwa index "$REF"

echo "📌 [2/9] Aligning reads with BWA..."
bwa mem "$REF" "$READ1" "$READ2" > "$SAM"

echo "📌 [3/9] Converting SAM to BAM..."
samtools view -S -b "$SAM" > "$BAM"

echo "📌 [4/9] Sorting BAM file..."
samtools sort "$BAM" -o "$SORT_BAM"

echo "📌 [5/9] Indexing BAM..."
samtools index "$SORT_BAM"

echo "📌 [6/9] Variant calling with BCFtools..."
bcftools mpileup -f "$REF" "$SORT_BAM" | bcftools call -mv -Oz -o "$VCF_GZ"

echo "📌 [7/9] Indexing VCF.gz..."
bcftools index "$VCF_GZ"

echo "📌 [8/9] Filtering variants with VCFtools..."
vcftools --gzvcf "$VCF_GZ" --recode --recode-INFO-all --out out_filtered --minQ 30 2> "$VCFTOOLS_LOG"

echo "📌 [9/9] Annotating with SnpEff..."
java -jar snpEff/snpEff.jar "$SNPEFF_GENOME" "$FILTERED_VCF" > "$SNPEFF_OUT"

echo "✅ Pipeline completed. Final output: $SNPEFF_OUT"

