#!/bin/bash
set -euo pipefail

# ==================== INPUTS ====================
SRR_ID=$1           # Example: SRR12345678
ORG="$2"            # Example: "Escherichia coli K12"
PREFIX="${SRR_ID}"  # Output file prefix

# ==================== CHECK ARGS ====================
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <SRR_ID> <Organism_Name>"
  exit 1
fi

# ==================== SETUP DIRS ====================
mkdir -p raw_data trimmed_data genome_data
cd raw_data

# ==================== DOWNLOAD SRA ====================
echo "[1/10] Downloading SRA: $SRR_ID"
fasterq-dump "$SRR_ID" --split-files -O .

READ1="${SRR_ID}_1.fastq"
READ2="${SRR_ID}_2.fastq"
cd ..

# ==================== FASTP TRIMMING ====================
echo "[2/10] Running fastp..."
fastp -i "raw_data/${READ1}" -I "raw_data/${READ2}" \
      -o "trimmed_data/${READ1}" -O "trimmed_data/${READ2}" \
      --detect_adapter_for_pe --thread 4 --html "${PREFIX}_fastp.html" --json "${PREFIX}_fastp.json"

# ==================== FETCH REFSEQ GENOME ====================
cd genome_data
echo "[3/10] Fetching reference genome for: $ORG"

ACC=$(esearch -db assembly -query "$ORG[Organism] AND latest[filter] AND complete genome[filter] AND \"reference genome\"[filter]" \
  | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | head -n 1)

if [[ -z "$ACC" ]]; then
  echo "❌ No matching RefSeq genome found for '$ORG'"
  exit 1
fi

BASE_NAME=$(basename "$ACC")
FA_GZ="${BASE_NAME}_genomic.fna.gz"
GB_GZ="${BASE_NAME}_genomic.gbff.gz"

echo "✅ Found genome: $BASE_NAME"
wget -q "${ACC}/${FA_GZ}" -O "$FA_GZ"
wget -q "${ACC}/${GB_GZ}" -O "$GB_GZ"
gunzip -f "$FA_GZ"
REF="${FA_GZ%.gz}"
cd ..

# ==================== VARIANT CALLING PIPELINE ====================
SAM="${PREFIX}.sam"
BAM="${PREFIX}.bam"
SORT_BAM="${PREFIX}_sorted.bam"
VCF="${PREFIX}.vcf.gz"
FILTERED_VCF="${PREFIX}_filtered.recode.vcf"

echo "[4/10] Indexing reference..."
bwa index "genome_data/$REF"

echo "[5/10] Aligning reads..."
bwa mem "genome_data/$REF" "trimmed_data/${READ1}" "trimmed_data/${READ2}" > "$SAM"

echo "[6/10] SAM -> BAM..."
samtools view -Sb "$SAM" > "$BAM"

echo "[7/10] Sorting BAM..."
samtools sort "$BAM" -o "$SORT_BAM"

echo "[8/10] Indexing BAM..."
samtools index "$SORT_BAM"

echo "[9/10] Variant calling..."
bcftools mpileup -f "genome_data/$REF" "$SORT_BAM" | bcftools call -mv -Oz -o "$VCF"
bcftools index "$VCF"

echo "[10/10] Filtering VCF..."
vcftools --gzvcf "$VCF" --recode --recode-INFO-all --out "${PREFIX}_filtered" --minQ 30 2> vcftools_error.log

echo "✅ All done!"
echo "  - Filtered VCF: $FILTERED_VCF"
echo "  - Reference: genome_data/$REF"
echo "  - GenBank: genome_data/$GB_GZ"

