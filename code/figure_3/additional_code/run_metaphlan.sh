#!/bin/bash

#SBATCH --account=open          # specify the account
#SBATCH --partition=open        # specify the partition
#SBATCH --job-name=metaphlan
#SBATCH --ntasks=4
#SBATCH --time=03:00:00
#SBATCH --mem=30GB

module load bowtie2

# Check if the required arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_fastq_R1> [<input_fastq_R2>]"
    echo "If single-end, provide only <input_fastq_R1>."
    exit 1
fi

# Input FASTQ files and database paths
INPUT_FASTQ_R1=$1
INPUT_FASTQ_R2=$2

# Extract the base name of the input FASTQ file (remove path and extension)
BASENAME=$(basename "$INPUT_FASTQ_R1" .fastq.gz)
BASENAME=${BASENAME%_R1}  # Remove _R1 if present

echo "${BASENAME}"

# Output directories
KNEADDATA_OUTPUT_DIR="/storage/home/kvm6065/scratch/kneaddata_output_${BASENAME}"
echo "${KNEADDATA_OUTPUT_DIR}"
mkdir -p "$KNEADDATA_OUTPUT_DIR"

kneaddata \
	--input1 "$INPUT_FASTQ_R1" \
	--input2 "$INPUT_FASTQ_R2" \
	--reference-db hdb/hg37dec_v0.1 \
	--trimmomatic Trimmomatic-0.39 \
	--output "$KNEADDATA_OUTPUT_DIR" \
	--threads 4

# Check if KneadData succeeded
if [ $? -ne 0 ]; then
    echo "KneadData preprocessing failed. Exiting."
    exit 1
fi

CLEANED_FASTQ="${KNEADDATA_OUTPUT_DIR}/${BASENAME}_R1_kneaddata_paired_1.fastq"
OUTPUT_FILE="${BASENAME}_metaphlan_profile.txt"
echo "$CLEANED_FASTQ"
echo "$OUTPUT_FILE"

# Run MetaPhlAn
metaphlan --force \
          --mpa3 \
	  --input_type fastq \
	  -t rel_ab_w_read_stats \
	  --index mpa_v30_CHOCOPhlAn_201901 \
	  --bowtie2db ./mpa_v30_CHOCOPhlAn_201901/ \
	  --nproc 4 \
	  $CLEANED_FASTQ \
	  -o "$OUTPUT_FILE"
