#!/usr/bin/env bash
set -euo pipefail

# ===== Settings =====
SRR="$1"            # csv with one column (SRR...)
THREADS=8                
ARIA_CONN=16             

# Check input
#if [[ ! -f "$CSV_FILE" ]]; then
#    echo "ERROR: CSV file not found: $CSV_FILE"
#    exit 1
#fi

#echo "Processing CSV: $CSV_FILE"
echo "Downloading $SRR..."

# ===== Main script =====
#while IFS=, read -r acc; do

    # Skip empty rows
#    [[ -z "${acc// }" ]] && continue

    # Skip header
#    if [[ "$acc" =~ ^(Run|run|Accession|accession)$ ]]; then
#        continue
#    fi

#    echo "=============================="
#    echo "Processing: $acc"
#    echo "=============================="

    # 1. Delete old files
    # rm -f "${acc}.sra" \
          # "${acc}_1.fastq.gz" \
          # "${acc}_2.fastq.gz" \
          # "${acc}.fastq.gz"

    # 2. Download .sra via Amazon S3
#    echo "[Download] ${acc}.sra"
    aria2c -s ${ARIA_CONN} -k 1M --continue=true \
           --split=${ARIA_CONN} \
           --file-allocation=none \
           "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}" \
           -d . \
           -o "${SRR}.sra"

    if [[ ! -s "${SRR}.sra" ]]; then
        echo "ERROR: download failed for $SRR"
        continue
    fi

    # 3. Fastq connversion
    echo "[fasterq-dump] ${SRR}"
    fasterq-dump --split-3 -e ${THREADS} "${SRR}.sra"

    # 4. Compression
    echo "[gzip]"
    pigz -p ${THREADS} "${SRR}"*.fastq

    echo "Done: $SRR"
    echo

#done < "$CSV_FILE"

echo "All finished."
