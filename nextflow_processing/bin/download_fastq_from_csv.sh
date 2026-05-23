#!/usr/bin/env bash
set -euo pipefail

# ===== НАСТРОЙКИ =====
SRR="$1"            # csv с одним столбцом (SRR...)
THREADS=8                # для fasterq-dump
ARIA_CONN=16             # соединений aria2c

# Проверка входного файла
#if [[ ! -f "$CSV_FILE" ]]; then
#    echo "ERROR: CSV file not found: $CSV_FILE"
#    exit 1
#fi

#echo "Processing CSV: $CSV_FILE"
echo "Downloading $SRR..."

# ===== ОСНОВНОЙ ЦИКЛ =====
#while IFS=, read -r acc; do

    # пропускаем пустые строки
#    [[ -z "${acc// }" ]] && continue

    # пропускаем header (если есть)
#    if [[ "$acc" =~ ^(Run|run|Accession|accession)$ ]]; then
#        continue
#    fi

#    echo "=============================="
#    echo "Processing: $acc"
#    echo "=============================="

    # 1. Удаляем старые файлы
    # rm -f "${acc}.sra" \
          # "${acc}_1.fastq.gz" \
          # "${acc}_2.fastq.gz" \
          # "${acc}.fastq.gz"

    # 2. Скачиваем .sra через Amazon S3
#    echo "[Download] ${acc}.sra"
    aria2c -s ${ARIA_CONN} -k 1M --continue=true \
           --split=${ARIA_CONN} \
           --file-allocation=none \
           "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}" \
           -d . \
           -o "${SRR}.sra"

    # Проверка скачивания
    if [[ ! -s "${SRR}.sra" ]]; then
        echo "ERROR: download failed for $SRR"
        continue
    fi

    # 3. Конвертация в FASTQ
    echo "[fasterq-dump] ${SRR}"
    fasterq-dump --split-3 -e ${THREADS} "${SRR}.sra"

    # 4. Сжатие
    echo "[gzip]"
    pigz -p ${THREADS} "${SRR}"*.fastq

    echo "Done: $SRR"
    echo

#done < "$CSV_FILE"

echo "All finished."
