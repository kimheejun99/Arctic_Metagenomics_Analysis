#!/bin/bash
# ==============================================================
# 01_run_kraken2.sh
# 62개 Arctic soil 샘플 Kraken2 종 동정 배치 실행
# conda env: ArcticVF
# ==============================================================
set -euo pipefail

DB="/mnt/program/db/plusPFP_16"          # Kraken2 DB 경로
FASTQ_DIR="."                             # 현재 디렉토리에 fastq.gz가 바로 있음
OUT_DIR="./kraken2_out"
LOG_DIR="./logs"
METADATA="./ids_62.txt"                   # sample-id 목록 (헤더 1줄 + 62개)
THREADS=16
CONFIDENCE=0.1                            # 이전 BEO 49샘플 파이프라인과 동일한 값으로 맞춰주세요

mkdir -p "$OUT_DIR" "$LOG_DIR"

SAMPLES=$(tail -n +2 "$METADATA" | cut -f1)
TOTAL=$(echo "$SAMPLES" | wc -l)
echo "총 ${TOTAL}개 샘플 처리 예정"

i=0
for SAMPLE in $SAMPLES; do
    i=$((i+1))
    R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
    SE="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
    KREPORT="${OUT_DIR}/${SAMPLE}.kreport"
    KOUT="${OUT_DIR}/${SAMPLE}.kraken.out"

    if [[ -s "$KREPORT" ]]; then
        echo "[$i/$TOTAL] ${SAMPLE}: 이미 완료됨 -> 스킵"
        continue
    fi

    echo "[$i/$TOTAL] ${SAMPLE}: Kraken2 실행 중..."

    if [[ -f "$R1" && -f "$R2" ]]; then
        # paired-end
        kraken2 \
            --db "$DB" --threads "$THREADS" --confidence "$CONFIDENCE" \
            --paired --gzip-compressed \
            --report "$KREPORT" --output "$KOUT" \
            "$R1" "$R2" \
            > "${LOG_DIR}/${SAMPLE}.log" 2>&1
    elif [[ -f "$SE" ]]; then
        # single-end
        kraken2 \
            --db "$DB" --threads "$THREADS" --confidence "$CONFIDENCE" \
            --gzip-compressed \
            --report "$KREPORT" --output "$KOUT" \
            "$SE" \
            > "${LOG_DIR}/${SAMPLE}.log" 2>&1
    else
        echo "[$i/$TOTAL] ${SAMPLE}: FASTQ 파일 없음 -> 스킵" | tee -a "${LOG_DIR}/missing_fastq.log"
        continue
    fi

    echo "[$i/$TOTAL] ${SAMPLE}: 완료"
done

echo "전체 Kraken2 실행 완료. kreport 파일: ${OUT_DIR}/*.kreport"
