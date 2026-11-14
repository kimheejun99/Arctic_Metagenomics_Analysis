#!/bin/bash
# ==============================================================
# V5_run_kraken_single.sh
# T3_run_kraken.sh(62개 배치 전용)의 핵심 로직을,
# run_pipeline.sh 루프에서 SRA ID 하나씩 호출할 수 있도록 재구성한 버전입니다.
# conda env: ArcticVF
# ==============================================================
set -uo pipefail

SAMPLE=$1

if [ -z "$SAMPLE" ]; then
    echo "사용법: $0 <SRA_RUN_ID>"
    exit 1
fi

DB="${KRAKEN2_DB:-/mnt/program/db/plusPFP_16}"   # 환경변수 KRAKEN2_DB로 override 가능
OUT_DIR="."
LOG_DIR="./logs"
THREADS="${KRAKEN2_THREADS:-16}"
CONFIDENCE="${KRAKEN2_CONFIDENCE:-0.1}"           # 이전 BEO 49샘플 파이프라인과 동일한 값

mkdir -p "$OUT_DIR" "$LOG_DIR"

R1="${SAMPLE}_1.fastq.gz"
R2="${SAMPLE}_2.fastq.gz"
SE="${SAMPLE}.fastq.gz"
KREPORT="${OUT_DIR}/${SAMPLE}.kreport"
KOUT="${OUT_DIR}/${SAMPLE}.kraken.out"

if [[ -s "$KREPORT" ]]; then
    echo "[V5] ${SAMPLE}: kreport가 이미 존재합니다 -> 스킵"
    exit 0
fi

echo "[V5] ${SAMPLE}: Kraken2 실행 중..."

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
    echo "[V5] ${SAMPLE}: FASTQ 파일 없음 -> 스킵" | tee -a "${LOG_DIR}/missing_fastq.log"
    exit 1
fi

echo "[V5] ${SAMPLE}: Kraken2 완료 -> ${KREPORT}"
