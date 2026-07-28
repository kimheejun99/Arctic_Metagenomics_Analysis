#!/bin/bash
# ==============================================================
# Step1_organize_samples.sh (타입별 폴더 버전)
#
# 현재 디렉토리의 flat한 파일들을 종류별로 정리:
#   {sample}.kreport, {sample}.kraken.out           -> ./kraken/
#   {sample}_annotated_parsed.tsv, _annotated.tsv,
#   _category_summary.tsv, _count.txt, _output.tsv  -> ./vfdb/
#   {sample}_args_output/                            -> ./arg/
#   fastq.gz는 이동하지 않고 현재 위치에 그대로 둠
#
# 사용법:
#   bash Step1_organize_samples.sh sample_list.txt
# ==============================================================
set -uo pipefail

SAMPLE_LIST="${1:-sample_list.txt}"

if [ ! -f "$SAMPLE_LIST" ]; then
    echo "오류: $SAMPLE_LIST 파일이 없습니다."
    exit 1
fi

mkdir -p ./kraken ./vfdb ./arg

while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    [[ "$SAMPLE" == "Run" ]] && continue

    echo "[Step1] ${SAMPLE} 정리 중..."
    moved=0

    # --- kraken ---
    for f in "${SAMPLE}.kreport" "${SAMPLE}.kraken.out"; do
        if [ -e "$f" ]; then
            mv "$f" ./kraken/
            moved=$((moved+1))
        fi
    done

    # --- vfdb ---
    for suffix in "_annotated_parsed.tsv" "_annotated.tsv" "_category_summary.tsv" "_count.txt" "_output.tsv"; do
        f="${SAMPLE}${suffix}"
        if [ -e "$f" ]; then
            mv "$f" ./vfdb/
            moved=$((moved+1))
        fi
    done

    # --- arg (args_output 폴더 통째로) ---
    if [ -d "${SAMPLE}_args_output" ]; then
        mv "${SAMPLE}_args_output" ./arg/
        moved=$((moved+1))
    fi

    if [ "$moved" -eq 0 ]; then
        echo "  [경고] ${SAMPLE} 관련 파일을 하나도 찾지 못했습니다."
    else
        echo "  -> ${moved}개 항목 정리 완료 (kraken/, vfdb/, arg/)"
    fi
done < "$SAMPLE_LIST"

echo "[Step1] 전체 완료."
