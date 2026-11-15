#!/bin/bash
# ==============================================================
# P5_aggregate_host_tracing.sh
# 전체 {SAMPLE}_host_tracing.tsv(P4 산출물)를 sample/group 정보와
# 함께 하나의 마스터 표로 합칩니다.
#
# 사용법: bash P5_aggregate_host_tracing.sh
# ※ eggnog_analysis/ 폴더 안에서 실행됩니다.
#    sample_list.txt / sample_group_mapping.csv 가 현재 폴더에 있어야 합니다.
# ==============================================================
set -euo pipefail

# sample_list.txt: 현재 폴더 우선, 없으면 상위 폴더
if [ -f "sample_list.txt" ]; then
    SAMPLE_LIST="sample_list.txt"
elif [ -f "../sample_list.txt" ]; then
    SAMPLE_LIST="../sample_list.txt"
else
    echo "❌ [P5] sample_list.txt를 현재 또는 상위 폴더에서 찾을 수 없습니다."
    exit 1
fi

GROUP_FILE="sample_group_mapping.csv"
OUT="host_tracing_master.tsv"

if [ ! -s "$GROUP_FILE" ]; then
    echo "❌ [P5] ${GROUP_FILE} 가 없습니다. Step7 make_group_mapping을 먼저 실행하세요."
    exit 1
fi

echo "[P5] sample_group_mapping.csv 헤더: $(head -1 "$GROUP_FILE")"

echo -e "read_id\tsample\tgroup\ttaxid\tKO\tKEGG_Pathway" > "$OUT"

N_MISSING=0
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    F="${SAMPLE}_host_tracing.tsv"
    if [ ! -s "$F" ]; then
        echo "⚠️  [P5] ${F} 없음 → 건너뜀 (P4를 먼저 실행하세요)"
        N_MISSING=$((N_MISSING + 1))
        continue
    fi
    GROUP=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $2}' "$GROUP_FILE")
    if [ -z "$GROUP" ]; then
        echo "⚠️  [P5] ${SAMPLE} 그룹 정보 없음 → UNKNOWN 처리"
        GROUP="UNKNOWN"
    fi
    tail -n +2 "$F" \
        | awk -F'\t' -v s="$SAMPLE" -v g="$GROUP" \
            '{print $1"\t"s"\t"g"\t"$2"\t"$3"\t"$4}' >> "$OUT"
done < "$SAMPLE_LIST"

N_TOTAL=$(( $(wc -l < "$OUT") - 1 ))
echo ""
echo "=== [P5] 완료 ==="
echo "  ${OUT}: 총 ${N_TOTAL}개 read-taxon-KO 레코드"
[ "$N_MISSING" -gt 0 ] && echo "  ⚠️  ${N_MISSING}개 샘플 누락 (P4 실행 확인 필요)"
