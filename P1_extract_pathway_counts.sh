#!/bin/bash
# ==============================================================
# P1_extract_pathway_counts.sh
# 샘플 1개의 *_emapper.emapper.annotations에서:
#   - KEGG_Pathway 전체 추출 (head로 안 자름)
#   - "Global and overview maps" 총괄 카테고리 제외
#   - 샘플 내 상대비율(%) 계산 (샘플마다 read 수가 달라서 raw count 비교 불가)
# 출력: {SAMPLE}_pathway_counts.tsv (long format: pathway_id, count, pct)
#
# 사용법: bash P1_extract_pathway_counts.sh <SAMPLE>
# ==============================================================
set -euo pipefail

SAMPLE=$1
IN_ANNOT="${SAMPLE}_emapper.emapper.annotations"
OUT_TSV="${SAMPLE}_pathway_counts.tsv"

if [[ -s "$OUT_TSV" ]]; then
    echo "[P1] ${SAMPLE}: ${OUT_TSV} 이미 존재 -> 스킵"
    exit 0
fi

if [ ! -s "$IN_ANNOT" ]; then
    echo "❌ [P1] ${IN_ANNOT} 가 없습니다."
    exit 1
fi

# KEGG_Pathway 컬럼 번호를 헤더 이름으로 자동 탐색 (하드코딩 안 함)
COL=$(awk -F'\t' '/^#query/{for(i=1;i<=NF;i++) if($i=="KEGG_Pathway") print i}' "$IN_ANNOT")
if [ -z "$COL" ]; then
    echo "❌ [P1] ${SAMPLE}: KEGG_Pathway 컬럼을 못 찾았습니다."
    exit 1
fi

# Global and overview maps (KEGG Mapper Reconstruct 화면에서 확인된 공식 목록)
EXCLUDE="01100|01110|01120|01130|01200|01210|01212|01220|01230|01232|01240|01250|01310|01320"

TMP="${SAMPLE}_pathway_raw.tmp"
grep -v '^#' "$IN_ANNOT" \
  | cut -f"$COL" \
  | grep -v '^-$' \
  | tr ',' '\n' \
  | grep '^map' \
  | egrep -v "map(${EXCLUDE})" \
  | sort | uniq -c | sed 's/^ *//' > "$TMP"

TOTAL=$(awk '{sum+=$1} END{print sum}' "$TMP")
if [ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ]; then
    echo "⚠️  [P1] ${SAMPLE}: 매칭된 pathway가 없습니다."
    TOTAL=1  # 0으로 나누기 방지
fi

echo -e "sample\tpathway_id\tcount\tpct" > "$OUT_TSV"
awk -v s="$SAMPLE" -v tot="$TOTAL" '{printf "%s\t%s\t%s\t%.4f\n", s, $2, $1, ($1/tot)*100}' "$TMP" >> "$OUT_TSV"

rm -f "$TMP"
N_PATHWAYS=$(($(wc -l < "$OUT_TSV") - 1))
echo "[P1] ${SAMPLE}: 완료. pathway ${N_PATHWAYS}개, 총 매칭 ${TOTAL}건 -> ${OUT_TSV}"
