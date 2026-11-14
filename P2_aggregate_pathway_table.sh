#!/bin/bash
# ==============================================================
# P2_aggregate_pathway_table.sh
# 전체 샘플의 P1 결과(long format)를 하나로 합쳐서
# "샘플 x pathway" wide table(count/pct 양쪽)을 만듭니다.
#
# 사용법: bash P2_aggregate_pathway_table.sh
# ※ 이 스크립트는 eggnog_analysis/ 폴더 안에서 실행됩니다.
#    sample_list.txt가 현재 폴더에 있거나 심볼릭 링크로 연결되어 있어야 합니다.
# ==============================================================
set -euo pipefail

# sample_list.txt: 현재 폴더 우선, 없으면 상위 폴더
if [ -f "sample_list.txt" ]; then
    SAMPLE_LIST="sample_list.txt"
elif [ -f "../sample_list.txt" ]; then
    SAMPLE_LIST="../sample_list.txt"
else
    echo "❌ [P2] sample_list.txt를 현재 폴더 또는 상위 폴더에서 찾을 수 없습니다."
    exit 1
fi

NAMES_FILE="kegg_pathway_names_clean.tsv"
OUT_WIDE_COUNT="pathway_abundance_wide_count.tsv"
OUT_WIDE_PCT="pathway_abundance_wide_pct.tsv"
OUT_NAMES="pathway_names_reference.tsv"

if [ ! -s "$NAMES_FILE" ]; then
    echo "[안내] ${NAMES_FILE} 없음 → KEGG REST API에서 자동 다운로드 중..."
    if curl -s --fail https://rest.kegg.jp/list/pathway -o "$NAMES_FILE" && [ -s "$NAMES_FILE" ]; then
        echo "[안내] ${NAMES_FILE} 다운로드 완료"
    else
        rm -f "$NAMES_FILE"
        echo "❌ ${NAMES_FILE} 자동 다운로드 실패 (인터넷 연결 확인)"
        echo "   수동: curl -s https://rest.kegg.jp/list/pathway > ${NAMES_FILE}"
        exit 1
    fi
fi

# ---- 1) P1 안 돌린 샘플 확인 및 실행 ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P1_SCRIPT=""
for cand in "${SCRIPT_DIR}/P1_extract_pathway_counts.sh" "../P1_extract_pathway_counts.sh"; do
    [ -f "$cand" ] && P1_SCRIPT="$cand" && break
done
if [ -z "$P1_SCRIPT" ]; then
    echo "❌ P1_extract_pathway_counts.sh 를 찾을 수 없습니다."
    exit 1
fi

echo "[P2] 각 샘플 pathway 집계 확인 중..."
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    bash "$P1_SCRIPT" "$SAMPLE" || echo "⚠️  [P2] ${SAMPLE} P1 실패 → 결과 없이 진행"
done < "$SAMPLE_LIST"

# ---- 2) long-format 합치기 ----
echo "[P2] 전체 결과 합치는 중..."
COMBINED="all_samples_pathway_long.tsv"
echo -e "sample\tpathway_id\tcount\tpct" > "$COMBINED"
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    F="${SAMPLE}_pathway_counts.tsv"
    [ -s "$F" ] && tail -n +2 "$F" >> "$COMBINED"
done < "$SAMPLE_LIST"

N_ROWS=$(( $(wc -l < "$COMBINED") - 1 ))
echo "[P2] 합쳐진 행 수: ${N_ROWS}"
if [ "$N_ROWS" -eq 0 ]; then
    echo "❌ [P2] 유효한 pathway 데이터가 없습니다. P1 결과를 확인하세요."
    exit 1
fi

# ---- 3) long → wide 변환 ----
python3 << 'PYEOF'
import pandas as pd, sys, os

df = pd.read_csv("all_samples_pathway_long.tsv", sep="\t")

wide_count = df.pivot_table(index="sample", columns="pathway_id", values="count", fill_value=0)
wide_count.to_csv("pathway_abundance_wide_count.tsv", sep="\t")
print(f"[P2] count wide table: {wide_count.shape[0]}샘플 × {wide_count.shape[1]}개 pathway")

wide_pct = df.pivot_table(index="sample", columns="pathway_id", values="pct", fill_value=0)
wide_pct.to_csv("pathway_abundance_wide_pct.tsv", sep="\t")
print(f"[P2] pct wide table:   {wide_pct.shape[0]}샘플 × {wide_pct.shape[1]}개 pathway")
PYEOF

# ---- 4) pathway 이름 참고 파일 ----
python3 << 'PYEOF'
import pandas as pd

wide = pd.read_csv("pathway_abundance_wide_count.tsv", sep="\t", index_col=0)
names = pd.read_csv("kegg_pathway_names_clean.tsv", sep="\t", header=None,
                    names=["pathway_id", "pathway_name"])
ref = pd.DataFrame({"pathway_id": wide.columns})
ref = ref.merge(names, on="pathway_id", how="left")
ref.to_csv("pathway_names_reference.tsv", sep="\t", index=False)
print(f"[P2] 이름 참고파일: {len(ref)}개 pathway")
PYEOF

echo ""
echo "=== [P2] 완료 ==="
echo "  pathway_abundance_wide_count.tsv : ANCOM-BC2 통계 입력용 (raw count)"
echo "  pathway_abundance_wide_pct.tsv   : 시각화/기술통계용 (%)"
echo "  pathway_names_reference.tsv      : pathway ID ↔ 이름"
