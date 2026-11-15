#!/bin/bash
# ==============================================================
# P1b_ko_abundance_table.sh
# MicrobiomeAnalyst(웹, Shotgun Data Profiling)의
# "A gene abundance table" 업로드용으로, KEGG_ko(개별 유전자) 단위
# 샘플 x KO count 표를 만듭니다. (P1/P2는 pathway 단위라 여기엔 안 맞음)
#
# 사용법: bash P1b_ko_abundance_table.sh
# 결과: ko_abundance_table.tsv  (MicrobiomeAnalyst 업로드용)
#       microbiomeanalyst_metadata.tsv (그룹 정보, 업로드용)
# ==============================================================
set -euo pipefail

SAMPLE_LIST="sample_list.txt"
GROUP_FILE="sample_group_mapping.csv"
LONG_OUT="all_samples_ko_long.tsv"

echo -e "sample\tko_id\tcount" > "$LONG_OUT"

while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    ANNOT="${SAMPLE}_emapper.emapper.annotations"
    if [ ! -s "$ANNOT" ]; then
        echo "⚠️  ${SAMPLE}: ${ANNOT} 없음 -> 스킵"
        continue
    fi

    KO_COL=$(awk -F'\t' '/^#query/{for(i=1;i<=NF;i++) if($i=="KEGG_ko") print i}' "$ANNOT")

    grep -v '^#' "$ANNOT" \
      | cut -f"$KO_COL" \
      | grep -v '^-$' \
      | tr ',' '\n' \
      | sed 's/ko://g' \
      | sort | uniq -c | sed 's/^ *//' \
      | awk -v s="$SAMPLE" '{print s"\t"$2"\t"$1}' >> "$LONG_OUT"

    echo "[P1b] ${SAMPLE} 완료"
done < "$SAMPLE_LIST"

# long -> wide (샘플 x KO)
python3 << 'PYEOF'
import pandas as pd
df = pd.read_csv("all_samples_ko_long.tsv", sep="\t")
wide = df.pivot_table(index="sample", columns="ko_id", values="count", fill_value=0).astype(int)
# MicrobiomeAnalyst는 보통 "행=유전자, 열=샘플" 형태를 기대하므로 전치
wide.T.to_csv("ko_abundance_table.tsv", sep="\t", index_label="#NAME")
print(f"[P1b] KO abundance table: {wide.shape[1]}개 KO x {wide.shape[0]}샘플")
PYEOF

# metadata 파일도 MicrobiomeAnalyst 형식에 맞게 별도 생성
awk -F',' 'NR==1{print "#SampleID\t"$2} NR>1{print $1"\t"$2}' "$GROUP_FILE" > microbiomeanalyst_metadata.tsv

echo ""
echo "=== [P1b] 완료 ==="
echo "  ko_abundance_table.tsv       : Abundance file로 업로드"
echo "  microbiomeanalyst_metadata.tsv : Metadata file로 업로드"
