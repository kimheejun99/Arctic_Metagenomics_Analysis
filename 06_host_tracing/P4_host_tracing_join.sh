#!/bin/bash
# ==============================================================
# P4_host_tracing_join.sh
# eggNOG-mapper 결과(ReadID+KO)와 기존 Kraken2 결과(ReadID+taxid)를
# ReadID 기준으로 join하여 "이 read는 taxon X이며 KO Y를 가진다"는
# host-tracing 결과를 만듭니다.
#
# eggNOG-mapper의 #query 컬럼(예: SRR3984960.66_3_149_-)에서
# FragGeneScan이 붙인 접미사(_유전자번호_길이_strand)를 제거하면
# 원본 ReadID(SRR3984960.66)가 복원되며, 이게 Kraken2 결과의
# ReadID와 동일한 형식입니다.
#
# Kraken2 파일(.kraken)이 수GB로 크기 때문에, awk 연관배열로
# 필요한 ReadID만 1회 스캔으로 찾아내는 방식을 씁니다(전체를
# 메모리에 다 올리지 않고, 필요한 lookup만 올림).
#
# 사용법: bash P4_host_tracing_join.sh <SAMPLE>
# 결과: {SAMPLE}_host_tracing.tsv (read_id, taxid, KO, pathway_ids)
# ==============================================================
set -euo pipefail

SAMPLE=$1
ANNOT="${SAMPLE}_emapper.emapper.annotations"
OUT="${SAMPLE}_host_tracing.tsv"

if [[ -s "$OUT" ]]; then
    echo "[P4] ${SAMPLE}: ${OUT} 이미 존재 -> 스킵"
    exit 0
fi
if [ ! -s "$ANNOT" ]; then
    echo "❌ [P4] ${ANNOT} 가 없습니다."
    exit 1
fi

# .kraken / .kraken.out 파일 탐색 (현재 폴더 → 상위 폴더, 두 확장자 모두 시도)
KRAKEN=""
for suffix in ".kraken" ".kraken.out"; do
    for dir in "." ".." "../kraken" "./kraken"; do
        candidate="${dir}/${SAMPLE}${suffix}"
        if [ -s "$candidate" ]; then
            KRAKEN="$candidate"
            echo "[P4] ${SAMPLE}: kraken 파일 발견 → $KRAKEN"
            break 2
        fi
    done
done

if [ -z "$KRAKEN" ]; then
    echo "❌ [P4] ${SAMPLE}.kraken(.out) 을 현재/상위/kraken 폴더 어디에서도 못 찾았습니다."
    exit 1
fi

# KEGG_ko, KEGG_Pathway 컬럼 번호를 헤더 이름으로 자동 탐색
read KO_COL PW_COL <<< $(awk -F'\t' '
/^#query/{
    for(i=1;i<=NF;i++){
        if($i=="KEGG_ko") ko=i
        if($i=="KEGG_Pathway") pw=i
    }
    print ko, pw
}' "$ANNOT")

if [ -z "$KO_COL" ] || [ -z "$PW_COL" ]; then
    echo "❌ [P4] ${SAMPLE}: KEGG_ko 또는 KEGG_Pathway 컬럼을 못 찾았습니다."
    exit 1
fi

# ---- 1) eggNOG 결과에서 원본 ReadID 복원 + KO/Pathway 추출 ----
TMP_EGG="${SAMPLE}_egg_readid.tmp"
awk -F'\t' -v ko="$KO_COL" -v pw="$PW_COL" '
!/^#/ && $ko != "-" {
    readid = $1
    sub(/_[0-9]+_[0-9]+_[+-]$/, "", readid)
    gsub("ko:", "", $ko)
    print readid "\t" $ko "\t" $pw
}' "$ANNOT" > "$TMP_EGG"

N_EGG=$(wc -l < "$TMP_EGG")
echo "[P4] ${SAMPLE}: eggNOG 결과에서 ${N_EGG}개 read-KO 쌍 추출됨"

# ---- 2) Kraken2 결과와 join (awk 연관배열, 1회 스캔) ----
# Kraken2 표준 출력 형식 가정: classified\tread_id\ttaxid\tlength\tlca_mapping
# (형식이 다르면 아래 $2, $3 부분을 실제 컬럼 위치로 수정하세요)
awk -F'\t' '
NR==FNR { ko[$1]=$2; pw[$1]=$3; next }
($2 in ko) { print $2 "\t" $3 "\t" ko[$2] "\t" pw[$2] }
' "$TMP_EGG" "$KRAKEN" > "${OUT}.tmp"

echo -e "read_id\ttaxid\tKO\tKEGG_Pathway" > "$OUT"
cat "${OUT}.tmp" >> "$OUT"
rm -f "$TMP_EGG" "${OUT}.tmp"

N_JOINED=$(($(wc -l < "$OUT") - 1))
echo "[P4] ${SAMPLE}: 완료. taxon 정보까지 연결된 read ${N_JOINED}개 -> ${OUT}"
echo "     (참고: join 안 된 read는 Kraken2가 unclassified로 판정했거나,"
echo "      원본 ReadID 형식이 예상과 다를 수 있습니다)"
