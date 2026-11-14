#!/bin/bash
# ==============================================================
# run_all_49_eggnog.sh
# sample_list.txt에 있는 49샘플 전체에 대해:
#   1) G1 (서브샘플링+FragGeneScan, 이미 된 샘플은 자동 스킵)
#   2) E2 (eggNOG-mapper, 이미 된 샘플은 자동 스킵)
# 를 순차 실행합니다.
#
# 검증된 옵션(--sensmode fast --dmnd_iterate no)이 이미 반영돼 있습니다
# (--iterate 켰을 때 샘플당 31일 -> 끄니 3.3시간으로 확인됨).
#
# 중간에 끊겨도 다시 실행하면, 이미 끝난 샘플/단계는 자동으로 건너뛰고
# 이어서 진행됩니다 (G1/E2 각각의 skip 로직 재사용).
#
# 사용법: nohup bash run_all_49_eggnog.sh > run_all_49.log 2>&1 &
#         disown
# ==============================================================
set -uo pipefail   # -e는 뺌: 한 샘플 실패해도 나머지는 계속 진행

SAMPLE_LIST="sample_list.txt"
THREADS="${THREADS:-24}"
DB_DIR="${EGGNOG_DB_DIR:-/mnt/program/db/eggnog}"
EMAPPER_BIN="$(conda info --base)/envs/eggnog/bin/emapper.py"

if [ ! -s "$SAMPLE_LIST" ]; then
    echo "❌ ${SAMPLE_LIST} 를 찾을 수 없습니다."
    exit 1
fi

N_TOTAL=$(grep -c . "$SAMPLE_LIST")
echo "[전체 실행] ${SAMPLE_LIST}에서 ${N_TOTAL}개 샘플 확인됨"
echo "[전체 실행] 시작 시각: $(date)"

N_DONE=0
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    N_DONE=$((N_DONE + 1))
    echo ""
    echo "=========================================="
    echo "[${N_DONE}/${N_TOTAL}] ${SAMPLE} 시작 ($(date))"
    echo "=========================================="

    # ---- 1) G1: 서브샘플링 + FragGeneScan (이미 됐으면 스킵) ----
    bash G1_prep_ghostkoala.sh "$SAMPLE"
    if [ $? -ne 0 ]; then
        echo "⚠️  [${SAMPLE}] G1 실패 -> 이 샘플 건너뛰고 계속 진행"
        continue
    fi

    IN_FAA="${SAMPLE}_ghostkoala_upload.faa"
    OUT_ANNOT="${SAMPLE}_emapper.emapper.annotations"

    # ---- 2) E2: eggNOG-mapper (이미 됐으면 스킵) ----
    if [[ -s "$OUT_ANNOT" ]]; then
        echo "[${SAMPLE}] eggNOG-mapper 결과 이미 있음 -> 스킵"
        continue
    fi

    if [ ! -s "$IN_FAA" ]; then
        echo "⚠️  [${SAMPLE}] ${IN_FAA} 가 없어서 eggNOG-mapper 건너뜀"
        continue
    fi

    T0=$(date +%s)
    "$EMAPPER_BIN" \
        -i "$IN_FAA" \
        -o "${SAMPLE}_emapper" \
        --data_dir "$DB_DIR" \
        --cpu "$THREADS" \
        --sensmode fast \
        --dmnd_iterate no
    T1=$(date +%s)
    ELAPSED=$((T1 - T0))

    if [[ -s "$OUT_ANNOT" ]]; then
        echo "[${SAMPLE}] 완료 ($((ELAPSED/60))분 소요)"
    else
        echo "⚠️  [${SAMPLE}] eggNOG-mapper 결과 파일이 안 만들어짐 (실패 가능성)"
    fi
done < "$SAMPLE_LIST"

echo ""
echo "[전체 실행] 종료 시각: $(date)"
echo "[전체 실행] 완료. 아래로 최종 결과 파일 개수 확인:"
ls -1 *_emapper.emapper.annotations 2>/dev/null | wc -l
echo "(49개가 나와야 전체 완료)"
