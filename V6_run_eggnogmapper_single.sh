#!/bin/bash
# ==============================================================
# V6_run_eggnogmapper_single.sh
#
# V5_run_kraken_single.sh와 같은 패턴으로, eggNOG-mapper를 run_pipeline.sh
# (Stage 0) 루프에서 SRA ID 하나씩 호출할 수 있게 만든 버전입니다.
#
# 직접 emapper.py 옵션을 다시 작성하지 않고, 이미 검증된
# run_all_49_eggnog.sh(G1 서브샘플링+FragGeneScan -> emapper.py 순서가
# 이미 확인된 스크립트)를 그대로 재사용합니다. sample_list.txt를
# 이 샘플 하나짜리로 임시 교체해서 호출한 뒤 원상복구합니다.
#
# conda env: eggnog
# ==============================================================
set -uo pipefail

SAMPLE=$1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "$SAMPLE" ]; then
    echo "사용법: $0 <SRA_RUN_ID>"
    exit 1
fi

OUT_ANNOT="${SAMPLE}_emapper.emapper.annotations"
if [[ -s "$OUT_ANNOT" ]]; then
    echo "[V6] ${SAMPLE}: emapper 결과가 이미 존재합니다 -> 스킵"
    exit 0
fi

echo "[V6] ${SAMPLE}: eggNOG-mapper 실행 중 (run_all_49_eggnog.sh 재사용)..."

# sample_list.txt를 건드리지 않기 위해, 이 샘플 하나짜리 임시 목록으로
# run_all_49_eggnog.sh를 호출합니다 (해당 스크립트가 sample_list.txt를
# 인자 없이 고정으로 읽는 구조이므로, 안전하게 백업 후 임시 교체합니다).
BACKUP_MADE=0
if [ -f sample_list.txt ]; then
    cp sample_list.txt sample_list.txt.v6backup
    BACKUP_MADE=1
fi
echo "$SAMPLE" > sample_list.txt

source ~/miniforge3/etc/profile.d/conda.sh
conda activate eggnog
bash "${SCRIPT_DIR}/run_all_49_eggnog.sh"
RESULT=$?

# 원상복구
if [ "$BACKUP_MADE" -eq 1 ]; then
    mv sample_list.txt.v6backup sample_list.txt
fi

if [ $RESULT -ne 0 ] || [ ! -s "$OUT_ANNOT" ]; then
    echo "[V6] ${SAMPLE}: 오류 - annotation 파일이 생성되지 않음"
    exit 1
fi

echo "[V6] ${SAMPLE}: eggNOG-mapper 완료 -> ${OUT_ANNOT}"
