#!/bin/bash
# ==============================================================
# run_pipeline.sh
#
# Step0_make_googlemap.py  -> 구글맵 시각화
# Step1_organize_samples.sh -> 샘플별 폴더 정리
# Step2_vf_arg_plots.py    -> VF/ARG 바그래프 + 히트맵
# Step3_host_tracing.py    -> Kraken 결과와 매칭해 host phylum 추적
# Step4_merge_kreports.py  -> kreport 병합, ASV/taxonomy/metadata 3파일 생성
#
# 사용법:
#   bash run_pipeline.sh sample_coords.txt sample_list.txt metadata.txt
#
#   sample_coords.txt : Run, latitude, longitude (Step0용, 헤더 필수)
#   sample_list.txt   : 샘플 ID 목록, 한 줄에 하나 (선택사항 -- 없으면 sample_coords.txt의
#                        Run 컬럼에서 자동 생성됩니다. 특정 서브셋만 돌리고 싶을 때만 직접 준비하세요)
#   metadata.txt      : sample-id, group (탭 구분, 없으면 "" 로 지정, Step2/3/4가 group 없이 진행)
#
# 환경변수:
#   GOOGLE_MAPS_API_KEY : Step0 구글맵에 쓸 API 키 (없으면 지도가 로드되지 않는 placeholder로 생성됨)
#   TAXDUMP             : Step3/4용 NCBI taxdump(nodes.dmp/names.dmp) 경로
#   JOBS                : Step3 병렬 처리 프로세스 수
# ==============================================================
set -uo pipefail

COORDS_FILE="${1:-sample_coords.txt}"
SAMPLE_LIST="${2:-sample_list.txt}"
METADATA="${3:-metadata.txt}"
TAXDUMP="${TAXDUMP:-/mnt/program/db/plusPFP_16}"
JOBS="${JOBS:-6}"

LOG_FILE="pipeline_execution_log_$(date +%Y%m%d_%H%M%S).txt"
echo "파이프라인 실행 로그: $(date)" > "$LOG_FILE"
echo "coords=$COORDS_FILE, samples=$SAMPLE_LIST, metadata=$METADATA" >> "$LOG_FILE"
echo "----------------------------------------" >> "$LOG_FILE"

# sample_list.txt가 없으면 sample_coords.txt의 Run 컬럼에서 자동 생성
# (직접 서브셋을 지정하고 싶으면 sample_list.txt를 미리 만들어두면 그걸 우선 사용합니다)
if [ ! -f "$SAMPLE_LIST" ]; then
    if [ -f "$COORDS_FILE" ]; then
        echo "[안내] ${SAMPLE_LIST} 없음 -> ${COORDS_FILE}에서 샘플 목록 자동 추출" | tee -a "$LOG_FILE"
        tail -n +2 "$COORDS_FILE" | cut -f1 > "$SAMPLE_LIST"
    else
        echo "오류: ${SAMPLE_LIST}도 ${COORDS_FILE}도 없습니다. 최소 하나는 있어야 샘플 목록을 알 수 있습니다." | tee -a "$LOG_FILE"
        exit 1
    fi
fi

run_step() {
    local desc="$1"; shift
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ${desc} 실행 중..." | tee -a "$LOG_FILE"
    "$@" >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 오류: ${desc} 실패 (로그 확인: $LOG_FILE)" | tee -a "$LOG_FILE"
        exit 1
    fi
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ${desc} 완료" | tee -a "$LOG_FILE"
    echo "----------------------------------------" >> "$LOG_FILE"
}

# Step0: 구글맵 (좌표 파일 있을 때만)
if [ -f "$COORDS_FILE" ]; then
    run_step "Step0_make_googlemap" python3 Step0_make_googlemap.py "$COORDS_FILE" --api_key "${GOOGLE_MAPS_API_KEY:-YOUR_GOOGLE_MAPS_API_KEY}"
else
    echo "[안내] ${COORDS_FILE} 없음 -> Step0 건너뜀" | tee -a "$LOG_FILE"
fi

# Step1: 샘플별 폴더 정리
run_step "Step1_organize_samples" bash Step1_organize_samples.sh "$SAMPLE_LIST"

# Step2: VF/ARG 바그래프 + 히트맵
run_step "Step2_vf_arg_plots" python3 Step2_vf_arg_plots.py \
    --sample_list "$SAMPLE_LIST" --metadata "$METADATA" \
    --vfdb_dir ./vfdb --arg_dir ./arg --outdir ./plots

# Step3: host tracing
run_step "Step3_host_tracing" python3 Step3_host_tracing.py \
    --sample_list "$SAMPLE_LIST" --metadata "$METADATA" \
    --kraken_dir ./kraken --vfdb_dir ./vfdb --arg_dir ./arg \
    --taxdump "$TAXDUMP" --jobs "$JOBS" --outdir ./host_tracing

# Step4: kreport 병합 -> 3개 txt
run_step "Step4_merge_kreports" python3 Step4_merge_kreports.py \
    --sample_list "$SAMPLE_LIST" --metadata "$METADATA" \
    --kraken_dir ./kraken --taxdump "$TAXDUMP" --outdir ./merged_tables

echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 전체 파이프라인 완료" | tee -a "$LOG_FILE"
