#!/bin/bash
# ==============================================================
# run_pipeline.sh  ─ Pipeline 1: 샘플별 원시 데이터 처리
#
# 역할:
#   V1(SRA 다운로드) → V5(Kraken2) → V2-1(fastq 카운트) →
#   V2-2(fastq→fasta) → V3(DIAMOND) → A2(ARGs-OAP) →
#   V6(eggNOG-mapper) → V4-1~3(VF 파싱) → V9(정리)
#   → result/ 생성 완료 후 Pipeline 2(run_pipeline_analysis.sh)로 자동 이어짐
#
# 사용법:
#   bash run_pipeline.sh biosample_sra_matches.txt
#
# 환경변수로 경로 override 가능:
#   VFDB_FASTA=/path/to/VFDB_setB_pro.fas
#   KRAKEN2_DB=/path/to/kraken2db
#   SAMPLE_COORDS=sample_coords.txt
#   SAMPLE_LIST_FILE=sample_list.txt
#   METADATA_FILE=metadata.txt
#   NEXT_PIPELINE=run_pipeline_analysis.sh  (Pipeline 2 스크립트 경로)
# ==============================================================
set -uo pipefail

# ──────────────────────────────────────────────
# 0) 경로 설정
# ──────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VFDB_FASTA="${VFDB_FASTA:-/mnt/program/db/VFDB/VFDB_setB_pro.fas}"
KRAKEN2_DB="${KRAKEN2_DB:-/mnt/program/db/plusPFP_16}"
SAMPLE_COORDS="${SAMPLE_COORDS:-sample_coords.txt}"
SAMPLE_LIST_FILE="${SAMPLE_LIST_FILE:-sample_list.txt}"
METADATA_FILE="${METADATA_FILE:-metadata.txt}"
NEXT_PIPELINE="${NEXT_PIPELINE:-${SCRIPT_DIR}/run_pipeline_analysis.sh}"

# ──────────────────────────────────────────────
# 1) 입력 파일 확인
# ──────────────────────────────────────────────
if [ $# -ne 1 ]; then
    echo "사용법: $0 biosample_sra_matches.txt"
    exit 1
fi
INPUT_FILE=$1

if [ ! -f "$INPUT_FILE" ]; then
    echo "오류: 입력 파일 '$INPUT_FILE'을 찾을 수 없습니다."
    exit 1
fi

# 필수 DB 존재 확인
if [ ! -f "$VFDB_FASTA" ]; then
    echo "오류: VFDB_FASTA 파일이 없습니다: $VFDB_FASTA"
    echo "  override: VFDB_FASTA=/your/path bash $0 ..."
    exit 1
fi
if [ ! -d "$KRAKEN2_DB" ]; then
    echo "오류: KRAKEN2_DB 디렉토리가 없습니다: $KRAKEN2_DB"
    echo "  override: KRAKEN2_DB=/your/path bash $0 ..."
    exit 1
fi

# 필수 스크립트 존재 확인 (실행 전 한 번에 검사)
PREFLIGHT_FAIL=0
for script in \
    V1_SRA_data_download.sh \
    V2-1_count_fastq.sh \
    V2-2_fq_to_fa_by_seqkit.sh \
    V3_run_diamond.sh \
    V4-1_blastx_tophit_annotator.py \
    V4-2_vf_annotation_parser.py \
    V4-3_categorize.py \
    V5_run_kraken_single.sh \
    V6_run_eggnogmapper_single.sh \
    A2_args_oap.sh \
    V9_clean.sh; do
    if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
        echo "[사전검사 실패] 필수 스크립트 없음: ${script}"
        PREFLIGHT_FAIL=1
    fi
done
if [ "$PREFLIGHT_FAIL" -eq 1 ]; then
    echo "위 스크립트를 이 파일과 같은 폴더에 두고 다시 실행하세요."
    exit 1
fi

# ──────────────────────────────────────────────
# 2) 로그 설정
# ──────────────────────────────────────────────
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline1_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Pipeline 1 시작: $(date)"
echo "입력 파일: $INPUT_FILE"
echo "VFDB_FASTA: $VFDB_FASTA"
echo "KRAKEN2_DB: $KRAKEN2_DB"
echo "로그 파일: $LOG_FILE"
echo "============================================================"

# ──────────────────────────────────────────────
# 3) 샘플별 처리
# ──────────────────────────────────────────────
FAIL_COUNT=0
SUCCESS_COUNT=0

tail -n +2 "$INPUT_FILE" | while read -r line; do
    SRA_RUN_ID=$(echo "$line" | awk '{print $NF}' | tr -d '\r ')

    if [ -z "$SRA_RUN_ID" ] || [ "$SRA_RUN_ID" = "No_SRA_IDs_found" ]; then
        echo "[$(date +%H:%M:%S)] SRA ID 없음, 건너뜀"
        continue
    fi

    echo ""
    echo "══════════════════════════════════════════"
    echo "[$(date +%H:%M:%S)] 처리 시작: $SRA_RUN_ID"
    echo "══════════════════════════════════════════"

    # 각 단계 실행 함수 (실패 시 해당 샘플만 skip)
    run_step() {
        local step_name="$1"; shift
        echo "[$(date +%H:%M:%S)] ${step_name} 실행 중..."
        if ! "$@"; then
            echo "[오류] ${step_name} 실패 → ${SRA_RUN_ID} 건너뜀"
            return 1
        fi
        echo "[$(date +%H:%M:%S)] ${step_name} 완료"
        return 0
    }

    (
        cd "$SCRIPT_DIR" || exit 1
        export KRAKEN2_DB
        export VFDB_FASTA

        run_step "V1_SRA_download"      bash "${SCRIPT_DIR}/V1_SRA_data_download.sh" "$SRA_RUN_ID"          || exit 1
        run_step "V5_kraken2"           bash "${SCRIPT_DIR}/V5_run_kraken_single.sh" "$SRA_RUN_ID" || exit 1
        run_step "V2-1_count_fastq"     bash "${SCRIPT_DIR}/V2-1_count_fastq.sh" "$SRA_RUN_ID"              || exit 1
        run_step "V2-2_fq_to_fa"        bash "${SCRIPT_DIR}/V2-2_fq_to_fa_by_seqkit.sh" "$SRA_RUN_ID"       || exit 1
        run_step "V3_diamond"           bash "${SCRIPT_DIR}/V3_run_diamond.sh" "$SRA_RUN_ID"                 || exit 1
        run_step "A2_args_oap"          bash "${SCRIPT_DIR}/A2_args_oap.sh" "$SRA_RUN_ID"                    || exit 1
        run_step "V6_eggnogmapper"      bash "${SCRIPT_DIR}/V6_run_eggnogmapper_single.sh" "$SRA_RUN_ID"     || exit 1
        run_step "V4-1_blastx_tophit"   python3 "${SCRIPT_DIR}/V4-1_blastx_tophit_annotator.py" "$SRA_RUN_ID" "$VFDB_FASTA" || exit 1
        run_step "V4-2_vf_parser"       python3 "${SCRIPT_DIR}/V4-2_vf_annotation_parser.py" "$SRA_RUN_ID"  || exit 1
        run_step "V4-3_categorize"      python3 "${SCRIPT_DIR}/V4-3_categorize.py" "$SRA_RUN_ID"             || exit 1
        run_step "V9_clean"             bash "${SCRIPT_DIR}/V9_clean.sh" "$SRA_RUN_ID"                       || exit 1
    )

    if [ $? -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] ✓ ${SRA_RUN_ID} 처리 완료"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "[$(date +%H:%M:%S)] ✗ ${SRA_RUN_ID} 실패 (일부 단계 오류)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    sleep 2
done

echo ""
echo "============================================================"
echo "Pipeline 1 완료: $(date)"
echo "  성공: ${SUCCESS_COUNT}개 / 실패: ${FAIL_COUNT}개"
echo "============================================================"

# ──────────────────────────────────────────────
# 4) Pipeline 2로 자동 이어 실행
# ──────────────────────────────────────────────
RESULT_DIR="${SCRIPT_DIR}/result"

if [ ! -f "$NEXT_PIPELINE" ]; then
    echo "[안내] ${NEXT_PIPELINE} 없음 → Pipeline 2 건너뜀"
    echo "  수동 실행: cd result && bash run_pipeline_analysis.sh sample_coords.txt sample_list.txt metadata.txt"
    exit 0
fi

if [ ! -d "$RESULT_DIR" ]; then
    echo "[안내] ${RESULT_DIR} 없음 → 처리된 샘플이 없을 수 있습니다. Pipeline 2 건너뜀"
    exit 0
fi

echo "[$(date +%H:%M:%S)] Pipeline 2로 이어서 실행..."
cp "$NEXT_PIPELINE" "$RESULT_DIR/run_pipeline_analysis.sh"

# 공통 파일 복사
for f in "$SAMPLE_COORDS" "$SAMPLE_LIST_FILE" "$METADATA_FILE"; do
    [ -f "${SCRIPT_DIR}/${f}" ] && cp "${SCRIPT_DIR}/${f}" "$RESULT_DIR/" && echo "  [복사] ${f}"
done

# 분석 스크립트 전체 복사
cp "${SCRIPT_DIR}"/*.py "${SCRIPT_DIR}"/*.R "$RESULT_DIR/" 2>/dev/null || true
for sh in Step1_organize_samples.sh Z1_organize_eggnog_outputs.sh \
          P1_extract_pathway_counts.sh P2_aggregate_pathway_table.sh \
          P4_host_tracing_join.sh P5_aggregate_host_tracing.sh; do
    [ -f "${SCRIPT_DIR}/${sh}" ] && cp "${SCRIPT_DIR}/${sh}" "$RESULT_DIR/"
done
echo "  [복사] 분석 스크립트 전체 → result/"

(
    cd "$RESULT_DIR" || exit 1
    bash run_pipeline_analysis.sh \
        "${SAMPLE_COORDS}" \
        "${SAMPLE_LIST_FILE}" \
        "${METADATA_FILE}"
)
