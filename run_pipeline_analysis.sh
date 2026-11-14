#!/bin/bash
# ==============================================================
# run_pipeline_analysis.sh  ─ Pipeline 2: 통합 분석 (Step 0~11)
#
# [Stage 1]  Step 0~4  : 좌표 시각화, VF/ARG 분석, host-tracing, taxonomy 병합
# [Stage 2]  Step 5~11 : eggNOG/KEGG pathway 기능 분석 및 통계 검증
#
# 사용법:
#   bash run_pipeline_analysis.sh sample_coords.txt sample_list.txt metadata.txt
#   bash run_pipeline_analysis.sh sample_coords.txt sample_list.txt metadata.txt 7
#   bash run_pipeline_analysis.sh "" sample_list.txt metadata.txt 5-10
#
# 환경변수:
#   TAXDUMP        : NCBI taxdump(nodes.dmp/names.dmp) 상위 경로
#   JOBS           : Step3 병렬 처리 수
#   CONDA_BASE     : miniforge/anaconda 루트 경로
#   EGGNOG_ENV_PY  : eggnog conda 환경의 python 경로
#   ANCOMBC_ENV_R  : r_ancombc conda 환경의 Rscript 경로
#   GOOGLE_MAPS_API_KEY : (선택) Step0 구글맵 API 키
# ==============================================================
set -uo pipefail

# ──────────────────────────────────────────────
# 0) 경로·파라미터 설정
# ──────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COORDS_FILE="${1:-sample_coords.txt}"
SAMPLE_LIST="${2:-sample_list.txt}"
METADATA="${3:-metadata.txt}"
ONLY_STEP="${4:-all}"

TAXDUMP="${TAXDUMP:-/mnt/program/db/plusPFP_16}"
JOBS="${JOBS:-6}"
CONDA_BASE="${CONDA_BASE:-/root/miniforge3}"
EGGNOG_PY="${EGGNOG_ENV_PY:-${CONDA_BASE}/envs/eggnog/bin/python}"
ANCOMBC_R="${ANCOMBC_ENV_R:-${CONDA_BASE}/envs/r_ancombc/bin/Rscript}"

# 출력 디렉토리 (모두 SCRIPT_DIR 기준 절대경로로 관리)
PLOTS_DIR="${SCRIPT_DIR}/plots"
KRAKEN_DIR="${SCRIPT_DIR}/kraken"
VFDB_DIR="${SCRIPT_DIR}/vfdb"
ARG_DIR="${SCRIPT_DIR}/arg"
HOST_TRACING_DIR="${SCRIPT_DIR}/host_tracing"
MERGED_DIR="${SCRIPT_DIR}/merged_tables"
EGGNOG_DIR="${SCRIPT_DIR}/eggnog_analysis"

# ──────────────────────────────────────────────
# 1) 로그 설정
# ──────────────────────────────────────────────
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_analysis_$(date +%Y%m%d_%H%M%S).log"

# stdout + stderr 모두 로그 파일에 기록 (tee로 터미널에도 출력)
exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Pipeline 2 시작: $(date)"
echo "  coords=$COORDS_FILE  samples=$SAMPLE_LIST"
echo "  metadata=$METADATA   only_step=$ONLY_STEP"
echo "  TAXDUMP=$TAXDUMP"
echo "  로그: $LOG_FILE"
echo "============================================================"

# ──────────────────────────────────────────────
# 2) sample_list.txt 자동 생성 (없는 경우)
# ──────────────────────────────────────────────
if [ ! -f "${SCRIPT_DIR}/${SAMPLE_LIST}" ] && [ ! -f "$SAMPLE_LIST" ]; then
    SAMPLE_LIST_FULL=""
else
    # 상대경로이면 SCRIPT_DIR 기준으로 절대화
    if [[ "$SAMPLE_LIST" != /* ]]; then
        SAMPLE_LIST_FULL="${SCRIPT_DIR}/${SAMPLE_LIST}"
    else
        SAMPLE_LIST_FULL="$SAMPLE_LIST"
    fi
fi

# coords에서 자동 추출 시도
if [ -z "${SAMPLE_LIST_FULL:-}" ] || [ ! -f "${SAMPLE_LIST_FULL}" ]; then
    COORDS_FULL="${SCRIPT_DIR}/${COORDS_FILE}"
    [[ "$COORDS_FILE" == /* ]] && COORDS_FULL="$COORDS_FILE"
    if [ -f "$COORDS_FULL" ]; then
        SAMPLE_LIST_FULL="${SCRIPT_DIR}/sample_list.txt"
        echo "[안내] sample_list.txt 없음 → ${COORDS_FILE}에서 자동 추출"
        tail -n +2 "$COORDS_FULL" | cut -f1 > "$SAMPLE_LIST_FULL"
    else
        echo "오류: sample_list.txt도 sample_coords.txt도 없습니다."
        exit 1
    fi
fi

# metadata 절대경로
if [[ "$METADATA" != /* ]]; then
    METADATA_FULL="${SCRIPT_DIR}/${METADATA}"
else
    METADATA_FULL="$METADATA"
fi

# coords 절대경로
if [[ "$COORDS_FILE" != /* ]]; then
    COORDS_FULL="${SCRIPT_DIR}/${COORDS_FILE}"
else
    COORDS_FULL="$COORDS_FILE"
fi

echo "[경로] sample_list: $SAMPLE_LIST_FULL"
echo "[경로] metadata:    $METADATA_FULL"
echo "[경로] coords:      $COORDS_FULL"

# ──────────────────────────────────────────────
# 3) 사전 검사 (Pre-flight check)
# ──────────────────────────────────────────────
echo ""
echo "──────────────────────────────────────────────"
echo "[사전검사] 필수 파일 및 스크립트 확인 중..."
echo "──────────────────────────────────────────────"
PREFLIGHT_FAIL=0

# metadata 검사
if [ ! -s "$METADATA_FULL" ]; then
    echo "[사전검사 실패] ${METADATA_FULL} 가 없거나 비어있습니다."
    echo "  형식(탭 구분): sample-id<TAB>group"
    echo "  예시:"
    echo '    printf "sample-id\tgroup\n" > metadata.txt'
    echo '    printf "SRR001\tAL\n" >> metadata.txt'
    echo '    printf "SRR002\tPERM\n" >> metadata.txt'
    PREFLIGHT_FAIL=1
else
    # 탭 구분자 확인
    if ! head -1 "$METADATA_FULL" | grep -qP '\t'; then
        echo "[사전검사 실패] ${METADATA_FULL}에 탭 구분자가 없습니다."
        echo "  cat -A ${METADATA_FULL} 로 확인하세요 (^I가 탭 문자)"
        PREFLIGHT_FAIL=1
    fi

    # sample_list의 모든 샘플이 metadata에 있는지 확인
    MISSING_IN_META=""
    while read -r s; do
        [ -z "$s" ] && continue
        grep -qP "^${s}\t" "$METADATA_FULL" || MISSING_IN_META="${MISSING_IN_META} ${s}"
    done < "$SAMPLE_LIST_FULL"
    if [ -n "$MISSING_IN_META" ]; then
        echo "[사전검사 실패] 아래 샘플이 ${METADATA_FULL}에 없습니다:"
        for s in $MISSING_IN_META; do echo "  - $s"; done
        PREFLIGHT_FAIL=1
    fi

    # AL/PERM 개수 확인 (경고만)
    N_AL=$(awk -F'\t' 'NR>1 && $2=="AL"' "$METADATA_FULL" | wc -l)
    N_PERM=$(awk -F'\t' 'NR>1 && $2=="PERM"' "$METADATA_FULL" | wc -l)
    echo "[사전검사] AL: ${N_AL}개, PERM: ${N_PERM}개"
    if [ "$N_AL" -lt 2 ] || [ "$N_PERM" -lt 2 ]; then
        echo "[사전검사 경고] AL/PERM 각각 2개 미만 → Step7(ANCOM-BC2)이 실패할 수 있습니다."
    fi
fi

# 필수 스크립트 확인
for req_script in \
    Step1_organize_samples.sh \
    Step2_vf_arg_plots.py \
    Step3_host_tracing.py \
    Step4_merge_kreports.py \
    VF_copies_per_cell.py \
    Z1_organize_eggnog_outputs.sh \
    reorder_heatmap.py \
    publication_boxplot.py \
    all_pathways_test.py \
    P1_extract_pathway_counts.sh \
    P2_aggregate_pathway_table.sh \
    P4_host_tracing_join.sh \
    P5_aggregate_host_tracing.sh \
    filter_ancombc2_by_category.py \
    ancombc2_volcano_plot.py \
    top_taxa_barplot.py \
    host_tracing_barplot.py \
    taxon_pathway_heatmap.py \
    taxon_pathway_heatmap_methanogens.py \
    verify_methanogenesis_taxon_source.py \
    batch_pathway_to_taxon.py \
    visualize_robust_pathways_taxon.py \
    check_sequencing_depth.py \
    host_tracing_normalized.py \
    A1_ancombc2_pathway.R \
    A2_ancombc2_balanced_subsample.R \
    A3_ancombc2_genus_corrected.R; do
    if [ ! -f "${SCRIPT_DIR}/${req_script}" ]; then
        echo "[사전검사 실패] 필수 스크립트 없음: ${req_script}"
        PREFLIGHT_FAIL=1
    fi
done

# taxdump 확인 (Step3/4에서 필요)
if [ ! -f "${TAXDUMP}/nodes.dmp" ] || [ ! -f "${TAXDUMP}/names.dmp" ]; then
    echo "[사전검사 실패] TAXDUMP 경로에 nodes.dmp / names.dmp 가 없습니다: $TAXDUMP"
    echo "  https://ftp.ncbi.nlm.nih.gov/pub/taxdump/ 에서 taxdump.tar.gz 다운로드 후 압축 해제"
    echo "  또는: TAXDUMP=/your/path bash $0 ..."
    PREFLIGHT_FAIL=1
fi

if [ "$PREFLIGHT_FAIL" -eq 1 ]; then
    echo ""
    echo "[사전검사 중단] 위 문제를 먼저 해결한 뒤 다시 실행하세요."
    exit 1
fi
echo "[사전검사 통과] 모든 파일·스크립트 확인 완료"
echo ""

# ──────────────────────────────────────────────
# 4) 헬퍼 함수
# ──────────────────────────────────────────────

# 이 step을 실행해야 하는지 판단
should_run_step() {
    local step_num="$1"
    if [ "$ONLY_STEP" = "all" ]; then return 0; fi
    if [[ "$ONLY_STEP" == *-* ]]; then
        local s="${ONLY_STEP%-*}" e="${ONLY_STEP#*-}"
        [ "$step_num" -ge "$s" ] && [ "$step_num" -le "$e" ] && return 0 || return 1
    fi
    [ "$ONLY_STEP" = "$step_num" ]
}

# step 실행 (실패 시 파이프라인 중단)
run_step() {
    local step_num="$1" desc="$2"; shift 2
    should_run_step "$step_num" || return 0

    local step_log="${LOG_DIR}/step${step_num}_${desc}_$(date +%Y%m%d_%H%M%S).log"
    echo ""
    echo "──────────────────────────────────────────────"
    echo "[$(date +%H:%M:%S)] Step${step_num} ${desc} 시작"
    echo "  로그: $step_log"
    echo "──────────────────────────────────────────────"

    if "$@" > "$step_log" 2>&1; then
        echo "[$(date +%H:%M:%S)] Step${step_num} ${desc} ✓ 완료"
    else
        echo "[$(date +%H:%M:%S)] Step${step_num} ${desc} ✗ 실패"
        echo "  → 오류 내용 (마지막 30줄):"
        tail -30 "$step_log" | sed 's/^/    /'
        echo "  → 전체 로그: $step_log"
        exit 1
    fi
}

# 지정 conda 환경에서 명령 실행
CONDA_BASE="${CONDA_BASE:-/root/miniforge3}"
run_in_env() {
    local env_name="$1"; shift
    bash -c "
        source '${CONDA_BASE}/etc/profile.d/conda.sh'
        conda activate '${env_name}'
        $*
    "
}

# 출력 디렉토리 사전 생성 (경로 못 찾음 방지)
ensure_dirs() {
    for d in "$@"; do
        mkdir -p "$d"
    done
}

# ──────────────────────────────────────────────
# 5) 출력 디렉토리 사전 생성
# ──────────────────────────────────────────────
ensure_dirs \
    "$PLOTS_DIR" \
    "$KRAKEN_DIR" \
    "$VFDB_DIR" \
    "$ARG_DIR" \
    "$HOST_TRACING_DIR" \
    "$MERGED_DIR" \
    "$EGGNOG_DIR" \
    "$LOG_DIR"

echo "[준비] 출력 디렉토리 생성 완료:"
echo "  plots/         → ${PLOTS_DIR}"
echo "  kraken/        → ${KRAKEN_DIR}"
echo "  vfdb/          → ${VFDB_DIR}"
echo "  arg/           → ${ARG_DIR}"
echo "  host_tracing/  → ${HOST_TRACING_DIR}"
echo "  merged_tables/ → ${MERGED_DIR}"
echo "  eggnog_analysis/ → ${EGGNOG_DIR}"
echo "  logs/          → ${LOG_DIR}"

# ══════════════════════════════════════════════
# Stage 1: VF / ARG / Taxonomy
# ══════════════════════════════════════════════

# Step 0: 구글맵 시각화 (좌표 파일 있을 때만)
if should_run_step 0; then
    if [ -f "$COORDS_FULL" ]; then
        run_step 0 "make_googlemap" \
            python3 "${SCRIPT_DIR}/Step0_make_googlemap.py" \
                "$COORDS_FULL" \
                --api_key "${GOOGLE_MAPS_API_KEY:-YOUR_GOOGLE_MAPS_API_KEY}" \
                --outdir "$SCRIPT_DIR"
    else
        echo "[Step0] ${COORDS_FULL} 없음 → 건너뜀"
    fi
fi

# Step 1: 샘플별 파일 폴더 정리 (kraken/ vfdb/ arg/)
run_step 1 "organize_samples" \
    bash "${SCRIPT_DIR}/Step1_organize_samples.sh" "$SAMPLE_LIST_FULL"

# Step 2: VF/ARG 바그래프 + 히트맵
run_step 2 "vf_arg_plots" \
    python3 "${SCRIPT_DIR}/Step2_vf_arg_plots.py" \
        --sample_list "$SAMPLE_LIST_FULL" \
        --metadata    "$METADATA_FULL" \
        --vfdb_dir    "$VFDB_DIR" \
        --arg_dir     "$ARG_DIR" \
        --outdir      "$PLOTS_DIR"

# Step 2 추가: VF copies-per-cell
run_step 2 "vf_copies_per_cell" \
    python3 "${SCRIPT_DIR}/VF_copies_per_cell.py" \
        --sample_list "$SAMPLE_LIST_FULL" \
        --metadata    "$METADATA_FULL" \
        --vfdb_dir    "$VFDB_DIR" \
        --arg_dir     "$ARG_DIR" \
        --outdir      "$PLOTS_DIR"

# Step 3: Host tracing (VF/ARG × Kraken2)
run_step 3 "host_tracing" \
    python3 "${SCRIPT_DIR}/Step3_host_tracing.py" \
        --sample_list "$SAMPLE_LIST_FULL" \
        --metadata    "$METADATA_FULL" \
        --kraken_dir  "$KRAKEN_DIR" \
        --vfdb_dir    "$VFDB_DIR" \
        --arg_dir     "$ARG_DIR" \
        --taxdump     "$TAXDUMP" \
        --jobs        "$JOBS" \
        --outdir      "$HOST_TRACING_DIR"

# Step 4: Kreport 병합 → ASV/taxonomy/metadata 테이블
run_step 4 "merge_kreports" \
    python3 "${SCRIPT_DIR}/Step4_merge_kreports.py" \
        --sample_list "$SAMPLE_LIST_FULL" \
        --metadata    "$METADATA_FULL" \
        --kraken_dir  "$KRAKEN_DIR" \
        --taxdump     "$TAXDUMP" \
        --outdir      "$MERGED_DIR"

# ══════════════════════════════════════════════
# Stage 2: eggNOG / KEGG Pathway 분석
# ══════════════════════════════════════════════

# Step 5: eggNOG-mapper 결과 정리 → eggnog_analysis/
run_step 5 "organize_eggnog_outputs" \
    bash -c "
        cd '${SCRIPT_DIR}'
        source '${CONDA_BASE}/etc/profile.d/conda.sh'
        conda activate eggnog
        bash '${SCRIPT_DIR}/Z1_organize_eggnog_outputs.sh'
    "

# Step 6a: KEGGaNOG 커뮤니티 완성도 분석
run_step 6 "run_kegganog" \
    bash -c "
        cd '${EGGNOG_DIR}'
        ls *_emapper.emapper.annotations 2>/dev/null > sample_list_for_kegganog.txt
        if [ ! -s sample_list_for_kegganog.txt ]; then
            echo '[오류] ${EGGNOG_DIR} 에 *.emapper.annotations 파일이 없습니다.'
            exit 1
        fi
        rm -rf kegganog_results
        '${CONDA_BASE}/envs/eggnog/bin/KEGGaNOG' -M -i sample_list_for_kegganog.txt -o kegganog_results
    "

# Step 6b: 히트맵 재정렬
run_step 6 "reorder_heatmap" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/reorder_heatmap.py'"

# Step 6c: 출판용 박스플롯
run_step 6 "publication_boxplot" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/publication_boxplot.py'"

# Step 6d: 전체 pathway 검정
run_step 6 "all_pathways_test" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/all_pathways_test.py'"

# Step 7a: Pathway abundance 집계 (P1 → P2)
run_step 7 "prep_pathway_abundance" \
    bash -c "
        cd '${EGGNOG_DIR}'
        # P1: 샘플별 pathway count
        while read -r s; do
            [ -z \"\$s\" ] && continue
            bash '${SCRIPT_DIR}/P1_extract_pathway_counts.sh' \"\$s\"
        done < '${SAMPLE_LIST_FULL}'
        # P2: wide table 생성
        # P2는 ../sample_list.txt를 기대하므로 심볼릭 링크로 연결
        ln -sf '${SAMPLE_LIST_FULL}' sample_list.txt 2>/dev/null || \
            cp '${SAMPLE_LIST_FULL}' sample_list.txt
        bash '${SCRIPT_DIR}/P2_aggregate_pathway_table.sh'
    "

# Step 7b: sample_group_mapping.csv 사전 생성 (Step 8~9에서도 필요하므로 run_step 밖에서 먼저 생성)
mkdir -p "${EGGNOG_DIR}"
{
    echo 'Sample,Group'
    tail -n +2 "${METADATA_FULL}" | awk -F'\t' '{print $1","$2}'
} > "${EGGNOG_DIR}/sample_group_mapping.csv"
echo "[준비] sample_group_mapping.csv 사전 생성 완료"

run_step 7 "make_group_mapping" \
    bash -c "
        cd '${EGGNOG_DIR}'
        echo 'Sample,Group' > sample_group_mapping.csv
        tail -n +2 '${METADATA_FULL}' | awk -F'\t' '{print \$1\",\"\$2}' >> sample_group_mapping.csv
        echo '[Step7] sample_group_mapping.csv 생성 완료'
    "

# Step 7c: ANCOM-BC2 pathway 분석
run_step 7 "ancombc2_pathway" \
    bash -c "cd '${EGGNOG_DIR}' && '${ANCOMBC_R}' '${SCRIPT_DIR}/A1_ancombc2_pathway.R'"

# Step 7d: 카테고리 필터링
run_step 7 "filter_ancombc2_category" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/filter_ancombc2_by_category.py'"

# Step 7e: Volcano plot
run_step 7 "ancombc2_volcano" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/ancombc2_volcano_plot.py'"

# Step 7f: 균형 서브샘플 robustness 검정 (20회)
run_step 7 "ancombc2_balanced_subsample" \
    bash -c "cd '${EGGNOG_DIR}' && '${ANCOMBC_R}' '${SCRIPT_DIR}/A2_ancombc2_balanced_subsample.R'"

# Step 8a: Host-tracing join (P4: eggNOG × Kraken2)
run_step 8 "host_tracing_join" \
    bash -c "
        cd '${EGGNOG_DIR}'
        # .kraken 파일을 현재 폴더에서도 찾을 수 있도록 kraken/ 폴더 심볼릭 링크
        [ ! -e kraken ] && ln -sf '${KRAKEN_DIR}' kraken
        while read -r s; do
            [ -z \"\$s\" ] && continue
            # .kraken 파일: P4가 현재 폴더 → 상위 폴더 순으로 찾으므로
            # eggnog_analysis/ 에 없으면 SCRIPT_DIR에서 복사 (원본 보존)
            if [ ! -s \"\${s}.kraken\" ] && [ -s '${SCRIPT_DIR}/'\"\${s}.kraken\" ]; then
                ln -sf '${SCRIPT_DIR}/'\"\${s}.kraken\" \"\${s}.kraken\"
            fi
            bash '${SCRIPT_DIR}/P4_host_tracing_join.sh' \"\$s\"
        done < '${SAMPLE_LIST_FULL}'
    "

# Step 8b: Host-tracing 집계 (P5)
run_step 8 "host_tracing_aggregate" \
    bash -c "
        cd '${EGGNOG_DIR}'
        # P5는 ../sample_list.txt를 기대
        ln -sf '${SAMPLE_LIST_FULL}' sample_list.txt 2>/dev/null || \
            cp '${SAMPLE_LIST_FULL}' sample_list.txt
        bash '${SCRIPT_DIR}/P5_aggregate_host_tracing.sh' || true
    "

# Step 8c: 시각화
run_step 8 "top_taxa_barplot" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/top_taxa_barplot.py'"

run_step 8 "host_tracing_barplot" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/host_tracing_barplot.py'"

run_step 8 "taxon_pathway_heatmap" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/taxon_pathway_heatmap.py'"

run_step 8 "taxon_pathway_heatmap_methanogens" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/taxon_pathway_heatmap_methanogens.py'"

# Step 9: 심화 taxon 귀속 검증
run_step 9 "verify_methanogenesis_taxon_source" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/verify_methanogenesis_taxon_source.py'"

run_step 9 "batch_pathway_to_taxon" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/batch_pathway_to_taxon.py'"

run_step 9 "visualize_robust_pathways_taxon" \
    bash -c "cd '${EGGNOG_DIR}' && '${EGGNOG_PY}' '${SCRIPT_DIR}/visualize_robust_pathways_taxon.py'"

# Step 10: 시퀀싱 깊이 교란요인 검증
run_step 10 "check_sequencing_depth" \
    bash -c "
        cd '${SCRIPT_DIR}'
        '${EGGNOG_PY}' '${SCRIPT_DIR}/check_sequencing_depth.py' \
            --sample_list '${SAMPLE_LIST_FULL}' \
            --outdir '${SCRIPT_DIR}'
    "

run_step 10 "host_tracing_normalized" \
    bash -c "
        # sequencing_depth_per_sample.csv를 eggnog_analysis/로 복사
        DEPTH_CSV='${SCRIPT_DIR}/sequencing_depth_per_sample.csv'
        if [ ! -s \"\$DEPTH_CSV\" ]; then
            echo '[오류] sequencing_depth_per_sample.csv 없음. Step10 check_sequencing_depth를 먼저 실행하세요.'
            exit 1
        fi
        cp \"\$DEPTH_CSV\" '${EGGNOG_DIR}/'
        cd '${EGGNOG_DIR}'
        '${EGGNOG_PY}' '${SCRIPT_DIR}/host_tracing_normalized.py'
    "

# Step 11: Genus 수준 ANCOM-BC2 최종 검증 (Step4 산출물 재사용)
run_step 11 "ancombc2_genus_final" \
    bash -c "
        ASV_TABLE='${MERGED_DIR}/1_ASV_table.txt'
        if [ ! -s \"\$ASV_TABLE\" ]; then
            echo '[오류] ${MERGED_DIR}/1_ASV_table.txt 없음. Step4가 먼저 완료되어야 합니다.'
            exit 1
        fi
        cp '${SCRIPT_DIR}/A3_ancombc2_genus_corrected.R' '${MERGED_DIR}/'
        cd '${MERGED_DIR}'
        '${ANCOMBC_R}' A3_ancombc2_genus_corrected.R
    "

# ──────────────────────────────────────────────
# 완료 요약
# ──────────────────────────────────────────────
echo ""
echo "============================================================"
echo "Pipeline 2 완료: $(date)  (STEP: ${ONLY_STEP})"
echo ""
echo "출력 위치:"
echo "  plots/           → ${PLOTS_DIR}"
echo "  host_tracing/    → ${HOST_TRACING_DIR}"
echo "  merged_tables/   → ${MERGED_DIR}"
echo "  eggnog_analysis/ → ${EGGNOG_DIR}"
echo "  logs/            → ${LOG_DIR}"
echo ""
echo "주요 산출물:"
echo "  [Step0]  sample_map.html"
echo "  [Step2]  plots/VF_*, ARG_*"
echo "  [Step3]  host_tracing/VF_host_tracing_by_group.png"
echo "  [Step4]  merged_tables/1_ASV_table.txt, 2_taxonomy_table.txt"
echo "  [Step6]  eggnog_analysis/kegganog_results/, heatmap_*.png"
echo "  [Step7]  eggnog_analysis/ancombc2_pathway_results.csv"
echo "  [Step8]  eggnog_analysis/host_tracing_master.tsv, *barplot*.png"
echo "  [Step11] merged_tables/ancombc2_genus_corrected_results.csv"
echo "============================================================"
