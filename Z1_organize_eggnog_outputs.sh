#!/bin/bash
# ==============================================================
# Z1_organize_eggnog_outputs.sh
# 전체 샘플 eggNOG-mapper 실행 완료 후, 관련 파일을 eggnog_analysis/로 정리.
# ※ SCRIPT_DIR(이 파일이 있는 폴더)에서 실행됩니다.
# ==============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# sample_list.txt 탐색: 현재 폴더 우선
if [ -f "${SCRIPT_DIR}/sample_list.txt" ]; then
    SAMPLE_LIST="${SCRIPT_DIR}/sample_list.txt"
else
    echo "❌ [Z1] sample_list.txt를 찾을 수 없습니다: ${SCRIPT_DIR}/sample_list.txt"
    exit 1
fi

DEST="${SCRIPT_DIR}/eggnog_analysis"

# ---- 0) 모든 샘플이 완료됐는지 확인 ----
N_DONE=$(ls -1 "${SCRIPT_DIR}"/*_emapper.emapper.annotations 2>/dev/null | wc -l)
N_TOTAL=$(grep -c . "$SAMPLE_LIST")
echo "[Z1] 완료된 eggNOG-mapper 결과: ${N_DONE} / ${N_TOTAL}"

if [ "$N_DONE" -lt "$N_TOTAL" ]; then
    echo "❌ [Z1] 아직 배치가 다 안 끝났습니다 (${N_DONE}/${N_TOTAL})."
    echo "   run_all_49.log 등 진행 중인 로그를 확인하고, 다 끝난 뒤 다시 실행하세요."
    exit 1
fi

echo "[Z1] 전체 완료 확인. 정리를 시작합니다..."
mkdir -p "$DEST"

# ---- 1) 테스트 파일 삭제 ----
echo "[Z1] 테스트 파일 삭제 중..."
rm -f "${SCRIPT_DIR}"/tiny_test* "${SCRIPT_DIR}"/medium_test* "${SCRIPT_DIR}"/large_test*
rm -rf "${SCRIPT_DIR}"/emappertmp_dmdn_*/

# ---- 2) 샘플별 eggNOG 파일 이동 ----
echo "[Z1] 샘플별 eggNOG 파일 이동 중..."
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    for pattern in \
        "${SAMPLE}_ghostkoala_upload.faa" \
        "${SAMPLE}_fgs.gff" \
        "${SAMPLE}_emapper.emapper.annotations" \
        "${SAMPLE}_emapper.emapper.hits" \
        "${SAMPLE}_emapper.emapper.seed_orthologs" \
        "${SAMPLE}_pathway_counts.tsv" \
        "${SAMPLE}_pathway_filtered.txt" \
        "${SAMPLE}_pathway_named.txt" \
        "${SAMPLE}_pathway_summary.txt" \
        "${SAMPLE}_host_tracing.tsv"; do
        f="${SCRIPT_DIR}/${pattern}"
        [ -f "$f" ] && mv "$f" "$DEST/"
    done
done < "$SAMPLE_LIST"

# ---- 3) 취합 산출물 이동 ----
echo "[Z1] 취합 결과물 이동 중..."
for f in \
    all_samples_pathway_long.tsv \
    pathway_abundance_wide_count.tsv \
    pathway_abundance_wide_pct.tsv \
    pathway_names_reference.tsv \
    kegg_pathway_names.tsv \
    kegg_pathway_names_clean.tsv \
    ancombc2_pathway_results.csv; do
    fp="${SCRIPT_DIR}/${f}"
    [ -f "$fp" ] && mv "$fp" "$DEST/"
done

# ---- 4) 스크립트/로그 사본 (원본 유지, 덮어쓰지 않음) ----
echo "[Z1] 관련 스크립트/로그 사본 중..."
mkdir -p "$DEST/scripts"
for script in \
    G1_prep_ghostkoala.sh E1_setup_eggnogmapper.sh \
    P1_extract_pathway_counts.sh P2_aggregate_pathway_table.sh \
    P4_host_tracing_join.sh A1_ancombc2_pathway.R; do
    fp="${SCRIPT_DIR}/${script}"
    [ -f "$fp" ] && cp -n "$fp" "$DEST/scripts/" 2>/dev/null || true
done
[ -f "${SCRIPT_DIR}/run_all_49.log" ] && \
    cp -n "${SCRIPT_DIR}/run_all_49.log" "$DEST/" 2>/dev/null || true

echo ""
echo "=== [Z1] 완료 ==="
echo "eggNOG 관련 파일 → ${DEST}/"
ls "$DEST" | head -20
echo "..."
echo "총 $(ls "$DEST" | wc -l)개 파일/폴더"
