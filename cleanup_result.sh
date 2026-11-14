#!/bin/bash
# ==============================================================
# cleanup_result.sh
# result/ 폴더를 깔끔하게 정리합니다.
#
# 정리 후 구조:
#   result/
#   ├── raw_data/         ← fastq.gz, fasta.gz (원시 데이터)
#   ├── plots/            ← Step2 시각화 (그대로)
#   ├── kraken/           ← Kraken2 결과 (그대로)
#   ├── vfdb/             ← VF 주석 결과 (그대로)
#   ├── arg/              ← ARG 결과 (그대로)
#   ├── host_tracing/     ← Host tracing 결과 (그대로)
#   ├── merged_tables/    ← ASV/taxonomy 테이블 (그대로)
#   ├── eggnog_analysis/  ← KEGG/eggNOG 결과 (그대로)
#   ├── logs/             ← 실행 로그 (그대로)
#   ├── misc_results/     ← 기타 결과 파일 (sample_map, depth 등)
#   └── *.sh *.py *.R     ← 스크립트 (그대로 — 경로 의존성 있음)
#
# 사용법: bash cleanup_result.sh
# 실행 위치: result/ 폴더 안에서 실행
# ==============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "result/ 정리 시작: $(date)"
echo "작업 위치: $SCRIPT_DIR"
echo "============================================================"

# ──────────────────────────────────────────────
# 1) raw_data/ — fastq.gz, fasta.gz
# ──────────────────────────────────────────────
echo ""
echo "[1] raw_data/ 생성 및 fastq/fasta 이동..."
mkdir -p raw_data

moved=0
for f in *.fastq.gz *.fasta.gz *.fastq *.fasta; do
    [ -f "$f" ] || continue
    mv "$f" raw_data/
    moved=$((moved+1))
done
echo "    → ${moved}개 파일 이동 완료"

# ──────────────────────────────────────────────
# 2) misc_results/ — 기타 결과 파일
# ──────────────────────────────────────────────
echo ""
echo "[2] misc_results/ 생성 및 기타 결과 이동..."
mkdir -p misc_results

moved=0
for f in \
    sample_map.html \
    sequencing_depth_AL_vs_PERM.png \
    sequencing_depth_per_sample.csv \
    pipeline_execution_log_*.txt \
    pipeline_analysis_*.log; do
    [ -f "$f" ] || continue
    mv "$f" misc_results/
    moved=$((moved+1))
done
echo "    → ${moved}개 파일 이동 완료"

# ──────────────────────────────────────────────
# 3) 불필요한 임시 파일 삭제
# ──────────────────────────────────────────────
echo ""
echo "[3] 임시 파일 정리..."
removed=0

# __pycache__
find . -maxdepth 3 -type d -name "__pycache__" | while read d; do
    rm -rf "$d"
    echo "    삭제: $d"
done

# .tmp 파일
for f in *.tmp *_raw.tmp; do
    [ -f "$f" ] || continue
    rm -f "$f"
    removed=$((removed+1))
done
echo "    → 임시 파일 ${removed}개 삭제"

# ──────────────────────────────────────────────
# 4) 완료 후 구조 출력
# ──────────────────────────────────────────────
echo ""
echo "============================================================"
echo "정리 완료: $(date)"
echo ""
echo "현재 result/ 구조:"
echo "  raw_data/         → $(ls raw_data/ 2>/dev/null | wc -l)개 파일 (fastq/fasta)"
echo "  plots/            → $(ls plots/ 2>/dev/null | wc -l)개 파일"
echo "  kraken/           → $(ls kraken/ 2>/dev/null | wc -l)개 파일"
echo "  vfdb/             → $(ls vfdb/ 2>/dev/null | wc -l)개 파일"
echo "  arg/              → $(ls arg/ 2>/dev/null | wc -l)개 항목"
echo "  host_tracing/     → $(ls host_tracing/ 2>/dev/null | wc -l)개 파일"
echo "  merged_tables/    → $(ls merged_tables/ 2>/dev/null | wc -l)개 파일"
echo "  eggnog_analysis/  → $(ls eggnog_analysis/ 2>/dev/null | wc -l)개 항목"
echo "  logs/             → $(ls logs/ 2>/dev/null | wc -l)개 파일"
echo "  misc_results/     → $(ls misc_results/ 2>/dev/null | wc -l)개 파일"
echo ""
echo "※ 스크립트(.sh/.py/.R)는 경로 의존성 때문에 그대로 유지됩니다."
echo "  재실행 시: bash run_pipeline_analysis.sh sample_coords.txt sample_list.txt metadata.txt"
echo "============================================================"
