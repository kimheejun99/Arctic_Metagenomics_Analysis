#!/bin/bash
# ==============================================================
# G1_prep_ghostkoala.sh
# GhostKOALA(웹, 수동 업로드)에 넣을 파일을 준비하는 단계.
# 이 스크립트가 하는 일(자동화 대상):
#   1) 49샘플 전부 동일 비율/동일 시드로 read 서브샘플링 (재현성 확보)
#   2) FragGeneScan(Rs)으로 assembly 없이 read에서 바로 유전자(protein) 예측
#   3) GhostKOALA 업로드 한도(서열 50만개 / 300MB) 안에 들어오는지 자동 확인
#
# 이 스크립트가 안 하는 일(직접 하셔야 함, 자동화 대상 아님):
#   - GhostKOALA 웹사이트에 파일 업로드
#   - 이메일 인증 클릭
#   - 결과 다운로드 후 KEGG Mapper Reconstruct에 붙여넣기
#
# 사용법: bash G1_prep_ghostkoala.sh <SRA_RUN_ID>
# 조정 가능(환경변수):
#   SUBSAMPLE_FRACTION (기본 0.01 = 1%)
#   SUBSAMPLE_SEED      (기본 100, 49샘플 전부 반드시 동일 값 사용)
#   THREADS             (기본 16)
# ==============================================================
set -euo pipefail

SAMPLE=$1
FRACTION="${SUBSAMPLE_FRACTION:-0.01}"
SEED="${SUBSAMPLE_SEED:-100}"
THREADS="${THREADS:-16}"

OUT_FAA="${SAMPLE}_ghostkoala_upload.faa"

# 이미 처리된 샘플이면 스킵 (기존 파이프라인과 동일한 idempotent 패턴)
if [[ -s "$OUT_FAA" ]]; then
    echo "[G1] ${SAMPLE}: ${OUT_FAA} 가 이미 존재합니다 -> 스킵"
    exit 0
fi

# ---- FragGeneScan(Rs) 준비 (없으면 최초 1회 설치) ----
# 먼저 지금 활성 환경에 그냥 바로 설치를 시도합니다 (구버전 fraggenescan은
# 순수 C/Perl 실행파일이라 보통 python 버전과 무관하게 잘 들어갑니다).
# 이것마저 실패할 때만 전용 환경으로 분리합니다.
FGS_BIN=""
for cand in FragGeneScanRs FragGeneScan run_FragGeneScan.pl; do
    if command -v "$cand" &> /dev/null; then FGS_BIN="$cand"; break; fi
done

if [ -z "$FGS_BIN" ]; then
    INSTALLER="conda"; command -v mamba &> /dev/null && INSTALLER="mamba"

    echo "[G1] 지금 활성 환경에 fraggenescan 설치를 바로 시도합니다..."
    if "$INSTALLER" install -y --override-channels -c bioconda -c conda-forge fraggenescan 2>/dev/null; then
        for cand in FragGeneScan run_FragGeneScan.pl; do
            if command -v "$cand" &> /dev/null; then FGS_BIN="$cand"; break; fi
        done
    fi

    # 그래도 안 되면 전용 환경(fgs_env)으로 분리 (진짜 버전 충돌이었을 경우의 대비책)
    if [ -z "$FGS_BIN" ]; then
        echo "[G1] 바로 설치는 실패 -> 전용 환경(fgs_env)으로 재시도합니다..."
        source ~/miniforge3/etc/profile.d/conda.sh
        if ! conda env list | grep -q "^fgs_env "; then
            if ! "$INSTALLER" create -y -n fgs_env --override-channels -c bioconda -c conda-forge fraggenescan; then
                echo "[G1] fraggenescan도 실패 -> fraggenescanrs로 재시도합니다..."
                conda env remove -y -n fgs_env 2>/dev/null || true
                if ! "$INSTALLER" create -y -n fgs_env --override-channels -c bioconda -c conda-forge fraggenescanrs; then
                    echo "❌ [G1] conda로 설치 실패했습니다."
                    echo "   수동 설치: cargo install fraggenescanrs"
                    echo "   (cargo 없으면: curl https://sh.rustup.rs -sSf | sh)"
                    exit 1
                fi
            fi
        fi
        conda activate fgs_env
        for cand in FragGeneScan FragGeneScanRs run_FragGeneScan.pl; do
            if command -v "$cand" &> /dev/null; then FGS_BIN="$cand"; break; fi
        done
    fi

    if [ -z "$FGS_BIN" ]; then
        echo "❌ [G1] 설치는 됐는데 실행파일을 못 찾았습니다. 'which FragGeneScan' 으로 직접 확인해보세요."
        exit 1
    fi
fi
echo "[G1] 사용할 gene prediction 도구: $FGS_BIN"

# ---- [버그 수정] conda로 설치된 FragGeneScan은 'train' 파라미터 폴더를
# 실행 위치(cwd) 기준 상대경로로 찾는데, conda 설치 구조상 실제로는
# $CONDA_PREFIX/share/... 안에 들어있어서 못 찾는 경우가 흔합니다.
# 실제 train 폴더를 찾아서 지금 작업 디렉토리에 심볼릭 링크로 연결합니다.
if [ "$FGS_BIN" = "FragGeneScan" ] && [ ! -e "./train" ]; then
    echo "[G1] FragGeneScan의 train 파라미터 폴더 위치를 찾는 중..."
    TRAIN_SRC=$(find "${CONDA_PREFIX:-/root/miniforge3}" -type d -name "train" -path "*fraggenescan*" 2>/dev/null | head -1)
    if [ -z "$TRAIN_SRC" ]; then
        # 위 패턴으로 못 찾으면, illumina_10 파일 자체를 기준으로 넓게 재검색
        HIT=$(find "${CONDA_PREFIX:-/root/miniforge3}" -type f -name "illumina_10" 2>/dev/null | head -1)
        [ -n "$HIT" ] && TRAIN_SRC=$(dirname "$HIT")
    fi
    if [ -n "$TRAIN_SRC" ]; then
        echo "[G1] train 폴더 발견: $TRAIN_SRC -> ./train 로 링크 생성"
        ln -s "$TRAIN_SRC" ./train
    else
        echo "⚠️  [G1] train 폴더를 자동으로 못 찾았습니다. 아래 명령으로 직접 찾아서"
        echo "    'ln -s <경로> ./train' 을 실행한 뒤 스크립트를 다시 돌려주세요:"
        echo "    find \$CONDA_PREFIX -iname 'illumina_10'"
    fi
fi

# ---- 1) 서브샘플링 (R1/R2 동일 seed로 페어 유지) ----
R1="${SAMPLE}_1.fastq.gz"; R2="${SAMPLE}_2.fastq.gz"; SE="${SAMPLE}.fastq.gz"
SUB1="${SAMPLE}_sub_1.fastq.gz"; SUB2="${SAMPLE}_sub_2.fastq.gz"; SUBSE="${SAMPLE}_sub.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    echo "[G1] ${SAMPLE}: paired-end ${FRACTION} 비율로 서브샘플링 중 (seed=${SEED})..."
    seqkit sample -p "$FRACTION" -s "$SEED" "$R1" -o "$SUB1"
    seqkit sample -p "$FRACTION" -s "$SEED" "$R2" -o "$SUB2"
    CAT_FASTQ="${SAMPLE}_sub_cat.fastq.gz"
    cat "$SUB1" "$SUB2" > "$CAT_FASTQ"
elif [[ -f "$SE" ]]; then
    echo "[G1] ${SAMPLE}: single-end ${FRACTION} 비율로 서브샘플링 중 (seed=${SEED})..."
    seqkit sample -p "$FRACTION" -s "$SEED" "$SE" -o "$SUBSE"
    CAT_FASTQ="$SUBSE"
else
    echo "❌ [G1] ${SAMPLE}: fastq 파일을 찾을 수 없습니다."
    exit 1
fi

# ---- fasta로 변환 (FragGeneScan 입력 형식) ----
SUB_FASTA="${SAMPLE}_sub.fasta"
seqkit fq2fa "$CAT_FASTQ" -o "$SUB_FASTA"
N_READS=$(grep -c '^>' "$SUB_FASTA")
echo "[G1] ${SAMPLE}: 서브샘플링된 read 수 = ${N_READS}"

# ---- 2) FragGeneScan(Rs)으로 유전자(protein) 예측 ----
# -w 0     : 완전한 유전체가 아니라 짧고 조각난 read라는 의미 (메타지노믹스 read에 필수)
# -t illumina_10 : Illumina 짧은 read용 학습 모델 (표준적으로 권장되는 기본값)
echo "[G1] ${SAMPLE}: FragGeneScan(Rs) 실행 중..."
"$FGS_BIN" -s "$SUB_FASTA" -o "${SAMPLE}_fgs" -w 0 -t illumina_10 -p "$THREADS"

mv "${SAMPLE}_fgs.faa" "${SAMPLE}_fgs_raw.faa"

# ---- [추가] FragGeneScan이 드물게 같은 ORF를 중복 출력하는 경우가 있어
# (GhostKOALA는 중복 헤더가 있으면 업로드 자체를 거부합니다), 최종 파일을
# 만들기 전에 헤더 기준 중복을 제거합니다.
awk '/^>/{h=$0; if (h in seen) skip=1; else {seen[h]=1; skip=0}} {if (!skip) print}' \
    "${SAMPLE}_fgs_raw.faa" > "$OUT_FAA"
N_DUP=$(( $(grep -c '^>' "${SAMPLE}_fgs_raw.faa") - $(grep -c '^>' "$OUT_FAA") ))
[ "$N_DUP" -gt 0 ] && echo "[G1] ${SAMPLE}: 중복 서열 ${N_DUP}개 제거함"
rm -f "${SAMPLE}_fgs_raw.faa"

# ---- 3) GhostKOALA 업로드 한도 확인 ----
N_SEQ=$(grep -c '^>' "$OUT_FAA" || echo 0)
SIZE_MB=$(du -m "$OUT_FAA" | cut -f1)

echo ""
echo "=== [G1] ${SAMPLE} 완료: ${OUT_FAA} ==="
echo "    예측된 protein 서열 수 : ${N_SEQ}  (GhostKOALA 한도: 500,000)"
echo "    파일 크기             : ${SIZE_MB} MB  (GhostKOALA 한도: 300 MB)"

if [ "$N_SEQ" -gt 500000 ] || [ "$SIZE_MB" -gt 300 ]; then
    echo "⚠️  [G1] 한도를 초과했습니다! SUBSAMPLE_FRACTION을 더 낮춰서(예: 0.005) 재실행하세요."
    echo "    예: SUBSAMPLE_FRACTION=0.005 bash G1_prep_ghostkoala.sh ${SAMPLE}"
else
    echo "✅  [G1] 한도 이내입니다. GhostKOALA(https://www.kegg.jp/ghostkoala/)에"
    echo "    이 파일을 직접 업로드하시면 됩니다."
fi

# ---- 중간 임시파일 정리 (최종 .faa만 남김) ----
rm -f "$SUB1" "$SUB2" "$SUBSE" "$CAT_FASTQ" "$SUB_FASTA" \
      "${SAMPLE}_fgs.ffn" "${SAMPLE}_fgs.out" 2>/dev/null || true

echo "[G1] ${SAMPLE}: 임시파일 정리 완료"
