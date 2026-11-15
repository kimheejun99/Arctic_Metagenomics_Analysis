#!/bin/bash

# 1. Conda 환경 활성화
# bash 내부에서 conda를 쓰기 위해 경로 재설정
source ~/miniforge3/etc/profile.d/conda.sh
conda activate args_oap

# 변수 설정
SRA_ID=$1
WORK_DIR="temp_args_${SRA_ID}"

echo "[A2] ARGs-OAP 분석 시작: ${SRA_ID}"

# 2. 임시 작업 디렉토리 생성 및 파일 복사 (원본 보호를 위해)
mkdir -p ${WORK_DIR}
# .fasta.gz 파일을 임시 폴더로 링크 (복사보다 빠름)
ln -s $(pwd)/${SRA_ID}_*.fasta.gz ${WORK_DIR}/

# 3. 확장자 변경 (args_oap는 종종 .fa를 선호함)
cd ${WORK_DIR}
for f in *.fasta.gz; do
    mv "$f" "${f%.fasta.gz}.fa.gz"
done
cd ..

# 4. ARGs-OAP 실행
# -i: 입력 디렉토리, -o: 결과 디렉토리, -f: 확장자(fa), -t: 스레드
args_oap stage_one -i ${WORK_DIR} -o "${SRA_ID}_args_output" -f fa -t 8

if [ $? -eq 0 ]; then
    args_oap stage_two -i "${SRA_ID}_args_output" -t 8
else
    echo "오류: ARGs-OAP stage_one 실패"
    exit 1
fi

# 5. 후정리
rm -rf ${WORK_DIR}

echo "[A2] ARGs-OAP 분석 완료"
