#!/bin/bash

#!/bin/bash
SRA_ID=$1

# 해당 ID로 시작하는 모든 fastq.gz 파일을 찾아 배열에 담음
FILES=$(ls ${SRA_ID}*.fastq.gz 2>/dev/null)

if [ -z "$FILES" ]; then
    echo "Error: No fastq files found for ${SRA_ID}"
    exit 1
fi

# 찾은 모든 파일(1개 혹은 2개 이상)에 대해 통계 생성
seqkit stats $FILES > "${SRA_ID}_count.txt"
