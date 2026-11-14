#!/bin/bash
sample="$1"

# 1. 파일 이름 후보군 설정
if [ -f "${sample}_1.fastq.gz" ]; then
    fq1="${sample}_1.fastq.gz"
    is_paired=true
elif [ -f "${sample}.fastq.gz" ]; then
    fq1="${sample}.fastq.gz"
    is_paired=false
else
    echo "❗ Error: '${sample}_1.fastq.gz' 또는 '${sample}.fastq.gz' 파일이 존재하지 않습니다."
    exit 1
fi

fq2="${sample}_2.fastq.gz"

# 2. 변환 로직
echo "▶ Converting $fq1 to Fasta..."
# 분석의 일관성을 위해 Single-end라도 결과물은 _1.fasta.gz로 만들어줍니다.
seqkit fq2fa "$fq1" -j 4 | gzip > "${sample}_1.fasta.gz"

if [ "$is_paired" = true ] && [ -f "$fq2" ]; then
    echo "▶ Converting paired file $fq2 ..."
    seqkit fq2fa "$fq2" -j 4 | gzip > "${sample}_2.fasta.gz"
    echo "✅ Done: ${sample}_1.fasta.gz and ${sample}_2.fasta.gz"
else
    echo "⚠ Notice: Single-end 모드로 변환을 완료했습니다."
    echo "✅ Done: ${sample}_1.fasta.gz"
fi
