import sys
import os
import pandas as pd
from Bio import SeqIO

def extract_vfg_id(entry):
    """VFG ID만 추출 (예: 'VFG037170(gb|WP_001081754)' → 'VFG037170')"""
    return entry.split('(')[0]

def load_fasta_headers(fasta_file):
    """FASTA 파일에서 VFG ID → 헤더 전체 매핑 딕셔너리 생성"""
    header_map = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        vfg_id = record.id.split('(')[0]
        header_map[vfg_id] = record.description
    return header_map

def main():
    if len(sys.argv) != 3:
        print("사용법: python blastx_tophit_annotator.py <diamond_output.tsv> <vfdb_fasta>")
        return

    input_file = sys.argv[1] + "_output.tsv"
    fasta_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"❌ 입력 파일이 존재하지 않습니다: {input_file}")
        return
    if not os.path.exists(fasta_file):
        print(f"❌ FASTA 파일이 존재하지 않습니다: {fasta_file}")
        return

    try:
        df = pd.read_csv(input_file, sep='\t', header=None)
    except Exception as e:
        print(f"❌ 파일 읽기 오류: {e}")
        return

    # 각 read에 대한 top hit 추출
    top_hits = df.groupby(0).first().reset_index()

    # VFG ID 추출
    top_hits['VFG_ID'] = top_hits[1].apply(extract_vfg_id)

    # VFDB FASTA 헤더 로딩
    vfg_headers = load_fasta_headers(fasta_file)

    # 주석 붙이기
    top_hits['Description'] = top_hits['VFG_ID'].map(vfg_headers)

    # 출력 파일명
    base = os.path.splitext(input_file)[0]
    output_file = sys.argv[1] + "_annotated.tsv"

    # 저장
    top_hits.drop(columns=["VFG_ID"]).to_csv(output_file, sep='\t', header=False, index=False)
    print(f"✅ 완료! 주석 결과 저장됨: {output_file}")

if __name__ == "__main__":
    main()

