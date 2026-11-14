import pandas as pd
import re
import sys
import os

def extract_info(annotation):
    # VF_Name
    vf_name_match = re.search(r"\[(.*?) \(VF(\d+)\)", annotation)
    vf_name = "VF" + vf_name_match.group(2) if vf_name_match else ""
    
    # VFC_ID
    vfc_id_match = re.search(r"\(VFC(\d+)\)", annotation)
    vfc_id = "VFC" + vfc_id_match.group(1) if vfc_id_match else ""

    # Category
    category_match = re.search(r"- (.*?) \(VFC", annotation)
    category = category_match.group(1).strip() if category_match else ""

    # Gene
    gene_match = re.search(r"\((\w+)\)", annotation)
    gene = gene_match.group(1) if gene_match else ""

    # Function (괄호 안의 유전자명 제거)
    function_match = re.search(r"\)\s+([^\[]+?)\s+\[", annotation)
    if function_match:
        raw_function = function_match.group(1).strip()
        function = re.sub(r"\s*\(.*?\)", "", raw_function).strip()
    else:
        function = ""

    # Strain
    strain_match = re.findall(r"\[(.*?)\]", annotation)
    strain = strain_match[-1] if strain_match else ""

    return pd.Series([vf_name, vfc_id, category, gene, function, strain],
                     index=["VF_Name", "VFC_ID", "Category", "Gene", "Function", "Strain"])

def main():
    if len(sys.argv) != 2:
        print("사용법: python parse_vf_info.py <input_output.txt>")
        return

    input_file = sys.argv[1] + "_annotated.tsv"
    if not os.path.exists(input_file):
        print(f"❌ 파일이 존재하지 않습니다: {input_file}")
        return

    # TSV 파일 읽기
    df = pd.read_csv(input_file, sep="\t", header=None)

    # 컬럼 이름 설정 및 필요한 컬럼 추출
    df.columns = ["SampleID", "VF_ID", "Identity", "Align_Length", "Mismatch", "Gap",
                  "Q_start", "Q_end", "S_start", "S_end", "Evalue", "BitScore", "Annotation"]
    df = df[["SampleID", "Annotation"]]

    # 정보 추출
    info_df = df["Annotation"].apply(extract_info)

    # 합치기
    final_df = pd.concat([df, info_df], axis=1)

    # 출력 파일 이름 생성
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}_parsed.tsv"

    # 저장
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"✅ 완료! 결과 파일: {output_file}")

if __name__ == "__main__":
    main()

