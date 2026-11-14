import pandas as pd
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("사용법: python category_summary.py <SampleID>")
        return

    sample_id = sys.argv[1]
    input_file = f"{sample_id}_annotated_parsed.tsv"

    if not os.path.exists(input_file):
        print(f"❌ 파일이 존재하지 않습니다: {input_file}")
        return

    # 파일 읽기 (TSV 형식)
    df = pd.read_csv(input_file, sep="\t")

    # 'Category'별 그룹화 후 Hit Count 계산
    category_summary = df.groupby("Category")["SampleID"].count().reset_index(name="Hit Count")

    # 전체 Hit 수 계산
    total_hits = category_summary["Hit Count"].sum()

    # 퍼센트 계산
    category_summary["Percentage"] = (category_summary["Hit Count"] / total_hits * 100).round(2)

    # 정렬
    category_summary = category_summary.sort_values(by="Hit Count", ascending=False)

    # 결과 저장 (TSV)
    output_file = f"{sample_id}_category_summary.tsv"
    category_summary.to_csv(output_file, sep="\t", index=False)

    # 출력
    print(category_summary)
    print(f"✅ 저장 완료: {output_file}")

if __name__ == "__main__":
    main()

