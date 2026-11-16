import pandas as pd

# ==============================================================
# host_tracing_normalized.py
# methanogenesis_host_tracing_barplot.py의 raw count 결과가 샘플별
# 시퀀싱 깊이 차이(특히 PERM이 평균 2.4배 더 많은 read를 가짐) 때문에
# 왜곡되지 않았는지, RPM(Reads Per Million) 정규화로 재확인합니다.
#
# 입력: host_tracing_master.tsv, sequencing_depth_per_sample.csv
#       (check_sequencing_depth.py가 만든 파일을 그대로 사용)
# 출력: raw count 순위 vs RPM 정규화 순위를 나란히 비교
#
# 사용법: python3 host_tracing_normalized.py
#   (host_tracing_master.tsv, sequencing_depth_per_sample.csv,
#    ../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt 가
#    모두 접근 가능한 위치에서 실행)
# ==============================================================

HOST_TRACING = "host_tracing_master.tsv"
DEPTH_CSV = "sequencing_depth_per_sample.csv"
TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"
OUT_PNG = "host_tracing_normalized_comparison.png"

MCR_KOS = ["K00399", "K00401", "K00402"]
TOP_N = 10

# ---- 1) 시퀀싱 깊이 정보 읽기 ----
depth = pd.read_csv(DEPTH_CSV)
depth_map = dict(zip(depth["sample_id"], depth["total_reads"]))
print(f"[정보] 시퀀싱 깊이 정보 있는 샘플: {len(depth_map)}개")

# ---- 2) host_tracing_master.tsv에서 mcr 유전자만 필터링 ----
df = pd.read_csv(HOST_TRACING, sep="\t", dtype=str)
df = df[df["taxid"] != "0"]
df = df[df["KO"].isin(MCR_KOS)]
print(f"[정보] mcr 유전자(K00399/401/402) 가진 read: {len(df)}개")

# ---- 3) 샘플별 총 read 수 매칭 (sample 컬럼이 샘플ID를 포함하는 형태일 수 있어 부분매칭) ----
def find_depth(sample_col_value):
    for srr, total in depth_map.items():
        if srr in str(sample_col_value):
            return total
    return None

df["total_reads_in_sample"] = df["sample"].apply(find_depth)
unmatched = df[df["total_reads_in_sample"].isna()]["sample"].unique()
if len(unmatched) > 0:
    print(f"[안내] 깊이 정보를 못 찾은 샘플: {list(unmatched)[:5]}")
df = df.dropna(subset=["total_reads_in_sample"])

# ---- 4) taxid x sample 별 read 개수 -> RPM 계산 ----
counts = df.groupby(["taxid", "sample", "total_reads_in_sample"]).size().reset_index(name="read_count")
counts["rpm"] = counts["read_count"] / counts["total_reads_in_sample"] * 1_000_000

# ---- 5) taxid별 합계: raw count 합 vs RPM 합 ----
raw_summary = df.groupby("taxid").size().rename("raw_count")
rpm_summary = counts.groupby("taxid")["rpm"].sum().rename("total_rpm")
summary = pd.concat([raw_summary, rpm_summary], axis=1).fillna(0)

# ---- 6) taxid 이름 매핑 ----
taxid_to_name = {}
try:
    with open(TAXID_NAME_FILE) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            taxid = parts[0]
            name = None
            for val in reversed(parts[1:]):
                val = val.strip()
                stripped = val.split("__", 1)[-1] if "__" in val else val
                if stripped and stripped not in ("-", "NA"):
                    name = stripped
                    break
            if name:
                taxid_to_name[taxid] = name
except FileNotFoundError:
    pass

missing = [t for t in summary.index if t not in taxid_to_name]
if missing:
    try:
        with open(NAMES_DMP) as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if taxid in missing and name_class == "scientific name":
                    taxid_to_name[taxid] = name
    except FileNotFoundError:
        pass

summary["name"] = summary.index.map(lambda t: taxid_to_name.get(t, f"taxid_{t}"))

# ---- 7) 두 기준으로 각각 순위 매겨서 비교 ----
by_raw = summary.sort_values("raw_count", ascending=False).head(TOP_N)
by_rpm = summary.sort_values("total_rpm", ascending=False).head(TOP_N)

print(f"\n=== Raw count 기준 상위 {TOP_N} ===")
print(by_raw[["name", "raw_count", "total_rpm"]].to_string())

print(f"\n=== RPM(깊이 정규화) 기준 상위 {TOP_N} ===")
print(by_rpm[["name", "raw_count", "total_rpm"]].to_string())

raw_top_taxa = set(by_raw.index)
rpm_top_taxa = set(by_rpm.index)
same_rank = raw_top_taxa == rpm_top_taxa
print(f"\n[핵심] Raw와 RPM top {TOP_N} taxon 구성이 동일한가: {same_rank}")
if not same_rank:
    print(f"  Raw에만 있음: {[taxid_to_name.get(t,t) for t in raw_top_taxa - rpm_top_taxa]}")
    print(f"  RPM에만 있음: {[taxid_to_name.get(t,t) for t in rpm_top_taxa - raw_top_taxa]}")

# ---- 8) 시각화: 나란히 비교 ----
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.barh(by_raw["name"][::-1], by_raw["raw_count"][::-1], color="#4C72B0", edgecolor="black")
ax1.set_title("Ranked by raw read count")
ax1.set_xlabel("Read count")

ax2.barh(by_rpm["name"][::-1], by_rpm["total_rpm"][::-1], color="#DD8452", edgecolor="black")
ax2.set_title("Ranked by RPM (depth-normalized)")
ax2.set_xlabel("Reads per million")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"\n완료: {OUT_PNG}")
