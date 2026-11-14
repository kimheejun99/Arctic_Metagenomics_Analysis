import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================
# top_taxa_barplot.py
# host_tracing_master.tsv 전체(특정 KO로 필터링하지 않음)에서
# read 개수가 가장 많은 상위 taxon들을 AL vs PERM 그룹별로 나눠
# 막대그래프로 보여줍니다.
#
# 사용법:
#   /root/miniforge3/envs/eggnog/bin/python top_taxa_barplot.py
# ==============================================================

IN_TSV = "host_tracing_master.tsv"
OUT_PNG = "top_taxa_by_group_barplot.png"
TOP_N = 20   # 10으로 바꾸고 싶으면 여기만 수정

TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"  # 상위계층(역/문/강 등) 보강용

# ---- 1) 읽기 (taxid=0 제외: Kraken2 unclassified) ----
df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
df = df[df["taxid"] != "0"]
print(f"[정보] 전체 레코드(unclassified 제외): {len(df)}개")

# ---- 2) 전체 read 수 기준 상위 N개 taxid 선정 ----
top_taxids = df["taxid"].value_counts().head(TOP_N).index.tolist()

# ---- 3) 그 taxid들에 대해서만 AL/PERM 개수 집계 ----
sub = df[df["taxid"].isin(top_taxids)]
pivot = sub.groupby(["taxid", "group"]).size().unstack(fill_value=0)
pivot = pivot.reindex(top_taxids)  # 원래 순위 순서 유지

for g in ["AL", "PERM"]:
    if g not in pivot.columns:
        pivot[g] = 0

print("\n[상위 taxid, AL/PERM 분포]")
print(pivot)

# ---- 4) 이름 매핑 ----
# ---- 4) 이름 매핑 ----
# 실제 형식: #TAXONOMY  Kingdom  Phylum  Class  Order  Family  Genus  Species
# (각 칸은 k__/p__/c__/o__/f__/g__/s__ 접두어가 붙어있음)
# taxid마다 "가장 구체적으로 채워진 칸"을 이름으로 사용 (species가 있으면 species,
# 없으면 genus, 그것도 없으면 그 위 단계...)
taxid_to_name = {}
try:
    with open(TAXID_NAME_FILE) as f:
        next(f)  # 헤더 줄(#TAXONOMY ...) 건너뜀
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            taxid = parts[0]
            name = None
            # Species(마지막 칸)부터 역순으로 훑어서 처음 채워진 값을 채택
            for val in reversed(parts[1:]):
                val = val.strip()
                # 접두어(k__/p__/c__/o__/f__/g__/s__)를 먼저 제거한 뒤 빈값 여부 판단
                stripped = val.split("__", 1)[-1] if "__" in val else val
                if stripped and stripped not in ("-", "NA"):
                    name = stripped
                    break
            if name:
                taxid_to_name[taxid] = name
    print(f"\n[정보] {TAXID_NAME_FILE} 에서 {len(taxid_to_name)}개 매핑 로드됨")
except FileNotFoundError:
    print(f"\n[안내] {TAXID_NAME_FILE} 을 못 찾아 taxid 숫자 그대로 표시합니다.")

# ---- 4b) 여전히 이름 없는 taxid는 names.dmp(NCBI 공식, 역/문/강 등 상위계층 포함)로 보강 ----
missing = [t for t in top_taxids if t not in taxid_to_name]
if missing:
    print(f"[정보] {len(missing)}개 taxid가 2_taxonomy_table.txt에 없어 names.dmp에서 추가 조회: {missing}")
    missing_set = set(missing)
    found = {}
    try:
        with open(NAMES_DMP) as f:
            for line in f:
                # 형식: tax_id\t|\tname_txt\t|\tunique_name\t|\tname_class\t|
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if taxid in missing_set and name_class == "scientific name":
                    found[taxid] = name
        taxid_to_name.update(found)
        print(f"[정보] names.dmp에서 {len(found)}개 추가로 찾음")
    except FileNotFoundError:
        print(f"[안내] {NAMES_DMP} 을 못 찾았습니다. NAMES_DMP 경로를 확인하세요.")

labels = []
for taxid in pivot.index:
    name = taxid_to_name.get(taxid)
    labels.append(f"{name}\n({taxid})" if name else f"taxid {taxid}")

# ---- 5) 그룹별 누적(stacked) 막대그래프 ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(13, 7))

x = range(len(pivot))
ax.bar(x, pivot["AL"], label="AL", color="#4C72B0", edgecolor="black", linewidth=0.5)
ax.bar(x, pivot["PERM"], bottom=pivot["AL"], label="PERM", color="#DD8452", edgecolor="black", linewidth=0.5)

# 총합 라벨
totals = pivot["AL"] + pivot["PERM"]
for i, v in enumerate(totals):
    ax.text(i, v + max(totals) * 0.01, str(v), ha="center", fontsize=9, fontweight="bold")

ax.set_xticks(list(x))
ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("Read count (host-tracing joined)")
ax.set_title(f"Top {TOP_N} taxa by read count in host-tracing results (AL vs PERM)")
ax.legend(title="Group")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"\n완료: {OUT_PNG}")
