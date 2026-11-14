import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================
# host_tracing_barplot.py
# host_tracing_master.tsv에서 메탄생성 핵심유전자(K00399/K00401/K00402,
# methyl-coenzyme M reductase)를 가진 read들이 어느 taxon에서 왔는지
# 막대그래프로 시각화합니다. Methanosarcina(taxid 2207)는 강조 색상으로
# 표시됩니다.
#
# 사용법:
#   /root/miniforge3/envs/eggnog/bin/python host_tracing_barplot.py
#
# 선택: taxonomy_table.txt 등 taxid->이름 매핑 파일이 있으면
#       TAXID_NAME_FILE 경로를 지정하면 막대 라벨이 실제 이름으로 나옵니다.
#       형식 가정: tab 구분, 컬럼에 taxid와 name이 있음 (아래 설정 확인)
# ==============================================================

IN_TSV = "host_tracing_master.tsv"
OUT_PNG = "methanogenesis_host_tracing_barplot.png"

TARGET_KOS = ["K00399", "K00401", "K00402"]  # methyl-coenzyme M reductase
TOP_N = 10
HIGHLIGHT_TAXID = "2207"  # Methanosarcina

# taxid -> 이름 매핑 파일 (선택, 없으면 taxid 숫자 그대로 표시)
TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"

# ---- 1) host_tracing_master.tsv 읽고 필터링 ----
df = pd.read_csv(IN_TSV, sep="\t", dtype=str)

pattern = "|".join(TARGET_KOS)
mask = df["KO"].str.contains(pattern, na=False) & (df["taxid"] != "0")
filtered = df[mask]

print(f"[정보] mcr 유전자(K00399/401/402) 매칭 read: {len(filtered)}개")

# ---- 2) taxid별 카운트, 상위 N개 ----
counts = filtered["taxid"].value_counts().head(TOP_N)
print("\n[상위 taxid]")
print(counts)

# ---- 3) 이름 매핑 시도 (파일 있으면) ----
# ---- 3) 이름 매핑 시도 (실제 형식: #TAXONOMY Kingdom~Species, k__/p__/... 접두어) ----
taxid_to_name = {}
try:
    with open(TAXID_NAME_FILE) as f:
        next(f)  # 헤더 줄 건너뜀
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
    print(f"\n[정보] {TAXID_NAME_FILE} 에서 {len(taxid_to_name)}개 taxid-이름 매핑 로드됨")
except FileNotFoundError:
    print(f"\n[안내] {TAXID_NAME_FILE} 을 못 찾아 taxid 숫자 그대로 표시합니다.")

# ---- 3b) 여전히 이름 없는 taxid는 names.dmp로 보강 ----
missing = [t for t in counts.index if t not in taxid_to_name]
if missing:
    missing_set = set(missing)
    found = {}
    try:
        with open(NAMES_DMP) as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if taxid in missing_set and name_class == "scientific name":
                    found[taxid] = name
        taxid_to_name.update(found)
        print(f"[정보] names.dmp에서 {len(found)}개 추가로 찾음")
    except FileNotFoundError:
        print(f"[안내] {NAMES_DMP} 을 못 찾았습니다.")

labels = []
for taxid in counts.index:
    name = taxid_to_name.get(taxid)
    if name:
        labels.append(f"{name}\n(taxid {taxid})")
    else:
        labels.append(f"taxid {taxid}\n(unidentified)")

# ---- 4) 막대그래프 (Methanosarcina 강조) ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(11, 6))

colors = ["#C0392B" if taxid == HIGHLIGHT_TAXID else "#8DA0CB" for taxid in counts.index]

bars = ax.bar(range(len(counts)), counts.values, color=colors, edgecolor="black", linewidth=0.5)

# 값 라벨 표시
for i, v in enumerate(counts.values):
    ax.text(i, v + 0.15, str(v), ha="center", fontsize=10, fontweight="bold")

ax.set_xticks(range(len(counts)))
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Read count (K00399/K00401/K00402 match)")
ax.set_title("Taxonomic origin of methyl-coenzyme M reductase reads\n(methanogenesis marker genes)")

# 범례: 강조색 설명
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#C0392B", edgecolor="black", label="Methanosarcina (methanogen)"),
    Patch(facecolor="#8DA0CB", edgecolor="black", label="Other taxa"),
]
ax.legend(handles=legend_elements, loc="upper right")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"\n완료: {OUT_PNG}")
