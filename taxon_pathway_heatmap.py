import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================
# taxon_pathway_heatmap.py
# host_tracing_master.tsv를 "taxon x pathway" 매트릭스로 재구성합니다.
# 지금까지의 그래프(taxon 순위만, 또는 기능 하나만)와 달리,
# "이 taxon이 어떤 기능들을 얼마나(%) 갖고 있는지"를 한 번에 보여줍니다.
#
# 값 = 그 taxon의 전체 read 중, 그 pathway를 가진 read의 비율(%)
#      (한 read가 여러 pathway를 가질 수 있어 가로 합이 100%는 아님)
#
# 사용법:
#   /root/miniforge3/envs/eggnog/bin/python taxon_pathway_heatmap.py
# ==============================================================

IN_TSV = "host_tracing_master.tsv"
OUT_PNG = "taxon_pathway_heatmap.png"

TOP_N_TAXA = 15
TOP_M_PATHWAYS = 20

TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"
PATHWAY_NAME_FILE = "kegg_pathway_names_clean.tsv"  # curl로 받아둔 KEGG 공식 이름 목록

# Global and overview maps (총괄 카테고리, 예전과 동일하게 제외)
EXCLUDE_PATHWAYS = {
    "map01100", "map01110", "map01120", "map01130", "map01200", "map01210",
    "map01212", "map01220", "map01230", "map01232", "map01240", "map01250",
    "map01310", "map01320",
    "01100", "01110", "01120", "01130", "01200", "01210",
    "01212", "01220", "01230", "01232", "01240", "01250", "01310", "01320",
}

# ---- 1) 데이터 읽기, taxid=0 제외 ----
df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
df = df[df["taxid"] != "0"]
df = df[df["KEGG_Pathway"].notna() & (df["KEGG_Pathway"] != "-")]
print(f"[정보] pathway 정보 있는 레코드: {len(df)}개")

# ---- 2) 상위 taxon 선정 (전체 read 수 기준, explode 전에) ----
top_taxa = df["taxid"].value_counts().head(TOP_N_TAXA).index.tolist()
df_top_taxa = df[df["taxid"].isin(top_taxa)].copy()

# ---- 3) KEGG_Pathway 컬럼 분리(explode) - 한 read가 여러 pathway를 가질 수 있음 ----
df_top_taxa["KEGG_Pathway"] = df_top_taxa["KEGG_Pathway"].str.split(",")
exploded = df_top_taxa.explode("KEGG_Pathway")
exploded["KEGG_Pathway"] = exploded["KEGG_Pathway"].str.strip()
exploded = exploded[~exploded["KEGG_Pathway"].isin(EXCLUDE_PATHWAYS)]
exploded = exploded[exploded["KEGG_Pathway"].str.startswith("map") | exploded["KEGG_Pathway"].str.match(r"^\d+$")]

# ---- 4) 상위 pathway 선정 (총괄 카테고리 제외 후, 상위 taxon들 내에서의 빈도 기준) ----
top_pathways = exploded["KEGG_Pathway"].value_counts().head(TOP_M_PATHWAYS).index.tolist()
print(f"[정보] 상위 taxon {len(top_taxa)}개, 상위 pathway {len(top_pathways)}개")

# ---- 5) 매트릭스 생성: taxon별 read 수(분모) 및 taxon x pathway count(분자) ----
taxon_total_reads = df_top_taxa.groupby("taxid")["read_id"].nunique()

sub = exploded[exploded["KEGG_Pathway"].isin(top_pathways)]
matrix_count = sub.groupby(["taxid", "KEGG_Pathway"])["read_id"].nunique().unstack(fill_value=0)
matrix_count = matrix_count.reindex(index=top_taxa, columns=top_pathways, fill_value=0)

matrix_pct = matrix_count.div(taxon_total_reads.reindex(top_taxa), axis=0) * 100

# ---- 6) taxid 이름 매핑 (2단계: taxonomy_table.txt -> names.dmp) ----
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

missing = [t for t in top_taxa if t not in taxid_to_name]
if missing:
    missing_set = set(missing)
    try:
        with open(NAMES_DMP) as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if taxid in missing_set and name_class == "scientific name":
                    taxid_to_name[taxid] = name
    except FileNotFoundError:
        pass

row_labels = [f"{taxid_to_name.get(t, 'taxid '+t)} ({t})" for t in top_taxa]

# ---- 7) pathway 이름 매핑 ----
pathway_to_name = {}
try:
    with open(PATHWAY_NAME_FILE) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                pathway_to_name[parts[0]] = parts[1]
except FileNotFoundError:
    print(f"[안내] {PATHWAY_NAME_FILE} 을 못 찾아 pathway ID 그대로 표시합니다.")

def pw_label(pid):
    key = pid if pid.startswith("map") else f"map{pid}"
    name = pathway_to_name.get(key, pathway_to_name.get(pid))
    return name if name else pid

col_labels = [pw_label(p) for p in top_pathways]

# ---- 8) 히트맵 ----
sns.set_style("white")
fig, ax = plt.subplots(figsize=(max(14, len(top_pathways) * 0.6), max(8, len(top_taxa) * 0.5)))

sns.heatmap(matrix_pct.values, annot=True, fmt=".0f", cmap="YlOrRd",
            xticklabels=col_labels, yticklabels=row_labels,
            cbar_kws={"label": "% of taxon's reads with this pathway"}, ax=ax)

ax.set_title(f"Taxon x KEGG Pathway matrix (top {TOP_N_TAXA} taxa, top {TOP_M_PATHWAYS} pathways)")
plt.xticks(rotation=60, ha="right", fontsize=8)
plt.yticks(fontsize=9)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)

print(f"\n완료: {OUT_PNG}")
print("\n[매트릭스 미리보기 (raw count)]")
print(matrix_count)
