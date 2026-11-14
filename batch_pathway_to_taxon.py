import pandas as pd

# ==============================================================
# batch_pathway_to_taxon.py
# ANCOM-BC2 표본수 균형 재검정(A2)에서 75% 이상 강건하게 유의했던
# 17개 pathway를 한 번에 확인하여, 각각 1위 taxon이 무엇인지 요약합니다.
# (pathway_to_taxon_lookup.py를 17번 따로 실행할 필요 없이 한번에 처리)
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python batch_pathway_to_taxon.py
# ==============================================================

IN_TSV = "host_tracing_master.tsv"
TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"

# A2_ancombc2_balanced_subsample.R 결과에서 75% 이상(15/20회+) 유의했던 17개
TARGET_PATHWAYS = {
    "00909": "Sesquiterpenoid and triterpenoid biosynthesis",
    "04744": "Phototransduction",
    "01062": "Biosynthesis of terpenoids and steroids",
    "00540": "Lipopolysaccharide biosynthesis",
    "05225": "Hepatocellular carcinoma",
    "00623": "Toluene degradation",
    "00905": "Brassinosteroid biosynthesis",
    "00073": "Cutin, suberine and wax biosynthesis",
    "00364": "Fluorobenzoate degradation",
    "00791": "Atrazine degradation",
    "00361": "Chlorocyclohexane and chlorobenzene degradation",
    "00622": "Xylene degradation",
    "00980": "Metabolism of xenobiotics by cytochrome P450",
    "00982": "Drug metabolism - cytochrome P450",
    "04630": "JAK-STAT signaling pathway",
    "04919": "Thyroid hormone signaling pathway",
    "05204": "Chemical carcinogenesis - DNA adducts",
}

# ---- 1) taxid 이름 매핑 ----
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
    print(f"[안내] {TAXID_NAME_FILE} 없음")

# ---- 2) 데이터 한 번만 읽기 ----
df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
df = df[df["taxid"] != "0"]
df = df[df["KEGG_Pathway"].notna() & (df["KEGG_Pathway"] != "-")]
df["pw_list"] = df["KEGG_Pathway"].apply(lambda s: [p.strip().replace("map", "") for p in s.split(",")])

# names.dmp는 필요할 때만 한 번에 로드 (성능)
all_taxids = set(df["taxid"].unique())
missing = all_taxids - set(taxid_to_name.keys())
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
        print(f"[안내] {NAMES_DMP} 없음")

# ---- 3) 17개 pathway 각각 요약 ----
results = []
for pw_id, pw_name in TARGET_PATHWAYS.items():
    mask = df["pw_list"].apply(lambda lst: pw_id in lst)
    sub = df[mask]
    if len(sub) == 0:
        results.append({"pathway_id": pw_id, "pathway_name": pw_name,
                         "total_reads": 0, "top1_taxon": "(없음)", "top1_pct": 0,
                         "n_distinct_taxa": 0})
        continue
    counts = sub["taxid"].value_counts()
    top1_taxid = counts.index[0]
    top1_name = taxid_to_name.get(top1_taxid, f"taxid_{top1_taxid}")
    top1_pct = round(100 * counts.iloc[0] / counts.sum(), 1)
    results.append({
        "pathway_id": pw_id, "pathway_name": pw_name,
        "total_reads": len(sub), "top1_taxon": top1_name, "top1_pct": top1_pct,
        "n_distinct_taxa": sub["taxid"].nunique(),
    })

result_df = pd.DataFrame(results)
result_df.to_csv("ancombc2_robust_pathways_taxon_summary.csv", index=False)

print(f"{'Pathway':<55} {'총read':>7} {'1위 taxon':<30} {'1위%':>6} {'taxon수':>7}")
print("-" * 115)
for _, r in result_df.iterrows():
    print(f"{r['pathway_name']:<55} {r['total_reads']:>7} {r['top1_taxon']:<30} {r['top1_pct']:>5}% {r['n_distinct_taxa']:>7}")

print(f"\n완료: ancombc2_robust_pathways_taxon_summary.csv")
