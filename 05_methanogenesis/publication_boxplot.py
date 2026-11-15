import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ==============================================================
# publication_boxplot.py
# KEGGaNOG의 merged_pathways.tsv에서 관심 pathway 몇 개만 골라
# AL vs PERM 그룹 boxplot(+개별 샘플 점)으로 비교. 논문 Figure용.
#
# 49개 샘플을 전부 펼쳐놓는 히트맵보다 훨씬 읽기 쉽고,
# 이 분야 논문에서 가장 흔히 쓰이는 2그룹 비교 형식입니다.
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python publication_boxplot.py
# ==============================================================

IN_TSV = "kegganog_results/merged_pathways.tsv"
OUT_PNG = "kegganog_results/methanogenesis_boxplot.png"

# ==============================================================
# 샘플 목록을 하드코딩하지 않고, metadata.txt(sample-id	group)에서
# 동적으로 읽어옵니다. 몇 개 샘플이 있든(2개든 49개든) 자동으로 맞춰집니다.
# metadata.txt는 이 스크립트를 실행하는 폴더(eggnog_analysis/)의
# 상위 폴더에 있다고 가정합니다 (run_pipeline.sh 구조와 동일).
# ==============================================================
import pandas as pd_meta_loader

_META_PATH = "../metadata.txt"
_meta_df = pd_meta_loader.read_csv(_META_PATH, sep="	")
_meta_df.columns = [c.strip() for c in _meta_df.columns]
AL_SAMPLES = _meta_df[_meta_df["group"] == "AL"]["sample-id"].astype(str).tolist()
PERM_SAMPLES = _meta_df[_meta_df["group"] == "PERM"]["sample-id"].astype(str).tolist()
print(f"[metadata.txt] AL {len(AL_SAMPLES)}개, PERM {len(PERM_SAMPLES)}개 로드됨")

# 히트맵에서 시각적으로 두드러졌던 pathway들
TARGET_PATHWAYS = [
    "Methanogenesis via CO2",
    "Methanogenesis via acetate",
    "Methanogenesis via dimethylamine",
    "Methanogenesis via methanol",
    "Methanogenesis via methylamine",
    "Methanogenesis via trimethylamine",
    "Methanol dehydrogenase",
    "Coenzyme M reduction to methane",
]

# ---- 1) 데이터 읽기 + 컬럼 자동 매칭 ----
df = pd.read_csv(IN_TSV, sep="\t")
pathway_col = df.columns[0]

def find_col(srr_id, columns):
    matches = [c for c in columns if srr_id in c]
    if len(matches) != 1:
        raise ValueError(f"{srr_id}: {len(matches)}개 매칭")
    return matches[0]

al_cols = [find_col(s, df.columns) for s in AL_SAMPLES]
perm_cols = [find_col(s, df.columns) for s in PERM_SAMPLES]

# ---- 2) long-format으로 변환 (pathway, group, completeness) ----
df_sub = df[df[pathway_col].isin(TARGET_PATHWAYS)].set_index(pathway_col)

records = []
for pw in TARGET_PATHWAYS:
    if pw not in df_sub.index:
        print(f"[skip] '{pw}' 가 데이터에 없어 건너뜁니다.")
        continue
    for c in al_cols:
        records.append({"Pathway": pw, "Group": "AL", "Completeness": df_sub.loc[pw, c]})
    for c in perm_cols:
        records.append({"Pathway": pw, "Group": "PERM", "Completeness": df_sub.loc[pw, c]})

long_df = pd.DataFrame(records)

# ---- 3) 그룹별 Mann-Whitney U test (비모수, 표본 작을 때 안전) ----
pvals = {}
for pw in TARGET_PATHWAYS:
    if pw not in df_sub.index:
        continue
    al_vals = df_sub.loc[pw, al_cols].values
    perm_vals = df_sub.loc[pw, perm_cols].values
    try:
        _, p = stats.mannwhitneyu(al_vals, perm_vals, alternative="two-sided")
    except ValueError:
        p = float("nan")
    pvals[pw] = p

# ---- 4) Boxplot + 개별 샘플 점(strip) 그리기 ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(10, 6))

sns.boxplot(data=long_df, x="Pathway", y="Completeness", hue="Group",
            palette={"AL": "#4C72B0", "PERM": "#DD8452"}, ax=ax,
            showfliers=False, width=0.6)
sns.stripplot(data=long_df, x="Pathway", y="Completeness", hue="Group",
              palette={"AL": "#2c4a70", "PERM": "#a05a2c"}, ax=ax,
              dodge=True, size=4, alpha=0.6, legend=False)

# p-value를 각 pathway 위에 표시
for i, pw in enumerate(TARGET_PATHWAYS):
    if pw in pvals and not pd.isna(pvals[pw]):
        p = pvals[pw]
        label = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
        ax.text(i, 1.05, label, ha="center", fontsize=9,
                fontweight="bold" if p < 0.05 else "normal")

ax.set_ylim(-0.05, 1.2)
ax.set_xticks(range(len(TARGET_PATHWAYS)))
ax.set_xticklabels(TARGET_PATHWAYS, rotation=35, ha="right")
ax.set_ylabel("KEGG Module Completeness")
ax.set_xlabel("")
ax.legend(title="Group", loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)

print(f"완료: {OUT_PNG}")
print("\n[Mann-Whitney U test 결과]")
for pw, p in pvals.items():
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {pw}: p={p:.4f} {sig}")
