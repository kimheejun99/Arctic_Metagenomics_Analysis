import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ==============================================================
# visualize_robust_pathways_taxon.py
# batch_pathway_to_taxon.py가 만든 ancombc2_robust_pathways_taxon_summary.csv를
# 시각화합니다. 두 개의 축(총 read 수, 1위 taxon 집중도%)으로 각 pathway가
# "한 종에 집중"됐는지 "수백 종에 분산"됐는지, 그리고 표본 자체가 너무
# 적어(노이즈) 신뢰하기 어려운지 한눈에 보여줍니다.
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python visualize_robust_pathways_taxon.py
# ==============================================================

IN_CSV = "ancombc2_robust_pathways_taxon_summary.csv"
OUT_PNG = "robust_pathways_taxon_pattern.png"
LOW_READ_THRESHOLD = 20  # 이 개수 미만이면 "표본 부족, 노이즈 의심"으로 강조

df = pd.read_csv(IN_CSV)
df = df.sort_values("total_reads", ascending=True)

sns.set_style("whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [1.3, 1]})

# ---- 왼쪽: 총 read 수 (log scale) - 노이즈 의심 항목 강조 ----
colors1 = ["#C0392B" if r < LOW_READ_THRESHOLD else "#4C72B0" for r in df["total_reads"]]
bars1 = ax1.barh(df["pathway_name"], df["total_reads"].clip(lower=0.5), color=colors1,
                  edgecolor="black", linewidth=0.5)
ax1.set_xscale("log")
ax1.axvline(LOW_READ_THRESHOLD, color="gray", linestyle="--", linewidth=1)
ax1.set_xlabel("Total reads (log scale)")
ax1.set_title(f"Total reads per pathway\n(red = under {LOW_READ_THRESHOLD}, likely noise)")
for i, (v, tot) in enumerate(zip(df["total_reads"], df["total_reads"])):
    ax1.text(max(v, 0.5) * 1.15, i, str(tot), va="center", fontsize=8)

# ---- 오른쪽: 1위 taxon 집중도(%) vs 전체 taxon 수 ----
scatter_colors = ["#C0392B" if r < LOW_READ_THRESHOLD else "#4C72B0" for r in df["total_reads"]]
sizes = np.clip(df["total_reads"] / 5, 20, 400)
sc = ax2.scatter(df["n_distinct_taxa"], df["top1_pct"], s=sizes, c=scatter_colors,
                  edgecolor="black", linewidth=0.6, alpha=0.8)
for _, row in df.iterrows():
    label = row["top1_taxon"] if row["total_reads"] >= LOW_READ_THRESHOLD else f"{row['top1_taxon']} (n={row['total_reads']})"
    ax2.annotate(label, (row["n_distinct_taxa"], row["top1_pct"]),
                 fontsize=7, xytext=(5, 3), textcoords="offset points")
ax2.set_xlabel("Number of distinct taxa detected")
ax2.set_ylabel("Top-1 taxon share (%)")
ax2.set_title("Concentration (high %, few taxa) vs\nDispersion (low %, many taxa)")
ax2.set_ylim(0, 105)

from matplotlib.patches import Patch
fig.legend(handles=[
    Patch(facecolor="#4C72B0", label=f"reads >= {LOW_READ_THRESHOLD} (해석 가능)"),
    Patch(facecolor="#C0392B", label=f"reads < {LOW_READ_THRESHOLD} (표본부족, 노이즈 의심)"),
], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
print(f"완료: {OUT_PNG}")
