import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================
# ancombc2_volcano_plot.py
# ancombc2_pathway_results.csv를 volcano plot으로 시각화합니다.
# x축: log fold change (양수=PERM에서 높음, 음수=AL에서 높음 -
#      실제 방향은 A1_ancombc2_pathway.R의 참조군 설정에 따름, 확인 필요)
# y축: -log10(q-value)
#
# 사용법:
#   /root/miniforge3/envs/r_ancombc/bin/python ancombc2_volcano_plot.py
#   (matplotlib/seaborn/pandas가 있는 환경이면 어디서든 가능;
#    r_ancombc 환경에 없으면 eggnog 환경에서 실행해도 무방)
# ==============================================================

IN_CSV = "ancombc2_pathway_results.csv"
OUT_PNG = "ancombc2_volcano_plot.png"

Q_THRESHOLD = 0.05
N_LABELS = 15          # 이름을 표시할 상위 유의 pathway 개수
HIGHLIGHT_PATHWAYS = {  # 특별히 강조해서 표시할 pathway (색을 다르게)
    "map00680": "Methane metabolism",
}

# ---- 1) 데이터 읽기 ----
df = pd.read_csv(IN_CSV)
df.columns = [c.strip().lower() for c in df.columns]
# 컬럼명이 pathway_id, log_fold_change, q_value, pathway_name 이라고 가정
required = {"pathway_id", "log_fold_change", "q_value", "pathway_name"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"필요한 컬럼이 없습니다: {missing}. 실제 컬럼: {list(df.columns)}")

df["q_value"] = df["q_value"].clip(lower=1e-300)  # log10(0) 방지
df["neg_log10_q"] = -np.log10(df["q_value"])

# ---- 2) 유의성 분류 ----
def classify(row):
    if row["pathway_id"] in HIGHLIGHT_PATHWAYS:
        return "highlight"
    if row["q_value"] < Q_THRESHOLD:
        return "significant"
    return "not significant"

df["category"] = df.apply(classify, axis=1)
n_sig = (df["q_value"] < Q_THRESHOLD).sum()
print(f"[정보] 전체 {len(df)}개 pathway 중 q<{Q_THRESHOLD} 유의: {n_sig}개")

# ---- 3) Volcano plot ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(11, 8))

palette = {"not significant": "#B0B0B0", "significant": "#4C72B0", "highlight": "#C0392B"}
sizes = {"not significant": 20, "significant": 30, "highlight": 90}

for cat in ["not significant", "significant", "highlight"]:
    sub = df[df["category"] == cat]
    ax.scatter(sub["log_fold_change"], sub["neg_log10_q"],
               c=palette[cat], s=sizes[cat],
               label=cat, alpha=0.7 if cat == "not significant" else 0.9,
               edgecolors="black" if cat == "highlight" else "none",
               linewidths=1.2 if cat == "highlight" else 0, zorder=3 if cat == "highlight" else 2)

# 유의성 기준선
ax.axhline(-np.log10(Q_THRESHOLD), color="gray", linestyle="--", linewidth=1)
ax.axvline(0, color="gray", linestyle="-", linewidth=0.8)

# 상위 유의 pathway 이름 라벨 (q-value 기준 상위 N개)
top_sig = df[df["q_value"] < Q_THRESHOLD].nsmallest(N_LABELS, "q_value")
for _, row in top_sig.iterrows():
    ax.annotate(row["pathway_name"], (row["log_fold_change"], row["neg_log10_q"]),
                fontsize=7, xytext=(4, 4), textcoords="offset points")

# 강조 pathway 라벨 (유의 여부와 무관하게 항상 표시)
for _, row in df[df["category"] == "highlight"].iterrows():
    ax.annotate(row["pathway_name"], (row["log_fold_change"], row["neg_log10_q"]),
                fontsize=10, fontweight="bold", color="#C0392B",
                xytext=(6, 6), textcoords="offset points")

ax.set_xlabel("log fold change")
ax.set_ylabel("-log10(q-value)")
ax.set_title(f"ANCOM-BC2: AL vs PERM (n={len(df)} pathways, q<{Q_THRESHOLD} highlighted)")
ax.legend(loc="upper right")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"\n완료: {OUT_PNG}")
