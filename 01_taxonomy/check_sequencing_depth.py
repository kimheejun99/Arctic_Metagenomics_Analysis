import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# ==============================================================
# check_sequencing_depth.py
# .kreport에서 총 read 수(시퀀싱 깊이)를 추출하여
# AL과 PERM 사이에 체계적인 차이가 있는지 확인합니다.
# ==============================================================

parser = argparse.ArgumentParser()
parser.add_argument("--sample_list", default=None)
parser.add_argument("--outdir", default=".")
args, _ = parser.parse_known_args()

OUT_PNG = os.path.join(args.outdir, "sequencing_depth_AL_vs_PERM.png")
OUT_CSV = os.path.join(args.outdir, "sequencing_depth_per_sample.csv")

# kreport 파일 탐색 경로 (우선순위 순)
KREPORT_SEARCH_DIRS = ["./kraken", "../kraken", ".", ".."]

# metadata.txt 탐색 (현재 폴더 → 상위 폴더)
_META_PATH = None
for _mp in ["metadata.txt", "../metadata.txt"]:
    if os.path.exists(_mp):
        _META_PATH = _mp
        break
if _META_PATH is None:
    print("[오류] metadata.txt를 현재/상위 폴더에서 찾을 수 없습니다.")
    sys.exit(1)

_meta_df = pd.read_csv(_META_PATH, sep="\t")
_meta_df.columns = [c.strip() for c in _meta_df.columns]
AL_SAMPLES = _meta_df[_meta_df["group"] == "AL"]["sample-id"].astype(str).tolist()
PERM_SAMPLES = _meta_df[_meta_df["group"] == "PERM"]["sample-id"].astype(str).tolist()
print(f"[metadata.txt] AL {len(AL_SAMPLES)}개, PERM {len(PERM_SAMPLES)}개 로드됨")

def find_kreport(srr):
    """여러 폴더에서 .kreport 파일 탐색"""
    for d in KREPORT_SEARCH_DIRS:
        p = os.path.join(d, f"{srr}.kreport")
        if os.path.exists(p):
            return p
    return None


def total_reads_from_kreport(path):
    """kreport의 최상위 두 행(Unclassified + Root)의 clade_reads 합 = 전체 read 수"""
    total = 0
    found_levels = set()
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            clade_reads, rank_code = cols[1], cols[3]
            if rank_code in ("U", "R") and rank_code not in found_levels:
                total += int(clade_reads)
                found_levels.add(rank_code)
    return total


rows = []
missing = []
for group, samples in [("AL", AL_SAMPLES), ("PERM", PERM_SAMPLES)]:
    for srr in samples:
        path = find_kreport(srr)
        if path is None:
            missing.append(srr)
            continue
        total = total_reads_from_kreport(path)
        rows.append({"sample_id": srr, "group": group, "total_reads": total})

if missing:
    print(f"[안내] .kreport를 못 찾은 샘플 {len(missing)}개: {missing}")

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)

print(f"\n[정보] 총 {len(df)}개 샘플 처리됨")
if df.empty or "group" not in df.columns or len(df) < 2:
    print("[안내] 데이터가 부족하여 통계 분석을 건너뜁니다.")
    sys.exit(0)
print(df.groupby("group")["total_reads"].describe())

# ---- Mann-Whitney U test ----
al_vals = df[df["group"] == "AL"]["total_reads"]
perm_vals = df[df["group"] == "PERM"]["total_reads"]
stat, p = mannwhitneyu(al_vals, perm_vals, alternative="two-sided")
print(f"\n[Mann-Whitney U test] AL vs PERM 총 read 수 차이: p = {p:.4f}")
if p < 0.05:
    print("  -> 유의한 차이 있음. 시퀀싱 깊이 차이가 다양성 비교에 혼란 요인일 수 있음.")
else:
    print("  -> 유의한 차이 없음. 시퀀싱 깊이는 AL/PERM 다양성 차이의 주 원인이 아닐 가능성이 높음.")

# ---- 시각화 ----
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot(data=df, x="group", y="total_reads", palette={"AL": "#4C72B0", "PERM": "#DD8452"}, ax=ax)
sns.stripplot(data=df, x="group", y="total_reads", color="black", alpha=0.5, size=4, ax=ax)
ax.set_ylabel("Total reads per sample")
ax.set_title(f"Sequencing depth: AL vs PERM (Mann-Whitney p={p:.4f})")
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"\n완료: {OUT_PNG}, {OUT_CSV}")
