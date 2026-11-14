import pandas as pd
import numpy as np
from scipy import stats

# ==============================================================
# all_pathways_test.py
# merged_pathways.tsv의 "전체" KEGG Module에 대해 AL vs PERM
# Mann-Whitney U test를 수행하고, FDR(Benjamini-Hochberg) 보정을
# 적용합니다. 특정 pathway를 미리 골라서 보는 게 아니라, 있는 대로
# 전부 검정해서 "우연히 골라서 유의했던 것"이라는 의심을 없애고,
# 놓쳤을 수 있는 다른 유의한 pathway도 찾아냅니다.
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python all_pathways_test.py
# 결과: kegganog_results/all_pathways_stats.csv (q-value 오름차순 정렬)
# ==============================================================

IN_TSV = "kegganog_results/merged_pathways.tsv"
OUT_CSV = "kegganog_results/all_pathways_stats.csv"

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


def benjamini_hochberg(pvals):
    """scipy 버전이 낮아 false_discovery_control이 없을 경우를 대비한
    수동 BH 보정 구현 (표준 공식 그대로)."""
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.clip(ranked_sorted, 0, 1)
    return qvals


df = pd.read_csv(IN_TSV, sep="\t")
pathway_col = df.columns[0]


def find_col(srr_id, columns):
    matches = [c for c in columns if srr_id in c]
    if len(matches) != 1:
        raise ValueError(f"{srr_id}: {len(matches)}개 매칭")
    return matches[0]


al_cols = [find_col(s, df.columns) for s in AL_SAMPLES]
perm_cols = [find_col(s, df.columns) for s in PERM_SAMPLES]

print(f"전체 pathway 수: {len(df)}, AL={len(al_cols)}개, PERM={len(perm_cols)}개")

results = []
for _, row in df.iterrows():
    pw = row[pathway_col]
    al_vals = row[al_cols].values.astype(float)
    perm_vals = row[perm_cols].values.astype(float)

    al_mean, perm_mean = al_vals.mean(), perm_vals.mean()

    if np.all(al_vals == al_vals[0]) and np.all(perm_vals == perm_vals[0]) and al_vals[0] == perm_vals[0]:
        p = 1.0
    else:
        try:
            _, p = stats.mannwhitneyu(al_vals, perm_vals, alternative="two-sided")
        except ValueError:
            p = 1.0

    results.append({
        "Pathway": pw,
        "AL_mean": al_mean,
        "PERM_mean": perm_mean,
        "Direction": "PERM > AL" if perm_mean > al_mean else ("AL > PERM" if al_mean > perm_mean else "동일"),
        "p_value": p,
    })

res_df = pd.DataFrame(results)

try:
    res_df["q_value"] = stats.false_discovery_control(res_df["p_value"].values, method="bh")
    print("scipy.stats.false_discovery_control 사용")
except AttributeError:
    res_df["q_value"] = benjamini_hochberg(res_df["p_value"].values)
    print("scipy 버전이 낮아 수동 BH 보정 사용")

res_df["Significant (q<0.05)"] = res_df["q_value"] < 0.05

res_df = res_df.sort_values("q_value")
res_df.to_csv(OUT_CSV, index=False)

n_sig = res_df["Significant (q<0.05)"].sum()
print(f"\n=== 완료: {OUT_CSV} ===")
print(f"전체 {len(res_df)}개 pathway 중 FDR<0.05 기준 유의: {n_sig}개\n")

print("[상위 20개, q-value 오름차순]")
print(res_df.head(20)[["Pathway", "Direction", "p_value", "q_value"]].to_string(index=False))

print("\n[메탄생성 관련 pathway들의 보정 후 결과 확인]")
methano_check = res_df[res_df["Pathway"].str.contains("Methanogenesis|Methanol dehydrogenase|Coenzyme M", case=False, na=False)]
print(methano_check[["Pathway", "Direction", "p_value", "q_value", "Significant (q<0.05)"]].to_string(index=False))
