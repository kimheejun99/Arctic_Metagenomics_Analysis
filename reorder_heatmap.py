import kegganog as kgn
import pandas as pd

# ==============================================================
# reorder_heatmap.py
# merged_pathways.tsv를 AL(36개) -> PERM(13개) 순서로 재배치해서
# 히트맵을 다시 그립니다.
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python reorder_heatmap.py
# ==============================================================

IN_TSV = "kegganog_results/merged_pathways.tsv"
OUT_PREFIX = "kegganog_results/heatmap_AL_PERM_sorted"

# sample_group_mapping.csv 기준 실제 그룹 목록
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

# ---- 1) 데이터 읽기 ----
df = pd.read_csv(IN_TSV, sep="\t")
print(f"원본 컬럼 수: {len(df.columns)}, 컬럼 예시: {list(df.columns[:3])}")

# ---- 2) 실제 컬럼명 자동 매칭 (SRR ID가 컬럼명에 포함돼 있는지로 찾음) ----
def find_col(srr_id, columns):
    matches = [c for c in columns if srr_id in c]
    if len(matches) != 1:
        raise ValueError(f"{srr_id}에 대해 컬럼을 {len(matches)}개 찾음(1개여야 함): {matches}")
    return matches[0]

pathway_col = df.columns[0]  # 첫 컬럼(pathway/module 이름)
al_cols = [find_col(s, df.columns) for s in AL_SAMPLES]
perm_cols = [find_col(s, df.columns) for s in PERM_SAMPLES]

print(f"AL 컬럼 {len(al_cols)}개, PERM 컬럼 {len(perm_cols)}개 매칭됨 (49개여야 함)")
assert len(al_cols) + len(perm_cols) == len(AL_SAMPLES) + len(PERM_SAMPLES), "49개가 안 맞습니다. 위 매칭 결과를 확인하세요."

# ---- 3) 재배치 ----
ordered_cols = [pathway_col] + al_cols + perm_cols
df_ordered = df[ordered_cols]

# ---- 4) 그룹 없이(-g 옵션의 레이아웃 버그 회피) 재시도, 넉넉한 figsize ----
kgnheat = kgn.heatmap(df_ordered, color="Blues", group=False, figsize=(26, 30))

# ---- 5) AL/PERM 경계에 구분선 + 라벨 추가 (어디까지가 AL/PERM인지 명확하게) ----
import matplotlib.pyplot as plt
import numpy as np

n_al = len(al_cols)
n_perm = len(perm_cols)

# kgnheat.ax가 단일 Axes가 아니라 여러 개로 나뉜 배열(Part1/2/3)일 수 있어
# 배열이면 각 축에 전부 적용, 단일 축이면 그대로 적용
axes = kgnheat.ax if isinstance(kgnheat.ax, (list, np.ndarray)) else [kgnheat.ax]

for ax in axes:
    ax.axvline(x=n_al, color="red", linewidth=2.5)
    ax.text(n_al / 2, -1.5, "AL (n=36)", ha="center", va="bottom",
            fontsize=14, fontweight="bold", color="black")
    ax.text(n_al + n_perm / 2, -1.5, "PERM (n=13)", ha="center", va="bottom",
            fontsize=14, fontweight="bold", color="black")

kgnheat.savefig(f"{OUT_PREFIX}.png", dpi=300)

print(f"완료: {OUT_PREFIX}.png (왼쪽=AL 36개, 오른쪽=PERM 13개)")
