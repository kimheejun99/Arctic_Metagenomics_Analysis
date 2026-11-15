#!/usr/bin/env python3
"""
VF_copies_per_cell.py

VF(virulence factor) 카테고리별 존재량을 "세포당 개수"로 정규화합니다.
ARGs-OAP가 args_output/metadata.txt 안에 계산해두는 nCell(단일카피
마커 유전자 기반 세포수 추정치)을 그대로 재사용하여, VF 카운트를
nCell로 나눕니다 (기존 VF2.py, VF_bar.py의 로직을 그대로 계승,
하드코딩된 49개 site A~G 대신 metadata.txt의 AL/PERM을 동적으로 읽음).

출력:
    VF_copies_per_cell_by_group_heatmap.png  (그룹별 Top10 히트맵)
    VF_copies_per_cell_overall_barplot.png   (전체 평균 Top10 막대그래프)
    VF_copies_per_cell_by_group.csv

사용법:
    python3 VF_copies_per_cell.py --sample_list sample_list.txt \
        --metadata metadata.txt --vfdb_dir ./vfdb --arg_dir ./arg \
        --outdir ./plots
"""

import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_list", default="sample_list.txt")
    p.add_argument("--metadata", default="metadata.txt", help="sample-id\\tgroup 형식")
    p.add_argument("--vfdb_dir", default="./vfdb", help="{sample}_annotated_parsed.tsv 들이 있는 폴더")
    p.add_argument("--arg_dir", default="./arg", help="{sample}_args_output/ 들이 있는 폴더")
    p.add_argument("--top_n", type=int, default=10)
    p.add_argument("--outdir", default="./plots")
    return p.parse_args()


def load_samples(sample_list_path):
    samples = []
    with open(sample_list_path) as f:
        for line in f:
            s = line.strip()
            if s and s.lower() != "run":
                samples.append(s)
    return samples


def load_group_map(metadata_path, samples):
    if not os.path.exists(metadata_path):
        return {s: "All" for s in samples}
    df = pd.read_csv(metadata_path, sep="\t")
    return dict(zip(df["sample-id"], df["group"]))


def load_ncell(arg_dir, sample):
    """ARGs-OAP가 args_output/ 안에 만드는 metadata.txt에서 nCell 값을 읽음."""
    meta_path = os.path.join(arg_dir, f"{sample}_args_output", "metadata.txt")
    if not os.path.exists(meta_path):
        return None
    try:
        meta_df = pd.read_csv(meta_path, sep=r"\s+", index_col="sample")
        if sample in meta_df.index:
            return meta_df.loc[sample, "nCell"]
    except Exception as e:
        print(f"[경고] {sample}: nCell 읽기 실패 - {e}")
    return None


def load_vf_copies_per_cell(vfdb_dir, arg_dir, sample):
    vf_path = os.path.join(vfdb_dir, f"{sample}_annotated_parsed.tsv")
    if not os.path.exists(vf_path):
        return None

    ncell = load_ncell(arg_dir, sample)
    if ncell is None or ncell == 0:
        print(f"[경고] {sample}: nCell 값을 찾을 수 없음 -> 스킵")
        return None

    vf_df = pd.read_csv(vf_path, sep="\t")
    if "Category" not in vf_df.columns:
        print(f"[경고] {sample}: Category 컬럼 없음 -> 스킵")
        return None

    raw_counts = vf_df["Category"].value_counts()
    return raw_counts / ncell


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    samples = load_samples(args.sample_list)
    group_map = load_group_map(args.metadata, samples)
    print(f"[VF_copies_per_cell] 대상 샘플 {len(samples)}개")

    per_sample = {}
    for s in samples:
        result = load_vf_copies_per_cell(args.vfdb_dir, args.arg_dir, s)
        if result is not None:
            per_sample[s] = result

    if not per_sample:
        print("[오류] 계산할 수 있는 샘플이 하나도 없습니다 (nCell 또는 VF 파일 확인 필요)")
        return

    combined_df = pd.concat(per_sample, axis=1).fillna(0.0)

    # ---- ① 그룹별(AL/PERM) 평균 히트맵 ----
    combined_t = combined_df.T
    combined_t["Group"] = combined_t.index.map(lambda s: group_map.get(s, "All"))
    group_avg = combined_t.groupby("Group").mean(numeric_only=True).T

    top_n_categories = group_avg.sum(axis=1).nlargest(args.top_n).index
    group_avg_top = group_avg.loc[top_n_categories]
    group_avg_top.to_csv(os.path.join(args.outdir, "VF_copies_per_cell_by_group.csv"))

    plt.figure(figsize=(10, 8))
    sns.heatmap(group_avg_top, annot=True, cmap="YlGnBu", fmt=".4f",
                linewidths=.5, cbar_kws={"label": "Average VF copies per cell"})
    plt.title(f"Top {args.top_n} VF Category by Group (copies per cell, absolute)", fontsize=14, pad=20)
    plt.xlabel("Groups")
    plt.ylabel("VF Categories")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "VF_copies_per_cell_by_group_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("[VF_copies_per_cell] 저장: VF_copies_per_cell_by_group_heatmap.png")

    # ---- ② 전체 평균 막대그래프 ----
    overall_mean = combined_df.mean(axis=1)
    top_overall = overall_mean.nlargest(args.top_n).sort_values(ascending=True)

    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.barh(top_overall.index, top_overall.values, color="#405B94", edgecolor="black", height=0.7)
    for bar in bars:
        width = bar.get_width()
        x_pos = width + (top_overall.max() * 0.015)
        y_pos = bar.get_y() + bar.get_height() / 2
        ax.text(x_pos, y_pos, f"{width:.4f}", va="center", ha="left", fontsize=11, color="#333333")

    ax.set_xlabel("Average VF copies per cell", fontsize=13, fontweight="bold", labelpad=15)
    ax.set_ylabel("VF Categories", fontsize=13, fontweight="bold", labelpad=15)
    plt.yticks(fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xlim(0, top_overall.max() * 1.15)
    plt.title(f"Top {args.top_n} Overall VF Category Composition ({len(per_sample)} Samples)",
              fontsize=15, fontweight="bold", pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "VF_copies_per_cell_overall_barplot.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("[VF_copies_per_cell] 저장: VF_copies_per_cell_overall_barplot.png")

    print("[VF_copies_per_cell] 완료.")


if __name__ == "__main__":
    main()
