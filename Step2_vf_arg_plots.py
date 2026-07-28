#!/usr/bin/env python3
"""
Step2_vf_arg_plots.py

샘플별 폴더(./{sample}/) 안의 VF(diamond) 결과와 ARGs-OAP 결과를 읽어
1) 전체 샘플 합산/평균 바그래프
2) 그룹별 히트맵
을 VF, ARG 각각 생성한다.

기대하는 파일 위치:
    ./{sample}/{sample}_annotated_parsed.tsv        (VF, Category 컬럼 사용)
    ./{sample}/{sample}_args_output/normalized_cell.type.txt   (ARG)

사용법:
    python3 Step2_vf_arg_plots.py \
        --sample_list sample_list.txt \
        --metadata metadata.txt \
        --base_dir . \
        --outdir ./plots
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_list", default="sample_list.txt")
    p.add_argument("--metadata", default="metadata.txt", help="sample-id\\tgroup 형식")
    p.add_argument("--vfdb_dir", default="./vfdb", help="{sample}_annotated_parsed.tsv 들이 있는 폴더")
    p.add_argument("--arg_dir", default="./arg", help="{sample}_args_output/ 들이 있는 폴더")
    p.add_argument("--arg_metric", default="normalized_cell",
                    choices=["normalized_cell", "normalized_16S", "tpm", "ppm", "rpkm",
                             "unnormalized_copy", "unnormalized_count"])
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


# ---------------- VF ----------------
def load_vf_category_pct(vfdb_dir, sample):
    f = os.path.join(vfdb_dir, f"{sample}_annotated_parsed.tsv")
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f, sep="\t")
    if "Category" not in df.columns:
        return None
    counts = df["Category"].value_counts()
    return counts, counts / counts.sum() * 100


def run_vf(args, samples, group_map):
    per_sample_counts, per_sample_pct = {}, {}
    for s in samples:
        result = load_vf_category_pct(args.vfdb_dir, s)
        if result is None:
            print(f"[Step2-VF][경고] {s}: annotated_parsed.tsv 없음 -> 스킵")
            continue
        per_sample_counts[s], per_sample_pct[s] = result

    total_counts = pd.DataFrame(per_sample_counts).fillna(0).sum(axis=1)
    total_pct = (total_counts / total_counts.sum() * 100).sort_values(ascending=False)
    top_categories = total_pct.head(args.top_n)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(top_categories.index[::-1], top_categories.values[::-1], color="#1f77b4")
    ax.set_title("Overall Total: VF Category Ratio", fontsize=13)
    ax.set_xlabel("Percentage of total reads (%)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "VF_total_category_barplot.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("[Step2-VF] 저장: VF_total_category_barplot.png")

    pct_df = pd.DataFrame(per_sample_pct).fillna(0)
    group_avg = {}
    for group in sorted(set(group_map.get(s, "All") for s in per_sample_pct)):
        members = [s for s in per_sample_pct if group_map.get(s, "All") == group]
        if members:
            group_avg[group] = pct_df[members].mean(axis=1)
    heat_df = pd.DataFrame(group_avg).loc[top_categories.index]

    fig, ax = plt.subplots(figsize=(7, 8))
    im = ax.imshow(heat_df.values, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(heat_df.columns))); ax.set_xticklabels(heat_df.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(heat_df.index))); ax.set_yticklabels(heat_df.index)
    for i in range(heat_df.shape[0]):
        for j in range(heat_df.shape[1]):
            ax.text(j, i, f"{heat_df.values[i,j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Mean relative abundance (%)")
    ax.set_title(f"Top {args.top_n} VF Category by Group")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "VF_group_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()
    heat_df.to_csv(os.path.join(args.outdir, "VF_group_relative_abundance_table.tsv"), sep="\t")
    print("[Step2-VF] 저장: VF_group_heatmap.png")


# ---------------- ARG ----------------
def load_arg_type_values(arg_dir, sample, metric):
    f = os.path.join(arg_dir, f"{sample}_args_output", f"{metric}.type.txt")
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f, sep="\t")
    if df.shape[1] < 2:
        return None
    s = df.set_index(df.columns[0])[df.columns[-1]]
    return pd.to_numeric(s, errors="coerce").fillna(0)


def run_arg(args, samples, group_map):
    per_sample_values = {}
    for s in samples:
        vals = load_arg_type_values(args.arg_dir, s, args.arg_metric)
        if vals is None:
            print(f"[Step2-ARG][경고] {s}: {args.arg_metric}.type.txt 없음 -> 스킵")
            continue
        per_sample_values[s] = vals

    value_df = pd.DataFrame(per_sample_values).fillna(0)
    overall_mean = value_df.mean(axis=1).sort_values(ascending=False)
    top_types = overall_mean.head(args.top_n)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(top_types.index[::-1], top_types.values[::-1], color="#ff7f0e")
    ax.set_title(f"Overall Mean: ARG Type Abundance ({args.arg_metric})", fontsize=13)
    ax.set_xlabel(f"Mean abundance ({args.arg_metric})")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f"ARG_total_type_barplot_{args.arg_metric}.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("[Step2-ARG] 저장: ARG_total_type_barplot.png")

    group_avg = {}
    for group in sorted(set(group_map.get(s, "All") for s in value_df.columns)):
        members = [s for s in value_df.columns if group_map.get(s, "All") == group]
        if members:
            group_avg[group] = value_df[members].mean(axis=1)
    heat_df = pd.DataFrame(group_avg).loc[top_types.index]

    fig, ax = plt.subplots(figsize=(7, 8))
    im = ax.imshow(heat_df.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(heat_df.columns))); ax.set_xticklabels(heat_df.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(heat_df.index))); ax.set_yticklabels(heat_df.index)
    for i in range(heat_df.shape[0]):
        for j in range(heat_df.shape[1]):
            ax.text(j, i, f"{heat_df.values[i,j]:.3f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label=f"abundance ({args.arg_metric})")
    ax.set_title(f"Top {args.top_n} ARG Type by Group")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f"ARG_group_heatmap_{args.arg_metric}.png"), dpi=300, bbox_inches="tight")
    plt.close()
    heat_df.to_csv(os.path.join(args.outdir, f"ARG_group_{args.arg_metric}_table.tsv"), sep="\t")
    print("[Step2-ARG] 저장: ARG_group_heatmap.png")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    samples = load_samples(args.sample_list)
    group_map = load_group_map(args.metadata, samples)
    print(f"[Step2] 대상 샘플 {len(samples)}개")

    run_vf(args, samples, group_map)
    run_arg(args, samples, group_map)
    print("[Step2] 완료.")


if __name__ == "__main__":
    main()
