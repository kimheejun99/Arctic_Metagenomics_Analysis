#!/usr/bin/env python3
"""
Step3_host_tracing.py

샘플별 폴더(./{sample}/) 안의 Kraken2 read-level 분류 결과(.kraken.out)와
VF(annotated_parsed.tsv) / ARG(blastout.filtered.txt) read-level 결과를
read ID로 조인해서, 카테고리/type별로 어떤 host phylum에서 검출됐는지 추적한다.

기대하는 파일 위치:
    ./{sample}/{sample}.kraken.out
    ./{sample}/{sample}_annotated_parsed.tsv
    ./{sample}/{sample}_args_output/blastout.filtered.txt

사용법:
    python3 Step3_host_tracing.py \
        --sample_list sample_list.txt \
        --metadata metadata.txt \
        --base_dir . \
        --taxdump /mnt/program/db/plusPFP_16 \
        --jobs 6 \
        --outdir ./host_tracing
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from multiprocessing import Pool

VALID_VF_CATEGORIES = {
    "Immune modulation", "Nutritional/Metabolic factor", "Adherence",
    "Effector delivery system", "Regulation", "Motility", "Biofilm",
    "Stress survival", "Exotoxin", "Others", "Exoenzyme",
    "Antimicrobial activity/Competitive advantage", "Invasion",
    "Post-translational modification"
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_list", default="sample_list.txt")
    p.add_argument("--metadata", default="metadata.txt")
    p.add_argument("--kraken_dir", default="./kraken")
    p.add_argument("--vfdb_dir", default="./vfdb")
    p.add_argument("--arg_dir", default="./arg")
    p.add_argument("--taxdump", default="/mnt/program/db/plusPFP_16")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--top_n_phyla", type=int, default=8)
    p.add_argument("--outdir", default="./host_tracing")
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


def load_nodes_names(taxdump):
    nodes_path = os.path.join(taxdump, "nodes.dmp")
    names_path = os.path.join(taxdump, "names.dmp")
    parent, rank = {}, {}
    with open(nodes_path) as f:
        for line in f:
            cols = [c.strip() for c in line.split("|")]
            parent[cols[0]] = cols[1]
            rank[cols[0]] = cols[2]
    name = {}
    with open(names_path) as f:
        for line in f:
            cols = [c.strip() for c in line.split("|")]
            if cols[3] == "scientific name":
                name[cols[0]] = cols[1]
    return parent, rank, name


def taxid_to_phylum(taxid, parent, rank, name, cache):
    if taxid in cache:
        return cache[taxid]
    current, seen = str(taxid), set()
    phylum = None
    while current in parent and current not in seen:
        seen.add(current)
        if rank.get(current) == "phylum":
            phylum = name.get(current, None)
            break
        if current == parent[current]:
            break
        current = parent[current]
    cache[taxid] = phylum
    return phylum


_WORKER = {}


def _init_worker(taxdump):
    parent, rank, name = load_nodes_names(taxdump)
    _WORKER["parent"], _WORKER["rank"], _WORKER["name"] = parent, rank, name
    _WORKER["cache"] = {}


def _extract_readid_arg(qseqid):
    return qseqid.split("@")[-1]


def _process_sample(args_tuple):
    sample, group, kraken_dir, vfdb_dir, arg_dir = args_tuple
    kraken_path = os.path.join(kraken_dir, f"{sample}.kraken.out")
    vf_path = os.path.join(vfdb_dir, f"{sample}_annotated_parsed.tsv")
    arg_path = os.path.join(arg_dir, f"{sample}_args_output", "blastout.filtered.txt")

    if not os.path.exists(kraken_path):
        return sample, group, None, None, "kraken.out 없음"

    # 필요한 read ID 수집 (VF + ARG)
    vf_df = pd.read_csv(vf_path, sep="\t") if os.path.exists(vf_path) else None
    arg_df = pd.read_csv(arg_path, sep="\t") if os.path.exists(arg_path) else None

    needed_ids = set()
    if vf_df is not None:
        needed_ids |= set(vf_df["SampleID"].astype(str))
    if arg_df is not None:
        needed_ids |= set(arg_df["qseqid"].astype(str).apply(_extract_readid_arg))

    readid_to_taxid = {}
    remaining = len(needed_ids)
    if remaining > 0:
        with open(kraken_path) as f:
            for line in f:
                if remaining == 0:
                    break
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3 or cols[0] != "C":
                    continue
                rid = cols[1]
                if rid in needed_ids:
                    readid_to_taxid[rid] = cols[2]
                    remaining -= 1

    parent, rank, name, cache = _WORKER["parent"], _WORKER["rank"], _WORKER["name"], _WORKER["cache"]

    vf_rows = []
    if vf_df is not None:
        for _, r in vf_df.iterrows():
            rid, cat = str(r["SampleID"]), r["Category"]
            taxid = readid_to_taxid.get(rid)
            phylum = taxid_to_phylum(taxid, parent, rank, name, cache) if taxid else None
            vf_rows.append({"Category": cat, "Phylum": phylum or "Unclassified"})

    arg_rows = []
    if arg_df is not None:
        for _, r in arg_df.iterrows():
            rid = _extract_readid_arg(str(r["qseqid"]))
            typ = r["type"]
            taxid = readid_to_taxid.get(rid)
            phylum = taxid_to_phylum(taxid, parent, rank, name, cache) if taxid else None
            arg_rows.append({"Type": typ, "Phylum": phylum or "Unclassified"})

    vf_out = pd.DataFrame(vf_rows)
    if len(vf_out):
        vf_out["group"] = group
        vf_out["sample"] = sample
    arg_out = pd.DataFrame(arg_rows)
    if len(arg_out):
        arg_out["group"] = group
        arg_out["sample"] = sample

    return sample, group, vf_out if len(vf_out) else None, arg_out if len(arg_out) else None, None


def make_stacked_plot(combined, group_col, cat_col, title_prefix, top_n_phyla, out_path):
    groups = sorted(combined[group_col].unique())
    top_phyla = (combined[combined["Phylum"] != "Unclassified"]["Phylum"]
                 .value_counts().head(top_n_phyla).index.tolist())

    def bucket(p):
        return "Unclassified" if p == "Unclassified" else (p if p in top_phyla else "Others")

    combined = combined.copy()
    combined["Phylum_bucket"] = combined["Phylum"].apply(bucket)
    phylum_order = top_phyla + ["Others", "Unclassified"]
    colors = cm.get_cmap("tab10", len(top_phyla))
    color_map = {p: colors(i) for i, p in enumerate(top_phyla)}
    color_map["Others"] = "#cccccc"
    color_map["Unclassified"] = "#dbe8c8"

    n_cats = {g: combined[combined[group_col] == g][cat_col].nunique() for g in groups}
    widths = [max(4, 0.4 * n_cats[g]) for g in groups]
    fig, axes = plt.subplots(1, len(groups), figsize=(sum(widths), 6.5), squeeze=False,
                              gridspec_kw={"width_ratios": widths})

    for i, g in enumerate(groups):
        ax = axes[0][i]
        sub = combined[combined[group_col] == g]
        ct = pd.crosstab(sub[cat_col], sub["Phylum_bucket"], normalize="index") * 100
        for p in phylum_order:
            if p not in ct.columns:
                ct[p] = 0
        ct = ct[phylum_order].sort_index()
        x = list(range(len(ct.index)))
        bottom = [0] * len(ct.index)
        for p in phylum_order:
            ax.bar(x, ct[p], bottom=bottom, color=color_map[p], label=p, edgecolor="white", linewidth=0.4, width=0.75)
            bottom = [b + v for b, v in zip(bottom, ct[p])]
        ax.set_title(f"Group {g} - {title_prefix}", fontsize=12, fontweight="bold")
        ax.set_ylabel("Host Contribution Ratio (%)")
        ax.set_ylim(0, 100)
        ax.set_xticks(x); ax.set_xticklabels(ct.index, rotation=60, ha="right", fontsize=8)
        ax.set_xlim(-0.6, len(x) - 0.4)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Host Phylum", bbox_to_anchor=(1.02, 0.9), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    samples = load_samples(args.sample_list)
    group_map = load_group_map(args.metadata, samples)
    tasks = [(s, group_map.get(s, "All"), args.kraken_dir, args.vfdb_dir, args.arg_dir) for s in samples]

    vf_all, arg_all = [], []
    with Pool(processes=args.jobs, initializer=_init_worker, initargs=(args.taxdump,)) as pool:
        for sample, group, vf_df, arg_df, err in pool.imap_unordered(_process_sample, tasks):
            if err:
                print(f"[Step3][경고] {sample}: {err} -> 스킵")
                continue
            if vf_df is not None:
                vf_all.append(vf_df)
            if arg_df is not None:
                arg_all.append(arg_df)
            print(f"[Step3] {sample} ({group}) 완료")

    if vf_all:
        vf_combined = pd.concat(vf_all, ignore_index=True)
        vf_combined.to_csv(os.path.join(args.outdir, "vf_host_tracing_raw.tsv"), sep="\t", index=False)
        vf_clean = vf_combined[vf_combined["Category"].isin(VALID_VF_CATEGORIES)]
        make_stacked_plot(vf_clean, "group", "Category", "VF Host Tracing", args.top_n_phyla,
                           os.path.join(args.outdir, "VF_host_tracing_by_group.png"))
        print("[Step3-VF] 저장 완료")

    if arg_all:
        arg_combined = pd.concat(arg_all, ignore_index=True)
        arg_combined.to_csv(os.path.join(args.outdir, "arg_host_tracing_raw.tsv"), sep="\t", index=False)
        make_stacked_plot(arg_combined, "group", "Type", "ARG Host Tracing", args.top_n_phyla,
                           os.path.join(args.outdir, "ARG_host_tracing_by_group.png"))
        print("[Step3-ARG] 저장 완료")

    print("[Step3] 완료.")


if __name__ == "__main__":
    main()
