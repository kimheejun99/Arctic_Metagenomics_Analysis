#!/usr/bin/env python3
"""
Step4_merge_kreports.py

샘플별 폴더(./{sample}/{sample}.kreport)를 병합하여
1_ASV_table.txt / 2_taxonomy_table.txt / 3_metadata.txt 3개 파일을 생성한다.
(Eukaryota 계통 필터링 + 동일 종명 taxid 중복 병합 포함)

사용법:
    python3 Step4_merge_kreports.py \
        --sample_list sample_list.txt \
        --metadata metadata.txt \
        --base_dir . \
        --taxdump /mnt/program/db/plusPFP_16 \
        --outdir ./merged_tables
"""

import argparse
import os
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_list", default="sample_list.txt")
    p.add_argument("--metadata", default="metadata.txt", help="sample-id\\tgroup (없으면 group 컬럼 비움)")
    p.add_argument("--kraken_dir", default="./kraken")
    p.add_argument("--taxdump", default="/mnt/program/db/plusPFP_16")
    p.add_argument("--level", default="S", choices=["S", "G", "F", "O", "C", "P", "D"])
    p.add_argument("--outdir", default="./merged_tables")
    return p.parse_args()


def load_samples(sample_list_path):
    samples = []
    with open(sample_list_path) as f:
        for line in f:
            s = line.strip()
            if s and s.lower() != "run":
                samples.append(s)
    return samples


def load_group_map(metadata_path):
    if not os.path.exists(metadata_path):
        return {}
    df = pd.read_csv(metadata_path, sep="\t")
    return dict(zip(df["sample-id"], df["group"]))


# ---------------- 1. kreport -> ASV table ----------------
def parse_kreport(path, level):
    result = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            clade_reads, rank_code, taxid = cols[1], cols[3], cols[4]
            if rank_code == level:
                result[taxid] = result.get(taxid, 0) + int(clade_reads)
    return result


def build_asv_table(kraken_dir, samples, level):
    combined = {}
    for s in samples:
        f = os.path.join(kraken_dir, f"{s}.kreport")
        if not os.path.exists(f):
            print(f"[Step4][경고] {f} 없음 -> 스킵")
            continue
        combined[s] = parse_kreport(f, level)
    asv = pd.DataFrame(combined).fillna(0).astype(int)
    asv.index.name = "#NAME"
    return asv


# ---------------- 2. taxonomy table ----------------
def load_nodes_names(taxdump):
    parent, rank = {}, {}
    with open(os.path.join(taxdump, "nodes.dmp")) as f:
        for line in f:
            cols = [c.strip() for c in line.split("|")]
            parent[cols[0]] = cols[1]
            rank[cols[0]] = cols[2]
    name = {}
    with open(os.path.join(taxdump, "names.dmp")) as f:
        for line in f:
            cols = [c.strip() for c in line.split("|")]
            if cols[3] == "scientific name":
                name[cols[0]] = cols[1]
    return parent, rank, name


LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
LEVEL_TO_COL = {"superkingdom": "Kingdom", "phylum": "Phylum", "class": "Class", "order": "Order",
                 "family": "Family", "genus": "Genus", "species": "Species"}


def lineage_for_taxid(taxid, parent, rank, name, cache):
    if taxid in cache:
        return cache[taxid]
    lineage = {lvl: "" for lvl in LEVELS}
    current, seen = str(taxid), set()
    while current in parent and current not in seen:
        seen.add(current)
        r = rank.get(current, "")
        if r in lineage:
            lineage[r] = name.get(current, "")
        if current == parent[current]:
            break
        current = parent[current]
    result = {LEVEL_TO_COL[lvl]: lineage[lvl] for lvl in LEVELS}
    cache[taxid] = result
    return result


def build_taxonomy_table(taxids, taxdump):
    cols = ["#TAXONOMY", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if not taxids:
        return pd.DataFrame(columns=cols)
    parent, rank, name = load_nodes_names(taxdump)
    cache = {}
    rows = []
    for t in taxids:
        lin = lineage_for_taxid(str(t), parent, rank, name, cache)
        rows.append({"#TAXONOMY": t, **lin})
    return pd.DataFrame(rows)[cols]


EUKARYOTA_KINGDOMS = {"Eukaryota", "Fungi", "Viridiplantae", "Metazoa", "Chordata",
                       "Sar", "Discoba", "Amoebozoa", "Rhodophyta", "Cryptophyceae", "Haptista"}


def filter_eukaryota(asv, tax):
    keep = tax[~tax["Kingdom"].isin(EUKARYOTA_KINGDOMS)]["#TAXONOMY"].astype(str).tolist()
    before = asv.shape[0]
    asv_f = asv.loc[asv.index.astype(str).isin(keep)]
    tax_f = tax[tax["#TAXONOMY"].astype(str).isin(keep)].reset_index(drop=True)
    print(f"[Step4] Eukaryota 필터링: {before} -> {asv_f.shape[0]} taxa")
    return asv_f, tax_f


def dedup_by_species(asv, tax):
    tax_idx = tax.set_index(tax["#TAXONOMY"].astype(str))
    species_of = tax_idx["Species"]
    named = species_of[species_of != ""]
    groups = named.groupby(named).groups
    dup_groups = {sp: list(ids) for sp, ids in groups.items() if len(ids) > 1}
    if not dup_groups:
        print("[Step4] 동일 종명 중복 없음")
        return asv, tax_idx.reset_index(drop=True)

    new_asv, new_tax, used = {}, [], set()
    for sp, ids in dup_groups.items():
        present = [i for i in ids if i in asv.index.astype(str).tolist()]
        if not present:
            continue
        merged = asv.loc[asv.index.astype(str).isin(present)].sum(axis=0)
        rep = present[0]
        new_asv[rep] = merged
        new_tax.append(tax_idx.loc[rep].copy())
        used.update(present)
        print(f"[Step4][Dedup] '{sp}' -> {present} 를 {rep} 로 병합")

    remaining = [i for i in asv.index.astype(str) if i not in used]
    asv_out = pd.concat([pd.DataFrame(new_asv).T if new_asv else pd.DataFrame(),
                          asv.loc[asv.index.astype(str).isin(remaining)]])
    tax_out = pd.concat([pd.DataFrame(new_tax) if new_tax else pd.DataFrame(),
                          tax_idx.loc[tax_idx.index.isin(remaining)]]).reset_index(drop=True)
    return asv_out, tax_out


# ---------------- 3. metadata ----------------
def get_classified_rate(kreport_path):
    try:
        with open(kreport_path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 6 and cols[3] == "U":
                    return round(100 - float(cols[0]), 4)
        return 100.0
    except (FileNotFoundError, ValueError):
        return None


def build_metadata_table(kraken_dir, samples, group_map):
    rows = []
    for s in samples:
        rate = get_classified_rate(os.path.join(kraken_dir, f"{s}.kreport"))
        rows.append({"#SampleID": s, "group": group_map.get(s, ""), "Classified_Rate": rate})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    samples = load_samples(args.sample_list)
    group_map = load_group_map(args.metadata)
    print(f"[Step4] 샘플 {len(samples)}개 대상 병합 시작")

    asv = build_asv_table(args.kraken_dir, samples, args.level)
    print(f"[Step4] 1차 ASV 매트릭스: {asv.shape[0]} taxa x {asv.shape[1]} 샘플")

    if asv.shape[1] == 0:
        print(f"[Step4][오류] kreport를 하나도 찾지 못했습니다.")
        print(f"  확인해주세요: {args.kraken_dir}/{{샘플ID}}.kreport 경로에 파일이 있는지")
        raise SystemExit(1)

    tax = build_taxonomy_table(asv.index.tolist(), args.taxdump)
    asv, tax = filter_eukaryota(asv, tax)
    asv, tax = dedup_by_species(asv, tax)
    tax = tax.rename(columns={tax.columns[0]: "#TAXONOMY"}) if tax.columns[0] != "#TAXONOMY" else tax

    asv.to_csv(os.path.join(args.outdir, "1_ASV_table.txt"), sep="\t")
    tax.to_csv(os.path.join(args.outdir, "2_taxonomy_table.txt"), sep="\t", index=False)
    meta = build_metadata_table(args.kraken_dir, samples, group_map)
    meta.to_csv(os.path.join(args.outdir, "3_metadata.txt"), sep="\t", index=False)

    print(f"[Step4] 완료: taxa {asv.shape[0]}개, 샘플 {asv.shape[1]}개")


if __name__ == "__main__":
    main()
