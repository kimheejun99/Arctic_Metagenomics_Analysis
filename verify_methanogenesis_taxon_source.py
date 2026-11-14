import pandas as pd

# ==============================================================
# verify_methanogenesis_taxon_source.py
# KEGG-Decoder(KEGGaNOG 내부 로직) 소스코드에서 직접 확인한, 메탄생성
# 8개 항목의 실제 구성 KO 목록을 이용해, host_tracing_master.tsv에서
# 각 KO가 실제로 어느 taxon(taxid)의 read에서 검출됐는지 확인합니다.
#
# 목적: "completeness가 높게 나온 게 한 종(예: Methanosarcina)에서
# 나온 것인지, 여러 종에 흩어진 것인지"를 직접 확인.
#
# 출처: https://github.com/bjtully/BioData/blob/master/KEGGDecoder/KEGG_decoder.py
#       의 methanogenesis() 함수 (2026-09 확인)
#
# 사용법: /root/miniforge3/envs/eggnog/bin/python verify_methanogenesis_taxon_source.py
# ==============================================================

IN_TSV = "host_tracing_master.tsv"
TAXID_NAME_FILE = "../merged_tables_excl_outlier_and_contam/2_taxonomy_table.txt"
NAMES_DMP = "/mnt/program/db/plusPFP_16/names.dmp"

# KEGG-Decoder 실제 소스코드에서 확인된 구성 KO (가중치는 표시만, 집계엔 안 씀)
PATHWAY_KOS = {
    "Methanogenesis via CO2": [
        "K00200", "K00201", "K00202", "K00203", "K00205", "K11261",
        "K00672", "K01499", "K13942", "K00320", "K00577", "K00578",
        "K00579", "K00580", "K00581", "K00582", "K00583", "K00584",
    ],
    "Methanogenesis via acetate": ["K00193", "K00194", "K00197"],
    "Methanogenesis via methanol": ["K14080", "K04480", "K14081"],
    "Methanogenesis via dimethylamine": ["K14082", "K16178"],
    "Methanogenesis via methylamine": ["K16178"],
    "Methanogenesis via trimethylamine": ["K14083"],
    "Coenzyme M reduction to methane": ["K00399", "K00401", "K00402"],
    "Coenzyme B/Coenzyme M regeneration": ["K03388", "K03389", "K03390", "K08264", "K08265"],
}

# ---- 1) taxid 이름 매핑 (2단계: taxonomy_table.txt -> names.dmp) ----
taxid_to_name = {}
try:
    with open(TAXID_NAME_FILE) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            taxid = parts[0]
            name = None
            for val in reversed(parts[1:]):
                val = val.strip()
                stripped = val.split("__", 1)[-1] if "__" in val else val
                if stripped and stripped not in ("-", "NA"):
                    name = stripped
                    break
            if name:
                taxid_to_name[taxid] = name
except FileNotFoundError:
    print(f"[안내] {TAXID_NAME_FILE} 없음, 1차 매핑 건너뜀")


def get_name(taxid, extra_cache):
    if taxid in taxid_to_name:
        return taxid_to_name[taxid]
    if taxid in extra_cache:
        return extra_cache[taxid]
    return f"taxid_{taxid}"


# ---- 2) host_tracing_master.tsv 읽기 ----
df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
df = df[df["taxid"] != "0"]

# ---- 3) 아직 이름 없는 taxid들을 한 번에 names.dmp에서 보강 ----
all_taxids_in_data = set(df["taxid"].unique())
missing = [t for t in all_taxids_in_data if t not in taxid_to_name]
extra_cache = {}
if missing:
    missing_set = set(missing)
    try:
        with open(NAMES_DMP) as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, name_class = parts[0], parts[1], parts[2], parts[3]
                if taxid in missing_set and name_class == "scientific name":
                    extra_cache[taxid] = name
    except FileNotFoundError:
        print(f"[안내] {NAMES_DMP} 없음, 2차 매핑 건너뜀")

# ---- 4) 각 pathway마다, 구성 KO별로 taxid 분포 확인 ----
for pathway, kos in PATHWAY_KOS.items():
    print(f"\n{'='*70}")
    print(f"[{pathway}]  (구성 KO: {', '.join(kos)})")
    print(f"{'='*70}")

    sub = df[df["KO"].isin(kos)]
    if len(sub) == 0:
        print("  -> 이 pathway를 구성하는 KO를 가진 read가 host_tracing_master.tsv에 없음")
        continue

    # KO별 + taxid별로 몇 개 read인지
    pivot = sub.groupby(["KO", "taxid"]).size().reset_index(name="read_count")
    pivot["taxon_name"] = pivot["taxid"].apply(lambda t: get_name(t, extra_cache))

    for ko in kos:
        ko_rows = pivot[pivot["KO"] == ko].sort_values("read_count", ascending=False)
        if len(ko_rows) == 0:
            print(f"  {ko}: 검출된 read 없음")
            continue
        top = ko_rows.head(3)
        detail = ", ".join(f"{r.taxon_name}({r.read_count})" for r in top.itertuples())
        n_distinct_taxa = ko_rows["taxid"].nunique()
        print(f"  {ko}: {n_distinct_taxa}개 taxon에서 검출, 상위 = {detail}")

    # 이 pathway 전체에서, taxon별로 "몇 개의 서로 다른 구성 KO를 커버하는지"
    coverage = sub.groupby("taxid")["KO"].nunique().sort_values(ascending=False)
    print(f"\n  [taxon별 커버리지: {len(kos)}개 구성 KO 중 몇 개를 갖고 있나]")
    for taxid, n_ko in coverage.head(5).items():
        name = get_name(taxid, extra_cache)
        print(f"    {name} ({taxid}): {n_ko}/{len(kos)} KO 보유")

print(f"\n{'='*70}")
print("해석: 위에서 '어느 한 taxon이 그 pathway 구성 KO 대부분(또는 전부)을")
print("혼자 커버'하면 -> completeness가 실제로 한 종의 능력을 반영한다는 뜻")
print("여러 taxon이 각각 한두 개씩만 나눠 가지면 -> 커뮤니티 차원의 '조각모음'")
print("일 뿐, 한 종이 이 경로를 완전히 갖췄다고 볼 근거는 약해집니다.")
