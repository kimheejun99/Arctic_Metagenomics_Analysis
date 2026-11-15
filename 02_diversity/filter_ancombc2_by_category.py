import pandas as pd
import requests
import re
import json
import os

# ==============================================================
# filter_ancombc2_by_category.py
# KEGG 공식 BRITE 계층(br08901)을 실시간으로 받아와서, ANCOM-BC2
# 결과에서 "Human Diseases", "Organismal Systems", "Drug Development"
# 카테고리(세균/고균 토양 메타지놈과 무관한, 사람/동물 특이적 pathway들)를
# 제외한 뒤 재검토합니다.
#
# 왜 이 3개만 빼는가: 이 3개 대분류는 KEGG 공식 정의상 전부 사람/동물
# 생리·질병·약물에 관한 것이라, 세균/고균에는 해당사항이 없습니다.
# 반대로 Metabolism/Genetic Information Processing/Environmental
# Information Processing/Cellular Processes는 세균 관련 pathway가
# 다수 섞여있어 통째로 빼면 안 됩니다(개별 검토 필요).
#
# 사용법:
#   /root/miniforge3/envs/eggnog/bin/python filter_ancombc2_by_category.py
# ==============================================================

IN_CSV = "ancombc2_pathway_results.csv"
OUT_CSV = "ancombc2_pathway_results_filtered.csv"
OUT_EXCLUDED_CSV = "ancombc2_excluded_pathways.csv"
BRITE_CACHE = "kegg_brite_br08901_cache.json"

EXCLUDE_CATEGORIES = {"Human Diseases", "Organismal Systems", "Drug Development"}
KEGG_BRITE_URL = "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901&format=json&filedir="

# ---- 1) KEGG 공식 카테고리 데이터 받아오기 (한 번 받으면 로컬에 캐싱) ----
if os.path.exists(BRITE_CACHE):
    print(f"[정보] 로컬 캐시({BRITE_CACHE})에서 KEGG BRITE 데이터 로드")
    with open(BRITE_CACHE, encoding="utf-8") as f:
        brite = json.load(f)
else:
    print("[정보] KEGG BRITE(br08901) 공식 카테고리 데이터 요청 중...")
    resp = requests.get(KEGG_BRITE_URL, timeout=30)
    resp.raise_for_status()
    brite = resp.json()
    with open(BRITE_CACHE, "w", encoding="utf-8") as f:
        json.dump(brite, f)
    print(f"[정보] 로컬 캐시로 저장됨: {BRITE_CACHE} (다음부터는 재다운로드 없이 재사용)")

# ---- 2) pathway_id -> 대분류(top category) 매핑 만들기 ----
pathway_to_category = {}

def walk(node, top_category):
    name = node.get("name", "")
    # "00680  Methane metabolism" 형태에서 앞의 숫자(pathway id) 추출
    m = re.match(r"^(\d{5})\s", name)
    if m:
        pathway_to_category[m.group(1)] = top_category
    for child in node.get("children", []):
        walk(child, top_category)

for top_node in brite.get("children", []):
    top_name = top_node["name"]
    for child in top_node.get("children", []):
        walk(child, top_name)

print(f"[정보] 총 {len(pathway_to_category)}개 pathway의 카테고리 매핑 완료")

# ---- 3) ANCOM-BC2 결과 읽고 pathway_id에서 숫자만 추출 ----
df = pd.read_csv(IN_CSV)
df.columns = [c.strip().lower() for c in df.columns]


def extract_id(pid):
    m = re.search(r"(\d{5})", str(pid))
    return m.group(1) if m else None


df["pw_num"] = df["pathway_id"].apply(extract_id)
df["category"] = df["pw_num"].map(pathway_to_category)

unmatched = df[df["category"].isna()]
if len(unmatched) > 0:
    print(f"[안내] 카테고리를 못 찾은 pathway {len(unmatched)}개 (그대로 유지됨): "
          f"{unmatched['pathway_id'].tolist()[:10]}")

# ---- 4) 필터링 ----
before_sig = (df["q_value"] < 0.05).sum()
excluded = df[df["category"].isin(EXCLUDE_CATEGORIES)]
kept = df[~df["category"].isin(EXCLUDE_CATEGORIES)]
after_sig = (kept["q_value"] < 0.05).sum()

print(f"\n=== 필터링 결과 ===")
print(f"전체 pathway: {len(df)}개")
print(f"제외됨 (Human Diseases/Organismal Systems/Drug Development): {len(excluded)}개")
print(f"남음: {len(kept)}개")
print(f"\n필터 전 유의(q<0.05): {before_sig}개")
print(f"필터 후 유의(q<0.05): {after_sig}개")

# 방향성(AL vs PERM) 재확인
kept_sig = kept[kept["q_value"] < 0.05]
n_perm_direction = (kept_sig["log_fold_change"] > 0).sum()
n_al_direction = (kept_sig["log_fold_change"] < 0).sum()
print(f"\n필터 후 유의 pathway 방향: PERM 우세 {n_perm_direction}개 / AL 우세 {n_al_direction}개")

kept.drop(columns=["pw_num", "category"]).to_csv(OUT_CSV, index=False)
excluded.drop(columns=["pw_num"]).to_csv(OUT_EXCLUDED_CSV, index=False)

print(f"\n완료: {OUT_CSV} (필터링된 결과)")
print(f"완료: {OUT_EXCLUDED_CSV} (제외된 항목 확인용)")
