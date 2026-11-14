# Arctic Metagenome Analysis (Stage 2)

circum-Arctic 메타지놈 1차 파이프라인
([arctic-metagenome-pipeline](https://github.com/kimheejun99/Arctic_Metagenomics_Pipeline))이
샘플별로 만들어낸 산출물(Kraken2 report, VFDB annotation, ARGs-OAP 결과 등)을
**여러 샘플에 걸쳐 취합·비교·시각화**하는 2차 분석 파이프라인입니다.

## Pipeline overview

```
sample_coords.txt (+ sample_list.txt) + metadata.txt
   │
   ├─ Step0_make_googlemap.py    좌표 기반 구글맵 시각화 (좌표 파일 있을 때만)
   ├─ Step1_organize_samples.sh  flat한 산출물을 종류별 폴더로 정리
   ├─ Step2_vf_arg_plots.py      VF/ARG 전체 평균 + 그룹별 비교 (바그래프, 히트맵)
   ├─ Step3_host_tracing.py      VF/ARG 검출 read의 host phylum 추적 (Kraken2 조인)
   └─ Step4_merge_kreports.py    kreport 병합 -> ASV/taxonomy/metadata 테이블
```

## Repository structure

```
.
├── environment.yml           # conda 환경 (ArcticAnalysis)
├── .gitignore
├── run_pipeline.sh           # 전체 파이프라인 오케스트레이션
├── Step0_make_googlemap.py
├── Step1_organize_samples.sh
├── Step2_vf_arg_plots.py
├── Step3_host_tracing.py
├── Step4_merge_kreports.py
└── sample_coords.txt.example
```

## Requirements

1차 파이프라인과 달리 DIAMOND, Kraken2 같은 생물정보학 전용 바이너리는 필요 없습니다.
순수 bash + Python(pandas, matplotlib)만으로 동작합니다.

```bash
conda env create -f environment.yml
conda activate ArcticAnalysis
```

| 도구 | 용도 |
|---|---|
| Python 3.10 | Step0, 2, 3, 4 |
| pandas | tsv/kreport 파싱, 그룹 집계 |
| matplotlib | Step2 바그래프/히트맵, Step3 stacked bar |
| bash | run_pipeline.sh, Step1 |

## Input files

| 파일 | 필수 여부 | 설명 |
|---|---|---|
| `sample_coords.txt` | 필수 (탭 구분, 헤더 필수: `Run  latitude  longitude`) | Step0 구글맵용. Run 컬럼이 곧 처리할 샘플 목록으로도 쓰입니다. |
| `sample_list.txt` | 선택 | 없으면 `sample_coords.txt`의 Run 컬럼에서 자동 생성됩니다. 특정 서브셋만 돌리고 싶을 때만 직접 준비하세요. |
| `metadata.txt` | 선택 (탭 구분: `sample-id  group`) | 없으면 전체 샘플이 "All" 그룹 하나로 처리됩니다. |

예시는 `sample_coords.txt.example` 참고 (실제로 쓰실 땐 `sample_coords.txt`로 이름을 바꿔서 채워 넣으세요).

## 샘플별 원본 파일 (1차 파이프라인 산출물, 현재 폴더에 flat하게 위치)

**Step1 실행 전**에는 현재 폴더에 아래처럼 흩어져 있어야 합니다.

```
./ (현재 위치)
├── sample_coords.txt
├── (sample_list.txt)
├── metadata.txt
│
├── SRR0000001_1.fastq.gz
├── SRR0000001_2.fastq.gz
├── SRR0000001.kreport
├── SRR0000001.kraken.out
├── SRR0000001_annotated_parsed.tsv    <- Category 컬럼 필요
├── SRR0000001_args_output/
│   ├── normalized_cell.type.txt
│   ├── normalized_cell.subtype.txt
│   ├── normalized_cell.gene.txt
│   └── blastout.filtered.txt
│
└── (다른 샘플들도 동일한 패턴)
```

Step1 실행 후에는 종류별 폴더(`./kraken/`, `./vfdb/`, `./arg/`)로 정리됩니다
(fastq.gz는 이동하지 않고 그대로 둡니다). Step2~4는 이 폴더들을 참조합니다.

## 외부 데이터베이스

Step3, Step4는 taxid → phylum/lineage 변환을 위해 **NCBI taxdump**(`nodes.dmp`, `names.dmp`)가 필요합니다.

- 이미 갖고 계신 Kraken2 DB(PlusPFP-16 등) 안에 `taxonomy/nodes.dmp`, `taxonomy/names.dmp`로 포함되어 있을 가능성이 높습니다.
- 없다면: https://ftp.ncbi.nlm.nih.gov/pub/taxdump/ 에서 `taxdump.tar.gz` 다운로드 후 압축 해제

## Usage

```bash
conda activate ArcticAnalysis
bash run_pipeline.sh sample_coords.txt sample_list.txt metadata.txt
# 또는 sample_list.txt 없이:
bash run_pipeline.sh sample_coords.txt "" metadata.txt
```

환경변수로 override 가능한 값:

```bash
GOOGLE_MAPS_API_KEY=your_key \
TAXDUMP=/mnt/program/db/plusPFP_16/taxonomy \
JOBS=8 \
bash run_pipeline.sh sample_coords.txt sample_list.txt metadata.txt
```

| 환경변수 | 기본값 | 설명 |
|---|---|---|
| `GOOGLE_MAPS_API_KEY` | `YOUR_GOOGLE_MAPS_API_KEY` (placeholder) | Step0 지도가 실제로 로드되려면 유효한 키 필요 |
| `TAXDUMP` | `/mnt/program/db/plusPFP_16` | Step3/4용 nodes.dmp/names.dmp 상위 경로 |
| `JOBS` | `6` | Step3 병렬 처리 프로세스 수 |

## Outputs

| Step | 산출물 |
|---|---|
| Step0 | `sample_map.html` |
| Step1 | `./kraken/`, `./vfdb/`, `./arg/` |
| Step2 | `./plots/VF_total_category_barplot.png`, `VF_group_heatmap.png`, `ARG_total_type_barplot_*.png`, `ARG_group_heatmap_*.png` 등 |
| Step3 | `./host_tracing/VF_host_tracing_by_group.png`, `ARG_host_tracing_by_group.png`, raw tsv |
| Step4 | `./merged_tables/1_ASV_table.txt`, `2_taxonomy_table.txt`, `3_metadata.txt` |

## Notes

- `Step0_make_googlemap.py`는 Google Maps JavaScript API를 사용하므로, 실제 지도를 보려면
  [Google Cloud Console](https://console.cloud.google.com/)에서 Maps JavaScript API 키를 발급받아야 합니다.
  키 없이 실행하면 HTML 파일은 만들어지지만 지도가 로드되지 않습니다.
- 모든 입력 파일은 `.txt` 확장자로 통일했습니다. 실제로는 tab-separated 텍스트이고
  코드도 확장자가 아니라 `sep="\t"`(또는 Step0는 구분자 자동감지 `sep=None`)로 직접 파싱하기 때문에
  확장자 자체는 동작에 영향을 주지 않습니다.

---

# Stage 3: KEGG Functional Analysis (eggNOG-mapper 기반)

Stage 2의 Step4가 만든 `./merged_tables/{1_ASV,2_taxonomy}_table.txt`를 이어받아,
eggNOG-mapper 기반 KEGG pathway 기능 분석 및 다각도 통계 검증을 수행합니다.

```
merged_tables/{1_ASV,2_taxonomy}_table.txt (Stage 2 Step4 산출물)
   │
   ├─ Step5  eggNOG-mapper 설치 (1회성)
   ├─ Step6  49샘플 배치 실행 + 결과 정리
   ├─ Step7  KEGGaNOG 커뮤니티 완성도 분석 (heatmap, boxplot, FDR 검정)
   ├─ Step8  통계적 검증 - ANCOM-BC2(존재량), 카테고리 필터링, 20회 재검정
   ├─ Step9  Host-tracing - read 단위 taxon-기능 직접 연결
   ├─ Step10 심화 taxon 귀속 검증 (KEGG-Decoder 실제 KO 기준)
   ├─ Step11 교란요인 검증 - 시퀀싱 깊이 정규화(RPM)
   └─ Step12 Genus 수준 ANCOM-BC2 최종 검증 (Step4 산출물 재사용)
```

## Usage

```bash
conda activate eggnog          # eggNOG-mapper용 별도 환경
bash run_pipeline_kegg.sh sample_list.txt metadata.txt

# 특정 Step만 실행
bash run_pipeline_kegg.sh sample_list.txt metadata.txt 9
```

## Requirements (추가)

Stage 2의 environment.yml에 더해, 아래 두 conda 환경이 추가로 필요합니다.

| 환경 | 용도 |
|---|---|
| `eggnog` (Python 3.11, eggNOG-mapper, KEGGaNOG, seaborn) | Step5~7, 9~11 |
| `r_ancombc` (R, bioconductor-ancombc, phyloseq) | Step8, 12 (ANCOM-BC2 통계) |

## Outputs (추가)

| Step | 산출물 |
|---|---|
| Step7 | `heatmap_AL_PERM_sorted.png`, `methanogenesis_boxplot.png`, `all_pathways_stats.csv` |
| Step8 | `ancombc2_pathway_results.csv`, `ancombc2_volcano_plot.png`, `ancombc2_balanced_summary.csv` |
| Step9 | `host_tracing_master.tsv`, `top_taxa_by_group_barplot.png`, `methanogenesis_host_tracing_barplot.png` |
| Step10 | `robust_pathways_taxon_pattern.png` |
| Step11 | `sequencing_depth_per_sample.csv`, `host_tracing_normalized_comparison.png` |
| Step12 | `merged_tables/ancombc2_genus_corrected_results.csv` |

## Notes (추가)

- Step5~11은 Stage 1(1차 파이프라인)의 원본 fastq/`.kraken` 파일이 있어야 합니다(Step9 host-tracing).
- Step12는 Stage 2의 taxonomy 데이터 정합성이 사전에 검증되어 있어야 신뢰할 수 있습니다.
  (taxonomy 계통정보 오류가 있으면 통계 결과가 왜곡될 수 있으므로, `2_taxonomy_table.txt`의
  Kingdom 컬럼에 `Bacteria/Archaea/Eukaryota/Viruses` 외의 값이 없는지 사전 점검을 권장합니다.)
