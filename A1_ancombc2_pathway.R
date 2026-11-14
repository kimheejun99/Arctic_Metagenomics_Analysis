# ==============================================================
# A1_ancombc2_pathway.R
# P2가 만든 pathway_abundance_wide_count.tsv(raw count)를 ANCOM-BC2로
# 검정하여, AL vs PERM 사이에 유의하게 다른 KEGG pathway를 찾습니다.
#
# 입력: pathway_abundance_wide_count.tsv (샘플 x pathway, raw count)
#       sample_group_mapping.csv (샘플 -> AL/PERM 그룹)
#       pathway_names_reference.tsv (pathway ID -> 이름)
#
# 사용법: Rscript A1_ancombc2_pathway.R
# ==============================================================

library(ANCOMBC)
library(phyloseq)
library(dplyr)

# ---- 1) 데이터 읽기 ----
count_wide <- read.delim("pathway_abundance_wide_count.tsv", row.names = 1,
                          check.names = FALSE)
groups <- read.csv("sample_group_mapping.csv")

# [확인 필요] sample_group_mapping.csv의 실제 컬럼명이 다르면 여기를 맞춰주세요.
# head(groups) 로 먼저 확인하시고, 아래 두 컬럼명이 실제와 다르면 수정하세요.
cat("sample_group_mapping.csv 컬럼:", paste(colnames(groups), collapse=", "), "\n")
SAMPLE_COL <- "Sample"   # 실제 컬럼명으로 수정
GROUP_COL  <- "Group"    # 실제 컬럼명으로 수정 (AL/PERM)

groups <- groups[, c(SAMPLE_COL, GROUP_COL)]
colnames(groups) <- c("sample", "group")
rownames(groups) <- groups$sample

# ---- 2) 샘플 순서 맞추기 (교집합만 사용) ----
common_samples <- intersect(rownames(count_wide), rownames(groups))
cat("공통 샘플 수:", length(common_samples), "/ 예상 49\n")
if (length(common_samples) < 49) {
  cat("⚠️  경고: 49개보다 적습니다. 샘플명 표기가 count_wide와 groups에서 다를 수 있습니다.\n")
}

count_wide <- count_wide[common_samples, ]
groups <- groups[common_samples, ]

# ---- 3) phyloseq 객체 구성 (rows=pathway, cols=sample 형태로 전치) ----
otu <- otu_table(t(as.matrix(count_wide)), taxa_are_rows = TRUE)
meta <- sample_data(groups)
physeq <- phyloseq(otu, meta)

# ---- 4) ANCOM-BC2 실행 (AL을 참조군으로) ----
out <- ancombc2(
  data = physeq,
  fix_formula = "group",
  p_adj_method = "fdr",
  group = "group",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05
)

res <- out$res

# ---- 5) 결과 정리 + pathway 이름 붙이기 ----
names_ref <- read.delim("pathway_names_reference.tsv")

# group 컬럼 이름은 ancombc2 결과에서 "lfc_group<참조군이 아닌 쪽 레벨>" 형태로 나옵니다.
# 실제 컬럼명은 str(res) 또는 colnames(res)로 확인 후 아래 grep 패턴이 맞는지 점검하세요.
cat("ANCOM-BC2 결과 컬럼:", paste(colnames(res), collapse=", "), "\n")

lfc_col <- grep("^lfc_group", colnames(res), value = TRUE)[1]
q_col   <- grep("^q_group",   colnames(res), value = TRUE)[1]

res_named <- res %>%
  select(taxon, all_of(lfc_col), all_of(q_col)) %>%
  rename(pathway_id = taxon, log_fold_change = all_of(lfc_col), q_value = all_of(q_col)) %>%
  left_join(names_ref, by = "pathway_id") %>%
  arrange(q_value)

write.csv(res_named, "ancombc2_pathway_results.csv", row.names = FALSE)

sig <- res_named %>% filter(q_value < 0.05)
cat("\n=== 유의한(q < 0.05) pathway 수:", nrow(sig), "/ 전체", nrow(res_named), "===\n")
print(sig)

cat("\n결과 저장됨: ancombc2_pathway_results.csv\n")
