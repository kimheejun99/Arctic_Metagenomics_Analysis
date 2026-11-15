# ==============================================================
# A2_ancombc2_balanced_subsample.R
# AL(36개)과 PERM(13개)의 표본 수 차이가 ANCOM-BC2 결과를 왜곡하는지
# 확인하기 위해, AL에서 13개를 무작위로 뽑아 PERM(13개)과 균형을 맞춘 뒤
# ANCOM-BC2를 재실행합니다. 우연한 추출의 영향을 배제하기 위해 20회 반복합니다.
#
# 출력: ancombc2_balanced_summary.csv
#   각 pathway가 20회 중 몇 번 유의(q<0.05)했는지, 그 방향(AL/PERM)이
#   얼마나 일관됐는지 정리한 표
#
# 사용법: /root/miniforge3/envs/r_ancombc/bin/Rscript A2_ancombc2_balanced_subsample.R
# ==============================================================

library(ANCOMBC)
library(phyloseq)
library(dplyr)

N_ITER <- 20        # 반복 횟수
SET_SEED_BASE <- 42 # 재현성을 위한 시드

# ---- 1) 데이터 읽기 (A1과 동일) ----
count_wide <- read.delim("pathway_abundance_wide_count.tsv", row.names = 1,
                          check.names = FALSE)
groups <- read.csv("sample_group_mapping.csv")

SAMPLE_COL <- "Sample"
GROUP_COL  <- "Group"
groups <- groups[, c(SAMPLE_COL, GROUP_COL)]
colnames(groups) <- c("sample", "group")
rownames(groups) <- groups$sample

common_samples <- intersect(rownames(count_wide), rownames(groups))
count_wide <- count_wide[common_samples, ]
groups <- groups[common_samples, ]

al_samples <- rownames(groups)[groups$group == "AL"]
perm_samples <- rownames(groups)[groups$group == "PERM"]
n_perm <- length(perm_samples)
cat("AL 전체:", length(al_samples), "개 / PERM 전체:", n_perm, "개\n")
cat("AL에서 매회", n_perm, "개씩 무작위 추출하여", N_ITER, "회 반복합니다.\n\n")

names_ref <- read.delim("pathway_names_reference.tsv")

# ---- 2) N_ITER회 반복 ----
all_results <- list()

for (i in 1:N_ITER) {
  set.seed(SET_SEED_BASE + i)
  al_sub <- sample(al_samples, n_perm)
  use_samples <- c(al_sub, perm_samples)

  sub_count <- count_wide[use_samples, ]
  sub_groups <- groups[use_samples, ]

  otu <- otu_table(t(as.matrix(sub_count)), taxa_are_rows = TRUE)
  meta <- sample_data(sub_groups)
  physeq <- phyloseq(otu, meta)

  out <- tryCatch({
    ancombc2(
      data = physeq, fix_formula = "group", p_adj_method = "fdr",
      group = "group", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05
    )
  }, error = function(e) {
    cat("  [반복", i, "] 에러 발생, 건너뜀:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(out)) next

  res <- out$res
  lfc_col <- grep("^lfc_group", colnames(res), value = TRUE)[1]
  q_col   <- grep("^q_group",   colnames(res), value = TRUE)[1]

  iter_res <- res %>%
    select(taxon, all_of(lfc_col), all_of(q_col)) %>%
    rename(pathway_id = taxon, log_fold_change = all_of(lfc_col), q_value = all_of(q_col)) %>%
    mutate(iter = i, significant = q_value < 0.05)

  all_results[[i]] <- iter_res
  cat("[반복", i, "/", N_ITER, "] 완료. 유의 pathway 수:", sum(iter_res$significant), "\n")
}

# ---- 3) 결과 취합: pathway별로 20회 중 몇 번 유의했는지 ----
combined <- bind_rows(all_results)

summary_df <- combined %>%
  group_by(pathway_id) %>%
  summarise(
    times_significant = sum(significant),
    times_tested = n(),
    pct_significant = round(100 * times_significant / times_tested, 1),
    mean_lfc = round(mean(log_fold_change, na.rm = TRUE), 3),
    direction = ifelse(mean_lfc > 0, "PERM", "AL")
  ) %>%
  left_join(names_ref, by = "pathway_id") %>%
  arrange(desc(times_significant))

write.csv(summary_df, "ancombc2_balanced_summary.csv", row.names = FALSE)

cat("\n=== 완료: ancombc2_balanced_summary.csv ===\n")
cat("\n[20회 중 15회(75%) 이상 유의했던 pathway - 강건한 결과]\n")
robust <- summary_df %>% filter(times_significant >= 15)
print(robust)

cat("\n[메탄생성 관련 pathway 확인]\n")
methane <- summary_df %>% filter(grepl("00680|00710|00720", pathway_id))
print(methane)
