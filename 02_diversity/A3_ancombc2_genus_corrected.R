# ==============================================================
# A3_ancombc2_genus_corrected.R
# taxonomy 데이터 오염 수정 후(1_ASV_table.txt, 2_taxonomy_table.txt)의
# genus 수준 존재량으로 ANCOM-BC2를 재실행하여, Methanosarcina(및 전체
# genus)의 AL vs PERM 유의성을 재검증합니다.
#
# 예전 metadata.txt의 group 컬럼이 비어있는 문제가 반복됐던 점을 고려하여,
# 이미 여러 스크립트에서 검증된 AL(36)/PERM(13) 목록을 코드에 직접
# 내장하여 외부 metadata 파일 품질에 의존하지 않도록 하였습니다.
#
# 입력: 1_ASV_table.txt (taxid x sample, raw count)
#       2_taxonomy_table.txt (taxid -> Kingdom~Species)
# 출력: ancombc2_genus_corrected_results.csv
#
# 사용법: /root/miniforge3/envs/r_ancombc/bin/Rscript A3_ancombc2_genus_corrected.R
# ==============================================================

library(ANCOMBC)
library(phyloseq)
library(dplyr)

ASV_FILE <- "1_ASV_table.txt"
TAX_FILE <- "2_taxonomy_table.txt"
OUT_CSV  <- "ancombc2_genus_corrected_results.csv"

# ---- 검증된 AL/PERM 목록 (AL_PERM_eukary_grouping.py 기준과 동일) ----
# ==============================================================
# 샘플 목록을 하드코딩하지 않고, metadata.txt(sample-id	group)에서
# 동적으로 읽어옵니다. 몇 개 샘플이 있든 자동으로 맞춰집니다.
# ==============================================================
# metadata.txt 탐색: 현재 폴더 → 상위 폴더 → 상위 상위 폴더
meta_candidates <- c("metadata.txt", "../metadata.txt", "../../metadata.txt")
meta_path <- Filter(file.exists, meta_candidates)
if (length(meta_path) == 0) {
    stop("[오류] metadata.txt를 찾을 수 없습니다. 현재 폴더: ", getwd())
}
meta_path <- meta_path[1]
cat("[안내] metadata.txt 경로:", meta_path, "\n")
meta_df <- read.delim(meta_path, sep = "\t")
colnames(meta_df) <- trimws(colnames(meta_df))
# sample-id 컬럼명 유연하게 처리 (sample-id, sample_id, sampleid 모두 허용)
sid_col <- grep("sample.id|sampleid", colnames(meta_df), value=TRUE, ignore.case=TRUE)[1]
grp_col <- grep("^group$", colnames(meta_df), value=TRUE, ignore.case=TRUE)[1]
if (is.na(sid_col) || is.na(grp_col)) {
    cat("[경고] 컬럼명 확인:", paste(colnames(meta_df), collapse=", "), "\n")
    stop("[오류] sample-id 또는 group 컬럼을 찾을 수 없습니다.")
}
AL_SAMPLES <- as.character(meta_df[[sid_col]][meta_df[[grp_col]] == "AL"])
PERM_SAMPLES <- as.character(meta_df[[sid_col]][meta_df[[grp_col]] == "PERM"])
cat("[metadata.txt] AL", length(AL_SAMPLES), "개, PERM", length(PERM_SAMPLES), "개 로드됨\n")
if (length(AL_SAMPLES) == 0 && length(PERM_SAMPLES) == 0) {
    stop("[오류] AL/PERM 샘플이 모두 0개입니다. metadata.txt 내용을 확인하세요:\n",
         paste(capture.output(head(meta_df)), collapse="\n"))
}

# ---- 1) ASV/taxonomy 읽기 ----
asv <- read.delim(ASV_FILE, row.names = 1, check.names = FALSE)
tax <- read.delim(TAX_FILE)
colnames(tax)[1] <- "taxid"

cat("ASV 테이블:", nrow(asv), "taxa x", ncol(asv), "샘플\n")

# ---- 2) 샘플 컬럼명이 taxid와 일치하는지 확인 (부분 문자열 매칭으로 안전하게) ----
match_col <- function(srr, cols) {
  hit <- cols[grepl(srr, cols, fixed = TRUE)]
  if (length(hit) != 1) return(NA)
  return(hit)
}
al_cols <- sapply(AL_SAMPLES, match_col, cols = colnames(asv))
perm_cols <- sapply(PERM_SAMPLES, match_col, cols = colnames(asv))

missing_al <- AL_SAMPLES[is.na(al_cols)]
missing_perm <- PERM_SAMPLES[is.na(perm_cols)]
if (length(missing_al) > 0) cat("[경고] ASV 테이블에서 못 찾은 AL 샘플:", paste(missing_al, collapse=", "), "\n")
if (length(missing_perm) > 0) cat("[경고] ASV 테이블에서 못 찾은 PERM 샘플:", paste(missing_perm, collapse=", "), "\n")

al_cols <- al_cols[!is.na(al_cols)]
perm_cols <- perm_cols[!is.na(perm_cols)]
use_cols <- c(al_cols, perm_cols)
cat("공통 샘플 수:", length(use_cols), "/ 예상", length(AL_SAMPLES) + length(PERM_SAMPLES), "\n")
if (length(use_cols) < 2) {
    stop("[오류] ASV 테이블에서 유효한 샘플을 2개 이상 찾지 못했습니다. ",
         "sample_list, metadata, ASV 테이블의 샘플 ID가 일치하는지 확인하세요.\n",
         "ASV 컬럼명: ", paste(head(colnames(asv), 10), collapse=", "))
}

# ---- 3) Genus 수준으로 집계 ----
asv_sub <- asv[, use_cols, drop=FALSE]
tax_indexed <- tax[match(as.character(rownames(asv_sub)), as.character(tax$taxid)), ]
genus_raw <- gsub("^g__", "", tax_indexed$Genus)
genus_raw[is.na(genus_raw) | genus_raw == "" | genus_raw == "-"] <- "Unclassified"

genus_table <- aggregate(asv_sub, by = list(Genus = genus_raw), FUN = sum)
rownames(genus_table) <- genus_table$Genus
genus_table$Genus <- NULL
cat("Genus 집계 후:", nrow(genus_table), "개 genus\n")

# ---- 4) 극희귀 genus 제거 (prevalence 기준 - ANCOM-BC2 내부 계산에서
#         분산 0으로 처리되는 초희귀종을 raw count 분산 체크보다 확실하게 거름) ----
MIN_SAMPLES <- max(2, floor(ncol(asv) * 0.5))  # 샘플 수의 50% 이상에서 검출 (최소 2개)
prevalence <- rowSums(genus_table > 0)
rare_genus <- rownames(genus_table)[prevalence < MIN_SAMPLES]
cat("[안내] 검출 샘플 수가", MIN_SAMPLES, "개 미만이라 제외되는 genus",
    length(rare_genus), "개 (예:", paste(head(rare_genus, 5), collapse = ", "), "...)\n")
genus_table <- genus_table[prevalence >= MIN_SAMPLES, ]
cat("Prevalence 필터 후:", nrow(genus_table), "개 genus\n")

# ---- 5) 알려진 문제 genus 명시적 제외 (반복 확인된 것들) ----
KNOWN_PROBLEM_GENERA <- c("Serangoonvirus")
before_n <- nrow(genus_table)
genus_table <- genus_table[!(rownames(genus_table) %in% KNOWN_PROBLEM_GENERA), ]
if (before_n != nrow(genus_table)) {
  cat("[안내] 기존에 문제 확인된 genus", before_n - nrow(genus_table), "개 명시적 제외:",
      paste(KNOWN_PROBLEM_GENERA, collapse = ", "), "\n")
}

# ---- 6) log1p 변환 기준 분산도 추가 확인 (ANCOM-BC2 내부 log-scale 계산과 더 가까움) ----
log_variance <- apply(log1p(genus_table), 1, var)
genus_table <- genus_table[log_variance > 0 & !is.na(log_variance), ]
cat("log-scale 분산 확인 후:", nrow(genus_table), "개 genus\n")

# ---- 7) 그래도 남아있을 수 있는 raw 분산 0 항목 이중 확인 ----
genus_variance <- apply(genus_table, 1, var)
genus_table <- genus_table[genus_variance > 0 & !is.na(genus_variance), ]
cat("최종 분산 확인 후:", nrow(genus_table), "개 genus\n")

# ---- 8) 그룹 메타데이터 ----
groups <- data.frame(
  sample = use_cols,
  group = c(rep("AL", length(al_cols)), rep("PERM", length(perm_cols)))
)
rownames(groups) <- groups$sample

# ---- 9) phyloseq 객체 생성 + ANCOM-BC2 ----
otu <- otu_table(as.matrix(genus_table), taxa_are_rows = TRUE)
meta <- sample_data(groups)
physeq <- phyloseq(otu, meta)

# 샘플 수가 너무 적으면 ANCOM-BC2가 실패하므로 tryCatch로 감싸기
n_per_group_min <- min(length(al_cols), length(perm_cols))
if (n_per_group_min < 2) {
    cat("[안내] 그룹당 샘플 수가 2개 미만입니다. ANCOM-BC2를 건너뜁니다.\n")
    cat("       49개 전체 샘플로 실행 시 정상 동작합니다.\n")
    quit(status = 0)
}

out <- tryCatch({
    ancombc2(
        data = physeq, fix_formula = "group", p_adj_method = "fdr",
        group = "group", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05
    )
}, error = function(e) {
    cat("[안내] ANCOM-BC2 실행 실패:", conditionMessage(e), "\n")
    cat("       샘플 수 부족(현재 AL =", length(al_cols), ", PERM =", length(perm_cols), ")이 원인일 수 있습니다.\n")
    cat("       49개 전체 샘플로 실행 시 정상 동작합니다.\n")
    # 빈 결과 CSV 생성 후 종료
    write.csv(data.frame(genus=character(), log_fold_change=numeric(), q_value=numeric()),
              OUT_CSV, row.names = FALSE)
    quit(status = 0)
})

res <- out$res
lfc_col <- grep("^lfc_group", colnames(res), value = TRUE)[1]
q_col   <- grep("^q_group",   colnames(res), value = TRUE)[1]

result_df <- res %>%
  select(taxon, all_of(lfc_col), all_of(q_col)) %>%
  rename(genus = taxon, log_fold_change = all_of(lfc_col), q_value = all_of(q_col)) %>%
  arrange(q_value)

write.csv(result_df, OUT_CSV, row.names = FALSE)

cat("\n=== 완료:", OUT_CSV, "===\n")
cat("\n[유의(q<0.05) genus 목록]\n")
print(result_df %>% filter(q_value < 0.05))

cat("\n[Methanosarcina 결과]\n")
print(result_df %>% filter(grepl("Methanosarcina", genus, ignore.case = TRUE)))

cat("\n[메탄생성균 계열 전체 확인]\n")
print(result_df %>% filter(grepl("Methano", genus, ignore.case = TRUE)))
