#Analysis of the CIBERSORTx ABSOLUTE mode RESULTS


# Job Parameters used for this run:
  
# Date: 2023-11-24 00:25:08
# Job type: Impute Cell Fractions
# Signature matrix file: LM22.update-gene-symbols.txt
# Mixture file: filtered_CIBERSORT_INPUT_RNASEQ_TUMOR_ONLY_522.txt
# Batch correction: enabled
# Batch correction mode: B-mode
# Disable quantile normalization: true
# Run mode (relative or absolute): absolute
# Permutations: 100

#Read the CIBERSORT ABSOLUTE mode results.csv into a df
cibersortx_ABSOLUTE_cell_proportions_tumor_only_df <- read.csv("cibersortx_impute_cell_fraction_ABS_mode_tumor_only.csv")
dim(cibersortx_ABSOLUTE_cell_proportions_tumor_only_df)
#[1] 411   2

#Format the submitter_ids so they match the ones in survival_info

cibersortx_ABSOLUTE_cell_proportions_tumor_only_df$submitter_id <- sub("\\.01A$", "", cibersortx_ABSOLUTE_cell_proportions_tumor_only_df$submitter_id)

cibersortx_ABSOLUTE_cell_proportions_tumor_only_df$submitter_id <- gsub("\\.", "-", cibersortx_ABSOLUTE_cell_proportions_tumor_only_df$submitter_id)


surv_info <- merged_survival_info_estimate_scores_gene_expression_limma[, c("submitter_id", "time", "event")]

# Read in file with 448 TCGA samples 

merged_survival_info_time_event_estimate_scores <- read.table("merged_survival_info_time_event_estimate_scores.txt", header = TRUE, sep = "\t")

# Filter to remove rows with NA values 

merged_survival_info_time_event_estimate_scores <- na.omit(merged_survival_info_time_event_estimate_scores, cols = "time")

dim(merged_survival_info_time_event_estimate_scores)
#[1] 439   6

surv_info_439 <- merged_survival_info_time_event_estimate_scores[, c("submitter_id", "time", "event")]


library(dplyr)

#Merge surv_info with cibersortx_ABSOLUTE_cell_proportions_tumor_only_df
cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info <- cibersortx_ABSOLUTE_cell_proportions_tumor_only_df %>%
  inner_join(surv_info, by = "submitter_id")

dim(cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info)
#[1] 302   4

cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info_439 <- cibersortx_ABSOLUTE_cell_proportions_tumor_only_df %>%
  inner_join(surv_info_439, by = "submitter_id")

dim(cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info_439)
#[1] 311   4



#Filter the df so that you can input into maxstat. 

#cibersort_ABSOLUTE_survival_info_filtered_df <- merged_cibersort_ABSOLUTE_survival_info_df[, c(1:6, ncol(merged_cibersort_ABSOLUTE_survival_info_df))]

library(survival)
library(maxstat)

# Create the Surv object
surv_obj <- with(cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info, Surv(time, event))

surv_obj_439 <- with(cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info_439, Surv(time, event))

# Find the optimal cutoff values using maxstat
maxstat_cibersort_ABSOLUTE_cutoff <- maxstat.test(surv_obj ~ Absolute.score, data = cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info, smethod = "LogRank")
#estimated cutpoint 1.559581

#optimal_cutoff_maxstat_cibersort_ABSOLUTE <- maxstat_cibersort_ABSOLUTE$statistic[1]
optimal_cutoff_maxstat_cibersort_ABSOLUTE <- 1.559581

maxstat_cibersort_ABSOLUTE_cutoff_439 <- maxstat.test(surv_obj_439 ~ Absolute.score, data = cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info_439, smethod = "LogRank")
# estimated cutpoint 1.559581
maxstat_cibersort_ABSOLUTE_cutoff_439 <- 1.559581

library(dplyr)
# Add column `ABS_TIL_levels` categorized into "high" vs "low"
cibersort_ABSOLUTE_TIL_levels <- cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info %>%
  mutate(ABS_TIL_levels = ifelse(Absolute.score > 1.559581, "high", "low"))

table(cibersort_ABSOLUTE_TIL_levels$ABS_TIL_levels)

# high  low 
# 266   36

cibersort_ABSOLUTE_TIL_levels_439 <- cibersortx_ABSOLUTE_cell_proportions_tumor_only_surv_info_439 %>%
  mutate(ABS_TIL_levels = ifelse(Absolute.score > 1.559581, "high", "low"))

table(cibersort_ABSOLUTE_TIL_levels_439$ABS_TIL_levels)

# high  low 
# 273   38

library(survminer)
# Fit Kaplan-Meier survival curves for the cibersort_ABSOLUTE_TIL_levels
surv_obj_311 <- with(cibersort_ABSOLUTE_TIL_levels_439, Surv(time, event))

km_fit_cibersort_ABSOLUTE_TIL_levels_439 <- survfit(surv_obj_439 ~ ABS_TIL_levels, data = cibersort_ABSOLUTE_TIL_levels_439)
ggsurvplot(km_fit_cibersort_ABSOLUTE_TIL_levels_439, data = cibersort_ABSOLUTE_TIL_levels_439, pval = TRUE, risk.table = TRUE)


surv_obj_302 <- with(cibersort_ABSOLUTE_TIL_levels, Surv(time, event))

km_fit_cibersort_ABSOLUTE_TIL_levels <- survfit(surv_obj_302 ~ ABS_TIL_levels, data = cibersort_ABSOLUTE_TIL_levels)

ggsurvplot(km_fit_cibersort_ABSOLUTE_TIL_levels, data = cibersort_ABSOLUTE_TIL_levels, pval = TRUE, risk.table = TRUE, dpi = 300)

