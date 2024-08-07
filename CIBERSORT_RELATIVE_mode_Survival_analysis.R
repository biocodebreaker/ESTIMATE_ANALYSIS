#Survival analysis of CIBERSORT results (LM22 Immune Cell Type expression values)


#Install maxstat and the required dependencies 

install.packages("exactRankTests")
install.packages("mvtnorm")

install.packages("/cloud/home/r1816512/ESTIMATE_ANALYSIS/maxstat_0.7-25.tar.gz", repos = NULL, type = "source", dependencies = TRUE)


# Read the .csv file into a data.frame
cibersortx_relative_mode_results_df <- read.csv("CIBERSORTx_RELATIVE_mode_Results_maxstat_survival_analysis.csv")

#Format the submitter_ids so that they match 
cibersortx_relative_mode_results_df$submitter_id <- sub("-(?!.*-).*", "", cibersortx_relative_mode_results_df$submitter_id, perl = TRUE)

#Merge survival info and cibersort results
merged_cibersort_survival_info_estimate_scores_df <- merge(merged_survival_info_censored_estimate_scores, cibersortx_relative_mode_results_df, by = "submitter_id", all.x = FALSE)


# Create the Surv object
cibersort_surv_obj <- with(merged_cibersort_survival_info_estimate_scores_df, Surv(time, event))


# Find the optimal cutoff values
# Find the optimal cutoff values
maxstat_cibersort <- maxstat.test(cibersort_surv_obj ~ merged_cibersort_survival_info_estimate_scores_df, data = merged_cibersort_survival_info_estimate_scores_df, smethod = "LogRank")


optimal_cutoff_stromal <- maxstat_stromal$statistic[1]