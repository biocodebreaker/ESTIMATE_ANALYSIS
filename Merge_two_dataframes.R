#Merge two data.frames

merged_survival_info_time_event_estimate_scores$rown <- row.names(merged_survival_info_time_event_estimate_scores)
stromal_immune_optimal_cutoff_categorized$rown <- row.names(stromal_immune_optimal_cutoff_categorized)

final_data <- merge(merged_survival_info_time_event_estimate_scores, stromal_immune_optimal_cutoff_categorized, by = "rown")
