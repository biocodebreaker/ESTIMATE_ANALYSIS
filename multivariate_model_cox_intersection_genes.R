# Prepare the data
complete_data <- merged_survival_info_estimate_scores_gene_expression_limma[!is.na(merged_survival_info_estimate_scores_gene_expression_limma$time), ]
complete_data <- complete_data[complete_data$time > 0, ]

# Extract the expression data for the genes in intersection_genes
intersection_genes_expression <- complete_data[, intersection_genes]

# Extract survival data
survival_data <- complete_data[, c("time", "event")]

# Create the Surv() object
surv_obj <- Surv(survival_data$time, survival_data$event)

# Fit the multivariate Cox proportional hazards model
multivariate_cox_model <- coxph(surv_obj ~ ., data = intersection_genes_expression)  # Only include predictor variables

# Display the model summary
summary(multivariate_cox_model)

# Extract the expression data for the 6 Stromal DEGs
Stromal_DEGs_expression <- complete_data[, c("OGN", "SFRP1", "THBS4", "CILP", "C7", "ADH1B")]

# Extract survival data
survival_data <- complete_data[, c("time", "event")]

# Create the Surv() object
surv_obj <- Surv(survival_data$time, survival_data$event)

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ OGN + SFRP1 + THBS4 + CILP + C7 + ADH1B, data = complete_data)

# Display the model summary
summary(cox_model)
