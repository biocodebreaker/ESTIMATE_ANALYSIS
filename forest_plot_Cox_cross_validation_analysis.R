# Cox proportional hazard regression analysis
#install.packages("forestplot")
#install.packages("glmnet")

library(forestplot)
library(glmnet)

# Exclude samples with missing survival data
complete_data <- merged_survival_info_estimate_scores_gene_expression_limma[!is.na(merged_survival_info_estimate_scores_gene_expression_limma$time), ]

# Remove samples with non-positive event times
complete_data <- complete_data[complete_data$time > 0, ]

# Extract the expression data for the 6 Stromal DEGs
Stromal_DEGs_expression <- complete_data[, c("OGN", "SFRP1", "THBS4", "CILP", "C7", "ADH1B")]

# Extract survival data
survival_data <- complete_data[, c("time", "event")]

# Set seed for reproducibility
set.seed(123)

# Function for cross-validation
cv_cox <- function(expr_data, surv_data, nfolds = 3, nrepeats = 1000) {
  n <- nrow(surv_data)
  foldid <- sample(rep(seq(nfolds), length.out = n))
  coefs <- matrix(NA, nrow = nrepeats, ncol = ncol(expr_data))
  
  for (i in 1:nrepeats) {
    fit <- cv.glmnet(x = as.matrix(expr_data), y = Surv(surv_data$time, surv_data$event), family = "cox", foldid = foldid, type.measure = "C")
    coefs[i, ] <- as.vector(coef(fit, s = "lambda.min")[-1]) # Convert to vector and exclude the intercept
  }
  
  return(coefs)
}

# Perform cross-validation
coefs <- cv_cox(Stromal_DEGs_expression, survival_data, nfolds = 3, nrepeats = 1000)


# Calculate the average coefficients
avg_coefs <- colMeans(coefs)

# Calculate the hazard ratios and their 95% confidence intervals
hr <- exp(avg_coefs)
lower_ci <- exp(avg_coefs - 1.96 * apply(coefs, 2, sd))
upper_ci <- exp(avg_coefs + 1.96 * apply(coefs, 2, sd))

# Combine the results into a data frame
hr_results <- data.frame(Gene = colnames(Stromal_DEGs_expression),
                         HR = hr,
                         LowerCI = lower_ci,
                         UpperCI = upper_ci)

# Create a forest plot

# Reverse the order of the rows to plot the genes from top to bottom
hr_results <- hr_results[nrow(hr_results):1, ]

# Generate the forest plot
forestplot(labeltext = list(hr_results$Gene, 
                            format(hr_results$HR, digits = 2),
                            format(hr_results$LowerCI, digits = 2),
                            format(hr_results$UpperCI, digits = 2)),
           mean = hr_results$HR,
           lower = hr_results$LowerCI,
           upper = hr_results$UpperCI,
           xlog = TRUE,
           is.summary = c(rep(FALSE, nrow(hr_results))),
           col = fpColors(box = "royalblue", lines = "darkblue", summary = "royalblue"))

