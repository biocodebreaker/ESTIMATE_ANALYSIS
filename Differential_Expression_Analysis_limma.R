# Read in the data

survival_info_estimate_scores_gene_expression_limma <- read.csv("survival_info_estimate_scores_gene_expression_limma.csv", header = TRUE)

dim(survival_info_estimate_scores_gene_expression_limma)
# [1]  381 9876

merged_survival_info_estimate_scores_gene_expression_limma <- survival_info_estimate_scores_gene_expression_limma


# Find sizes of stromal_group and immune_group
table(merged_survival_info_estimate_scores_gene_expression_limma$Stromal_group)

# high  low 
# 42  339

# table(merged_survival_info_estimate_scores_gene_expression_limma$Immune_group)

# high  low 
# 65  316 

## create the design matrix for both stromal and immune groups

Stromal_design <- data.frame(
  Stromallow = as.numeric(merged_survival_info_estimate_scores_gene_expression_limma$Stromal_group == "low"),
  Stromalhigh = as.numeric(merged_survival_info_estimate_scores_gene_expression_limma$Stromal_group == "high")
)
rownames(Stromal_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id

Immune_design <- data.frame(
  Immunelow = as.numeric(merged_survival_info_estimate_scores_gene_expression_limma$Immune_group == "low"),
  Immunehigh = as.numeric(merged_survival_info_estimate_scores_gene_expression_limma$Immune_group == "high")
)
rownames(Immune_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id


# Create a contrast matrix for the Stromal_group
Stromal_contrasts <- makeContrasts(high_vs_low = Stromalhigh - Stromallow, levels = Stromal_design)

# Create a contrast matrix for the Immune_group
Immune_contrasts <- makeContrasts(high_vs_low = Immunehigh - Immunelow, levels = Immune_design)

#set the rownames of the design matrices to be the same as the submitter_id column
#in the merged_survival_info_estimate_scores_gene_expression_limma

#rownames(Stromal_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id
#rownames(Immune_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id

#make sure that the gene expression data also has the submitter_id as its rownames
gene_expression_data <- merged_survival_info_estimate_scores_gene_expression_limma[, 12:ncol(merged_survival_info_estimate_scores_gene_expression_limma)]
rownames(gene_expression_data) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id

#The lmFit function expects the number of columns in the data to match the number of rows in the design matrix.

# Transpose the gene expression data
transposed_gene_expression_data <- t(gene_expression_data)

# Run lmFit for the Stromal_group
Stromal_fit <- lmFit(transposed_gene_expression_data, Stromal_design)

# Apply contrasts and compute empirical Bayes statistics
Stromal_fit <- contrasts.fit(Stromal_fit, Stromal_contrasts)
Stromal_fit <- eBayes(Stromal_fit)

# Linear modeling and empirical Bayes moderation for Immune_group
Immune_fit <- lmFit(transposed_gene_expression_data, Immune_design)
Immune_fit <- contrasts.fit(Immune_fit, Immune_contrasts)
Immune_fit <- eBayes(Immune_fit)

# Find top DEGs for Stromal_group
TCGA_Stromal_DEGs <- topTable(Stromal_fit, n = Inf, adjust.method = "fdr")

# Find top DEGs for Immune_group
TCGA_Immune_DEGs <- topTable(Immune_fit, n = Inf, adjust.method = "fdr")

# Get all Significant Immune DEGs
TCGA_Significant_Immune_DEGs_adj.P.Val_0.05 <- TCGA_Immune_DEGs[TCGA_Immune_DEGs$adj.P.Val < 0.05, ]

# Get all Significant Stromal DEGs

TCGA_Significant_Stromal_DEGs_adj.P.Val_0.05 <- TCGA_Stromal_DEGs[TCGA_Stromal_DEGs$adj.P.Val < 0.05, ]

#genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

#genes_to_check %in% rownames(TCGA_Significant_Stromal_DEGs_adj.P.Val_0.05)

#genes_to_check %in% rownames(TCGA_Significant_Immune_DEGs_adj.P.Val_0.05)

# Survival Analysis of TCGA Data

surv_obj <- Surv(merged_survival_info_estimate_scores_gene_expression_limma$time, merged_survival_info_estimate_scores_gene_expression_limma$event)


# Initialize univariate_results data.frame

TCGA_Significant_intersection_genes_univariate_results <- list()

# Perform univariate Cox regression for each gene in TCGA_Significant_intersection_genes
for (gene in TCGA_Significant_intersection_genes) {
  cox_model <- coxph(surv_obj ~ get(gene), data = merged_survival_info_estimate_scores_gene_expression_limma)
  result <- summary(cox_model)
  TCGA_Significant_intersection_genes_univariate_results[[gene]] <- result
}

# Multivariate analysis 
dim(TCGA_univariate_significant_genes_expression_data)
# [1] 381 878

TCGA_univariate_significant_genes_multivariate_results <- coxph(surv_obj ~ ., data = TCGA_univariate_significant_genes_expression_data[, 1:439])

#length(TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs)
#[1] 44

TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs <- intersect(TCGA_univariate_significant_genes, GSE84426_GSE84433_univariate_significant_genes)

TCGA_GSE84426_GSE84433_significant_intersection_DEGs <- intersect(TCGA_Significant_intersection_genes, GSE84426_GSE84433_intersection_genes)

length(TCGA_GSE84426_GSE84433_significant_intersection_DEGs)
#[1] 1253

# MULTIVARIATE ANALYSIS OF THE 44 TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs


TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs_expression_data <- merged_survival_info_estimate_scores_gene_expression_limma[, TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs]

surv_obj <- Surv(merged_survival_info_estimate_scores_gene_expression_limma$time, merged_survival_info_estimate_scores_gene_expression_limma$event)

#multivariate_cox_model <- coxph(surv_obj ~ ., data = TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs_expression_data)

TCGA_GSE84426_GSE84433_multivariate_results <- coxph(surv_obj ~ ., data = TCGA_GSE84426_GSE84433_univariate_significant_intersection_DEGs_expression_data)



# Display the summary results for the first 40 genes
for (gene in intersection_genes[60:90]) {
  cat("Results for", gene, "\n")
  print(univariate_results[[gene]])
}




# Filter DEGs based on FDR-adjusted p-value and log fold-change thresholds
Stromal_DEGs <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.0001 & abs(Stromal_topTable$logFC) > 4, ]
Stromal_DEGs_p_value_0.05_logFC_1 <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) > 1, ]

Immune_DEGs <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) > 1, ]

# Overexpressed genes in the Immune group (positive logFC values)
Immune_overexpressed <- Immune_DEGs[Immune_DEGs$logFC > 0, ]

# Underexpressed genes in the Immune group (negative logFC values)
Immune_underexpressed <- Immune_DEGs[Immune_DEGs$logFC < 0, ]


#Find overlapping Genes between Stromal and Immune Groups 

#Find underExpressed Stromal DEGs (0)
UnderExpressed_Stromal_DEGs <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) < -1.5, ]

#Find overExpressed Stromal DEGs (863)
Overexpressed_Stromal_DEGs <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) > 1.5, ]

#Overexpressed Immune DEGs (325)
OverExpressed_Immune_DEGs <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) > 1.5, ]

#UnderExpressed Immune DEGs (0)
UnderExpressed_Immune_DEGs <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) < -1.5, ]

#Generate a heatmap with unsupervised hierarchical clustering for the Stromal_group DEGs
library(pheatmap)
# Extract gene expression values for the selected DEGs
Stromal_DEGs_expression <- transposed_gene_expression_data[rownames(transposed_gene_expression_data) %in% rownames(Stromal_DEGs), ]
Immune_DEGs_expression <- transposed_gene_expression_data[rownames(transposed_gene_expression_data) %in% rownames(Immune_DEGs), ]

# Scale the expression values by gene (row)
Stromal_DEGs_scaled <- t(scale(t(Stromal_DEGs_expression)))
Immune_DEGs_scaled <- t(scale(t(Immune_DEGs_expression)))

# Create a heatmap with hierarchical clustering

#Produce Heatmap of the DEGs with gene names
pheatmap(Stromal_DEGs_scaled, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         scale = "none", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         annotation_row = data.frame(Gene = rownames(Stromal_DEGs_scaled), row.names = rownames(Stromal_DEGs_scaled)))

