#convert the Ensembl gene IDs to gene symbols

#gene_id to gene_name mapping


# Remove version information from gene_id
Gene_ID_Gene_Name_Gencode_v36_df$gene_id <- gsub("\\..*", "", Gene_ID_Gene_Name_Gencode_v36_df$gene_id)

# Create a named vector for gene_id to gene_name mapping
gene_id_to_gene_name <- setNames(Gene_ID_Gene_Name_Gencode_v36_df$gene_name, Gene_ID_Gene_Name_Gencode_v36_df$gene_id)

# Map the gene_ids in transposed_rnaCounts_voom_Normalized_df_modified to gene names
colnames(transposed_rnaCounts_voom_Normalized_df_modified)[-1] <- gene_id_to_gene_name[colnames(transposed_rnaCounts_voom_Normalized_df_modified)[-1]]

# Check the first five rows and columns of the modified data.frame
transposed_rnaCounts_voom_Normalized_df_modified[1:5, 1:5]

# Find the common gene symbols
common_gene_symbols <- intersect(colnames(transposed_rnaCounts_voom_Normalized_df_modified), gene_symbols_estimate_common_genes_9883_df)

#Filter to only keep the gene_symbols_estimate_common_genes_9883_df
filtered_rnaCounts_voom_Normalized_df <- transposed_rnaCounts_voom_Normalized_df_modified[, c("submitter_id", common_gene_symbols)]


# Create a list of Ensembl gene IDs
ensembl_ids <- colnames(transposed_rnaCounts_voom_Normalized_df_modified)[-1]

# Connect to the Ensembl database
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Convert Ensembl gene IDs to gene symbols
gene_symbols <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                      filters="ensembl_gene_id",
                      values=ensembl_ids,
                      mart=ensembl)

# Create a named vector to map Ensembl gene IDs to gene symbols
id_to_symbol <- setNames(gene_symbols$external_gene_name, gene_symbols$ensembl_gene_id)

# Replace Ensembl gene IDs with gene symbols in the data.frame
colnames(transposed_rnaCounts_voom_Normalized_df_modified)[-1] <- id_to_symbol[colnames(transposed_rnaCounts_voom_Normalized_df_modified)[-1]]


# Filtering out rows with missing values in time and status
filtered_merged_data <- merged_data[!is.na(merged_data$time) & !is.na(merged_data$event), ]

# Extracting the new time, status, and x variables
filtered_time <- filtered_merged_data$time
filtered_status <- filtered_merged_data$event
filtered_x <- log2_gene_expression_data[rownames(log2_gene_expression_data) %in% rownames(filtered_merged_data), ]

#Load workspace file 
load("/cloud/home/r1816512/ESTIMATE_ANALYSIS/TCGA-STAD.RData.gz")

#Install R package rbsurv https://bioconductor.org/packages/release/bioc/html/rbsurv.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rbsurv")

library(rbsurv)

#Run rbsurv
rbsurv(filtered_time, filtered_status, transposed_filtered_x, gene.ID=rownames(transposed_filtered_x), method="efron", n.iter=10, n.fold=3)

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# convert Ensembl gene IDs to gene symbols
ensembl_ids <- results$gene.list
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
               filters = "ensembl_gene_id", 
               values = ensembl_ids, 
               mart = ensembl)

# extract gene symbols from the result
gene_symbols <- genes$external_gene_name