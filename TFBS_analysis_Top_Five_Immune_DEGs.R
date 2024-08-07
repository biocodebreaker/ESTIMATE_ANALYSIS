#Transcription factor binding site (TFBS) analysis using JASPAR2022 and TFBSTools

#install TFBSTools

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TFBSTools")


# install JASPAR2022
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("JASPAR2022")

library(Biostrings)
library(TFBSTools)
library(JASPAR2022)

#Load the biomaRt library
library(biomaRt)

#Specify the species and type of profiles you want to retrieve. 
#In this case, we are interested in human profiles ("species" = 9606) and SELEX type
opts <- list()
opts[["species"]] <- 9606
opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE


#Get the JASPAR matrix set based on the specified options

PFMatrixList <- getMatrixSet(JASPAR2022, opts)

#scan for TFBSs in the promoter regions of Immune_DEGs
# Get the top_five Immune_DEGs 
top_five_Immune_DEGs <- Immune_DEGs[order(Immune_DEGs$adj.P.Val),][1:5, ]


#extract the gene symbols from the Immune_DEGs data.frame
top_five_Immune_DEGs_symbols <- rownames(top_five_Immune_DEGs)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#map the gene symbols to Ensembl gene IDs

top5_Immune_DEGs_ensembl_gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                       filters = "hgnc_symbol",
                                       values = top_five_Immune_DEGs_symbols,
                                       mart = ensembl)


#retrieve the promoter sequences for your Immune_DEGs

#Retrieve the upstream and downstream sequences separately

top5_Immune_DEGs_upstream_sequences <- getSequence(id = top5_Immune_DEGs_ensembl_gene_ids$ensembl_gene_id, type = "ensembl_gene_id", seqType = "gene_flank", upstream = 2000, mart = ensembl)

top5_Immune_DEGs_downstream_sequences <- getSequence(id = top5_Immune_DEGs_ensembl_gene_ids$ensembl_gene_id, type = "ensembl_gene_id", seqType = "gene_flank", downstream = 2000, mart = ensembl)

# make sure the gene IDs are sorted in the same order for both data frames

top5_Immune_DEGs_upstream_sequences <- top5_Immune_DEGs_upstream_sequences[order(top5_Immune_DEGs_upstream_sequences$ensembl_gene_id), ]
top5_Immune_DEGs_downstream_sequences <- top5_Immune_DEGs_downstream_sequences[order(top5_Immune_DEGs_downstream_sequences$ensembl_gene_id), ]

#Combine the upstream and downstream sequences


combined_promoter_sequences_upstream_downstream_top5_Immune_DEGs <- data.frame(ensembl_gene_id = top5_Immune_DEGs_upstream_sequences$ensembl_gene_id, sequence = paste(top5_Immune_DEGs_upstream_sequences$gene_flank, top5_Immune_DEGs_downstream_sequences$gene_flank, sep = ""))

top5_Immune_DEGs_promoter_sequences_set <- DNAStringSet(combined_promoter_sequences_upstream_downstream_top5_Immune_DEGs$sequence)

#find all functions in TFBSTools

#tfbs_tools_functions <- ls(getNamespace("TFBSTools"))


PWMMatrixList <- lapply(PFMatrixList, toPWM)

predicted_TFBS_list_top5_Immune_DEGs <- lapply(PWMMatrixList, function(pwm) {
  lapply(top5_Immune_DEGs_promoter_sequences_set, function(seq) {
    TFBSTools::searchSeq(pwm, seq, min.score = "90%")
  })
})


# Extract the transcription factors from predicted_TFBS_list

predicted_TFs_top5_Immune_DEGs <- sapply(predicted_TFBS_list_top5_Immune_DEGs, function(x) x[[1]]@pattern@name)

# Print the predicted transcription factors
print(predicted_TFs_top5_Immune_DEGs)


#Print each PWM alone 

predicted_TFBS_list_top5_Immune_DEGs$MA0002.1

#TFBS analysis done for each ensembl gene ID separately 

#Map the gene symbols to Ensembl gene IDs

LPXN_Immune_DEGs_ensembl_gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = "LPXN", mart = ensembl)

# Fetch the upstream sequence for the gene LPXN
LPXN_Immune_DEG_upstream_sequences <- getSequence(id = "ENSG00000110031", type = "ensembl_gene_id", seqType = "gene_flank", upstream = 2000, mart = ensembl)

LPXN_Immune_DEG_downstream_sequences <- getSequence(id = "ENSG00000110031", type = "ensembl_gene_id", seqType = "gene_flank", downstream = 2000, mart = ensembl)

#Combine upstream and downstream promoter sequences for LPXN
LPXN_Immune_DEG_promoter_sequences <- data.frame(ensembl_gene_id = "ENSG00000110031", sequence = paste(LPXN_Immune_DEG_upstream_sequences$gene_flank, LPXN_Immune_DEG_downstream_sequences$gene_flank, sep = ""))

#Make a DNAStringSet object for LPXN
LPXN_Immune_DEG_promoter_sequences_set <- DNAStringSet(LPXN_Immune_DEG_promoter_sequences$sequence)

#Predict TFs for LPXN using promoter sequences set

predicted_TFBS_list_LPXN_Immune_DEG <- lapply(PWMMatrixList, function(pwm) {
  lapply(LPXN_Immune_DEG_promoter_sequences_set, function(seq) {
    TFBSTools::searchSeq(pwm, seq, min.score = "90%")
  })
})

#List all the PWMs for the predicted TFs
predicted_TFs_LPXN_Immune_DEG <- sapply(predicted_TFBS_list_LPXN_Immune_DEG, function(x) x[[1]]@pattern@name)

#Print the PWM 

predicted_TFBS_list_LPXN_Immune_DEG$MA0002.1


top5_Immune_DEGs_info_attributes <- getBM(attributes = desired_attributes,
                                      filters = "ensembl_gene_id",
                                      values = top5_Immune_DEGs_ensembl_gene_ids$ensembl_gene_id,
                                      mart = ensembl)
