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

#scan for TFBSs in the promoter regions of Stromal_DEGs
#and prognostic gene signature (29 genes).

#extract the gene symbols from the Stromal_DEGs data.frame
stromal_DEG_symbols <- rownames(Stromal_DEGs)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#map the gene symbols to Ensembl gene IDs
Stromal_gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "hgnc_symbol",
                  values = stromal_DEG_symbols,
                  mart = ensembl)

#retrieve the promoter sequences for your Stromal_DEGs
#Retrieve the upstream and downstream sequences separately

stromal_upstream_sequences <- getSequence(id = gene_ids$ensembl_gene_id, type = "ensembl_gene_id", seqType = "coding_gene_flank", upstream = 2000, mart = ensembl)
stromal_downstream_sequences <- getSequence(id = gene_ids$ensembl_gene_id, type = "ensembl_gene_id", seqType = "coding_gene_flank", downstream = 2000, mart = ensembl)


# make sure the gene IDs are sorted in the same order for both data frames

stromal_upstream_sequences <- stromal_upstream_sequences[order(stromal_upstream_sequences$ensembl_gene_id), ]
stromal_downstream_sequences <- stromal_downstream_sequences[order(stromal_downstream_sequences$ensembl_gene_id), ]

#Combine the upstream and downstream sequences to create the promoter_sequences data frame 

stromal_promoter_sequences <- data.frame(ensembl_gene_id = stromal_upstream_sequences$ensembl_gene_id, sequence = paste(stromal_upstream_sequences$coding_gene_flank, stromal_downstream_sequences$coding_gene_flank, sep = ""))

stromal_promoter_sequences_set <- DNAStringSet(stromal_promoter_sequences$sequence)


PWMMatrixList <- lapply(PFMatrixList, toPWM)

stromal_predicted_TFBS_list <- lapply(PWMMatrixList, function(pwm) {
  lapply(stromal_promoter_sequences_set, function(seq) {
    TFBSTools::searchSeq(pwm, seq, min.score = "90%")
  })
})


# Extract the transcription factors from predicted_TFBS_list
stromal_predicted_TFs <- sapply(stromal_predicted_TFBS_list, function(x) x[[1]]@pattern@name)

# Print the predicted transcription factors
print(stromal_predicted_TFs)


#get attributes of interest from biomart
biomaRt_attributes_of_interest <- c("ensembl_gene_id", "external_gene_name",  "chromosome_name", "start_position", "end_position", "transcription_start_site", "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end")

stromal_gene_info_attributes <- getBM(attributes = biomaRt_attributes_of_interest,
                                      filters = "ensembl_gene_id",
                                      values = stromal_gene_ids$ensembl_gene_id,
                                      mart = ensembl)
#attributes from page sequences
sequences_attributes <- listAttributes(ensembl, page = "sequences")
sequences_attributes <- sequences_attributes[sequences_attributes$name %in% biomaRt_attributes_of_interest, "name"]

