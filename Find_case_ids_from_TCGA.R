# Use the GDC Data Transfer Tool Client and Manifest file to download the .tsv files

#  wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip

# ./gdc-client download -m gdc_manifest.2023-12-11.txt

# Download case_ids for the TCGA-CESC project 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicDataCommons")

#Load the GenomicDataCommons library which has the results_all() function
library(GenomicDataCommons)

#Load dplyr which has the %>% operator
library(dplyr)


# Find cases for the TCGA-CESC project
cases_CESC <- cases() %>%
  filter(project.project_id == "TCGA-CESC") %>%
  results_all()

# Extract case_ids

case_ids_CESC <- cases_CESC$id

length (case_ids_CESC)
# [1] 307

# Download clinical data using case_ids

TCGA_CESC_clinical_data <- gdc_clinical(case_ids_CESC)

# Load required libraries
library(dplyr)

# Extract the relevant data.frames from clinical_data

TCGA_CESC_clinical_main_data <- TCGA_CESC_clinical_data$main
TCGA_CESC_demographic_data <- TCGA_CESC_clinical_data$demographic
TCGA_CESC_diagnoses_data <- TCGA_CESC_clinical_data$diagnoses

# Remove suffixes from submitter_id values in demographic and diagnoses data.frames

TCGA_CESC_demographic_data$submitter_id <- gsub("_demographic", "", TCGA_CESC_demographic_data$submitter_id)

TCGA_CESC_diagnoses_data$submitter_id <- gsub("_diagnosis", "", TCGA_CESC_diagnoses_data$submitter_id)


# Merge the dataframes based on submitter_id

TCGA_CESC_merged_data <- TCGA_CESC_clinical_main_data %>%
  left_join(TCGA_CESC_demographic_data, by = "submitter_id") %>%
  left_join(TCGA_CESC_diagnoses_data, by = "submitter_id")


# Extract submitter_id, days_to_death, days_to_last_follow_up, and vital_status

TCGA_CESC_survival_info <- TCGA_CESC_merged_data %>%
  dplyr::select(submitter_id, days_to_death, days_to_last_follow_up, vital_status)


#Calculate the survival time using the days_to_death and days_to_last_follow_up columns.
#If a patient is dead, use the days_to_death as the survival_time; 
#otherwise, if the patient is alive, use the days_to_last_follow_up.

TCGA_CESC_survival_time <- TCGA_CESC_survival_info %>%
  mutate(time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up))


#Set the event status based on vital_status, using 1 for "Dead" and 0 for "Alive"

TCGA_CESC_survival_time_censored <- TCGA_CESC_survival_time %>%
  mutate(event = ifelse(vital_status == "Dead", 1, 0))

# Extract the GeneSymbols from GeneCards_MSigDB_Oxidative_Stress_genes_9570
selected_genesymbols <- GeneCards_MSigDB_Oxidative_Stress_genes_9570$GeneSymbol

# Subset CIBERSORT_INPUT_TCGA_CESC_304 based on selected GeneSymbols
TCGA_CESC_304_GeneCards_MSigDB_Oxidative_Stress_genes_8712 <- CIBERSORT_INPUT_TCGA_CESC_304[CIBERSORT_INPUT_TCGA_CESC_304$GeneSymbol %in% selected_genesymbols, ]

dim(TCGA_CESC_304_GeneCards_MSigDB_Oxidative_Stress_genes_8712)
#[1] 8712  305

#GDC TCGA-STAD SUMMARY STATS

gender_counts <-metaMatrix.RNA %>% group_by(gender) %>% tally()

sample_type_counts <- metaMatrix.RNA %>% group_by(sample_type) %>% tally()

vital_status_counts <- metaMatrix.RNA %>% group_by(vital_status) %>% tally()

gender_counts$category <- c("gender","gender")

colnames(gender_counts)[1] <- "value"

sample_type_counts$category <- c("sample_type","sample_type")

colnames(sample_type_counts)[1] <- "value"

vital_status_counts$category <- c("vital_status","vital_status", "vital_status")

colnames(vital_status_counts)[1] <- "value"

combined_counts <-rbind(gender_counts, sample_type_counts, vital_status_counts)

combined_counts <- combined_counts[,c("category","value", "n")]

combined_counts

#Convert to markdown table 

knitr::kable(combined_counts)

summary_stats <- metaMatrix.RNA[,c("age_at_diagnosis", "days_to_death", "days_to_last_follow_up")] %>% summary()


summary_stats_df <- as.data.frame(summary_stats)
knitr::kable(summary_stats_df)



# get vector of file names (full names of tsv files)


files <- list.files("/cloud/home/r1816512/CERVICAL_CANCER/TCGA-CESC", pattern="\\.tsv$", full.names = TRUE, recursive = T)

# Read the tab-delimited file into a data.frame
TCGA_CESC_Sample_Sheet_309 <- read.delim("TCGA_CESC_Sample_Sheet_309.txt", header = TRUE, stringsAsFactors = FALSE)

# TCGA_CESC_gdc_manifest <- read.delim("gdc_manifest.2023-12-11.txt", header = TRUE, stringsAsFactors = FALSE)

# create vector of target columns
# targ_cols <- c("gene_id", "gene_name", "gene_type","fpkm_unstranded")

# targ_cols <- c("gene_id", "gene_name", "gene_type","tpm_unstranded")

targ_cols <- c("gene_id", "gene_name", "gene_type","unstranded")


# Load the data.table package
library(data.table)

# read each file, selecting the columns of interest
result <- lapply(files, \(f) fread(f,header=F, skip=6, select = c(1,2,3,4),col.names = targ_cols))

# set the names of this list
result <- setNames(result, basename(files))


# rowbind these all into one big data table
result <- rbindlist(result, idcol  = "filename")

# merge the meta matrix with result, and drop the filename column
result <- result[TCGA_CESC_Sample_Sheet_309, on="filename"][,filename:=NULL]

# swing the result wide, using dcast
result <- dcast(result, gene_id+gene_name+gene_type~Sample_ID, value.var = "unstranded")

# TCGA_CESC_RNASEQ_60660_FPKM_Unstranded_309 <- result

# TCGA_CESC_RNASEQ_60660_TPM_Unstranded_309 <- result

TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309 <- result




#Download case_ids for the TCGA-STAD project 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicDataCommons")


#Load dplyr which has the %>% operator
library(dplyr)

#Load the GenomicDataCommons library which has the results_all() function
library(GenomicDataCommons)

# Find cases for the TCGA-STAD project
cases_stad <- cases() %>%
  filter(project.project_id == "TCGA-STAD") %>%
  results_all()


# Find cases for the TCGA-CESC project
cases_CESC <- cases() %>%
  filter(project.project_id == "TCGA-CESC") %>%
  results_all()

# Extract case_ids
case_ids <- cases_stad$id

case_ids_CESC <- cases_CESC$id

length (case_ids_CESC)
# [1] 307

# Check if the number of cases is 448
cat("Number of cases in TCGA-STAD project:", length(case_ids), "\n")

# Download clinical data using case_ids
clinical_data <- gdc_clinical(case_ids)

TCGA_CESC_clinical_data <- gdc_clinical(case_ids_CESC)

# Check the structure of the clinical_data object which is a GDCClinicalList
str(clinical_data)

# Load required libraries
library(dplyr)

# Extract the relevant data.frames from clinical_data
main_data <- clinical_data$main
demographic_data <- clinical_data$demographic
diagnoses_data <- clinical_data$diagnoses

TCGA_CESC_clinical_main_data <- TCGA_CESC_clinical_data$main

TCGA_CESC_demographic_data <- TCGA_CESC_clinical_data$demographic
TCGA_CESC_diagnoses_data <- TCGA_CESC_clinical_data$diagnoses


# Remove suffixes from submitter_id values in demographic and diagnoses data.frames
demographic_data$submitter_id <- gsub("_demographic", "", demographic_data$submitter_id)
diagnoses_data$submitter_id <- gsub("_diagnosis", "", diagnoses_data$submitter_id)

TCGA_CESC_demographic_data$submitter_id <- gsub("_demographic", "", TCGA_CESC_demographic_data$submitter_id)

TCGA_CESC_diagnoses_data$submitter_id <- gsub("_diagnosis", "", TCGA_CESC_diagnoses_data$submitter_id)


# Merge the dataframes based on submitter_id
merged_data <- main_data %>%
  left_join(demographic_data, by = "submitter_id") %>%
  left_join(diagnoses_data, by = "submitter_id")


TCGA_CESC_merged_data <- TCGA_CESC_clinical_main_data %>%
  left_join(TCGA_CESC_demographic_data, by = "submitter_id") %>%
  left_join(TCGA_CESC_diagnoses_data, by = "submitter_id")


# Extract submitter_id, days_to_death, days_to_last_follow_up, and vital_status
survival_info <- merged_data %>%
  dplyr::select(submitter_id, days_to_death, days_to_last_follow_up, vital_status)


TCGA_CESC_survival_info <- TCGA_CESC_merged_data %>%
  dplyr::select(submitter_id, days_to_death, days_to_last_follow_up, vital_status)


# Install package GDCRNATools


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GDCRNATools")

library(GDCRNATools)

# Differential expression analysis and Survival Analysis

# Fetch metadata for TCGA-CESC RNAseq data using GDC API
TCGA_CESC_metaMatrix <- gdcParseMetadata(project.id='TCGA-CESC', data.type='RNAseq')

#Rename metaMatrix to dataframe 
TCGA_CESC_metadata_file_df <- data.frame(case_id=rownames(TCGA_CESC_metaMatrix), TCGA_CESC_metaMatrix, row.names=NULL)

#Merge RNAseq by raw counts (unstranded)
TCGA_CESC_rnaCounts <- gdcRNAMerge(metadata=TCGA_CESC_metaMatrix, path='/cloud/home/r1816512/CERVICAL_CANCER/TCGA-CESC', data.type='RNAseq')

TCGA_CESC_rnaCounts_df <- as.data.frame(TCGA_CESC_rnaCounts)

## Normalize the raw rnaCounts using voom transformation 

TCGA_CESC_rnaCounts_voom_Normalized <- gdcVoomNormalization(counts = TCGA_CESC_rnaCounts_df, filter = FALSE)

TCGA_CESC_rnaCounts_voom_Normalized_df <- as.data.frame(TCGA_CESC_rnaCounts_voom_Normalized)

#transpose data.frame
transposed_rnaCounts_voom_Normalized_df <- t(rnaCounts_voom_Normalized_df)

#Add submitter_id to header (column 1)
transposed_rnaCounts_voom_Normalized_df <- data.frame(submitter_id = rownames(transposed_rnaCounts_voom_Normalized_df), transposed_rnaCounts_voom_Normalized_df)
rownames(transposed_rnaCounts_voom_Normalized_df) <- NULL

#Remove the "-01" suffix to make submitter_ids consistent 
transposed_rnaCounts_voom_Normalized_df_modified <- transposed_rnaCounts_voom_Normalized_df

transposed_rnaCounts_voom_Normalized_df_modified$submitter_id <- gsub("-01$", "", transposed_rnaCounts_voom_Normalized_df$submitter_id)

#not all the gene symbols in gene_symbols_estimate_common_genes_9883_df 
#are present in the transposed_rnaCounts_voom_Normalized_df_modified data frame.
#use the intersect function to find the common gene symbols between the two data frames

common_gene_symbols <- intersect(colnames(transposed_rnaCounts_voom_Normalized_df_modified), gene_symbols_estimate_common_genes_9883_df)

#Merge the data.frames
merged_survival_info_estimate_scores_gene_expression_limma <- merge(merged_survival_info_censored_estimate_scores_categorized, filtered_rnaCounts_voom_Normalized_df, by = "submitter_id")


#gene symbols in gene_symbols_estimate_common_genes_9883_df that were not found in transposed_rnaCounts_voom_Normalized_df_modified
missing_genes <- gene_symbols_estimate_common_genes_9883_df[!gene_symbols_estimate_common_genes_9883_df %in% colnames(transposed_rnaCounts_voom_Normalized_df_modified)]
missing_genes


## perform differential gene expression analysis using limma
DEG_analysis_Limma <- gdcDEAnalysis(counts = rnaCounts, 
                                    group      = metaMatrix.RNA$sample_type, 
                                    comparison = 'PrimaryTumor-SolidTissueNormal', 
                                    method     = 'limma',
                                    filter=TRUE)

## Generate Volcano plot and BarPlot
gdcVolcanoPlot(DEG_analysis_Limma)

gdcBarPlot(deg = DEG_analysis_Limma, angle = 45, data.type = 'RNAseq')

## filter out all the significant differentially expressed genes based on gene type.

DE_protein_coding <- gdcDEReport(deg = DEG_analysis_Limma, gene.type = 'protein_coding')

DE_ALL <- gdcDEReport(deg = DEG_analysis_Limma, gene.type = 'all')


##A total of 2288 genes were found differentially expressed between PrimaryTumor and SolidTissueNormal with statistical significance.

dim(DE_ALL)

##A total of 2031 protein coding genes were found differentially expressed between PrimaryTumor and SolidTissueNormal with statistical significance.
dim(DE_protein_coding) 


### Functional Enrichment Analysis 
Functional_EnrichOutput <- gdcEnrichAnalysis(gene = rownames(DEG_analysis_Limma), simplify = TRUE)

## Kaplan Meier survival analysis 

KM_Survival_Analysis_Output <- gdcSurvivalAnalysis(gene = rownames(DEG_analysis_Limma), 
                                                   method   = 'KM', 
                                                   rna.expr = rnaExpr, 
                                                   metadata = metaMatrix, 
                                                   sep      = 'median')



# Sort the output of univariate survival analysis based on pvalue
KM_Survival_Analysis_Output$pValue <- as.numeric(KM_Survival_Analysis_Output$pValue)

sorted_survOutput<- KM_Survival_Analysis_Output[order(KM_Survival_Analysis_Output$pValue),]

sorted_survOutput[1:5,]




#Calculate the survival time using the days_to_death and days_to_last_follow_up columns.
#If a patient is dead, use the days_to_death as the survival_time; 
#otherwise, if the patient is alive, use the days_to_last_follow_up.
survival_info_calculated <- survival_info %>%
  mutate(time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up))

TCGA_CESC_survival_time <- TCGA_CESC_survival_info %>%
  mutate(time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up))


#Set the event status based on vital_status, using 1 for "Dead" and 0 for "Alive"

TCGA_CESC_survival_time_censored <- TCGA_CESC_survival_time %>%
  mutate(event = ifelse(vital_status == "Dead", 1, 0))

# Remove unnecessary columns

#survival_info_time_event <- survival_info_censored %>%
 # select(submitter_id, time, event)

# Print survival information
#print(survival_info_time_event)

#Replace periods with hyphens in submitter_id values so that they are consistent

transposed_estimate_scores_448_df$submitter_id <- gsub("\\.01A", "", transposed_estimate_scores_448_df$submitter_id)
transposed_estimate_scores_448_df$submitter_id <- gsub("\\.", "-", transposed_estimate_scores_448_df$submitter_id)

library(dplyr)
library("survival")
library(survminer)
# survival_info_time_event contains survival information for 443 patients
#transposed_estimate_scores_448_df contains estimate scores for 448 patients. 
#only patients (submitter_ids) with both survival information and estimate score data available will be included.
#merge the data.frames using inner_join so that only matching submitter_ids in both data.frames will be kept.

library(dplyr)

# Read in TCGA survival info
survival_info_time_event <- read.csv("survival_info_time_event.csv", header = TRUE)

dim(survival_info_time_event)
#[1] 443   3

# Prepare the data: filter out rows with NA or non-positive `time` values
survival_info_time_event_filtered <- survival_info_time_event %>%
  filter(!is.na(time) & time > 0)

dim(survival_info_time_event_filtered)
dim(survival_info_time_event_filtered)
#[1] 412   3

# read in estimate scores
transposed_estimate_scores_448_df <- read.delim("transposed_estimate_scores_448.txt", header = TRUE)

dim(transposed_estimate_scores_448_df)
#[1] 448   4

# Replace periods with hyphens in `transposed_estimate_scores_448_df`
transposed_estimate_scores_448_df$submitter_id <- gsub("\\.", "-", transposed_estimate_scores_448_df$submitter_id)

# Remove the last part after the last hyphen (e.g., '-01A') using base R's sub function
transposed_estimate_scores_448_df$submitter_id <- sub("-[^-]+$", "", transposed_estimate_scores_448_df$submitter_id)

# Merge the data frames by `submitter_id`
merged_survival_info_time_event_estimate_scores_filtered <- survival_info_time_event_filtered %>%
  inner_join(transposed_estimate_scores_448_df, by = "submitter_id")

dim(merged_survival_info_time_event_estimate_scores_filtered)
dim(merged_survival_info_time_event_estimate_scores_filtered)
#[1] 417   6

# Merge the data frames by `submitter_id`
merged_survival_info_time_event_estimate_scores <- survival_info_time_event %>%
  inner_join(transposed_estimate_scores_448_df, by = "submitter_id")

dim(merged_survival_info_time_event_estimate_scores)
#[1] 448   6

# Merge the data frames by `submitter_id`
merged_survival_info_time_event_estimate_scores_noSTN <- survival_info_time_event %>%
  inner_join(transposed_estimate_scores_412_noSTN, by = "submitter_id")

# Prepare the data: filter out rows with NA or non-positive `time` values
merged_survival_info_time_event_estimate_scores_noSTN_filtered <- merged_survival_info_time_event_estimate_scores_noSTN %>%
  filter(!is.na(time) & time > 0)

dim(merged_survival_info_time_event_estimate_scores_noSTN_filtered)
#[1] 383   6

# Prepare the data: filter out rows with NA or non-positive `time` values
merged_survival_info_time_event_estimate_scores_filtered <- merged_survival_info_time_event_estimate_scores %>%
  filter(!is.na(time) & time > 0)


# Prepare the data: filter out rows with NA or non-positive `time` values
merged_survival_info_time_event_estimate_scores_filtered <- merged_survival_info_time_event_estimate_scores %>%
  filter(!is.na(time) & time > 0)

dim(merged_survival_info_time_event_estimate_scores_filtered)
#[1] 417   6


#determine stromal and immune score optimal cutoff using survminer package

stromal_immune_optimal_cutoff <- surv_cutpoint(merged_survival_info_time_event_estimate_scores, time = "time", event = "event",
                                               variables = c("StromalScore", "ImmuneScore"))


summary(stromal_immune_optimal_cutoff)

plot(stromal_immune_optimal_cutoff, "StromalScore", palette = "npg")

plot(stromal_immune_optimal_cutoff, "ImmuneScore", palette = "npg")

#Categorize stromal and immune scores into low and high categories 
stromal_immune_optimal_cutoff_categorized <- surv_categorize(stromal_immune_optimal_cutoff)


# Fit survival curves and visualize
library("survival")
library(survminer)
fit <- survfit(Surv(time, event) ~StromalScore, data = stromal_immune_optimal_cutoff_categorized )
ggsurvplot(fit, data = stromal_immune_optimal_cutoff_categorized, risk.table = TRUE, pval = TRUE, dpi = 300)


#Perform survival analysis using survival library and maxstat package. 

library(survival)
library(maxstat)

# Create the Surv object using with()
surv_obj <- with(merged_survival_info_censored_estimate_scores, Surv(time, event))

# Find the optimal cutoff values
maxstat_stromal <- maxstat.test(surv_obj ~ merged_survival_info_censored_estimate_scores$StromalScore, data = merged_survival_info_censored_estimate_scores, smethod = "LogRank")
optimal_cutoff_stromal <- maxstat_stromal$statistic[1]

maxstat_immune <- maxstat.test(surv_obj ~ merged_survival_info_censored_estimate_scores$ImmuneScore, data = merged_survival_info_censored_estimate_scores, smethod = "LogRank")
optimal_cutoff_immune <- maxstat_immune$statistic[1]

# Add StromalScore_group and ImmuneScore_group columns
merged_survival_info_censored_estimate_scores_categorized <- merged_survival_info_time_event_estimate_scores %>%
  mutate(StromalScore_group = ifelse(StromalScore > 1323.487, "high", "low"),
         ImmuneScore_group = ifelse(ImmuneScore > 1957.582, "high", "low"))

# Create the Surv object
surv_obj <- with(merged_survival_info_censored_estimate_scores_categorized, Surv(time, event))

# Fit Kaplan-Meier survival curves for the StromalScore group
km_fit_stromal <- survfit(surv_obj ~ StromalScore_group, data = merged_survival_info_censored_estimate_scores_categorized)
ggsurvplot(km_fit_stromal, data = merged_survival_info_censored_estimate_scores_categorized, pval = TRUE, risk.table = TRUE)

# Fit Kaplan-Meier survival curves for the Immune group
km_fit_immune <- survfit(surv_obj ~ Immune_group, data = merged_survival_info_censored_estimate_scores_categorized)

# Plot Kaplan-Meier survival curves for the Immune group
ggsurvplot(km_fit_immune, data = merged_survival_info_censored_estimate_scores_categorized, pval = TRUE, risk.table = TRUE)



# Make submitter_ids in both data.frames consistent 
rownames(transposed_estimate_common_genes_9883) <- gsub("\\.", "-", rownames(transposed_estimate_common_genes_9883))
rownames(transposed_estimate_common_genes_9883) <- substr(rownames(transposed_estimate_common_genes_9883), 1, nchar(rownames(transposed_estimate_common_genes_9883)) - 4)

# Filter the data.frame to keep only rows with submitter_id ending in -01A
#Keep gene expression data only for "Primary Tumor" that have the "-01A" suffix

transposed_estimate_common_genes_9883_primary_tumor <- transposed_estimate_common_genes_9883[grep("-01A$", rownames(transposed_estimate_common_genes_9883)),]

#perform an inner join of the two data.frames, keeping only the rows that have a matching submitter_id in both data.frames.
merged_surv_info_estimate_scores_gene_expression_primary_tumor_9883_genes <- merge(merged_survival_info_censored_estimate_scores_categorized, transposed_estimate_common_genes_9883_primary_tumor, by = "submitter_id", all = FALSE)




# Write the result data frame to a tab-delimited text file
# write.table(Merged_GDC_RNASEQ_60660_FPKM_Unstranded, "Merged_GDC_RNASEQ_60660_FPKM_Unstranded.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# write.table(TCGA_CESC_RNASEQ_60660_FPKM_Unstranded_309, "TCGA_CESC_RNASEQ_60660_FPKM_Unstranded_309.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# write.table(TCGA_CESC_RNASEQ_60660_TPM_Unstranded_309, "TCGA_CESC_RNASEQ_60660_TPM_Unstranded_309.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


# GTEx Analysis

# Load the data.table package
library(data.table)

# Reading the first 24 rows from the GTEx gene read counts file

# GTEx_v8_RNASeQ_gene_tpm_24 <- read.delim(file="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip=2, nrows = 24)

GTEx_v8_RNASeQ_gene_read_counts_24 <- read.delim(file="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2, nrows = 24)


# Sample identifiers to keep

GTEx_SAMPIDs <- c(
  "GTEX-N7MT-0826-SM-EV794",
  "GTEX-OHPN-2226-SM-E9TI7",
  "GTEX-P78B-2326-SM-EZ6KO",
  "GTEX-P78B-2426-SM-EWRM2",
  "GTEX-PLZ4-2226-SM-EZ6KS",
  "GTEX-Q2AG-2326-SM-EZ6KY",
  "GTEX-S32W-1526-SM-4AD6Z",
  "GTEX-S32W-1626-SM-4AD6G",
  "GTEX-S341-1126-SM-4AD6T",
  "GTEX-S341-1326-SM-4AD72",
  "GTEX-S4UY-1426-SM-4AD6Y",
  "GTEX-T5JW-0726-SM-4DM6D",
  "GTEX-T5JW-0826-SM-EYYVD",
  "GTEX-T6MO-1426-SM-4DM73",
  "GTEX-TML8-0726-SM-4DXTT",
  "GTEX-TSE9-2726-SM-4DXSQ",
  "GTEX-TSE9-2826-SM-4DXTF",
  "GTEX-U3ZN-1626-SM-4DXTZ",
  "GTEX-ZPIC-1326-SM-DO91Y"
)


# Read only the GTEx_SAMPIDs columns from the file into the data.frame

# GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs <- fread(file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2, select = c("Name", "Description", GTEx_SAMPIDs))

GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs <- fread(file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip = 2, select = c("Name", "Description", GTEx_SAMPIDs))

dim(GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs)
#[1] 56200    21

class(GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs)
#[1] "data.table" "data.frame"

# Drop the "Name" column
# GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs[, Name := NULL]

GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs[, Name := NULL]


# Sort the data.frame by the "Description" column
# GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs <- GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs[order(Description)]

GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs <- GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs[order(Description)]


# Rename the "Description" column to "GeneSymbol"
# setnames(GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs, old = "Description", new = "GeneSymbol")

setnames(GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs, old = "Description", new = "GeneSymbol")

# Extract the "GeneSymbol" column from GeneCards_MSigDB_Oxidative_Stress_genes_9570

# Read the file into a data.frame

GeneCards_MSigDB_Oxidative_Stress_genes_9570 <- read.delim("GeneCards_MSigDB_Oxidative_Stress_genes_9570.txt", header = TRUE, stringsAsFactors = FALSE)


selected_genesymbols <- GeneCards_MSigDB_Oxidative_Stress_genes_9570$GeneSymbol

# Filter rows in GTEx_v8_RNASeQ_gene_tpm_CESC_SAMPIDs where "GeneSymbol" is in selected_genesymbols

GTEx_CESC_OS_gene_read_counts_8311 <- GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs[GTEx_v8_RNASeQ_gene_read_counts_CESC_SAMPIDs$GeneSymbol %in% selected_genesymbols, ]

dim(GTEx_CESC_OS_gene_read_counts_8311)
# [1] 8311   20

# Drop the "gene_id" and "gene_type" columns
TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309[, gene_id := NULL]
TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309[, gene_type := NULL]


# Rename the "gene_name" column to "GeneSymbol"
setnames(TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309, old = "gene_name", new = "GeneSymbol")

# Sort the data.frame by the "GeneSymbol" column

TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309 <- TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309[order(GeneSymbol)]

# Filter rows in TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309 where "GeneSymbol" is in selected_genesymbols

TCGA_CESC_OS_gene_read_counts_8712 <- TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309[TCGA_CESC_RNASEQ_60660_gene_read_counts_unstranded_309$GeneSymbol %in% selected_genesymbols, ]

dim(TCGA_CESC_OS_gene_read_counts_8712)
# [1] 8712  310

class(TCGA_CESC_OS_gene_read_counts_8712)
# [1] "data.table" "data.frame"


# Remove duplicates from TCGA data frame
TCGA_CESC_OS_gene_read_counts_dedup_8690 <- TCGA_CESC_OS_gene_read_counts_8712[!duplicated(TCGA_CESC_OS_gene_read_counts_8712$GeneSymbol), ]

dim(TCGA_CESC_OS_gene_read_counts_dedup_8690)
# [1] 8253  310

# Remove duplicates from GTEx data frame
GTEx_CESC_OS_gene_read_counts_dedup_8264 <- GTEx_CESC_OS_gene_read_counts_8311[!duplicated(GTEx_CESC_OS_gene_read_counts_8311$GeneSymbol), ]
dim(GTEx_CESC_OS_gene_read_counts_dedup_8264)
# [1] 8264  20

# Find common gene symbols

common_genes <- intersect(TCGA_CESC_OS_gene_read_counts_dedup_8690$GeneSymbol, GTEx_CESC_OS_gene_read_counts_dedup_8264$GeneSymbol)

length(common_genes)
length(common_genes)
#[1] 8253

# Subsetting TCGA dataset
TCGA_CESC_OS_gene_read_counts_dedup_8253 <- TCGA_CESC_OS_gene_read_counts_dedup_8690[TCGA_CESC_OS_gene_read_counts_dedup_8690$GeneSymbol %in% common_genes, ]

# Subsetting GTEx dataset
GTEx_CESC_OS_gene_read_counts_dedup_8253 <- GTEx_CESC_OS_gene_read_counts_dedup_8264[GTEx_CESC_OS_gene_read_counts_dedup_8264$GeneSymbol %in% common_genes, ]


# Specify the columns to transfer
columns_to_transfer_names <- c("TCGA-MY-A5BF-11A", "TCGA-FU-A3EO-11A", "TCGA-HM-A3JJ-11A")

# Create a data.frame with the specified columns
columns_to_transfer <- TCGA_CESC_OS_gene_read_counts_dedup_8253[, ..columns_to_transfer_names, with = FALSE]

# Check the dimensions of the columns_to_transfer data.frame
dim(columns_to_transfer)

# Drop the specified columns from TCGA_CESC_OS_gene_read_counts_dedup_8253
TCGA_CESC_OS_gene_read_counts_dedup_8253[, c("TCGA-MY-A5BF-11A", "TCGA-FU-A3EO-11A", "TCGA-HM-A3JJ-11A") := NULL]


TCGA_CESC_OS_gene_read_counts_dedup_8253_307 <- TCGA_CESC_OS_gene_read_counts_dedup_8253

# Check the dimensions of the updated data.frame
dim(TCGA_CESC_OS_gene_read_counts_dedup_8253_307)
# [1] 8253  307

# Merge columns_to_transfer with GTEx_CESC_OS_gene_read_counts_dedup_8253
GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3 <- cbind(GTEx_CESC_OS_gene_read_counts_dedup_8253, columns_to_transfer)

# Check the dimensions of the merged data
dim(GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3)
# [1] 8253   23

GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3_limma <- GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3

GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3_limma[, GeneSymbol := NULL]

TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal <- cbind(TCGA_CESC_OS_gene_read_counts_dedup_8253_307, GTEx_CESC_OS_gene_read_counts_dedup_8253_TCGA_3_limma)
dim(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal)
#[1] 8253  329

class(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal)
#[1] "data.table" "data.frame"

# Convert data.table to data.frame
TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal <- as.data.frame(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal)

#set rownames to GeneSymbol
rownames(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal) <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal$GeneSymbol
TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal$GeneSymbol <- NULL

# Ensure all columns are numeric
TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal[] <- lapply(
  TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal,
  as.numeric
)


# Remove rows/genes with very few read counts before runnning voom

keep <- rowSums(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal) > 20

TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_filtered <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal[keep,]

TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910 <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_filtered

dim(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910)
#[1] 7910  328

# Use limma R package, to identify differentially expressed genes between 
# TCGA_tpm_CESC_Oxidative_Stress_genes_8253_307 which contains 306 tumor samples and 
# GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_TCGA_3 which contains 22 normal samples.

library(limma)

# Create the design matrix
# sample_groups <- factor(rep(c("Tumor", "Normal"), each = 32))

sample_groups <- factor(rep(c("Tumor", "Normal"), c(306, 22)))


# Convert the sample_groups to a factor
sample_groups <- factor(sample_groups, levels = c("Tumor", "Normal"))

# Create a data frame for the design matrix
design_matrix <- model.matrix(~0 + sample_groups)

# Set appropriate column names in the design matrix
colnames(design_matrix) <- c("Tumor", "Normal")

# Normalize using voom 

v <- voom(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910, design_matrix)

# Fit the linear model
fit <- lmFit(v, design_matrix)


# Define the contrast matrix
contrast_matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design_matrix)

# Apply the contrast to the fit
fit_contrast <- contrasts.fit(fit, contrast_matrix)

# Perform empirical Bayes moderation of the standard errors
fit_eBayes <- eBayes(fit_contrast, robust = TRUE)

# Warning message: Zero sample variances detected, have been offset away from zero


# When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and for prioritizing genes that are biologically important. 
# treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods.

fit_treat <- treat(fit_contrast, fc = 1.5)

fit_treat_robust <- treat(fit_contrast, fc = 1.5, robust = TRUE) 
# Warning message: 89 very small variances detected, have been offset away from zero 


# Get the results of differential expression analysis
TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs <- topTable(fit_eBayes, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "fdr", sort.by = "P")

TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust <- topTreat(fit_treat_robust, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "fdr", sort.by = "P")

# write.table(TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs, file = "TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs.txt", sep = "\t", row.names = TRUE)

# write.table(TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust, file = "TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust.txt", sep = "\t", row.names = TRUE)write.table(TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust, file = "TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust.txt", sep = "\t", row.names = TRUE)

# write.table(TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust, file = "TCGA_GTEx_v8_tpm_CESC_OS_DEGs_8253_fit_treat_robust.txt", sep = "\t", row.names = TRUE)


write.table(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs, file = "TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910_DEGs.txt", sep = "\t", row.names = TRUE)
 

#Filter to include genes with both positive and negative logFC values

TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs_filtered <- TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs[
  TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$adj.P.Val < 0.001 & abs(TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$logFC) > 4
  | TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$logFC < -4,
]

#Filter to find underexpressed (logFC < 0) and significant genes (adj.P.Val < 0.05)

Bcell_spec_expr_TN_Pairs_LM10_199_DEGs_underexpressed <- Bcell_spec_expr_TN_Pairs_LM10_199_DEGs[
  Bcell_spec_expr_TN_Pairs_LM10_199_DEGs$logFC < 0 & Bcell_spec_expr_TN_Pairs_LM10_199_DEGs$adj.P.Val < 0.05,
]

# Filter underexpressed DEGs (negative logFC < -4 and adj.P.Val < 0.0001)
underexpressed_DEGs_logFC4_adj.P.Val_0.0001 <- TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs[
  TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$logFC < -4 & TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$adj.P.Val < 0.0001,
]



# Filter overexpressed DEGs (positive logFC > 4 and adj.P.Val < 0.0001)
overexpressed_DEGs_logFC4_adj.P.Val_0.0001 <- TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs[
  TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$logFC > 4 & TCGA_GTEx_v8_tpm_CESC_Oxidative_Stress_genes_8253_DEGs$adj.P.Val < 0.0001,
]

# Combine both sets of DEGs
combined_DEGs <- rbind(underexpressed_DEGs, overexpressed_DEGs)

# Extract gene names of the combined DEGs
combined_DEGs_names <- rownames(combined_DEGs)
# Extract gene expression values of the combined DEGs from the main data frame
combined_DEGs_expression <- Bcell_spec_expr_TN_Pairs_LM10_199[combined_DEGs_names ,]

combined_DEGs_scaled <- scale(combined_DEGs_expression)

# Create a heatmap with hierarchical clustering

# Produce Heatmap of the combined DEGs with gene names
pheatmap(combined_DEGs_scaled, 
         color = colorRampPalette(c("blue", "black", "yellow"))(100), 
         scale = "none", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         annotation_row = data.frame(Gene = rownames(combined_DEGs_scaled), row.names = rownames(combined_DEGs_scaled)))


# TCGA CESC SURVIVAL ANALYSIS 
# Combine expression (gene read counts) used in limma DEG analysis with surv info (time, event)

library(dplyr)

TCGA_CESC_survival_time_censored <- readRDS("TCGA_CESC_survival_time_censored.rds")

TCGA_CESC_survival_time_censored <- as.data.frame(TCGA_CESC_survival_time_censored)
 
# Add a column "GeneSymbol" containing original row names
 
TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910$GeneSymbol <- rownames(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910)
 
 # Set row names to row numbers
 
rownames(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910) <- NULL
 
# Reorder the columns so that "GeneSymbol" is the first column

TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910 <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910 %>%
  select(GeneSymbol, everything())
 
# keep only columns 1 to 307 so that GTEX and TCGA Normal are removed and TCGA Tumors (Primary and Metastatic) remain
TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors <- TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910[, 1:307]

dim(TCGA_CESC_OS_expression_gene_read_counts_limma_306_Tumors)
#[1] 7910  307

# Prepare the data
complete_data <- merged_survival_info_estimate_scores_gene_expression_limma[!is.na(merged_survival_info_estimate_scores_gene_expression_limma$time), ]
complete_data <- complete_data[complete_data$time > 0, ]

# Extract survival data
survival_data <- complete_data[, c("time", "event")]

# Create the Surv() object
surv_obj <- Surv(survival_data$time, survival_data$event)

# Initialize an empty list to store results
univariate_results <- list()

# Perform univariate Cox regression for each gene in intersection_genes
for (gene in intersection_genes) {
  cox_model <- coxph(surv_obj ~ get(gene), data = complete_data)
  result <- summary(cox_model)
  univariate_results[[gene]] <- result
}


# Display the summary results for the first 40 genes
for (gene in intersection_genes[60:90]) {
  cat("Results for", gene, "\n")
  print(univariate_results[[gene]])
}

# Define the significance level
significance_level <- 0.05

# Extract significant genes
TCGA_univariate_significant_genes <- names(TCGA_Significant_intersection_genes_univariate_results)[sapply(TCGA_Significant_intersection_genes_univariate_results, function(result) {
  p_value <- result$coefficients["get(gene)", "Pr(>|z|)"]
  p_value < significance_level
})]

# Display the list of significant genes
print(TCGA_univariate_significant_genes)

# Extract significant genes
GSE84426_GSE84433_univariate_significant_genes <- names(GSE84426_GSE84433_intersection_genes_univariate_results)[sapply(GSE84426_GSE84433_intersection_genes_univariate_results, function(result) {
  p_value <- result$coefficients["get(gene)", "Pr(>|z|)"]
  p_value < significance_level
})]

# Display the list of significant genes
print(GSE84426_GSE84433_univariate_significant_genes)

length(GSE84426_GSE84433_univariate_significant_genes)
# [1] 454

# Combine gene expressions into a formula using reformulate
gene_formula <- reformulate(TCGA_Significant_intersection_genes, response = "surv_obj")

# Fit the multivariate Cox regression model
multivariate_cox_model <- coxph(gene_formula, data = merged_survival_info_estimate_scores_gene_expression_limma)

# Summarize the results
multivariate_results <- summary(multivariate_cox_model)

# Print the results
print(multivariate_results)
