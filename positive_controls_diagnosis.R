

library(dplyr)
library(tidyr)



################################
### INPUT AND CLEAN ALLELE DATA
################################

INPUT_FOLDER <- "MULE_NextSEQ02_120625_RESULTS_v0.1.8_FILTERED/" 
data <- read.csv(paste0(INPUT_FOLDER, "/allele_data_global_max_0_filtered.csv"))


min_allele_read_count <- 10
maf_filter <- 0.02


clean_data <- function(data, maf_filter, min_allele_read_count) {

  # Identify sampleID with >=50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
  good_sampleID <- data %>%
    group_by(sampleID, locus) %>%
    summarize(reads = sum(reads)) %>%
    
    # Filter rows where reads >= 100
    filter(reads >= 100) %>%   
    
    # Filter NIDA where n_loci >= 50
    summarise(n_loci = n_distinct(locus)) %>%
    filter(n_loci >= 50) %>%
    
    pull(sampleID)
  
  # Keep good quality samples
  data <- data[data$sampleID %in% good_sampleID, ]
  
  # Keep samples with > 10000 reads: ANDRÉS
  read_counts_above_10K_reads <- data %>%
    group_by(sampleID) %>%
    summarise(total_reads = sum(reads)) 
  
  read_counts_above_10K_reads <- read_counts_above_10K_reads[read_counts_above_10K_reads$total_reads > 10000,]
  data <- data[data$sampleID %in% read_counts_above_10K_reads$sampleID, ]
  
  # Keep pool 1A
  data <- data[grepl("-1A$", data$locus), ]
  
  # Ignore masking, turn it into ref (.)
  data$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", data$pseudo_cigar) # Remove masking
  data$pseudo_cigar <- ifelse(data$pseudo_cigar == "" | is.na(data$pseudo_cigar), ".", data$pseudo_cigar) # If empty, add "." since it was reference
  
  # Aggregate unmasked sequences
  data_agg <- data %>%
    group_by(sampleID, locus, pseudo_cigar) %>%
    summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus))
  
  # Create allele column
  data_agg$allele <- paste0(data_agg$locus, "__", data_agg$pseudo_cigar)
  
  # Remove indels
  data_agg <- data_agg[!grepl("I=", data_agg$allele), ] # Remove alleles with I (insertion)
  data_agg <- data_agg[!grepl("D=", data_agg$allele), ] # Remove alleles with D (deletion)
  
  # Apply MAF filter
  data_agg <- data_agg[data_agg$norm.reads.locus > maf_filter, ] # From contaminants study
  
  # Remove alleles with low read counts
  data_agg <- data_agg[data_agg$reads > min_allele_read_count, ]
  
  # Remove pseudo_cigar column
  data_agg <- data_agg %>% select(-pseudo_cigar)
  
  # remove undetermined
  data_agg <- data_agg[!grepl("Undetermined", data_agg$sampleID, ignore.case = T),]
  
  return(data_agg)
}

# clean data
data_agg <- clean_data(data, maf_filter, min_allele_read_count) 




################################
### INPUT AND CLEAN TRUTH CONTROL DATA
################################

truth <- read.csv("pool1A_truth.csv")

# Ignore masking, turn it into ref (.)
truth$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", truth$pseudo_cigar) # Remove masking
truth$pseudo_cigar <- ifelse(truth$pseudo_cigar == "" | is.na(truth$pseudo_cigar), ".", truth$pseudo_cigar) # If empty, add "." since it was reference

# Create allele column
truth$allele <- paste0(truth$locus, "__", truth$pseudo_cigar)



################################
### EXTRACT CONTROLS SEPARATELY
################################

# controls
controls_3d7 <- data_agg[grepl("3d7", data_agg$sampleID, ignore.case = T),]
controls_dd2 <- data_agg[grepl("dd2", data_agg$sampleID, ignore.case = T),]

# extract field data
field <- data_agg[!grepl("3d7|dd2", data_agg$sampleID, ignore.case = TRUE), ]

# controls truth genotype
truth_3d7 <- truth[grepl("3d7", truth$Strain, ignore.case = T),]
truth_dd2 <- truth[grepl("dd2", truth$Strain, ignore.case = T),]



################################
### EXTRACT NON-REF ALLELES
################################

nonref_3d7 <- setdiff(controls_3d7$allele, truth_3d7$allele)
nonref_dd2 <- setdiff(controls_dd2$allele, truth_dd2$allele)

controls_3d7_nonref_alleles <- controls_3d7[controls_3d7$allele %in% nonref_3d7, ]
controls_dd2_nonref_alleles <- controls_dd2[controls_dd2$allele %in% nonref_dd2, ]



################################
### CREATE NONREF TABLES
################################

# checking within run or unknown. no cross contam here. would need a db to compare with beforehand
controls_3d7_nonref_alleles$source <- ifelse(controls_3d7_nonref_alleles$allele %in% field$allele, "field", "unknown")
controls_dd2_nonref_alleles$source <- ifelse(controls_dd2_nonref_alleles$allele %in% field$allele, "field", "unknown")

ALL_nonref_alleles <- do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0, 
                                           list(controls_3d7_nonref_alleles, controls_dd2_nonref_alleles)))



################################
### CREATE SUMMARY
################################

SUMMARY <- ALL_nonref_alleles %>%
  group_by(sampleID) %>%
  summarise(
    fixed_alleles = sum(norm.reads.locus == 1, na.rm = TRUE),
    source_field = sum(source == "field", na.rm = TRUE),
    source_unknown = sum(source == "unknown", na.rm = TRUE)
  )

# add the clean controls just in case
all_controls <- unique(c(controls_3d7$sampleID, controls_dd2$sampleID))

# Ensure all_controls are included and fill in missing rows with 0s
SUMMARY_COMPLETE <- tibble(sampleID = all_controls) %>%
  left_join(SUMMARY, by = "sampleID") %>%
  replace_na(list(
    fixed_alleles = 0,
    source_unknown = 0,
    source_field = 0
  )) 



################################
### OUTPUTS
################################

write.csv(SUMMARY_COMPLETE, "summary.csv", row.names = F)
write.csv(ALL_nonref_alleles, "all_nonref_alleles.csv", row.names = F)
