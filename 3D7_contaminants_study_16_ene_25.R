
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(vegan)
library(tibble)
library(reshape2)


############################################################
####### PARAMETERS AND INPUTS #######--------------------
############################################################

# provice folder name containing all *.FILTERED genomic data
main_dir <- "../amp_performance/ALL_SEQ_RESULTS_FILTERED"

# import genomic data
filtered_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

list_of_dfs <- list()

for (dir in filtered_dirs) {
  
  print(paste0("Importing: ", dir))
  
  folder_name <- gsub("_RESULTS.*", "", basename(dir))
  
  csv_file <- file.path(dir, "allele_data_global_max_0_filtered.csv")
  
  if (file.exists(csv_file)) {
    
    df <- read.csv(csv_file)
    df$run <- folder_name
    list_of_dfs[[folder_name]] <- df
  }
}

#merge genomic data into one df
merged_dfs <- bind_rows(list_of_dfs)


# identify sampleID with >=50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
good_sampleID <- merged_dfs %>%
  
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  
  # Filter rows where reads >= 100
  filter(reads >= 100) %>%   
  
  # Filter NIDA where n_loci >= 50
  summarise(n_loci = n_distinct(locus)) %>%
  filter(n_loci >= 50) %>%
  
  pull(sampleID)


# Keep good quality samples
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% good_sampleID, ]


# keep 1A amps
merged_dfs <- merged_dfs[grepl("-1A$", merged_dfs$locus),]


#create allele column
merged_dfs$allele <- paste0(merged_dfs$locus, "__", merged_dfs$pseudo_cigar)

# add run name to sampleID
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)


############################### REMOVE MASKING
#IGNORE MASKING
merged_dfs$pseudo_cigar <-  gsub("\\d+\\+[^N]*N", "", merged_dfs$pseudo_cigar)
merged_dfs$pseudo_cigar <- ifelse(merged_dfs$pseudo_cigar == "" | is.na(merged_dfs$pseudo_cigar), ".", merged_dfs$pseudo_cigar)

#create allele column
merged_dfs$allele <- paste0(merged_dfs$locus, "__", merged_dfs$pseudo_cigar)

merged_dfs <- merged_dfs %>%
  select(sampleID, locus, allele,reads, run)

#recalc reads
merged_dfs <- merged_dfs %>%
  group_by(run, sampleID, locus, allele) %>%  # Group by sampleID, locus, and allele
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")  # Sum reads and drop grouping

#recalc norm.reads.locus
merged_dfs <- merged_dfs %>%
  group_by(sampleID, locus) %>%  # Group by sampleID and locus
  mutate(norm.reads.locus = reads / sum(reads, na.rm = TRUE)) %>%  # Calculate the proportion
  ungroup() 

# remove icae1 run (cARLA); exclude bachmi
#merged_dfs <- merged_dfs[!grepl("ICAE_NextSeq01|MDACD_NextSeq02|ASINT_NextSeq04", merged_dfs$sampleID, ignore.case = TRUE),]


# merged_dfs <- merged_dfs[!grepl("I=", merged_dfs$allele),] #remove alleles with I (insertion)
# merged_dfs <- merged_dfs[!grepl("D=", merged_dfs$allele),] #remove alleles with D (deletion)
# merged_dfs <- merged_dfs[!merged_dfs$reads < 10,] #remove alleles with low read counts

controls_data <- merged_dfs[grepl("3d7", merged_dfs$sampleID, ignore.case = TRUE) &
                              !grepl("plus|hb3|dd2", merged_dfs$sampleID, ignore.case = TRUE),] 
                              #,]& !grepl("000000000", merged_dfs$sampleID, ignore.case = TRUE),] # remove CISM runs


TOTAL_RUNS <- length(unique(controls_data$run))
TOTAL_CONTROLS <- length(unique(controls_data$sampleID))



############################################################
###############3 CONTAMINANTS EDA ###################------------
############################################################

CONTAMINANTS <- controls_data[!grepl("\\.", controls_data$allele),]

CONTAMINANTS_RUNS <- length(unique(CONTAMINANTS$run))
CONTAMINANTS_CONTROLS <- length(unique(CONTAMINANTS$sampleID))

print(paste(CONTAMINANTS_RUNS, "out of", TOTAL_RUNS, "runs had contaminants"))
print(paste(CONTAMINANTS_CONTROLS, "out of", TOTAL_CONTROLS, "controls across runs had contaminants"))


# *** ARE CONTAMINANT ALLELES REPEATED ACROSS RUNS? ***

ALLELE_COUNT <- table(CONTAMINANTS$allele)
allele_df <- as.data.frame(ALLELE_COUNT)
colnames(allele_df) <- c("allele", "count")

allele_df$allele <- factor(allele_df$allele, levels = allele_df$allele[order(allele_df$count)])

# Create a ggplot histogram
contam_alleles <- ggplot(allele_df[allele_df$count > 1,], aes(x = allele, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(
    x = "Allele",
    y = "Count",
    title = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()

contam_alleles

ggsave("contam_alleles.png", contam_alleles, dpi = 300, height = 20, width = 15, bg = "white")



LOCI_COUNT <- table(CONTAMINANTS$locus)
loci_df <- as.data.frame(LOCI_COUNT)
colnames(loci_df) <- c("locus", "count")

loci_df$locus <- factor(loci_df$locus, levels = loci_df$locus[order(loci_df$count)])

# Create a ggplot histogram
contam_loci <- ggplot(loci_df, aes(x = locus, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(
    x = "Locus",
    y = "# Contaminant Alleles",
    title = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()

contam_loci

ggsave("contam_loci.png", contam_loci, dpi = 300, height = 20, width = 15, bg = "white")



# *** NUMBER OF CONTAMINANT ALLELES PER CONTROL ***-------------------------

n_contams_per_run <- CONTAMINANTS %>%
  group_by(run, sampleID) %>%
  summarise(n_contams = length(unique(allele)))

n_contams_per_run

n_contams_per_run <- n_contams_per_run %>%
  mutate(sampleID = factor(sampleID, levels = sampleID[order(-n_contams)]))

set.seed(42069)

num_colors <- length(unique(n_contams_per_run$run))
color_palette <- rgb(runif(num_colors), runif(num_colors), runif(num_colors))

n_contams_per_run <- n_contams_per_run %>%
  mutate(sampleID = factor(sampleID, levels = sampleID[order(-n_contams)]))

contams1  <- ggplot(n_contams_per_run, aes(x = sampleID, y = n_contams, fill = run)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Sample ID",
    y = " # Contaminant Alleles",
    fill = "Run"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank()
  ) +
  guides(fill = FALSE)
  #guides(fill = guide_legend(ncol = 2))

contams1

ggsave("contams1.png", contams1, dpi = 300, height = 10, width = 10, bg = "white")


#check the pair with the same high number of contaminant alleles. is it the same? if so, likely mislabelling issue
CONTAMINANTS[CONTAMINANTS$sampleID == "N3D7100KA_S7__BOH22_Nextseq01",]$allele %in% CONTAMINANTS[CONTAMINANTS$sampleID == "N3D710kA_S56__BOH22_Nextseq01",]$allele

#do these come from the cism contaminant strain?
CONTAMINANTS[CONTAMINANTS$sampleID == "N3D7-10K_S7_L001__240530_M07977_0028_000000000-LBCV5",]$allele %in% CONTAMINANTS[CONTAMINANTS$sampleID == "N3D710kA_S56__BOH22_Nextseq01",]$allele

##
#visualize similaritites of samples with many contaminants through heatmap to track potential origin
SAMPLES_WITH_MANY_CONTAMS <- n_contams_per_run[n_contams_per_run$n_contams > 0,]$sampleID # SELECT MINIMUM AMOUNF OF CONTAMINANTS

# Createreshape2# Create a binary matrix for allele presence/absence
allele_matrix <- CONTAMINANTS[CONTAMINANTS$sampleID %in% SAMPLES_WITH_MANY_CONTAMS,] %>%
  select(sampleID, locus) %>%
  distinct() %>%
  mutate(present = 1) %>%
  spread(key = locus, value = present, fill = 0) %>%
  column_to_rownames("sampleID")

# Calculate the Jaccard distance between sampleIDs
dist_matrix <- vegdist(allele_matrix, method = "jaccard")

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Get the order of samples based on clustering
ordered_samples <- hc$labels[hc$order]

# Convert distance matrix into a data frame for ggplot
dist_df <- as.data.frame(as.matrix(dist_matrix))
dist_df$SampleID1 <- rownames(dist_df)
dist_df <- melt(dist_df, id.vars = "SampleID1")
dist_df <- dist_df %>%
  rename(SampleID2 = variable, distance = value) %>%
  mutate(
    SampleID1 = factor(SampleID1, levels = ordered_samples),
    SampleID2 = factor(SampleID2, levels = ordered_samples)
  )

# Plot heatmap using ggplot2
hist_jaccard_all <- ggplot(dist_df, aes(x = SampleID1, y = SampleID2, fill = distance)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high = "blue", limits = c(0, 1)) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "Jaccard's Distance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )

hist_jaccard_all

ggsave("hist_jaccard_all.png", hist_jaccard_all, dpi = 300, height = 12, width = 15, bg = "white")


##
#visualize similaritites of samples with many contaminants through heatmap to track potential origin
SAMPLES_WITH_MANY_CONTAMS <- n_contams_per_run[n_contams_per_run$n_contams > 10,]$sampleID # SELECT MINIMUM AMOUNF OF CONTAMINANTS

# Createreshape2# Create a binary matrix for allele presence/absence
allele_matrix <- CONTAMINANTS[CONTAMINANTS$sampleID %in% SAMPLES_WITH_MANY_CONTAMS,] %>%
  select(sampleID, locus) %>%
  distinct() %>%
  mutate(present = 1) %>%
  spread(key = locus, value = present, fill = 0) %>%
  column_to_rownames("sampleID")

# Calculate the Jaccard distance between sampleIDs
dist_matrix <- vegdist(allele_matrix, method = "jaccard")

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Get the order of samples based on clustering
ordered_samples <- hc$labels[hc$order]

# Convert distance matrix into a data frame for ggplot
dist_df <- as.data.frame(as.matrix(dist_matrix))
dist_df$SampleID1 <- rownames(dist_df)
dist_df <- melt(dist_df, id.vars = "SampleID1")
dist_df <- dist_df %>%
  rename(SampleID2 = variable, distance = value) %>%
  mutate(
    SampleID1 = factor(SampleID1, levels = ordered_samples),
    SampleID2 = factor(SampleID2, levels = ordered_samples)
  )

# Plot heatmap using ggplot2
hist_jaccard_10plus <- ggplot(dist_df, aes(x = SampleID1, y = SampleID2, fill = distance)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high = "blue", limits = c(0, 1)) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "Jaccard's Distance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )

hist_jaccard_10plus

ggsave("hist_jaccard_10plus.png", hist_jaccard_10plus, dpi = 300, height = 9, width = 11, bg = "white")


# *** FREQUENCY DISTRIBUTION OF CONTAMINANT ALLELES *** -------------------

set.seed(420)

num_colors <- length(unique(CONTAMINANTS$sampleID))

color_palette <- rgb(runif(num_colors), runif(num_colors), runif(num_colors))

ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = sampleID)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity", color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "In-sample allele frequency",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())+
  guides(fill = guide_legend(ncol = 1))

contams2 <- ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = sampleID)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "In-sample allele frequency",
    y = "Non-reference alleles"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())+
  guides(fill = guide_legend(ncol = 1))+
  theme(
    legend.position = "none"
  ) +
  facet_wrap(~run, scales = "free_y")

contams2

ggsave("contams2.png", contams2, dpi = 300, height = 10, width = 17, bg = "white")


# controls with contaminant alleles that have a freq = 1
CONTAMINATED_CONTROLS_FREQ1 <- length(unique(CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]$sampleID))
RUNS_WITH_CONTAMINATED_CONTROLS_FREQ1  <- length(unique(CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]$run))

print(paste("There were", CONTAMINATED_CONTROLS_FREQ1, "3D7 controls with at least 1 fixed contaminant allele (norm.reads.locus = 1) across", RUNS_WITH_CONTAMINATED_CONTROLS_FREQ1, "runs"))



# *** CONTAMINANT ALLELES IN FIELD SAMPLES FROM EACH RUN?  *** -------------------

# remove controls from original data
field_data <- merged_dfs[!merged_dfs$sampleID %in% CONTAMINANTS$sampleID,]

SAMPLES <- unique(CONTAMINANTS$sampleID)

contam_procedence_results <- list()

for (sample in SAMPLES) {

  subsample_CONTAMINANTS_data <- CONTAMINANTS[CONTAMINANTS$sampleID == sample,]

  run <- strsplit(sample, "__")[[1]][2]

  subsample_filed_data <- field_data[field_data$run == run,]
  
  comparison_result <- unique(subsample_CONTAMINANTS_data$allele) %in% unique(subsample_filed_data$allele)
  
  percentage_contams_in_field_samples <- mean(comparison_result) * 100

  n_contams_in_control <- length(comparison_result)

  contam_procedence_results[[sample]] <- data.frame(
    run = run,
    sampleID = sample,
    n_contams_in_control = n_contams_in_control,
    percentage_contams_in_field_samples_from_run = percentage_contams_in_field_samples
  )
}

contam_procedence_results <- do.call(rbind, contam_procedence_results)

contam_procedence_results$sampleID <- sub("__.*", "", contam_procedence_results$sampleID)
rownames(contam_procedence_results) <- NULL

print(contam_procedence_results)

write.csv(contam_procedence_results, "contam_procedence_results.csv")


# ggplot(contam_procedence_results, aes(x = percentage_contams_in_field_samples_from_run)) +
#   geom_histogram(bins = 50, color = "black", fill = "blue", alpha = 0.7) +
#   labs(
#     title = "",
#     x = "Percentage of Contaminants in Field Samples",
#     y = "# Controls"
#   ) +
#   theme_minimal()

# % contaminants in field samples
contams3 <- ggplot(contam_procedence_results, aes(x = sampleID, y = percentage_contams_in_field_samples_from_run, fill = run)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "3D7 Controls",
    y = "% Contaminants in Field Samples",
    fill = "Run"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank()
  )+
  facet_wrap(~run, scales = "free_x")+
  guides(fill = FALSE)

contams3

ggsave("contams3.png", contams3, dpi = 300, height = 15, width = 15, bg = "white")


# number of contaminants in field samples
contam_procedence_results$n_contams_in_field_samples <- contam_procedence_results$n_contams_in_control * (contam_procedence_results$percentage_contams_in_field_samples_from_run/100)

contams4 <- ggplot(contam_procedence_results, aes(x = sampleID, y = n_contams_in_field_samples, fill = run)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "3D7 Controls",
    y = "# Contaminants in Field Samples",
    fill = "Run"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank()
  )+
  facet_wrap(~run, scales = "free_x")+
  guides(fill = FALSE)

contams4

ggsave("contams4.png", contams4, dpi = 300, height = 15, width = 15, bg = "white")


# *** MISSING (REF) ALLELES IN 3D7 CONTROLS  *** -------------------

all_ref_alleles <- unique(controls_data[grepl("_.$", controls_data$allele),]$locus)
n_expected_alleles <- length(all_ref_alleles) # number of expected alleles for 3d7 controls  == number of loci

allele_counts <- controls_data %>%
  filter(grepl("_.$", allele)) %>%
  group_by(sampleID) %>%
  summarise(n_correctly_sequenced_ref_alleles = length(unique(allele)))

allele_counts$missing_alleles <-  n_expected_alleles - allele_counts$n_correctly_sequenced_ref_alleles

allele_counts <- allele_counts %>%
  separate(sampleID, into = c("sampleID", "run"), sep = "__")


missing <- ggplot(allele_counts, aes(x = sampleID, y = missing_alleles, fill = run)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "3D7 Controls",
    y = "# Missing Alleles",
    fill = "Run"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank()
  )+
  facet_wrap(~run, scales = "free_x")+
  guides(fill = FALSE)

missing

ggsave("missing.png", missing, dpi = 300, height = 15, width = 15, bg = "white")


### PROPOSED THRESHOLDS

#remove alleles with freq of 1 (possible mislabelling of controls)
CONTAMINANTS_less_than_1 <- CONTAMINANTS[CONTAMINANTS$norm.reads.locus< 1,]

# interpret as for 50%, 50% of contaminants appear at a freq of n;
contam_thresholds <- as.data.frame(t(t(quantile(CONTAMINANTS_less_than_1$norm.reads.locus, probs= c(0.5, 0.75, 0.9, 0.95, 0.99)))))

colnames(contam_thresholds)[1] <- "MAF_threshold"

write.csv(contam_thresholds, "comtam_thresholds.csv")


thresholds <- ggplot(CONTAMINANTS_less_than_1, aes(x = norm.reads.locus)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity", color = "black", fill = "#69b3a2") +
  geom_vline(data = as.data.frame(contam_thresholds), 
             aes(xintercept = MAF_threshold, color = factor(rownames(contam_thresholds))), 
             linetype = "solid", size = 1) +
  scale_color_manual(name = "Thresholds", values = c("red", "blue", "green", "purple", "orange")) +
  labs(
    title = "",
    x = "In-sample allele frequency",
    y = "Contaminant alleles"
  ) +
  theme_minimal() +
  guides(color = guide_legend(title = "Thresholds")) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

thresholds

ggsave("thresholds.png", thresholds, height = 5, width = 8, dpi = 300, bg = "white")


# missing alleles
ggplot(allele_counts, aes(x = missing_alleles)) +
  geom_histogram(
    binwidth = 1,  # Adjust binwidth for your data
    fill = "#69b3a2",  # Custom fill color
    color = "black",   # Border color
    alpha = 0.8        # Transparency
  ) +
  labs(
    title = "",
    x = "Contaminant Alleles",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +  # Clean minimal theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center and bold title
    axis.title = element_text()
  )

#interpret backwards (when 0.1, 90% of samples have n missing alleles)
quantile(allele_counts$missing_alleles, probs= c(0.1, 0.25, 0.5, 0.75, 0.9))

