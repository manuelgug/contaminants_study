# Contamination Analysis Script for Genomic Data

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(vegan)
library(tibble)
library(reshape2)
library(stringr)
library(forcats)
library(ggpubr)

# ============================================================================
# 0) PARAMETERS AND DATA IMPORT
# ============================================================================

# Specify folder containing all filtered genomic data
main_dir <- "../amp_performance/ALL_SEQ_RESULTS_FILTERED"

# Import genomic data from all subdirectories
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

# Merge all genomic data into one dataframe
merged_dfs <- bind_rows(list_of_dfs)

# Add run name to sampleID to avoid aggregating samples with identical names
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)

# ============================================================================
# 1) QUALITY FILTERING
# ============================================================================

# Identify high-quality samples (>=50 loci with >= 100 reads)
good_sampleID <- merged_dfs %>%
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads), .groups = "drop") %>%
  filter(reads >= 100) %>%
  group_by(sampleID) %>%
  summarise(n_loci = n_distinct(locus), .groups = "drop") %>%
  filter(n_loci >= 50) %>%
  pull(sampleID)

# Keep only good quality samples
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% good_sampleID, ]

# Keep samples with > 10,000 reads
read_counts_above_10K <- merged_dfs %>%
  group_by(sampleID) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  filter(total_reads > 10000)

merged_dfs <- merged_dfs[merged_dfs$sampleID %in% read_counts_above_10K$sampleID, ]

# Keep only pool 1A
merged_dfs <- merged_dfs[grepl("-1A$", merged_dfs$locus), ]

# ============================================================================
# 2) DATA PROCESSING
# ============================================================================

# Convert masking to reference
merged_dfs$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", merged_dfs$pseudo_cigar)
merged_dfs$pseudo_cigar <- ifelse(merged_dfs$pseudo_cigar == "" | is.na(merged_dfs$pseudo_cigar), 
                                  ".", merged_dfs$pseudo_cigar)

# Aggregate unmasked sequences
merged_dfs_agg <- merged_dfs %>%
  group_by(sampleID, locus, pseudo_cigar, run, Category) %>%
  summarise(reads = sum(reads), 
            norm.reads.locus = sum(norm.reads.locus), 
            .groups = "drop")

# Remove indels
merged_dfs_agg <- merged_dfs_agg[!grepl("I=|D=", merged_dfs_agg$pseudo_cigar), ]

# Apply 1% MAF filter
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$norm.reads.locus > 0.01, ]

# Remove specific pools
merged_dfs_agg <- merged_dfs_agg[!grepl("-1AB$|1B2", merged_dfs_agg$locus), ]

# Remove known contaminated runs
merged_dfs_agg <- merged_dfs_agg[!merged_dfs_agg$run %in% 
                                   c("ASINT_NextSeq01", "ICAE_NextSeq01_demux_"), ]

# Remove alleles with less than 10 reads
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$reads > 10, ]

# ============================================================================
# 3) CREATE USEFUL VARIABLES
# ============================================================================

# Create pool variable
merged_dfs_agg$pool <- str_extract(merged_dfs_agg$locus, "[^-]+$")

# Create allele column
merged_dfs_agg$allele <- paste0(merged_dfs_agg$locus, "__", merged_dfs_agg$pseudo_cigar)

# Add lab variable
merged_dfs_agg$lab <- ifelse(grepl("000000000", merged_dfs_agg$run), "CISM", "ISG")

# ============================================================================
# 4) EXTRACT AND PROCESS CONTROL DATA
# ============================================================================

# Extract 3D7 controls
controls_data <- merged_dfs_agg[grepl("3d7", merged_dfs_agg$sampleID, ignore.case = TRUE) &
                                  !grepl("plus|hb3|dd2", merged_dfs_agg$sampleID, ignore.case = TRUE), ]

# Remove known mislabeled controls
mislabeled_controls <- c("N3D7-10K_S7_L001__240530_M07977_0028_000000000-LBCV5",
                         "N3D7100KA_S7__BOH22_Nextseq01",
                         "N3D710kA_S56__BOH22_Nextseq01",
                         "N3D7-10K_S7_L001__231129_M07977_0018_000000000-KHHTK")

controls_data <- controls_data[!controls_data$sampleID %in% mislabeled_controls, ]

# Check control read counts
controls_data %>%
  group_by(sampleID) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  arrange(total_reads)

# Assign parasitemia categories
controls_data$parasitemia <- ifelse(grepl("100K|10K", controls_data$sampleID, ignore.case = TRUE), 
                                    "High", "Low")

# Calculate summary statistics
TOTAL_RUNS <- length(unique(controls_data$run))
TOTAL_CONTROLS <- length(unique(controls_data$sampleID))

total_runs_per_lab <- table(distinct(controls_data[c("run", "lab")])$lab)
total_controls_per_lab <- table(distinct(controls_data[c("sampleID", "lab")])$lab)

print("Total runs per lab:")
print(total_runs_per_lab)
print("Total controls per lab:")
print(total_controls_per_lab)

# ============================================================================
# 5) CONTAMINATION ANALYSIS
# ============================================================================

# Identify contaminants (non-reference alleles)
CONTAMINANTS <- controls_data[!grepl("\\.", controls_data$allele), ]

CONTAMINANTS_RUNS <- length(unique(CONTAMINANTS$run))
CONTAMINANTS_CONTROLS <- length(unique(CONTAMINANTS$sampleID))

print(paste(CONTAMINANTS_RUNS, "out of", TOTAL_RUNS, "runs had contaminants"))
print(paste(CONTAMINANTS_CONTROLS, "out of", TOTAL_CONTROLS, "controls had contaminants"))

# Contamination by lab
contam_runs_per_lab <- table(distinct(CONTAMINANTS[c("run", "lab")])$lab)
contam_controls_per_lab <- table(distinct(CONTAMINANTS[c("sampleID", "lab")])$lab)

print("Contamination rates by lab:")
print(contam_runs_per_lab / total_runs_per_lab)
print(contam_controls_per_lab / total_controls_per_lab)

# Create contamination summary
contam_summ <- CONTAMINANTS %>%
  group_by(lab, parasitemia) %>%
  summarise(nonref_alleles = length(unique(allele)),
            nonref_loci = length(unique(locus)),
            controls_nonref = length(unique(sampleID)),
            .groups = "drop")

control_summ <- controls_data %>%
  group_by(lab, parasitemia) %>%
  summarise(controls_total = length(unique(sampleID)), .groups = "drop")

all_summ <- merge(control_summ, contam_summ, by = c("lab", "parasitemia"))
all_summ$perc_controls_nonref <- round(all_summ$controls_nonref / all_summ$controls_total * 100, 1)

# Save summary
write.csv(all_summ, "summary_contaminants.csv", row.names = FALSE)

# ============================================================================
# 6) VISUALIZATION: CONTAMINATION PATTERNS
# ============================================================================

# Plot 1: Contamination by parasitemia and lab
CONTAMINANTS_SUMMARY_2 <- CONTAMINANTS %>%
  group_by(lab, pool, parasitemia, run, sampleID) %>%
  summarise(unique_nonref = length(unique(allele)), .groups = "drop")

parasitemia_plots2 <- ggplot(CONTAMINANTS_SUMMARY_2, aes(x = pool, y = unique_nonref, fill = parasitemia)) +
  geom_boxplot() +
  facet_wrap(~lab) +
  scale_fill_manual(values = c("High" = "orchid", "Low" = "lightblue")) +
  labs(x = "Pool", y = "Unique Non-Ref Alleles", fill = "Parasitemia") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  stat_compare_means(aes(group = parasitemia), method = "kruskal.test", 
                     label = "p.signif", size = 5)

ggsave("parasitemia_plots2.png", parasitemia_plots2, dpi = 300, height = 6, width = 9, bg = "white")

# Plot 2: Combined parasitemia comparison
parasitemia_combined_plot <- ggplot(CONTAMINANTS_SUMMARY_2, aes(x = pool, y = unique_nonref, fill = lab)) +
  geom_boxplot() +
  facet_wrap(~parasitemia) +
  scale_fill_manual(values = c("ISG" = "orange2", "CISM" = "green3")) +
  labs(x = "Pool", y = "Unique Non-Ref Alleles", fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  stat_compare_means(aes(group = lab), method = "kruskal.test", 
                     label = "p.signif", size = 5)

ggsave("parasitemia_combined_plot.png", parasitemia_combined_plot, dpi = 300, height = 6, width = 9, bg = "white")

# Plot 3: Allele counts
ALLELE_COUNT <- CONTAMINANTS %>%
  group_by(lab, pool, parasitemia, allele) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) %>%
  mutate(allele = fct_reorder(allele, -as.numeric(factor(pool))))

allele_counts_plot <- ggplot(ALLELE_COUNT, aes(x = allele, y = count, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_grid(lab ~ parasitemia, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +
  labs(x = "Non-Reference Allele", y = "Count", fill = "Pool") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  coord_flip()

ggsave("allele_counts_plot.png", allele_counts_plot, dpi = 300, height = 5, width = 10, bg = "white")

# ============================================================================
# 7) CONTAMINATION SOURCE ANALYSIS
# ============================================================================

# Save contaminant data
CONTAMINANTS_output <- CONTAMINANTS %>%
  separate(sampleID, into = c("sampleID", "run"), sep = "__", extra = "merge") %>%
  select(-allele) %>%
  arrange(lab, run, pool, sampleID, Category)

write.csv(CONTAMINANTS_output, "CONTAMINANTS_output.csv", row.names = FALSE)

# Number of contaminants per control
n_contams_per_run <- CONTAMINANTS %>%
  group_by(run, sampleID, lab, pool, parasitemia) %>%
  summarise(n_contams = length(unique(allele)), .groups = "drop")

# Plot contamination per control
n_contams_per_run$sampleID_plot <- gsub("__.*", "", n_contams_per_run$sampleID)

contams1 <- ggplot(n_contams_per_run, aes(x = sampleID, y = n_contams, fill = pool)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(lab ~ parasitemia, scales = "free", space = "free") +
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +
  labs(x = "Controls", y = "Non-reference alleles", fill = "Pool") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  coord_flip()

ggsave("contams1.png", contams1, dpi = 300, height = 12, width = 15, bg = "white")

# ============================================================================
# 8) JACCARD DISTANCE ANALYSIS
# ============================================================================

# Calculate Jaccard distance for samples with contamination
SAMPLES_WITH_CONTAMS <- n_contams_per_run %>%
  filter(n_contams > 0) %>%
  pull(sampleID)

# Jaccard distance calculation
calculate_jaccard_distance <- function(contaminant_data, sample_ids) {
  dist_df_all <- contaminant_data %>%
    filter(sampleID %in% sample_ids) %>%
    group_by(pool) %>%
    do({
      pool_data <- .
      allele_matrix <- pool_data %>%
        select(sampleID, locus) %>%
        distinct() %>%
        mutate(present = 1) %>%
        pivot_wider(names_from = locus, values_from = present, values_fill = 0) %>%
        column_to_rownames("sampleID")
      
      if (nrow(allele_matrix) > 1) {
        dist_matrix <- vegdist(allele_matrix, method = "jaccard")
        df <- as.data.frame(as.matrix(dist_matrix))
        df$SampleID1 <- rownames(df)
        df_long <- df %>%
          pivot_longer(cols = -SampleID1, names_to = "SampleID2", values_to = "distance") %>%
          mutate(SampleID1 = factor(SampleID1, levels = rownames(allele_matrix)),
                 SampleID2 = factor(SampleID2, levels = rownames(allele_matrix)))
        df_long
      } else {
        data.frame()
      }
    }) %>%
    ungroup()
  
  return(dist_df_all)
}

# Calculate and plot Jaccard distances
dist_df_all <- calculate_jaccard_distance(CONTAMINANTS, SAMPLES_WITH_CONTAMS)

if (nrow(dist_df_all) > 0) {
  hist_jaccard_all <- ggplot(dist_df_all, aes(x = SampleID1, y = SampleID2, fill = distance)) +
    geom_tile() +
    scale_fill_gradient(low = "orange", high = "blue", limits = c(0, 1)) +
    labs(x = "", y = "", fill = "Jaccard's Distance") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
          axis.title.x = element_text(size = 14),
          strip.text = element_text(size = 17, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    facet_wrap(~pool, scales = "free", ncol = 2)
  
  ggsave("jaccard_distance_heatmap.png", hist_jaccard_all, dpi = 300, height = 10, width = 12, bg = "white")
}

# ============================================================================
# 9) FREQUENCY DISTRIBUTION ANALYSIS
# ============================================================================

# Plot frequency distribution of contaminants
contams2 <- ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = pool)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
  facet_grid(lab ~ parasitemia, scales = "free", space = "free") +
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +
  labs(x = "In-sample allele frequency", y = "Non-reference alleles", fill = "Pool") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave("frequency_distribution.png", contams2, dpi = 300, height = 10, width = 17, bg = "white")

# ============================================================================
# 10) CONTAMINATION SOURCE TRACKING
# ============================================================================

# Remove controls from original data for field sample comparison
field_data <- merged_dfs_agg[!merged_dfs_agg$sampleID %in% CONTAMINANTS$sampleID, ]
field_data <- field_data[field_data$run %in% unique(controls_data$run), ]

# Analyze contamination sources
analyze_contamination_sources <- function(contaminant_data, field_data) {
  SAMPLES <- unique(contaminant_data$sampleID)
  results <- list()
  
  for (sample in SAMPLES) {
    subsample_contaminants <- contaminant_data[contaminant_data$sampleID == sample, ]
    run <- strsplit(sample, "__")[[1]][2]
    subsample_field <- field_data[field_data$run == run, ]
    
    # Find contaminants in field samples from same run
    contaminants_in_field <- intersect(subsample_contaminants$allele, subsample_field$allele)
    contaminants_not_in_field <- setdiff(subsample_contaminants$allele, subsample_field$allele)
    
    # Cross-contamination from other runs
    cross_contams <- intersect(contaminants_not_in_field, field_data$allele)
    cross_contam_runs <- field_data[field_data$allele %in% contaminants_not_in_field, ]$run
    
    # Unknown origin contaminants
    unknown_contams <- setdiff(contaminants_not_in_field, field_data$allele)
    
    results[[sample]] <- data.frame(
      run = run,
      sampleID = sample,
      n_nonref_in_control = length(unique(subsample_contaminants$allele)),
      n_nonref_in_field = length(unique(contaminants_in_field)),
      n_cross_nonref = length(unique(cross_contams)),
      n_unknown_nonref = length(unique(unknown_contams))
    )
  }
  
  return(do.call(rbind, results))
}

# Analyze and save results
contam_source_results <- analyze_contamination_sources(CONTAMINANTS, field_data)

# Add parasitemia and lab information
contam_source_results <- left_join(contam_source_results, 
                                   unique(CONTAMINANTS[c("sampleID", "parasitemia", "lab")]), 
                                   by = "sampleID")

contam_source_results$sampleID <- sub("__.*", "", contam_source_results$sampleID)
contam_source_results <- contam_source_results %>%
  arrange(lab, parasitemia, run)

write.csv(contam_source_results, "contamination_source_analysis.csv", row.names = FALSE)


# Reshape data to long format for ggplot
contam_procedence_long <- contam_source_results %>%
  select(sampleID, run, lab, parasitemia, n_nonref_in_field, n_cross_nonref, n_unknown_nonref) %>%
  pivot_longer(cols = c(n_nonref_in_field, n_cross_nonref, n_unknown_nonref),
               names_to = "nonref_type",
               values_to = "count")


contam_procedence_long$sampleID_plot <- paste0(contam_procedence_long$sampleID, "__", contam_procedence_long$run)

contam_procedence_long <- contam_procedence_long %>%
  mutate(sampleID_plot = fct_reorder(sampleID_plot, -as.numeric(factor(run))))

# Create stacked bar plot
contams3 <- ggplot(contam_procedence_long, aes(x = sampleID_plot, y = count, fill = nonref_type)) +
  geom_bar(stat = "identity", position = "stack") +  
  facet_grid(lab ~ parasitemia, scales = "free", space="free") + 
  scale_fill_manual(values = c( "#e31a1c", "#1f78b4", "black")) +  
  labs(
    x = "Control",
    y = "Non-reference alleles",
    fill = "Source",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),
    axis.title.x = element_text(size = 14),        
    axis.title.y = element_text(size = 14),     
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  ) +
  coord_flip()

contams3



ggsave("contams3210.png", contams3, dpi = 300, height = 12, width = 20, bg = "white")




# ============================================================================
# 11) THRESHOLD ANALYSIS
# ============================================================================

# Calculate contamination thresholds (excluding frequency = 1)
CONTAMINANTS_less_than_1 <- CONTAMINANTS[CONTAMINANTS$norm.reads.locus < 1, ]

# Calculate thresholds by lab, parasitemia, and pool
contam_thresholds_pools <- CONTAMINANTS_less_than_1 %>%
  group_by(lab, parasitemia, pool) %>%
  summarise(`50%` = quantile(norm.reads.locus, probs = 0.5, na.rm = TRUE),
            `75%` = quantile(norm.reads.locus, probs = 0.75, na.rm = TRUE),
            `90%` = quantile(norm.reads.locus, probs = 0.9, na.rm = TRUE),
            `95%` = quantile(norm.reads.locus, probs = 0.95, na.rm = TRUE),
            `99%` = quantile(norm.reads.locus, probs = 0.99, na.rm = TRUE),
            .groups = "drop")

write.csv(contam_thresholds_pools, "contamination_thresholds_by_pool.csv", row.names = FALSE)

# Plot thresholds
contam_thresholds_long <- contam_thresholds_pools %>%
  pivot_longer(cols = `50%`:`99%`, names_to = "percentile", values_to = "MAF_threshold")

thresholds_plot <- ggplot(CONTAMINANTS_less_than_1, aes(x = norm.reads.locus)) +
  geom_histogram(bins = 100, alpha = 0.5, fill = "steelblue", color = "white") +
  geom_vline(data = contam_thresholds_long,
             aes(xintercept = MAF_threshold, color = percentile),
             linetype = "solid", size = 1) +
  facet_grid(lab ~ parasitemia) +
  scale_color_manual(name = "Thresholds", values = c("red", "blue", "green", "purple", "orange")) +
  labs(x = "In-sample allele frequency", y = "Non-reference alleles") +
  theme_minimal() +
  guides(color = guide_legend(title = "% Non-ref alleles\neliminated")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 17, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave("contamination_thresholds.png", thresholds_plot, height = 8, width = 12, dpi = 300, bg = "white")

# ============================================================================
# 12) FINAL SUMMARY
# ============================================================================

print("Analysis complete!")
print(paste("Total runs analyzed:", TOTAL_RUNS))
print(paste("Total controls analyzed:", TOTAL_CONTROLS))
print(paste("Runs with contamination:", CONTAMINANTS_RUNS))
print(paste("Controls with contamination:", CONTAMINANTS_CONTROLS))
print(paste("Contamination rate:", round(CONTAMINANTS_CONTROLS/TOTAL_CONTROLS*100, 1), "%"))

# Save session info
sessionInfo()

