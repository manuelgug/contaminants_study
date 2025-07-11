library(dplyr)
library(tidyr)
library(data.table)

# Read truth data once (assumed to be the same across runs)
truth <- read.csv("pool1A_truth.csv")
truth$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", truth$pseudo_cigar)
truth$pseudo_cigar <- ifelse(truth$pseudo_cigar == "" | is.na(truth$pseudo_cigar), ".", truth$pseudo_cigar)
truth$allele <- paste0(truth$locus, "__", truth$pseudo_cigar)
truth <- unique(truth)

# Get all subdirectories
input_dirs <- list.dirs("../../amp_performance/ALL_SEQ_RESULTS_FILTERED/", recursive = FALSE)

# Parameters
min_allele_read_count <- 10
maf_filter <- 0.02

# Cleaning function
clean_data <- function(data, maf_filter, min_allele_read_count) {
  good_sampleID <- data %>%
    group_by(sampleID, locus) %>%
    summarize(reads = sum(reads), .groups = "drop") %>%
    filter(reads >= 100) %>%
    group_by(sampleID) %>%
    summarise(n_loci = n_distinct(locus), .groups = "drop") %>%
    filter(n_loci >= 50) %>%
    pull(sampleID)
  
  data <- data[data$sampleID %in% good_sampleID, ]
  
  read_counts <- data %>%
    group_by(sampleID) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    filter(total_reads > 10000)
  
  data <- data[data$sampleID %in% read_counts$sampleID, ]
  data <- data[grepl("-1A$", data$locus), ]
  data$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", data$pseudo_cigar)
  data$pseudo_cigar <- ifelse(data$pseudo_cigar == "" | is.na(data$pseudo_cigar), ".", data$pseudo_cigar)
  
  data_agg <- data %>%
    group_by(sampleID, locus, pseudo_cigar) %>%
    summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus), .groups = "drop")
  
  data_agg$allele <- paste0(data_agg$locus, "__", data_agg$pseudo_cigar)
  data_agg <- data_agg[!grepl("I=|D=", data_agg$allele), ]
  data_agg <- data_agg[data_agg$norm.reads.locus > maf_filter, ]
  data_agg <- data_agg[data_agg$reads > min_allele_read_count, ]
  data_agg <- data_agg %>% select(-pseudo_cigar)
  data_agg <- data_agg[!grepl("Undetermined", data_agg$sampleID, ignore.case = TRUE), ]
  
  return(data_agg)
}

# Helper for counting strain references
count_strains <- function(x) {
  sum(grepl("3D7", x, ignore.case = TRUE),
      grepl("DD2", x, ignore.case = TRUE),
      grepl("HB3", x, ignore.case = TRUE))
}

# Loop through each input folder
for (INPUT_FOLDER in input_dirs) {
  file_path <- file.path(INPUT_FOLDER, "allele_data_global_max_0_filtered.csv")
  if (!file.exists(file_path)) next
  
  message("Processing: ", INPUT_FOLDER)
  data <- read.csv(file_path)
  data_agg <- clean_data(data, maf_filter, min_allele_read_count)
  
  # Controls and field
  controls_3d7 <- data_agg[grepl("3d7", data_agg$sampleID, ignore.case = TRUE), ]
  controls_dd2 <- data_agg[grepl("dd2", data_agg$sampleID, ignore.case = TRUE), ]
  field <- data_agg[!grepl("3d7|dd2", data_agg$sampleID, ignore.case = TRUE), ]
  
  # Control truth
  truth_3d7 <- truth[grepl("3d7", truth$Strain, ignore.case = TRUE), ]
  truth_dd2 <- truth[grepl("dd2$", truth$Strain, ignore.case = TRUE), ]
  
  total_alleles_3d7 <- controls_3d7 %>% group_by(sampleID) %>% summarise(alleles_sequenced = n_distinct(allele), .groups = "drop")
  total_alleles_dd2 <- controls_dd2 %>% group_by(sampleID) %>% summarise(alleles_sequenced = n_distinct(allele), .groups = "drop")
  All_total_alleles <- rbind(total_alleles_3d7, total_alleles_dd2)
  
  # Nonref
  nonref_3d7 <- setdiff(controls_3d7$allele, truth_3d7$allele)
  nonref_dd2 <- setdiff(controls_dd2$allele, truth_dd2$allele)
  
  controls_3d7_nonref <- controls_3d7[controls_3d7$allele %in% nonref_3d7, ]
  controls_dd2_nonref <- controls_dd2[controls_dd2$allele %in% nonref_dd2, ]
  
  # Missing
  missing_3d7 <- controls_3d7 %>%
    group_by(sampleID) %>%
    summarise(missing_alleles = sum(!truth_3d7$allele %in% allele), .groups = "drop")
  
  missing_dd2 <- controls_dd2 %>%
    group_by(sampleID) %>%
    summarise(missing_alleles = sum(!truth_dd2$allele %in% allele), .groups = "drop")
  
  missing_all <- rbind(missing_3d7, missing_dd2)
  
  # Successful
  successful_3d7 <- controls_3d7 %>%
    group_by(sampleID) %>%
    summarise(successful_alleles = sum(truth_3d7$allele %in% allele),
              expected_alleles = length(unique(truth_3d7$allele)),
              .groups = "drop")
  
  successful_dd2 <- controls_dd2 %>%
    group_by(sampleID) %>%
    summarise(successful_alleles = sum(truth_dd2$allele %in% allele),
              expected_alleles = length(unique(truth_dd2$allele)),
              .groups = "drop")
  
  successful_all <- rbind(successful_3d7, successful_dd2)
  summary_all <- merge(All_total_alleles, successful_all, by = "sampleID")
  summary_all <- merge(summary_all, missing_all, by = "sampleID")
  
  # Source of nonref
  controls_3d7_nonref$source <- ifelse(controls_3d7_nonref$allele %in% field$allele, "field", "unknown")
  controls_dd2_nonref$source <- ifelse(controls_dd2_nonref$allele %in% field$allele, "field", "unknown")
  
  all_nonref <- rbindlist(list(controls_3d7_nonref, controls_dd2_nonref), use.names = TRUE, fill = TRUE)
  
  summary <- all_nonref %>%
    group_by(sampleID) %>%
    summarise(
      fixed_nonref_alleles = sum(norm.reads.locus == 1, na.rm = TRUE),
      nonref_source_field = sum(source == "field", na.rm = TRUE),
      nonref_source_unknown = sum(source == "unknown", na.rm = TRUE),
      .groups = "drop"
    )
  
  all_controls <- unique(c(controls_3d7$sampleID, controls_dd2$sampleID))
  
  SUMMARY_COMPLETE <- tibble(sampleID = all_controls) %>%
    left_join(summary, by = "sampleID") %>%
    replace_na(list(
      fixed_nonref_alleles = 0,
      nonref_source_field = 0,
      nonref_source_unknown = 0
    ))
  
  # Filter mixes
  SUMMARY_COMPLETE <- SUMMARY_COMPLETE %>% filter(sapply(sampleID, count_strains) == 1)
  all_nonref <- all_nonref %>% filter(sapply(sampleID, count_strains) == 1)
  
  # Add run info
  run_name <- basename(INPUT_FOLDER)
  SUMMARY_COMPLETE$run <- run_name
  all_nonref$run <- run_name
  
  # Merge summary
  SUMMARY_COMPLETE <- merge(SUMMARY_COMPLETE, summary_all, by = "sampleID")
  SUMMARY_COMPLETE <- SUMMARY_COMPLETE %>% select(sampleID, run, expected_alleles, alleles_sequenced,
                                                  successful_alleles, missing_alleles,
                                                  fixed_nonref_alleles, nonref_source_field, nonref_source_unknown)
  all_nonref <- all_nonref %>% select(sampleID, run, everything())
  
  # Output
  write.csv(SUMMARY_COMPLETE, paste0("summary_", run_name, ".csv"), row.names = FALSE)
  write.csv(all_nonref, paste0("all_nonref_alleles_", run_name, ".csv"), row.names = FALSE)
}


##############3


# K-means Clustering Analysis for UMAP Data
# Load required libraries
library(dplyr)
library(uwot)
library(ggplot2)
library(factoextra)
library(cluster)
library(gridExtra)

# ========================================
# STEP 1: Data Preparation (from your code)
# ========================================

# Define the folder where the summary files are stored
summary_folder <- "."

# List all files starting with "summary_" and ending in .csv
summary_files <- list.files(path = summary_folder, pattern = "^summary_.*\\.csv$", full.names = TRUE)

# Read and combine all summary files into a single data frame
all_summaries <- lapply(summary_files, read.csv) 

# Remove empty data frames
all_summaries <- all_summaries[sapply(all_summaries, nrow) > 0]

# Combine all non-empty summary data frames
all_summaries <- bind_rows(all_summaries)

# Optionally save the combined summary
write.csv(all_summaries, file = paste0(summary_folder, "ALL_summary_combined.csv"), row.names = FALSE)

# ========================================
# STEP 2: UMAP Analysis (from your code)
# ========================================

# Select the numeric feature columns from 'alleles_sequenced' onward
feature_data <- all_summaries %>%
  select(alleles_sequenced:nonref_source_unknown)

# Run UMAP
set.seed(123)  # for reproducibility
umap_result <- umap(feature_data, n_neighbors = 15, min_dist = 0.1)

# Create a data frame with UMAP results and add metadata
umap_df <- data.frame(
  sampleID = all_summaries$sampleID,
  run = all_summaries$run,
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2]
)

# ========================================
# STEP 3: K-means Clustering Analysis
# ========================================

# 3.1: Determine optimal number of clusters using multiple methods
print("=== DETERMINING OPTIMAL NUMBER OF CLUSTERS ===")

# Method 1: Elbow Method (Within Sum of Squares)
elbow_plot <- fviz_nbclust(umap_result, kmeans, method = "wss", k.max = 10) +
  ggtitle("Elbow Method for Optimal k") +
  theme_minimal()

# Method 2: Silhouette Method
silhouette_plot <- fviz_nbclust(umap_result, kmeans, method = "silhouette", k.max = 10) +
  ggtitle("Silhouette Method for Optimal k") +
  theme_minimal()

# Method 3: Gap Statistic
gap_plot <- fviz_nbclust(umap_result, kmeans, method = "gap_stat", k.max = 10) +
  ggtitle("Gap Statistic Method for Optimal k") +
  theme_minimal()

# Display all three plots
grid.arrange(elbow_plot, silhouette_plot, gap_plot, ncol = 2)

# 3.2: Calculate silhouette scores for different k values
k_range <- 2:20
silhouette_scores <- sapply(k_range, function(k) {
  kmeans_temp <- kmeans(umap_result, centers = k, nstart = 25)
  sil <- silhouette(kmeans_temp$cluster, dist(umap_result))
  mean(sil[, 3])
})

# Print silhouette scores
sil_df <- data.frame(k = k_range, silhouette_score = silhouette_scores)
print("Silhouette Scores for different k values:")
print(sil_df)

# Find optimal k based on highest silhouette score
optimal_k <- k_range[which.max(silhouette_scores)]
print(paste("Optimal k based on silhouette score:", optimal_k))

# 3.3: Perform k-means clustering with optimal k
set.seed(123)
kmeans_result <- kmeans(umap_result, centers = optimal_k, nstart = 25)

# Add cluster assignments to the dataframe
umap_df$cluster <- as.factor(kmeans_result$cluster)

# Print cluster summary
print("=== CLUSTER SUMMARY ===")
print(table(umap_df$cluster))
print(paste("Total samples:", nrow(umap_df)))
print(paste("Number of clusters:", optimal_k))
print(paste("Within-cluster sum of squares:", round(kmeans_result$tot.withinss, 2)))
print(paste("Between-cluster sum of squares:", round(kmeans_result$betweenss, 2)))
print(paste("Total sum of squares:", round(kmeans_result$totss, 2)))

# ========================================
# STEP 4: Visualization
# ========================================

# 4.1: Plot clusters on UMAP
cluster_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = paste("K-means Clustering (k =", optimal_k, ")"),
       subtitle = "Clusters overlaid on UMAP projection",
       x = "UMAP1", y = "UMAP2", color = "Cluster") +
  theme(legend.position = "right")

print(cluster_plot)

# 4.2: Compare with original run groupings
run_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = run)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Original Run Groupings",
       subtitle = "Colored by experimental run",
       x = "UMAP1", y = "UMAP2", color = "Run") +
  theme(legend.position = "right")

# Display both plots side by side
grid.arrange(cluster_plot, run_plot, ncol = 1)

# # 4.3: Create a comparison plot showing both cluster and run
# comparison_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
#   geom_point(aes(color = cluster, shape = run), size = 2, alpha = 0.8) +
#   scale_color_viridis_d() +
#   theme_minimal() +
#   labs(title = "Clusters vs. Runs",
#        subtitle = "Color = Cluster, Shape = Run",
#        x = "UMAP1", y = "UMAP2") +
#   theme(legend.position = "right")
# 
# print(comparison_plot)

# ========================================
# STEP 5: Cluster Characterization
# ========================================

# 5.1: Analyze cluster composition by run
print("=== CLUSTER COMPOSITION BY RUN ===")
cluster_run_table <- table(umap_df$cluster, umap_df$run)
print(cluster_run_table)

# Calculate proportions
cluster_run_prop <- prop.table(cluster_run_table, margin = 1)
print("Proportions within each cluster:")
print(round(cluster_run_prop, 3))

# 5.2: Statistical summary of original features by cluster
print("=== FEATURE MEANS BY CLUSTER ===")
feature_data_with_clusters <- feature_data
feature_data_with_clusters$cluster <- umap_df$cluster

cluster_means <- feature_data_with_clusters %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(cluster_means)

# 5.3: Create a heatmap of cluster characteristics
# Scale the features for better visualization
scaled_features <- scale(feature_data)
scaled_df <- data.frame(scaled_features, cluster = umap_df$cluster)

# Step 1: Compute cluster profiles
cluster_profiles <- scaled_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Step 2: Compute row order using hierarchical clustering
mat_for_clustering <- as.matrix(cluster_profiles %>% select(-cluster))
row_order <- hclust(dist(mat_for_clustering))$order

# Step 3: Reorder the factor levels of cluster by clustering order
cluster_profiles$cluster <- factor(cluster_profiles$cluster,
                                   levels = cluster_profiles$cluster[row_order])

# Step 4: Convert to long format for heatmap
cluster_profiles_long <- melt(cluster_profiles, id.vars = "cluster")

# Step 5: Plot heatmap with legend on the left
heatmap_plot <- ggplot(cluster_profiles_long, aes(x = variable, y = cluster, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "left"
  ) +
  labs(
    title = "Cluster Profiles Heatmap",
    subtitle = "Standardized feature means by cluster (rows sorted by similarity)",
    x = "Features", y = "Cluster", fill = "Z-score"
  )

print(heatmap_plot)

ggsave("heatmap.png", heatmap_plot, bg = "white", width = 8, height = 5)

# ========================================
# STEP 6: Save Results
# ========================================

# Save the results
umap_df_with_features <- cbind(umap_df, feature_data)
write.csv(umap_df_with_features, "clustered_results.csv", row.names = FALSE)

# Save cluster centers
cluster_centers <- data.frame(
  cluster = 1:optimal_k,
  UMAP1_center = kmeans_result$centers[, 1],
  UMAP2_center = kmeans_result$centers[, 2]
)
write.csv(cluster_centers, "cluster_centers.csv", row.names = FALSE)

print("=== ANALYSIS COMPLETE ===")
print("Results saved to:")
print("- clustered_results.csv (samples with cluster assignments)")
print("- cluster_centers.csv (cluster center coordinates)")
print(paste("- Optimal number of clusters:", optimal_k))
print(paste("- Silhouette score:", round(max(silhouette_scores), 3)))
