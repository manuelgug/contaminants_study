
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


####### 0) PARAMETERS AND INPUTS #######--------------------

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


########## 0) ADD RUN NAME TO SAMPLEID TO AVOID AGGREGATING SAMPLES WITH IDENTIAL NAME BUT IN DIFF RUNS ################-------
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)


########## 1) KEEP QUALITY SAMPLES ################-------
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


# keep samples with > 10000 reads: ANDRÉS
read_counts_above_10K_reads <- merged_dfs %>% group_by(sampleID) %>% summarise(total_reads = sum(reads)) 
read_counts_above_10K_reads <- read_counts_above_10K_reads[read_counts_above_10K_reads$total_reads > 10000,]
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% read_counts_above_10K_reads$sampleID, ]


### keep pool 1A
merged_dfs <- merged_dfs[grepl("-1A$", merged_dfs$locus),]


########## 2) CONVERT MASKING INTO REF ################-------
# ignore masking, turn it to ref (.)
merged_dfs$pseudo_cigar <-  gsub("\\d+\\+[^N]*N", "", merged_dfs$pseudo_cigar) #remove masking
merged_dfs$pseudo_cigar <- ifelse(merged_dfs$pseudo_cigar == "" | is.na(merged_dfs$pseudo_cigar), ".", merged_dfs$pseudo_cigar) # if empty, add "." since it was reference

#aggregate unmasked sequences
merged_dfs_agg <- merged_dfs %>% group_by(sampleID, locus, pseudo_cigar, run, Category) %>% summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus))


########## 3) REMOVE INDELS ################-------
# remove indels
merged_dfs_agg <- merged_dfs_agg[!grepl("I=", merged_dfs_agg$pseudo_cigar),] #remove alleles with I (insertion)
merged_dfs_agg <- merged_dfs_agg[!grepl("D=", merged_dfs_agg$pseudo_cigar),] #remove alleles with D (deletion)


########### 4) 1% MAF FILTER
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$norm.reads.locus > 0.01,]


########### 5) REMOVE POOLS 1AB AND 1B2
merged_dfs_agg <- merged_dfs_agg[!grepl("-1AB$|1B2", merged_dfs_agg$locus),]


########### 6) REMOVE KNOWN CONTAMINATED RUNS
merged_dfs_agg <- merged_dfs_agg[!merged_dfs_agg$run %in% c("ASINT_NextSeq01", "ICAE_NextSeq01_demux_"), ]


########### 7) REMOVE ALLELES WITH LESS THAN 10 READS
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$reads > 10,]



########## 2) CREATE/MODIFY USEFUL VARIABLES ################-------

#make pool variable
merged_dfs_agg$pool <-  str_extract(merged_dfs_agg$locus, "[^-]+$")

#create allele column
merged_dfs_agg$allele <- paste0(merged_dfs_agg$locus, "__", merged_dfs_agg$pseudo_cigar)

# add lab variable
merged_dfs_agg$lab <- ifelse(grepl("000000000", merged_dfs_agg$run), "CISM", "ISG")


########## 5) EXTRACT CONTROL DATA ################-------
controls_data <- merged_dfs_agg[grepl("3d7", merged_dfs_agg$sampleID, ignore.case = TRUE) &
                              !grepl("plus|hb3|dd2", merged_dfs_agg$sampleID, ignore.case = TRUE),] 
                              #,]& !grepl("000000000", merged_dfs_agg$sampleID, ignore.case = TRUE),] # remove CISM runs

#remove known mislabelled controls:
mislabelled_controls <- c("N3D7-10K_S7_L001__240530_M07977_0028_000000000-LBCV5", "N3D7100KA_S7__BOH22_Nextseq01", "N3D710kA_S56__BOH22_Nextseq01")
controls_data <- controls_data[!controls_data$sampleID %in% mislabelled_controls,]

controls_data %>% group_by(sampleID) %>% summarise(total_reads= sum(reads)) %>% arrange(total_reads) # check if all controls have > 10000 reads

#asign parasitemia categories. 10K and 100K = high; 1K and 100 = low.
controls_data$parasitemia <- ifelse(grepl("100K|10K", controls_data$sampleID, ignore.case = T), "High","Low")


TOTAL_RUNS <- length(unique(controls_data$run))
TOTAL_CONTROLS <- length(unique(controls_data$sampleID))

# total runs and controls per lab
total_runs_per_lab <- table(distinct(controls_data[c("run", "lab")])$lab); total_runs_per_lab
total_controls_per_lab <- table(distinct(controls_data[c("sampleID", "lab")])$lab); total_controls_per_lab


############################################################
###############3 CONTAMINANTS EDA ###################------------
############################################################

CONTAMINANTS <- controls_data[!grepl("\\.", controls_data$allele),]

CONTAMINANTS_RUNS <- length(unique(CONTAMINANTS$run))
CONTAMINANTS_CONTROLS <- length(unique(CONTAMINANTS$sampleID))

print(paste(CONTAMINANTS_RUNS, "out of", TOTAL_RUNS, "runs had contaminants"))
print(paste(CONTAMINANTS_CONTROLS, "out of", TOTAL_CONTROLS, "controls across runs had contaminants"))

# by lab
contam_runs_per_lab <- table(distinct(CONTAMINANTS[c("run", "lab")])$lab); contam_runs_per_lab
contam_controls_per_lab<- table(distinct(CONTAMINANTS[c("sampleID", "lab")])$lab); contam_controls_per_lab


contam_runs_per_lab/total_runs_per_lab #percentage contaminated runs per pool
contam_controls_per_lab/total_controls_per_lab #percentage contaminated runs per pool




# *** ARE CONTAMINANT ALLELES REPEATED ACROSS RUNS? ***

# CONTAMINANTNS_SUMMARY <- CONTAMINANTS %>% group_by(lab, pool, parasitemia) %>% summarise(unique_nonref =length(unique(allele)))
# 
# parasitemia_plots <- ggplot(CONTAMINANTNS_SUMMARY, aes(x = pool, y = unique_nonref, fill = parasitemia)) +
#   geom_bar(stat = "identity", position = "dodge") +  # Dodged bars for comparison
#   facet_wrap(~lab) +  # Facet by lab
#   scale_fill_manual(values = c("High_Parasitemia" = "orange2", "Low_Parasitemia" = "green3")) +  # Custom colors
#   labs(
#     x = "Pool",
#     y = "Unique Non-Ref Alleles",
#     fill = "",
#     title = ""
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 15, angle = 0, hjust = 1),
#     axis.title.x = element_text(size = 14),        
#     strip.text = element_text(size = 17, face = "bold"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
#   )
# 
# parasitemia_plots
# 
# ggsave("parasitemia_plots.png", parasitemia_plots, dpi = 300, height = 5, width = 8, bg = "white")


CONTAMINANTNS_SUMMARY_2 <- CONTAMINANTS %>% group_by(lab, pool, parasitemia, run, sampleID) %>% summarise(unique_nonref =length(unique(allele)))

parasitemia_plots2 <- ggplot(CONTAMINANTNS_SUMMARY_2, aes(x = pool, y = unique_nonref, fill = parasitemia)) +
  geom_boxplot() +  # Dodged bars for comparison
  facet_wrap(~lab) +  # Facet by lab
  scale_fill_manual(values = c("High" = "orchid", "Low" = "lightblue")) +  # Custom colors
  labs(
    x = "Pool",
    y = "Unique Non-Ref Alleles",
    fill = "Parasitemia",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 1),
    axis.title.x = element_text(size = 14),        
    axis.title.y = element_text(size = 14),   
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  )+
  stat_compare_means(aes(group = parasitemia),
                     method = "kruskal.test",
                     label = "p.signif",
                     label.y = max(CONTAMINANTNS_SUMMARY_2$unique_nonref) * 1.05,
                     size = 5)


parasitemia_plots2

ggsave("parasitemia_plots2.png", parasitemia_plots2, dpi = 300, height = 6, width = 9, bg = "white")


parasitemia_plots3 <- ggplot(CONTAMINANTNS_SUMMARY_2[CONTAMINANTNS_SUMMARY_2$parasitemia == "High",], aes(x = pool, y = unique_nonref, fill = lab)) +
  geom_boxplot() +  # Dodged bars for comparison
  facet_wrap(~parasitemia) +  # Facet by lab
  scale_fill_manual(values = c("ISG" = "orange2", "CISM" = "green3")) +  # Custom colors
  labs(
    x = "Pool",
    y = "Unique Non-Ref Alleles",
    fill = "Lab",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 1),
    axis.title.x = element_text(size = 14),        
    axis.title.y = element_text(size = 14),   
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  )+
  stat_compare_means(aes(group = lab),
                     method = "kruskal.test",
                     label = "p.signif",
                     label.y = max(CONTAMINANTNS_SUMMARY_2[CONTAMINANTNS_SUMMARY_2$parasitemia == "High",]$unique_nonref) * 1.05,
                     size = 5)


parasitemia_plots3

ggsave("parasitemia_plots3.png", parasitemia_plots3, dpi = 300, height = 6, width = 7, bg = "white")



ALLELE_COUNT <-  CONTAMINANTS %>% group_by(lab, pool, parasitemia, allele) %>% summarise(count =length(allele))

length(unique(ALLELE_COUNT$allele))
length(unique(CONTAMINANTS$locus))

ALLELE_COUNT_above1 <- ALLELE_COUNT[ALLELE_COUNT$count > 1,]

ALLELE_COUNT_above1 <- ALLELE_COUNT_above1 %>%
  mutate(allele = fct_reorder(allele, -as.numeric(factor(pool))))

allele_counts_plot <- ggplot(ALLELE_COUNT_above1, aes(x = allele, y = count, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Dodged bars for comparison
  facet_grid(lab ~ parasitemia, scales = "free_y", space="free_y") +  # Facet by lab and pool
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +  # Custom colors
  labs(
    x = "Non-Reference Allele",
    y = "Count",
    fill = "Pool",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),    
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  )+
  coord_flip()

allele_counts_plot

ggsave("allele_counts_plot.png", allele_counts_plot, dpi = 300, height = 5, width = 10, bg = "white")


# ALLELE_COUNT <- table(CONTAMINANTS$allele)
# allele_df <- as.data.frame(ALLELE_COUNT)
# colnames(allele_df) <- c("allele", "count")
# 
# allele_df$allele <- factor(allele_df$allele, levels = allele_df$allele[order(allele_df$count)])
# 
# 
# #put the pool variable
# allele_df <- left_join(allele_df, controls_data[c("allele", "pool")], by = "allele") %>% distinct()
# 
# allele_df <- allele_df %>%
#   group_by(pool) %>%
#   arrange(pool, desc(count), .by_group = TRUE)
# 
# # Create a ggplot histogram
# contam_alleles <- ggplot(
#   ALLELE_COUNT %>% 
#     filter(count > 1) %>% 
#     mutate(allele = reorder(allele, count)),  # Reorder by count
#   aes(x = allele, y = count)
# ) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "white") +
#   labs(
#     x = "Allele",
#     y = "Count",
#     title = ""
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 15, angle = 90, hjust = 1),
#     axis.title.x = element_text(size = 14),        
#     strip.text = element_text(size = 17, face = "bold"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
#   ) +
#   coord_flip() +
#   facet_wrap(~pool, scales = "free", ncol = 1)  # Use "free_y" to allow different y-axis scales per facet
# 
# 
# contam_alleles
# 
# ggsave("contam_alleles2.png", contam_alleles, dpi = 300, height = 14, width = 27, bg = "white")



# save contaminant allele df
CONTAMINANTS_output <- CONTAMINANTS %>%
  separate(sampleID, into = c("sampleID", "run"), sep = "__", extra = "merge") %>%
  select(-allele)

CONTAMINANTS_output <- CONTAMINANTS_output %>% arrange(lab, run, pool, sampleID, Category)

write.csv(CONTAMINANTS_output, "CONTAMINANTS_output.csv", row.names = F)


LOCI_COUNT <-  CONTAMINANTS %>% group_by(lab, pool, parasitemia, locus) %>% summarise(count =length(unique(allele)))

length(unique(LOCI_COUNT$locus))

LOCI_COUNT <- LOCI_COUNT %>%
  mutate(locus = fct_reorder(locus, -as.numeric(factor(pool))))

loci_counts_plot <- ggplot(LOCI_COUNT, aes(x = locus, y = count, fill = pool)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Dodged bars for comparison
  facet_grid(lab ~ parasitemia, scales = "free_y", space="free_y") +  # Facet by lab and pool
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +  # Custom colors
  labs(
    x = "Non-reference allele",
    y = "Count",
    fill = "Pool",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 14),       
    axis.title.y = element_text(size = 14),    
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  )+
  coord_flip()

loci_counts_plot


# LOCI_COUNT <- table(CONTAMINANTS$locus)
# loci_df <- as.data.frame(LOCI_COUNT)
# colnames(loci_df) <- c("locus", "count")
# 
# loci_df$locus <- factor(loci_df$locus, levels = loci_df$locus[order(loci_df$count)])
# 
# loci_df <- left_join(loci_df, controls_data[c("locus", "pool")], by = "locus") %>% distinct()
# 
# # Create a ggplot histogram
# contam_loci <- ggplot( loci_df %>% 
#                          filter(count > 1) %>% 
#                          mutate(locus = reorder(locus, count)),  # Reorder by count
#                        aes(x = locus, y = count)) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "white") +
#   labs(
#     x = "Locus",
#     y = "Non-reference alleles",
#     title = ""
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
#     axis.title.x = element_text(size = 14),        
#     strip.text = element_text(size = 17, face = "bold"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
#   ) +
#   coord_flip()+
#   facet_wrap(~pool, scales = "free", ncol = 2) 
# 
# contam_loci
# 
# ggsave("contam_loci2.png", contam_loci, dpi = 300, height = 15, width = 15, bg = "white")



# *** NUMBER OF CONTAMINANT ALLELES PER CONTROL ***-------------------------

n_contams_per_run <- CONTAMINANTS %>%
  group_by(run, sampleID, lab, pool, parasitemia) %>%
  summarise(n_contams = length(unique(allele)))

n_contams_per_run

# n_contams_per_run <- n_contams_per_run %>%
#   mutate(sampleID = factor(sampleID, levels = sampleID[order(-n_contams)]))

set.seed(42069)

num_colors <- length(unique(n_contams_per_run$run))
color_palette <- rgb(runif(num_colors), runif(num_colors), runif(num_colors))

#remove run name from sample id. some controls aare gonna appear stacked because they are named the same, it`s fine
n_contams_per_run$sampleID_plot <- gsub("__.*", "", n_contams_per_run$sampleID)

#order for plot
n_contams_per_run <- n_contams_per_run %>%
  group_by(run) %>%
  arrange(run, desc(n_contams), sampleID_plot) %>%
  mutate(sampleID_plot = factor(sampleID_plot, levels = unique(sampleID_plot))) %>%
  ungroup()

n_contams_per_run <- n_contams_per_run %>%
  mutate(sampleID = fct_reorder(sampleID, -as.numeric(factor(run))))


#plot
contams1  <- ggplot(n_contams_per_run,
                    aes(x = sampleID, y = n_contams, fill = pool)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(lab ~ parasitemia, scales = "free", space="free") +  # Facet by lab and pool
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +  # Custom colors
  labs(
    x = "Controls",
    y = "Non-reference alleles",
    fill = "Pool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),        
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  )+
  guides(fill = guide_legend(ncol = 1))+
  coord_flip()


contams1

ggsave("contams12.png", contams1, dpi = 300, height = 12, width = 15, bg = "white")



##### ?) JACCARD'S DISTANCE ###########33-------------------------

# 1. Get the sampleIDs with many contaminants
SAMPLES_WITH_MANY_CONTAMS <- n_contams_per_run %>%
  filter(n_contams > 0) %>%
  pull(sampleID)

# 2. For each pool (within the selected samples), calculate the Jaccard distance
dist_df_all <- CONTAMINANTS %>%
  filter(sampleID %in% SAMPLES_WITH_MANY_CONTAMS) %>%
  group_by(pool) %>%
  do({
    pool_data <- .
    # 3. Create the binary allele matrix for this pool:
    allele_matrix <- pool_data %>%
      select(sampleID, locus) %>%
      distinct() %>%
      mutate(present = 1) %>%
      spread(key = locus, value = present, fill = 0)
    
    # Set sampleID as row names (assumed unique within a pool)
    allele_matrix <- allele_matrix %>% column_to_rownames("sampleID")
    
    # 4. Calculate the Jaccard distance for this pool
    dist_matrix <- vegdist(allele_matrix, method = "jaccard")
    
    # 5. Convert the distance matrix to a long data frame
    df <- as.data.frame(as.matrix(dist_matrix))
    df$SampleID1 <- rownames(df)
    df_long <- melt(df, id.vars = "SampleID1")
    df_long <- df_long %>%
      rename(SampleID2 = variable, distance = value) %>%
      mutate(
        SampleID1 = factor(SampleID1, levels = rownames(as.matrix(dist_matrix))),
        SampleID2 = factor(SampleID2, levels = rownames(as.matrix(dist_matrix)))
      )
    df_long
  }) %>%
  ungroup()

dist_df_all <- dist_df_all %>%
  group_by(pool) %>%
  mutate(n_samples = n_distinct(as.character(SampleID1))) %>% 
  ungroup() %>%
  # Remove self-comparison rows only when there's exactly one sample in that pool
  filter(!(n_samples == 1 & SampleID1 == SampleID2)) %>%
  select(-n_samples)

# 6. Plot the heatmap faceted by pool
hist_jaccard_all <- ggplot(dist_df_all, aes(x = SampleID1, y = SampleID2, fill = distance)) +
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
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.title.x = element_text(size = 14),        
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  ) +
  facet_wrap(~ pool, scales = "free", ncol = 2)

# Display the plot
hist_jaccard_all

ggsave("hist_jaccard_all3.png", hist_jaccard_all, dpi = 300, height = 20, width = 20, bg = "white")




##
# 1. Get the sampleIDs with many contaminants
SAMPLES_WITH_MANY_CONTAMS <- n_contams_per_run %>%
  filter(n_contams > 1) %>%
  pull(sampleID)

# 2. For each pool (within the selected samples), calculate the Jaccard distance
dist_df_all <- CONTAMINANTS %>%
  filter(sampleID %in% SAMPLES_WITH_MANY_CONTAMS) %>%
  group_by(pool) %>%
  do({
    pool_data <- .
    # 3. Create the binary allele matrix for this pool:
    allele_matrix <- pool_data %>%
      select(sampleID, locus) %>%
      distinct() %>%
      mutate(present = 1) %>%
      spread(key = locus, value = present, fill = 0)
    
    # Set sampleID as row names (assumed unique within a pool)
    allele_matrix <- allele_matrix %>% column_to_rownames("sampleID")
    
    # 4. Calculate the Jaccard distance for this pool
    dist_matrix <- vegdist(allele_matrix, method = "jaccard")
    
    # 5. Convert the distance matrix to a long data frame
    df <- as.data.frame(as.matrix(dist_matrix))
    df$SampleID1 <- rownames(df)
    df_long <- melt(df, id.vars = "SampleID1")
    df_long <- df_long %>%
      rename(SampleID2 = variable, distance = value) %>%
      mutate(
        SampleID1 = factor(SampleID1, levels = rownames(as.matrix(dist_matrix))),
        SampleID2 = factor(SampleID2, levels = rownames(as.matrix(dist_matrix)))
      )
    df_long
  }) %>%
  ungroup()

dist_df_all <- dist_df_all %>%
  group_by(pool) %>%
  mutate(n_samples = n_distinct(as.character(SampleID1))) %>% 
  ungroup() %>%
  # Remove self-comparison rows only when there's exactly one sample in that pool
  filter(!(n_samples == 1 & SampleID1 == SampleID2)) %>%
  select(-n_samples)

# 6. Plot the heatmap faceted by pool
hist_jaccard_10plus <- ggplot(dist_df_all, aes(x = SampleID1, y = SampleID2, fill = distance)) +
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
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.title.x = element_text(size = 14),        
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  ) +
  facet_wrap(~ pool, scales = "free", ncol = 2)

# Display the plot
hist_jaccard_10plus

ggsave("hist_jaccard_10plus3.png", hist_jaccard_10plus, dpi = 300, height = 20, width = 20, bg = "white")





# *** FREQUENCY DISTRIBUTION OF CONTAMINANT ALLELES *** -------------------

set.seed(420)

num_colors <- length(unique(CONTAMINANTS$sampleID))

color_palette <- rgb(runif(num_colors), runif(num_colors), runif(num_colors))

# ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = sampleID)) +
#   geom_histogram(bins = 100, alpha = 0.5, position = "identity", color = "black") +
#   scale_fill_manual(values = color_palette) +
#   labs(
#     title = "",
#     x = "In-sample allele frequency",
#     y = "Count"
#   ) +
#   theme_minimal() +
#   theme(legend.title = element_blank())+
#   guides(fill = guide_legend(ncol = 1))

contams2 <- ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = pool)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
  facet_grid(lab ~ parasitemia, scales = "free", space="free") +  # Facet by lab and pool
  scale_fill_manual(values = c("1A" = "red", "1B" = "blue", "2" = "green")) +  # Custom colors
  labs(
    title = "",
    x = "In-sample allele frequency",
    y = "Non-reference alleles",
    fill = "Pool"
  ) +
  theme_minimal() +
  #theme(legend.title = element_blank())+
  #guides(fill = guide_legend(ncol = 1))+
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),
    axis.title.x = element_text(size = 14),        
    axis.title.y = element_text(size = 14),     
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  ) 

contams2

ggsave("contams23.png", contams2, dpi = 300, height = 10, width = 17, bg = "white")



# contams2_pool <- ggplot(CONTAMINANTS, aes(x = norm.reads.locus, fill = sampleID)) +
#   geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
#   scale_fill_manual(values = color_palette) +
#   labs(
#     title = "",
#     x = "In-sample allele frequency",
#     y = "Non-reference alleles"
#   ) +
#   theme_minimal() +
#   theme(legend.title = element_blank())+
#   guides(fill = guide_legend(ncol = 1))+
#   theme(
#     axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
#     axis.title.x = element_text(size = 10),        
#     strip.text = element_text(size = 9, face = "bold"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     legend.position = "none"
#   ) +
#   facet_wrap(~pool, scales = "free_y", ncol = 2)
# 
# contams2_pool
# 
# ggsave("contams22_pool.png", contams2_pool, dpi = 300, height = 7, width = 7, bg = "white")


# controls with contaminant alleles that have a freq = 1
CONTAMINATED_CONTROLS_FREQ1 <- length(unique(CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]$sampleID))
RUNS_WITH_CONTAMINATED_CONTROLS_FREQ1  <- length(unique(CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]$run))

print(paste("There were", CONTAMINATED_CONTROLS_FREQ1, " controls with at least 1 fixed contaminant allele (norm.reads.locus = 1) across", RUNS_WITH_CONTAMINATED_CONTROLS_FREQ1, "runs"))



# *** CONTAMINANT ALLELES IN FIELD SAMPLES FROM EACH RUN?  *** -------------------

# remove controls from original data AND keep runs from which the controls used in the study are taken (for instance, nanna's run is not used because it doesn't distinguish controls, etc.)
field_data <- merged_dfs_agg[!merged_dfs_agg$sampleID %in% CONTAMINANTS$sampleID,]
field_data <- field_data[field_data$run %in% unique(controls_data$run),]

SAMPLES <- unique(CONTAMINANTS$sampleID)


## A) for all pools together

contam_procedence_results <- list()

for (sample in SAMPLES) {
  
  subsample_CONTAMINANTS_data <- CONTAMINANTS[CONTAMINANTS$sampleID == sample,]
  
  run <- strsplit(sample, "__")[[1]][2]
  
  subsample_filed_data <- field_data[field_data$run == run,]
  
  #contaminants found in field samples from same run 
  contaminants_in_field <- intersect(subsample_CONTAMINANTS_data$allele, subsample_filed_data$allele)
  
  #store the contaminants not found in field samples from same run 
  contaminants_not_in_field <- setdiff(subsample_CONTAMINANTS_data$allele, subsample_filed_data$allele)
  
  
  # # of cross contamination contaminants (contaminants not found in same run but found in others) and thee runs they belong to
  n_cross_comtams <- intersect(contaminants_not_in_field, field_data$allele)
  cross_contam_runs <- field_data[field_data$allele %in% contaminants_not_in_field,]$run
  
  # contaminants not in same and other runs. oriigin unknown
  n_contams_unknown <- setdiff(contaminants_not_in_field, field_data$allele)
  
  
  contam_procedence_results[[sample]] <- data.frame(
    run = run,
    sampleID = sample,
    n_nonref_in_control = length(unique(subsample_CONTAMINANTS_data$allele)),
    n_nonref_in_field = length(unique(contaminants_in_field)),
    nonref_in_field = I(list(contaminants_in_field)),
    n_nonref_not_in_field = length(unique(contaminants_not_in_field)),
    nonref_not_in_field = I(list(contaminants_not_in_field)),
    n_cross_nonref = length(unique(n_cross_comtams)),
    cross_nonref_runs = I(list(unique(cross_contam_runs))),
    n_unknown_nonref = length(unique(n_contams_unknown))  
  )
}

contam_procedence_results <- do.call(rbind, contam_procedence_results)

# add parasitemia classification
contam_procedence_results <- left_join(contam_procedence_results, unique(CONTAMINANTS[c("sampleID", "parasitemia", "lab")]), by= "sampleID")

contam_procedence_results$sampleID <- sub("__.*", "", contam_procedence_results$sampleID)
rownames(contam_procedence_results) <- NULL

# Ensure list-columns (e.g., vector columns) are stored as character strings
contam_procedence_results <- contam_procedence_results %>%
  mutate(across(where(is.list), ~sapply(., toString)))  # Converts list elements to comma-separated strings

# Replace blanks (empty values) with NA for proper CSV formatting
contam_procedence_results[contam_procedence_results == ""] <- NA

contam_procedence_results <- contam_procedence_results %>% 
  select(run, sampleID, lab, parasitemia, n_nonref_in_control, n_nonref_in_field, n_cross_nonref, cross_nonref_runs, n_unknown_nonref) %>%
  arrange(lab, parasitemia, run)

# Write to CSV
write.csv(contam_procedence_results, "contam_procedence_results2_all_pools_together.csv", row.names = FALSE, na = "")



# ### B) for each pool separately
# 
# # Create an empty list to store results from each pool
# all_results <- list()
# 
# unique_pools <- unique(CONTAMINANTS$pool)
# 
# for(pool in unique_pools) {
#   
#   CONTAMINANTS_pool <- CONTAMINANTS[CONTAMINANTS$pool == pool,]
#   field_data_pool <- field_data[field_data$pool == pool,]
#   
#   contam_procedence <- list()
#   
#   for (sample in SAMPLES) {
#     
#     subsample_CONTAMINANTS_data <- CONTAMINANTS_pool[CONTAMINANTS_pool$sampleID == sample,]
#     
#     run <- strsplit(sample, "__")[[1]][2]
#     
#     subsample_filed_data <- field_data_pool[field_data_pool$run == run,]
#     
#     #contaminants found in field samples from same run 
#     contaminants_in_field <- intersect(subsample_CONTAMINANTS_data$allele, subsample_filed_data$allele)
#     
#     #store the contaminants not found in field samples from same run 
#     contaminants_not_in_field <- setdiff(subsample_CONTAMINANTS_data$allele, subsample_filed_data$allele)
#     
#     # # of cross contamination contaminants (contaminants not found in same run but found in others) and thee runs they belong to
#     n_cross_comtams <- intersect(contaminants_not_in_field, field_data$allele)
#     cross_contam_runs <- field_data[field_data$allele %in% contaminants_not_in_field,]$run
#     
#     # contaminants not in same and other runs. oriigin unknown
#     n_contams_unknown <- setdiff(contaminants_not_in_field, field_data$allele)
#     
#     
#     contam_procedence[[sample]] <- data.frame(
#       run = run,
#       sampleID = sample,
#       n_nonref_in_control = length(unique(subsample_CONTAMINANTS_data$allele)),
#       n_nonref_in_field = length(unique(contaminants_in_field)),
#       nonref_in_field = I(list(contaminants_in_field)),
#       n_nonref_not_in_field = length(unique(contaminants_not_in_field)),
#       nonref_not_in_field = I(list(contaminants_not_in_field)),
#       n_cross_nonref = length(unique(n_cross_comtams)),
#       cross_nonref_runs = I(list(unique(cross_contam_runs))),
#       n_unknown_nonref = length(unique(n_contams_unknown))  
#     )
#   }
#   
#   # Combine the results for the current pool into one data frame
#   pool_results <- do.call(rbind, contam_procedence)
#   
#   # Remove the part after "__" from sampleID if desired
#   pool_results$sampleID <- sub("__.*", "", pool_results$sampleID)
#   rownames(pool_results) <- NULL
#   
#   # Convert list-columns to comma-separated strings for CSV output
#   pool_results <- pool_results %>%
#     mutate(across(where(is.list), ~ sapply(., toString)))
#   
#   # Replace blank strings with NA for proper CSV formatting
#   pool_results[pool_results == ""] <- NA
#   
#   # Select columns in desired order (optional)
#   pool_results <- pool_results %>%
#     select(run, sampleID, n_nonref_in_control, n_nonref_in_field, n_cross_nonref, cross_nonref_runs, n_unknown_nonref)
#   
#   # Add a new column 'pool' to indicate the current pool
#   pool_results$pool <- pool
#   
#   # Store the pool results in the overall list
#   all_results[[pool]] <- pool_results
#   
# }
# 
# # Combine results from all pools into a single data frame
# contam_procedence_results_POOLS <- do.call(rbind, all_results)
#                                            
# # Write final results to a CSV file
# write.csv(contam_procedence_results_POOLS, "contam_procedence_results2_separate_pools.csv", row.names = FALSE, na = "")
# 
# 

# ggplot(contam_procedence_results, aes(x = percentage_contams_in_field_samples_from_run)) +
#   geom_histogram(bins = 50, color = "black", fill = "blue", alpha = 0.7) +
#   labs(
#     title = "",
#     x = "Percentage of Contaminants in Field Samples",
#     y = "# Controls"
#   ) +
#   theme_minimal()


# # Calculate percentages based on n_nonref_in_control (100%)
# contam_procedence_results_percentages <- contam_procedence_results %>%
#   mutate(
#     perc_nonref_in_field = (n_nonref_in_field / n_nonref_in_control) * 100,
#     perc_cross_nonref = (n_cross_nonref / n_nonref_in_control) * 100,
#     perc_unknown_nonref = (n_unknown_nonref / n_nonref_in_control) * 100
#   )
# 
# 
# # Reshape data to long format for ggplot
# contam_procedence_long <- contam_procedence_results_percentages %>%
#   select(sampleID, run, perc_nonref_in_field, perc_cross_nonref, perc_unknown_nonref) %>%
#   pivot_longer(cols = c(perc_nonref_in_field, perc_cross_nonref, perc_unknown_nonref),
#                names_to = "nonref_type",
#                values_to = "count")
# 
# contam_procedence_long$nonref_type <- gsub("perc_", "",contam_procedence_long$nonref_type)


# Reshape data to long format for ggplot
contam_procedence_long <- contam_procedence_results %>%
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
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
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

ggsave("contams3210.png", contams3, dpi = 300, height = 12, width = 15, bg = "white")


# # per pool
# 
# # Calculate percentages based on n_nonref_in_control (100%)
# contam_procedence_results_percentages <- contam_procedence_results_POOLS %>%
#   mutate(
#     perc_nonref_in_field = (n_nonref_in_field / n_nonref_in_control) * 100,
#     perc_cross_nonref = (n_cross_nonref / n_nonref_in_control) * 100,
#     perc_unknown_nonref = (n_unknown_nonref / n_nonref_in_control) * 100
#   )
# 
# # contam_procedence_results_percentages$sampleID <- paste0(contam_procedence_results_percentages$sampleID, "_", seq_len(nrow(contam_procedence_results_percentages)))
# 
# # Reshape data to long format for ggplot
# contam_procedence_long <- contam_procedence_results_percentages %>%
#   select(sampleID, run, perc_nonref_in_field, perc_cross_nonref, perc_unknown_nonref, pool) %>%
#   pivot_longer(cols = c(perc_nonref_in_field, perc_cross_nonref, perc_unknown_nonref),
#                names_to = "nonref_type",
#                values_to = "count")
# 
# contam_procedence_long$nonref_type <- gsub("perc_", "",contam_procedence_long$nonref_type)
# 
# contams3_pools <- ggplot(contam_procedence_long, aes(x = sampleID, y = count, fill = nonref_type)) +
#   geom_bar(stat = "identity", position = "stack") +  # Stacked bars
#   facet_grid(run ~ pool, scales = "free", space = "free") +  # Facet by both 'pool' and 'run'
#   scale_fill_manual(values = c("#e31a1c", "#1f78b4",  "black")) +  # Custom colors
#   labs(
#     x = "Control",
#     y = "% Non-reference alleles",
#     fill = "Source",
#     title = ""
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     strip.text = element_text(size = 8, face = "bold"),  # Format facet labels
#     strip.text.y = element_text(angle = 0, hjust = 0)
#     #panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
#   )
# 
# contams3_pools
# 
# ggsave("contams32_pools.png", contams3_pools, dpi = 300, height = 14, width = 20, bg = "white")

# # *** MISSING (REF) ALLELES IN 3D7 CONTROLS  *** -------------------
# 
# all_ref_alleles <- unique(controls_data[grepl("_.$", controls_data$allele),]$locus)
# n_expected_alleles <- length(all_ref_alleles) # number of expected alleles for 3d7 controls  == number of loci
# 
# allele_counts <- controls_data %>%
#   filter(grepl("_.$", allele)) %>%
#   group_by(sampleID) %>%
#   summarise(n_correctly_sequenced_ref_alleles = length(unique(allele)))
# 
# allele_counts$missing_alleles <-  n_expected_alleles - allele_counts$n_correctly_sequenced_ref_alleles
# 
# allele_counts <- allele_counts %>%
#   separate(sampleID, into = c("sampleID", "run"), sep = "__")
# 
# 
# missing <- ggplot(allele_counts, aes(x = sampleID, y = missing_alleles, fill = run)) +
#   geom_bar(stat = "identity", color = "black") +
#   scale_fill_manual(values = color_palette) +
#   labs(
#     title = "",
#     x = "3D7 Controls",
#     y = "# Missing Alleles",
#     fill = "Run"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.title = element_blank()
#   )+
#   facet_wrap(~run, scales = "free_x")+
#   guides(fill = FALSE)
# 
# missing
# 
# ggsave("missing.png", missing, dpi = 300, height = 15, width = 15, bg = "white")



### PROPOSED THRESHOLDS------------

#remove alleles with freq of 1 (possible mislabelling of controls)
CONTAMINANTS_less_than_1 <- CONTAMINANTS[CONTAMINANTS$norm.reads.locus< 1,]


# all pools together

## interpret as for 50%, 50% of contaminants appear at a freq of n;
# contam_thresholds <- as.data.frame(t(t(quantile(CONTAMINANTS_less_than_1$norm.reads.locus, probs= c(0.5, 0.75, 0.9, 0.95, 0.99)))))
# 
# colnames(contam_thresholds)[1] <- "MAF_threshold"
# 
# contam_thresholds
# 
# write.csv(contam_thresholds, "comtam_thresholds2.csv")
# 
# 
# thresholds <- ggplot(CONTAMINANTS_less_than_1, aes(x = norm.reads.locus)) +
#   geom_histogram(bins = 100, alpha = 0.5, position = "identity", fill = "steelblue", color = "white") +
#   geom_vline(data = as.data.frame(contam_thresholds), 
#              aes(xintercept = MAF_threshold, color = factor(rownames(contam_thresholds))), 
#              linetype = "solid", size = 1) +
#   scale_color_manual(name = "Thresholds", values = c("red", "blue", "green", "purple", "orange")) +
#   labs(
#     title = "",
#     x = "In-sample allele frequency",
#     y = "Non-reference alleles"
#   ) +
#   theme_minimal() +
#   guides(color = guide_legend(title = "% Non-ref alleles\n    eliminated")) +
#   theme(
#     legend.position = "right",
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 8)
#   )
# 
# thresholds
# 
# ggsave("thresholds2.png", thresholds, height = 5, width = 8, dpi = 300, bg = "white")


## by pool
contam_thresholds_pools <- CONTAMINANTS_less_than_1 %>%
  group_by(lab, parasitemia, pool) %>%
  summarise(
    `50%` = quantile(norm.reads.locus, probs = 0.5, na.rm = TRUE),
    `75%` = quantile(norm.reads.locus, probs = 0.75, na.rm = TRUE),
    `90%` = quantile(norm.reads.locus, probs = 0.9, na.rm = TRUE),
    `95%` = quantile(norm.reads.locus, probs = 0.95, na.rm = TRUE),
    `99%` = quantile(norm.reads.locus, probs = 0.99, na.rm = TRUE)
  )

contam_thresholds_pools

contam_thresholds_pools_OUT <- as.data.frame(t(contam_thresholds_pools))

colnames(contam_thresholds_pools_OUT) <- contam_thresholds_pools_OUT[1,]
contam_thresholds_pools_OUT <- contam_thresholds_pools_OUT[-1,]

contam_thresholds_pools_OUT

write.csv(contam_thresholds_pools_OUT, "comtam_thresholds2_pools.csv", row.names = T)


# Convert contam_thresholds_pools to long format for better ggplot mapping
contam_thresholds_long <- contam_thresholds_pools %>%
  pivot_longer(cols = `50%`:`99%`, names_to = "percentile", values_to = "MAF_threshold")

# Plot faceted by pool
thresholds <- ggplot(CONTAMINANTS_less_than_1, aes(x = norm.reads.locus)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity", fill = "steelblue", color = "white") +
  geom_vline(data = contam_thresholds_long, 
             aes(xintercept = MAF_threshold, color = percentile), 
             linetype = "solid", size = 1) +
  facet_grid(lab~parasitemia) +  # Facet by pool
  scale_color_manual(name = "Thresholds", values = c("red", "blue", "green", "purple", "orange")) +
  labs(
    title = "",
    x = "In-sample allele frequency",
    y = "Non-reference alleles"
  ) +
  theme_minimal() +
  
  guides(color = guide_legend(title = "% Non-ref alleles\n    eliminated")) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),
    axis.title.x = element_text(size = 14),        
    axis.title.y = element_text(size = 14),     
    strip.text = element_text(size = 17, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.text = element_text(size = 12),   # Increases legend text size
    legend.title = element_text(size = 14) 
  )

thresholds

ggsave("thresholds23_pools.png", thresholds, height = 8, width = 12, dpi = 300, bg = "white")




### EXPLORE FIXED CONTAMINANTS --------
fixed_contams <- CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]

fixed_contams$allele <-  gsub(".*__", "", fixed_contams$allele)
fixed_contams$sampleID <-  gsub("__.*$", "", fixed_contams$sampleID)

fixed_contams_summary <- fixed_contams %>% 
  group_by(sampleID, run, lab, parasitemia, pool, locus) %>% 
  summarise(
    #n_nonref_alleles_fixed = n_distinct(allele),  # Count unique alleles
    unique_alleles = toString(unique(allele)),    # Convert list to a comma-separated string
    .groups = "drop"  # Ungroup after summarization
  ) %>%
  arrange(lab, locus, unique_alleles, run, pool, parasitemia)

fixed_contams_summary


write.csv(fixed_contams_summary, "fixed_contams_summary.csv", row.names = F)

#asv for andres
fixed_contams <- CONTAMINANTS[CONTAMINANTS$norm.reads.locus == 1,]
fixed_contams_asv<- merge(fixed_contams, merged_dfs[c("sampleID", "locus", "run", "Category", "reads", "pseudo_cigar","norm.reads.locus", "asv")], by = c("sampleID", "locus", "run", "Category", "reads", "pseudo_cigar","norm.reads.locus"))

fixed_contams_asv <- fixed_contams_asv %>% arrange(lab, locus, pseudo_cigar, run, pool, parasitemia)

write.csv(fixed_contams_asv, "fixed_contams_summary_ASV.csv", row.names = F)


#align the sequences from amplicon Pf3D7_07_v3-403801-404052-1B
library(Biostrings)
library(DECIPHER)

amp <- fixed_contams_asv[fixed_contams_asv$locus == "Pf3D7_07_v3-403801-404052-1B",]

seqs <- DNAStringSet(amp$asv)
names(seqs) <- paste(amp$sampleID, amp$allele)

# Align sequences
aligned_seqs <- AlignSeqs(seqs)

# View the aligned sequences
BrowseSeqs(aligned_seqs)


# # missing alleles
# ggplot(allele_counts, aes(x = missing_alleles)) +
#   geom_histogram(
#     binwidth = 1,  # Adjust binwidth for your data
#     fill = "#69b3a2",  # Custom fill color
#     color = "black",   # Border color
#     alpha = 0.8        # Transparency
#   ) +
#   labs(
#     title = "",
#     x = "Contaminant Alleles",
#     y = "Frequency"
#   ) +
#   theme_minimal(base_size = 14) +  # Clean minimal theme
#   theme(
#     plot.title = element_text(hjust = 0.5),  # Center and bold title
#     axis.title = element_text()
#   )
# 
# #interpret backwards (when 0.1, 90% of samples have n missing alleles)
# quantile(allele_counts$missing_alleles, probs= c(0.1, 0.25, 0.5, 0.75, 0.9))

