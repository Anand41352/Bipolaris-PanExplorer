#!/usr/bin/env Rscript
# ============================================================
# Pangenome visualization and rarefaction analysis
# Dataset: OrthoFinder Results_BS19_v1_fastree_FINAL
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  stop("Install UpSetR before running this script.")
}

# Run from the directory containing Results_BS19_v1_fastree_FINAL,
# or edit results_dir to the appropriate path.
results_dir <- "Results_BS19_v1_fastree_FINAL"
orthogroup_file <- file.path(results_dir, "Orthogroups",
                             "Orthogroups.GeneCount.tsv")
unassigned_file <- file.path(results_dir, "Orthogroups",
                             "Orthogroups_UnassignedGenes.tsv")
species_tree_file <- file.path(results_dir, "Species_Tree",
                               "SpeciesTree_rooted_node_labels.txt")

prefix <- "BS19"
softcore_min <- 17
shell_min <- 2
Nperm <- 200
set.seed(1234)

if (!file.exists(orthogroup_file)) stop("Missing: ", orthogroup_file)

# Load orthogroup gene counts and remove the summary Total column.
df <- read.delim(orthogroup_file, row.names = 1,
                 check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(df) > 0 && grepl("^Total$", colnames(df)[ncol(df]))) {
  df <- df[, -ncol(df), drop = FALSE]
}
df[] <- lapply(df, function(x) as.integer(x > 0))

# Add unassigned genes as individual presence/absence rows.
if (file.exists(unassigned_file)) {
  df_uniq <- read.delim(unassigned_file, row.names = 1,
                        check.names = FALSE, stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
  missing_cols <- setdiff(colnames(df), colnames(df_uniq))
  for (nm in missing_cols) df_uniq[[nm]] <- NA
  df_uniq <- df_uniq[, colnames(df), drop = FALSE]
  df_uniq[] <- lapply(df_uniq, function(x) {
    as.integer(!is.na(x) & trimws(as.character(x)) != "")
  })
  df_combined <- rbind(df, df_uniq)
} else {
  warning("Unassigned-gene file not found; assigned orthogroups only.")
  df_combined <- df
}

df_combined[is.na(df_combined)] <- 0
df_combined <- as.data.frame(df_combined, check.names = FALSE)
n_genomes <- ncol(df_combined)
if (n_genomes != 19) warning("Expected 19 genomes; found ", n_genomes, ".")

# Categorize orthogroups: core, soft-core, shell and cloud.
presence_counts <- rowSums(df_combined != 0)
category <- case_when(
  presence_counts == n_genomes ~ "Core",
  presence_counts >= softcore_min ~ "Softcore",
  presence_counts >= shell_min &
    presence_counts <= (n_genomes - 3) ~ "Shell",
  presence_counts == 1 ~ "Cloud",
  TRUE ~ "Shell"
)
cat_levels <- c("Core", "Softcore", "Shell", "Cloud")
category <- factor(category, levels = cat_levels)

og_info <- tibble(
  Orthogroup = rownames(df_combined),
  N_present = presence_counts,
  Category = category
)

category_colours <- c(
  Core = "#F8766D", Softcore = "#7CAE00",
  Shell = "#00BFC4", Cloud = "#C77CFF"
)

# A. UpSet plot: genome-level intersections of orthogroups.
pdf(paste0(prefix, "_FigA_UpSet_Genome_Intersections.pdf"),
    width = 12, height = 7)
UpSetR::upset(
  as.data.frame(df_combined),
  nsets = ncol(df_combined),
  order.by = "freq",
  mainbar.y.label = "Number of orthogroups",
  sets.x.label = "Orthogroups per genome"
)
dev.off()

# B. Composition of the super-pangenome.
pie_df <- og_info %>%
  count(Category, .drop = FALSE) %>%
  mutate(Percent = round(n / sum(n) * 100, 2))

pdf(paste0(prefix, "_FigB_Pie_SuperPangenome.pdf"),
    width = 6, height = 6)
print(
  ggplot(pie_df, aes(x = "", y = n, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") + theme_void() +
    scale_fill_manual(values = category_colours, drop = FALSE) +
    geom_text(aes(label = paste0(Percent, "%")),
              position = position_stack(vjust = 0.5), size = 4) +
    ggtitle("Composition of the Super-pangenome")
)
dev.off()

# C. Number of orthogroups in each category per genome.
df_long <- df_combined %>%
  rownames_to_column("Orthogroup") %>%
  pivot_longer(-Orthogroup, names_to = "Genome",
               values_to = "Presence") %>%
  left_join(og_info, by = "Orthogroup") %>%
  filter(Presence == 1)

per_genome_counts <- df_long %>%
  group_by(Genome, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  complete(Genome, Category = factor(cat_levels, levels = cat_levels),
           fill = list(Count = 0))

write.csv(per_genome_counts,
          paste0(prefix, "_Table_PerGenome_Category_Counts.csv"),
          row.names = FALSE)

pdf(paste0(prefix, "_FigC_PerGenome_Category_Counts.pdf"),
    width = 11, height = 6)
print(
  ggplot(per_genome_counts, aes(x = Genome, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = category_colours, drop = FALSE) +
    theme_minimal() + xlab("") + ylab("Number of orthogroups") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Core, soft-core, shell and cloud orthogroups per genome")
)
dev.off()

# D. Export orthogroup lists by category.
dir.create("Category_Exports", showWarnings = FALSE)
for (cat_name in cat_levels) {
  writeLines(
    og_info$Orthogroup[og_info$Category == cat_name],
    file.path("Category_Exports", paste0(cat_name, "_OGs.txt"))
  )
}
write.csv(og_info, paste0(prefix, "_Table_Orthogroup_Categories.csv"),
          row.names = FALSE)

# E. Rarefaction: 200 random genome permutations.
genomes <- colnames(df_combined)

pan_core_one_order <- function(order_of_genomes) {
  seen <- rep(FALSE, nrow(df_combined))
  tibble(
    k = seq_along(order_of_genomes),
    Pan = map_int(seq_along(order_of_genomes), function(i) {
      seen <<- seen | (df_combined[, order_of_genomes[i]] == 1)
      sum(seen)
    }),
    Core = map_int(seq_along(order_of_genomes), function(i) {
      sum(rowSums(df_combined[, order_of_genomes[seq_len(i)],
                              drop = FALSE]) == i)
    })
  )
}

sim_results <- map_dfr(seq_len(Nperm), function(p) {
  pc <- pan_core_one_order(sample(genomes))
  pc$perm <- p
  pc
})

write.csv(sim_results, paste0(prefix, "_Table_PanCore_Rarefaction.csv"),
          row.names = FALSE)

sim_summary <- sim_results %>%
  group_by(k) %>%
  summarise(
    Pan_mean = mean(Pan),
    Pan_lower = quantile(Pan, 0.025),
    Pan_upper = quantile(Pan, 0.975),
    Core_mean = mean(Core),
    Core_lower = quantile(Core, 0.025),
    Core_upper = quantile(Core, 0.975),
    .groups = "drop"
  )
write.csv(sim_summary,
          paste0(prefix, "_Table_PanCore_Rarefaction_Summary.csv"),
          row.names = FALSE)

pdf(paste0(prefix, "_FigE_Rarefaction_PanCore.pdf"),
    width = 7, height = 5)
print(
  ggplot(sim_results, aes(x = factor(k))) +
    geom_boxplot(aes(y = Pan), color = "red", fill = "pink",
                 alpha = 0.6, outlier.size = 0.5) +
    geom_boxplot(aes(y = Core), color = "blue", fill = "lightblue",
                 alpha = 0.6, outlier.size = 0.5) +
    xlab("Number of genomes sampled") +
    ylab("Number of orthogroups") +
    ggtitle("Rarefaction of pan- and core-orthogroups") +
    theme_minimal()
)
dev.off()

# F. Model fitting.
# The model-fitting code was not present in the supplied script:
# it was replaced by the note "(same as your version, left unchanged
# for brevity)". Insert the original Heaps' law and exponential-decay
# fitting code here; no equations or fitting choices are invented.

# G. Presence/absence heatmap.
binary_matrix <- as.matrix(df_combined)
storage.mode(binary_matrix) <- "numeric"
row_order <- order(rowMeans(binary_matrix), decreasing = TRUE)
binary_matrix <- binary_matrix[row_order, , drop = FALSE]

ann_row <- data.frame(
  Category = og_info$Category[
    match(rownames(binary_matrix), og_info$Orthogroup)
  ]
)
rownames(ann_row) <- rownames(binary_matrix)

# Read and validate the species tree. pheatmap itself does not draw a
# phylogenetic tree; this file is checked here but not overlaid on heatmap.
if (file.exists(species_tree_file) && requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(species_tree_file)
  missing_tips <- setdiff(colnames(binary_matrix), tree$tip.label)
  if (length(missing_tips) > 0) {
    warning("Species tree is missing matrix labels: ",
            paste(missing_tips, collapse = ", "))
  }
} else {
  warning("Species tree not read (file missing or package ape unavailable).")
}

pdf(paste0(prefix, "_FigG_Heatmap_PresenceAbsence.pdf"),
    width = 9, height = 10)
pheatmap(
  binary_matrix,
  color = c("#3288bd", "#9e0142"),
  breaks = c(-0.01, 0.5, 1.01),
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 90,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  annotation_colors = list(Category = category_colours),
  main = paste0("Gene presence/absence (", ncol(binary_matrix), " genomes)")
)
dev.off()

message("Pangenome analysis plots and tables completed.")
