setwd("/home/naren/Desktop/ANAND/PAN_GENOME/ANNOTATION/GO_ENRICHMENT_FILE_PREPARATION/final_go_sheet")

# ============================================================================
# GO Term Annotation Statistics Analysis for Bipolaris Pangenome
# Fixed Version with Correct Percentages and No Grids
# ============================================================================

# Load required libraries
library(tidyverse)
library(RColorBrewer)
library(ggVennDiagram)
library(patchwork)
library(scales)
library(ggpubr)
library(UpSetR)
library(viridis)
library(grid)  # For grid.text function
library(magick)  # For image merging

# Set working directory
setwd("~/Desktop/ANAND/PAN_GENOME/ANNOTATION/GO_ENRICHMENT_FILE_PREPARATION/final_go_sheet")

# ============================================================================
# 1. DEFINE TOTAL GENE COUNTS PER CLASS
# ============================================================================

# Total genes per pangenome class (including those without GO annotations)
TOTAL_GENES <- data.frame(
  Gene_Class = c("CORE", "SOFTCORE", "SHELL", "CLOUD"),
  Total_Genes = c(10323, 1444, 3350, 1864)
)

TOTAL_ALL_GENES <- sum(TOTAL_GENES$Total_Genes)  # 16981

# ============================================================================
# 2. DATA LOADING AND PREPROCESSING
# ============================================================================

# Read merged data (with all unique GO terms per gene)
merged_go <- read.csv("/home/naren/Desktop/ANAND/PAN_GENOME/ANNOTATION/GO_ENRICHMENT_FILE_PREPARATION/final_go_sheet/ALL_SOURCES_MERGED_GO_comined.csv", stringsAsFactors = FALSE)

# Read individual source files
eggnog <- read.csv("eggnog_go_labeled.csv", stringsAsFactors = FALSE)
interpro <- read.csv("InterPro_GO_labeled.csv", stringsAsFactors = FALSE)
deepfri <- read.csv("/home/naren/Desktop/ANAND/PAN_GENOME/ANNOTATION/GO_ENRICHMENT_FILE_PREPARATION/final_go_sheet/DeepFRI_GO_0.5_filtered.csv", stringsAsFactors = FALSE)
swissprot <- read.csv("/home/naren/Desktop/ANAND/PAN_GENOME/ANNOTATION/GO_ENRICHMENT_FILE_PREPARATION/final_go_sheet/SWISSPROT_GO_labeled.csv", stringsAsFactors = FALSE)

# Standardize column names
colnames(eggnog)[1:2] <- c("GeneID", "GO_terms")
colnames(interpro)[1:2] <- c("GeneID", "GO_terms")
colnames(deepfri)[1:2] <- c("GeneID", "GO_terms")
colnames(swissprot)[1:2] <- c("GeneID", "GO_terms")
colnames(merged_go)[1:2] <- c("GeneID", "GO_terms")

# Function to process GO data
process_go_data <- function(df, source_name) {
  df %>%
    mutate(
      GO_terms = ifelse(is.na(GO_terms) | GO_terms == "", NA, GO_terms),
      GO_count = ifelse(is.na(GO_terms), 0,
                        sapply(strsplit(as.character(GO_terms), "[,;]"), length)),
      Source = source_name,
      Has_Annotation = !is.na(GO_terms) & GO_terms != ""
    )
}

# Process all datasets
eggnog_proc <- process_go_data(eggnog, "EggNOG")
interpro_proc <- process_go_data(interpro, "InterProScan")
deepfri_proc <- process_go_data(deepfri, "DeepFRI")
swissprot_proc <- process_go_data(swissprot, "SwissProt")
merged_proc <- process_go_data(merged_go, "Merged")

# ============================================================================
# 3. CREATE COMPREHENSIVE DATASET
# ============================================================================

# Get all unique genes
all_genes <- unique(merged_go$GeneID)

# Create comprehensive dataset
comprehensive_data <- data.frame(GeneID = all_genes) %>%
  left_join(merged_go %>% select(GeneID, Gene_Class, GO_terms), by = "GeneID") %>%
  mutate(
    # Track which sources have annotations for each gene
    EggNOG = GeneID %in% eggnog_proc$GeneID[eggnog_proc$Has_Annotation],
    InterProScan = GeneID %in% interpro_proc$GeneID[interpro_proc$Has_Annotation],
    DeepFRI = GeneID %in% deepfri_proc$GeneID[deepfri_proc$Has_Annotation],
    SwissProt = GeneID %in% swissprot_proc$GeneID[swissprot_proc$Has_Annotation],

    # Count number of sources
    N_Sources = EggNOG + InterProScan + DeepFRI + SwissProt,

    # Total unique GO terms (from merged file)
    Total_GO_Terms = ifelse(is.na(GO_terms) | GO_terms == "", 0,
                            sapply(strsplit(as.character(GO_terms), "[,;]"), length))
  )

# Add individual source GO counts
add_source_go_count <- function(gene_id, source_df) {
  idx <- which(source_df$GeneID == gene_id)
  if(length(idx) == 0) return(0)
  return(source_df$GO_count[idx[1]])
}

comprehensive_data <- comprehensive_data %>%
  mutate(
    GO_EggNOG = sapply(GeneID, function(x) add_source_go_count(x, eggnog_proc)),
    GO_InterProScan = sapply(GeneID, function(x) add_source_go_count(x, interpro_proc)),
    GO_DeepFRI = sapply(GeneID, function(x) add_source_go_count(x, deepfri_proc)),
    GO_SwissProt = sapply(GeneID, function(x) add_source_go_count(x, swissprot_proc))
  )

# ============================================================================
# 4. SUMMARY STATISTICS (Using Correct Total Counts)
# ============================================================================

# Overall statistics by source (percentage based on TOTAL_ALL_GENES)
summary_stats <- data.frame(
  Source = c("EggNOG", "InterProScan", "DeepFRI", "SwissProt", "Merged_Unique"),
  Total_Genes_Annotated = c(
    sum(comprehensive_data$EggNOG),
    sum(comprehensive_data$InterProScan),
    sum(comprehensive_data$DeepFRI),
    sum(comprehensive_data$SwissProt),
    sum(comprehensive_data$Total_GO_Terms > 0)
  ),
  Percentage_Annotated = c(
    sum(comprehensive_data$EggNOG)/TOTAL_ALL_GENES*100,
    sum(comprehensive_data$InterProScan)/TOTAL_ALL_GENES*100,
    sum(comprehensive_data$DeepFRI)/TOTAL_ALL_GENES*100,
    sum(comprehensive_data$SwissProt)/TOTAL_ALL_GENES*100,
    sum(comprehensive_data$Total_GO_Terms > 0)/TOTAL_ALL_GENES*100
  ),
  Mean_GO_per_Gene = c(
    mean(comprehensive_data$GO_EggNOG[comprehensive_data$EggNOG]),
    mean(comprehensive_data$GO_InterProScan[comprehensive_data$InterProScan]),
    mean(comprehensive_data$GO_DeepFRI[comprehensive_data$DeepFRI]),
    mean(comprehensive_data$GO_SwissProt[comprehensive_data$SwissProt]),
    mean(comprehensive_data$Total_GO_Terms[comprehensive_data$Total_GO_Terms > 0])
  ),
  Median_GO_per_Gene = c(
    median(comprehensive_data$GO_EggNOG[comprehensive_data$EggNOG]),
    median(comprehensive_data$GO_InterProScan[comprehensive_data$InterProScan]),
    median(comprehensive_data$GO_DeepFRI[comprehensive_data$DeepFRI]),
    median(comprehensive_data$GO_SwissProt[comprehensive_data$SwissProt]),
    median(comprehensive_data$Total_GO_Terms[comprehensive_data$Total_GO_Terms > 0])
  ),
  Max_GO_per_Gene = c(
    max(comprehensive_data$GO_EggNOG[comprehensive_data$EggNOG]),
    max(comprehensive_data$GO_InterProScan[comprehensive_data$InterProScan]),
    max(comprehensive_data$GO_DeepFRI[comprehensive_data$DeepFRI]),
    max(comprehensive_data$GO_SwissProt[comprehensive_data$SwissProt]),
    max(comprehensive_data$Total_GO_Terms)
  ),
  SD_GO_per_Gene = c(
    sd(comprehensive_data$GO_EggNOG[comprehensive_data$EggNOG]),
    sd(comprehensive_data$GO_InterProScan[comprehensive_data$InterProScan]),
    sd(comprehensive_data$GO_DeepFRI[comprehensive_data$DeepFRI]),
    sd(comprehensive_data$GO_SwissProt[comprehensive_data$SwissProt]),
    sd(comprehensive_data$Total_GO_Terms[comprehensive_data$Total_GO_Terms > 0])
  )
)

# Source overlap statistics (based on genes in CSV only)
overlap_stats <- comprehensive_data %>%
  group_by(N_Sources) %>%
  summarise(Count = n(), Percentage = n()/nrow(comprehensive_data)*100, .groups = "drop")

# Statistics by Gene Class (using correct totals)
class_stats <- comprehensive_data %>%
  filter(!is.na(Gene_Class)) %>%
  group_by(Gene_Class) %>%
  summarise(
    Genes_in_CSV = n(),
    Genes_with_GO = sum(Total_GO_Terms > 0),
    Mean_Sources = mean(N_Sources),
    Mean_GO_Terms = mean(Total_GO_Terms[Total_GO_Terms > 0]),
    Median_GO_Terms = median(Total_GO_Terms[Total_GO_Terms > 0]),
    Max_GO_Terms = max(Total_GO_Terms),
    .groups = "drop"
  ) %>%
  left_join(TOTAL_GENES, by = "Gene_Class") %>%
  mutate(
    Percentage_Annotated = Genes_with_GO/Total_Genes*100
  ) %>%
  select(Gene_Class, Total_Genes, Genes_in_CSV, Genes_with_GO, Percentage_Annotated,
         Mean_Sources, Mean_GO_Terms, Median_GO_Terms, Max_GO_Terms)

# Print and save statistics
cat("=== SUMMARY STATISTICS BY SOURCE ===\n")
print(summary_stats)
write.csv(summary_stats, "GO_annotation_summary_statistics.csv", row.names = FALSE)

cat("\n=== GENES ANNOTATED BY N SOURCES ===\n")
print(overlap_stats)
write.csv(overlap_stats, "GO_source_overlap_statistics.csv", row.names = FALSE)

cat("\n=== STATISTICS BY GENE CLASS ===\n")
print(class_stats)
write.csv(class_stats, "GO_statistics_by_gene_class.csv", row.names = FALSE)

# Save comprehensive dataset
write.csv(comprehensive_data, "Comprehensive_GO_Annotation_Data.csv", row.names = FALSE)

# ============================================================================
# 5. VISUALIZATION SETUP
# ============================================================================

# Create output directories
dir.create("GO_Analysis_Plots", showWarnings = FALSE)
dir.create("GO_Analysis_Plots/with_titles", showWarnings = FALSE)
dir.create("GO_Analysis_Plots/without_titles", showWarnings = FALSE)

# Define color palettes
source_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
names(source_colors) <- c("EggNOG", "InterProScan", "DeepFRI", "SwissProt")

class_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
names(class_colors) <- c("CORE", "SHELL", "SOFTCORE", "CLOUD")

# ============================================================================
# 6. FIGURE 1: ENHANCED UPSET PLOT
# ============================================================================

# Prepare data for UpSet plot
upset_data <- comprehensive_data %>%
  select(EggNOG, InterProScan, DeepFRI, SwissProt) %>%
  mutate(across(everything(), as.numeric))

# Custom UpSet plot with axis label changes
# Order sets: SwissProt (smallest) -> EggNOG -> InterProScan -> DeepFRI (largest)
upset_data_ordered <- upset_data %>%
  select(SwissProt, EggNOG, InterProScan, DeepFRI)

# PDF version with title
pdf("GO_Analysis_Plots/with_titles/01_Enhanced_UpSet_Plot.pdf", width = 14, height = 9)
upset(upset_data_ordered,
      sets = c("SwissProt", "EggNOG", "InterProScan", "DeepFRI"),
      keep.order = TRUE,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#2C3E50",
      sets.bar.color = c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB"),
      matrix.color = "#E74C3C",
      matrix.dot.alpha = 0.6,
      point.size = 4.5,
      line.size = 1.5,
      text.scale = c(2, 1.8, 1.5, 1.5, 2, 1.8),
      number.angles = 0,
      shade.color = "#F8F9FA",
      shade.alpha = 0.3,
      set_size.show = TRUE,
      set_size.scale_max = max(colSums(upset_data_ordered)) * 1.1,
      mb.ratio = c(0.65, 0.35))
grid.text("Source Annotation Overlap - UpSet Plot Analysis",
          x = 0.5, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("No. of Proteins", x = 0.15, y = 0.75, rot = 90,
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# PNG version with title
png("GO_Analysis_Plots/with_titles/01_Enhanced_UpSet_Plot.png",
    width = 14, height = 9, units = "in", res = 400)
upset(upset_data_ordered,
      sets = c("SwissProt", "EggNOG", "InterProScan", "DeepFRI"),
      keep.order = TRUE,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#2C3E50",
      sets.bar.color = c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB"),
      matrix.color = "#E74C3C",
      matrix.dot.alpha = 0.6,
      point.size = 4.5,
      line.size = 1.5,
      text.scale = c(2, 1.8, 1.5, 1.5, 2, 1.8),
      number.angles = 0,
      shade.color = "#F8F9FA",
      shade.alpha = 0.3,
      set_size.show = TRUE,
      set_size.scale_max = max(colSums(upset_data_ordered)) * 1.1,
      mb.ratio = c(0.65, 0.35))
grid.text("Source Annotation Overlap - UpSet Plot Analysis",
          x = 0.5, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("No. of Proteins", x = 0.15, y = 0.75, rot = 90,
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# PDF version without title
pdf("GO_Analysis_Plots/without_titles/01_Enhanced_UpSet_Plot_no_title.pdf", width = 14, height = 9)
upset(upset_data_ordered,
      sets = c("SwissProt", "EggNOG", "InterProScan", "DeepFRI"),
      keep.order = TRUE,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#2C3E50",
      sets.bar.color = c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB"),
      matrix.color = "#E74C3C",
      matrix.dot.alpha = 0.6,
      point.size = 4.5,
      line.size = 1.5,
      text.scale = c(2, 1.8, 1.5, 1.5, 2, 1.8),
      number.angles = 0,
      shade.color = "#F8F9FA",
      shade.alpha = 0.3,
      set_size.show = TRUE,
      set_size.scale_max = max(colSums(upset_data_ordered)) * 1.1,
      mb.ratio = c(0.65, 0.35))
#grid.text("No. of Proteins", x = 0.15, y = 0.75, rot = 90,
  #        gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# PNG version without title (high resolution for merging)
png("GO_Analysis_Plots/without_titles/01_Enhanced_UpSet_Plot_no_title.png",
    width = 14, height = 9, units = "in", res = 1200)
upset(upset_data_ordered,
      sets = c("SwissProt", "EggNOG", "InterProScan", "DeepFRI"),
      keep.order = TRUE,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#2C3E50",
      sets.bar.color = c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB"),
      matrix.color = "#E74C3C",
      matrix.dot.alpha = 0.6,
      point.size = 4.5,
      line.size = 1.5,
      text.scale = c(2, 1.8, 1.5, 1.5, 2, 1.8),
      number.angles = 0,
      shade.color = "#F8F9FA",
      shade.alpha = 0.3,
      set_size.show = TRUE,
      set_size.scale_max = max(colSums(upset_data_ordered)) * 1.1,
      mb.ratio = c(0.65, 0.35))
grid.text("No. of Proteins", x = 0.15, y = 0.75, rot = 90,
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# ============================================================================
# 7. FIGURE 2: VENN DIAGRAM
# ============================================================================

genes_list <- list(
  EggNOG = comprehensive_data$GeneID[comprehensive_data$EggNOG],
  InterProScan = comprehensive_data$GeneID[comprehensive_data$InterProScan],
  DeepFRI = comprehensive_data$GeneID[comprehensive_data$DeepFRI],
  SwissProt = comprehensive_data$GeneID[comprehensive_data$SwissProt]
)

# With title
p_venn_genes <- ggVennDiagram(genes_list,
                              label_alpha = 0,
                              label_size = 5,
                              edge_size = 1.5,
                              category.names = c("EggNOG", "InterProScan", "DeepFRI", "SwissProt")) +
  scale_fill_gradient(low = "#E8F5E9", high = "#1B5E20", name = "Gene Count") +
  scale_color_manual(values = rep("black", 4)) +
  labs(title = "Gene Annotation Overlap Across Sources",
       subtitle = "Venn Diagram showing shared and unique gene annotations") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12)
  )

ggsave("GO_Analysis_Plots/with_titles/02_Venn_Diagram_Genes.pdf",
       p_venn_genes, width = 12, height = 9, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/02_Venn_Diagram_Genes.png",
       p_venn_genes, width = 12, height = 9, dpi = 300)

# Without title
p_venn_genes_no_title <- ggVennDiagram(genes_list,
                                       label_alpha = 0,
                                       label_size = 5,
                                       edge_size = 1.5,
                                       category.names = c("EggNOG", "InterProScan", "DeepFRI", "SwissProt")) +
  scale_fill_gradient(low = "#E8F5E9", high = "#1B5E20", name = "Gene Count") +
  scale_color_manual(values = rep("black", 4)) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12)
  )

ggsave("GO_Analysis_Plots/without_titles/02_Venn_Diagram_Genes_no_title.pdf",
       p_venn_genes_no_title, width = 12, height = 9, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/02_Venn_Diagram_Genes_no_title.png",
       p_venn_genes_no_title, width = 12, height = 9, dpi = 300)

# ============================================================================
# 8. FIGURE 3: BAR PLOT - GENES BY N SOURCES (NO GRIDS)
# ============================================================================

# With title
p_n_sources <- ggplot(overlap_stats, aes(x = factor(N_Sources), y = Percentage)) +
  geom_bar(stat = "identity", aes(fill = factor(N_Sources)),
           color = "black", linewidth = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", Percentage,
                                format(Count, big.mark = ","))),
            vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26")) +
  labs(title = "Distribution of Genes by Number of Annotation Sources",
       subtitle = "How many sources contributed annotations to each gene",
       x = "Number of Sources with Annotations",
       y = "Percentage of Genes (%)") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/with_titles/03_Genes_by_N_Sources.pdf",
       p_n_sources, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/03_Genes_by_N_Sources.png",
       p_n_sources, width = 11, height = 8, dpi = 300)

# Without title
p_n_sources_no_title <- ggplot(overlap_stats, aes(x = factor(N_Sources), y = Percentage)) +
  geom_bar(stat = "identity", aes(fill = factor(N_Sources)),
           color = "black", linewidth = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", Percentage,
                                format(Count, big.mark = ","))),
            vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26")) +
  labs(x = "Number of Sources with Annotations",
       y = "Percentage of Genes (%)") +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/without_titles/03_Genes_by_N_Sources_no_title.pdf",
       p_n_sources_no_title, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/03_Genes_by_N_Sources_no_title.png",
       p_n_sources_no_title, width = 11, height = 8, dpi = 300)

# ============================================================================
# 9. FIGURE 4: BOX PLOT - GO BY GENE CLASS (NO GRIDS)
# ============================================================================

go_class_data <- comprehensive_data %>%
  filter(!is.na(Gene_Class) & Total_GO_Terms > 0)

stats_by_class <- go_class_data %>%
  group_by(Gene_Class) %>%
  summarise(
    mean_go = mean(Total_GO_Terms),
    median_go = median(Total_GO_Terms),
    .groups = "drop"
  )

# With title - avoid label overlap by adjusting y position
p_boxplot_class <- ggplot(go_class_data,
                          aes(x = Gene_Class, y = Total_GO_Terms, fill = Gene_Class)) +
  geom_violin(alpha = 0.3, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = 21,
               outlier.size = 2, outlier.alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "#E74C3C", color = "black", stroke = 1.5) +
  geom_text(data = stats_by_class %>% filter(Gene_Class != "CLOUD"),
            aes(x = Gene_Class, y = max(go_class_data$Total_GO_Terms) * 1.05,
                label = sprintf("μ=%.1f\nM=%.0f", mean_go, median_go)),
            vjust = 0, size = 4, fontface = "bold") +
  geom_text(data = stats_by_class %>% filter(Gene_Class == "CLOUD"),
            aes(x = Gene_Class, y = max(go_class_data$Total_GO_Terms) * 0.85,
                label = sprintf("μ=%.1f\nM=%.0f", mean_go, median_go)),
            vjust = 1, size = 4, fontface = "bold") +
  scale_fill_manual(values = class_colors) +
  labs(title = "Distribution of Unique GO Terms per Protein Across Pangenome Classes",
       subtitle = "Red diamond = mean, μ = mean, M = median",
       x = "Pangenome Category",
       y = "Number of Unique GO Terms per Protein") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/with_titles/04_GO_Distribution_by_Gene_Class.pdf",
       p_boxplot_class, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/04_GO_Distribution_by_Gene_Class.png",
       p_boxplot_class, width = 11, height = 8, dpi = 300)

# Without title
p_boxplot_class_no_title <- ggplot(go_class_data,
                                   aes(x = Gene_Class, y = Total_GO_Terms, fill = Gene_Class)) +
  geom_violin(alpha = 0.3, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = 21,
               outlier.size = 2, outlier.alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "#E74C3C", color = "black", stroke = 1.5) +
  scale_fill_manual(values = class_colors) +
  labs(x = "Pangenome Category",
       y = "Number of Unique GO Terms per Protein") +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/without_titles/04_GO_Distribution_by_Gene_Class_no_title.pdf",
       p_boxplot_class_no_title, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/04_GO_Distribution_by_Gene_Class_no_title.png",
       p_boxplot_class_no_title, width = 11, height = 8, dpi = 300)

# ============================================================================
# 10. FIGURE 5: BOX PLOT - GO BY SOURCE (NO GRIDS)
# ============================================================================

go_source_long <- rbind(
  data.frame(Source = "EggNOG",
             GO_Count = comprehensive_data$GO_EggNOG[comprehensive_data$EggNOG]),
  data.frame(Source = "InterProScan",
             GO_Count = comprehensive_data$GO_InterProScan[comprehensive_data$InterProScan]),
  data.frame(Source = "DeepFRI",
             GO_Count = comprehensive_data$GO_DeepFRI[comprehensive_data$DeepFRI]),
  data.frame(Source = "SwissProt",
             GO_Count = comprehensive_data$GO_SwissProt[comprehensive_data$SwissProt])
) %>%
  mutate(Source = factor(Source, levels = c("EggNOG", "InterProScan", "DeepFRI", "SwissProt")))

stats_by_source <- go_source_long %>%
  group_by(Source) %>%
  summarise(mean_go = mean(GO_Count), median_go = median(GO_Count), .groups = "drop")

# With title - avoid label overlap by adjusting positions
p_boxplot_source <- ggplot(go_source_long, aes(x = Source, y = GO_Count, fill = Source)) +
  geom_violin(alpha = 0.3, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = 21,
               outlier.size = 2, outlier.alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "#E74C3C", color = "black", stroke = 1.5) +
  geom_text(data = stats_by_source %>% filter(Source != "DeepFRI"),
            aes(x = Source, y = max(go_source_long$GO_Count) * 1.05,
                label = sprintf("μ=%.1f\nM=%.0f", mean_go, median_go)),
            vjust = 0, size = 4, fontface = "bold") +
  geom_text(data = stats_by_source %>% filter(Source == "DeepFRI"),
            aes(x = Source, y = max(go_source_long$GO_Count) * 0.85,
                label = sprintf("μ=%.1f\nM=%.0f", mean_go, median_go)),
            vjust = 1, size = 4, fontface = "bold") +
  scale_fill_manual(values = source_colors) +
  labs(title = "Distribution of GO Terms per Protein Across Annotation Sources",
       subtitle = "Red diamond = mean, μ = mean, M = median",
       x = "Annotation Source",
       y = "Number of GO Terms per Protein") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/with_titles/05_GO_Distribution_by_Source.pdf",
       p_boxplot_source, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/05_GO_Distribution_by_Source.png",
       p_boxplot_source, width = 11, height = 8, dpi = 300)

# Without title
p_boxplot_source_no_title <- ggplot(go_source_long, aes(x = Source, y = GO_Count, fill = Source)) +
  geom_violin(alpha = 0.3, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = 21,
               outlier.size = 2, outlier.alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "#E74C3C", color = "black", stroke = 1.5) +
  scale_fill_manual(values = source_colors) +
  labs(x = "Annotation Source",
       y = "Number of GO Terms per Protein") +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/without_titles/05_GO_Distribution_by_Source_no_title.pdf",
       p_boxplot_source_no_title, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/05_GO_Distribution_by_Source_no_title.png",
       p_boxplot_source_no_title, width = 11, height = 8, dpi = 300)

# ============================================================================
# 11. FIGURE 6: STACKED BAR - SOURCE COVERAGE BY CLASS (NO GRIDS)
# ============================================================================

source_class_data <- comprehensive_data %>%
  filter(!is.na(Gene_Class)) %>%
  group_by(Gene_Class, N_Sources) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Gene_Class) %>%
  mutate(Percentage = Count/sum(Count)*100,
         Total = sum(Count))

# With title
p_stacked <- ggplot(source_class_data,
                    aes(x = Gene_Class, y = Percentage, fill = factor(N_Sources))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  geom_text(data = source_class_data %>% group_by(Gene_Class) %>% slice(1),
            aes(x = Gene_Class, y = 102, label = sprintf("n=%s", format(Total, big.mark = ","))),
            inherit.aes = FALSE, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26"),
                    name = "Number of\nSources",
                    labels = c("1", "2", "3", "4")) +
  labs(title = "Annotation Source Coverage Across Gene Classes",
       subtitle = "Stacked bar showing proportion of genes annotated by different numbers of sources",
       x = "Gene Class",
       y = "Percentage of Genes (%)") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("GO_Analysis_Plots/with_titles/06_Source_Coverage_by_Class.pdf",
       p_stacked, width = 12, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/06_Source_Coverage_by_Class.png",
       p_stacked, width = 12, height = 8, dpi = 300)

# Without title
p_stacked_no_title <- ggplot(source_class_data,
                             aes(x = Gene_Class, y = Percentage, fill = factor(N_Sources))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26"),
                    name = "Number of\nSources",
                    labels = c("1", "2", "3", "4")) +
  labs(x = "Gene Class",
       y = "Percentage of Genes (%)") +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("GO_Analysis_Plots/without_titles/06_Source_Coverage_by_Class_no_title.pdf",
       p_stacked_no_title, width = 12, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/06_Source_Coverage_by_Class_no_title.png",
       p_stacked_no_title, width = 12, height = 8, dpi = 300)

# ============================================================================
# 12. FIGURE 7: ANNOTATION COVERAGE BAR PLOT (NO GRIDS)
# ============================================================================

coverage_data <- summary_stats %>%
  filter(Source != "Merged_Unique") %>%
  arrange(desc(Percentage_Annotated))

# With title - thinner bars (width = 0.4)
p_coverage <- ggplot(coverage_data,
                     aes(x = reorder(Source, Percentage_Annotated),
                         y = Percentage_Annotated)) +
  geom_bar(stat = "identity", aes(fill = Source), color = "black", linewidth = 0.8, width = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%\n(%s proteins)",
                                Percentage_Annotated,
                                format(Total_Genes_Annotated, big.mark = ","))),
            hjust = -0.1, size = 5, fontface = "bold") +
  scale_fill_manual(values = source_colors) +
  coord_flip() +
  labs(title = "Annotation Coverage by Source",
       subtitle = "Percentage of total proteins annotated by each source",
       x = "Annotation Source",
       y = "Percentage of Proteins Annotated (%)") +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, 100))

ggsave("GO_Analysis_Plots/with_titles/07_Annotation_Coverage.pdf",
       p_coverage, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/07_Annotation_Coverage.png",
       p_coverage, width = 11, height = 8, dpi = 300)

# Without title - thinner bars
p_coverage_no_title <- ggplot(coverage_data,
                              aes(x = reorder(Source, Percentage_Annotated),
                                  y = Percentage_Annotated)) +
  geom_bar(stat = "identity", aes(fill = Source), color = "black", linewidth = 0.8, width = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%\n(%s proteins)",
                                Percentage_Annotated,
                                format(Total_Genes_Annotated, big.mark = ","))),
            hjust = -0.1, size = 5, fontface = "bold") +
  scale_fill_manual(values = source_colors) +
  coord_flip() +
  labs(x = "Annotation Source",
       y = "Percentage of Proteins Annotated (%)") +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, 100))

ggsave("GO_Analysis_Plots/without_titles/07_Annotation_Coverage_no_title.pdf",
       p_coverage_no_title, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/07_Annotation_Coverage_no_title.png",
       p_coverage_no_title, width = 11, height = 8, dpi = 300)

# High resolution version (1200 DPI) for merging
ggsave("GO_Analysis_Plots/without_titles/07_Annotation_Coverage_no_title_1200dpi.png",
       p_coverage_no_title, width = 11, height = 8, dpi = 1200)

# ============================================================================
# 13. FIGURE 8: HEATMAP - SOURCE × CLASS (NO GRIDS)
# ============================================================================

source_class_go <- rbind(
  comprehensive_data %>% filter(EggNOG & !is.na(Gene_Class)) %>%
    select(Gene_Class, GO_EggNOG) %>%
    rename(GO_Count = GO_EggNOG) %>%
    mutate(Source = "EggNOG"),
  comprehensive_data %>% filter(InterProScan & !is.na(Gene_Class)) %>%
    select(Gene_Class, GO_InterProScan) %>%
    rename(GO_Count = GO_InterProScan) %>%
    mutate(Source = "InterProScan"),
  comprehensive_data %>% filter(DeepFRI & !is.na(Gene_Class)) %>%
    select(Gene_Class, GO_DeepFRI) %>%
    rename(GO_Count = GO_DeepFRI) %>%
    mutate(Source = "DeepFRI"),
  comprehensive_data %>% filter(SwissProt & !is.na(Gene_Class)) %>%
    select(Gene_Class, GO_SwissProt) %>%
    rename(GO_Count = GO_SwissProt) %>%
    mutate(Source = "SwissProt")
)

heatmap_data <- source_class_go %>%
  group_by(Source, Gene_Class) %>%
  summarise(Mean_GO = mean(GO_Count, na.rm = TRUE), .groups = "drop")

# With title
p_heatmap <- ggplot(heatmap_data, aes(x = Source, y = Gene_Class, fill = Mean_GO)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1f", Mean_GO)),
            size = 6, fontface = "bold", color = "white") +
  scale_fill_gradientn(
    colors = c("#3288BD", "#66C2A5", "#FFFFBF", "#FC8D59", "#D53E4F"),
    name = "Mean GO\nTerms"
  ) +
  labs(title = "Mean GO Terms per Gene: Source × Gene Class Heatmap",
       subtitle = "Average number of GO terms assigned by each source to different gene classes",
       x = "Annotation Source",
       y = "Gene Class") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/with_titles/08_Heatmap_Source_Class.pdf",
       p_heatmap, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/08_Heatmap_Source_Class.png",
       p_heatmap, width = 11, height = 8, dpi = 300)

# Without title
p_heatmap_no_title <- ggplot(heatmap_data, aes(x = Source, y = Gene_Class, fill = Mean_GO)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1f", Mean_GO)),
            size = 6, fontface = "bold", color = "white") +
  scale_fill_gradientn(
    colors = c("#3288BD", "#66C2A5", "#FFFFBF", "#FC8D59", "#D53E4F"),
    name = "Mean GO\nTerms"
  ) +
  labs(x = "Annotation Source",
       y = "Gene Class") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text = element_text(color = "black", size = 13, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/without_titles/08_Heatmap_Source_Class_no_title.pdf",
       p_heatmap_no_title, width = 11, height = 8, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/08_Heatmap_Source_Class_no_title.png",
       p_heatmap_no_title, width = 11, height = 8, dpi = 300)

# ============================================================================
# 14. FIGURE 9: FACETED BOX PLOT (NO GRIDS)
# ============================================================================

# With title
p_faceted <- ggplot(source_class_go, aes(x = Gene_Class, y = GO_Count, fill = Gene_Class)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, outlier.alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5,
               fill = "red", color = "darkred") +
  facet_wrap(~Source, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = class_colors) +
  labs(title = "GO Term Distribution by Source and Gene Class",
       subtitle = "Faceted view showing how each source annotates different gene classes",
       x = "Gene Class",
       y = "Number of GO Terms per Gene") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#34495E", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 13, color = "white"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/with_titles/09_Faceted_Source_Class.pdf",
       p_faceted, width = 14, height = 11, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/09_Faceted_Source_Class.png",
       p_faceted, width = 14, height = 11, dpi = 300)

# Without title
p_faceted_no_title <- ggplot(source_class_go, aes(x = Gene_Class, y = GO_Count, fill = Gene_Class)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, outlier.alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5,
               fill = "red", color = "darkred") +
  facet_wrap(~Source, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = class_colors) +
  labs(x = "Gene Class",
       y = "Number of GO Terms per Gene") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#34495E", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 13, color = "white"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("GO_Analysis_Plots/without_titles/09_Faceted_Source_Class_no_title.pdf",
       p_faceted_no_title, width = 14, height = 11, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/09_Faceted_Source_Class_no_title.png",
       p_faceted_no_title, width = 14, height = 11, dpi = 300)

# ============================================================================
# 15. FIGURE 10: STATISTICS COMPARISON (NO GRIDS)
# ============================================================================

stats_comparison <- summary_stats %>%
  filter(Source != "Merged_Unique") %>%
  select(Source, Mean_GO_per_Gene, Median_GO_per_Gene, Max_GO_per_Gene) %>%
  pivot_longer(cols = -Source, names_to = "Statistic", values_to = "Value") %>%
  mutate(Statistic = factor(Statistic,
                            levels = c("Mean_GO_per_Gene", "Median_GO_per_Gene", "Max_GO_per_Gene"),
                            labels = c("Mean", "Median", "Maximum")))

# With title
p_stats_comparison <- ggplot(stats_comparison,
                             aes(x = Source, y = Value, fill = Source)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Value)),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  facet_wrap(~Statistic, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = source_colors) +
  labs(title = "GO Terms per Gene: Statistical Comparison Across Sources",
       subtitle = "Mean, Median, and Maximum GO terms assigned per gene by each source",
       x = "Annotation Source",
       y = "Number of GO Terms") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#34495E", color = "black"),
    strip.text = element_text(face = "bold", size = 13, color = "white"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/with_titles/10_Statistics_Comparison.pdf",
       p_stats_comparison, width = 14, height = 7, dpi = 300)
ggsave("GO_Analysis_Plots/with_titles/10_Statistics_Comparison.png",
       p_stats_comparison, width = 14, height = 7, dpi = 300)

# Without title
p_stats_comparison_no_title <- ggplot(stats_comparison,
                                      aes(x = Source, y = Value, fill = Source)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Value)),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  facet_wrap(~Statistic, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = source_colors) +
  labs(x = "Annotation Source",
       y = "Number of GO Terms") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#34495E", color = "black"),
    strip.text = element_text(face = "bold", size = 13, color = "white"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("GO_Analysis_Plots/without_titles/10_Statistics_Comparison_no_title.pdf",
       p_stats_comparison_no_title, width = 14, height = 7, dpi = 300)
ggsave("GO_Analysis_Plots/without_titles/10_Statistics_Comparison_no_title.png",
       p_stats_comparison_no_title, width = 14, height = 7, dpi = 300)

# ============================================================================
# 16. SUMMARY REPORT
# ============================================================================

sink("GO_Annotation_Analysis_Report.txt")
cat(strrep("=", 80), "\n")
cat("GO ANNOTATION STATISTICS REPORT\n")
cat("Bipolaris Pangenome - Four Source Comparison\n")
cat(strrep("=", 80), "\n\n")

cat("ANALYSIS DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. OVERALL DATASET STATISTICS\n")
cat(strrep("-", 80), "\n")
cat("Total genes in pangenome:", format(TOTAL_ALL_GENES, big.mark = ","), "\n")
cat("  - CORE:", format(TOTAL_GENES$Total_Genes[TOTAL_GENES$Gene_Class == "CORE"], big.mark = ","), "\n")
cat("  - SOFTCORE:", format(TOTAL_GENES$Total_Genes[TOTAL_GENES$Gene_Class == "SOFTCORE"], big.mark = ","), "\n")
cat("  - SHELL:", format(TOTAL_GENES$Total_Genes[TOTAL_GENES$Gene_Class == "SHELL"], big.mark = ","), "\n")
cat("  - CLOUD:", format(TOTAL_GENES$Total_Genes[TOTAL_GENES$Gene_Class == "CLOUD"], big.mark = ","), "\n\n")

cat("Genes with GO annotations (in CSV):", format(nrow(comprehensive_data), big.mark = ","), "\n")
cat("Genes with at least 1 source:", format(sum(comprehensive_data$N_Sources > 0), big.mark = ","), "\n")
cat("Overall annotation rate:", sprintf("%.2f%%\n\n",
                                        nrow(comprehensive_data)/TOTAL_ALL_GENES*100))

cat("2. ANNOTATION COVERAGE BY SOURCE\n")
cat(strrep("-", 80), "\n")
print(summary_stats, row.names = FALSE)
cat("\n")

cat("3. SOURCE OVERLAP STATISTICS\n")
cat(strrep("-", 80), "\n")
print(overlap_stats, row.names = FALSE)
cat("\n")

cat("4. STATISTICS BY GENE CLASS\n")
cat(strrep("-", 80), "\n")
print(class_stats, row.names = FALSE)
cat("\n")

cat("5. KEY FINDINGS\n")
cat(strrep("-", 80), "\n")
cat("- Best coverage source:",
    summary_stats$Source[which.max(summary_stats$Percentage_Annotated[1:4])],
    sprintf("(%.1f%%)\n", max(summary_stats$Percentage_Annotated[1:4])))
cat("- Highest mean GO/gene:",
    summary_stats$Source[which.max(summary_stats$Mean_GO_per_Gene[1:4])],
    sprintf("(%.1f terms)\n", max(summary_stats$Mean_GO_per_Gene[1:4])))
cat("- Genes with all 4 sources:",
    format(overlap_stats$Count[overlap_stats$N_Sources == 4], big.mark = ","),
    sprintf("(%.1f%% of genes in CSV)\n", overlap_stats$Percentage[overlap_stats$N_Sources == 4]))

cat("\n6. OUTPUT FILES GENERATED\n")
cat(strrep("-", 80), "\n")
cat("CSV Files:\n")
cat("  - GO_annotation_summary_statistics.csv\n")
cat("  - GO_source_overlap_statistics.csv\n")
cat("  - GO_statistics_by_gene_class.csv\n")
cat("  - Comprehensive_GO_Annotation_Data.csv\n\n")

cat("Plots (PDF + PNG):\n")
cat("  WITH TITLES (GO_Analysis_Plots/with_titles/):\n")
cat("    01. Enhanced UpSet Plot (source combinations)\n")
cat("    02. Venn Diagram (gene overlap)\n")
cat("    03. Bar plot (genes by N sources)\n")
cat("    04. Box plot (GO by gene class)\n")
cat("    05. Box plot (GO by source)\n")
cat("    06. Stacked bar (source coverage by class)\n")
cat("    07. Bar plot (annotation coverage)\n")
cat("    08. Heatmap (source × class)\n")
cat("    09. Faceted box plot (source and class)\n")
cat("    10. Statistics comparison (mean/median/max)\n\n")

cat("  WITHOUT TITLES (GO_Analysis_Plots/without_titles/):\n")
cat("    All 10 figures above also generated without titles\n\n")

cat(strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("All outputs saved in:\n")
cat("  - GO_Analysis_Plots/with_titles/\n")
cat("  - GO_Analysis_Plots/without_titles/\n")
cat(strrep("=", 80), "\n")
sink()

# ============================================================================
# 17. PRINT COMPLETION MESSAGE
# ============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("✓ ANALYSIS COMPLETE!\n")
cat(strrep("=", 80), "\n\n")
cat("Generated Files:\n")
cat("  - 4 CSV summary tables\n")
cat("  - 20 publication-quality plots (10 with titles + 10 without)\n")
cat("  - 1 comprehensive text report\n\n")
cat("Output directories:\n")
cat("  - GO_Analysis_Plots/with_titles/\n")
cat("  - GO_Analysis_Plots/without_titles/\n\n")
cat("Key Statistics:\n")
cat(sprintf("  - Total genes in pangenome: %s\n", format(TOTAL_ALL_GENES, big.mark = ",")))
cat(sprintf("  - Genes with GO annotations: %s (%.1f%%)\n",
            format(nrow(comprehensive_data), big.mark = ","),
            nrow(comprehensive_data)/TOTAL_ALL_GENES*100))
cat(sprintf("  - Best coverage: %s (%.1f%%)\n",
            summary_stats$Source[which.max(summary_stats$Percentage_Annotated[1:4])],
            max(summary_stats$Percentage_Annotated[1:4])))
cat(sprintf("  - Highest mean GO/gene: %s (%.1f terms)\n",
            summary_stats$Source[which.max(summary_stats$Mean_GO_per_Gene[1:4])],
            max(summary_stats$Mean_GO_per_Gene[1:4])))
cat(strrep("=", 80), "\n")


