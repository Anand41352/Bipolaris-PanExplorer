```r
#!/usr/bin/env Rscript

# ============================================================================
# GENOME COMPARTMENTALIZATION ANALYSIS
# Bipolaris sorokiniana – LK93 strain
#
# Statistical analysis:
#   1. Separate 5' and 3' FIR comparisons
#   2. Two-sided Wilcoxon rank-sum tests
#   3. Benjamini-Hochberg correction across all 10 tests
#      (5 functional categories x 2 FIR directions)
#   4. Category-specific background excluding focal category genes
#
# Figure:
#   - Original colours retained
#   - Original panel arrangement retained
#   - Original legend styling retained
#   - Original output dimensions retained (24 x 19.2 inches)
#   - Resolution: 600 dpi
# ============================================================================


# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================

library(ggplot2)
library(hexbin)
library(data.table)
library(dplyr)
library(cowplot)
library(stringr)
library(grid)


# ============================================================================
# SET WORKING DIRECTORY
# ============================================================================

setwd("~/Desktop/ANAND/PAN_GENOME/COMPARTMENTALIZATION/3_RUN_FOR_STRAIN/2_LK93/RUN")

cat("\n")
cat("============================================================\n")
cat("GENOME COMPARTMENTALIZATION ANALYSIS – LK93\n")
cat("============================================================\n\n")


# ============================================================================
# STEP 1: PARSE GFF3 FILE
# ============================================================================

cat("STEP 1: Parsing GFF3 file...\n")

gff <- fread(
  "LK93_clean.gff3",
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  comment.char = "#"
)

# Retain valid GFF3 records
gff <- gff[!is.na(V3)]

# Assign GFF3 column names
setnames(
  gff,
  c(
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes"
  )
)

# Extract gene features only
gff_genes <- gff[feature == "gene"]

# Extract gene IDs from the attributes column
gff_genes[, gene_id := str_extract(attributes, "ID=[^;]+")]
gff_genes[, gene_id := sub("^ID=", "", gene_id)]

# Remove records without gene IDs
gff_genes <- gff_genes[!is.na(gene_id) & gene_id != ""]

# Retain required columns and convert coordinates to integers
gff_genes <- gff_genes[, .(
  seqname,
  gene_id,
  start = as.integer(start),
  end = as.integer(end)
)]

# Remove invalid coordinates
gff_genes <- gff_genes[
  !is.na(start) &
  !is.na(end) &
  start > 0 &
  end >= start
]

# Check for duplicated gene IDs
if (anyDuplicated(gff_genes$gene_id)) {
  stop(
    "Duplicated gene IDs detected in the GFF3 gene records. ",
    "Please resolve these before calculating FIRs."
  )
}

# Order genes by chromosome and genomic position
setorder(gff_genes, seqname, start)

cat(sprintf("  Total genes loaded: %d\n", nrow(gff_genes)))
cat(sprintf(
  "  Number of chromosomes: %d\n",
  uniqueN(gff_genes$seqname)
))

cat("\nFirst three gene records:\n")
print(head(gff_genes, 3))

if (nrow(gff_genes) == 0) {
  stop("No valid gene features found in the GFF3 file.")
}

cat("\n")


# ============================================================================
# STEP 2: CALCULATE FLANKING INTERGENIC REGIONS (FIRs)
# ============================================================================

cat("STEP 2: Calculating FIRs...\n")

fir_list <- list()

for (chr in unique(gff_genes$seqname)) {

  chr_genes <- copy(gff_genes[seqname == chr])

  setorder(chr_genes, start)

  # At least three genes are required to retain genes with
  # both upstream and downstream neighbouring genes.
  if (nrow(chr_genes) < 3) {
    next
  }

  # Distance to the preceding gene:
  # Genomic 5' flanking intergenic region
  chr_genes[, fir_5prime :=
              start - data.table::shift(end, n = 1)]

  # Distance to the following gene:
  # Genomic 3' flanking intergenic region
  chr_genes[, fir_3prime :=
              data.table::shift(start, n = 1, type = "lead") - end]

  chr_genes[, chromosome := chr]

  fir_list[[chr]] <- chr_genes
}

# Combine chromosome-specific results
fir_data <- rbindlist(
  fir_list,
  use.names = TRUE,
  fill = TRUE
)

# Exclude chromosome-terminal genes and loci without valid
# upstream or downstream intergenic regions.
fir_data <- fir_data[
  !is.na(fir_5prime) &
  !is.na(fir_3prime)
]

# Exclude overlapping or adjacent gene configurations
# with zero or negative intergenic distances.
fir_data <- fir_data[
  fir_5prime > 0 &
  fir_3prime > 0
]

# Ensure one record per gene
if (anyDuplicated(fir_data$gene_id)) {
  stop("Duplicated gene IDs detected in the FIR dataset.")
}

if (nrow(fir_data) == 0) {
  stop("No genes with valid 5' and 3' FIR values were retained.")
}

cat(sprintf(
  "  Genes with valid 5' and 3' FIRs: %d\n",
  nrow(fir_data)
))

cat(sprintf(
  "  Genes excluded from FIR analysis: %d\n",
  nrow(gff_genes) - nrow(fir_data)
))

# Genome-wide median FIR values on log10 scale
GENOME_MEDIAN_X <- median(log10(fir_data$fir_5prime))
GENOME_MEDIAN_Y <- median(log10(fir_data$fir_3prime))

cat(sprintf(
  "  Genome median 5' FIR (log10): %.4f\n",
  GENOME_MEDIAN_X
))

cat(sprintf(
  "  Genome median 3' FIR (log10): %.4f\n",
  GENOME_MEDIAN_Y
))

cat("\n")


# ============================================================================
# STEP 3: LOAD FUNCTIONAL GENE CATEGORY FILES
# ============================================================================

cat("STEP 3: Loading gene category files...\n")

load_category <- function(filename, category_name) {

  if (!file.exists(filename)) {
    stop(
      paste(
        "Required category file not found:",
        filename
      )
    )
  }

  genes <- scan(
    filename,
    what = character(),
    quiet = TRUE
  )

  # Remove empty entries and duplicate gene IDs
  genes <- unique(genes[genes != ""])

  cat(sprintf(
    "  %s: %d unique gene IDs loaded\n",
    category_name,
    length(genes)
  ))

  return(genes)
}


category_files <- list(
  "Effectors" = "EFFECTOR_LK93.txt",
  "CAZymes" = "CAZYMES_LK93.txt",
  "PHI" = "PHI_LK93.txt",
  "Positive_Selection" = "POSITIVE_SELECTION_LK93.txt",
  "CAFE" = "CAFE_LK93.txt"
)


gene_categories <- list()

for (category_name in names(category_files)) {

  gene_categories[[category_name]] <- load_category(
    category_files[[category_name]],
    category_name
  )
}

cat("\n")


# ============================================================================
# STEP 4: SEPARATE STATISTICAL TESTING OF 5' AND 3' FIR DISTRIBUTIONS
# ============================================================================

cat("STEP 4: Performing separate 5' and 3' FIR statistical tests...\n")

stats_list <- list()

for (category_name in names(gene_categories)) {

  cat("\n")
  cat("  Analysing:", category_name, "\n")

  # Obtain unique genes belonging to the focal category
  cat_genes <- unique(gene_categories[[category_name]])

  # Retain only genes present in the valid FIR dataset
  cat_genes_in_data <- intersect(
    cat_genes,
    fir_data$gene_id
  )

  # Category-specific background:
  # Exclude genes belonging to the focal category.
  #
  # Genes belonging to other functional categories are retained
  # unless they also belong to the focal category.
  other_genes <- setdiff(
    fir_data$gene_id,
    cat_genes_in_data
  )

  # Check the number of genes available for comparison
  if (length(cat_genes_in_data) < 5) {

    warning(
      paste(
        "Skipping", category_name,
        "- fewer than five category genes have valid FIR values."
      )
    )

    next
  }

  if (length(other_genes) < 5) {

    warning(
      paste(
        "Skipping", category_name,
        "- fewer than five background genes are available."
      )
    )

    next
  }

  # --------------------------------------------------------------------------
  # Extract 5' FIR values
  # --------------------------------------------------------------------------

  cat_fir_5 <- fir_data[
    gene_id %in% cat_genes_in_data,
    fir_5prime
  ]

  other_fir_5 <- fir_data[
    gene_id %in% other_genes,
    fir_5prime
  ]

  # --------------------------------------------------------------------------
  # Extract 3' FIR values
  # --------------------------------------------------------------------------

  cat_fir_3 <- fir_data[
    gene_id %in% cat_genes_in_data,
    fir_3prime
  ]

  other_fir_3 <- fir_data[
    gene_id %in% other_genes,
    fir_3prime
  ]

  # --------------------------------------------------------------------------
  # Separate two-sided Wilcoxon rank-sum tests
  # --------------------------------------------------------------------------

  # 5' FIR comparison
  test_5 <- wilcox.test(
    x = cat_fir_5,
    y = other_fir_5,
    alternative = "two.sided",
    exact = FALSE
  )

  # 3' FIR comparison
  test_3 <- wilcox.test(
    x = cat_fir_3,
    y = other_fir_3,
    alternative = "two.sided",
    exact = FALSE
  )

  # --------------------------------------------------------------------------
  # Calculate medians and fold changes
  # --------------------------------------------------------------------------

  med_cat_5 <- median(cat_fir_5)
  med_other_5 <- median(other_fir_5)

  med_cat_3 <- median(cat_fir_3)
  med_other_3 <- median(other_fir_3)

  fold_5 <- med_cat_5 / med_other_5
  fold_3 <- med_cat_3 / med_other_3

  # --------------------------------------------------------------------------
  # Store category-specific results
  # --------------------------------------------------------------------------

  stats_list[[category_name]] <- data.table(

    Category = category_name,

    N_genes_category = length(cat_genes_in_data),
    N_genes_background = length(other_genes),

    Median_5FIR_category = med_cat_5,
    Median_5FIR_other = med_other_5,
    Fold_5FIR = fold_5,
    P_value_5FIR = test_5$p.value,

    Median_3FIR_category = med_cat_3,
    Median_3FIR_other = med_other_3,
    Fold_3FIR = fold_3,
    P_value_3FIR = test_3$p.value
  )

  cat(sprintf(
    "    Category genes retained: %d\n",
    length(cat_genes_in_data)
  ))

  cat(sprintf(
    "    Background genes: %d\n",
    length(other_genes)
  ))

  cat(sprintf(
    "    5' FIR: median category = %.2f; median background = %.2f; raw P = %.4e\n",
    med_cat_5,
    med_other_5,
    test_5$p.value
  ))

  cat(sprintf(
    "    3' FIR: median category = %.2f; median background = %.2f; raw P = %.4e\n",
    med_cat_3,
    med_other_3,
    test_3$p.value
  ))
}


# Combine all category-specific results
if (length(stats_list) == 0) {
  stop("No functional categories could be analysed.")
}

stats_results <- rbindlist(
  stats_list,
  use.names = TRUE,
  fill = TRUE
)


# ============================================================================
# STEP 4A: BENJAMINI-HOCHBERG CORRECTION ACROSS ALL 10 TESTS
# ============================================================================

cat("\n")
cat("Applying Benjamini-Hochberg correction...\n")

# Collect all raw P-values:
# First five 5' FIR tests, followed by five 3' FIR tests.
all_pvalues <- c(
  stats_results$P_value_5FIR,
  stats_results$P_value_3FIR
)

# Apply BH correction across the complete family of tests
all_bh <- p.adjust(
  all_pvalues,
  method = "BH"
)

n_categories <- nrow(stats_results)

# Assign adjusted P-values back to the corresponding columns
stats_results[
  ,
  P_value_BH_5FIR := all_bh[seq_len(n_categories)]
]

stats_results[
  ,
  P_value_BH_3FIR :=
    all_bh[n_categories + seq_len(n_categories)]
]


# ============================================================================
# STEP 4B: INTERPRET STATISTICAL RESULTS
# ============================================================================

# Interpretation is based on adjusted P-value and median fold change.

stats_results[
  ,
  Interpretation_5FIR :=
    fifelse(
      P_value_BH_5FIR < 0.05 & Fold_5FIR > 1,
      "Significantly larger 5' FIR",
      fifelse(
        P_value_BH_5FIR < 0.05 & Fold_5FIR < 1,
        "Significantly smaller 5' FIR",
        "Not significant"
      )
    )
]

stats_results[
  ,
  Interpretation_3FIR :=
    fifelse(
      P_value_BH_3FIR < 0.05 & Fold_3FIR > 1,
      "Significantly larger 3' FIR",
      fifelse(
        P_value_BH_3FIR < 0.05 & Fold_3FIR < 1,
        "Significantly smaller 3' FIR",
        "Not significant"
      )
    )
]


# Reorder columns for the final output
setcolorder(
  stats_results,
  c(
    "Category",
    "N_genes_category",
    "N_genes_background",
    "Median_5FIR_category",
    "Median_5FIR_other",
    "Fold_5FIR",
    "P_value_5FIR",
    "P_value_BH_5FIR",
    "Interpretation_5FIR",
    "Median_3FIR_category",
    "Median_3FIR_other",
    "Fold_3FIR",
    "P_value_3FIR",
    "P_value_BH_3FIR",
    "Interpretation_3FIR"
  )
)


# Save revised statistical results
fwrite(
  stats_results,
  "Compartmentalization_Statistics_SEPARATE_FIR.csv"
)

cat("\n")
cat("============================================================\n")
cat("SEPARATE FIR STATISTICAL RESULTS\n")
cat("============================================================\n")

print(stats_results)

cat("\nStatistical analysis completed.\n")
cat("Results saved as:\n")
cat("Compartmentalization_Statistics_SEPARATE_FIR.csv\n\n")


# ============================================================================
# STEP 5: DEFINE MANUSCRIPT FIGURE THEME AND PLOT FUNCTIONS
# ============================================================================

cat("STEP 5: Creating manuscript figure templates...\n")


# Original figure theme retained
base_theme <- theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(
      color = "black",
      size = 18
    ),

    axis.title = element_blank(),

    plot.title = element_text(
      face = "bold",
      size = 26,
      hjust = 0.5
    ),

    legend.position = "none"
  )


# --------------------------------------------------------------------------
# STEP 5A: BACKGROUND PLOT – ALL GENES
# --------------------------------------------------------------------------

create_background_plot <- function(fir_data) {

  plot_data <- fir_data[
    ,
    .(
      log5 = log10(fir_5prime),
      log3 = log10(fir_3prime)
    )
  ]

  ggplot(
    plot_data,
    aes(x = log5, y = log3)
  ) +

    geom_hex(
      bins = 30,
      aes(fill = after_stat(count)),
      color = NA
    ) +

    scale_fill_gradient(
      low = "lightgrey",
      high = "black",
      trans = "sqrt"
    ) +

    geom_vline(
      xintercept = GENOME_MEDIAN_X,
      color = "black",
      linetype = "dashed",
      linewidth = 1.0
    ) +

    geom_hline(
      yintercept = GENOME_MEDIAN_Y,
      color = "black",
      linetype = "dashed",
      linewidth = 1.0
    ) +

    labs(
      title = "All Genes"
    ) +

    base_theme
}


# --------------------------------------------------------------------------
# STEP 5B: FUNCTIONAL CATEGORY PLOT
# --------------------------------------------------------------------------
# Each category plot displays:
#   - Background genes as grey-to-black hexbin density
#   - Focal category genes as blue points
#   - Category density contours in blue
#   - Category median lines in red
#   - Genome-wide median lines as black dashed lines
#   - Separate BH-adjusted P-values for 5' and 3' FIRs
# --------------------------------------------------------------------------

create_compartment_plot <- function(
  fir_data,
  gene_ids,
  category_name,
  p_5,
  p_3
) {

  # Prepare plotting data
  plot_data <- fir_data[
    ,
    .(
      gene_id,
      fir_5prime,
      fir_3prime
    )
  ]

  # Assign genes to focal category or background
  plot_data[
    ,
    group := ifelse(
      gene_id %in% unique(gene_ids),
      "Category",
      "Other"
    )
  ]

  # Log10 transformation
  plot_data[
    ,
    log5 := log10(fir_5prime)
  ]

  plot_data[
    ,
    log3 := log10(fir_3prime)
  ]

  # Separate background and category genes
  other_data <- plot_data[group == "Other"]

  category_data <- plot_data[group == "Category"]

  if (nrow(category_data) == 0) {
    stop(
      paste(
        "No valid FIR records available for plotting:",
        category_name
      )
    )
  }

  # Category-specific median positions
  med_x <- median(
    category_data$log5,
    na.rm = TRUE
  )

  med_y <- median(
    category_data$log3,
    na.rm = TRUE
  )

  # Format separate BH-adjusted P-values
  p_5_text <- formatC(
    p_5,
    format = "e",
    digits = 2
  )

  p_3_text <- formatC(
    p_3,
    format = "e",
    digits = 2
  )

  exact_p_text <- paste0(
    "5' FIR BH-p = ", p_5_text,
    "\n",
    "3' FIR BH-p = ", p_3_text
  )

  # Generate plot
  ggplot(
    other_data,
    aes(x = log5, y = log3)
  ) +

    geom_hex(
      bins = 30,
      aes(fill = after_stat(count)),
      color = NA,
      alpha = 0.8
    ) +

    scale_fill_gradient(
      low = "lightgrey",
      high = "black",
      trans = "sqrt"
    ) +

    # Genome-wide median reference lines
    geom_vline(
      xintercept = GENOME_MEDIAN_X,
      color = "black",
      linetype = "dashed",
      linewidth = 1.0
    ) +

    geom_hline(
      yintercept = GENOME_MEDIAN_Y,
      color = "black",
      linetype = "dashed",
      linewidth = 1.0
    ) +

    # Focal category genes
    geom_point(
      data = category_data,
      aes(x = log5, y = log3),
      colour = "dodgerblue",
      alpha = 0.5,
      size = 2
    ) +

    # Density contours of focal category genes
    stat_density_2d(
      data = category_data,
      aes(x = log5, y = log3),
      colour = "blue",
      linewidth = 0.8,
      bins = 4
    ) +

    # Category-specific median reference lines
    geom_vline(
      xintercept = med_x,
      color = "red",
      linewidth = 1.2
    ) +

    geom_hline(
      yintercept = med_y,
      color = "red",
      linewidth = 1.2
    ) +

    labs(
      title = category_name
    ) +

    base_theme +

    # Separate adjusted P-values for the two FIR directions
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = exact_p_text,
      hjust = -0.05,
      vjust = 1.1,
      fontface = "bold.italic",
      color = "darkred",
      size = 6.0
    )
}


# ============================================================================
# STEP 6: GENERATE INDIVIDUAL PANELS AND ASSEMBLE FIGURE
# ============================================================================

cat("STEP 6: Generating figure panels...\n")

plot_list <- list()

# Panel A: Genome-wide background
plot_list[[1]] <- create_background_plot(fir_data)


# Original panel order retained
categories_order <- c(
  "Effectors",
  "CAZymes",
  "PHI",
  "CAFE",
  "Positive_Selection"
)

display_names <- c(
  "Effectors",
  "CAZymes",
  "PHI",
  "CAFE",
  "Positive Selection"
)


# Generate functional category panels
for (i in seq_along(categories_order)) {

  cat_id <- categories_order[i]

  # Check that statistical results are available
  if (!cat_id %in% stats_results$Category) {
    stop(
      paste(
        "Statistical results are missing for category:",
        cat_id
      )
    )
  }

  # Retrieve separate BH-adjusted P-values
  p_5 <- stats_results[
    Category == cat_id,
    P_value_BH_5FIR
  ]

  p_3 <- stats_results[
    Category == cat_id,
    P_value_BH_3FIR
  ]

  # Generate category plot
  plot_list[[i + 1]] <- create_compartment_plot(

    fir_data = fir_data,

    gene_ids = gene_categories[[cat_id]],

    category_name = display_names[i],

    p_5 = p_5,

    p_3 = p_3
  )

  cat(sprintf(
    "  Generated panel: %s\n",
    display_names[i]
  ))
}


# Assemble six panels in the original 3-column x 2-row layout
main_grid <- cowplot::plot_grid(

  plotlist = plot_list,

  labels = LETTERS[1:6],

  label_size = 28,

  ncol = 3,

  nrow = 2,

  align = "hv"
)


# --------------------------------------------------------------------------
# ADD COMMON X- AND Y-AXIS LABELS
# --------------------------------------------------------------------------

x_label_plot <- ggplot() +

  labs(
    x = expression(
      bold(
        log[10] * "(5' Flanking Intergenic Region, bp)"
      )
    )
  ) +

  theme_void() +

  theme(
    axis.title.x = element_text(
      face = "bold",
      size = 22,
      margin = margin(t = 10, b = 5)
    )
  )


y_label_plot <- ggplot() +

  labs(
    y = expression(
      bold(
        log[10] * "(3' Flanking Intergenic Region, bp)"
      )
    )
  ) +

  theme_void() +

  theme(
    axis.title.y = element_text(
      face = "bold",
      size = 22,
      angle = 90,
      margin = margin(r = 10)
    )
  )


# Combine y-axis label with main grid
grid_with_y <- cowplot::plot_grid(

  y_label_plot,

  main_grid,

  rel_widths = c(0.03, 0.97),

  align = "v"
)


# Combine x-axis label with main grid
complete_grid <- cowplot::plot_grid(

  grid_with_y,

  x_label_plot,

  rel_heights = c(0.95, 0.05),

  ncol = 1
)


# ============================================================================
# STEP 7: GENERATE HORIZONTAL UNIFIED LEGEND
# ============================================================================
# Original legend colours, line types, key dimensions and positioning retained.
# ============================================================================

cat("STEP 7: Creating unified horizontal legend...\n")

legend_dummy_data <- data.table(
  x = c(1, 2),
  y = c(1, 2)
)


legend_plot <- ggplot(legend_dummy_data) +

  # Hexbin gene-density legend
  geom_hex(
    data = fir_data,

    aes(
      x = log10(fir_5prime),
      y = log10(fir_3prime),
      fill = after_stat(count)
    ),

    color = NA
  ) +

  scale_fill_gradient(

    low = "lightgrey",

    high = "black",

    name = "Gene Count",

    trans = "sqrt",

    guide = guide_colorbar(

      barwidth = unit(10, "cm"),

      barheight = unit(0.7, "cm"),

      title.theme = element_text(
        face = "bold",
        size = 18
      ),

      label.theme = element_text(
        size = 16
      ),

      title.position = "top",

      title.hjust = 0.5,

      order = 2
    )
  ) +

  # Category median legend key
  geom_vline(
    aes(
      xintercept = NA_real_,
      color = "Category Median"
    ),
    linewidth = 1.2
  ) +

  # Genome median legend key
  geom_vline(
    aes(
      xintercept = NA_real_,
      color = "Genome Median"
    ),
    linetype = "dashed",
    linewidth = 1.0
  ) +

  scale_color_manual(

    name = NULL,

    values = c(
      "Category Median" = "red",
      "Genome Median" = "black"
    ),

    guide = guide_legend(

      title = NULL,

      order = 1,

      keywidth = unit(3.0, "cm"),

      override.aes = list(

        linetype = c(
          "solid",
          "dashed"
        ),

        linewidth = c(
          1.2,
          1.2
        )
      )
    )
  ) +

  theme_classic() +

  theme(

    legend.position = "bottom",

    legend.box = "horizontal",

    legend.text = element_text(
      size = 16,
      face = "bold"
    ),

    legend.spacing.x = unit(1.0, "cm"),

    legend.background = element_blank(),

    legend.box.background = element_blank()
  )


# Extract combined legends
combined_legends <- cowplot::get_legend(legend_plot)


# ============================================================================
# STEP 8: ASSEMBLE AND SAVE FINAL MANUSCRIPT FIGURE
# ============================================================================

cat("STEP 8: Assembling final manuscript figure...\n")

final_figure <- cowplot::plot_grid(

  complete_grid,

  combined_legends,

  ncol = 1,

  rel_heights = c(1, 0.08)
)


# Save final figure
# Original dimensions and resolution retained:
# Width  = 24 inches
# Height = 19.2 inches
# DPI    = 600
# Background = white

ggsave(

  filename = "Figure_Combined_Manuscript_600dpi_NEW.png",

  plot = final_figure,

  width = 24,

  height = 19.2,

  dpi = 600,

  bg = "white"
)


# ============================================================================
# STEP 9: SAVE FINAL ANALYSIS SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nOutput files:\n")

cat("1. Compartmentalization_Statistics_SEPARATE_FIR.csv\n")

cat("2. Figure_Combined_Manuscript_600dpi.png\n")

cat("\nStatistical approach:\n")

cat("- Separate 5' and 3' FIR comparisons\n")

cat("- Two-sided Wilcoxon rank-sum tests\n")

cat("- Category-specific background excluding focal category genes\n")

cat("- BH correction across all 10 comparisons\n")

cat("- Overlapping functional categories analysed independently\n")

cat("\nFigure specifications:\n")

cat("- Original figure colours retained\n")

cat("- Original panel arrangement retained\n")

cat("- Original horizontal legend retained\n")

cat("- Figure dimensions: 24 x 19.2 inches\n")

cat("- Resolution: 600 dpi\n")

cat("- White background\n")

cat("\nAll analyses and figure generation finished.\n")

cat("============================================================\n")
```
