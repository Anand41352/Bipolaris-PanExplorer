# Bipolaris-PanExplorer

## An Interactive Pangenome Dashboard for *Bipolaris sorokiniana*

Bipolaris PanExplorer** is an interactive **R Shiny** application developed to explore the pangenome, functional annotation, virulence-associated genes, evolutionary dynamics, and protein–protein interaction networks of **19 globally distributed *Bipolaris sorokiniana* genomes**.

The dashboard accompanies the manuscript:

> **Integrated Pangenomic and Systems Biology Analyses Reveal the Genomic Basis of Virulence and Adaptation in  Bipolaris sorokiniana***

---

# Authors

**Anand Kumar Shukla¹²**, **Narendra Kadoo¹²***

¹ Biochemical Sciences Division, CSIR–National Chemical Laboratory (CSIR-NCL), Pune 411008, India

² Academy of Scientific and Innovative Research (AcSIR), Ghaziabad 201002, India

**Corresponding author**

Dr. Narendra Kadoo

Email: **ny.kadoo.ncl@csir.res.in**

---

# Live Dashboard

The dashboard is freely available online:

## https://pangenomebipolaris.shinyapps.io/bipolaris_panexplorer/

---

# Dashboard Overview

The dashboard integrates comparative genomics and systems biology into an interactive web interface.

Available modules include:

- **Home**
  - Project overview
  - World map showing strain origins
  - Genome assembly statistics
  - BUSCO completeness assessment

- **Pangenome Summary**
  - Core, soft-core, shell, and cloud genome composition
  - Gene category statistics
  - Orthogroup summary

- **Presence–Absence Matrix**
  - Interactive visualization of orthogroup distribution across all strains

- **Orthogroup Explorer**
  - Orthogroup search
  - Protein sequence retrieval
  - FASTA download

- **Functional Annotation**
  - EggNOG
  - InterProScan
  - Swiss-Prot
  - NCBI NR annotations

- **Virulence Factors**
  - PHI-base annotations
  - CAZyme annotations
  - CAZyme subclass distribution
  - Effector predictions

- **Gene Family Evolution**
  - CAFE expansion and contraction analysis
  - Expanded and contracted orthogroups
  - Heatmap visualization

- **Positive Selection**
  - Core genes under positive selection
  - Functional annotation of positively selected genes

- **Protein–Protein Interaction Network**
  - Core interactome
  - Hub genes
  - Community networks
  - Topological analysis

- **Summary**
  - Downloadable summary tables
  - Publication figures

---

# Repository Structure

```
Bipolaris_PanExplorer/
│
├── app.R
├── pangenome_summary.txt
├── orthogroup_summary.tsv
├── presence_absence_matrix.csv
├── gene_category_summary.csv
├── busco_summary.csv
├── strain_metadata.csv
├── strain_summary.txt
│
├── annotation/
│   ├── annotation_summary.txt
│   ├── eggnog_annotations.csv
│   ├── interproscan_annotations.tsv
│   ├── nr_annotations.tsv
│   └── swissprot_annotations.csv
│
├── cafe/
│   ├── cafe_expanded_contracted_annotations.csv
│   ├── cafe_gene_family_categories.txt
│   ├── cafe_heatmap_data.csv
│   └── cafe_strain_summary.csv
│
├── orthogroups/
│   └── all_orthogroup_sequences.fa
│
├── ppi/
│   ├── Community1_Edges.csv
│   ├── Community1_Nodes.csv
│   ├── Community2_Edges.csv
│   ├── Community2_Nodes.csv
│   ├── Community3_Edges.csv
│   ├── Community3_Nodes.csv
│   ├── Community4_Edges.csv
│   ├── Community4_Nodes.csv
│   ├── Top10_Edges.csv
│   ├── Top10_Nodes.csv
│   ├── core_ppi_network.svg
│   └── core_interactome/
│
├── selection_pressure/
│   ├── core_positive_selection.csv
│   ├── positively_selected_gene_annotations.csv
│   └── selection_method_comparison.png
│
├── summary/
│
└── virulence_category/
    ├── cazyme_annotations.csv
    ├── cazyme_subclass_distribution.csv
    ├── cazyme_orthogroups.txt
    ├── effector_annotations.csv
    ├── effector_summary.txt
    ├── phi_annotations.tsv
    └── phi_phenotype_distribution.csv
```

---

# Installation

Clone the repository

```bash
git clone https://github.com/Anand41352/Bipolaris-PanExplorer.git

cd Bipolaris-PanExplorer
```

---

# Requirements

- R (≥ 4.3)
- RStudio (recommended)

Required R packages include

```r
shiny
shinydashboard
DT
dplyr
data.table
ggplot2
plotly
leaflet
igraph
visNetwork
tidyr
stringr
RColorBrewer
shinyWidgets
shinycssloaders
```

Install missing packages

```r
packages <- c(
  "shiny",
  "shinydashboard",
  "DT",
  "dplyr",
  "data.table",
  "ggplot2",
  "plotly",
  "leaflet",
  "igraph",
  "visNetwork",
  "tidyr",
  "stringr",
  "RColorBrewer",
  "shinyWidgets",
  "shinycssloaders"
)

install.packages(
  setdiff(packages, rownames(installed.packages()))
)
```

---

# Running the Dashboard

Open R or RStudio.

Set the project directory

```r
setwd("Bipolaris-PanExplorer")
```

Launch the application

```r
shiny::runApp()
```

The application will automatically load all required datasets from the repository.

---

# Updating the Online Dashboard

The dashboard is deployed using **shinyapps.io**.

To deploy a new version:

```r
library(rsconnect)

rsconnect::deployApp(
  "path/to/Bipolaris-PanExplorer"
)
```

---

# Data Included

The repository contains

- Pangenome summary statistics
- Orthogroup summary
- Presence–absence matrix
- Functional annotations
- EggNOG annotations
- InterProScan annotations
- Swiss-Prot annotations
- NCBI NR annotations
- CAZyme annotations
- PHI-base annotations
- Effector annotations
- Gene family evolution (CAFE)
- Positive selection analyses
- Protein–protein interaction network
- Network community datasets
- Orthogroup protein sequences
- Summary tables used in the manuscript

---

# Citation

If you use **Bipolaris PanExplorer** or the associated datasets, please cite:

Shukla AK, Kadoo N.**

*Integrated Pangenomic and Systems Biology Analyses Reveal the Genomic Basis of Virulence and Adaptation in Bipolaris sorokiniana.*

---

# License

This repository is distributed under the **MIT License**.

---

# Contact

# Contact

**Anand Kumar Shukla**  
Biochemical Sciences Division, CSIR–National Chemical Laboratory (CSIR-NCL)  
Pune 411008, Maharashtra, India  
Email: **anand.shukla7066@gmail.com**

**Corresponding Author**

**Dr. Narendra Kadoo**  
Biochemical Sciences Division, CSIR–National Chemical Laboratory (CSIR-NCL)  
Academy of Scientific and Innovative Research (AcSIR)  
Email: **ny.kadoo.ncl@csir.res.in**
