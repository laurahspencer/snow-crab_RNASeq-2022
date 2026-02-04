# 🦀 Snow Crab Ocean Acidification Gene Expression Explorer

An interactive Shiny application for visualizing **gene expression patterns**, **differential expression**, and **PCA structure** across ocean acidification (OA) treatments in juvenile snow crab (*Chionoecetes opilio*).  

Developed for exploring RNA-seq results from multi-treatment OA exposure experiments.

---

## 🔧 Features

- **Interactive filtering** by:
  - Gene ID, gene name, & protein name from [snow crab reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016584305.1/), or UniProt SPID, UniProt ID   
- **Faceted boxplots** showing variance stabilization normalized counts, or those counts further adjusted to control for surrogate variable  
- **DEG highlighting** with comparison labels (e.g., “DEG in ambient vs. severe (88 d)”)
- **Principal Component Analysis (PCA)**:
  - Run using prcomp() on user-supplied gene lists or on genes currently filtered in boxplots 
  - Custom selection of treatments to include  
- Plots built with `ggplot2`  

---

## 📁 Folder Structure

Place all files in a single directory:

snow_crab_expression_explorer/
├── app.R # Main Shiny app 
├── vsd_treatment.rds # DESeq2 variance-stabilized counts
├── Y_adj.rds # SV-regressed expression matrix
├── sample_info.rds # Metadata with sample IDs and treatments
├── snow_annot.rds # Annotation table (gene_id, protein, UniProt, etc.); Uniprot annotations derived from Diamond Blast  
├── deg_info.rds # DEG metadata (gene_id + comparison)
└── README.md # This file

## How to use 

### ▶️ Running locally in R or RStudio:

```r
# Install required packages if not yet installed
install.packages(c("shiny", "tidyverse", "stringr", "DESeq2", "shinyjs"))

# Launch the app - specifying your actual path! 
shiny::runApp("path/to/snow_crab_expression_explorer")

### Access on Shiny Server: 