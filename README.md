# This repo is for the manuscript:

### _Short-term mechanisms, long-term consequences: molecular effects of ocean acidification on juvenile snow crab_
Laura H Spencer1,3, Ingrid Spies3, Jennifer L. Gardner2, Steven Roberts1, W. Christopher Long2

1. University of Washington, School of Aquatic and Fishery Sciences, 1122 NE Boat St, Seattle, WA 98195, USA
2. Kodiak Laboratory, Alaska Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, 301 Research Court, Kodiak, AK 99615, USA
3. Alaska Fisheries Science Center, National Marine Fisheries Service, National Oceanic and Atmospheric Administration, 7600 Sand Point Way NE, Seattle, WA 98115, USA

This repository contains the RNA-seq analysis workflow for a 2022 ocean acidification experiment in juvenile snow crab (_Chionoecetes opilio_). The goal of this project was to characterize short- and long-term transcriptional responses to reduced pH in an OA-tolerant species.

* This is a complementary study to the [Long et al. 2026](https://doi.org/10.1016/j.jembe.2025.152153).
* Raw data is available on NCBI SRA under [BioProject PRJNA1419231](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1419231).
* The manuscript has been submitted to Journal of Experimental Biology and is available as a pre-print on BioRxiv, ID: BIORXIV/2026/703865 (link TBD).
* Expression data is available to explore! Check out this interactive app, [Snow Crab Ocean Acidification Gene Expression Explorer](https://i0o16k-laura0h0spencer.shinyapps.io/snow\_crab\_expression\_explorer/), where you can search for genes by name, protein name, GO term, Uniprot SPID and visualize expression by treatment via boxplots, and PCA.

---

## Project overview

- *Species*: Snow crab (_Chionoecetes opilio_)  
- *Experiment:* Controlled ocean acidification exposure for ~400 days at pH 7.5 (severe OA), pH 7.8 (moderate OA), and ambient/control (pH 8.0)  
- *Tissue:* Homogenized whole=bodies  
- *Approach:* Bulk RNA-seq with differential expression analysis, using the \[snow crab draft genome, GCA\_016584305.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA\_016584305.1/) for STAR alignment and gene-level quantificaton  
- *Primary contrasts:* OA treatments vs. ambient at 8 hours of exposure, and OA treatments at 88 days (no control at 88 days).  

---

## Repository overview
snow-crab_RNASeq-2022/  
├── data/ # Sample metadata, experimental design files  
├── fastqc/ # Raw and trimmed read QC summaries (FastQC, MultiQC)  
├── manuscript/ # Manuscript files (some content excluded from GitHub)  
├── notebooks/ # R analysis notebooks  
├── references/ # Annotation files and reference publications  
├── results/ # Outputs from alignment, differential expression, enrichment  
├── scripts/ # SLURM and utility scripts  
├── snow_crab_expression_explorer/ # Shiny app files  
└── README.md  
