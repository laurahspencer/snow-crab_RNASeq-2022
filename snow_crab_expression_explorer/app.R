#==============================================================
# Snow Crab Gene Expression Explorer (with DEG and GO filters)
#==============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
options(repos = c(
    BiocManager::repositories(version = "3.18"),
    CRAN = "https://cran.rstudio.com"
))

library(shiny)
library(tidyverse)
library(stringr)
library(shinyjs)
library(ggpubr)

#-------------------------------
# Load data
#-------------------------------
vsd_matrix  <- readRDS("vsd_matrix.rds")
Y_adj       <- readRDS("Y_adj.rds")
sample.info <- readRDS("sample_info.rds")
snow.annot  <- readRDS("snow_annot.rds")
deg.info    <- readRDS("deg_info.rds")

#-------------------------------
# Defaults
#-------------------------------
default_genes <- c(
    "GWK47_042301","GWK47_031192","GWK47_039354",
    "GWK47_019251","GWK47_049593","GWK47_003819"
)

treatment_levels <- c(
    "ambient", "moderate_short", "severe_short",
    "moderate_long", "severe_long"
)

treatment_labels <- c(
    "ambient"         = "Control, 0-hr",
    "moderate_short"  = "Moderate OA (pH 7.8), 8-hr exposure",
    "severe_short"    = "Severe OA (pH 7.5), 8-hr exposure",
    "moderate_long"   = "Moderate OA (pH 7.8), 88-d exposure",
    "severe_long"     = "Severe OA (pH 7.5), 88-d exposure"
)

color_scale <- c(
    "ambient"        = "green",
    "moderate_short" = "yellow",
    "moderate_long"  = "#ffc043",
    "severe_short"   = "red",
    "severe_long"    = "darkred"
)

#-------------------------------
# UI
#-------------------------------
ui <- fluidPage(
    useShinyjs(),
    
    tags$head(tags$style(HTML("
    .nav-tabs > li > a {
      font-size: 18px !important;
      font-weight: bold !important;
    }
    .form-group {
      margin-bottom: 6px !important;
    }
  "))),
    
    titlePanel("Snow Crab Ocean Acidification Gene Expression Explorer"),
    
    sidebarLayout(
        sidebarPanel(
            
            checkboxGroupInput(
                "treatments_pca",
                "Select treatments to include:",
                choices  = setNames(treatment_levels, treatment_labels),
                selected = treatment_levels
            ),
            
            checkboxInput("show_deg_only", "Show only DEGs", FALSE),
            
            selectInput(
                "dataset",
                "Select dataset:",
                choices = c(
                    "Normalized counts (VSD)" = "vsd",
                    "SVA-regressed counts (Y_adj)" = "Y_adj"
                ),
                selected = "vsd"
            ),
            
            tags$p(
                tags$strong("Add one or more filters using the below boxes"),
                style = "font-size:16px; margin-bottom:4px;"
            ),
            
            textInput(
                "gene_id",
                "Filter by gene_id(s):",
                value = paste(default_genes, collapse = ", "),
                placeholder = "e.g. ca7"
            ),
            
            tags$p(
                HTML(
                    'Find gene_id&#39;s by differential expression at this link:
           <a href="https://docs.google.com/spreadsheets/d/1wZ32GN-zT5lpuEgUaaFvXcN_qEDqaaRvaAMcDDT74yc/edit?usp=sharing"
              target="_blank">
              <strong>list of DEGs, Top PC Loadings</strong>
           </a>'
                ),
                style = "font-size:15px; margin-top:2px;"
            ),
            br(),
            
            textInput(
                "gene_name",
                "Filter by gene_name (regex):",
                "",
                placeholder = "e.g. ca7"
            ),
            
            textInput(
                "protein",
                "Filter by protein (regex):",
                "",
                placeholder = "e.g. heat shock"
            ),
            
            textInput(
                "spid",
                "Filter by UniProt SPID (regex or list):",
                "",
                placeholder = "e.g. P43166"
            ),
            
            textInput(
                "go_id",
                "Filter by GO Term ID(s):",
                "",
                placeholder = "e.g. GO:0008150 GO:0003674"
            ),
            
            selectInput(
                "go_category",
                "Select GO Category:",
                choices = c(
                    "Any" = "any",
                    "Biological Process (BP)" = "BP",
                    "Molecular Function (MF)" = "MF",
                    "Cellular Component (CC)" = "CC"
                ),
                selected = "any"
            ),
            
            textInput(
                "go_term",
                "Filter by GO term description (regex):",
                "",
                placeholder = "e.g. oxidoreductase, metabolism, membrane"
            ),
            
            checkboxInput("wrap_labels", "Wrap long labels", FALSE),
            numericInput("ncol", "Facet columns:", 3, min = 1, max = 10),
            
            actionButton("go", "Update Boxplots"),
            hr(),
            
            h4("PCA by Gene List"),
            checkboxInput("use_filtered_genes", "Use genes from current boxplot filters", FALSE),
            helpText("Paste a list of gene IDs (one per line or comma-separated):"),
            
            textAreaInput("gene_vector", "Gene IDs for PCA:", "", rows = 6),
            actionButton("run_pca", "Run PCA")
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "Boxplots",
                    h4(textOutput("n_matches"),
                       style = "margin-top:10px;margin-bottom:5px;font-weight:bold;"),
                    plotOutput("expression_plot", height = "900px")
                ),
                
                tabPanel(
                    "PCA",
                    plotOutput("pca_plot", height = "550px"),
                    textOutput("pca_status"),
                    h5("Variance Explained by Principal Components:"),
                    tableOutput("pca_variance")
                )
            )
        )
    )
)

#-------------------------------
# SERVER
#-------------------------------
server <- function(input, output, session) {
    
    observe({ shinyjs::toggleState("gene_vector", !input$use_filtered_genes) })
    
    observeEvent(TRUE, {
        shinyjs::click("go")
    }, once = TRUE)
    
    #-------------------------------
    # Boxplot data
    #-------------------------------
    filtered_data <- eventReactive(input$go, {
        
        expr_mat <- if (input$dataset == "vsd") vsd_matrix else Y_adj
        expr_df  <- expr_mat %>% as.data.frame() %>% rownames_to_column("gene_id")
        
        snow.annot.clean <- snow.annot %>%
            select(
                gene_id, gene_name, protein, spid, gene_uni, protein_names,
                evalue, gene_ontology_i_ds,
                gene_ontology_biological_process,
                gene_ontology_molecular_function,
                gene_ontology_cellular_component,
                gene_ontology_go
            ) %>%
            group_by(gene_id) %>%
            slice_min(evalue, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            distinct(gene_id, .keep_all = TRUE)
        
        ids   <- str_split(input$gene_id, "\\s+|,")[[1]] %>% str_trim() %>% discard(~ .x == "")
        spids <- str_split(input$spid, "\\s+|,")[[1]] %>% str_trim() %>% discard(~ .x == "")
        goids <- str_split(input$go_id, "\\s+|,")[[1]] %>% str_trim() %>% discard(~ .x == "")
        
        filtered_ids <- snow.annot.clean %>%
            filter(
                (length(ids) == 0 | gene_id %in% ids),
                (input$gene_name == "" | grepl(input$gene_name, gene_name, ignore.case = TRUE)),
                (input$protein == "" | grepl(input$protein, protein, ignore.case = TRUE)),
                (length(spids) == 0 | spid %in% spids)
            )
        
        if (length(goids) > 0) {
            filtered_ids <- filtered_ids %>%
                filter(str_detect(gene_ontology_i_ds, paste(goids, collapse = "|")))
        }
        
        if (input$go_term != "") {
            filtered_ids <- filtered_ids %>%
                filter(
                    case_when(
                        input$go_category == "BP" ~ str_detect(gene_ontology_biological_process, regex(input$go_term, TRUE)),
                        input$go_category == "MF" ~ str_detect(gene_ontology_molecular_function, regex(input$go_term, TRUE)),
                        input$go_category == "CC" ~ str_detect(gene_ontology_cellular_component, regex(input$go_term, TRUE)),
                        TRUE ~ str_detect(gene_ontology_go, regex(input$go_term, TRUE))
                    )
                )
        }
        
        filtered_ids <- filtered_ids %>% pull(gene_id) %>% unique()
        
        if (input$show_deg_only) {
            filtered_ids <- intersect(filtered_ids, unique(deg.info$gene_id))
        }
        
        expr_df %>%
            filter(gene_id %in% filtered_ids) %>%
            pivot_longer(starts_with("S")) %>%
            left_join(sample.info[c("sampleID", "treatment")], by = c("name" = "sampleID")) %>%
            filter(treatment %in% input$treatments_pca) %>%
            mutate(treatment = factor(treatment, levels = treatment_levels)) %>%
            left_join(snow.annot.clean, by = "gene_id") %>%
            left_join(deg.info, by = "gene_id") %>%
            mutate(
                DEG_flag = ifelse(!is.na(comparison), "DEG", "non-DEG"),
                label = ifelse(
                    !is.na(comparison),
                    paste0(protein, "\n(", gene_name, " - ", gene_id, ")\nDEG in ", comparison),
                    paste0(protein, "\n(", gene_name, " - ", gene_id, ")")
                ),
                label = if (input$wrap_labels) str_wrap(label, 25) else label
            )
        
    }, ignoreInit = TRUE)
    
    output$n_matches <- renderText({
        fd <- filtered_data()
        if (is.null(fd) || nrow(fd) == 0)
            "No results yet — set filters and click “Update Boxplots”."
        else
            paste("Number of genes matching filters:", n_distinct(fd$gene_id))
    })
    
    output$expression_plot <- renderPlot({
        fd <- filtered_data()
        validate(need(nrow(fd) > 0, "No data to plot."))
        
        ggplot(fd, aes(treatment, value, fill = treatment, alpha = DEG_flag)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, size = 1.2) +
            geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray40") +
            facet_wrap(~label, scales = "free", ncol = input$ncol) +
            scale_fill_manual(values = color_scale, labels = treatment_labels) +
            scale_alpha_manual(values = c("DEG" = 1, "non-DEG" = 0.4), guide = "none") +
            theme_minimal(base_size = 16) +
            theme(
                axis.text.x  = element_blank(),
                axis.title.x = element_blank(),
                legend.position = "top"
            ) +
            ylab("Expression value") +
            ggtitle("Gene Expression by Treatment")
    })
    
    #-------------------------------
    # PCA logic
    #-------------------------------
    pca_data <- eventReactive(input$run_pca, {
        
        expr_mat <- if (input$dataset == "vsd") vsd_matrix else Y_adj
        
        genes <- if (input$use_filtered_genes) {
            req(filtered_data())
            filtered_data() %>% distinct(gene_id) %>% pull(gene_id)
        } else {
            req(input$gene_vector)
            str_split(input$gene_vector, "\\s+|,")[[1]] %>%
                str_trim() %>%
                discard(~ .x == "")
        }
        
        genes <- intersect(rownames(expr_mat), genes)
        validate(need(length(genes) >= 2, "Need at least two valid genes for PCA"))
        
        mat <- expr_mat[genes, , drop = FALSE] %>% t()
        
        samples_to_keep <- sample.info %>%
            filter(treatment %in% input$treatments_pca) %>%
            pull(sampleID)
        
        mat <- mat[rownames(mat) %in% samples_to_keep, , drop = FALSE]
        
        pca <- prcomp(mat, scale. = TRUE)
        
        pc_data <- as.data.frame(pca$x[, 1:2]) %>%
            rownames_to_column("sampleID") %>%
            left_join(sample.info[c("sampleID", "treatment")], by = "sampleID")
        
        var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)
        
        list(
            pca = pca,
            data = pc_data,
            variance = tibble(
                PC = paste0("PC", seq_along(var_exp)),
                Variance_Explained = paste0(var_exp, "%")
            )
        )
    })
    
    output$pca_plot <- renderPlot({
        req(pca_data())
        
        df <- pca_data()$data
        var_exp <- pca_data()$variance$Variance_Explained
        
        ggscatter(
            df,
            x = "PC1",
            y = "PC2",
            col = "treatment",
            size = 3.5,
            alpha = 0.85,
            ellipse = FALSE,
            star.plot = TRUE
        ) +
            theme_minimal(base_size = 16) +
            ggtitle("PCA of Selected Genes") +
            xlab(paste0("PC1 (", var_exp[1], ")")) +
            ylab(paste0("PC2 (", var_exp[2], ")")) +
            scale_color_manual(
                name   = "Treatment",
                values = color_scale,
                labels = treatment_labels
            ) +
            theme(
                legend.position = "right",
                legend.text  = element_text(size = 9),
                legend.title = element_text(size = 10)
            )
    })
    
    output$pca_status <- renderText({
        src <- if (input$use_filtered_genes)
            "current boxplot filters"
        else
            "manual gene list"
        paste("PCA computed using genes from", src)
    })
    
    output$pca_variance <- renderTable({
        req(pca_data())
        head(pca_data()$variance, 5)
    })
}

#-------------------------------
# Run App
#-------------------------------
shinyApp(ui, server)