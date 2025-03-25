# ---------------------------------------------------------------
# Synaptic Gene Enrichment Analysis using Enrichr
# Author: Anny Silva Adri
# Date: 03/25/2025
# ---------------------------------------------------------------

# ==== 1. Install (if necessary) and load required packages ====
if (!require("enrichR")) install.packages("enrichR", dependencies = TRUE)

library(enrichR)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(readxl)
library(ggplot2)

# ==== 2. Set working directory ====
# (Change this path to your local environment if needed)
setwd("/Users/annysilvaadri/Documents/Transcriptoma certo/Meus peipers/Paper 1/Novas análises Synapse/Metanalise PBMC/teste synapse com reactome e bp")

# ==== 3. Load input gene list ====
input_data <- read_excel("~/Documents/Transcriptoma certo/Metanalise/Sangue bulk novo/sig_notconverted.xlsx")

# ==== 4. Define enrichment function ====
search_synaptic_genes_enrichr <- function(genes) {
  
  # Relevant databases for synaptic biology
  relevant_dbs <- c(
    "GO_Biological_Process_2025",
    "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025",
    "SynGO_2024",
    "KEGG_2021_Human",
    "WikiPathways_2024_Human",
    "Reactome_Pathways_2024",
    "MSigDB_Hallmark_2024",
    "Human_Phenotype_Ontology",
    "Jensen_COMPARTMENTS"
  )
  
  # Check which of them are available in Enrichr
  available_dbs <- enrichR::listEnrichrDbs()$libraryName
  valid_dbs <- relevant_dbs[relevant_dbs %in% available_dbs]
  
  if (length(valid_dbs) == 0) {
    stop("No valid databases found in Enrichr.")
  }
  
  # Keyword to identify synapse-related terms
  keyword <- "synap"
  
  # Initialize result containers
  gene_counts <- data.frame(Gene = character(), Database = character(), stringsAsFactors = FALSE)
  enrichment_tables <- list()
  
  # Loop through each valid database
  for (db in valid_dbs) {
    cat("Querying database:", db, "\n")
    
    enrichment <- tryCatch({
      enrichR::enrichr(genes, database = db)
    }, error = function(e) NULL)
    
    if (!is.null(enrichment) && db %in% names(enrichment)) {
      results <- enrichment[[db]]
      
      filtered <- results %>%
        filter(str_detect(Term, regex(keyword, ignore_case = TRUE))) %>%
        dplyr::select(Term, Genes, Adjusted.P.value) %>%
        separate_rows(Genes, sep = ";") %>%
        mutate(Genes = trimws(Genes), Database = db) %>%
        distinct()
      
      gene_counts <- bind_rows(gene_counts, filtered %>% select(Gene = Genes, Database))
      
      enrichment_tables[[db]] <- filtered %>%
        select(Database, Term = Term, Genes, Adj_P = Adjusted.P.value)
    }
  }
  
  # Count in how many databases each gene appeared (only if in 2+)
  gene_summary <- gene_counts %>%
    distinct(Gene, Database) %>%
    group_by(Gene) %>%
    summarise(Num_Databases = n(), .groups = "drop") %>%
    filter(Num_Databases >= 2)
  
  enrichment_tables[["Synaptic_Genes"]] <- gene_summary
  
  return(enrichment_tables)
}

# ==== 5. Run enrichment analysis ====
input_genes <- input_data$Genes
enrichment_result <- search_synaptic_genes_enrichr(input_genes)

# ==== 6. Save output as Excel (multiple sheets) ====
write.xlsx(enrichment_result, "synapse_enrichment_results.xlsx", rowNames = FALSE)

# ==== 7. Filter original data to retain only synaptic DEGs ====
synaptic_genes <- enrichment_result[["Synaptic_Genes"]]$Gene

filtered_sig <- input_data %>%
  filter(Genes %in% synaptic_genes) %>%
  select(Genes, metap, metafc)

# Save filtered DEGs
write.xlsx(filtered_sig, "synapse_genes_filtered_min2dbs.xlsx", rowNames = FALSE)

