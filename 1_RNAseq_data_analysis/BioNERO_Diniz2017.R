library(dplyr)
library(DESeq2)
library(SummarizedExperiment)
library(BioNERO)

set.seed(137)

################################################################################
## Set up gene annotation data

lta <- read.csv("locus_tag_association.csv", sep = "\t")
Sce_annotations <- read.csv("Sce_annotations.csv", sep = "\t", skipNul = TRUE)
colnames(Sce_annotations) <- c("gene_name", "interpro", "ensembl", "GO")

annotations <- merge(lta, Sce_annotations, by = "gene_name", all.x = TRUE)
annotations <- annotations[,-1]

################################################################################
## Load Diniz et al. 2017 data

countData <- read.csv("counts_Diniz2017.csv", header = TRUE, row.names = 1)

metadata_Diniz2017 <- data.frame(sample_original = colnames(countData), 
                                 condition = c("4h","4h",
                                               "1h","1h",
                                               "0h","0h",
                                               "4h","4h",
                                               "1h","1h",
                                               "0h","0h",
                                               "4h","4h",
                                               "1h","1h",
                                               "0h","0h"))

metadata_Diniz2017$sample <- c("4hRNA9_2","4hRNA9_1",
                               "1hRNA8_2","1hRNA8_1",
                               "0hRNA7_2","0hRNA7_1",
                               "4hRNA6_2","4hRNA6_1",
                               "1hRNA5_2","1hRNA5_1",
                               "0hRNA4_2","0hRNA4_1",
                               "4hRNA3_2","4hRNA3_1",
                               "1hRNA2_2","1hRNA2_1",
                               "0hRNA1_2","0hRNA1_1")


# Run DESeq
col_sel = names(countData)     # Get all but first column name

mdata <- countData %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.", "", name)   # Get the shorter treatment names
  )

meta_df <- data.frame(Sample = names(countData)) %>%
  mutate(
    Type = gsub("-.", "", metadata_Diniz2017$condition)
  ) 

dds <- DESeqDataSetFromMatrix(round(countData),
                              meta_df,
                              design = ~ Type)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

# Construct SummarizedExperiment object
Diniz2017 <- SummarizedExperiment(assays=list(counts=normalized_counts),
                                  colData=metadata_Diniz2017$condition)

################################################################################
## Data processing and EDA

# One-step preprocessing
Diniz2017_final_exp <- exp_preprocess(Diniz2017,
                                      NA_rm = TRUE,               # Replace missing values
                                      remove_nonexpressed = TRUE, # Remove non-expressed genes
                                      Zk_filtering = TRUE,        # Remove outliers
                                      zk = -2,
                                      method = "median",
                                      remove_confounders = TRUE,  # Adjust for confounding artifacts
                                      min_exp = 10,
                                      cor_method = "pearson",
                                      vstransform = FALSE)          # Apply DESeq2's variance stabilizing transformation

# Exploratory data analysis
plot_heatmap(Diniz2017_final_exp,
             type = "samplecor",
             show_rownames = FALSE)

plot_heatmap(Diniz2017_final_exp[1:50, ],
             type = "expr",
             show_rownames = FALSE,
             show_colnames = FALSE)

plot_PCA(Diniz2017_final_exp)

################################################################################
# Infer gene coexpression networks

Diniz2017_sft <- SFT_fit(Diniz2017_final_exp, 
                         net_type = "signed hybrid", 
                         cor_method = "pearson")

Diniz2017_power <- Diniz2017_sft$power

Diniz2017_sft$plot

Diniz2017_net <- exp2gcn(Diniz2017_final_exp, 
                         net_type = "signed hybrid", 
                         SFTpower = Diniz2017_power,
                         cor_method = "pearson")

plot_dendro_and_colors(Diniz2017_net)

plot_eigengene_network(Diniz2017_net)

plot_ngenes_per_module(Diniz2017_net)

module_stability(Diniz2017_final_exp, Diniz2017_net, nRuns = 5)

Diniz2017_MEtrait <- module_trait_cor(exp = Diniz2017_final_exp, MEs = Diniz2017_net$MEs)

plot_module_trait_cor(Diniz2017_MEtrait)

Diniz2017_enrichment <- module_enrichment(net = Diniz2017_net,
                                             background_genes = rownames(Diniz2017_final_exp),
                                             annotation = annotations)

Diniz2017_hubs <- get_hubs_gcn(Diniz2017_final_exp, Diniz2017_net)

## Analyse specific modules
plot_expression_profile(exp = Diniz2017_final_exp,
                        net = Diniz2017_net,
                        plot_module = TRUE,
                        modulename = "darkturquoise"
)

Diniz2017_edges <- get_edge_list(Diniz2017_net, module="darkturquoise")

Diniz2017_edges_filtered <- get_edge_list(Diniz2017_net, 
                                          module = "darkturquoise", 
                                          filter = TRUE)

plot_gcn(edgelist_gcn = Diniz2017_edges_filtered,
         net = Diniz2017_net,
         color_by = "module",
         hubs = Diniz2017_hubs)

################################################################################
GO_terms_to_genes <- read.csv("kmarx_locus_tag_GOs_terms_clean.csv", sep = ",", fileEncoding = "UTF-8")
GO_terms_to_genes <- GO_terms_to_genes[,-3]
GO_terms_to_genes <- GO_terms_to_genes[,c(2,1)]

GO_terms_to_names <- read.csv("kmarx_locus_tag_GOs_terms_clean.csv", sep = ",", fileEncoding = "UTF-8")
GO_terms_to_names <- GO_terms_to_names[,-1]
#GO_terms_to_names <- GO_terms_to_names[,c(2,1)]

# Get the unique colours (WGCNA module IDs) from the second column of your data
colours <- unique(Diniz2017_net$genes_and_modules$Modules)
#colours <- c("coral1","cyan","ivory","steelblue","darkturquoise")

# Define a function to perform GO enrichment analysis for a given WGCNA module
perform_enrichment <- function(module_colour, module_df) {
  # Subset the data to include only the genes belonging to the current module
  module_genes <- module_df$Genes[module_df$Modules == module_colour]
  # Use enricher to perform GO enrichment analysis for the current module
  result <- enricher(gene = module_genes,
                     pvalueCutoff = 0.05,
                     #pAdjustMethod = "bonferroni",
                     TERM2GENE = GO_terms_to_genes,
                     TERM2NAME = GO_terms_to_names)
  # Return the enrichment results
  return(result)
}

# Perform GO enrichment analysis for each module and store the results in a list
enrichment_results <- list()

for (colour in colours) {
  enrichment_results[[colour]] <- perform_enrichment(colour, Diniz2017_net$genes_and_modules)
}

# 0h
dotplot(enrichment_results[["coral1"]], showCategory = 20, title = "Module coral1")
# 1h
dotplot(enrichment_results[["cyan"]], showCategory = 20, title = "Module cyan")
dotplot(enrichment_results[["ivory"]], showCategory = 20, title = "Module ivory")
dotplot(enrichment_results[["steelblue"]], showCategory = 20, title = "Module steelblue")
# 4h
dotplot(enrichment_results[["darkturquoise"]], showCategory = 20, title = "Module darkturquoise")
dotplot(enrichment_results[["darkseagreen4"]], showCategory = 20, title = "Module darkseagreen4")

################################################################################
## Infer gene regulatory networks

# Load TRN
TRN <- read.csv("TRN.csv", header = TRUE, sep = ";")
colnames(TRN) <- c("regulator","gene_name_target","direction","evidence","gene_name")
TRN <- merge(lta, TRN, by = "gene_name", all.x = TRUE)
TRN <- na.omit(TRN)
unique_regulators = unique(TRN$gene_id)

Diniz2017_grn <- exp2grn(exp = Diniz2017_final_exp,
                         regulators = unique_regulators,
                         nTrees = 10)

Diniz2017_grn_hubs <- get_hubs_grn(Diniz2017_grn)

plot_grn(Diniz2017_grn)
