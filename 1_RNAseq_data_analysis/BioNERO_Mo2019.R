library(dplyr)
library(DESeq2)
library(SummarizedExperiment)
library(BioNERO)
library(clusterProfiler)

set.seed(137)

################################################################################
## Set up gene annotation data

lta <- read.csv("locus_tag_association.csv", sep = "\t")
Sce_annotations <- read.csv("Sce_annotations.csv", sep = "\t", skipNul = TRUE)
colnames(Sce_annotations) <- c("gene_name", "interpro", "ensembl", "GO")

annotations <- merge(lta, Sce_annotations, by = "gene_name", all.x = TRUE)
#annotations <- annotations[,-1]
annotations <- annotations[,c(1,3,4,5)]

################################################################################
## Load Mo et al. 2019 data
countData <- read.delim("counts_Mo2019.csv", header = TRUE, sep = '\t')
#countData <- merge(lta, countData, by = "gene_name", all.x = TRUE)
#countData <- countData[,-1]
#rownames(countData) <- countData$gene_id
rownames(countData) <- countData$gene_name
countData <- countData[,-1]
countData <- na.omit(countData)

metadata_Mo2019 <- data.frame(sample_original = colnames(countData), 
                              condition = c("0%-KM","0%-KM","0%-KM",
                                            "0%-100d","0%-100d","0%-100d",
                                            "4%-KM","4%-KM","4%-KM",
                                            "4%-100d","4%-100d","4%-100d",
                                            "6%-KM","6%-KM","6%-KM",
                                            "6%-100d","6%-100d","6%-100d"))

metadata_Mo2019$sample <- c("0%-KM_1","0%-KM_2","0%-KM_3",
                            "0%-100d_1","0%-100d_2","0%-100d_3",
                            "4%-KM_1","4%-KM_2","4%-KM_3",
                            "4%-100d_1","4%-100d_2","4%-100d_3",
                            "6%-KM_1","6%-KM_2","6%-KM_3",
                            "6%-100d_1","6%-100d_2","6%-100d_3")

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
    Type = gsub("-.", "", metadata_Mo2019$condition)
  ) 

dds <- DESeqDataSetFromMatrix(round(countData),
                              meta_df,
                              design = ~ Type)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)

# Construct SummarizedExperiment object
Mo2019 <- SummarizedExperiment(assays=list(counts=normalized_counts),
                                  colData=metadata_Mo2019$condition)

################################################################################
## Data processing and EDA

# One-step preprocessing
Mo2019_final_exp <- exp_preprocess(Mo2019,
                                   NA_rm = TRUE,               # Replace missing values
                                   remove_nonexpressed = TRUE, # Remove non-expressed genes
                                   Zk_filtering = TRUE,        # Remove outliers
                                   zk = -2,
                                   method = "median",
                                   remove_confounders = TRUE,  # Adjust for confounding artifacts
                                   min_exp = 10,
                                   cor_method = "pearson",
                                   vstransform = FALSE)        # Apply DESeq2's variance stabilizing transformation

# Exploratory data analysis
plot_heatmap(Mo2019_final_exp,
             type = "samplecor",
             show_rownames = FALSE)

plot_heatmap(Mo2019_final_exp[1:50, ],
             type = "expr",
             show_rownames = FALSE,
             show_colnames = FALSE)

plot_PCA(Mo2019_final_exp)

################################################################################
# Infer gene coexpression networks

Mo2019_sft <- SFT_fit(Mo2019_final_exp, 
                      net_type = "signed hybrid", 
                      cor_method = "pearson")

Mo2019_power <- Mo2019_sft$power

Mo2019_sft$plot

Mo2019_net <- exp2gcn(Mo2019_final_exp, 
                      net_type = "signed hybrid", 
                      SFTpower = Mo2019_power,
                      cor_method = "pearson")

plot_dendro_and_colors(Mo2019_net)

plot_eigengene_network(Mo2019_net)

plot_ngenes_per_module(Mo2019_net)

module_stability(Mo2019_final_exp, Mo2019_net, nRuns = 5)

Mo2019_MEtrait <- module_trait_cor(exp = Mo2019_final_exp, MEs = Mo2019_net$MEs)

plot_module_trait_cor(Mo2019_MEtrait)


Mo2019_enrichment <- module_enrichment(net = Mo2019_net,
                                       background_genes = rownames(Mo2019_final_exp),
                                       annotation = annotations)

Mo2019_hubs <- get_hubs_gcn(Mo2019_final_exp, Mo2019_net)

## Analyse specific modules
plot_expression_profile(exp = Mo2019_final_exp,
                        net = Mo2019_net,
                        plot_module = TRUE,
                        modulename = "brown4"
)

Mo2019_edges <- get_edge_list(Mo2019_net, module="brown4")

Mo2019_edges_filtered <- get_edge_list(Mo2019_net, 
                                       module = "brown4", 
                                       filter = TRUE)

plot_gcn(edgelist_gcn = Mo2019_edges_filtered,
         net = Mo2019_net,
         color_by = "module",
         hubs = Mo2019_hubs)

################################################################################
GO_terms_to_genes <- read.csv("kmarx_gene_names_GOs_terms_clean.csv", sep = ",", fileEncoding = "UTF-8")
GO_terms_to_genes <- GO_terms_to_genes[,-3]
GO_terms_to_genes <- GO_terms_to_genes[,c(2,1)]

GO_terms_to_names <- read.csv("kmarx_gene_names_GOs_terms_clean.csv", sep = ",", fileEncoding = "UTF-8")
GO_terms_to_names <- GO_terms_to_names[,-1]
#GO_terms_to_names <- GO_terms_to_names[,c(2,1)]

# Get the unique colours (WGCNA module IDs) from the second column of your data
colours <- unique(Mo2019_net$genes_and_modules$Modules)
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
  enrichment_results[[colour]] <- perform_enrichment(colour, Mo2019_net$genes_and_modules)
}
# 0%-KM
#dotplot(enrichment_results[["salmon4"]], showCategory = 20, title = "Module salmon4")
# 4%-KM
#dotplot(enrichment_results[["floralwhite"]], showCategory = 20, title = "Module floralwhite")
dotplot(enrichment_results[["antiquewhite4"]], showCategory = 20, title = "Module antiquewhite4")
dotplot(enrichment_results[["midnightblue"]], showCategory = 20, title = "Module midnightblue")
# 6%-KM
#dotplot(enrichment_results[["grey"]], showCategory = 20, title = "Module grey")
#dotplot(enrichment_results[["grey60"]], showCategory = 20, title = "Module grey60")
dotplot(enrichment_results[["cyan"]], showCategory = 20, title = "Module cyan")
dotplot(enrichment_results[["bisque4"]], showCategory = 20, title = "Module bisque4")
# 0%-100d
#dotplot(enrichment_results[["coral2"]], showCategory = 20, title = "Module coral2")
dotplot(enrichment_results[["cyan"]], showCategory = 20, title = "Module cyan")
# 4%-100d
dotplot(enrichment_results[["darkred"]], showCategory = 20, title = "Module darkred")
#dotplot(enrichment_results[["steelblue"]], showCategory = 20, title = "Module steelblue")
dotplot(enrichment_results[["darkgrey"]], showCategory = 20, title = "Module darkgrey")
dotplot(enrichment_results[["green"]], showCategory = 20, title = "Module green")
#dotplot(enrichment_results[["lightyellow"]], showCategory = 20, title = "Module lightyellow")
# 6%-100d
#dotplot(enrichment_results[["brown4"]], showCategory = 20, title = "Module brown4")
dotplot(enrichment_results[["honeydew1"]], showCategory = 20, title = "Module honeydew1")
#dotplot(enrichment_results[["darkslateblue"]], showCategory = 20, title = "Module darkslateblue")
dotplot(enrichment_results[["plum1"]], showCategory = 20, title = "Module plum1")

################################################################################
## Infer gene regulatory networks

# Load TRN
TRN <- read.csv("TRN.csv", header = TRUE, sep = ";")
colnames(TRN) <- c("regulator","gene_name_target","direction","evidence","gene_name")
#TRN <- merge(lta, TRN, by = "gene_name", all.x = TRUE)
#TRN <- na.omit(TRN)
unique_regulators = unique(TRN$gene_name)

Mo2019_grn <- exp2grn(exp = Mo2019_final_exp,
                      regulators = unique_regulators,
                      nTrees = 10)

Mo2019_grn_hubs <- get_hubs_grn(Mo2019_grn)

plot_grn(Mo2019_grn)
