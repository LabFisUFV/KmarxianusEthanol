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
annotations <- annotations[,c(1,3,4,5)]

################################################################################
## Load Diniz et al. 2017 data

countData <- read.csv("counts_Diniz2017.csv", header = TRUE)

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

################################################################################
## Load Mo et al. 2019 data
countData <- read.delim("counts_Mo2019.csv", header = TRUE, sep = '\t')
countData <- merge(lta, countData, by = "gene_name", all.x = TRUE)
countData <- countData[,-1]
rownames(countData) <- countData$gene_id
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

################################################################################
## Network comparison

num_rows <- nrow(Diniz2017_final_exp)
sample_indices <- sample(1:num_rows, 3826, replace = FALSE)
Diniz2017_final_exp2 <- Diniz2017_final_exp[sample_indices, ]

Kmarx_list <- list(set1 = Diniz2017_final_exp2, set2 = Mo2019_final_exp)

cons_sft <- consensus_SFT_fit(Kmarx_list, 
                              setLabels = c("Diniz2017", "Mo2018"),
                              cor_method = "pearson")

powers <- cons_sft$power

cons_sft$plot

consensus <- consensus_modules(Kmarx_list, 
                               power = powers, 
                               cor_method = "pearson")

consensus_trait <- consensus_trait_cor(consensus, 
                                       cor_method = "pearson")

################################################################################
## Module preservation

kmarx_og <- read.csv("kmarx_og.csv", sep = "\t")

ortho_exp <- exp_genes2orthogroups(Kmarx_list, kmarx_og, summarize = "mean")


