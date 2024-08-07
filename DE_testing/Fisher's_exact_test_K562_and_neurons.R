library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(dplyr)

#Want to use the same input dfs as I did for the gRNA analysis, which are the two dfs below;
#note that these do not include NTCs

K562_results = read.table("~/Documents/Shendure_Lab/CRISPRa_QTL/gRNA_analyses/K562_results_all_targets.csv",
                                  sep = ",",
                                  header = TRUE) #397 rows

neuron_results = read.table("~/Documents/Shendure_Lab/CRISPRa_QTL/gRNA_analyses/neuron_results_all_targets.csv",
                                  sep = ",",
                                  header = TRUE) #383 rows

#is the log2FC greater than 0?
K562_results$log2FC_greater_than_0 = NA
K562_results$log2FC_greater_than_0[K562_results$avg_log2FC >0 ] = "yes"
K562_results$log2FC_greater_than_0[K562_results$avg_log2FC <0 ] = "no"

neuron_results$log2FC_greater_than_0 = NA
neuron_results$log2FC_greater_than_0[neuron_results$avg_log2FC >0 ] = "yes"
neuron_results$log2FC_greater_than_0[neuron_results$avg_log2FC <0 ] = "no"

#now want to get the NTCs - this let's me do the Fisher's test
#getting NTCs - this code is from when I did the NTC correlations
K562_target_genes = unique(K562_results$target_gene)
length(K562_target_genes)
#K562_target_genes

neuron_target_genes = unique(neuron_results$target_gene)
length(neuron_target_genes)
#neuron_target_genes

length(neuron_target_genes %in% K562_target_genes)
#all neuron targets are in K562 targets, so use the neuron_target_genes_list; this is 58 genes

#filter aggregated results for NTCs only - K562 dataset
K562_aggr_results = read.table("~/Documents/Shendure_Lab/CRISPRa_QTL/DE_results_primary_targets/K562_aggregated_results.txt",
                               sep = "\t",
                               header = TRUE)

#head(K562_aggr_results)
K562_aggr_results_NTCs_primary_targets = K562_aggr_results %>%
  dplyr::filter(grepl('NTC', target_guide)) %>%
  subset(target_gene %in% neuron_target_genes) 

#K562_aggr_results_NTCs_primary_targets = K562_aggr_results_NTCs_primary_targets[sample(nrow(K562_aggr_results_NTCs_primary_targets), 300), ]

head(K562_aggr_results_NTCs_primary_targets)
dim(K562_aggr_results_NTCs_primary_targets)
length(unique(K562_aggr_results_NTCs_primary_targets$target_gene))

#filter aggregated results for NTCs only - neuron dataset
neuron_aggr_results = read.table("~/Documents/Shendure_Lab/CRISPRa_QTL/DE_results_primary_targets/neuron_aggregated_results.txt",
                                 sep = "\t",
                                 header = TRUE)

#head(neuron_aggr_results)
neuron_aggr_results_NTCs_primary_targets = neuron_aggr_results %>%
  dplyr::filter(grepl('NTC', target_guide)) %>%
  subset(target_gene %in% neuron_target_genes)

#neuron_aggr_results_NTCs_primary_targets = neuron_aggr_results_NTCs_primary_targets[sample(nrow(neuron_aggr_results_NTCs_primary_targets), 300), ]

head(neuron_aggr_results_NTCs_primary_targets)
dim(neuron_aggr_results_NTCs_primary_targets)
length(unique(neuron_aggr_results_NTCs_primary_targets$target_gene))

#Now have the correct input dfs

#fold change less than or greater than 0?
K562_aggr_results_NTCs_primary_targets$log2FC_greater_than_0 = NA
K562_aggr_results_NTCs_primary_targets$log2FC_greater_than_0[K562_aggr_results_NTCs_primary_targets$avg_log2FC >0 ] = "yes"
K562_aggr_results_NTCs_primary_targets$log2FC_greater_than_0[K562_aggr_results_NTCs_primary_targets$avg_log2FC <0 ] = "no"

neuron_aggr_results_NTCs_primary_targets$log2FC_greater_than_0 = NA
neuron_aggr_results_NTCs_primary_targets$log2FC_greater_than_0[neuron_aggr_results_NTCs_primary_targets$avg_log2FC >0 ] = "yes"
neuron_aggr_results_NTCs_primary_targets$log2FC_greater_than_0[neuron_aggr_results_NTCs_primary_targets$avg_log2FC <0 ] = "no"

#now should be able to combine everything into one df
#K562
K562_results$targeting_guide = "targeting gRNA"
K562_aggr_results_NTCs_primary_targets$targeting_guide = "non-targeting gRNA"

K562_primary_results_short = K562_results[,c(12,13)]
K562_aggr_results_NTCs_primary_targets_short = K562_aggr_results_NTCs_primary_targets[,c(10,11)]
K562_targeting_and_NTC = rbind(K562_primary_results_short, K562_aggr_results_NTCs_primary_targets_short)  

#neurons
neuron_results$targeting_guide = "targeting gRNA"
neuron_aggr_results_NTCs_primary_targets$targeting_guide = "non-targeting gRNA"
neuron_primary_results_short = neuron_results[,c(12,13)]
neuron_aggr_results_NTCs_primary_targets_short = neuron_aggr_results_NTCs_primary_targets[,c(10,11)]
neuron_targeting_and_NTC = rbind(neuron_primary_results_short, neuron_aggr_results_NTCs_primary_targets_short)  

#Fischer's test
fisher.test(table(K562_targeting_and_NTC))
fisher.test(table(neuron_targeting_and_NTC))

table(K562_targeting_and_NTC)
table(neuron_targeting_and_NTC)

#leftover old code
