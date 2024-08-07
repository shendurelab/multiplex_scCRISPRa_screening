# CLI script to test guide effects, original from Troy & Diego on April 15th, 2022

library(tidyverse); library(repr); library(Matrix)
library(Seurat); library(broom); library(ggridges)
library(ggrepel); library(patchwork); library(data.table)
library(methods)


# Pass to this script (AS CLI)
partition <- commandArgs(trailingOnly=TRUE)[1]
# partition <- '1' # for testing
# k is the parition of guides to test
# Rscript 1mb_DE_Script_snakemake.R k

# Change this to the directory you want it to run in.
setwd("/net/shendure/vol10/projects/troym/ResQTL/nobackup/promoter_pilot/Neuron_1Mb_DE/")

##Be sure to precompute the required datasets using 
#"1mb_DE_DB_Pre_Compute_aggr.R" and save them to this directory 
#to run

# read the previously computed datasets we need
pilot_subset <- readRDS('nobackup/pilot_seurat_object.Rds')
All_1MB_Neighbouring_Genes_2_Test <- readRDS('nobackup/neighboring_genes.Rds')
FBC_DF <- readRDS('nobackup/filtered_FBC_DF.Rds') # using filtered version
FBC_DF_filt <- readRDS('nobackup/filtered_FBC_DF_umi5.Rds')
pilot_subset <- readRDS('nobackup/pilot_seurat_object.Rds')
ngene_lookup <- readRDS('nobackup/ngene_lookup.Rds')

# only keep gRNA tests with the partition the script is using
gRNA_Test_IDs <- (readRDS('nobackup/gRNA_test_ids.Rds') %>% 
	filter(a_partition==partition))$guide_ID
	
##Loop to iterate through the tests

results <- rbindlist(lapply(1:length(gRNA_Test_IDs), function(i) {
	gRNA_Test_ID <- gRNA_Test_IDs[i]
	# print the gRNA being tested
	print(paste0('starting guide ', i, ': ', gRNA_Test_ID))

	# #Isolate gRNA associated with that test
	
	gRNA <- ngene_lookup[[gRNA_Test_ID]]$gRNA_Name

	# #Isolate target gene associated with that test

	target_gene <- ngene_lookup[[gRNA_Test_ID]]$hgnc_symbol

	# Filter for cells with a given guide

	gRNA_FBC_DF <- FBC_DF_filt %>%
		filter(Gene %in% gRNA) 
	gRNA_Cells <- unique(gRNA_FBC_DF$Cell)

	# This makes control all other cells
	# exclude cells in the cells to exclude)
	Control_Cells_DF <- FBC_DF %>% 
		filter(!(Cell %in% gRNA_Cells))
	Control_Cells <- unique(Control_Cells_DF$Cell)

	# print(paste0('starting testing'))
	Result <- FindMarkers(pilot_subset, ident.1 = gRNA_Cells, 
	     ident.2 = Control_Cells, min.pct = 0, min.cells.feature = 0,
		 min.cells.group = 0, features = target_gene, logfc.threshold = 0)

	# print(paste0('finished testing'))
	return(as.data.frame(Result) %>%
	       mutate(target_guide=gRNA, target_gene=target_gene,
	              n_cells=length(gRNA_Cells),
	              n_control_cells=length(Control_Cells)))
}))

# writing results to disk for this partition
results_DF <- as.data.frame(results)
results_DF %>% write_csv(paste0('nobackup/', partition, '_results.txt'))
