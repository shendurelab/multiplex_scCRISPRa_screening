library(tidyverse); library(repr); library(Matrix)
library(Seurat); library(broom); library(ggridges)
library(ggrepel); library(patchwork); library(data.table)
library(methods)

##change this to the directory you want snakemake to run in.
setwd("/net/shendure/vol10/projects/troym/ResQTL/nobackup/promoter_pilot/K562_1Mb_DE")

#Change these directories to the filtered feature barcode matrix you want to analyze (e.g. aggregated as below of for each lane).
#This will read in the data and save them as .Rds objects to save time. 

data_dir <- paste0('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/')

##Making seurat object from FBC matrix
CRISPRaQTL_Pilot_data <- Read10X(data.dir = data_dir)
pilot <- CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`, min.cells = 4,
	min.features = 200)

##Filtering for high quality cells
pilot[["percent.mt"]] <- PercentageFeatureSet(pilot, pattern = "^MT-")
pilot_subset <- NormalizeData(subset(pilot, subset = percent.mt < 10 &
	nCount_RNA > 4000))

##Reading in FBC matrix
matrix_dir <- paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")


mat <- readMM(file = matrix.path)
feature.names <- read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names <- read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1


##Indexing the matrix for only high quality cells in the seurat
# object (x is the index vector of high quality cell names)
x <- rownames(pilot_subset@meta.data); filtered_fbc_mat <- mat[ ,x]

##Creating a gRNA-specific FBC Matrix
Ft_Ref <- read_csv(paste0("/net/shendure/vol10/projects/ResQTL/",
	"nobackup/10X_single_cell_pilot_expt_Nov_2021/10X_pilot_1_",
	"cellranger_count_SNV_fixed/outs/crispr_analysis/feature_reference.csv"))  %>%
    rename(gRNA = id)
x <- Ft_Ref$gRNA; gRNA_Mat <- filtered_fbc_mat[x, ]



##Calculating the number of cells with each perturbation to only
# input guides with >1 cell
N_Cells_With_Pert_Mat <- rowSums(gRNA_Mat > 5)
N_Cells_With_Pert_DF <- as.data.frame(as.table(N_Cells_With_Pert_Mat))
N_Cells_With_Pert_DF <- N_Cells_With_Pert_DF %>% rename(Gene = Var1) %>%
    rename(N_Cells_With_Perturbation = Freq)
No_Cell_gRNAs <- N_Cells_With_Pert_DF %>%
	filter(N_Cells_With_Perturbation < 2) %>% select(Gene)
No_Cell_gRNAs <- No_Cell_gRNAs$Gene


##Creating a dataframe version of the full filtered FBC
transposed_mat <- t(filtered_fbc_mat)
FBC_DF <- tidy(transposed_mat) %>% rename(Cell = row) %>%
    rename(Gene = column) %>% rename(UMI_Count = value)

##Reading in the list of all gRNA-1mb neighbouring gene pairings for testing
All_1MB_Neighbouring_Genes_2_Test <- read_csv(paste0("/net/shendure/vol10/",
	"projects/troym/ResQTL/nobackup/promoter_pilot/K562_Pilot_gRNA_",
	"Coordinates/CRISPRaQTL_K562_pilot_hg38_gRNA_All_distinct_neighbouring",
	"_genes_1MB_to_test.csv"))
All_1MB_Neighbouring_Genes_2_Test <- All_1MB_Neighbouring_Genes_2_Test %>%
    group_by(gRNA_Name) %>% mutate(Test_ID = row_number()) %>%
    unite("Test_ID",  c("gRNA_Name", "hgnc_symbol", "Test_ID"), remove = FALSE)

##Filter the list to only test genes that are detected
Detected_Genes <- rownames(pilot_subset) 
All_1MB_Neighbouring_Genes_2_Test <- All_1MB_Neighbouring_Genes_2_Test %>%
    filter(hgnc_symbol %in% Detected_Genes)

##Making a unique gRNA-target_gene list to itereate through for DE testing
gRNA_Test_IDs <- All_1MB_Neighbouring_Genes_2_Test %>%
    filter(!gRNA_Name %in% No_Cell_gRNAs)
gRNA_Test_IDs <- gRNA_Test_IDs$Test_ID


##Saving everything to the nobackup folder
saveRDS(pilot_subset, 'nobackup/pilot_seurat_object.Rds')
saveRDS(All_1MB_Neighbouring_Genes_2_Test, 'nobackup/neighboring_genes.Rds')
saveRDS(FBC_DF, 'nobackup/filtered_FBC_DF.Rds')
saveRDS(FBC_DF %>% filter(UMI_Count > 5), 'nobackup/filtered_FBC_DF_umi5.Rds')

# saving this a as a data frame that includes
# partitions of guides so that we can easily parallelize.
n_jobs <- 300 # this will then expect to split things over 300 jobs
partitioned_gRNAs <- data.frame(guide_ID=gRNA_Test_IDs) %>%
    mutate(a_partition=cut_number(n=n_jobs, row_number(), labels=c(1:n_jobs)))

saveRDS(partitioned_gRNAs, 'nobackup/gRNA_test_ids.Rds')


#pre-filter and make a look-up table to speed things
#up.
ngene_lookup <- split(All_1MB_Neighbouring_Genes_2_Test, All_1MB_Neighbouring_Genes_2_Test$Test_ID)
saveRDS(ngene_lookup, 'nobackup/ngene_lookup.Rds')

# ONLY RUN THIS WHEN SNAKEMAKE IS FINISHED. This aggregates all the output after snakemake has run
af <- list.files('nobackup/', pattern='*_results.txt')
all <- rbindlist(lapply(af, function(p) { fread(paste0('nobackup/', p)) }))
all %>% write.table("aggregated_results.txt", quote=F, row.names=F, sep='\t')
saveRDS(all, "aggregated_results.Rds")
