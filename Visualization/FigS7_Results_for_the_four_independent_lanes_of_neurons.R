#### Figure S6 | Results for four independent 10X lanes from neuron screen.

# A) The four 10X lanes profiled consisted of two lanes with dCas9-VPH neurons that were sorted on the top 40% of GFP expression in these cells, and two lanes that were not on the top 40% of GFP expression in these cells. The cells that were not sorted were still 100% GFP+. Following QC and gRNA assignment we identified an average of 7.71, 7.91, 4.55, and 4.39 gRNAs/cell for the four different 10X lanes profiled (median 7.71, 7.91, 4.55, and 4.39 gRNAs per cell). Note that sorting neurons on the top 40% of GFP expression boosted the median and mean gRNAs/cell ~2 fold. Also note that piggyBac integrations per cell distribution is not well-modeled by a standard Poisson distribution and is better approximated by an exponential function.
# B) Multiplexing multiple perturbations per cell yielded an average of 218.7, 189.9, 118.1, and 114.4 cells/gRNA for the four different 10X lanes profiled (median 166, 146, 98, and 96 cells/gRNA).
# C) QQ-plots displaying observed vs. expected p-value distributions for targeting (blue) and NTC (downsampled) populations across the four different 10X lanes profiled.
# D) QQ-plots for targeting tests against their intended/programmed target (blue) compared to targeting tests of all other genes with 1Mb of each gRNA (pink) and NTCs (gray downsampled) across the four different 10X lanes profiled.
# E) Matrix correlation plot displaying the Pearson correlations of the log2(fold change) of target gene expression values for programmed targets across the four different 10X lanes profiled.
# F) Violin plot displaying the log2(fold change) of target gene expression values for programmed targets for neurons that were sorted on the top 40% GFP expression (sorted) and neurons that were not sorted (not sorted).

##Loading required libraries
library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)
library(reshape2)

#### A) gRNAs per cell and B) cells per gRNA for each of the four 10x lanes profiled
#Start by making a Seurat object, filtering out cells with high mitochondrial reads and low RNA counts, and then building appropriate matrices and dataframes to plot the data in the desired histogram format. 

#For these plots, the axes are different for the sorted vs. not sorted populations, so am splitting into two for loops

lane_IDs_sorted = c("lane_1_sorted", "lane_2_sorted")
lane_IDs_not_sorted = c("lane_3_not_sorted", "lane_4_not_sorted")

for (lane in lane_IDs_sorted){
  print(c("current lane being analyzed:", lane))
  data_dir = paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_count_output_", lane, "/outs/filtered_feature_bc_matrix/")
  CRISPRaQTL_Pilot_data = Read10X(data.dir = data_dir)
  CRISPRaQTL_Pilot_Seurat_Object = CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`)
  
  #Checking how much mito DNA we have in the Seurat object
  CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]] = PercentageFeatureSet(CRISPRaQTL_Pilot_Seurat_Object, pattern = "^MT-")
  CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]]
  CRISPRaQTL_Pilot_Seurat_Object_Subset <- subset(CRISPRaQTL_Pilot_Seurat_Object, subset = percent.mt < 17 & nCount_RNA > 1500)
  print(CRISPRaQTL_Pilot_Seurat_Object_Subset)
  
  #Reading in Feature Barcode Matrix
  
  barcode.path <- paste0(data_dir, "barcodes.tsv.gz")
  features.path <- paste0(data_dir, "features.tsv.gz")
  matrix.path <- paste0(data_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  #Indexing the matrix for only good cells in the seurat object
  x = rownames(CRISPRaQTL_Pilot_Seurat_Object_Subset@meta.data)
  Seurat_Filtered_FBC_Mat <- mat[,x]
  
  #Creating a dataframe version of the full filtered FBC 
  transposed_mat <- t(Seurat_Filtered_FBC_Mat)
  FBC_DF <- tidy(transposed_mat)  %>% 
    rename(Cell = row) %>% 
    rename(Gene = column) %>% 
    rename(UMI_Count = value)
  #Making a gRNA matrix; note this feature reference is the same for all samples/lanes
  Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id)
  x = Ft_Ref$gRNA
  gRNA_Mat <- Seurat_Filtered_FBC_Mat[x,]
  
  #Creating a gRNA UMI DF for visualization
  gRNA_Mat <- t(as.matrix(gRNA_Mat))
  gRNA_DF <- as.data.frame(gRNA_Mat)
  gRNA_DF <- rownames_to_column(gRNA_DF, "Cell")
  
  #Now we can calculate the single-cell gRNA capture rate (MOI)
  MOI_Mat <- rowSums(gRNA_Mat > 5)
  MOI_DF <- as.data.frame(as.table(MOI_Mat))
  MOI_DF <- MOI_DF %>%  
    rename(Cell = Var1) %>% 
    rename(MOI = Freq)  
  print(summary(MOI_DF))
  
  #Now we can calculate the assignment rate (proportion of cells with 1 or more gRNA)
  Assignment_Rate <- 1 - (sum(MOI_DF == 0 ) / length(MOI_DF$MOI))
  print(Assignment_Rate)
  
  #Creating an MOI plot
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(MOI_DF, aes(x = MOI)) +
          geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 1) +  
          theme_classic() +
          xlim(-1,41) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(axis.text = element_text(family="Arial", colour = "black", size = 24)))
  
  #Same plot as above without text labels for final figure
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(MOI_DF, aes(x = MOI)) +
          geom_histogram(colour = "#A7A9AC", fill = "#A7A9AC", binwidth = 1) +  
          theme_classic() +
          xlim(-1,41) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          theme(axis.text = element_blank()) +
          labs(title = "", x = "", y = "")) 
  
  #Calculating the number of cells with each perturbation
  N_Cells_With_Pert_Mat <- colSums(gRNA_Mat > 5)
  N_Cells_With_Pert_DF <- as.data.frame(as.table(N_Cells_With_Pert_Mat))
  N_Cells_With_Pert_DF <- N_Cells_With_Pert_DF %>%  
    rename(Gene = Var1) %>% 
    rename(N_Cells_With_Perturbation = Freq)  
  print(summary(N_Cells_With_Pert_DF))
  
  #Making a plot of the number of cells per perturbation
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
          geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 10) +  
          theme_classic() +
          xlim(-10,601) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(axis.text = element_text(family="Arial", colour = "black", size = 24)))
  
  #Same plot as above without axis text for final figure
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
          geom_histogram(colour = "#A7A9AC", fill = "#A7A9AC", binwidth = 10) +  
          theme_classic() +
          xlim(-10,601) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.text = element_blank()) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = ""))
}

#now for the two unsorted lanes
for (lane in lane_IDs_not_sorted){
  print(c("current lane being analyzed:", lane))
  data_dir = paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_count_output_", lane, "/outs/filtered_feature_bc_matrix/")
  CRISPRaQTL_Pilot_data = Read10X(data.dir = data_dir)
  CRISPRaQTL_Pilot_Seurat_Object = CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`)
  
  #Checking how much mito DNA we have in the Seurat object
  CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]] = PercentageFeatureSet(CRISPRaQTL_Pilot_Seurat_Object, pattern = "^MT-")
  CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]]
  CRISPRaQTL_Pilot_Seurat_Object_Subset <- subset(CRISPRaQTL_Pilot_Seurat_Object, subset = percent.mt < 17 & nCount_RNA > 1500)
  print(CRISPRaQTL_Pilot_Seurat_Object_Subset)
  
  #Reading in Feature Barcode Matrix
  
  barcode.path <- paste0(data_dir, "barcodes.tsv.gz")
  features.path <- paste0(data_dir, "features.tsv.gz")
  matrix.path <- paste0(data_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  #Indexing the matrix for only good cells in the seurat object
  x = rownames(CRISPRaQTL_Pilot_Seurat_Object_Subset@meta.data)
  Seurat_Filtered_FBC_Mat <- mat[,x]
  
  #Creating a dataframe version of the full filtered FBC 
  transposed_mat <- t(Seurat_Filtered_FBC_Mat)
  FBC_DF <- tidy(transposed_mat)  %>% 
    rename(Cell = row) %>% 
    rename(Gene = column) %>% 
    rename(UMI_Count = value)
  #Making a gRNA matrix; note this feature reference is the same for all samples/lanes
  Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id)
  x = Ft_Ref$gRNA
  gRNA_Mat <- Seurat_Filtered_FBC_Mat[x,]
  
  #Creating a gRNA UMI DF for visualization
  gRNA_Mat <- t(as.matrix(gRNA_Mat))
  gRNA_DF <- as.data.frame(gRNA_Mat)
  gRNA_DF <- rownames_to_column(gRNA_DF, "Cell")
  
  #Now we can calculate the single-cell gRNA capture rate (MOI)
  MOI_Mat <- rowSums(gRNA_Mat > 5)
  MOI_DF <- as.data.frame(as.table(MOI_Mat))
  MOI_DF <- MOI_DF %>%  
    rename(Cell = Var1) %>% 
    rename(MOI = Freq)  
  print(summary(MOI_DF))
  
  #Now we can calculate the assignment rate (proportion of cells with 1 or more gRNA)
  Assignment_Rate <- 1 - (sum(MOI_DF == 0 ) / length(MOI_DF$MOI))
  print(Assignment_Rate)
  
  #Creating an MOI plot
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(MOI_DF, aes(x = MOI)) +
          geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 1) +  
          theme_classic() +
          xlim(-1,31) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(axis.text = element_text(family="Arial", colour = "black", size = 24)))
  
  #Same plot as above without text labels for final figure
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(MOI_DF, aes(x = MOI)) +
          geom_histogram(colour = "#A7A9AC", fill = "#A7A9AC", binwidth = 1) +  
          theme_classic() +
          xlim(-1,31) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          theme(axis.text = element_blank()) +
          labs(title = "", x = "", y = "")) 
  
  #Calculating the number of cells with each perturbation
  N_Cells_With_Pert_Mat <- colSums(gRNA_Mat > 5)
  N_Cells_With_Pert_DF <- as.data.frame(as.table(N_Cells_With_Pert_Mat))
  N_Cells_With_Pert_DF <- N_Cells_With_Pert_DF %>%  
    rename(Gene = Var1) %>% 
    rename(N_Cells_With_Perturbation = Freq)  
  print(summary(N_Cells_With_Pert_DF))
  
  #Making a plot of the number of cells per perturbation
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
          geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 10) +  
          theme_classic() +
          xlim(-10,401) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(axis.text = element_text(family="Arial", colour = "black", size = 24)))
  
  #Same plot as above without axis text for final figure
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
          geom_histogram(colour = "#A7A9AC", fill = "#A7A9AC", binwidth = 10) +  
          theme_classic() +
          xlim(-10,401) +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
          theme(axis.text = element_blank()) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = ""))
}

#### C) QQ plots for all four lanes (left to right: sorted lane 1, sorted lane 2, not sorted lane 3, not sorted lane 4)

#For the QQ plot visualization, we start out by reading in the DE results which have the DE testing significance values we need

lanes = c("lane_1", 
          "lane_2", 
          "lane_3", 
          "lane_4")

for (lane in lanes){
  print(c("current lane being analyzed:", lane))
  Full_Results = read_tsv(paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/1MB_DE_aggregated_results_", lane, ".txt")) %>%
    mutate(plot_pval = p_val+2.225074e-308) %>% 
    mutate(plot_pval = -log10(plot_pval)) %>% 
    filter(pct.1 > 0.0019) %>%
    filter(pct.1 > 0.0019) %>%
    unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)
  
  #Read in detected neighbouring genes DF with additional metadata
  Detc_Neighbouring_Genes = readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
    unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
    distinct(Test_ID, .keep_all = TRUE) 
  Full_Results = Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")
  Targeting_Results = Full_Results %>% 
    filter(Targeting == "Targeting")
  set.seed(12)
  NTC_Results = Full_Results %>% 
    filter(Targeting == "Non-targeting") 
  NTC_Results %>% 
    arrange(p_val)
  
  #Randomly downsampling the NTCs then recombining for plotting
  NTC_Results <- NTC_Results %>%   
    slice_sample(n = length(Targeting_Results$p_val)) 
  DS_Full_Results <- rbind(NTC_Results, Targeting_Results) 
  
  #QQ plot comparing p-value distributions
  options(repr.plot.width=12, repr.plot.height=8)
  
  print(ggplot(DS_Full_Results, aes(sample = plot_pval, colour = Targeting)) +
          stat_qq(size = 2) +
          scale_color_manual(values=c("#999999", "#56B4E9")) +
          theme_classic() +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "-log10(expected P-values)", y = "-log10(observed P-values)") +
          theme(text = element_text(family="Arial", colour = "black", size = 30)) +
          theme(axis.text = element_text(colour = "black")) +
          xlim(0,4))
  
  ##same QQ plot as above without axis text for final figures
  options(repr.plot.width=12, repr.plot.height=8)
  print(ggplot(DS_Full_Results, aes(sample = plot_pval, colour = Targeting)) +
          stat_qq(size = 2) +
          scale_color_manual(values=c("#999999", "#56B4E9")) +
          theme_classic() +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(text = element_text(family="Arial", colour = "black", size = 30)) +
          theme(axis.text = element_blank()) +
          xlim(0,4)+
          theme(legend.position = "none"))
  
  ##Now making QQ plot with alternate target DE results
  Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
  Primary_Target_Tests = Ft_Ref$Test_ID
  Primary_Targeting_Results = Targeting_Results %>% 
    filter(Test_ID %in% Primary_Target_Tests)
  Primary_Targeting_Results$Primary_Targeting = "Primary_Target"
  Other_Targeting_Results = Targeting_Results %>% 
    filter(!Test_ID %in% Primary_Target_Tests)
  Other_Targeting_Results$Primary_Targeting <- "Alternate_Target"
  NTC_Results$Primary_Targeting = "Non-targeting"
  DS_Full_Results_Primary_Targets = rbind(NTC_Results, Other_Targeting_Results, Primary_Targeting_Results) 
  
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(DS_Full_Results_Primary_Targets, aes(sample = plot_pval, colour = Primary_Targeting)) +
          stat_qq(size = 2) +
          scale_color_manual(values=c("#EC008C", "#999999", "#56B4E9")) +
          theme_classic() +
          theme(axis.line = element_line(colour = 'black', size = 0.4)) +
          theme(axis.ticks.y = element_line(colour = "black", size = 0.4)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "-log10(expected P-values)", y = "-log10(observed P-values)") +
          theme(text = element_text(family="Arial", colour = "black", size = 30)) +
          theme(axis.text = element_text(colour = "black")) +
          xlim(0,4))
  
  ##Same QQ plot as above without axis text for final figures
  options(repr.plot.width=7, repr.plot.height=6)
  print(ggplot(DS_Full_Results_Primary_Targets, aes(sample = plot_pval, colour = Primary_Targeting)) +
          stat_qq(size = 2) +
          scale_color_manual(values=c("#EC008C", "#999999", "#56B4E9")) +
          theme_classic() +
          theme(axis.line = element_line(colour = 'black', size = 0.7)) +
          theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
          theme(axis.ticks.length=unit(.25, "cm")) +
          labs(title = "", x = "", y = "") +
          theme(legend.position = "none") +
          theme(text = element_text(family="Arial", colour = "black", size = 30)) +
          theme(axis.text = element_blank()) +
          xlim(0,4))
}

#### E) Matrix correlation plot displaying the Pearson correlations of the log2(fold change) of target gene expression values for programmed targets across the four different 10X lanes profiled.

##Read in detected neighbouring genes DF with additional metadata
Detc_Neighbouring_Genes = readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
  unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 

lanes = c("lane_1", 
          "lane_2", 
          "lane_3", 
          "lane_4")

for (i in 1:length(lanes)){
  Full_Results <- read_tsv(paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/1MB_DE_aggregated_results_", lanes[i], ".txt"))%>%
    mutate(plot_pval = p_val+2.225074e-308) %>% 
    mutate(plot_pval = -log10(plot_pval))  %>% 
    filter(pct.1 > 0.0019) %>%
    filter(pct.1 > 0.0019) %>%
    unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)
  Full_Results <- Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")
  Targeting_Results <- Full_Results %>% 
    filter(Targeting == "Targeting")
  
  Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
  Primary_Target_Tests <- Ft_Ref$Test_ID
  Primary_Targeting_Results <- Targeting_Results %>% 
    filter(Test_ID %in% Primary_Target_Tests)
  Primary_Targeting_Results$Primary_Targeting <- "Primary_Target"
  #paste0("desired_cols_", lane) = Primary_Targeting_Results[,c(2,6)]
  assign(paste0("desired_columns_", lanes[i]), Primary_Targeting_Results[,c(2,6)])
}

# Put results in one dataframe
log2fc_df = merge(desired_columns_lane_1, desired_columns_lane_2, by = "Test_ID")
log2fc_df = merge(log2fc_df, desired_columns_lane_3, by = "Test_ID")
log2fc_df = merge(log2fc_df, desired_columns_lane_4, by = "Test_ID")
colnames(log2fc_df) = c("Test_ID", "sorted_lane_1", "sorted_lane_2", "not_sorted_lane_1", "not_sorted_lane_2")
rownames(log2fc_df) = log2fc_df[,1]
log2fc_df = log2fc_df[,-c(1)]

cormat <- round(cor(log2fc_df),2)
cormat
melted_cormat <- melt(cormat)
head(melted_cormat)


# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Correlation matrix visualization
ggplot(data = melted_cormat, aes(Var1, -c(Var2), fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient(low = "black", high = "#56B4E9", 
                      limit = c(-1,1), space = "Lab", 
                      name="Lane\nCorrelation") +
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+
  coord_fixed()


ggplot(data = melted_cormat, aes(Var1, -c(Var2), fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient(low = "black", high = "#56B4E9", 
                      limit = c(-1,1), space = "Lab", 
                      name="Lane\nCorrelation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 20, hjust = 1), 
        axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 20, hjust = 1),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+
  coord_fixed()

#### F) Violin plot displaying the log2(fold change) of target gene expression values for programmed targets for neurons that were sorted prior to expression profiling and not sorted prior to expression profiling (lanes 1&2 vs. lanes 3&4)

##Read in detected neighbouring genes DF with additional metadata
Detc_Neighbouring_Genes = readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
  unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 

lanes = c("lane_1", 
          "lane_2", 
          "lane_3", 
          "lane_4")

results_list = list()

for (i in 1:length(lanes)){
  Full_Results <- read_tsv(paste0("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/1MB_DE_aggregated_results_", lanes[i], ".txt")) %>%
    mutate(plot_pval = p_val+2.225074e-308) %>% 
    mutate(plot_pval = -log10(plot_pval))  %>% 
    filter(pct.1 > 0.0019) %>%
    filter(pct.1 > 0.0019) %>%
    unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)
  Full_Results <- Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")
  Targeting_Results <- Full_Results %>% 
    filter(Targeting == "Targeting")
  
  Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
  Primary_Target_Tests <- Ft_Ref$Test_ID
  Primary_Targeting_Results <- Targeting_Results %>% 
    filter(Test_ID %in% Primary_Target_Tests)
  Primary_Targeting_Results$Primary_Targeting <- "Primary_Target"
  results_list[[i]] = Primary_Targeting_Results[,c(2)]
}

lane_1_results = results_list[[1]]
lane_2_results = results_list[[2]]
lane_3_results = results_list[[3]]
lane_4_results = results_list[[4]]

lane_1_results$sort = "sorted"
lane_2_results$sort = "sorted"
lane_3_results$sort = "not_sorted"
lane_4_results$sort = "not_sorted"

log2fc_df = rbind(lane_1_results, 
                  lane_2_results, 
                  lane_3_results, 
                  lane_4_results)

print(c(head(log2fc_df), tail(log2fc_df)))

# violin plot, first with labels, then without for paper figure
ggplot(log2fc_df, aes(x=sort, y = avg_log2FC, fill = sort)) + 
  geom_violin() +
  theme_classic()+
  scale_fill_manual(values = c("#56B4E9","steelblue4"))

ggplot(log2fc_df, aes(x=sort, y = avg_log2FC, fill = sort)) + 
  geom_violin() +
  theme_classic()+
  scale_fill_manual(values = c("#56B4E9","steelblue4"))+
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank())+
  coord_fixed()+
  theme(legend.position="none")















