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

##Reading in data and creating a seaurat object

data_dir <- '/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/'
list.files(data_dir)

CRISPRaQTL_Pilot_data <- Read10X(data.dir = data_dir)
CRISPRaQTL_Pilot_Seurat_Object = CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`, min.cells = 100, min.features = 200)

CRISPRaQTL_Pilot_Seurat_Object

##Checking how much mito DNA we have in the Seurat object

CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(CRISPRaQTL_Pilot_Seurat_Object, pattern = "^MT-")

##Filtering for high quality cells 

CRISPRaQTL_Pilot_Seurat_Object_Subset <- subset(CRISPRaQTL_Pilot_Seurat_Object, subset = percent.mt < 17 & nCount_RNA > 1500)

CRISPRaQTL_Pilot_Seurat_Object_Subset

##Reading in Feature Barcode Matrix

matrix_dir = "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

##Indexing the matrix for only good cells in the seurat object

x <- rownames(CRISPRaQTL_Pilot_Seurat_Object_Subset@meta.data)

Seurat_Filtered_FBC_Mat <- mat[ ,x]

head(Seurat_Filtered_FBC_Mat)

##Creating a dataframe version of the full filtered FBC 

transposed_mat <- t(Seurat_Filtered_FBC_Mat)

FBC_DF <- tidy(transposed_mat)  %>% 
    rename(Cell = row) %>% 
    rename(Gene = column) %>% 
    rename(UMI_Count = value)

head(FBC_DF)

##Creating a gRNA Matrix 

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id)
    
x <- Ft_Ref$gRNA

gRNA_Mat <- Seurat_Filtered_FBC_Mat[x, ]

gRNA_Mat

##Creating a gRNA UMI DF for visualization

gRNA_Mat <- t(as.matrix(gRNA_Mat))

gRNA_DF <- as.data.frame(gRNA_Mat)

gRNA_DF <- rownames_to_column(gRNA_DF, "Cell")

gRNA_DF

##Calculating MOI

MOI_Mat <- rowSums(gRNA_Mat > 5)

MOI_DF <- as.data.frame(as.table(MOI_Mat))

MOI_DF <- MOI_DF %>%  
    rename(Cell = Var1) %>% 
    rename(MOI = Freq)  

summary(MOI_DF)

##Calculate assignment rate (the proportion of cells with 1 or more gRNA)

Assignment_Rate <- 1 - (sum(MOI_DF == 0 ) / length(MOI_DF$MOI))

Assignment_Rate


##Creating an MOI plot

options(repr.plot.width=7, repr.plot.height=6)

ggplot(MOI_DF, aes(x = MOI)) +
  geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 1) +  
  theme_classic() +
  xlim(-1,50) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) 


##Same plot as above without text labels for final figure

options(repr.plot.width=7, repr.plot.height=6)

ggplot(MOI_DF, aes(x = MOI)) +
  geom_histogram(colour = "#D1D4D4", fill = "#D1D4D4", binwidth = 1) +  
  theme_classic() +
  xlim(-1,50) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text = element_blank()) +
  labs(title = "", x = "", y = "") 


##Calculating the number of cells with each perturbation

N_Cells_With_Pert_Mat <- colSums(gRNA_Mat > 5)

N_Cells_With_Pert_DF <- as.data.frame(as.table(N_Cells_With_Pert_Mat))

N_Cells_With_Pert_DF <- N_Cells_With_Pert_DF %>%  
    rename(Gene = Var1) %>% 
    rename(N_Cells_With_Perturbation = Freq)  

summary(N_Cells_With_Pert_DF)

N_Cells_With_Pert_DF

##Then making a plot of the number of cells per perturbation

options(repr.plot.width=7, repr.plot.height=6)

ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
  geom_histogram(colour = "#56B4E9", fill = "#56B4E9", binwidth = 30) +  
  theme_classic() +
  xlim(-10,2500) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) 

##Same plot as above without axis text for final figure
options(repr.plot.width=7, repr.plot.height=6)


ggplot(N_Cells_With_Pert_DF, aes(x = N_Cells_With_Perturbation)) +
  geom_histogram(colour = "#D1D4D4", fill = "#D1D4D4", binwidth = 30) +  
  theme_classic() +
  xlim(-10,2500) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") 

