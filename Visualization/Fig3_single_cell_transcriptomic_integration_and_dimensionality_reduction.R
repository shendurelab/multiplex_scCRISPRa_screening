## Single-cell transcriptomic integration and dimensionality reduction of CRISPRa dataset and Treutlein timecourse sc dataset

# In this script, I make two Seurat objects for the two different datasets from Cellranger output files, and then integrate them using the Seurat v3 integration feature

#Required libraries

library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)

# read in Treutlein single-cell data, from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10632
# downloaded the matrices_timecourse.tar.gz files
# original names: 'counts.mtx.gz''features.tsv.gz''meta.tsv.gz', renamed to matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz to be compatible with the Seurat pipeline
data_dir = '/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/Treutlein_single_cell_data/'
list.files(data_dir)
timecourse_data <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1)
timecourse_seurat_object = CreateSeuratObject(counts = timecourse_data, min.cells = 4, min.features = 200)
timecourse_seurat_object #29554 cells, 17525 genes

##Visualize percent MT vs. n genes in Treutlein dataset (which has already been QC filtered)

timecourse_seurat_object[["percent.mt"]] <- PercentageFeatureSet(timecourse_seurat_object, pattern = "^MT-")

options(repr.plot.width=10, repr.plot.height=8)

plot1 <- FeatureScatter(timecourse_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = c("#56B4E9"))
plot2 <- FeatureScatter(timecourse_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = c("#56B4E9"))
plot1 + 
  theme(text = element_text(family="Arial", colour = "black", size = 24)) + theme(axis.text = element_text(family="Arial", colour = "black", size = 24))+
  geom_vline(xintercept = 1000)+
  geom_hline(yintercept=10)+ geom_hex(bins=100)+ggtitle("Treutlein dataset")

plot2 + 
  theme(text = element_text(family="Arial", colour = "black", size = 24)) + theme(axis.text = element_text(family="Arial", colour = "black", size = 24))+
  geom_vline(xintercept = 1000)+
  geom_hex(bins=100)

#add in meta data provided by Treutlein lab (same source website as listed above)
meta_data = read.table("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/Treutlein_single_cell_data/metadata.tsv",
                       row.names = 1)

#just making seurat object name shorter for simplicity
Treutlein_obj = timecourse_seurat_object
Treutlein_obj@meta.data = cbind(Treutlein_obj@meta.data, meta_data)

#reading in the CRISPRa dataset
data_dir <- '/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/'
list.files(data_dir)

QTL_neurons_data <- Read10X(data.dir = data_dir)
QTL_neurons = CreateSeuratObject(counts = QTL_neurons_data$`Gene Expression`, min.cells = 4, min.features = 200)
QTL_neurons

QTL_neurons[["percent.mt"]] <- PercentageFeatureSet(QTL_neurons, pattern = "^MT-")
#QC cutoffs of <17% mitochondrial reads and >1500 UMIs/cell
QTL_neurons <- subset(QTL_neurons, subset = percent.mt < 17 & nCount_RNA > 1500)
plot1 <- FeatureScatter(QTL_neurons, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = c("#56B4E9"))
plot2 <- FeatureScatter(QTL_neurons, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = c("#56B4E9"))
plot1 + 
  theme(text = element_text(family="Arial", colour = "black", size = 24)) + theme(axis.text = element_text(family="Arial", colour = "black", size = 24))+
  geom_vline(xintercept = 1500)+
  geom_hline(yintercept=17)+ geom_hex(bins=100)+ggtitle("CRISPRa dataset")

plot2 + 
  theme(text = element_text(family="Arial", colour = "black", size = 24)) + theme(axis.text = element_text(family="Arial", colour = "black", size = 24))+
  geom_vline(xintercept = 1500)+
  geom_hex(bins=100)

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

##Indexing the matrix for only high quality cells in the seurat object (x is the index vector of high quality cell names)
x <- rownames(QTL_neurons@meta.data)
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
Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
  rename(gRNA = id)
x <- Ft_Ref$gRNA
gRNA_Mat <- Seurat_Filtered_FBC_Mat[x, ]
head(gRNA_Mat)

gRNA_names = rownames(tail(Seurat_Filtered_FBC_Mat, n = 493))
head(gRNA_names)
length(gRNA_names)

gRNA_matrix_only = tail(Seurat_Filtered_FBC_Mat, n = 493)
head(gRNA_matrix_only)

for (num in 1:length(gRNA_names)){ 
  print(gRNA_names[num])
  gRNA_ID = gRNA_names[num]
  new_column = as.vector(gRNA_matrix_only[gRNA_ID,])
  QTL_neurons[[paste0(gRNA_ID)]] = new_column
}

# want to filter out cells that have fewer than 5 gRNA UMIs
QTL_neurons@meta.data[,c(5:497)][QTL_neurons@meta.data[,c(5:497)] <6 ] <- 0

#checking that cells with less than 5 gRNAs are filtered out
ggplot() + 
  geom_histogram(data = QTL_neurons@meta.data, 
                 aes(x=SCN2A_241_promoter	), bins = 100) +
  xlim(1,10)

#when we integrate the data, we only want to use genes that are in both datasets 
genes_intersect = intersect(rownames(Treutlein_obj), rownames(QTL_neurons))
Treutlein_obj_subset <- Treutlein_obj[rownames(Treutlein_obj) %in% genes_intersect, ] 
QTL_neurons_subset <- QTL_neurons[rownames(QTL_neurons) %in% genes_intersect, ] 
Treutlein_obj_subset
QTL_neurons_subset

#We are really only interested in d1, d2, d5, d14, d28, and d35 of the Treutlein dataset, so we will omit the other timepoints
Treutlein_cells_to_use = Treutlein_obj_subset@meta.data[(Treutlein_obj_subset@meta.data$timepoint != "h0" &  
                                                           Treutlein_obj_subset@meta.data$timepoint != "h6/12"),]

Treutlein_cells = rownames(Treutlein_cells_to_use)
Treutlein_subset_final <- subset(Treutlein_obj_subset, cells = Treutlein_cells)
Treutlein_subset = Treutlein_subset_final
table(Treutlein_subset@meta.data$timepoint)

#Want to downsample the CRISPRa dataset so it doesn't overtake the entire dataset
random_5000_cells = sample(rownames(QTL_neurons_subset@meta.data), 5000)
print(head(random_5000_cells))

cells.use <- intersect(colnames(QTL_neurons_subset), random_5000_cells)
QTL_neurons_subset_final <- subset(QTL_neurons_subset, cells = cells.use)
QTL_neurons_subset_final

#Preprocess CRISPRa data
QTL_neurons_subset_final <- NormalizeData(QTL_neurons_subset_final, normalization.method = "LogNormalize", scale.factor = 10000)
QTL_neurons_subset_final <- FindVariableFeatures(QTL_neurons_subset_final, selection.method = "vst", nfeatures = 2000)

#Preprocess Treutlein data
Treutlein_obj_subset <- NormalizeData(Treutlein_subset, normalization.method = "LogNormalize", scale.factor = 10000)
Treutlein_obj_subset <- FindVariableFeatures(Treutlein_obj_subset, selection.method = "vst", nfeatures = 2000)

#Seurat v3 integration (from here: https://www.cell.com/action/showPdf?pii=S0092-8674%2819%2930559-8)
#Note: this next step takes a long time (depending on your compute environment, can take anywhere from 10 minutes+)
features <- SelectIntegrationFeatures(object.list = list(Treutlein_obj_subset, QTL_neurons_subset_final))
neuron.anchors <- FindIntegrationAnchors(object.list = list(Treutlein_obj_subset, QTL_neurons_subset_final), 
                                         anchor.features = features)
neuron.combined <- IntegrateData(anchorset = neuron.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(neuron.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
neuron.combined <- ScaleData(neuron.combined, verbose = FALSE)
neuron.combined <- RunPCA(neuron.combined, npcs = 30, verbose = FALSE)
neuron.combined <- RunUMAP(neuron.combined, reduction = "pca", dims = 1:30)

#At this point, it's a good idea to save this seurat object so you don't need to do all the preprocessing steps again
saveRDS(neuron.combined,
        "path_to_where_you_want_to_save_here/processed_seurat_object.rds")

table(neuron.combined@meta.data$timepoint)
neuron.combined@meta.data$timepoint[is.na(neuron.combined@meta.data$timepoint)] <- "d19_QTL"
table(neuron.combined@meta.data$timepoint)

neuron.combined$timepoint <- factor(neuron.combined$timepoint, 
                                    levels = c("d1", "d2", "d5", "w2", "d19_QTL", "w4", "w5"))

p1 <- DimPlot(neuron.combined, reduction = "umap", group.by = "timepoint", size = 0.01)
p2 <- DimPlot(neuron.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1

#getting UMAP coordinates to facilitate nicer visualization of the UMAP
UMAP_coord = as.data.frame(neuron.combined[["umap"]]@cell.embeddings)
UMAP_coord_with_metadata = cbind(UMAP_coord, neuron.combined@meta.data)

#head(UMAP_coord_with_metadata)

UMAP_coord_with_metadata$day <- with(UMAP_coord_with_metadata, ifelse(timepoint == "d1", 1,
                                                                      ifelse(timepoint == "d2", 2,
                                                                             ifelse(timepoint == "d2", 2,
                                                                                    ifelse(timepoint == "d5", 5,
                                                                                           ifelse(timepoint == "w2", 14,
                                                                                                  ifelse(timepoint == "w4", 28,
                                                                                                         ifelse(timepoint == "w5", 35, "d19_QTL")))))))) #day16 should really be d19 to align to Treutlein timescale

UMAP_coord_with_metadata_treutlein = UMAP_coord_with_metadata[UMAP_coord_with_metadata$day != "d19_QTL",]
UMAP_coord_with_metadata_treutlein$day = as.numeric(UMAP_coord_with_metadata_treutlein$day)
UMAP_coord_with_metadata_QTL = UMAP_coord_with_metadata[UMAP_coord_with_metadata$day == "d19_QTL",]

UMAP_coord_with_metadata_treutlein$day = as.factor(UMAP_coord_with_metadata_treutlein$day)

#This next function is how we made the plot for Figure 2d
ggplot() +
  geom_point(data =
               UMAP_coord_with_metadata_treutlein,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_with_metadata_treutlein,
             aes(x = UMAP_1,
                 y = UMAP_2,
                 color = day),
             stroke = 0,
             size = 0.5)+
  theme_void()+
  scale_colour_grey(start = 1, end = 0)+
  geom_point(data =
               UMAP_coord_with_metadata_QTL,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_with_metadata_QTL,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "#56B4E9",
             stroke = 0,
             size = 0.5)








