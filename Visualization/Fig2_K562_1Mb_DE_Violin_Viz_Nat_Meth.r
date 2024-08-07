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
library(scales)

##Reading in data and creating a seaurat object

data_dir <- '/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/'
list.files(data_dir)

CRISPRaQTL_Pilot_data <- Read10X(data.dir = data_dir)
CRISPRaQTL_Pilot_Seurat_Object = CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`, min.cells = 4, min.features = 200)

CRISPRaQTL_Pilot_Seurat_Object

##Checking how much mito DNA we have in the Seurat object

CRISPRaQTL_Pilot_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(CRISPRaQTL_Pilot_Seurat_Object, pattern = "^MT-")

##Filtering for high quality cells 

CRISPRaQTL_Pilot_Seurat_Object_Subset <- subset(CRISPRaQTL_Pilot_Seurat_Object, subset = percent.mt < 10 & nCount_RNA > 4000)

CRISPRaQTL_Pilot_Seurat_Object_Subset

##Reading in Feature Barcode Matrix

matrix_dir = "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/"
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

##Normalize the data with seurats default paramters before statistical testing

CRISPRaQTL_Pilot_Seurat_Object_Subset <- NormalizeData(CRISPRaQTL_Pilot_Seurat_Object_Subset)

##Creating reference dataframes with gRNA-gene and gene-gRNA mappings

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id)
    

##Making violin for ANXA1 enhancer hit 

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("chr9.871_489_Gasperini_enhancer")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "ANXA1", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

Control_Cells_DS <- sample(Control_Cells, 178)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'ANXA1', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,4.25), breaks = c(0,1,2,3,4)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for second ANXA1 enhancer hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("chr9.871_490_Gasperini_enhancer")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "ANXA1", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'ANXA1', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,4.25), breaks = c(0,1,2,3,4)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for TSPAN5 enhancer hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("chr4.2290_458_Gasperini_enhancer")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "TSPAN5", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'TSPAN5', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,2.5), breaks = c(0,1,2)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for TMSB4X enhancer hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("chrX.232_495_Gasperini_enhancer")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "TMSB4X", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'TMSB4X', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 4.5), breaks = c(0,1,2,3,4)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for CCND2 TSS hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("CCND2_54_TSS_pos_ctrl")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "CCND2", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'CCND2', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 3.25), breaks = c(0,1,2,3)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for ANK2 promoter hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("ANK2_2_promoter")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "ANK2", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'ANK2', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#56B4E9", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 2.5), breaks = c(0,1,2)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for BCL11A promoter hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("BCL11A_41_promoter")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "BCL11A", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'BCL11A', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#56B4E9", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 2.5), breaks = c(0,1,2)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for FOXP1 enhancer hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("FOXP1_148_Flashfry_promoter")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "FOXP1", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'FOXP1', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 2.5), breaks = c(0,1,2)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

##Making violin for GNB2 TSS positive control hit gRNA

##isolate cells with guide 

gRNA_FBC_DF <- FBC_DF  %>% 
    filter(Gene %in% c("GNB2_180_TSS_pos_ctrl")) %>% 
    filter(UMI_Count > 5) %>% 
    select(Cell)

gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

    # This makes control cells all cells that dont have the targeting guide
        Control_Cells_DF <- FBC_DF %>% 
        filter(!Cell %in% gRNA_Cells)

    Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for DE testing

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"

##Test for differential expression for a gene of interest

Result <- FindMarkers(CRISPRaQTL_Pilot_Seurat_Object_Subset, ident.1 = "gRNA Cells", 
                 ident.2 = "Control Cells", min.pct = 0, features = "GNB2", logfc.threshold = 0)

print(paste0('finished testing'))
    print(head(Result))
    return(as.data.frame(Result) %>%
           mutate(n_cells=length(gRNA_Cells),
                  n_control_cells=length(Control_Cells)))


##Downsample the control cells to have the same number of cells as those that have the targeting gRNA for plotting

set.seed(26)

##Make a violin plot

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "gRNA Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells_DS) <- "Control Cells Downsample"


options(repr.plot.width=4, repr.plot.height=6)

p = VlnPlot(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, features = 'GNB2', idents = c("Control Cells Downsample", "gRNA Cells"),
       pt.size = 0)
p$layers[[1]]$aes_params$size = 0
p$layers[[1]]$aes_params$fill = NA

p + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0), limits = c(0.0001, 3.05), breaks = c(0,1,2,3)) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 0))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))
