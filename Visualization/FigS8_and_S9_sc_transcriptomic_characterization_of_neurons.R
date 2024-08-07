# Figure S7 | Single-cell transcriptomic characterization of iPSC-derived neurons in the neuron screen.

# A) Expression feature plots of canonical pluripotency markers NANOG, POU5F1, KLF4, FBXO15, and PODXL.
# B) Expression feature plots of pan-neuronal markers MAP2, RBFOX3, MAPT, ANK3, and NCAM1.
# C) Expression feature plots of central nervous system marker genes LHX9, GPM6A, and POU4F1.
# D) Expression feature plots of peripheral nervous system marker genes PHOX2B and PRPH.
# E) Expression feature plots of cortical excitatory neuron markers HOMER1, CUX1, and SLC17A7.
# F) Expression feature plots of GABAergic neuron marker genes GAD1 and GAD2.

## Load required libraries

library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)

## Read in combined dataset containing both CRISPRa data and Treutlein neuron timecourse data
## This dataset is already QC filtered, processed, and ready for visualization

neuron.combined = readRDS("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/Seurat_object_for_sc_data_visualization/processed_seurat_obj_with_gRNA_data.rds")
neuron.combined

## Subset combined dataset to look at individual datasets

neuron.combined_QTL = neuron.combined[,neuron.combined@meta.data$timepoint == "d19_QTL"]
neuron.combined_Treutlein = neuron.combined[,neuron.combined@meta.data$timepoint != "d19_QTL"]
neuron.combined_QTL
neuron.combined_Treutlein

## Get the UMAP coordinates as data points
## this is so we can plot the UMAP as a geom_point() instead of using a built in Seurat visualization function (nicer for visualization)

UMAP_coord_QTL = as.data.frame(neuron.combined_QTL[["umap"]]@cell.embeddings)
UMAP_coord_QTL_with_metadata = cbind(UMAP_coord_QTL, neuron.combined_QTL@meta.data)

UMAP_coord_Treutlein = as.data.frame(neuron.combined_Treutlein[["umap"]]@cell.embeddings)
UMAP_coord_Treutlein_with_metadata = cbind(UMAP_coord_Treutlein, neuron.combined_Treutlein@meta.data)

## Changing timepoint to a numerical day value 

table(UMAP_coord_Treutlein_with_metadata$timepoint)
UMAP_coord_Treutlein_with_metadata$day <- with(UMAP_coord_Treutlein_with_metadata, ifelse(timepoint == "d1", 1,
                                                                                          ifelse(timepoint == "d2", 2,
                                                                                                 ifelse(timepoint == "d2", 2,
                                                                                                        ifelse(timepoint == "d5", 5,
                                                                                                               ifelse(timepoint == "w2", 14,
                                                                                                                      ifelse(timepoint == "w4", 28,
                                                                                                                             ifelse(timepoint == "w5", 35, "other"))))))))

UMAP_coord_Treutlein_with_metadata$day = as.numeric(UMAP_coord_Treutlein_with_metadata$day)

## Can see how many cells there are per timepoint
table(UMAP_coord_Treutlein_with_metadata$day)

## Plotting the UMAP (both datasets) as a geom_point plot - nice for visualization

UMAP_coord_Treutlein_with_metadata$day = as.factor(UMAP_coord_Treutlein_with_metadata$day)

ggplot() +
  geom_point(data =
               UMAP_coord_Treutlein_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_Treutlein_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2,
                 color = day),
             stroke = 0,
             size = 0.5)+
  theme_void()+
  scale_colour_grey(start = 1, end = 0)+
  # a second geom_point to overlay our data on top of the Treutlein data
  geom_point(data =
               UMAP_coord_QTL_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_QTL_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "#56B4E9",
             stroke = 0,
             size = 0.5)+
  xlim(11,-11)

# without legend for actual paper figure
ggplot() +
  geom_point(data =
               UMAP_coord_Treutlein_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_Treutlein_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2,
                 color = day),
             stroke = 0,
             size = 0.5)+
  theme_void()+
  scale_colour_grey(start = 1, end = 0)+
  geom_point(data =
               UMAP_coord_QTL_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_QTL_with_metadata,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "#56B4E9",
             stroke = 0,
             size = 0.5)+
  xlim(11,-11)+ theme(legend.position = "none")

#save png of figure at dpi=300
#ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_QTL_on_Treutlein_data.png", 
#dpi=300)

## making a dataframe with just the counts of each gene, and transposing it to have cells as rows, and genes as columns

counts.df <- as.data.frame(neuron.combined_QTL@assays$RNA@counts)
counts.df_transpose = t(counts.df)

## adding in all the metadata to the dataframe we just made above so we can plot this data in a geom_point

UMAP_coord_QTL_with_RNA_counts = cbind(UMAP_coord_QTL_with_metadata, counts.df_transpose)
head(UMAP_coord_QTL_with_RNA_counts)

## Now we can plot the expression of genes any way we want using geom_point
## but I actually think the Seurat FeaturePlot() function is better for visualization (better contrast)
## FeaturePlot() visualizations are in the next cell and onwards
ggplot() +
  geom_point(data =
               UMAP_coord_QTL_with_RNA_counts,
             aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data =
               UMAP_coord_QTL_with_RNA_counts,
             aes(x = UMAP_1,
                 y = UMAP_2,
                 color = PRPH), #IMPORTANT: this is where you pick the gene to look at
             stroke = 0,
             size = 0.5)+
  theme_void()+
  scale_color_gradient(low = "black", high = "#56B4E9", 
                       limit = c(-1,1), space = "Lab", 
                       name="Normalized\nExpression")+
  xlim(11,-11) #reversing the x axis scale is necessary to flip the UMAP 180 degrees

# The FeaturePlot() function is what I use for all the supplemental figure plots, all plots made below. 
# A) Expression feature plots of canonical pluripotency markers NANOG, POU5F1, KLF4, FBXO15, and PODXL.

pluripotency_markers = c("NANOG", "POU5F1", "KLF4", "FBXO15", "PODXL")

for (marker in pluripotency_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

# B) Expression feature plots of pan-neuronal markers MAP2, RBFOX3, MAPT, ANK3, and NCAM1.

pan_neuronal_markers = c("MAP2", "RBFOX3", "MAPT", "ANK3", "NCAM1")

for (marker in pan_neuronal_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

# C) Expression feature plots of central nervous system marker genes LHX9, GPM6A, and POU4F1.

CNS_markers = c("LHX9", "GPM6A", "POU4F1")

for (marker in CNS_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

# D) Expression feature plots of peripheral nervous system marker genes PHOX2B and PRPH.

PNS_markers = c("PHOX2B", "PRPH")

for (marker in PNS_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

# E) Expression feature plots of cortical excitatory neuron markers HOMER1, CUX1, and SLC17A7.

CEN_markers = c("HOMER1", "CUX1", "SLC17A7")

for (marker in CEN_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

# F) Expression feature plots of GABAergic neuron marker genes GAD1 and GAD2. 

GABAergic_markers = c("GAD1", "GAD2")

for (marker in GABAergic_markers){
  print(c("current marker:", marker))
  print(FeaturePlot(neuron.combined_QTL, features = c(marker), 
                    pt.size = 1, 
                    cols = c("#E6E6E6", "black"),
                    order = TRUE)+
          xlim(11,-11)+
          theme_void()+
          labs(title=""))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/UMAP_marker_gene_plots/", marker, ".png"), 
  #dpi=300)
}

#### Figure S8 | Distribution of CRISPRa gRNAs in single-cell neuron transcriptome data.

# A) Cells harboring specific CRISPRa gRNAs (dark blue) overlaid onto the NGN2-induced neuron differentiation transcriptome data from Lin et al. No readily apparent spatial enrichment of gRNAs is observed in UMAP plots. Note that the CRISPRa dataset was randomly downsampled to 5000 cells for all UMAP comparison analyses.

#read in the single-cell data as well as the gRNA identities

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
x <- rownames(neuron.combined_QTL@meta.data)
Seurat_Filtered_FBC_Mat <- mat[ ,x]

##Creating a dataframe version of the full filtered FBC 
transposed_mat <- t(Seurat_Filtered_FBC_Mat)
FBC_DF <- tidy(transposed_mat)  %>% 
  rename(Cell = row) %>% 
  rename(Gene = column) %>% 
  rename(UMI_Count = value)

##Creating a gRNA Matrix 
Ft_Ref <- read_csv("/net/shendure/vol10/projects/ResQTL/nobackup/10X_single_cell_pilot_expt_Nov_2021/10X_pilot_1_cellranger_count_SNV_fixed/outs/crispr_analysis/feature_reference.csv")  %>% 
  rename(gRNA = id)
x <- Ft_Ref$gRNA
gRNA_Mat <- Seurat_Filtered_FBC_Mat[x, ]

gRNA_names = rownames(tail(Seurat_Filtered_FBC_Mat, n = 493))
head(gRNA_names)
length(gRNA_names)

UMAP_coord_QTL = as.data.frame(neuron.combined_QTL[["umap"]]@cell.embeddings)
UMAP_coord_QTL_with_metadata = cbind(UMAP_coord_QTL, neuron.combined_QTL@meta.data)

UMAP_coord_Treutlein = as.data.frame(neuron.combined_Treutlein[["umap"]]@cell.embeddings)
UMAP_coord_Treutlein_with_metadata = cbind(UMAP_coord_Treutlein, neuron.combined_Treutlein@meta.data)

df_base = UMAP_coord_QTL_with_metadata[,c(1:18)]
df_base

## gRNA UMAPs on top of Treutlein data

for (num in 1:length(gRNA_names)){ #for all the gRNAs
  print(gRNA_names[num])
  current_gRNA = gRNA_names[num]
  df_OI = UMAP_coord_QTL_with_metadata[,c(current_gRNA)] 
  final_df = cbind(df_base, df_OI)
  print(sum(final_df$df_OI))
  
  UMAP_coord_QTL_with_metadata_no_zeros = final_df[df_OI != 0, ] #column 19 will be the gRNA of interest
  print(sum(UMAP_coord_QTL_with_metadata_no_zeros$df_OI))
  print(nrow(UMAP_coord_QTL_with_metadata_no_zeros))
  print(ggplot()+
          geom_point(data =
                       UMAP_coord_Treutlein_with_metadata,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "#CACACA",
                     stroke = 0,
                     size = 1) +
          geom_point(data =
                       UMAP_coord_Treutlein_with_metadata,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "#E6E6E6",
                     stroke = 0,
                     size = 0.75) +
          geom_point(data =
                       UMAP_coord_QTL_with_metadata_no_zeros,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "black",
                     stroke = 0,
                     size = 2.75)+
          geom_point(data =
                       UMAP_coord_QTL_with_metadata_no_zeros,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "steelblue",
                     stroke = 0,
                     size = 2.5)+
          ggtitle(gRNA_names[num])+
          theme_void()+ 
          labs(title = current_gRNA,
               subtitle = paste0(nrow(UMAP_coord_QTL_with_metadata_no_zeros), " gRNA cells shown"))+
          theme(plot.title = element_text(family="Arial", size = 28, face = "bold"),
                plot.subtitle = element_text(family="Arial", size = 26))+
          xlim(11,-11))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNAs_on_UMAPs/lighter_background/",gRNA_names[num],".png"),
  #dpi = 300,
  #height = 8, width = 8)
}

# looking at the neuron hits

neuron_hit_list = c("ZC3HAV1_393_TSS_pos_ctrl",
                    "TBR1_313_promoter",
                    "TBR1_316_Flashfry_promoter",
                    "TCF4_329_promoter",
                    "CCND2_54_TSS_pos_ctrl",
                    "TCF4_326_promoter",
                    "ANO5_31_TSS_pos_ctrl",
                    "DNMT3B_83_TSS_pos_ctrl",
                    "TCF4_328_promoter",
                    "FOXP1_91_promoter",
                    "TCF4_378_Flashfry_promoter",
                    "BCL11A_33_promoter",
                    "BCL11A_41_promoter",
                    "CCNE2_55_TSS_pos_ctrl",
                    "FOXP1_148_Flashfry_promoter",
                    "TCF4_330_promoter",
                    "TCF4_348_Flashfry_promoter")

for (num in 1:length(neuron_hit_list)){ 
  print(neuron_hit_list[num])
  current_gRNA = neuron_hit_list[num]
  df_OI = UMAP_coord_QTL_with_metadata[,c(current_gRNA)] 
  final_df = cbind(df_base, df_OI)
  print(sum(final_df$df_OI))
  
  UMAP_coord_QTL_with_metadata_no_zeros = final_df[df_OI != 0, ] #column 19 will be the gRNA of interest
  print(sum(UMAP_coord_QTL_with_metadata_no_zeros$df_OI))
  print(nrow(UMAP_coord_QTL_with_metadata_no_zeros))
  print(ggplot()+
          geom_point(data =
                       UMAP_coord_Treutlein_with_metadata,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "#E6E6E6",
                     stroke = 0,
                     size = 1) +
          geom_point(data =
                       UMAP_coord_Treutlein_with_metadata,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "#E6E6E6",
                     stroke = 0,
                     size = 0.75) +
          geom_point(data =
                       UMAP_coord_QTL_with_metadata_no_zeros,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "black",
                     stroke = 0,
                     size = 2.75)+
          geom_point(data =
                       UMAP_coord_QTL_with_metadata_no_zeros,
                     aes(x = UMAP_1,
                         y = UMAP_2),
                     color = "steelblue",
                     stroke = 0,
                     size = 2.5)+
          ggtitle(gRNA_names[num])+
          theme_void()+ 
          labs(title = current_gRNA,
               subtitle = paste0(nrow(UMAP_coord_QTL_with_metadata_no_zeros), " gRNA cells shown"))+
          theme(plot.title = element_text(family="Arial", size = 28, face = "bold"),
                plot.subtitle = element_text(family="Arial", size = 26))+
          xlim(11,-11))
  #ggsave(paste0("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNAs_on_UMAPs/neuron_hits/lighter_background/",neuron_hit_list[num],".png"),
  #dpi = 300,
  #height = 8, width = 8)
}







