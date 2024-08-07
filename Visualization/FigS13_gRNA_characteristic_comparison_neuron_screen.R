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

# Data frame objects and lists I need for this analysis: 
# 1) Neuron promoter hits (list)
# 2) Neuron aggregated results 1 Mb (data frame)
# 3) Neuron aggregated results, promoter hits only (data frame)

#### 1. Do hit gRNAs biasly target highly or lowly expressed genes?
#all non-hit gRNAs contain the data for genes that were targeted in the screen, but weren't hits, which I think is correct in terms of what we want to look at 

#filter all DE testing results for all the target genes in the neuron hits df;
#that way, should have all gRNAs against all intended and unintended targets (just a few of these):
Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/cellranger_input_files/feature_reference.csv")  %>% 
  rename(target_guide = id) %>%
  rename(target_gene = target_gene_name)

neuron_results = read.table("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/aggregated_results.txt",
                            sep = "\t",
                            header = TRUE)
#tail(neuron_results)
#neuron_hits
neuron_hits = read.table("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Hits.csv",
                         sep = ",",
                         header = TRUE)
neuron_hit_list = neuron_hits$target_guide
length(neuron_hit_list)
#head(neuron_hits)

#Filter aggregated results by target genes in neuron_hits and get rid of NTCs (don't want for this analysis)
#Exceptions for neuron: chr6.2134_473_Gasperini_NTC (enhancer) and the two ASIC1 enhancers, CALM3, RPS18 and WWC3 (promoters)
#neuron_results_targets_only = neuron_results %>%
#filter(target_gene %in% neuron_hits$target_gene | target_gene %in% Ft_Ref$target_gene_name) #this is so I can get all targets
#neuron_results_targets_only = neuron_results %>%
#    filter(target_gene %in% Ft_Ref$target_gene | target_guide %in% Ft_Ref$target_guide) #this is so I can get all the intended targets, will add back unintended targets later

#merge keeps all the matching gRNA/target_gene pairs with the DE results data
neuron_results_targets_only = merge(neuron_results, Ft_Ref) #NTCs get filtered out during this step
#neuron_results_targets_only

#filter out the exceptions: CALM3
CALM3_hit = neuron_results %>%
  filter(grepl('PPP5D1_238_TSS_pos_ctrl', target_guide), grepl("CALM3", target_gene))
promoter_exceptions = CALM3_hit
#promoter_exceptions

#enhancer_exceptions
neuron_results_targets_only = neuron_results_targets_only[,c(3,4,5,6,7,1,2,8,9)]
#head(neuron_results_targets_only)
#head(promoter_exceptions)
#head(enhancer_exceptions)
#adding the promoter and enhancer exceptions back
neuron_results_all_targets = rbind(neuron_results_targets_only, promoter_exceptions)
neuron_results_all_targets

#Now want to add in if each row is a hit or not - can get this information from neuron_hits
neuron_results_all_targets$hit[neuron_results_all_targets$target_guide %in% neuron_hit_list] = "yes"
neuron_results_all_targets$hit[!neuron_results_all_targets$target_guide %in% neuron_hit_list] = "no"
table(neuron_results_all_targets$hit)

#Need to manually fix the exceptions (chr12.1559_411_Gasperini_enhancer, chr12.1559_412_Gasperini_enhancer,
#PPP5D1_238_TSS_pos_ctrl, SYNGAP1_295_Flashfry_promoter, FOXP1_146_Flashfry_promoter)
neuron_results_all_targets$hit[neuron_results_all_targets$target_guide == "PPP5D1_238_TSS_pos_ctrl" &
                                 neuron_results_all_targets$target_gene == "PPP5D1"] = "no"

table(neuron_results_all_targets$hit) #this is correct now
head(neuron_results_all_targets)
#neuron_results_all_targets
write.csv(neuron_results_all_targets,
          "/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/neuron_results_all_targets.csv")

#Now, neuron_results_targets_only has all the genes targeted in the screen, including the unexpected targets (for both promoters and enhancers)
#I'm missing 2 gRNA-target pairs
#The two missing are: chr4.2290_457_Gasperini_enhancer, and SPZ1_288_TSS_pos_ctrl - these probably didn't have enough coverage so were filtered out during the QC steps

neuron_results_all_targets
neuron_results_all_targets_promoters = neuron_results_all_targets
#then make dfs for just the hits
neuron_results_hits_promoters = neuron_results_all_targets_promoters[neuron_results_all_targets_promoters$hit == "yes",]
neuron_results_hits_promoters

#Now have df for promoter gRNAs in neurons - can use this for all those plots for the supp figure
#downsample the non-hits:
non_hits_neurons = neuron_results_all_targets_promoters[neuron_results_all_targets_promoters$hit == "no",]
hits_neurons = neuron_results_all_targets_promoters[neuron_results_all_targets_promoters$hit == "yes",]

non_hits_neurons_DS = non_hits_neurons[sample(nrow(non_hits_neurons), 17), ]
neuron_results_all_targets_promoters_DS = rbind(hits_neurons, non_hits_neurons_DS)
neuron_results_all_targets_promoters_DS

## B) Comparison of the percentage of cells expressing the target gene of gRNAs that resulted in an eFDR<0.1 versus gRNAs that resulted in an eFDR>0.1

#violin plots - neuron promoters
options(repr.plot.width=4, repr.plot.height=6)
ggplot(neuron_results_all_targets_promoters_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#just plots for figure without axis labels, etc.
options(repr.plot.width=4, repr.plot.height=6)
ggplot(neuron_results_all_targets_promoters_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  #theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_blank())
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/pct.2_neurons.png",
       width = 4,
       height = 6,
       dpi=300)

# D) Number of cells harboring each gRNA for gRNAs that resulted in an eFDR<0.1 versus gRNAs that resulted in an eFDR>0.1

#violin plots - neuron promoters
options(repr.plot.width=4, repr.plot.height=6)
ggplot(neuron_results_all_targets_promoters_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#just plots for figure without axis labels, etc.
options(repr.plot.width=4, repr.plot.height=6)
ggplot(neuron_results_all_targets_promoters_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  #theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_blank())
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/n_cells_neurons.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(neuron_results_all_targets_promoters_DS$n_cells, neuron_results_all_targets_promoters_DS$hit == "yes"),
            subset(neuron_results_all_targets_promoters_DS$n_cells, neuron_results_all_targets_promoters_DS$hit == "no"))

#wilcoxon rank sum test
wilcox.test(subset(neuron_results_all_targets_promoters_DS$pct.2, neuron_results_all_targets_promoters_DS$hit == "yes"),
            subset(neuron_results_all_targets_promoters_DS$pct.2, neuron_results_all_targets_promoters_DS$hit == "no"))

## Looking at normalized expression values
##Reading in data and creating a seaurat object

data_dir <- '/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/sequencing_data/cellranger_aggr_output/outs/count/filtered_feature_bc_matrix/'
list.files(data_dir)

CRISPRaQTL_Pilot_data <- Read10X(data.dir = data_dir)
CRISPRaQTL_Pilot_Seurat_Object = CreateSeuratObject(counts = CRISPRaQTL_Pilot_data$`Gene Expression`, min.cells = 4, min.features = 200)

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

##Normalize the data with seurats default paramters before statistical testing

CRISPRaQTL_Pilot_Seurat_Object_Subset <- NormalizeData(CRISPRaQTL_Pilot_Seurat_Object_Subset)

##Creating reference dataframes with gRNA-gene and gene-gRNA mappings

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/cellranger_input_files/feature_reference.csv")  %>% 
  rename(gRNA = id)

#need to get target genes from the promoter and enhancer dataframes
neuron_promoter_hit_list = unique(neuron_results_hits_promoters$target_guide)
length(neuron_promoter_hit_list)
neuron_promoter_hit_list

tail(FBC_DF)
#doing neuron promoters first
gRNA_FBC_DF <- FBC_DF  %>% 
  filter(Gene %in% neuron_promoter_hit_list) %>% 
  filter(UMI_Count > 5) %>% 
  select(Cell)
head(gRNA_FBC_DF)
gRNA_Cells <- gRNA_FBC_DF$Cell 

##Isolate control cells

# This makes control cells all cells that dont have the targeting guide
Control_Cells_DF <- FBC_DF %>% 
  filter(!Cell %in% gRNA_Cells)

Control_Cells <- unique(Control_Cells_DF$Cell)

##Assign those cells to identities in Seurat object for plotting

Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = Control_Cells) <- "Control Cells"
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "promoter gRNA Cells"
length(gRNA_Cells)
length(Control_Cells)

#Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset)
levels(x = CRISPRaQTL_Pilot_Seurat_Object_Subset)

average_exprsn = as.data.frame(AverageExpression(
  CRISPRaQTL_Pilot_Seurat_Object_Subset,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
))

dim(average_exprsn)
head(average_exprsn)

#promoters - neuron
ft_ref_promoters = Ft_Ref %>%
  filter(!grepl("Gasperini|NTC", gRNA))
#ft_ref_promoters

target_gene_list_promoters = unique(ft_ref_promoters$target_gene_name)

target_gene_list_promoters = append(target_gene_list_promoters, c('CALM3'))
target_gene_list_promoters = target_gene_list_promoters[target_gene_list_promoters != "PPP5D1"]

target_gene_list_promoters

average_exprsn = subset(average_exprsn, rownames(average_exprsn) %in% target_gene_list_promoters)
#average_exprsn
#rownames(average_exprsn)

average_exprsn_hits = subset(average_exprsn, rownames(average_exprsn) %in% neuron_results_hits_promoters$target_gene)
#average_exprsn_hits

#now can make violin plots comparing the control cells in the hits vs. non hit dataframes! 
average_exprsn_non_hits = subset(average_exprsn, !rownames(average_exprsn) %in% neuron_results_hits_promoters$target_gene)
#average_exprsn_non_hits

#now, make a new df including just control cells
average_exprsn_hits$hit = "yes"
hit_control_cell_exprsn = average_exprsn_hits[,c(2,3)]
dim(hit_control_cell_exprsn)
#hit_control_cell_exprsn

average_exprsn_non_hits$hit = "no"
non_hit_control_cell_exprsn = average_exprsn_non_hits[,c(2,3)]
dim(non_hit_control_cell_exprsn)
#non_hit_control_cell_exprsn

all_data = rbind(hit_control_cell_exprsn, non_hit_control_cell_exprsn)
head(all_data)
tail(all_data)

#now can make a violin plot
#violin plots - log10 scale K562
options(repr.plot.width=6, repr.plot.height=6)
ggplot(all_data, aes(x = hit, y = RNA.Control.Cells, fill = hit)) + 
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_y_continuous(trans='log10') +  
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

options(repr.plot.width=6, repr.plot.height=6)
ggplot(all_data, aes(x = hit, y = RNA.Control.Cells	, fill = hit)) + 
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_y_continuous(trans='log10') +  
  scale_x_discrete(labels = c('non hit gRNAs', 'hit gRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  #theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_blank())
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/norm_exprsn_neurons.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(all_data$RNA.Control.Cells, all_data$hit == "yes"),
            subset(all_data$RNA.Control.Cells, all_data$hit == "no"))

#### Looking at GC content

#neuron promoters
ft_ref_short = Ft_Ref[,c(1,5)]
head(ft_ref_short)

#neuron_hit_list = neuron_hits$target_guide
#neuron_hit_list

gRNA_info_hits = subset(ft_ref_short, ft_ref_short$gRNA %in% neuron_results_hits_promoters$target_guide)
gRNA_info_hits

gRNA_info_non_hits = subset(ft_ref_short, !ft_ref_short$gRNA %in% neuron_results_hits_promoters$target_guide)
#gRNA_info_non_hits = subset(ft_ref_short, !ft_ref_short$gRNA %in% neuron_enhancer_hits)

#gRNA_info_non_hits
#take out NTCs 

gRNA_info_non_hits = gRNA_info_non_hits %>%
  filter(!grepl('NTC', gRNA))
gRNA_info_non_hits

#downsample to 48
gRNA_info_non_hits = gRNA_info_non_hits[sample(nrow(gRNA_info_non_hits), 17), ]
gRNA_info_non_hits

#Calculate GC content of hit vs. non-hit gRNAs
#Add column with GC content (%)
gRNA_info_hits$GC_content = ((str_count(gRNA_info_hits$sequence, "C") + str_count(gRNA_info_hits$sequence, "G"))/19)*100
mean(gRNA_info_hits$GC_content)
#gRNA_info_hits

gRNA_info_non_hits$GC_content = ((str_count(gRNA_info_non_hits$sequence, "C") + str_count(gRNA_info_non_hits$sequence, "G"))/19)*100
#gRNA_info_non_hits
mean(gRNA_info_non_hits$GC_content)

gRNA_info_hits$hit = "yes"
gRNA_info_non_hits$hit = "no"

gRNA_info_rbound = rbind(gRNA_info_hits, gRNA_info_non_hits)
#gRNA_info_rbound

options(repr.plot.width=6, repr.plot.height=6)
ggplot(gRNA_info_rbound, aes(x = hit, y = GC_content, fill = hit)) + 
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "GC content of gRNA (%)") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  scale_x_discrete(labels = c('non hit \n gRNAs', 'hit \ngRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 18)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text.x = element_text(family="Arial", colour = "black", size = 18))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#without axis and tick titles
ggplot(gRNA_info_rbound, aes(x = hit, y = GC_content, fill = hit)) + 
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  geom_violin(size = 0.5, trim = FALSE) +
  theme_bw()+
  #scale_y_continuous(trans='log10')+
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm"), legend.position = "none") +
  labs(title = "", x = "", y = "") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  #scale_y_continuous(trans='log10', labels = scales::number_format(accuracy = 1.0), limits = c(0.0001,10), breaks = c(0.1,1,10)) +
  #scale_x_discrete(labels = c('non hit \n gRNAs', 'hit \ngRNAs'))+
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 24)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_blank())
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/GC_content_neuron.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(gRNA_info_rbound$GC_content, gRNA_info_rbound$hit == "yes"),
            subset(gRNA_info_rbound$GC_content, gRNA_info_rbound$hit == "no"))


















