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
# 1) K562 promoter hits (list)
# 2) K562 enhancer hits (list)
# 3) K562 aggregated results 1 Mb (data frame)
# 4) K562 aggregated results, promoter hits only (data frame)
# 6) K562 aggregated results, enhancer hits only (data frame)

#### 1. Do hit gRNAs biasly target highly or lowly expressed genes?
#all non-hit gRNAs contain the data for genes that were targeted in the screen, but weren't hits

#filter all DE testing results for all the target genes in the K562 hits df;
#that way, should have all gRNAs against all intended and unintended targets (just a few of these):
Ft_Ref = read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
  rename(target_guide = id) %>%
  rename(target_gene = target_gene_name)

K562_results = read.table("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/aggregated_results.txt",
                          sep = "\t",
                          header = TRUE)
#tail(K562_results)
#K562_hits
K562_hits = read.table("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv",
                       sep = ",",
                       header = TRUE)
K562_hit_list = K562_hits$target_guide
length(K562_hit_list)
#head(K562_hits)

#Filter aggregated results by target genes in K562_hits and get rid of NTCs (don't want for this analysis)
#Exceptions for K562: chr6.2134_473_Gasperini_NTC (enhancer) and the two ASIC1 enhancers, CALM3, RPS18 and WWC3 (promoters)
#K562_results_targets_only = K562_results %>%
#filter(target_gene %in% K562_hits$target_gene | target_gene %in% Ft_Ref$target_gene_name) #this is so I can get all targets
#K562_results_targets_only = K562_results %>%
#    filter(target_gene %in% Ft_Ref$target_gene | target_guide %in% Ft_Ref$target_guide) #this is so I can get all the intended targets, will add back unintended targets later

#merge keeps all the matching gRNA/target_gene pairs with the DE results data
K562_results_targets_only = merge(K562_results, Ft_Ref) #NTCs get filtered out during this step
#K562_results_targets_only

#filter out the exceptions: CALM3, RPS18, and WWC3 promoter hits
CALM3_hit = K562_results %>%
  filter(grepl('PPP5D1_238_TSS_pos_ctrl', target_guide), grepl("CALM3", target_gene))
RPS18_hit = K562_results %>%
  filter(grepl('SYNGAP1_295_Flashfry_promoter', target_guide), grepl("RPS18", target_gene))
WWC3_hit = K562_results %>%
  filter(grepl('FOXP1_146_Flashfry_promoter', target_guide), grepl("WWC3", target_gene))
promoter_exceptions = rbind(CALM3_hit,
                            RPS18_hit,
                            WWC3_hit)
#promoter_exceptions

#take out the three non-specific NTC enhancer hits before I filter out NTCs:
chr6.2134_473_Gasperini_NTC_hit = K562_results %>%
  filter(grepl('chr6.2134_473_Gasperini_NTC', target_guide), grepl("HMGA1", target_gene))
chr12.1559_411_Gasperini_enhancer_hit = K562_results %>%
  filter(grepl('chr12.1559_411_Gasperini_enhancer', target_guide), grepl("ASIC1", target_gene))
chr12.1559_412_Gasperini_enhancer_hit = K562_results %>%
  filter(grepl('chr12.1559_412_Gasperini_enhancer', target_guide), grepl("ASIC1", target_gene))

enhancer_exceptions = rbind(chr6.2134_473_Gasperini_NTC_hit,
                            chr12.1559_411_Gasperini_enhancer_hit,
                            chr12.1559_412_Gasperini_enhancer_hit)
#enhancer_exceptions
K562_results_targets_only = K562_results_targets_only[,c(3,4,5,6,7,1,2,8,9)]
#head(K562_results_targets_only)
#head(promoter_exceptions)
#head(enhancer_exceptions)
#adding the promoter and enhancer exceptions back
K562_results_all_targets = rbind(K562_results_targets_only, promoter_exceptions, enhancer_exceptions)
K562_results_all_targets

#Now want to add in if each row is a hit or not - can get this information from K562_hits
K562_results_all_targets$hit[K562_results_all_targets$target_guide %in% K562_hit_list] = "yes"
K562_results_all_targets$hit[!K562_results_all_targets$target_guide %in% K562_hit_list] = "no"
table(K562_results_all_targets$hit)

#Need to manually fix the exceptions (chr12.1559_411_Gasperini_enhancer, chr12.1559_412_Gasperini_enhancer,
#PPP5D1_238_TSS_pos_ctrl, SYNGAP1_295_Flashfry_promoter, FOXP1_146_Flashfry_promoter)
K562_results_all_targets$hit[K562_results_all_targets$target_guide == "PPP5D1_238_TSS_pos_ctrl" &
                               K562_results_all_targets$target_gene == "PPP5D1" |
                               K562_results_all_targets$target_guide == "SYNGAP1_295_Flashfry_promoter" &
                               K562_results_all_targets$target_gene == "SYNGAP1"|
                               K562_results_all_targets$target_guide == "FOXP1_146_Flashfry_promoter" &
                               K562_results_all_targets$target_gene == "FOXP1"|
                               K562_results_all_targets$target_guide == "chr12.1559_411_Gasperini_enhancer" &
                               K562_results_all_targets$target_gene == "TUBA1A"|
                               K562_results_all_targets$target_guide == "chr12.1559_412_Gasperini_enhancer" &
                               K562_results_all_targets$target_gene == "TUBA1A"] = "no"

table(K562_results_all_targets$hit) 

#K562_results_all_targets
write.csv(K562_results_all_targets,
          "/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/K562_results_all_targets.csv")

#Now, K562_results_targets_only has all the genes targeted in the screen, including the unexpected targets (for both promoters and enhancers)
#I'm missing 2 gRNA-target pairs
#The two missing are: chr4.2290_457_Gasperini_enhancer, and SPZ1_288_TSS_pos_ctrl - these must have been filtered out during the QC steps

#now, split the above df into promoters and enhancers - can do this by the string "Gasperini"
K562_results_all_targets_promoters = K562_results_all_targets %>%
  filter(!grepl('*Gasperini*', target_guide)) 
K562_results_all_targets_promoters
length(unique(K562_results_all_targets_promoters$target_guide))

K562_results_all_targets_enhancers = K562_results_all_targets %>%
  filter(grepl('*Gasperini*', target_guide)) 

K562_results_all_targets_enhancers
length(unique(K562_results_all_targets_enhancers$target_guide))

#then make dfs for just the hits
K562_results_hits_promoters = K562_results_all_targets_promoters[K562_results_all_targets_promoters$hit == "yes",]
K562_results_hits_enhancers = K562_results_all_targets_enhancers[K562_results_all_targets_enhancers$hit == "yes",]
K562_results_hits_promoters
K562_results_hits_enhancers

#Now have df for promoter and enhancers in K562 cells - can use this for all those plots for this supp figure
#Now have df for promoter gRNAs in K562 - can use this for all those plots for the supp figure
#downsample the non-hits:
non_hits_K562 = K562_results_all_targets_promoters[K562_results_all_targets_promoters$hit == "no",]
hits_K562 = K562_results_all_targets_promoters[K562_results_all_targets_promoters$hit == "yes",]

non_hits_K562_DS = non_hits_K562[sample(nrow(non_hits_K562), 48), ]
K562_results_all_targets_promoters_DS = rbind(hits_K562, non_hits_K562_DS)
K562_results_all_targets_promoters_DS

non_hits_K562 = K562_results_all_targets_enhancers[K562_results_all_targets_enhancers$hit == "no",]
hits_K562 = K562_results_all_targets_enhancers[K562_results_all_targets_enhancers$hit == "yes",]

non_hits_K562_DS = non_hits_K562[sample(nrow(non_hits_K562), 12), ]
K562_results_all_targets_enhancers_DS = rbind(hits_K562, non_hits_K562_DS)
K562_results_all_targets_enhancers_DS

#violin plots - K562 promoters
options(repr.plot.width=4, repr.plot.height=6)
ggplot(K562_results_all_targets_promoters_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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

#violin plots - K562 enhancers
options(repr.plot.width=4, repr.plot.height=6)
ggplot(K562_results_all_targets_enhancers_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggplot(K562_results_all_targets_promoters_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/pct.2_K562_promoters.png",
       width = 4,
       height = 6,
       dpi=300)

ggplot(K562_results_all_targets_enhancers_DS, aes(x = hit, y = pct.2, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/pct.2_K562_enhancers.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test to test for significance
wilcox.test(subset(K562_results_all_targets_promoters_DS$pct.2, K562_results_all_targets_promoters_DS$hit == "yes"),
            subset(K562_results_all_targets_promoters_DS$pct.2, K562_results_all_targets_promoters_DS$hit == "no"))

wilcox.test(subset(K562_results_all_targets_enhancers_DS$pct.2, K562_results_all_targets_enhancers_DS$hit == "yes"),
            subset(K562_results_all_targets_enhancers_DS$pct.2, K562_results_all_targets_enhancers_DS$hit == "no"))

#violin plots - K562 promoters
options(repr.plot.width=4, repr.plot.height=6)
ggplot(K562_results_all_targets_promoters_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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

#violin plots - K562 enhancers
options(repr.plot.width=4, repr.plot.height=6)
ggplot(K562_results_all_targets_enhancers_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggplot(K562_results_all_targets_promoters_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/n_cells_promoter_hits_K562.png",
       width = 4,
       height = 6,
       dpi=300)

ggplot(K562_results_all_targets_enhancers_DS, aes(x = hit, y = n_cells, fill = hit)) + scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/n_cells_enhancer_hits_K562.png",
       width = 4,
       height = 6,
       dpi=300)

wilcox.test(subset(K562_results_all_targets_promoters_DS$n_cells, K562_results_all_targets_enhancers_DS$hit == "yes"),
            subset(K562_results_all_targets_promoters_DS$n_cells, K562_results_all_targets_enhancers_DS$hit == "no"))

wilcox.test(subset(K562_results_all_targets_enhancers_DS$n_cells, K562_results_all_targets_enhancers_DS$hit == "yes"),
            subset(K562_results_all_targets_enhancers_DS$n_cells, K562_results_all_targets_enhancers_DS$hit == "no"))

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

#need to get target genes from the promoter and enhancer dataframes
K562_promoter_hit_list = unique(K562_results_hits_promoters$target_guide)
length(K562_promoter_hit_list)
K562_promoter_hit_list

tail(FBC_DF)
#doing K562 promoters first
gRNA_FBC_DF <- FBC_DF  %>% 
  filter(Gene %in% K562_promoter_hit_list) %>% 
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

#promoters - K562
ft_ref_promoters = Ft_Ref %>%
  filter(!grepl("Gasperini|NTC", gRNA))
#ft_ref_promoters

ft_ref_enhancers = Ft_Ref %>%
  filter(grepl("Gasperini_enhancer|chr6.2134_473_Gasperini_NTC", gRNA))
#ft_ref_enhancers

target_gene_list_promoters = unique(ft_ref_promoters$target_gene_name)
target_gene_list_enhancers = unique(ft_ref_enhancers$target_gene_name)

target_gene_list_promoters = append(target_gene_list_promoters, c('WWC3', 'CALM3', 'RPS18'))
target_gene_list_promoters = target_gene_list_promoters[target_gene_list_promoters != "PPP5D1"]

target_gene_list_enhancers = append(target_gene_list_enhancers, "HMGA1")
target_gene_list_enhancers = target_gene_list_enhancers[target_gene_list_enhancers != "TUBA1A"]
target_gene_list_enhancers = target_gene_list_enhancers[target_gene_list_enhancers != "Non-Targeting"]

target_gene_list_promoters
target_gene_list_enhancers

average_exprsn = subset(average_exprsn, rownames(average_exprsn) %in% target_gene_list_promoters)
#average_exprsn
#rownames(average_exprsn)

average_exprsn_hits = subset(average_exprsn, rownames(average_exprsn) %in% K562_results_hits_promoters$target_gene)
#average_exprsn_hits

#now can make violin plots comparing the control cells in the hits vs. non hit dataframes! 
average_exprsn_non_hits = subset(average_exprsn, !rownames(average_exprsn) %in% K562_results_hits_promoters$target_gene)
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

#violin plots 
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/norm_exprsn_K562_promoters.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(all_data$RNA.Control.Cells, all_data$hit == "yes"),
            subset(all_data$RNA.Control.Cells, all_data$hit == "no"))

#now, do the enhancers
#need to get target genes from the promoter and enhancer dataframes
K562_enhancer_hit_list = unique(K562_results_hits_enhancers$target_guide)
length(K562_enhancer_hit_list)
K562_enhancer_hit_list

tail(FBC_DF)
#doing K562 promoters first
gRNA_FBC_DF <- FBC_DF  %>% 
  filter(Gene %in% K562_enhancer_hit_list) %>% 
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
Idents(object = CRISPRaQTL_Pilot_Seurat_Object_Subset, cells = gRNA_Cells) <- "enhancer gRNA Cells"
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

average_exprsn = subset(average_exprsn, rownames(average_exprsn) %in% target_gene_list_enhancers)
#average_exprsn
#rownames(average_exprsn)

average_exprsn_hits = subset(average_exprsn, rownames(average_exprsn) %in% K562_results_hits_enhancers$target_gene)
#average_exprsn_hits

#now can make violin plots comparing the control cells in the hits vs. non hit dataframes! 
average_exprsn_non_hits = subset(average_exprsn, !rownames(average_exprsn) %in% K562_results_hits_enhancers$target_gene)
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

#not downsampling bc the data sizes are small (7 and 18 rows)
#non_hit_control_cell_exprsn_DS = non_hit_control_cell_exprsn[sample(nrow(non_hit_control_cell_exprsn), 7), ]
#dim(non_hit_control_cell_exprsn_DS)

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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/norm_exprsn_K562_enhancers.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(all_data$RNA.Control.Cells, all_data$hit == "yes"),
            subset(all_data$RNA.Control.Cells, all_data$hit == "no"))


## Looking at GC content
#K562 promoters
ft_ref_short = Ft_Ref[,c(1,5)]
head(ft_ref_short)

#K562_hit_list = K562_hits$target_guide
#K562_hit_list

gRNA_info_hits = subset(ft_ref_short, ft_ref_short$gRNA %in% K562_results_hits_promoters$target_guide)
gRNA_info_hits

gRNA_info_non_hits = subset(ft_ref_short, !ft_ref_short$gRNA %in% K562_results_hits_promoters$target_guide)
#gRNA_info_non_hits = subset(ft_ref_short, !ft_ref_short$gRNA %in% K562_enhancer_hits)

#gRNA_info_non_hits
#take out NTCs 

gRNA_info_non_hits = gRNA_info_non_hits %>%
  filter(!grepl('NTC', gRNA))
gRNA_info_non_hits

#downsample to 48
gRNA_info_non_hits = gRNA_info_non_hits[sample(nrow(gRNA_info_non_hits), 48), ]
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
ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/gRNA_analysis_Feb_2023/GC_content_K562_promoters.png",
       width = 4,
       height = 6,
       dpi=300)

#wilcoxon rank sum test
wilcox.test(subset(gRNA_info_rbound$GC_content, gRNA_info_rbound$hit == "yes"),
            subset(gRNA_info_rbound$GC_content, gRNA_info_rbound$hit == "no"))





