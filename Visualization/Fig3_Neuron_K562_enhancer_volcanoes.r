##Loading libraries

library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)
library(corrplot)
library(RColorBrewer)
library(scales)

##K562 volcano 

options(repr.plot.width=7, repr.plot.height=6)

##Read in DE results

Full_Results <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/aggregated_results.txt")%>%
  mutate(plot_pval = p_val+2.225074e-308) %>% 
  mutate(plot_pval = -log10(plot_pval))  %>% 
  unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)

Full_Results %>% 
    arrange(p_val) 


##Read in detected neighbouring genes DF with additional metadata

Detc_Neighbouring_Genes <- readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 

##Left joining

Full_Results <- Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")

##Splitting out targeting and NTC results

Targeting_Results <- Full_Results %>% 
    filter(Targeting == "Targeting")

Targeting_Results %>% 
arrange(p_val)

##NTC results

set.seed(12)

NTC_Results <- Full_Results %>% 
    filter(Targeting == "Non-targeting") 

NTC_Results %>% 
arrange(p_val)

##Isolating just primary targeting results

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
 
Primary_Target_Tests <- Ft_Ref$Test_ID

Primary_Targeting_Results <- Targeting_Results %>% 
 filter(Test_ID %in% Primary_Target_Tests) 

##Isolating just enhancer results

En_Results_Primary_Targets <- Primary_Targeting_Results %>% 
    filter(grepl("chr", Test_ID))

##Downsampling the NTCs then recombining for plotting

##NTC results

set.seed(11)

En_NTC_Results <- NTC_Results %>% 
 filter(target_gene %in% En_Results_Primary_Targets$target_gene)

En_NTC_Results <- En_NTC_Results %>%   
    slice_sample(n = length(En_Results_Primary_Targets$p_val)) 

En_NTC_Results

##Combining for plotting

En_DS_Full_Results <- rbind(En_Results_Primary_Targets, En_NTC_Results)


##Volcano time

ggplot(En_DS_Full_Results, aes(x = avg_log2FC, y = plot_pval, colour = Targeting)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(legend.position = "none") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank()) +
  ylim(0,80) +
  xlim(-0.5,1.5)

##Neuron volcano 

options(repr.plot.width=7, repr.plot.height=6)

##Read in DE results

Full_Results <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/aggregated_results.txt")%>%
  mutate(plot_pval = p_val+2.225074e-308) %>% 
  mutate(plot_pval = -log10(plot_pval))  %>% 
  unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)

Full_Results %>% 
    arrange(p_val) 


##Read in detected neighbouring genes DF with additional metadata

Detc_Neighbouring_Genes <- readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 

##Left joining

Full_Results <- Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")

##Splitting out targeting and NTC results

Targeting_Results <- Full_Results %>% 
    filter(Targeting == "Targeting")

Targeting_Results %>% 
arrange(p_val)

##NTC results

set.seed(12)

NTC_Results <- Full_Results %>% 
    filter(Targeting == "Non-targeting") 

NTC_Results %>% 
arrange(p_val)

##Isolating just primary targeting results

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
 
Primary_Target_Tests <- Ft_Ref$Test_ID

Primary_Targeting_Results <- Targeting_Results %>% 
 filter(Test_ID %in% Primary_Target_Tests) 

##Isolating just enhancer results

En_Results_Primary_Targets <- Primary_Targeting_Results %>% 
    filter(grepl("chr", Test_ID))

##Downsampling the NTCs then recombining for plotting

##NTC results

set.seed(11)

En_NTC_Results <- NTC_Results %>% 
 filter(target_gene %in% En_Results_Primary_Targets$target_gene)

En_NTC_Results <- En_NTC_Results %>%   
    slice_sample(n = length(En_Results_Primary_Targets$p_val)) 

En_NTC_Results

##Combining for plotting

En_DS_Full_Results <- rbind(En_Results_Primary_Targets, En_NTC_Results)


##Volcano time

ggplot(En_DS_Full_Results, aes(x = avg_log2FC, y = plot_pval, colour = Targeting)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(legend.position = "none") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank()) +
  ylim(0,80) +
  xlim(-0.5,1.5)
