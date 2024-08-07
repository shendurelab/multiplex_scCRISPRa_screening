#### Figure S9 | Comparison of K562 CRISPRa-QTL hits vs. neuron CRISPRa-QTL hits. 

# A) Venn diagram showing number of overlapping promoter-targeting gRNA hits (left) and enhancer-targeting gRNA hits (right) between the K562 and neuron CRISPRa screens. 
# B) Correlation plots of log2 fold changes of TSS positive control targeting gRNAs (left) and enhancer targeting gRNAs (right) between the K562 and neuron CRISPRa screens. 
# C) Strong CRISPRaQTL scores for targeting gRNAs are enriched for proximity to their target gene, NTCs are not. 
# D) Same plot as in Fig. S7C, with the y-axis zoomed in from 0 to 50. 

library(tidyverse)
library(repr)
library(Matrix)
library(Seurat)
library(broom)
library(ggridges)
library(ggrepel)
library(patchwork)
library(data.table)
library(VennDiagram)
library(RColorBrewer)

#### A) Venn diagram of TSS+promoter (top) and enhancer (bottom) hit overlap between cell types
## Input data tables to use for now, but substitute later with the final tables (once Troy has them)

K562_results = read.table("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Updated_input_files_011823/K562_Primary_Targeting_Results.csv",
                          sep = ",",
                          header = TRUE)

neuron_results = read.table("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Updated_input_files_011823/Neuron_Primary_Targeting_Results.csv",
                            sep = ",",
                            header = TRUE)

K562_hits_table = read.table("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Updated_input_files_011823/K562_Hits.csv",
                             sep = ",",
                             header = TRUE)
neuron_hits_table = read.table("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Updated_input_files_011823/Neuron_Hits.csv",
                               sep = ",",
                               header = TRUE)

K562_hits = K562_hits_table$target_guide
neuron_hits = neuron_hits_table$target_guide

length(K562_hits)
length(neuron_hits)

K562_promoter_hits = K562_hits[grepl("*promoter*|*TSS*", K562_hits)]
K562_promoter_hits
length(K562_promoter_hits)

neuron_promoter_hits = neuron_hits[grepl("*promoter*|*TSS*", neuron_hits)]
neuron_promoter_hits
length(neuron_promoter_hits)

setwd("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/venn_diagrams_v2")
getwd()

venn.diagram(
  x = list(K562_promoter_hits, neuron_promoter_hits),
  category.names = c("" , ""),
  filename = 'promoter_venn_diagramm.png',
  output=TRUE,
  lwd = 1,
  cex = .6,
  fontface = "bold",
  fontfamily = "Arial",
  fill = c(alpha("#56B4E9",0.9), alpha('steelblue',0.9)),
  col=c("#56B4E9","steelblue"))

# with no numbers
venn.diagram(
  x = list(K562_promoter_hits, neuron_promoter_hits),
  category.names = c("" , ""),
  filename = 'promoter_venn_diagramm_no_numbers.png',
  output=TRUE,
  lwd = 1,
  cex = 0,
  fontface = "bold",
  fontfamily = "Arial",
  fill = c(alpha("#56B4E9",0.9), alpha('steelblue',0.9)),
  col=c("#56B4E9","steelblue"))

## Enhancer hit overlaps; chr6.2134_473_Gasperini_NTC was an enhancer-gene pair hit in our CRISPRa screen

K562_enhancer_hits = K562_hits[grepl("*enhancer|chr6.2134_473_Gasperini_NTC*", K562_hits)]
K562_enhancer_hits
length(K562_enhancer_hits)

neuron_enhancer_hits = neuron_hits[grepl("*enhancer*", neuron_hits)]
neuron_enhancer_hits
length(neuron_enhancer_hits)

venn.diagram(
  x = list(K562_enhancer_hits, neuron_enhancer_hits),
  category.names = c("" , ""),
  filename = 'enhancer_venn_diagramm.png',
  output=TRUE,
  lwd = 1,
  cex = .6,
  fontface = "bold",
  fontfamily = "Arial",
  fill = c(alpha("#56B4E9",0.9), alpha('steelblue',0.9)),
  col=c("#56B4E9","steelblue"))

# with no numbers
venn.diagram(
  x = list(K562_enhancer_hits, K562_enhancer_hits),
  category.names = c("" , ""),
  filename = 'enhancer_venn_diagramm_no_numbers.png',
  output=TRUE,
  lwd = 1,
  cex = 0,
  fontface = "bold",
  fontfamily = "Arial",
  fill = c(alpha("#56B4E9",0.9), alpha('steelblue',0.9)),
  col=c("#56B4E9","steelblue"))

#### B) Correlation in log2FC between neurons and K562 for TSS controls

K562_TSS_controls = K562_results[grepl("*TSS*", K562_results$target_guide),]
nrow(K562_TSS_controls)

neuron_TSS_controls = neuron_results[grepl("*TSS*", neuron_results$target_guide),]
nrow(neuron_TSS_controls)

## Make a merged df
K562_TSS_controls_log2fc = K562_TSS_controls[,c(2,7)]
colnames(K562_TSS_controls_log2fc) = c("K562_avg_log2FC", "target_guide")

neuron_TSS_controls_log2fc = neuron_TSS_controls[,c(2,7)]
colnames(neuron_TSS_controls_log2fc) = c("neuron_avg_log2FC", "target_guide")

K562_and_neuron_log2fc = merge(K562_TSS_controls_log2fc, neuron_TSS_controls_log2fc)
nrow(K562_and_neuron_log2fc)

## Corelation plot for log2FC

lm_eqn <- function(df){
  m <- lm(neuron_avg_log2FC ~ K562_avg_log2FC, df);
  eq <- substitute(~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}

ggplot(K562_and_neuron_log2fc, aes(x = K562_avg_log2FC, y = neuron_avg_log2FC, color = "steelblue4", label = target_guide)) +
  geom_hline(yintercept=0, color = "darkgrey", size = 1)+
  geom_vline(xintercept=0, color = "darkgrey", size = 1)+
  geom_point(size = 3) +
  scale_color_manual(values = c("steelblue4")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "K562 gRNAs (log2FC)", y = "Neuron gRNAs (log2FC)") +
  theme(text = element_text(family="Arial", colour = "black", size = 20)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position="none")+
  geom_smooth(method = "lm", se=FALSE, size = 0)+
  xlim(-0.3,1.5)+ylim(-0.3,1.5)+ 
  geom_text(x = 0.1, y = 1.25, label = lm_eqn(K562_and_neuron_log2fc), parse = TRUE, size = 8)+
  geom_smooth(method = "lm", se = FALSE)


ggplot(K562_and_neuron_log2fc, aes(x = K562_avg_log2FC, y = neuron_avg_log2FC, color = "steelblue4", label = target_guide)) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", size = 1)+
  geom_vline(xintercept=0, linetype="dashed", color = "darkgrey", size = 1)+
  geom_point(size = 3) +
  scale_color_manual(values = c("steelblue4")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "K562 gRNAs (log2FC)", y = "Neuron gRNAs (log2FC)") +
  theme(text = element_text(family="Arial", colour = "black", size = 20)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position="none")+
  geom_smooth(method = "lm", se=FALSE, size = 0)+
  xlim(-0.3,1.5)+ylim(-0.3,1.5)+ 
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank())+
  coord_fixed()+
  theme(legend.position="none")+
  geom_smooth(method = "lm", se = FALSE)
#ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Fig_S7/log2FC_correlation_plot_TSS_controls.png", 
#dpi=300)

#### C) Correlation in log2FC between neurons and K562 for enhancers

K562_enhancers = K562_results[grepl("*enhancer|chr6.2134_473_Gasperini_NTC*", K562_results$target_guide),]
nrow(K562_enhancers)

neuron_enhancers = neuron_results[grepl("*enhancer*", neuron_results$target_guide),]
nrow(neuron_enhancers)

## Make a merged df
K562_enhancers_log2fc = K562_enhancers[,c(2,7)]
colnames(K562_enhancers_log2fc) = c("K562_avg_log2FC", "target_guide")

neuron_enhancers_log2fc = neuron_enhancers[,c(2,7)]
colnames(neuron_enhancers_log2fc) = c("neuron_avg_log2FC", "target_guide")

K562_and_neuron_log2fc_enhancers = merge(K562_enhancers_log2fc, neuron_enhancers_log2fc)
nrow(K562_and_neuron_log2fc_enhancers)
K562_and_neuron_log2fc_enhancers

## Corelation plot for enhancers (log2FC)

lm_eqn <- function(df){
  m <- lm(neuron_avg_log2FC ~ K562_avg_log2FC, df);
  eq <- substitute(~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(-sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}


ggplot(K562_and_neuron_log2fc_enhancers, aes(x = K562_avg_log2FC, y = neuron_avg_log2FC, color = "#56B4E9", label = target_guide)) +
  geom_hline(yintercept=0, color = "darkgrey", size = 1)+
  geom_vline(xintercept=0, color = "darkgrey", size = 1)+  
  geom_point(size = 3) +
  scale_color_manual(values = c("#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "K562 gRNAs (log2FC)", y = "Neuron gRNAs (log2FC)") +
  theme(text = element_text(family="Arial", colour = "black", size = 20)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position="none")+
  geom_smooth(method = "lm", se=FALSE, size = 0)+
  xlim(-0.2,1.1)+ylim(-0.2,1.1)+ 
  geom_text(x = 0, y = 0.8, label = lm_eqn(K562_and_neuron_log2fc_enhancers), parse = TRUE, size = 8)+
  geom_smooth(method = "lm", se = FALSE)


ggplot(K562_and_neuron_log2fc_enhancers, aes(x = K562_avg_log2FC, y = neuron_avg_log2FC, color = "#56B4E9", label = target_guide)) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", size = 1)+
  geom_vline(xintercept=0, linetype="dashed", color = "darkgrey", size = 1)+  
  geom_point(size = 3) +
  scale_color_manual(values = c("#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "K562 gRNAs (log2FC)", y = "Neuron gRNAs (log2FC)") +
  theme(text = element_text(family="Arial", colour = "black", size = 20)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position="none")+
  geom_smooth(method = "lm", se=FALSE, size = 0)+
  xlim(-0.2,1.1)+ylim(-0.2,1.1)+ 
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank())+
  coord_fixed()+
  theme(legend.position="none")+
  geom_smooth(method = "lm", se = FALSE)
#ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Fig_S7/log2FC_correlation_plot_enhancers.png", 
#dpi=300)

#### C) Strong CRISPRaQTL scores for targeting gRNAs are enriched for proximity to their target gene, NTCs are not. D) Same plot as in Fig. S7C, with the y-axis zoomed in from 0 to 50. 

Full_Results = read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/aggregated_results.txt")%>%
  mutate(plot_pval = p_val+2.225074e-308) %>% 
  mutate(plot_pval = -log10(plot_pval))  %>% 
  filter(pct.1 > 0.0019) %>%
  filter(pct.1 > 0.0019) %>%
  unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)

#Read in detected neighbouring genes DF with additional metadata
Detc_Neighbouring_Genes = readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
  unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 
Full_Results = Full_Results %>% 
  left_join(Detc_Neighbouring_Genes, by = "Test_ID")

#Splitting out targeting and NTC results
Targeting_Results = Full_Results %>% 
  filter(Targeting == "Targeting")
set.seed(11)
NTC_Results = Full_Results %>% 
  filter(Targeting == "Non-targeting") 
NTC_Results <- NTC_Results %>%   
  slice_sample(n = length(Targeting_Results$p_val)) 
DS_Full_Results <- rbind(NTC_Results, Targeting_Results) 

##Target gene distance vs. p values for upstream gRNAs

DS_Full_Results_Rev = rbind(Targeting_Results, NTC_Results)
DS_Full_Results_Rev_US = DS_Full_Results_Rev %>% 
  filter(target_gene_distance < 0)

##Target gene distance vs. p-values for upstream guides
options(repr.plot.width=14, repr.plot.height=4)
ggplot(DS_Full_Results_Rev_US, aes(x = target_gene_distance, y = plot_pval, colour = Targeting)) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  geom_point(size = 1.5) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "Target gene distance (bp)", y = "-log10(observed P-values)") +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text = element_text(colour = "black"))

##Same plot as above without axis labels for final figure 
options(repr.plot.width=58, repr.plot.height=6)
DS_Full_Results_Rev_US <- DS_Full_Results_Rev %>% 
  filter(target_gene_distance < 0)
options(repr.plot.width=58, repr.plot.height=6)
ggplot(DS_Full_Results_Rev_US, aes(x = target_gene_distance, y = plot_pval, colour = Targeting)) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  geom_point(size = 6.5) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1.4)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.4)) +
  theme(axis.ticks.length=unit(.5, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "none")
#ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Fig_S7/neuron_target_gene_distance.png", 
#height = 4,
#width = 18,
#dpi=300)

##Same plot as above without axis labels and y-axis clipped at 50
options(repr.plot.width=58, repr.plot.height=6)
DS_Full_Results_Rev_US <- DS_Full_Results_Rev %>% 
  filter(target_gene_distance < 0)
ggplot(DS_Full_Results_Rev_US, aes(x = target_gene_distance, y = plot_pval, colour = Targeting)) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  geom_point(size = 6.5) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1.4)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.4)) +
  theme(axis.ticks.length=unit(.5, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "none") +
  ylim(0,50)
#ggsave("/net/shendure/vol1/home/chardonf/ResQTL/Paper_Figures/Fig_S7/neuron_target_gene_distance_y_axis_to_50.png", 
#height = 4,
#width = 18,
#dpi=300)
