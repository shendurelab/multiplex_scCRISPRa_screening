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

##Read in DE results

Full_Results <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/aggregated_results.txt")%>%
  mutate(plot_pval = p_val+2.225074e-308) %>% 
  mutate(plot_pval = -log10(plot_pval))  %>% 
  filter(pct.1 > 0.0019) %>% 
  filter(pct.2 > 0.0019) %>%
  unite("Test_ID",  c("target_guide", "target_gene"), remove = FALSE)

Full_Results %>% 
    arrange(p_val) 


##Read in detected neighbouring genes DF with additional metadata

Detc_Neighbouring_Genes <- readRDS('/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/nobackup/neighboring_genes.Rds') %>% 
unite("Test_ID",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Test_ID, .keep_all = TRUE) 

Detc_Neighbouring_Genes

Full_Results <- Full_Results %>% 
    left_join(Detc_Neighbouring_Genes, by = "Test_ID")

Full_Results

##Splitting out targeting and NTC results

Targeting_Results <- Full_Results %>% 
    filter(Targeting == "Targeting")

Targeting_Results %>% 
arrange(p_val)

##NTC results

set.seed(22)

NTC_Results <- Full_Results %>% 
    filter(Targeting == "Non-targeting") 

NTC_Results %>% 
arrange(p_val)

##Downsampling the NTCs then recombining for plotting

##NTC results

set.seed(22)

NTC_Results <- NTC_Results %>%   
    slice_sample(n = length(Targeting_Results$p_val)) 

##Combining for plotting

DS_Full_Results <- rbind(NTC_Results, Targeting_Results) 

##QQ plot comparing p-value distributions

options(repr.plot.width=12, repr.plot.height=8)

revlog_trans <- function(base = exp(1)){
  ## Define the desired transformation.
  trans <- function(x){
    log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}

ggplot(DS_Full_Results, aes(sample = plot_pval, colour = Targeting)) +
  stat_qq(size = 3) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "-log10(expected P-values)", y = "-log10(observed P-values)") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_text(colour = "black")) +
  xlim(0,4)

##Same QQ plot as above without axis text for final figures
options(repr.plot.width=12, repr.plot.height=8)

revlog_trans <- function(base = exp(1)){
  ## Define the desired transformation.
  trans <- function(x){
    log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}

ggplot(DS_Full_Results, aes(sample = plot_pval, colour = Targeting)) +
  stat_qq(size = 3) +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "", y = "") +
  theme(text = element_text(family="Arial", colour = "black", size = 30)) +
  theme(axis.text = element_blank()) +
  xlim(0,4)

##Looking at full results zoomed in

ggplot(DS_Full_Results, aes(sample = plot_pval, colour = Targeting)) +
  stat_qq(size = 2.5) +
  theme_classic() +
  scale_color_manual(values=c("#999999", "#56B4E9")) +
  theme(axis.line = element_line(colour = 'black', size = 0.4)) +
  theme(axis.ticks.y = element_line(colour = "black", size = 0.4)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  labs(title = "", x = "Theoretical p-value distribution", y = "-log10(p-value)") +
  ylim(0,50) +
  theme(text = element_text(family="Helvetica", colour = "black", size = 30)) +
  theme(axis.text = element_text(colour = "black")) +
  xlim(0,3) 
  #scale_y_continuous(trans = revlog_trans(base = 10))

##Target gene distiance vs. p values for upstream gRNAs

DS_Full_Results_Rev <- rbind(Targeting_Results, NTC_Results)

DS_Full_Results_Rev_US <- DS_Full_Results_Rev %>% 
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
