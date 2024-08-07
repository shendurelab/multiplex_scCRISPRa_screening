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

set.seed(21)

NTC_Results <- Full_Results %>% 
    filter(Targeting == "Non-targeting") 

NTC_Results %>% 
arrange(p_val)

##Making a dataset indicate if it was the primary/programmed/intended target gene by filtering using the Ft_Ref file

Ft_Ref <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/cellranger_input_files/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)
 
Primary_Target_Tests <- Ft_Ref$Test_ID

Primary_Target_Tests

##Filtering for the primary target tests

Primary_Targeting_Results <- Targeting_Results %>% 
 filter(Test_ID %in% Primary_Target_Tests)

Primary_Targeting_Results$Primary_Targeting <- "Primary_Target"

Other_Targeting_Results <- Targeting_Results %>% 
 filter(!Test_ID %in% Primary_Target_Tests)

Other_Targeting_Results$Primary_Targeting <- "Alternate_Target"

NTC_Results$Primary_Targeting <- "Non-targeting"

DS_Full_Results_Primary_Targets <- rbind(NTC_Results, Other_Targeting_Results, Primary_Targeting_Results) 

##Wirting targeting results to file for gviz

write_csv(Primary_Targeting_Results, "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Primary_Targeting_Results.csv")



##Heatmap showing degree of down or upregulation accross primary target tests 

options(repr.plot.width=58, repr.plot.height=1.25)


ggplot(Primary_Targeting_Results, aes(reorder(Test_ID, avg_log2FC), y = 1, fill= avg_log2FC)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("#999999", "white", "dodgerblue2"),
                       values = scales::rescale(c(-0.5, -0.25, 0, 0.05, 2.0)))
  

##Setting a 10% empirical FDR 

Total_NTC_Tests <- Full_Results %>% 
  filter(Targeting == "Non-targeting") 

NTC_Pvals <- Total_NTC_Tests$p_val

Total_NTC_Tests <- length(Total_NTC_Tests$Targeting)


Targeting_Results <- Full_Results %>% 
  filter(Targeting == "Targeting")

Targeting_Results <- Targeting_Results %>%
  rowwise() %>% 
  mutate(n_lower_NTC_pvals = sum(NTC_Pvals < p_val))


##Defining empirical FDR 

Targeting_Results <- Targeting_Results %>%
  mutate(Empirical_p_val = ((n_lower_NTC_pvals+1)/(Total_NTC_Tests+1)))


##BH correcting the pvals

BH_Corrected_Empirical_p_val = p.adjust(Targeting_Results$Empirical_p_val, method = "BH")

Targeting_Results$BH_Corrected_Empirical_p_val <- BH_Corrected_Empirical_p_val

Hits <- Targeting_Results %>% 
  filter(BH_Corrected_Empirical_p_val < 0.1)

Enhancer_Hits <- Hits %>% 
  filter(grepl("chr", Test_ID))

Primary_Targeting_Results <- Targeting_Results %>% 
 filter(Test_ID %in% Primary_Target_Tests)

Hits

##Writing hits to a file 

write_csv(Hits, "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")



##Determining if any gRNAs produced upregulation of more than one gene

Hits %>% 
  group_by(target_guide) %>% 
  filter(n()>1)

##Determining how many promoter hits upregulated their expected target gene

Hits %>%
  filter(!grepl("chr", Test_ID)) %>% 
  filter(Test_ID %in% Primary_Targeting_Results$Test_ID)

##Determining hot many hits upregulated an alternate target gene

Hits %>% 
  filter(!grepl("chr", Test_ID)) %>%
  filter(!Test_ID %in% Primary_Targeting_Results$Test_ID)

##Determining how many enhancer hits there are

Enhancer_Hits <- Hits %>% 
  filter(grepl("chr", Test_ID))

Enhancer_Hits

##Determining how many enhancer hits upregulated their expected target gene

Enhancer_Hits %>% 
  filter(Test_ID %in% Primary_Targeting_Results$Test_ID)

##Determining how many enhancer hits upregulated an alternate gene

Enhancer_Hits %>% 
  filter(!Test_ID %in% Primary_Targeting_Results$Test_ID)

##BH corrected empircal FDR < 0.1 significance threshold heatmap

Sig_Primary_Targeting_Results <- Primary_Targeting_Results

Sig_Primary_Targeting_Results$BH_Corrected_Empirical_p_val[Sig_Primary_Targeting_Results$BH_Corrected_Empirical_p_val < 0.1] <- 0.1
Sig_Primary_Targeting_Results$BH_Corrected_Empirical_p_val[Sig_Primary_Targeting_Results$BH_Corrected_Empirical_p_val > 0.1] <- 0


options(repr.plot.width=58, repr.plot.height=1.25)

ggplot(Sig_Primary_Targeting_Results, aes(reorder(Test_ID, avg_log2FC), y = 1, fill= BH_Corrected_Empirical_p_val)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("white", "black"))

##P-value < e-10 significance threshold heatmap

Sig_Primary_Targeting_Results <- Primary_Targeting_Results

Sig_Primary_Targeting_Results$plot_pval[Sig_Primary_Targeting_Results$plot_pval < 10] <- 0
Sig_Primary_Targeting_Results$plot_pval[Sig_Primary_Targeting_Results$plot_pval > 10] <- 1

options(repr.plot.width=58, repr.plot.height=1.25)

ggplot(Sig_Primary_Targeting_Results, aes(reorder(Test_ID, avg_log2FC), y = 1, fill= plot_pval)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("white", "black"))
