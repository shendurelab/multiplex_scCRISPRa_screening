##Loading required libraries

library(tidyverse)


K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

K562_Hits

##Count enhancer hits

Enhancer_Hits <- K562_Hits %>% 
  filter(grepl("chr", Test_ID))

Enhancer_Hits

TSS_Hits <- K562_Hits %>% 
  filter(grepl("TSS", Test_ID))

TSS_Hits

NDD_Promoter_Hits <- K562_Hits %>% 
  filter(!grepl("TSS", Test_ID)) %>% 
  filter(!grepl("chr", Test_ID))

NDD_Promoter_Hits

unique(NDD_Promoter_Hits$target_gene)

# Create Data on gRNA categories for pie chart

data <- data.frame(
  gRNA_class=as.factor(c("NDD promoter (n=40)", "TSS positive control (n=8)", "K562 candidate enhancer (n=12)")),
  value=c(40, 8, 12))

data$gRNA_class <- factor(data$gRNA_class, levels = c("NDD promoter (n=40)", "TSS positive control (n=8)", "K562 candidate enhancer (n=12)"))
data$gRNA_class

# Piechart
ggplot(data, aes(x="", y=value, fill=gRNA_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = c("#E69F00", "#00B06B", "#56B4E9", "#56B4E9")) +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels

# Create Data on gRNA categories for pie chart

data <- data.frame(
  gRNA_class=as.factor(c("Flashfry pipeline (n = 19)", "Horlbeck Library (n=29)", "K562 candidate enhancer hit (n=11)", "K562 candidate enhancer hit (n=1)")),
  value=c(19, 29, 1, 11))

data$gRNA_class <- factor(data$gRNA_class, levels = c("Flashfry pipeline (n = 19)", "Horlbeck Library (n=29)", "K562 candidate enhancer hit (n=11)", "K562 candidate enhancer hit (n=1)"))
data$gRNA_class

# Piechart
ggplot(data, aes(x="", y=value, fill=gRNA_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = c("Black", "#999999", "#56B4E9", "#56B4E9")) +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels

CRISPRaQTL_Pilot_gRNAs <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/CRISPRaQTL_pilot_guides.csv") %>% 
  rename(sequence = Spacer_Last_19) %>% 
  select(sequence, Source)

CRISPRaQTL_Pilot_gRNAs$sequence  <- toupper(CRISPRaQTL_Pilot_gRNAs$sequence)

CRISPRaQTL_Pilot_gRNAs 

Ft_Ref <- read_csv("/net/shendure/vol10/projects/ResQTL/nobackup/10X_single_cell_pilot_expt_Nov_2021/10X_pilot_1_cellranger_count_SNV_fixed/outs/crispr_analysis/feature_reference.csv")  %>% 
    rename(gRNA = id) %>% 
    filter(!target_gene_name == "Non-Targeting") %>% 
    unite("Test_ID",  c("gRNA", "target_gene_name"), remove = FALSE)

Ft_Ref

Ft_Ref <- Ft_Ref %>%
  left_join(CRISPRaQTL_Pilot_gRNAs, by = "sequence") %>% 
  rename(target_guide = gRNA)

Ft_Ref

Hits_Source <- K562_Hits %>% 
     left_join(Ft_Ref, by = "target_guide")

Hits_Source

