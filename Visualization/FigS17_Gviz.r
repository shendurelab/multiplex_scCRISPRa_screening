##Loading required libraries

library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(tidyverse)

##ANXA1 K562

symbol = "ANXA1"
gen = "hg38"
chr = "chr9"
from = 73140000
to = 73180000


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chr9.871") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,40), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")


overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/ENCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,15))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 


##TMEM56 K562

symbol = "TMEM56"
gen = "hg38"
chr = "chr1"
from = 95109054
to = 95152908


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = "TLCD4",
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= "TLCD4"),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chr1.7358") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,35), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")

overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,5))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/NCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2))

##TSMB4X K562


symbol = "TMSB4X"
gen = "hg38"
chr = "chrX"
from = 12952580
to = 12980000


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chrX.232") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,30), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")


overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/ENCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,20))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 

##ASIC1 K562

symbol = "ASIC1"
gen = "hg38"
chr = "chr12"
from = 50035092
to = 50090115


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == "TUBA1A") %>% 
  filter(Target_Site == "chr12.1559") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- K562_Hits %>% 
  filter(target_gene == symbol)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,25), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")

overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/ENCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,20))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 

##ANXA1 neuron

symbol = "ANXA1"
gen = "hg38"
chr = "chr9"
from = 73140000
to = 73180000


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chr9.871") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "white", legend=FALSE, ylim = c(0,40), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")

overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_ATAC-seq_Song_2019_hg38.50bp_trim.merged.srt.nodup.no_chrM_MT.filtered.bw"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_H3K27ac_CUT+RUN_Song_2019_hg38.srt.no_dups.filtered.bw"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,15))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 

##TMEM56 neuron

symbol = "TMEM56"
gen = "hg38"
chr = "chr1"
from = 95109054
to = 95152908


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = "TLCD4",
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= "TLCD4"),
                                    biomart = bm)  


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chr1.7358") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits


Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,35), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_ATAC-seq_Song_2019_hg38.50bp_trim.merged.srt.nodup.no_chrM_MT.filtered.bw"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,5))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_H3K27ac_CUT+RUN_Song_2019_hg38.srt.no_dups.filtered.bw"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2))

##TSMB4X neuron

symbol = "TMSB4X"
gen = "hg38"
chr = "chrX"
from = 12952580
to = 12980000


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == symbol) %>% 
  filter(Target_Site == "chrX.232") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol) 

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% K562_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,30), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_ATAC-seq_Song_2019_hg38.50bp_trim.merged.srt.nodup.no_chrM_MT.filtered.bw"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_H3K27ac_CUT+RUN_Song_2019_hg38.srt.no_dups.filtered.bw"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,20))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 

##ASIC1 neuron

symbol = "ASIC1"
gen = "hg38"
chr = "chr12"
from = 50035092
to = 50090115

##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, hgnc_symbol= symbol),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/Gasperini_664_Enhancer_Gene_Pairs_hg38_coords_0505_2022.csv", col_names = TRUE) %>% 
  dplyr::rename(Gene = target_gene_short, Start = start.candidate_enhancer, End = stop.candidate_enhancer) %>% 
  filter(Gene == "TUBA1A") %>% 
  filter(Target_Site == "chr12.1559") 
  

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements") 

##Then reading in the gRNA mappings and function data and generating the gRNA function score track

gRNA_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) %>%
  arrange(gRNA_Name) 

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

K562_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/differential_gene_expression_testing_outputs/K562_1Mb_DE/K562_Hits.csv")

Hits <- K562_Hits %>% 
  filter(target_gene == symbol)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% K562_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "white", legend=FALSE, ylim = c(0,25), col.axis = "black", background = "white", fontcolor.title = "black", name = "-log10(P-values)", background.title = "white", col.frame = "white", fill = "white")


overlapTracks <- OverlayTrack(trackList = list(Hit_Track), legend= TRUE, background.title = "white", col.frame = "white", fill = "white")

##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_ATAC-seq_Song_2019_hg38.50bp_trim.merged.srt.nodup.no_chrM_MT.filtered.bw"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,10))
class(ATAC_track)

ATAC_track



##Adding a H3k27Ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_H3K27ac_CUT+RUN_Song_2019_hg38.srt.no_dups.filtered.bw"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "H3K27AC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,20))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")
displayPars(H3k27ac) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.2)) 
