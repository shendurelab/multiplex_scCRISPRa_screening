##Loading libraries 

library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(tidyverse)

##Making a visualization of te effects of each gRNA at a given locus (ANK2 K562)

symbol = "ANK2"
gen = "hg38"
chr = "chr4"
from = 112700000
to = 113500000


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track

gtrack <- GenomeAxisTrack()


##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(with_refseq_mrna=TRUE, ensembl_transcript_id= c("ENST00000672356", "ENST00000671951", "ENST00000671727", "ENST00000672240", "ENST00000672411")),
                                    biomart = bm) 


##Making a track for the location of each regulatory element

RE_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gencode.v38.basic.coding.transcripts.500bp_promoters.CRISPRa_QTL_PILOT.txt", col_names = FALSE) %>% 
  dplyr::rename(Chromosome = X1, Start = X2, End = X3, Gene_Transcript = X4, Strand = X6) %>% 
  separate(Gene_Transcript, into = c("Gene", "Transcripts"), sep= ";") %>% 
  filter(Gene == symbol)

start <- RE_Coordinates$Start
end <- RE_Coordinates$End

REtrack <- AnnotationTrack(start = start,
                    end = end, chromosome = chr, genome = gen, 
                    name = "Candidate Regulatory Elements")

##Then reading in the gRNA mappings and functional data and generating the gRNA score track

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

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,308))
Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,308))

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track, Hit_Track), legend= TRUE)


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



##Adding a H3K27ac track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/ENCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,5))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(overlapTracks) <- list(col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.3))  

##Making a visualization of te effects of each gRNA at a given locus (CHD8 K562)

symbol = "CHD8"
gen = "hg38"
chr = "chr14"
from = 21380000
to = 21465000


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

RE_Coordinates <- read_tsv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/gRNA_design/gencode.v38.basic.coding.transcripts.500bp_promoters.CRISPRa_QTL_PILOT.txt", col_names = FALSE) %>% 
  dplyr::rename(Chromosome = X1, Start = X2, End = X3, Gene_Transcript = X4, Strand = X6) %>% 
  separate(Gene_Transcript, into = c("Gene", "Transcripts"), sep= ";") %>% 
  filter(Gene == symbol)

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

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,16))
Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,16))

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track, Hit_Track), legend= TRUE)


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
plotTracks(ATAC_track, from = from, to = to)



##Adding a H3K27AC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/K562_data/K562_ATAC_H3K27ac/ENCFF381NDD.bigWig"
H3k27ac <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "#56B4E9", 
                    fill.histogram = "#56B4E9", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,15))
class(H3k27ac)

H3k27ac


##Changing the display features

displayPars(overlapTracks) <- list(col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.15)) 
