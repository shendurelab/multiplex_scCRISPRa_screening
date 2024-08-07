##Loading required libraries

library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(tidyverse)

##Making a visualization of the effects of each gRNA at each target locus (TCF4 iPSC-derived neurons)

symbol = "TCF4"
gen = "hg38"
chr = "chr18"
from = 55187714
to = 55670157


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track 

gtrack <- GenomeAxisTrack()


##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(ensembl_transcript_id= c("ENST00000561831", "ENST00000637169", "ENST00000643689", "ENST00000544241", "ENST00000561992", "ENST00000543082", "ENST00000354452", "ENST00000398339")),
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

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

Neuron_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% Neuron_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% Neuron_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,39))
Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,39))

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track, Hit_Track), legend= TRUE)


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

displayPars(overlapTracks) <- list(col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.3)) 


##Making a visualization of the effects of each gRNA at each target locus (TBR1 iPSC-derived neurons)

symbol = "TBR1"
gen = "hg38"
chr = "chr2"
from = 161413904
to = 161426264


##First loading a biomart dataset from ensembl

bm=biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

##Creating a genome axis track and ideogram

gtrack <- GenomeAxisTrack()

##Then constructing biomart gene annotation track for a given gene

biomTrack <- BiomartGeneRegionTrack(genome = gen, chromosome = chr, 
                                    name = "ENSEMBL",
                                    symbol = symbol,
                                    filter = list(ensembl_transcript_id= c("ENST00000389554", "ENST00000410035")),
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

gRNA_DE_Results <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Primary_Targeting_Results.csv") 

gRNA_Data <- gRNA_DE_Results %>% 
  filter(target_gene == symbol)

##Filtering for hits or non-hits and coloring 

Neuron_Hits <- read_csv("/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/differential_gene_expression_testing_outputs/Neuron_1Mb_DE/Neuron_Hits.csv")

Hits <- gRNA_Data %>% 
  filter(target_guide %in% Neuron_Hits$target_guide)
Hits$Hit <- "Hit"

Non_Hits <- gRNA_Data %>% 
  filter(!target_guide %in% Neuron_Hits$target_guide)
Non_Hits$Hit <- "Non-Hit" 

##Create hit and non-hit tracks and overlay them

df1_label = "Non-Hit"
df2_label = "Hit"
df1 <- Non_Hits
df2 <- Hits

Hit_Track <- DataTrack(data = Hits$plot_pval, start = Hits$start, end = Hits$end, chromosome = chr, genome = gen, groups = factor(df2_label,levels = c(df1_label, df2_label)), col = "#56B4E9", legend=FALSE, ylim = c(0,160))
Non_Hit_Track <- DataTrack(data = Non_Hits$plot_pval, start = Non_Hits$start, end = Non_Hits$end, chromosome = chr, genome = gen, groups = factor(df1_label,levels = c(df1_label, df2_label)), col = "#999999", legend=FALSE, ylim = c(0,160))

overlapTracks <- OverlayTrack(trackList = list(Non_Hit_Track, Hit_Track), legend= TRUE)


##Adding an ENCODE ATAC track


bamFile <- "/net/shendure/vol10/www/content/members/CRISPRa_QTL_website/public/data/neuron_data/Neuron_ATAC_H3K27ac/WTC11_NGN2_iPSC_7-8wk_ExN_ATAC-seq_Song_2019_hg38.50bp_trim.merged.srt.nodup.no_chrM_MT.filtered.bw"
ATAC_track <- DataTrack(range = bamFile, genome = gen, type = "l", 
                     name = "ATAC", window = -1, 
                     chromosome = chr, type = "hist", 
                    col.histogram = "black", 
                    fill.histogram = "black", background.title = "white", 
                    fontcolor.title = "black", ylim = c(0,15))
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

displayPars(overlapTracks) <- list(col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(REtrack) <- list(fill = "darkgray", col.axis = "black", background.title = "white", fontcolor.title = "black")
displayPars(biomTrack) <- list(background.title = "white", fontcolor.title = "black", col = "#999999", fill = "#999999", col.line = "#999999", utr3 = "#999999", utr5 = "#999999", protein_coding = "#999999")
displayPars(ATAC_track) <- list(background.title = "white", fontcolor.title = "black")


##Plotting the tracks

pdf("", width = 8, height = 4)
plotTracks(list(overlapTracks, REtrack, H3k27ac, ATAC_track, biomTrack), stackHeight = 0.3, transcriptAnnotation = "symbol", from = from, to = to, sizes=c(0.5,0.1,0.1,0.1,0.3)) 

