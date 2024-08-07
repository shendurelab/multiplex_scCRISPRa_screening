##Change working directory

setwd("/net/shendure/vol10/projects/troym/ResQTL/nobackup/promoter_pilot/K562_Pilot_gRNA_Coordinates/")

getwd()

##Loading libraries and BSgenomes

library(tidyverse)
library(BSgenome)
library(AnnotationHub)
library(corrplot)
library(ggridges)
library(ggrepel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(Biostrings)

available.genomes()
installed.genomes()

##Mapping a list of gRNAs to genomic coordinates. 

##First read in gRNAs as FASTA

x <- readDNAStringSet("/net/shendure/vol10/projects/troym/ResQTL/nobackup/promoter_pilot/K562_Pilot_gRNA_Coordinates/CRISPRaQTL_K562_Pilot_gRNAs.fa", "fasta")
x


writeHits <- function(seqname, matches, strand, file="", append=FALSE)
{
if (file.exists(file) && !append)
warning("existing file ", file, " will be overwritten with 'append=FALSE'")
if (!file.exists(file) && append)
warning("new file ", file, " will have no header with 'append=TRUE'")
hits <- data.frame(seqname=rep.int(seqname, length(matches)),
start=start(matches),
end=end(matches),
strand=rep.int(strand, length(matches)),
patternID=names(matches),
check.names=FALSE)
write.table(hits, file=file, append=append, quote=FALSE, sep="\t",
row.names=FALSE, col.names=!append)
}
runAnalysis1 <- function(dict0, outfile="")
{
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
seqnames <- seqnames(genome)
seqnames_in1string <- paste(seqnames, collapse=", ")
cat("Target:", metadata(genome)$genome,
"chromosomes", seqnames_in1string, "\n")
append <- FALSE
for (seqname in seqnames) {
subject <- genome[[seqname]]
cat(">>> Finding all hits in chromosome", seqname, "...\n")
for (i in seq_len(length(dict0))) {
patternID <- names(dict0)[i]
pattern <- dict0[[i]]
plus_matches <- matchPattern(pattern, subject)
names(plus_matches) <- rep.int(patternID, length(plus_matches))
writeHits(seqname, plus_matches, "+", file=outfile, append=append)
append <- TRUE
rcpattern <- reverseComplement(pattern)
minus_matches <- matchPattern(rcpattern, subject)
names(minus_matches) <- rep.int(patternID, length(minus_matches))
writeHits(seqname, minus_matches, "-", file=outfile, append=append)
}
cat(">>> DONE\n")
}
}

runAnalysis1(x, outfile="/net/shendure/vol10/projects/troym/ResQTL/nobackup/promoter_pilot/K562_Pilot_gRNA_Coordinates/CRISPRaQTL_K562_pilot_hg38_gRNA_Coordinates.txt")

##This works for one set of cooridnates

mart <- useEnsembl(biomart = 'genes', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 90)
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start","end")
values <- list(9, 72644043, 73644043)
all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)

## Reading in our list of gRNA coordinates 

gRNA_Coordinates <- read_tsv("CRISPRaQTL_K562_pilot_hg38_gRNA_Coordinates.txt") %>% 
  dplyr::rename(Chromosome = seqname, gRNA_Name = patternID) %>% 
  separate(Chromosome, into = c("Chr", "Chromosome"),
            sep= 3) 

##Adding a unique identifier for each gRNA mapping, and calculating 1mb before and after each gRNA

gRNA_Coordinates <- gRNA_Coordinates %>% 
  group_by(gRNA_Name) %>%
  mutate(Unique_gRNA_mapping_number = row_number()) %>% 
  unite("Unique_gRNA_mapping_name",  c("gRNA_Name", "Unique_gRNA_mapping_number"), remove = FALSE)

gRNA_Coordinates <- gRNA_Coordinates %>% 
  mutate(Start_1mb_Window = start - 1000000) %>%
  mutate(End_1mb_Window = end + 1000000) 


##Now loopin through to find all genes within 1mb of all guides 

Genes_1MB <- NULL;
for (i in c(gRNA_Coordinates$Unique_gRNA_mapping_name)) {
print(paste0('Finding genes near ', i))
Chromosome <- gRNA_Coordinates %>% 
  filter(Unique_gRNA_mapping_name == i)
Chromosome <- Chromosome$Chromosome
Start_1mb_Window <- gRNA_Coordinates %>% 
  filter(Unique_gRNA_mapping_name == i)
Start_1mb_Window <- Start_1mb_Window$Start_1mb_Window
End_1mb_Window <- gRNA_Coordinates %>% 
  filter(Unique_gRNA_mapping_name == i)
End_1mb_Window <- End_1mb_Window$End_1mb_Window
mart <- useMart("ensembl", host="mar2017.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start","end")
values <- list(Chromosome, Start_1mb_Window, End_1mb_Window)
all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes <- add_row(all.genes, ensembl_gene_id = NA)
all.genes$Unique_gRNA_mapping_name <- i
print(all.genes)
Genes_1MB <- rbind(Genes_1MB, all.genes)
} 

##Writing results to file

write_csv(Genes_1MB, "CRISPRaQTL_K562_pilot_hg38_gRNA_neighbouring_genes_1MB2.csv")




##Reading in results

Genes_1MB <- read_csv("CRISPRaQTL_K562_pilot_hg38_gRNA_neighbouring_genes_1MB.csv")


##left joining with gRNA coordinates 

Genes_1MB <- Genes_1MB %>% 
  left_join(gRNA_Coordinates, by = "Unique_gRNA_mapping_name") 

##Rename columns to make them more explicit

Genes_1MB <- Genes_1MB %>% 
  dplyr::rename(target_gene_strand = strand.x, gRNA_strand = strand.y)

##Calculating distance of each gene to gRNA

Sense_Genes_1MB <- Genes_1MB %>% 
  filter(target_gene_strand == "1") %>% 
  mutate(target_gene_distance = start - start_position)

Antisense_Genes_1MB <- Genes_1MB %>% 
  filter(target_gene_strand == "-1") %>% 
  mutate(target_gene_distance = end_position - start)

Genes_1MB <- rbind(Sense_Genes_1MB, Antisense_Genes_1MB) %>% 
  arrange(Unique_gRNA_mapping_name)

##Filter to only keep upstream genes within 500kb (or arbitrary window)

Genes_1MB <- Genes_1MB %>% 
  filter(target_gene_distance < 1000001)

##Only keeping distinct genes by HGNC symbol (what Seurat uses for testing) and writing to file


Distinct_Genes_1MB <- Genes_1MB %>% 
  group_by(Unique_gRNA_mapping_name) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  filter(!hgnc_symbol == "NA")



Distinct_Genes_1MB

##Writing results to file

write_csv(Distinct_Genes_1MB, "CRISPRaQTL_K562_pilot_hg38_gRNA_distinct_neighbouring_genes_1MB2.csv") 


Distinct_Genes_1MB <- read_csv("CRISPRaQTL_K562_pilot_hg38_gRNA_distinct_neighbouring_genes_1MB.csv")

Distinct_Genes_1MB

##Now I want to pair all potentially targeted genes with those same NTCs 

##First find the list of the 50 NTCs that are truly NTCs

NTC_gRNAs <- read_csv("/net/shendure/vol10/projects/ResQTL/nobackup/10X_single_cell_pilot_expt_Nov_2021/10X_pilot_1_cellranger_count_SNV_fixed/outs/crispr_analysis/feature_reference.csv") %>% 
  dplyr::rename(gRNA_Name = id) %>% 
  filter(grepl("Horlbeck_NTC", gRNA_Name))

##Now writing a loop to artificially pair each of the 50 NTCs with all genes within 1mb of a targeting guide. 

##First making a modifiable version of the list of genes within 1 megabase as an output for the NTC assignment loop

Distinct_Genes_1MB_NTC <- Distinct_Genes_1MB

Genes_1MB_NTC <- NULL;
for (i in c(NTC_gRNAs$gRNA_Name)) {
print(paste0('Artificially assigning target genes to ', i))
Distinct_Genes_1MB_NTC$gRNA_Name <- i
print(Distinct_Genes_1MB_NTC)
Genes_1MB_NTC <- rbind(Distinct_Genes_1MB_NTC, Genes_1MB_NTC)
}

Distinct_Genes_1MB_NTC <- Genes_1MB_NTC


Distinct_Genes_1MB_NTC


Distinct_Genes_1MB_NTC <- Distinct_Genes_1MB_NTC %>% 
  unite("Unique_gRNA_mapping_name",  c("gRNA_Name", "hgnc_symbol"), remove = FALSE) %>% 
  distinct(Unique_gRNA_mapping_name, .keep_all = TRUE)

Distinct_Genes_1MB_NTC




##Writing results to file

write_csv(Distinct_Genes_1MB_NTC, "CRISPRaQTL_K562_pilot_hg38_gRNA_distinct_neighbouring_genes_1MB_NTC_pairings.csv")

##Combining the targeting genes and NTC gene lists within 1MB for testing 

Distinct_Genes_1MB$Targeting <- "Targeting"
Distinct_Genes_1MB_NTC$Targeting <- "Non-targeting"

All_Distinct_1MB_gene_gRNA_pairings_to_test <- rbind(Distinct_Genes_1MB, Distinct_Genes_1MB_NTC)

write_csv(All_Distinct_1MB_gene_gRNA_pairings_to_test, "CRISPRaQTL_K562_pilot_hg38_gRNA_All_distinct_neighbouring_genes_1MB_to_test.csv")

##Reading in full list of tests 

All_Distinct_1MB_gene_gRNA_pairings_to_test <- read_csv("CRISPRaQTL_K562_pilot_hg38_gRNA_All_distinct_neighbouring_genes_1MB_to_test.csv")

##How many targeting guides are there?

Distinct_Genes_1MB %>%
  ungroup() %>% 
  distinct(gRNA_Name)

##How many distinct mappings? Should be less than or equal to the number of gRNA coordinates

Distinct_Genes_1MB %>% 
  distinct(Unique_gRNA_mapping_name)

##Note there are 443 unique guides, exactly what we would expect with 50 true NTCs, many of which map to multiple regions in the genome - the true set of 50 NTCs will have to be the new NTCs for this which couldn't be done with the 10x analysis. 

gRNA_Coordinates %>% 
  group_by(gRNA_Name) %>% 
  filter(n()>1)

##See how many genes a given gRNA is typically paired with

N_Distinct_1MB_Genes_Per_Spacer <- Distinct_Genes_1MB %>% 
  dplyr::count(gRNA_Name)



N_Distinct_1MB_Genes_Per_Spacer %>% 
arrange(-n)

##See how many genes a given gRNA is typically paired with

options(repr.plot.width=7, repr.plot.height=6)

ggplot(N_Distinct_1MB_Genes_Per_Spacer, aes(x = n)) +
   geom_histogram(color="#999999", fill="#999999", binwidth = 0.5) +  
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  labs(title = "", x = "Number of genes within 1 Mb of each gRNA", y = "Count") +
  theme(text = element_text(family="Arial", colour = "black", size = 24)) 

##Visualization of gRNA distance to neighbouring genes

ggplot(Distinct_Genes_1MB, aes(x = target_gene_distance)) +
   geom_histogram(color="#999999", fill="#999999", bins = 100) +  
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 0.7)) +
  theme(axis.ticks = element_line(colour = "black", size = 0.7)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text = element_text(family="Arial", colour = "black", size = 24)) +
  labs(title = "", x = "Distance of gRNAs to neighbouring genes (bp)", y = "Count") +
  theme(text = element_text(family="Arial", colour = "black", size = 24))  
