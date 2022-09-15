#!/usr/bin/env Rscript

## PART ONE: IMPORT SNV AND GWAS FILES, MAKE BED FILES
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs()

print(args)
#inputs
significant_SNVs = args[8]
GWAS = args[9]
chain_38_19 = args[10]
condition = args[11]
#outputs
(SNVoutput = paste0(
  "ENH_",
  condition,
  "_hg19.bed"))
(formattedSNVoutput = paste0(
  "ENH_",
  condition,
  "_hg19.csv"))
(GWAS_SNVoutput =  paste0(
  "GWAS_",
  condition,
  "_SNV_merge_hg19.bed"))
  

# bed.dir = "bed"
# dir.create(file.path(".", bed.dir))

#import liftOver chain file
ch_38_19 = import.chain(chain_38_19)

#IMPORT SNV VARIANT FILES AND CONVERT TO HG19
#BRAIN
(significant_SNVs_hg38 <- makeGRangesFromDataFrame(
  data.table::fread(significant_SNVs, sep=","),
  keep.extra.columns = T))

seqlevelsStyle(significant_SNVs_hg38) = "UCSC"
significant_SNVs_hg19 <- liftOver(significant_SNVs_hg38, ch_38_19) %>%  unlist() %>% 
  as_tibble() 

# head(significant_SNVs_hg19)
print("significant_SNVs_hg19")
(significant_SNVs_hg19_bed<- significant_SNVs_hg19 %>% 
  dplyr::select(seqnames, start, end))

#print out bed file

export.bed(object=GenomicRanges::makeGRangesFromDataFrame(significant_SNVs_hg19),
           con= SNVoutput) 
data.table::fwrite(significant_SNVs_hg19, file=formattedSNVoutput)

#IMPORT GWAS FILE AND CONVERT TO HG19
(GWAS_clumped_hg19 <- GenomicRanges::makeGRangesFromDataFrame(
     data.table::fread(GWAS, sep="\t"),
    keep.extra.columns = F, start.field = "POS", end.field = "POS", seqnames.field = "CHR"
))

seqlevelsStyle(GWAS_clumped_hg19)="UCSC"

GWAS_clumped_hg19_bed <- as_tibble(GWAS_clumped_hg19) %>% 
  dplyr::select(seqnames, start, end)


### PART 2: MERGE SNV AND GWAS datasets, then merge


(SNV_GWAS <- GenomicRanges::makeGRangesFromDataFrame(
    rbind(significant_SNVs_hg19_bed,GWAS_clumped_hg19_bed)))


export.bed(object=SNV_GWAS, 
            con=GWAS_SNVoutput )