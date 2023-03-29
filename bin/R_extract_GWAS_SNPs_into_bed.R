#!/usr/bin/env Rscript

## PART ONE: IMPORT SNV AND GWAS FILES, MAKE BED FILES
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs()

print(args)
#inputs 

print(args[9])

GWAS_QC_nodups <- fread(file = args[8])
collected_bed_files_for_enhancers <- read_lines(args[9])
clumped_SNPs <-fread(file = args[10], select=c("CHR","SNP"))
condition = args[11]

#outputs
(SNPs_to_extract_out <- paste0(condition,"_clumped_GWAS_SNPs_plus_those_in_bed_files.bed"))
clumped_GWAS_out <- paste0(condition,"_clumped_GWAS_QC_nodups.tsv.gz")





#read each SNP list and join to GWAS
totalbed = data.frame()
for(bedfile in collected_bed_files_for_enhancers){
  print(bedfile)
  (bed = fread(file = bedfile, select=c(1:4), col.names = c("seqnames","start","end","name")) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T))
  seqlevelsStyle(x = bed) <- "UCSC"
  (bed=as_tibble(bed))
  
  (totalbed = rbind(totalbed,bed))
  
}
sample_n(totalbed, 20)

#remove internal range dupls
(totalbed = makeGRangesFromDataFrame(totalbed, keep.extra.columns = T) %>% 
    GenomicRanges::reduce(drop.empty.ranges = T))

#extract GWAS SNPs wihtin ranges
print("GWAS head")
print(head(GWAS_QC_nodups))
(GWAS_QC_nodups_GR = makeGRangesFromDataFrame(GWAS_QC_nodups, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(GWAS_QC_nodups_GR) <- "UCSC"
(full_PGC_GWAS_overlap_beds = subsetByOverlaps(x = GWAS_QC_nodups_GR, ranges = totalbed, type="any"))

##clumped GWAS
clumped_GWAS_QC_nodups <- GWAS_QC_nodups %>% dplyr::filter(SNP %in% clumped_SNPs$SNP)
fwrite(clumped_GWAS_QC_nodups, file=clumped_GWAS_out, sep = "\t", compress = "gzip")
(clumped_GWAS_QC_nodups_GR = makeGRangesFromDataFrame(clumped_GWAS_QC_nodups, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(clumped_GWAS_QC_nodups_GR) <- "UCSC"

(SNPs_to_extract = rbind(
  as_tibble(full_PGC_GWAS_overlap_beds),
  as_tibble(clumped_GWAS_QC_nodups_GR)
))


(SNPs_to_extract <- SNPs_to_extract %>% 
        group_by(SNP) %>% slice_head(n=1) %>% ungroup())

#make bed
SNPs_to_extract %>% select(seqnames,  start, SNP) %>%
  mutate(POS=start, score=".", strand=".") %>%
  relocate(seqnames, start, POS, SNP, score, strand) %>%
  fwrite(file=SNPs_to_extract_out, sep="\t", col.names=F)
