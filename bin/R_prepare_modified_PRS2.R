#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
# library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

(condition = args[8])
clumped_SNPs = args[9]
merged_GWAS = args[10]



# #OUTPUT
#split scores
Original_NoEPs_base_outfilename = paste0(condition,  "_", Sys.Date(),"_originalGWAS_SNPs_clumped_NoEPs_overlap.tsv")
NonOverlapOriginal_TS_EPs_base_outfilename = paste0(condition,  "_", Sys.Date(),"_TS_EPs_GWAS_REC_clumped_NonOverlapOriginal.tsv")


clumped_SNPs <- fread(clumped_SNPs, select=c("CHR","SNP"))

merged_GWAS <- fread(merged_GWAS, select=c("CHR","BP","SNP","A1","A2","P","OR","tidytype"))

clumped_merged_GWAS <- merged_GWAS %>%
    dplyr::filter(SNP %in% clumped_SNPs$SNP)

Original_NoEPs_base <- clumped_merged_GWAS %>%
    dplyr::filter(tidytype=="original") %>%
    dplyr::select(-tidytype)

NonOverlapOriginal_TS_EPs_base <- clumped_merged_GWAS %>%
    dplyr::filter(tidytype!="original") %>%
    dplyr::select(-tidytype)

fwrite(x = Original_NoEPs_base, file = Original_NoEPs_base_outfilename, sep="\t")
fwrite(x = NonOverlapOriginal_TS_EPs_base, file = NonOverlapOriginal_TS_EPs_base_outfilename, sep="\t") 