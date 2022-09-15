#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

(condition = args[8])
full_GWAS = args[9]
annotated_ORs = args[10]
# GW_LD_blocks = args[11]



#OUTPUT
original_base_outfilename = paste0(condition, "_", Sys.Date(), "_original_dedup_GWAS.tsv")
modified_base_outfilename = paste0(condition, "_", Sys.Date(), "_substituted_GWAS.tsv")
merged_GWAS_outfilename = paste0(condition, "_", Sys.Date(), "_merged_GWAS.tsv")
all_TS_EPs_base_outfilename = paste0(condition, "_", Sys.Date(), "_all_TS_EPs_associations.tsv")
all_TS_EPs_ZEROP_outfilename = paste0(condition, "_", Sys.Date(), "_all_TS_EPs_associations_Pdivide250.tsv")
tissue_EPeQTL_base_outfilename = paste0(condition, "_", Sys.Date(), "_only_tissue_EPeQTL_associations.tsv")
tissue_facet_base_outfilename = paste0(condition, "_", Sys.Date(), "_only_tissue_facet_associations.tsv")


original_base<-fread(full_GWAS ) 
#deduplicate by best P using data.table, much faster
(original_base<-original_base[original_base[, .I[which.min(P)], by=SNP]$V1] %>%
    dplyr::select(CHR,BP=POS,SNP,A1,A2,P,OR) %>% as_tibble())
##export original, deduplicated GWAS in right format for PRSice:
fwrite(x= original_base, file = original_base_outfilename, sep="\t")

#import annotated_ORs and format like GWAS
(annotated_OR_E_Ps <- fread(annotated_ORs) %>% 
    dplyr::select(CHR=seqnames,BP=start,SNP=ID,A1,A2=REF,P=P_recessive,OR=OR_recessive, tidytype) %>% 
    as_tibble())

#export all formatted E-Ps
all_TS_EPs <- annotated_OR_E_Ps %>% 
  dplyr::select(-tidytype)
fwrite(x= all_TS_EPs, file = all_TS_EPs_base_outfilename, sep="\t")
all_TS_EPs_zerop <- annotated_OR_E_Ps %>% 
    mutate(P=P/250) %>%
    dplyr::select(-tidytype)
fwrite(x= all_TS_EPs_zerop, file = all_TS_EPs_ZEROP_outfilename, sep="\t")

#export formatted E-P_eQTLs
(EP_eQTL <- annotated_OR_E_Ps %>% dplyr::filter(tidytype=="ESpos_Contact_overlap_eQTLs") %>% 
    dplyr::select(-tidytype))
fwrite(x = EP_eQTL, file= tissue_EPeQTL_base_outfilename, sep="\t")

#export formatted TS facet E-Ps
(facet <- annotated_OR_E_Ps %>% dplyr::filter(tidytype=="neuron_facet_not_eQTLs") %>% 
    dplyr::select(-tidytype))
fwrite(x = facet, file = tissue_facet_base_outfilename, sep="\t")

#obtain substituted GWAS file
(significant_E_Ps <-
    annotated_OR_E_Ps %>% dplyr::select(-tidytype) %>% dplyr::filter(P<0.05))
(modified_base<- 
    rbind(
      #remove an equal number of random SNPs from original base
      original_base[-sample(nrow(significant_E_Ps)),],
      #and instead, add SNPs of interest
      #setNames(
      significant_E_Ps 
      #names(original_base))
    ) %>% as.data.table())

#remove duplications
#deduplicate by best P using data.table, much faster
(modified_base<-modified_base[modified_base[, .I[which.min(P)], by=SNP]$V1] )
#export substituted GWAS file
fwrite(x = modified_base, file = modified_base_outfilename, sep="\t")


#merge annotated_OR_E_Ps and original_base, map to LD blocks, clump, separate 2 lists
#merge annotated_OR_E_Ps and original_base
(merged_GWAS<- 
    rbind(
      original_base %>% mutate(tidytype="original"),
      annotated_OR_E_Ps 
    ) %>% 
    as.data.table())

#remove duplications
#deduplicate by best P using data.table, much faster
(merged_GWAS<-merged_GWAS[merged_GWAS[, .I[which.min(P)], by=SNP]$V1] )
#export substituted GWAS file
fwrite(x = merged_GWAS, file = merged_GWAS_outfilename, sep="\t")


