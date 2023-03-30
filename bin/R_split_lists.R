#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

    # input:
    # R_split_lists.R "${ENH_list}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_EPWAS} ${EP_ES_gene_brain_exp}
(ENH_list_cohort = args[8])
clumped_SNPs = args[9]
noclump_residual_GWAS_compartment = args[10]
noclump_EPWAS = args[11]
EP_ES_gene_brain_exp = args[12]
(ES_multiplier = as.numeric(args[13]))

pdivby = 100000



# #OUTPUT
clumped_residual_GWAS_compartment_out = paste0(ENH_list_cohort, "_clumped_residual_GWAS_compartment.tsv.gz")
clumped_EPWAS_out = paste0(ENH_list_cohort,"_X_", ES_multiplier,"_clumped_EPWAS.tsv.gz")




#import clumped SNP list
(clumped_SNPs <- fread(clumped_SNPs, select=c("CHR","SNP")))

############################################################ import ES and exp data ############################################################
#import ES per enh
#  [1] "chr"                                      "start"                                    "end"                                     
#  [4] "enh"                                      "score"                                    "strand"                                  
#  [7] "log_TS_FANTOM_enh_tpm_1_4"                "log_tissue_FANTOM_enh_tpm_1_2"            "log_max_ES_perEnh_contact_1_3"           
# [10] "maxESperEnh_contact_X_TSEnhFantomExp_1_7" "highestES_gene_contact_ensemblID"         "enh_GRB_gal"                             
# [13] "enh_GRB_mus"                              "numGenes_perEnh_contact"                  "ENH_with_contact"          
(EP_ES_gene_brain_exp_info<-fread(EP_ES_gene_brain_exp)%>% 
    ############################################################ SCALE ES         ############################
    # mutate(measure1=scales::rescale(log_max_ES_perEnh_contact_1_3, to=c(1,10))) %>% 
    # elog_max_ES_perEnh_contact_X_10
    mutate(measure1= log_max_ES_perEnh_contact_1_3 * ES_multiplier)%>%
    mutate(measure2= log_TS_FANTOM_enh_tpm_1_4 * ES_multiplier)%>%
    # dplyr::filter(brain_exp_more_than_brain_median==1) %>% # N = 28100
    # dplyr::filter(brain_exp_more_than_brain_median==1 & brain_exp_more_than_other_tissues==1) %>% # N = 9176
    # dplyr::filter(brain_exp_tissue_specific==1) %>% # N = 7157
    makeGRangesFromDataFrame(keep.extra.columns = T))
seqlevelsStyle(EP_ES_gene_brain_exp_info) = "NCBI"






noclump_EPWAS_granges <- fread(noclump_EPWAS, #    // CHR	SNP	POS	A1	A2	P	OR
                                   select=c("CHR","POS","SNP","A1","A2","P","OR")) %>%
  select(chr=CHR, start=POS, end=POS, SNP, A1, A2, P, OR) %>%
  makeGRangesFromDataFrame( keep.extra.columns = T)
seqlevelsStyle(noclump_EPWAS_granges) = "NCBI"

#annotate overlapping SNPs
ES_annotated_overlaps<- mergeByOverlaps(query = noclump_EPWAS_granges,
                                        subject = EP_ES_gene_brain_exp_info, select=c("all"))
ES_annotated_overlaps$EP_ES_gene_brain_exp_info <- NULL
(ES_annotated_overlaps<-as_tibble(ES_annotated_overlaps) %>% 
    select(SNP, measure1, measure2) )


############################################################ CAN MULTIPLY P TO RESTORE ORIGINAL VALUE ############################
(clumped_EPWAS <- as_tibble(fread(noclump_EPWAS, select=c("CHR","POS","SNP","A1","A2","P","OR"))) %>%
    dplyr::mutate(P=P*pdivby) %>%
    dplyr::filter(SNP %in% clumped_SNPs$SNP) %>% left_join(ES_annotated_overlaps) %>% 
    mutate_at(vars(measure1:measure2), ~replace(., is.na(.), 1)) 
    )

#multiply OR by ES for overlapping SNPs - only if measure1 or 2 are != 1
clumped_EPWAS$OR_by_measure1 <- ifelse(
    test= clumped_EPWAS$measure1 == 1, 
    yes = clumped_EPWAS$OR, 
    no  = exp(log(clumped_EPWAS$OR) * (clumped_EPWAS$measure1))
)
clumped_EPWAS$OR_by_measure2 <- ifelse(
    test= clumped_EPWAS$measure2 == 1, 
    yes = clumped_EPWAS$OR, 
    no  = exp(log(clumped_EPWAS$OR) * (clumped_EPWAS$measure2))
)


clumped_EPWAS <- clumped_EPWAS %>%
    dplyr::select(CHR,POS,SNP,A1,A2,P,OR, OR_by_measure1, OR_by_measure2, measure1, measure2) %>%
    group_by(SNP) %>% slice_max(OR,with_ties = F) %>% ungroup()

fwrite(x =clumped_EPWAS , file = clumped_EPWAS_out, sep="\t", compress="gzip") 


#residual compartment
clumped_residual_GWAS_compartment <- as_tibble(fread(noclump_residual_GWAS_compartment)) %>% dplyr::select("CHR","POS","SNP","A1","A2","P","OR") %>%#CHR	SNP	POS	A1	A2	P	OR
    group_by(SNP) %>% slice_min(P, with_ties=F) %>% ungroup()%>%
    dplyr::filter(SNP %in% clumped_SNPs$SNP) %>%
    group_by(SNP) %>% slice_max(OR,with_ties = F) %>% ungroup()

fwrite(x =clumped_residual_GWAS_compartment , file = clumped_residual_GWAS_compartment_out, sep="\t", compress="gzip")