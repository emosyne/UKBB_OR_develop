#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

(condition = args[8])
PLINK_ORs_recessive = args[9]
PLINK_ORs_standard = args[10]
processed_ENH_SNP_lists_hg19 = args[11]
GWAS = args[12]
hg38ToHg19_chain = import.chain(args[13])
GW_LD_blocks = args[14]


(results_per_snp <- rbind(
  cbind(data.table::fread(PLINK_ORs_recessive), recessive=rep("recessive")),
  cbind(data.table::fread(PLINK_ORs_standard), recessive=rep("standardGLM"))
))

#OUTPUT
annotated_ORs_outfilename = paste0(condition, "_", Sys.Date(), "_annotated_ORs.csv")
dir.create(file.path(".", "figs"))
ORfigure = paste0("figs/PLINK_ORs_",condition, "_", Sys.Date(),"_UKBB.png")

## Annotate with  INFO scores, MAF


##if dir exists, annotate info scores and MAF
if (dir.exists(file.path("/mnt/biggley/home/emanuele/2HH/biobank/genetics/imputed_data/qualityscores/"))){
  # list all files in folder: one per chr
  files <- dir("/mnt/biggley/home/emanuele/2HH/biobank/genetics/imputed_data/qualityscores/", pattern = "^ukb_mfi_chr", full.names = TRUE, ignore.case = TRUE)
  results_per_snp2 <-results_per_snp
  for (file in files) {
    print(file)
    infoscores <- data.table::fread(file, sep="\t", 
                                    select = c(2,6,7,8), col.names = c("SNP", "MAF","minor_allele","INFO"), 
                                    fill = T, showProgress = F)
    results_per_snp2<-left_join(results_per_snp2, infoscores, by=c("ID"="SNP"))
  }
  
  results_per_snp2<-results_per_snp2 %>%  pivot_longer(
    cols = starts_with("MAF"),
    names_to = "pippo",
    values_to = "MAF",
    values_drop_na = T
  ) %>%  pivot_longer(
    cols = starts_with("INFO"),
    names_to = "pippo2",
    values_to = "INFO",
    values_drop_na = T
  ) %>%  pivot_longer(
    cols = starts_with("minor_allele"),
    names_to = "pippo3",
    values_to = "minor_allele",
    values_drop_na = T
  ) %>% dplyr::select(-starts_with(match = "pippo"))
} else {
  #if no info, just create NA values
  results_per_snp2 <- results_per_snp %>%
    mutate(MAF=NA, INFO=NA)
}


results_per_snp2


# filter unique SNVs and then unique E-P_eQTLs by joining BB ORs back to E-P and eQTL SNV data 


results_per_snp3 <- as_tibble(results_per_snp2, .name_repair = "universal") %>% 
  dplyr::mutate(seqnames=.CHROM,start=POS, end=POS, .CHROM=NULL, POS=NULL) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
seqlevelsStyle(results_per_snp3) = "UCSC"
results_per_snp3


#NEED TO MERGE ORS WITH SNV DATA
significant_E_Ps_hg19 <- data.table::fread(processed_ENH_SNP_lists_hg19)
table(significant_E_Ps_hg19$tidytype)

#FIRST MERGE THOSE OVERLAPPING GTEX EQTLS
significant_SNVs_withinGTEx_hg38 <- significant_E_Ps_hg19 %>% 
  dplyr::filter(tidytype=="ESpos_Contact_overlap_eQTLs") %>% 
  #make granges with GTEx SNV position
  dplyr::select(-start ,-end, start=GTEx_SNV_pos_hg38) %>% 
  dplyr::mutate(end=start) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

seqlevelsStyle(significant_SNVs_withinGTEx_hg38) = "UCSC"

significant_SNVs_withinGTEx_hg19 <- liftOver(significant_SNVs_withinGTEx_hg38, hg38ToHg19_chain) %>%  
  unlist() %>%  as_tibble()



# (as_tibble(results_per_snp3))

significant_SNVs_withinGTEx_hg19 <- 
  inner_join(
    as_tibble(results_per_snp3) %>% dplyr::select(-width, -strand) ,
    as_tibble(significant_SNVs_withinGTEx_hg19)%>% dplyr::select(-width, -strand) ,
    by = c("seqnames", "start", "end")
  ) 



#THEN MERGE THOSE not overlapping eqtls
significant_E_Ps_notGTEx_hg19 <- significant_E_Ps_hg19 %>% 
  dplyr::filter(tidytype!="ESpos_Contact_overlap_eQTLs") %>% 
  #make granges with enh position
  dplyr::select(-GTEx_SNV_pos_hg38) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
seqlevelsStyle(significant_E_Ps_notGTEx_hg19) = "UCSC"
#find overlaps between ENH and ORs
overlaps<-mergeByOverlaps(query = results_per_snp3, subject = significant_E_Ps_notGTEx_hg19) 

significant_unique_E_Ps_notGTEx_hg19 <- 
  cbind(
    as_tibble(overlaps$results_per_snp3),
    as_tibble(overlaps$significant_E_Ps_notGTEx_hg19)
  ) %>% 
  #select cols in other significant dataset to merge
  select(names(significant_SNVs_withinGTEx_hg19)) %>%
  #keep only one per OR (SNP)
  group_by(ID,enh,ensembl_gene_id, recessive) %>% slice_max(effect_size, n=1, with_ties = F) %>%  ungroup()


dim(significant_SNVs_withinGTEx_hg19)
dim(significant_unique_E_Ps_notGTEx_hg19)


#### MERGE the 2 lists
table(names(significant_SNVs_withinGTEx_hg19)==names(significant_unique_E_Ps_notGTEx_hg19))
significant_SNVs <-
  rbind (
    significant_SNVs_withinGTEx_hg19,
    #add the SNVs in enhancers, minus the ones already in E_P_eQTLs
    significant_unique_E_Ps_notGTEx_hg19 %>% dplyr::filter(!ID %in% significant_SNVs_withinGTEx_hg19$ID)
  ) 
table(duplicated(significant_SNVs$ID))


#now deal with duplicated rows due to multiple genes associated with one SNV:
significant_SNVs2<-significant_SNVs %>% 
  #select(seqnames,start,snp,enh,ensembl_gene_id, GTEx_variant_id, overlapGTEx_eQTLs, effect_size, OR_2_p, hgnc_symbol) %>% 
  #group by SNV vars (all but gene vars)
  group_by(ID, recessive) %>%
  #merge all gene ids for each SNV
  mutate(ensembl_gene_id = paste0(ensembl_gene_id, collapse = ","),
         hgnc_symbol = paste0(hgnc_symbol, collapse = ","))  %>% 
  #now keep one row per SNV
  slice_max(effect_size, with_ties = F) %>% 
  ungroup() %>% arrange(P) 
significant_SNVs2


# #what about same E-P, but separate SNPs?
# results_per_snp2 %>%  group_by(enh, ensembl_gene_id) %>% 
#   select(seqnames,start,ID,enh,ensembl_gene_id, GTEx_variant_id, overlapGTEx_eQTLs, effect_size, P) %>% 
#   filter(n()>1) %>% arrange(enh) 
# (results_per_snp2<- results_per_snp2 %>%  group_by(enh, ensembl_gene_id) %>% 
#   slice_min(P, with_ties = F) %>% arrange(P) )

head(significant_SNVs2[c("ID","enh","hgnc_symbol", "recessive")])
table(significant_SNVs2$tidytype)

# # #export UCSC track in bed format
# # names(results_per_snp3)
# # data.table::fwrite(x =
# #   results_per_snp3 %>% select ("seqnames", "end", "snp", "OR_2", "OR_2_p" ,"enh","hgnc_symbol") %>% 
# #     mutate(start=end-1, OR_scaled_p=round(scales::rescale(-log10(OR_2_p), to = c(0, 1000))), CHR=paste0("chr",seqnames), seqnames=NULL, strand=".") %>%
# #     relocate(CHR, start, end, snp, OR_scaled_p, strand, OR_2, OR_2_p, enh, hgnc_symbol)
# #   , file = "~/2HH/GWAS_summary_results/brainSNVs_22May23_hg19_UCSC.bed", sep = "\t", row.names = F, col.names = F, scipen = 999
# # )



# # import GWAS_LD block overlap data in hg19, by condition
# ##LD block work from 1000 genomes - self calculated on EUR pop

##read all block files and merge them into one variable

blocks_GR_hg19<-data.table::fread(file = GW_LD_blocks, select=c("CHR","BP1", "BP2")) %>% 
  mutate(seqnames=CHR, start=BP1, end=BP2, BP1=NULL, BP2=NULL, CHR=NULL) %>% 
  makeGRangesFromDataFrame(ignore.strand = T)
seqlevelsStyle(blocks_GR_hg19) = "UCSC"  

## Import PGC wave 3 v2 schizophrenia GWAS results (hg19) and annotate SNPs with LD blocks
#The liftOver function will create a GRangesList.
GWAS_hg19 <- data.table::fread(file = GWAS, select=c("CHR","SNP",	"POS", "OR", "P"), sep = "\t",  fill = T) %>% 
  mutate(seqnames=CHR, CHR=NULL, start=POS, POS=NULL, end=start) %>% 
  relocate(seqnames, start, end) %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
seqlevelsStyle(GWAS_hg19) = "UCSC"  


#####   MAP GWAS SNPs to LD BLOCK     ###====

###### FIND overlaps between GWAS and LD blocks

#for those overlapping a block, annotate block start and end, called start and end, and save SNP location as SNP_position
overlaps<-mergeByOverlaps(query = GWAS_hg19, subject = blocks_GR_hg19) 
GWAS_block_merge_hg19 <- 
  cbind(
    as_tibble(overlaps$blocks_GR_hg19),
    as_tibble(overlaps$GWAS_hg19) %>% mutate(SNP_position=start, start=NULL, end=NULL, seqnames=NULL, width=NULL, strand=NULL)
  )

#add SNPs not matching LD blocks
no_overlaps<-subsetByOverlaps(x = GWAS_hg19, ranges  = blocks_GR_hg19, invert = T) %>% 
  as_tibble() %>% mutate(SNP_position=start)
GWAS_block_merge_hg19 <- rbind(GWAS_block_merge_hg19,no_overlaps)

rm(GWAS_hg19)
rm(overlaps)

#remove duplications
GWAS_block_merge_hg19 <- GWAS_block_merge_hg19[!duplicated(GWAS_block_merge_hg19$SNP), ]

#check for dups
table(duplicated(GWAS_block_merge_hg19$SNP))
# GWAS_block_merge_hg19
GWAS_block_merge_hg19<-makeGRangesFromDataFrame(GWAS_block_merge_hg19,  keep.extra.columns = T)
seqlevelsStyle(GWAS_block_merge_hg19) = "UCSC" 
# sample(GWAS_block_merge_hg19)
#remove internal range duplications (overlaps). first sort by increasing p
GWAS_block_merge_hg19 <- sort(GWAS_block_merge_hg19, by = ~ P, decreasing = F)
GWAS_block_merge_hg19 <- GWAS_block_merge_hg19[unique(findOverlaps(GWAS_block_merge_hg19, type = "any", select = "first"))]
GWAS_block_merge_hg19 <- GWAS_block_merge_hg19 %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)  
seqlevelsStyle(GWAS_block_merge_hg19) = "UCSC"


#annotate GWAS SCZ overlap including significance P VALUE and EXPORT E-P_eQTL list
significant_SNVs2<-makeGRangesFromDataFrame(significant_SNVs2,  keep.extra.columns = T)
seqlevelsStyle(significant_SNVs2) = "UCSC"

# significant_SNVs2
# GWAS_block_merge_hg19

#for those overlapping a block, annotate block start and end, called start and end, and save SNP location as SNP_position
overlaps<-mergeByOverlaps(query = significant_SNVs2, subject = GWAS_block_merge_hg19) 

results_per_snp4 <- 
  cbind(
    as_tibble(overlaps$significant_SNVs2),
    as_tibble(overlaps$GWAS_block_merge_hg19) %>% 
      mutate( block_start=start, start=NULL, GWAS_SNP_position = SNP_position, SNP_position = NULL,
              block_end=end, end=NULL, seqnames=NULL, width=NULL, strand=NULL,
              GWAS_SNP=SNP, SNP=NULL, GWAS_block_top_OR=OR, OR=NULL, GWAS_block_top_p=P,
              P=NULL)
  )

table(duplicated(results_per_snp4$snp))
#add SNPs not matching LD blocks
no_overlaps<-subsetByOverlaps(x = significant_SNVs2, ranges  = GWAS_block_merge_hg19, invert = T) %>% as_tibble() %>%
  mutate(block_start=NA, GWAS_SNP_position=NA, block_end=NA,GWAS_SNP=NA,GWAS_block_top_OR=NA,
         GWAS_block_top_p=NA)
names(results_per_snp4)[names(results_per_snp4) != names(no_overlaps)]
results_per_snp4 <- rbind(results_per_snp4,no_overlaps)
results_per_snp4<-results_per_snp4 %>% arrange(P)
table(duplicated(results_per_snp4[results_per_snp4$recessive=="recessive",]$ID))



# #export SNPs of interest
#remove chr as in base file
results_per_snp4$seqnames<-gsub(pattern = 'chr', replacement = "", x = results_per_snp4$seqnames)
names(results_per_snp4)
results_per_snp4 <- results_per_snp4 %>% 
  select(seqnames,start,end,ID, tidytype,
         recessive,
         REF, ALT, A1, INFO, MAF,
         enh,GTEx_variant_id,ensembl_gene_id,
         OR, LOG.OR._SE, P,OBS_CT,
         #  OR_1,OR_1_lowCI,OR_1_high_CI,OR_1_p,
         #  OR_2, OR_2_lowCI, OR_2_high_CI,OR_2_p,FDR,
         #  OR_A_AA, OR_A_AA_lowCI, OR_A_AA_high_CI, OR_A_AA_p,
         #  OR_AA, OR_AA_lowCI, OR_AA_high_CI,OR_AA_p,
         hgnc_symbol,
         GWAS_SNP_position, block_start, block_end,     GWAS_SNP, GWAS_block_top_OR,GWAS_block_top_p
  )

table(results_per_snp4$tidytype)

#calculate OR CIs
results_per_snp4$BETA <- log(results_per_snp4$OR)
results_per_snp4$B_lower <- results_per_snp4$BETA-1.96*results_per_snp4$LOG.OR._SE
results_per_snp4$B_upper <- results_per_snp4$BETA+1.96*results_per_snp4$LOG.OR._SE
results_per_snp4$OR_CI_LOWER<-exp(results_per_snp4$B_lower)
results_per_snp4$OR_CI_UPPER<-exp(results_per_snp4$B_upper)
results_per_snp4 <- results_per_snp4 %>% arrange(ID)



#pivot wider (one row per SNP)
names(results_per_snp4)
(results_per_snp4 <- as_tibble(results_per_snp4) %>% 
  pivot_wider(names_from=recessive, 
              values_from = c(OR:OBS_CT,BETA:OR_CI_UPPER)) )

data.table::fwrite(x = results_per_snp4, file = annotated_ORs_outfilename)

#plot OR figure

plot <- results_per_snp4 %>%
  dplyr::filter((P_recessive<0.01 | P_standardGLM<0.01)) %>%
  mutate(SNP=factor(ID), SNP = fct_reorder(SNP, OR_recessive))

# str(plot)
(max_plot = max(plot$OR_CI_UPPER_recessive) + 15)

png(filename = ORfigure, width = 4000, height=2800, res = 250)
ggplot(  data = plot , aes(x=SNP ) ) + 
  scale_color_gradient(low = "navy blue", high = "white")+
  geom_pointrange(
    aes(y=OR_recessive,
        ymin = OR_CI_LOWER_recessive, 
        ymax = OR_CI_UPPER_recessive,
        col=P_recessive)) +   
  geom_pointrange( 
    color="yellow", fill="white", alpha=0.25,
    aes(
      y=OR_standardGLM,
      ymin = OR_CI_LOWER_standardGLM, 
      ymax = OR_CI_UPPER_standardGLM#,         col=OR_1_p
    )
  ) +
  ylab("log10(OR)")+ 
  xlab('')+
  scale_y_log10(limits = c(0.2, max_plot))+  
  geom_hline(yintercept =1, linetype=2)+ 
  coord_flip() +
  ggnewscale::new_scale_color() +
  geom_text( 
    aes(label=paste(enh,"gene",hgnc_symbol, "INFO", round(INFO,3), "MAF",round(MAF,2), "A1", A1), #, "chr", seqnames
        y=0, colour=factor(tidytype)), 
    size=3, hjust = 0) + scale_colour_manual(values=c("#000000", "#FF5733")) +
  geom_text( aes(label=paste("GWAS SNP: ",GWAS_SNP," OR: ",GWAS_block_top_OR,"| p ",GWAS_block_top_p), 
                 y=max_plot, hjust = "inward"), size=3) + 
  #theme(legend.position = "none") + 
  labs(title =  paste(condition, " OR in Biobank for significant E-Ps."),
       subtitle= "Coords are hg19. 
       Red text are neural facet E-Ps, black text are E-P_eQTLs.
       The blue/dark blue pointranges are the recessive allele OR with covariates, in yellow standard additive GLM model")

dev.off()








# #GO term analysis of lists
# ```{r}
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(data.table)
# library(tidyverse)


# # import brain e-QTLS and E-P pairs from our method
# #only keep gene ids
# (brain_eqtl <- fread("/mnt/biggley/home/emanuele/2HH/GTEx_eQTL/GTEx_Analysis_v8_QTLs_brain_meta_avgSlope_p0.05.tsv.gz") %>% 
#     dplyr::select(ensembl_gene_id=gene_id) %>% distinct())

# (neural_facet_any_contact <- 
#   read.csv("/mnt/biggley/home/emanuele/2HH/E_P_files/annotated_neural_facet_E_Ps_hg19_1Jun22.csv") %>% 
#     dplyr::select(ensembl_gene_id) %>% distinct())

# (E_P_eQTL <- fread("E_P_files/June22_significant_SNVs_allEPpairs_ESpos_someContact_allTotScores_anyTissue_overlapeQTLs_BRAIN_and_neuron_facet_E_Ps_hg38_uniqueEnhGenes.csv") %>% filter(type=="ESpos_someContact_overlap_eQTL") %>% 
#     dplyr::select(ensembl_gene_id) %>% distinct() )

# (neuron_facet_E_Ps_with_contact <- 
#     fread("E_P_files/June22_significant_SNVs_allEPpairs_ESpos_someContact_allTotScores_anyTissue_overlapeQTLs_BRAIN_and_neuron_facet_E_Ps_hg38_uniqueEnhGenes.csv") %>%
#     filter(neuron_facet_E_Ps_with_contact=="neuron_facet_E_Ps_with_contact")%>% 
#     dplyr::select(ensembl_gene_id) %>% distinct())


# # GO term analysis and plotting
# gene_lists <- c("brain_eqtl", "neural_facet_any_contact", "E_P_eQTL", "neuron_facet_E_Ps_with_contact")
# gene_list_N <- length(gene_lists)
  
# Enrich_go_function <- function (gene_list) {
#   enrichGO(gene = gene_list,
#                 OrgDb = "org.Hs.eg.db",
#                 keyType = "ENSEMBL",
#                 ont = "ALL",
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)
# }


# i=1
# for (gene_list in gene_lists) { 
#   print (gene_list ) 
#   ego <- Enrich_go_function (unlist(get(gene_list)))
#   (nam <- paste("fig", i, sep = ""))
#   print(nam)
#   assign(nam, 
#          dotplot(ego, showCategory=25, title=paste("GO term analysis",gene_list,"genes, N=",nrow(get(gene_list))))
#   )
#   i<-i+1
  
#   }


# png("figs/GO_plot_17Jun22.png", width = 25, height = 15, units = 'in', res = 300)
# fig1 + fig2 + fig3 + fig4
# dev.off()
# ```
