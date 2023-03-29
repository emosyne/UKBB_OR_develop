#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

# R_annotate_ORs.R ${condition} ${PLINK_ORs_recessive} ${PLINK_ORs_dominant} ${PLINK_ORs_additive} ${full_GWAS} ${SNP_frq} 

(condition = args[8])
PLINK_ORs_recessive = args[9]
PLINK_ORs_dominant = args[10]
PLINK_ORs_additive = args[11]
full_GWAS = args[12]
SNP_frq = args[13]


(results_per_snp <- rbind(
  data.table::fread(PLINK_ORs_recessive, na.strings = ".", header = T, col.names = c("CHR", "start", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "OBS_CT", "OR", "Beta_SE", "L95", "U95", "Z_STAT", "P", "ERRCODE" )),
  data.table::fread(PLINK_ORs_additive, na.strings = ".", header = T, col.names = c("CHR", "start", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "OBS_CT", "OR", "Beta_SE", "L95", "U95", "Z_STAT", "P", "ERRCODE" )),
  data.table::fread(PLINK_ORs_dominant, na.strings = ".", header = T, col.names = c("CHR", "start", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "OBS_CT", "OR", "Beta_SE", "L95", "U95", "Z_STAT", "P", "ERRCODE" ))
  ))

#OUTPUT
annotated_ORs_outfilename = paste0(condition, "_", PLINK_ORs_recessive, "_", Sys.Date(), "_annotated_ORs.csv")
dir.create(file.path(".", "figs"))
ORfigure = paste0("figs/PLINK_ORs_", condition, "_", PLINK_ORs_recessive, "_", Sys.Date(),"_UKBB.png")

## Annotate with  INFO scores, MAF


(results_per_snp2 <- results_per_snp %>% dplyr::select(-ALT,-FIRTH,-Z_STAT,-OBS_CT) %>% 
    left_join(frq) %>% group_by(SNP, TEST) %>% slice_head(n=1) %>% ungroup())





# filter unique SNVs and then unique E-P_eQTLs by joining BB ORs back to E-P and eQTL SNV data 


(results_per_snp3 <- as_tibble(results_per_snp2, .name_repair = "universal") %>% 
  #remove error SNPs and low MAF
  dplyr::filter(MAF>0.01, is.na(ERRCODE)) %>% #
  dplyr::mutate(seqnames=CHR, end=start, CHR=NULL, ERRCODE=NULL, NCHROBS=NULL, Beta_SE=NULL) )
  


(results_per_snp_GR = results_per_snp3  %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T))
seqlevelsStyle(results_per_snp_GR) = "UCSC"
results_per_snp_GR


## now annotate with GWAS
(gwas = data.table::fread(GWAS) %>% 
    #make fields equivalent to results
    dplyr::select(seqnames=CHR, start=POS, end=POS, SNP, A1, A2, OR, P) 
    # dplyr::mutate(TEST="PGC", L95=NA, U95=NA, MAF=NA )
    )
(gwas_GR = makeGRangesFromDataFrame(gwas, keep.extra.columns = T))
seqlevelsStyle(gwas_GR) = "UCSC"
gwas_GR

(hits = unlist(nearest(results_per_snp_GR, gwas_GR, select="arbitrary")))

(edited_GWAS = gwas[hits,] %>% dplyr::select(closest_PGC_OR=OR,closest_PGC_P=P,closest_PGC_SNP=SNP))

(merge = as_tibble(cbind(results_per_snp3,edited_GWAS)))



#plot OR figure

(plot=merge %>%
  dplyr::select(SNP, A1, TEST,OR, L95,U95,P, closest_PGC_OR, closest_PGC_P, closest_PGC_SNP) %>% 
  pivot_wider(id_cols = c(SNP, A1,closest_PGC_OR, closest_PGC_P, closest_PGC_SNP), names_from = TEST, values_from = c(OR,L95,U95,P)) %>% 
  dplyr::filter((P_DOM<0.005 )) %>%#| P_REC<0.001
  dplyr::mutate(SNP=factor(SNP), SNP = fct_reorder(SNP, OR_DOM, .na_rm = T))
)
plot[grep(x = plot$SNP, pattern = "rs3438"),]

colours = values=MetBrewer::met.brewer("Austria", 4)

(add_dominant = ggplot(  data = plot , aes(x=SNP ) ) + 
    scale_color_gradientn(colours = c(colours[4],"gray","gray"),
                          values = scales::rescale(c(0,0.005,1)), breaks = c(0,1),
                          guide = "colorbar", limits=c(0,1),name = "BH-adjusted p value",)+
    geom_point( 
      alpha=1,
      aes(
        y=closest_PGC_OR,
        col=closest_PGC_P
      ))+
    geom_pointrange(
      color=colours[1], fill="white", alpha=0.5,
      aes(
        y=OR_ADD,
        ymin = L95_ADD,
        ymax = U95_ADD#,         col=OR_1_p
      )
    ) +
    geom_pointrange( 
      color=colours[2], fill="white", alpha=0.6,
      aes(
        y=OR_DOM,
        ymin = L95_DOM, 
        ymax = U95_DOM#,         col=OR_1_p
      )
    ) +
    ylab("log10(OR)")+ 
    xlab('')+
    scale_y_log10()+  
    geom_hline(yintercept =1, linetype=2)+ 
    coord_flip() +
    theme_bw()+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 6), 
          title = element_text(size = 8)) + 
    ylab("log-transformed OR (p<0.005 for dominant) for schizophrenia")+
    labs(title =  paste0(condition, " ORs for enhancer-based SNPs."),
         subtitle = "Blue: dominant model. Red: additive model.\nSmall dots: closest PGC SNP - yellow if p<0.005"))

ggsave(
  filename = paste0("figs/PLINK_ORs_", condition, "_add_dom_", Sys.Date(),"_UKBB.pdf"), 
  add_dominant,  
  width = 5, height = 6, dpi = "print", device = "pdf", scale = 1.2)



(plot=merge %>%
    dplyr::select(SNP, A1, TEST,OR, L95,U95,P, closest_PGC_OR, closest_PGC_P, closest_PGC_SNP) %>% 
    pivot_wider(id_cols = c(SNP, A1,closest_PGC_OR, closest_PGC_P, closest_PGC_SNP), names_from = TEST, values_from = c(OR,L95,U95,P)) %>% 
    dplyr::filter((P_REC<0.005 )) %>%#| P_REC<0.001
    dplyr::mutate(SNP=factor(SNP), SNP = fct_reorder(SNP, OR_REC, .na_rm = T))
)
plot[grep(x = plot$SNP, pattern = "rs3438"),]
(recessive = ggplot(  data = plot , aes(x=SNP ) ) + 
    scale_color_gradientn(colours = c(colours[4],"gray","gray"),
                         values = scales::rescale(c(0,0.005,1)), breaks = c(0,1),
                         guide = "colorbar", limits=c(0,1),name = "BH-adjusted p value",)+
    geom_point( 
      # color="black", alpha=0.5,
      aes(
        y=closest_PGC_OR,
        col=closest_PGC_P
      ))+
    geom_pointrange(
      color=colours[3], fill="white", alpha=0.7,
      aes(y=OR_REC,
          ymin = L95_REC,
          ymax = U95_REC#,col=P_REC
          )) +
    geom_pointrange(
      color=colours[1], fill="white", alpha=0.5,
      aes(
        y=OR_ADD,
        ymin = L95_ADD,
        ymax = U95_ADD#,         col=OR_1_p
      )
    ) +
    # geom_pointrange( 
    #   color=colours[2], fill="white", alpha=0.5,
    #   aes(
    #     y=OR_DOM,
    #     ymin = L95_DOM, 
    #     ymax = U95_DOM#,         col=OR_1_p
    #   )
    # ) +
    ylab("log10(OR)")+ 
    xlab('')+
    scale_y_log10()+  
    geom_hline(yintercept =1, linetype=2)+ 
    coord_flip() +
    theme_bw()+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 6), 
          title = element_text(size = 8)) + 
    ylab("log-transformed OR (p<0.005 for recessive) for schizophrenia")+
    labs(title =  paste0(condition, " ORs for enhancer-based SNPs."),
         subtitle = "Green: recessive model. Red: additive model.\nSmall dots: closest PGC SNP - yellow if p<0.005"))


ggsave(
  filename = paste0("figs/PLINK_ORs_", condition, "_add_REC_", Sys.Date(),"_UKBB.pdf"), 
  recessive,  
  width = 5, height = 6, dpi = "print", device = "pdf", scale = 1.2)


#produce bed file for pipeline
#make bed
results_per_snp3 %>% select(seqnames,  start,end, SNP) %>%
  mutate(score=".", strand=".") %>%
  relocate(seqnames, start, end, SNP, score, strand) %>%
  fwrite(file="EP_WAS.bed", sep="\t", col.names=F)

# produce file in HCM GWAS format
#SNP     A1      A2      Z       N       FRQ     P       POS     CHR     BETA    SE
data.table::fread(PLINK_ORs_recessive, na.strings = ".", header = T, col.names = c("CHR", "POS", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "N", "OR", "SE", "L95", "U95", "Z", "P", "ERRCODE" )) %>% 
  dplyr::mutate(BETA=log(OR), OR=NULL) %>% 
  dplyr::select(SNP,A1,A2,Z,N, P,POS,CHR,BETA,SE) %>% 
  data.table::fwrite(sep = "\t", file = "UKBB_ENH_associations_REC.tsv.gz")
  
data.table::fread(PLINK_ORs_dominant, na.strings = ".", header = T, col.names = c("CHR", "POS", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "N", "OR", "SE", "L95", "U95", "Z", "P", "ERRCODE" )) %>% 
  dplyr::mutate(BETA=log(OR), OR=NULL) %>% 
  dplyr::select(SNP,A1,A2,Z,N, P,POS,CHR,BETA,SE) %>% 
  data.table::fwrite(sep = "\t", file = "UKBB_ENH_associations_DOM.tsv.gz")

data.table::fread(PLINK_ORs_additive, na.strings = ".", header = T, col.names = c("CHR", "POS", "SNP", "A2", "ALT", "A1", "FIRTH", "TEST", "N", "OR", "SE", "L95", "U95", "Z", "P", "ERRCODE" )) %>% 
  dplyr::mutate(BETA=log(OR), OR=NULL) %>% 
  dplyr::select(SNP,A1,A2,Z,N, P,POS,CHR,BETA,SE) %>% 
  data.table::fwrite(sep = "\t", file = "UKBB_ENH_associations_ADD.tsv.gz")