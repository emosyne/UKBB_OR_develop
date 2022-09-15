#!/usr/bin/env Rscript

library(tidyverse)

#INPUT
args = commandArgs()

print(args)
# ${meta} ${original_PRS_prsice} ${original_NoEPsOverlap_prsice} ${TS_partition_PRS_prsice} ${TS_EPs_prsice} 
(condition = args[8])
original_PRS_prsice = args[9]
original_NoEPsOverlap_prsice = args[10]
TS_partition_PRS_prsice = args[11]
TS_EPs_prsice = args[12]



#OUTPUT
PRS_comparison_figure_path = paste0(condition, "_", Sys.Date(), "_PRS_comparison_figure.png")


# Compare r2 BETWEEN PRSs:

original<- cbind(
  data.table::fread(original_PRS_prsice),
  dataset="Original_GWAS_PRS_additive"
  )


original_NoEPsOverlap_prsice<- cbind(
  data.table::fread(original_NoEPsOverlap_prsice),
  dataset="Original_NoEPsOverlap_clumped"
  )

TS <- cbind(
  data.table::fread(TS_partition_PRS_prsice),
  dataset="PRS_TS_partition_only_recessive"
  )

TS_EPs_no_overlap_clumped<- cbind(
  data.table::fread(TS_EPs_prsice),
  dataset="TS_EPs_no_overlap_clumped"
  )

#include both bests in thresholds
all_PRS <- rbind(original,
                original_NoEPsOverlap_prsice,
                 TS, TS_EPs_no_overlap_clumped) %>% 
  filter(Threshold %in% c(5e-08,0.00010005, 0.01, 0.0218359, 0.0267362, 0.0555001, 0.0697, 0.1021 ,0.2024, 0.3003, 0.4981, 1)) 

# table(all_PRS$dataset)

# library(MetBrewer)
dodge <- position_dodge(1)
plot
png(filename = PRS_comparison_figure_path,
    width = 3000, height=2000, res = 300)
ggplot(all_PRS, aes(x=factor(Threshold), y=R2.adj, label=paste0("Adj R2=",round(R2.adj,3),",\n N=",Num_SNP), 
                     fill=factor(dataset))) +
  geom_dotplot(binaxis='y', position = "dodge", stackdir='center') +  scale_colour_brewer(palette = "Set1")+ #scale_fill_manual(values=met.brewer("Lakota", 7))+
  ggtitle("Adjusted pseudo-R2 calculated by PRSice at several thresholds", 
  subtitle = paste("Original GWAS PRS, vs modified PRSs for", condition) )+
  xlab(label = "PRSice p-value threshold")+labs(fill='PRS') + 
  ggrepel::geom_text_repel(all_PRS, mapping=aes(segment.color="red"),position = dodge,min.segment.length = 0)+ 
  theme(legend.position="bottom")

dev.off()




# # PLOT PPV BY SNP CASE/CONTROL

# ```{r plot as case-control}
# library(tidyverse)
# #extract genotypes at SNP
# geno_at_snp <- BEDMatrix::BEDMatrix("~/2HH/2HHworkfiles/2022_06_10_BRAIN_SNVs_includingNeuronFacet_hg19_allchr_brainSNVsonly.bed", simple_names = T) %>% 
#   as.matrix() %>% as_tibble(rownames = "FID") %>%  dplyr::select(FID,rs34380086)
# head(geno_at_snp)

# #extract PRS per id
# individual_PRS<-   data.table::fread("~/2HH/2HHworkfiles/PRSice_output/2022_06_01_ORIG_BRAIN_PGC3_SCZ_wave3_public.v2.GRCh37_noline5328757_clump_p1_1_500kb_r20.1_NOclump_keepAmbig.best", select = c("FID", "PRS"), colClasses = c(FID="character",PRS="numeric")) %>% mutate(original_PRS=PRS, PRS=NULL)
# individual_PRS

# #import scz diagnosis
# scz<-data.table::fread("~/2HH/biobank/genetics/scz_newDec21.pheno", select=c("FID", "scz"), colClasses = c(FID="character",scz="numeric"))

# #join all
# individual_PRS_geno_rs34380086<- right_join(geno_at_snp, individual_PRS, by=c("FID")) %>% left_join(scz, by=c("FID")) %>% 
#   remove_missing() %>% mutate(rs34380086_caco=factor(ifelse(rs34380086==2,1,0)), rs34380086=NULL) %>% 
#   mutate(scz=factor(ifelse(scz==2,1,0))) %>% 
#   ##mean PRS by genotype
#   #group_by(rs34380086_caco) %>%    mutate(PRS_mean_by_rs34380086_caco=mean(original_PRS)) %>% ungroup()%>% 
#   #generate quantiles
#   mutate(quantile = factor(ntile(original_PRS, 10))) %>% 
#   #generate PRS cutoff 
#   group_by(quantile) %>% mutate(cutoff=min(original_PRS)) %>% #dplyr::slice_min(order_by =original_PRS,n =  1) %>%
#   ungroup() 
# table(individual_PRS_geno_rs34380086$quantile, individual_PRS_geno_rs34380086$cutoff)
# individual_PRS_geno_rs34380086<-individual_PRS_geno_rs34380086 %>% 
#   mutate(
#     qm10=ifelse(test = original_PRS>-0.0126497359, yes = 1, no = 0),
#     qm20=ifelse(test = original_PRS>-0.0125646899, 1, 0),
#     qm50=ifelse(test = original_PRS>-0.0125363065, 1, 0),
#     qp50=ifelse(test = original_PRS>-0.012529286, 1, 0),
#     qp90=ifelse(test = original_PRS>-0.0124936337, 1, 0)
#   ) %>% arrange(-original_PRS)
# sample_n(individual_PRS_geno_rs34380086, size = 10)
# table(individual_PRS_geno_rs34380086$qm10, individual_PRS_geno_rs34380086$quantile)

# gmodels::CrossTable(individual_PRS_geno_rs34380086$scz)
# gmodels::CrossTable(individual_PRS_geno_rs34380086$rs34380086_caco,individual_PRS_geno_rs34380086$scz, prop.r = T, prop.c = T, prop.chisq = F, prop.t = F, expected = T)


# #calculate PPVs
# (qm10<-individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm10,scz) %>% table() %>% as.matrix())
# qm20<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm20,scz) %>% table() %>% as.matrix())
# qm50<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm50,scz) %>% table() %>% as.matrix())
# obs<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#         select(scz) %>% table() %>% as.matrix())
# qp50<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qp50,scz) %>% table() %>% as.matrix())
# qp90<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qp90,scz) %>% table() %>% as.matrix())
# PPVs_rs34380086neg<-rbind(
#     qm10=qm10[2,2]*100/(qm10[2,2]+qm10[2,1]),
#     qm20=qm20[2,2]*100/(qm20[2,2]+qm20[2,1]),
#     qm50=qm50[2,2]*100/(qm50[2,2]+qm50[2,1]),
#     obs=obs[2,1]*100/sum(obs),
#     qp50=qp50[2,2]*100/(qp50[2,2]+qp50[2,1]),
#     qp90=qp90[2,2]*100/(qp90[2,2]+qp90[2,1])
# ) %>% as_tibble(rownames = NA) %>%  dplyr::rename("PPVs_rs34380086neg"=`1`)
# qm10_SNP<-individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#          select(qm10,scz) %>% table() %>% as.matrix()
# qm20_SNP<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#          select(qm20,scz) %>% table() %>% as.matrix())
# qm50_SNP<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#          select(qm50,scz) %>% table() %>% as.matrix())
# (obs_SNP<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#         select(scz) %>% table() %>% as.matrix()))
# (qp50_SNP<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#          select(qp50,scz) %>% table() %>% as.matrix()))
# qp90_SNP<-(individual_PRS_geno_rs34380086 %>% dplyr::filter(rs34380086_caco=="1") %>% 
#          select(qp90,scz) %>% table() %>% as.matrix())
# PPVs_rs34380086pos<-rbind(
#     qm10=qm10_SNP[1,2]*100/(qm10_SNP[1,1]+qm10_SNP[1,2]),
#     qm20=qm20_SNP[2,2]*100/(qm20_SNP[2,2]+qm20_SNP[2,1]),
#     qm50=qm50_SNP[2,2]*100/(qm50_SNP[2,2]+qm50_SNP[2,1]),
#     obs=obs_SNP[2,1]*100/sum(obs_SNP),
#     qp50=qp50_SNP[2,2]*100/(qp50_SNP[2,2]+qp50_SNP[2,1]),
#     qp90=qp90_SNP[2,2]*100/(qp90_SNP[2,2]+qp90_SNP[2,1])
# ) %>% as_tibble(rownames = NA) %>% dplyr::rename("PPVs_rs34380086pos"="1")


# #plot similar to Davies et al, PPV over PRS quantiles
# png(filename = paste0("figs/",Sys.Date(),
#                       "_PPV_SCZ_by_PRS_cutoff_Davies_et_al.png"), 
#     width = 3500, height=2000, res = 300)
# cbind(PPVs_rs34380086neg,PPVs_rs34380086pos) %>% rownames_to_column(var="quantile") %>% pivot_longer(cols = c("PPVs_rs34380086neg", "PPVs_rs34380086pos")) %>% mutate(value=round(value, 2)) %>% 
#   ggplot(aes(x=factor(quantile), y=value, label=value, fill=factor(name))) +
#   geom_dotplot(binaxis='y', position = "dodge", stackdir='center') +  
#   ggthemes::scale_fill_colorblind()+
#   ggtitle("PPVs (y axis) for SCZ based on various cutoffs of SCZ_PRS.", 
#           subtitle = "Colours differentiate values from the PPVs_rs34380086pos cohort (orange) versus values estimated from the PPVs_rs34380086 negatives (black)")+
#   xlab(label = "Subgroups by different percentile PRS cutoff")+
#   labs(fill='rs34380086 status') + 
#   ggrepel::geom_text_repel(mapping=aes(segment.color="red"))
# dev.off()

# #write plink pheno file
# data.table::fwrite(
#   data.frame(FID = individual_PRS_geno_rs34380086[individual_PRS_geno_rs34380086$rs34380086==2,]$FID, IID = individual_PRS_geno_rs34380086[individual_PRS_geno_rs34380086$rs34380086==2,]$FID),
#   sep = "\t", file = "~/2HH/biobank/genetics/rs34380086_cases.pheno")


# #plot as density plots
# # Use semi-transparent fill
# ggplot(individual_PRS_geno_rs34380086, aes(x=original_PRS, fill=rs34380086_caco)) +
#   geom_density(alpha=0.2)

# #plot as quantile plots
# individual_PRS_geno_rs34380086  %>% 
#   ggplot(  aes(x=reorder(snp, (OR_2)))) + 
#   geom_boxplot(    aes(y=original_PRS,x=quantile, colour=rs34380086_caco))
         
         

# ```


