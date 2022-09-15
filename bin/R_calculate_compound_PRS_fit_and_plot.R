#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
require(rms)

#INPUT
args = commandArgs()


print(args)


(condition = args[8])
diagnosis = fread(args[9], select=c(1,3))
colnames(diagnosis) = c("FID", "dx")
TS_EPs_summary = args[10]
# TS_EPs_prsice = args[11]
TS_EPs_best_individual_PRS_scores = 
  fread(args[12]) %>% dplyr::rename(PRS_TS_EPs = PRS)
original_NoEPsOverlap_summary = args[13]
# original_NoEPsOverlap_prsice = args[14]
original_NoEPsOverlap_best_individual_PRS_scores = fread(args[15])  %>% 
  dplyr::rename(PRS_original_noOverlap = PRS)
original_summary = args[16]
original_best = fread(args[17]) %>% dplyr::rename(PRS_original = PRS)



#OUTPUT
analysis_output_txt = paste0(condition, "_", Sys.Date(),"_model_fit_measures.txt")
PRS_QUANTILES_PLOT  = paste0(condition, "_", Sys.Date(),"_PRS_QUANTILES_PLOT.png")


#import thresholds
(TS_EPs_threshold = as.numeric(
  fread(TS_EPs_summary) %>% dplyr::select(Threshold)
  ))
(original_no_over_threshold = as.numeric(
  fread(original_NoEPsOverlap_summary) %>% dplyr::select(Threshold)
  ))
(original_threshold = as.numeric(
  fread(original_summary) %>% dplyr::select(Threshold)
  ))


#create total PRS score
(total_PRS <- original_best %>%
    left_join(original_NoEPsOverlap_best_individual_PRS_scores) %>% 
    left_join(TS_EPs_best_individual_PRS_scores) %>% 
    left_join(diagnosis) %>% 
    mutate(dx=factor(dx), IID=factor(IID)) %>% 
    select(-FID, -In_Regression) %>% 
    remove_missing())
# hist(total_PRS$PRS_TS_EPs)
scaled_total_PRS <- total_PRS
scaled_total_PRS[,c(2,3,4)] <-  as_tibble(scale(total_PRS[,c(2,3,4)]))+10
    
(scaled_total_PRS<-scaled_total_PRS %>% 
  mutate(total_PRS = PRS_TS_EPs + PRS_original_noOverlap) %>% 
    remove_missing() %>% 
    #generate quantiles
    mutate(base_quantile = factor(ntile(PRS_original, n = 10))) %>% 
    mutate(base_noover_quantile = factor(ntile(PRS_original_noOverlap, n = 10))) %>% 
    mutate(PRS_TS_EPs_quantile = factor(ntile(PRS_TS_EPs, 10))) %>% 
    mutate(total_PRS_quantile = factor(ntile(total_PRS, 10))) )

# str(scaled_total_PRS)



#measures of model fit 
# Start writing to an output file
sink(analysis_output_txt)
(mod0 <- lrm(dx ~ PRS_original, data = scaled_total_PRS))
(mod <- lrm(dx ~ PRS_original_noOverlap, data = scaled_total_PRS))
(mod1 <- lrm(dx ~ PRS_original_noOverlap + PRS_TS_EPs, data = scaled_total_PRS))
(mod2 <- lrm(dx ~ PRS_original_noOverlap * PRS_TS_EPs, data = scaled_total_PRS))
(mod3 <- lrm(dx ~ total_PRS, data = scaled_total_PRS))
# Stop writing to the file
sink()

##plot ORs
# base_ORs
summary(logistic<-glm(formula = dx ~ base_quantile, data = scaled_total_PRS, family = binomial, na.action = "na.omit"))
(base_ORs<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
colnames(base_ORs) <- c("quantile", "OR", "LCI", "UCI")
base_ORs$type=paste0("PRS_base_thresh_",original_threshold)

#base_no overlap
summary(logistic<-glm(formula = dx ~ base_noover_quantile, data = scaled_total_PRS, family = binomial, na.action = "na.omit"))
(base_noover_ORs<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
colnames(base_noover_ORs) <- c("quantile", "OR", "LCI", "UCI")
base_noover_ORs$type=paste0("PRS_base_no_overlap_thresh_",original_no_over_threshold)


summary(logistic<-glm(formula = dx ~ total_PRS_quantile, data = scaled_total_PRS, family = binomial, na.action = "na.omit"))
(total_PRS_ORs<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
colnames(total_PRS_ORs) <- c("quantile", "OR", "LCI", "UCI")
total_PRS_ORs$type=paste0("PRS_base_plus_TS_EPs_thresh_",TS_EPs_threshold)
total_PRS_ORs

(tot_ORs <- 
    rbind(base_ORs, base_noover_ORs, total_PRS_ORs) %>% 
    separate(type, c("type", "threshold"), sep = "_thresh_") %>% 
    filter(quantile!="(Intercept)") %>% 
    mutate(quantile=gsub(perl = T, pattern = "[a-z,A-Z,_]", replacement="", x=quantile)) %>% 
    mutate(quantile=as.numeric(quantile)) %>% 
    pivot_wider(names_from = type, values_from = c(OR,LCI,UCI, threshold ), id_cols = quantile))
tot_ORs<-rbind(
  tot_ORs,
  rep(x = 1, times=ncol(tot_ORs))
)

png(filename = PRS_QUANTILES_PLOT, 
    width = 1500, height=2000, res = 200)
ggplot(  data = tot_ORs , aes(x=reorder(quantile, (quantile)))) + 
  scale_color_gradient(low = "navy blue", high = "white")+
  geom_pointrange(
    color="black", alpha=0.3,
    aes(y=OR_PRS_base,
        ymin = LCI_PRS_base, 
        ymax = UCI_PRS_base
    ), position=position_dodge(width=1)) +   
  geom_pointrange(
    color="green", alpha=0.5,
    aes(y=OR_PRS_base_no_overlap,
        ymin = LCI_PRS_base_no_overlap, 
        ymax = UCI_PRS_base_no_overlap#,        col=OR_2_p
    ), position=position_dodge(width=1)) +   
  geom_pointrange( 
    color="red", alpha=0.5,
    aes(
      y=OR_PRS_base_plus_TS_EPs,
      ymin = LCI_PRS_base_plus_TS_EPs, 
      ymax = UCI_PRS_base_plus_TS_EPs#,         col=OR_1_p
    ), position=position_dodge(width=-1)
  ) +
  ylab("OR for dx by PRS decile")+ 
  xlab('PRS decile')+ylim(0.5,NA)+scale_y_continuous(breaks = seq(0, 100, by = 1))+
  #scale_y_log10()+ 
  geom_hline(yintercept = 1, linetype=2)+
  #theme(legend.position = "none") + 
  labs(title =  "dx OR in Biobank for PRS decile",
       subtitle= paste("In black, base, thresh = ",tot_ORs$threshold_PRS_base[1],
       "in green, non-overlapping base PRS (PGC), C+T threshold=", tot_ORs$threshold_PRS_base_no_overlap[1],
       "In red, base PRS + TS_EPs PRS, C+T p=1.",tot_ORs$threshold_PRS_base_plus_TS_EPs[1]) )
dev.off()  





# # PLOT PPV 
# scaled_total_PRS<-scaled_total_PRS %>% 
#   #generate PRS cutoff 
#   group_by(base_quantile) %>% mutate(base_quantile_cutoff=min(PRS)) %>% #dplyr::slice_min(order_by =original_PRS,n =  1) %>%
#   ungroup() %>% 
#   #generate PRS cutoff 
#   group_by(total_PRS_quantile) %>% mutate(total_PRS_quantile_cutoff=min(total_PRS)) %>% #dplyr::slice_min(order_by =original_PRS,n =  1) %>%
#   ungroup() 
# table(scaled_total_PRS$base_quantile_cutoff)
# table(scaled_total_PRS$total_PRS_quantile_cutoff)
# (scaled_total_PRS_2<-scaled_total_PRS %>% 
#   mutate(
#     qm10=ifelse(test = PRS<8.7203522242309, yes = 1, no = 0),
#     qm20=ifelse(test = PRS<9.15800755903283, 1, 0),
#     qm50=ifelse(test = PRS<9.99869005614688, 1, 0),
#     qm60=ifelse(test = PRS<10.2510892227546, 1, 0),
#     qp90=ifelse(test = PRS>11.2814484781398, 1, 0)
#   ) %>% 
#     mutate(
#       qm10_tot=ifelse(test = total_PRS<18.2962887232296, yes = 1, no = 0),
#       qm20_tot=ifelse(test = total_PRS<18.9245936262333, 1, 0),
#       qm50_tot=ifelse(test = total_PRS<20.0343609910009, 1, 0),
#       qm60_tot=ifelse(test = total_PRS<20.3608376468006, 1, 0),
#       qp90_tot=ifelse(test = total_PRS>21.7174666341562, 1, 0)
#     ) 
# )
# sample_n(scaled_total_PRS_2[c(1,2,7,15:22)], size = 10)
# table(scaled_total_PRS_2$qp90_tot, scaled_total_PRS_2$total_PRS_quantile)


# #calculate PPVs
# (qm10<-scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#     select(qm10,dx) %>% table() %>% as.matrix())
# qm20<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm20,dx) %>% table() %>% as.matrix())
# qm50<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm50,dx) %>% table() %>% as.matrix())
# obs<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#         select(dx) %>% table() %>% as.matrix())
# qm60<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm60,dx) %>% table() %>% as.matrix())
# qp90<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qp90,dx) %>% table() %>% as.matrix())
# (PPVs_original<-rbind(
#   qm10=qm10[2,2]*100/(qm10[2,2]+qm10[2,1]),
#   qm20=qm20[2,2]*100/(qm20[2,2]+qm20[2,1]),
#   qm50=qm50[2,2]*100/(qm50[2,2]+qm50[2,1]),
#   obs=obs[2,1]*100/sum(obs),
#   qm60=qm60[2,2]*100/(qm60[2,2]+qm60[2,1]),
#   qp90=qp90[2,2]*100/(qp90[2,2]+qp90[2,1])) %>% as_tibble(rownames = "quant") %>% 
#   mutate(type="original_nooverlap"))
# colnames(PPVs_original) = c("threshold", "PPV_at_threshold", "type")
# PPVs_original

# (qm10<-scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#     select(qm10_tot,dx) %>% table() %>% as.matrix())
# qm20<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm20_tot,dx) %>% table() %>% as.matrix())
# qm50<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm50_tot,dx) %>% table() %>% as.matrix())
# obs<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#         select(dx) %>% table() %>% as.matrix())
# qm60<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qm60_tot,dx) %>% table() %>% as.matrix())
# qp90<-(scaled_total_PRS_2 %>% #dplyr::filter(rs34380086_caco=="0") %>% 
#          select(qp90_tot,dx) %>% table() %>% as.matrix())
# (PPVs_tot<-rbind(
#   qm10=qm10[2,2]*100/(qm10[2,2]+qm10[2,1]),
#   qm20=qm20[2,2]*100/(qm20[2,2]+qm20[2,1]),
#   qm50=qm50[2,2]*100/(qm50[2,2]+qm50[2,1]),
#   obs=obs[2,1]*100/sum(obs),
#   qm60=qm60[2,2]*100/(qm60[2,2]+qm60[2,1]),
#   qp90=qp90[2,2]*100/(qp90[2,2]+qp90[2,1])
# ) %>% as_tibble(rownames = "quant") %>% 
#   mutate(type="original_plus_TS_EPs"))
# colnames(PPVs_tot) = c("threshold", "PPV_at_threshold", "type")


# (plot <- rbind(PPVs_original,PPVs_tot) %>% 
#     pivot_wider(names_from = c("type"), values_from = PPV_at_threshold))

# #plot similar to Davies et al, PPV over PRS quantiles
# png(filename = paste0("figs/",Sys.Date(),
#                       "_PPV_dx_by_PRS_cutoff_original_vs_partitionedPRS.png"), 
#     width = 2000, height=2000, res = 300)
# ggplot(data = plot, aes(x=factor(threshold))) +
#   geom_dotplot(aes(y=original_nooverlap), binaxis='y', position = "dodge")  +
#   geom_dotplot(aes(y=original_plus_TS_EPs), binaxis='y', position = "dodge", 
#                stackdir='center', fill="red", colour="red")  +
#   ggthemes::scale_fill_colorblind()+ylim(0, NA)+
#   ggtitle("PPVs (y axis) for dx based on various cutoffs of dx_PRS.", 
#            subtitle = "In black, PPV for the original PRS, in red partitioned PRS")+
#   xlab(label = "Subgroups by different percentile PRS cutoff")+
#   ylab(label = "PPV")+scale_y_continuous(breaks = seq(0, 100, by = 0.05))
# dev.off()
