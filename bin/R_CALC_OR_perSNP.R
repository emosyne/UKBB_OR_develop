#!/usr/bin/env Rscript

## PART ONE: IMPORT SNV AND GWAS FILES, MAKE BED FILES
library(tidyverse)
# library(GenomicRanges)
# library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

(extracted_allchrom_genotypes_SNVs_hg19 = args[8])
(phenofile = args[9])

#OUTPUT
(ORs.name = paste0(sub('\\.pheno', '', phenofile), "_odds_ratios_UKBB.csv"))


#IMPORT SNV GENOTYPES ONLY
m1 <- BEDMatrix::BEDMatrix(extracted_allchrom_genotypes_SNVs_hg19, simple_names = T)
m1<-as.matrix(m1) 

sparsem1<-Matrix::Matrix(m1, sparse = T)
sparsem1[1:8,1:10]
dim(sparsem1)

rm(m1)
# gc()

#remove participants with names starting with minus (withdrawn)
#(sparsem1<-sparsem1[-grepl(pattern = "^-", x = rownames(sparsem1),  perl=TRUE), ])
sparsem1<-sparsem1[,c("rs34380086")]
#sparsem1[,1:2]

dim(sparsem1)


# Calculate OR by genotype and make quantile plots ====
#merge burden info with diagnosis and calculate OR
(dx <- data.table::fread(phenofile,
                        select = c(2,3), colClasses = c("IID"="character"), col.names=c("IID", "dx")) )


#per snp
results_per_snp <- data.frame(
  snp=NA_character_,
  OR_1=NA_integer_,
  OR_1_lowCI=NA_integer_,
  OR_1_high_CI=NA_integer_,
  OR_1_p=NA_integer_,
  OR_2=NA_integer_,
  OR_2_lowCI=NA_integer_,
  OR_2_high_CI=NA_integer_,
  OR_2_p=NA_integer_,
  OR_A_AA=NA_integer_,
  OR_A_AA_lowCI=NA_integer_,
  OR_A_AA_high_CI=NA_integer_,
  OR_A_AA_p=NA_integer_,
  OR_AA=NA_integer_,
  OR_AA_lowCI=NA_integer_,
  OR_AA_high_CI=NA_integer_,
  OR_AA_p=NA_integer_
)
#for(i in 1: length(colnames(sparsem1))){
for(i in 1: 1){
  print (paste(i," out of ", length(colnames(sparsem1))))
  
  snp <- "rs34380086"
  
  print(snp)
  SNP_perIID_genotypes<-as_tibble(as.matrix(sparsem1)[,i], rownames="IID") %>% left_join(dx) %>% 
    mutate(genotype=factor(value,levels = c(0,1,2)), dx = factor(dx, levels = c(1,2), labels = c("no dx", "dx")), value=NULL)
  
  #glm
  SNP_perIID_genotypes$genotype_A_AA<-0
  SNP_perIID_genotypes$genotype_AA<-0
  SNP_perIID_genotypes$genotype_A_AA<-ifelse(SNP_perIID_genotypes$genotype==1 | SNP_perIID_genotypes$genotype==2, yes = 1, 0)
  SNP_perIID_genotypes$genotype_AA<-ifelse(SNP_perIID_genotypes$genotype==2, yes = 1, 0)
  sample_n(SNP_perIID_genotypes, size=10)
  
  #0,1,2
  summary(m12<-glm(formula = dx~genotype, family = binomial, data = SNP_perIID_genotypes))
  #A, AA
  summary(mAAA<-glm(formula = dx~genotype_A_AA+genotype_AA, family = binomial, data = SNP_perIID_genotypes))
  
  (OR12<-exp(cbind(coef(m12), confint(m12))))
  (ORAA<-exp(cbind(coef(mAAA), confint(mAAA))))
  
  print(OR12)
  print(OR12[2,])
  print(OR12[3,])
  results_per_snp[i,]<-c(snp, 
                         OR12[2,],coef(summary(m12))[2,4],  OR12[3,], coef(summary(m12))[3,4],
                         ORAA[2,],coef(summary(mAAA))[2,4],  ORAA[3,], coef(summary(mAAA))[3,4]
                         )
  print(results_per_snp[i,])
}


# saveRDS(results_per_snp, file = ORs.name, 
#         compress = TRUE)

data.table::fwrite(results_per_snp, file = ORs.name)
