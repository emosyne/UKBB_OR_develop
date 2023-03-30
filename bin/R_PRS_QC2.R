#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

#INPUT
args = commandArgs()

print(args)


hetfile = args[8]
target_bimfile = args[9]
HCM_GWAS = args[10]
(cohort = args[11])

#OUTPUT
het_valid_out = paste0(cohort,"_het_valid_out_vs_HCM_GWAS.sample")
restranded_recoded_cohort_bim = paste0(cohort,"_a1_cohort_bim_vs_HCM_GWAS")
mismatching_SNPs = paste0(cohort,"_mismatching_SNPs_vs_HCM_GWAS")
# annotated_mismatching_file = paste0(cohort,"_HCM_GWAS_vs_UKBB_annotated_mismatching.csv")


#Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. 

# Read in file
dat <- fread(hetfile)

# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)]

#remove "minus" subjects (withdrawn)
valid<-dplyr::filter(valid, !grepl("^-",FID))


# print FID and IID for valid samples (HET and not minus)
fwrite(valid[,c("FID","IID")], het_valid_out, sep="\t")



### Mismatching SNPs

## 1. Load the cohort_bim file, the summary statistic and the QC SNP list into R

# Read in QCed cohort_bim file
cohort_bim <- read.table(target_bimfile)
colnames(cohort_bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
cohort_bim
# Read in the GWAS data
GWAS <- fread(file=HCM_GWAS, select = c("CHR", "SNP", "POS", "A1", "A2",  "P"))
GWAS <- GWAS %>% rename(BP=POS)

  # read.table(gzfile(HCM_GWAS),
  #            header = T,
  #            stringsAsFactors = F, 
  #            sep="\t")


# Change all alleles to upper case for easy comparison
GWAS$A1 <- toupper(GWAS$A1)
GWAS$A2 <- toupper(GWAS$A2)
cohort_bim$B.A1 <- toupper(cohort_bim$B.A1)
cohort_bim$B.A2 <- toupper(cohort_bim$B.A2)

## 2. Identify SNPs that require strand flipping 

# Merge summary statistic with target
print("head info")
info <- merge(cohort_bim, GWAS, by = c("SNP", "CHR", "BP"))
head(info)

# Function for finding the complementary allele
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}
# Get SNPs that have the same alleles across base and target
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
# Identify SNPs that are complementary between base and target
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Update the complementary alleles in the cohort_bim file
# This allow us to match the allele in subsequent analysis
complement.snps <- cohort_bim$SNP %in% info.complement$SNP
cohort_bim[complement.snps,]$B.A1 <-
  sapply(cohort_bim[complement.snps,]$B.A1, complement)
cohort_bim[complement.snps,]$B.A2 <-
  sapply(cohort_bim[complement.snps,]$B.A2, complement)

## 3. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)

# identify SNPs that need recoding
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
# Update the recode SNPs
recode.snps <- cohort_bim$SNP %in% info.recode$SNP
tmp <- cohort_bim[recode.snps,]$B.A1
cohort_bim[recode.snps,]$B.A1 <- cohort_bim[recode.snps,]$B.A2
cohort_bim[recode.snps,]$B.A2 <- tmp

# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
com.snps <- cohort_bim$SNP %in% info.crecode$SNP
tmp <- cohort_bim[com.snps,]$B.A1
cohort_bim[com.snps,]$B.A1 <- as.character(sapply(cohort_bim[com.snps,]$B.A2, complement))
cohort_bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))

# Output updated cohort_bim file
write.table(
  cohort_bim[,c("SNP", "B.A1")],
  restranded_recoded_cohort_bim,
  quote = F,
  row.names = F,
  col.names = F,
  sep="\t"
)

##4. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)

mismatch <-
  cohort_bim$SNP[!(cohort_bim$SNP %in% info.match$SNP |
              cohort_bim$SNP %in% info.complement$SNP | 
              cohort_bim$SNP %in% info.recode$SNP |
              cohort_bim$SNP %in% info.crecode$SNP)]
write.table(
  mismatch,
  mismatching_SNPs,
  quote = F,
  row.names = F,
  col.names = F
)

# print("gen annotated_mismatching")
# annotated_mismatching = info %>% mutate(mismatch = 
#     ifelse(test=(SNP %in% mismatch), "mismatch","match")
#   ) %>% select(-C.A1,-C.A2)
# print(" head annotated_mismatching")
# head(annotated_mismatching)
# #SNP CHR        BP CM B.A1 B.A2 A1 A2      P      OR tidytype mismatch
# colnames(annotated_mismatching) = c("SNP",	"CHR",	"BP",	"CM",	"UKBB.A1",	"UKBB.A2",	"PGC_A1",	"PGC_A2",	"P",	"OR", "tidytype", "mismatch")

# fwrite(annotated_mismatching, file = annotated_mismatching_file)
