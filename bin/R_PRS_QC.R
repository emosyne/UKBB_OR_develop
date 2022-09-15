#!/usr/bin/env Rscript
library(data.table)

#INPUT
args = commandArgs()

print(args)

(condition = args[8])
hetfile = args[9]
bimfile = args[10]
GWAS = args[11]


#OUTPUT
het_valid_out = paste0(condition, "_het_valid_out.sample")
restranded_recoded_BIM = paste0(condition, "_a1_bim")
mismatching_SNPs = paste0(condition, "_mismatching_SNPs")


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

## 1. Load the bim file, the summary statistic and the QC SNP list into R

# Read in QCed bim file
bim <- read.table(bimfile)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
bim
# Read in the GWAS data
GWAS <-
  read.table(gzfile(GWAS),
             header = T,
             stringsAsFactors = F, 
             sep="\t")
colnames(GWAS) <- c("CHR", "SNP", "BP","A1", "A2",  "P", "OR")
# Change all alleles to upper case for easy comparison
GWAS$A1 <- toupper(GWAS$A1)
GWAS$A2 <- toupper(GWAS$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)

## 2. Identify SNPs that require strand flipping 

# Merge summary statistic with target
info <- merge(bim, GWAS, by = c("SNP", "CHR", "BP"))

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
# Update the complementary alleles in the bim file
# This allow us to match the allele in subsequent analysis
complement.snps <- bim$SNP %in% info.complement$SNP
bim[complement.snps,]$B.A1 <-
  sapply(bim[complement.snps,]$B.A1, complement)
bim[complement.snps,]$B.A2 <-
  sapply(bim[complement.snps,]$B.A2, complement)

## 3. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)

# identify SNPs that need recoding
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
# Update the recode SNPs
recode.snps <- bim$SNP %in% info.recode$SNP
tmp <- bim[recode.snps,]$B.A1
bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
bim[recode.snps,]$B.A2 <- tmp

# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
com.snps <- bim$SNP %in% info.crecode$SNP
tmp <- bim[com.snps,]$B.A1
bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))

# Output updated bim file
write.table(
  bim[,c("SNP", "B.A1")],
  restranded_recoded_BIM,
  quote = F,
  row.names = F,
  col.names = F,
  sep="\t"
)

##4. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)

mismatch <-
  bim$SNP[!(bim$SNP %in% info.match$SNP |
              bim$SNP %in% info.complement$SNP | 
              bim$SNP %in% info.recode$SNP |
              bim$SNP %in% info.crecode$SNP)]
write.table(
  mismatch,
  mismatching_SNPs,
  quote = F,
  row.names = F,
  col.names = F
)