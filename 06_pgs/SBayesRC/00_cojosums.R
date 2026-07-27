# library 


library(data.table)


ppd_f<- fread("03_GSEM/gsem_model/PPD_resid_F_results.txt") 


## Make COJO-style table
## COJO canonical columns: SNP A1 A2 freq b se p N
cojo <- ppd_f[, .(
  SNP,                # rsID
  A1,                 # effect allele
  A2,                 # other allele
  freq = MAF,         # allele frequency for A1 (assumes MAF is for A1)
  b    = est,         # effect size
  se   = SE,          # standard error
  p    = Pval_Estimate         
)]



# Basic QC
summary(dat$freq)
summary(dat$b)
summary(dat$se)
summary(dat$p)
summary(dat$N)


Neff <- median(ppd_f$Neff)
cojo[, N := Neff]

fwrite(cojo, "SBayesRC/COJO/ppd_f_sumstats_cojo.ma")  



######################################################################################################

mdd_f<- fread("03_GSEM/gsem_model/MDD_F_shared_results.txt") 

## COJO canonical columns: SNP A1 A2 freq b se p N
cojo <- mdd_f[, .(
  SNP,                # rsID
  A1,                 # effect allele
  A2,                 # other allele
  freq = MAF,         # allele frequency for A1 (assumes MAF is for A1)
  b    = est,         # effect size
  se   = SE,          # standard error
  p    = Pval_Estimate         
)]

Neff <- median(mdd_f$Neff)
cojo[, N := Neff]

fwrite(cojo, "SBayesRC/COJO/mdd_f_sumstats_cojo.ma")  

###########################################################################################
ppd_org<- fread("GWAS_QC/ppd.qc.sumstats.gz") 

## COJO canonical columns: SNP A1 A2 freq b se p N
cojo <- ppd_org[, .(
  SNP,                # rsID
  A1,                 # effect allele
  A2,                 # other allele
  freq = MAF,         # allele frequency for A1 (assumes MAF is for A1)
  b    = BETA,         # effect size
  se   = SE,          # standard error
  p    = P,
  N         
)]



fwrite(cojo, "SBayesRC/COJO/ppd_orginal_sumstats_cojo.ma")  


####################################################################################
mdd_f_org<- fread("GWAS_QC/mdd_f.qc.sumstats.gz") 

## COJO canonical columns: SNP A1 A2 freq b se p N
cojo <- mdd_f_org[, .(
  SNP,                # rsID
  A1,                 # effect allele
  A2,                 # other allele
  freq = MAF,         # allele frequency for A1 (assumes MAF is for A1)
  b    = BETA,         # effect size
  se   = SE,          # standard error
  p    = P,
  N         
)]



fwrite(cojo, "SBayesRC/COJO/mdd_f_orginal_sumstats_cojo.ma")  
































































