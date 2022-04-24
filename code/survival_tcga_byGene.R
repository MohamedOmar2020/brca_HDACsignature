## Clean working space
rm(list = ls())

library(survminer)
library(survival)
library(patchwork)
library(readxl)

#######################
# Goal: refine the initial 125 genes to 60 genes (ET60) based on survival in the Metabric dataset

`%!in%` <- Negate(`%in%`)
################
# Load the METABRIC Expr_tcgaession and Pheno_tcgatype data
Expr_tcga <- read.delim("./Data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt")
Pheno_tcga <- read.delim("./Data/brca_tcga/data_clinical_patient.txt")
Pheno_tcga <- Pheno_tcga[-c(1:4), ]

#################
## get the ET124 genes
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
ET124 <- All$`ET-125`

ET124[ET124=='GPR56'] <- 'ADGRG1'
ET124[ET124=='FAM116B'] <- 'DENND6B'
#ET124[ET124=='FAM46B'] <- 'TENT5B'
ET124[ET124=='HSA011916'] <- 'CTDNEP1'

## get the ET60 genes
ET60 <- All$`ET-60`
ET60 <- ET60[!is.na(ET60)]

###########################
## Annotation
head(rownames(Expr_tcga))
Expr_tcga <- Expr_tcga[!duplicated(Expr_tcga$Hugo_Symbol), ]
rownames(Expr_tcga) <- Expr_tcga$Hugo_Symbol
Expr_tcga$Hugo_Symbol <- NULL
Expr_tcga$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga)))
sel <- which(apply(Expr_tcga, 1, function(x) all(is.finite(x)) ))
Expr_tcga <- Expr_tcga[sel, ]
Expr_tcga <- Expr_tcga[!is.na(rownames(Expr_tcga)),]
dim(Expr_tcga)

range(Expr_tcga)
#Expr_tcga <- log2(Expr_tcga + 6)
#Expr_tcga <- t(scale(t(Expr_tcga), center = TRUE, scale = TRUE))


## subset the Expr_tcgaession to ET124 genes
Expr_tcga <- Expr_tcga[rownames(Expr_tcga) %in% ET124, ]

# fix the column names
colnames(Expr_tcga) <- gsub('\\.01', '', colnames(Expr_tcga))

####################################
### Modify the Pheno_tcgatype

Pheno_tcga$X.Patient.Identifier <- gsub("\\-", "\\.", Pheno_tcga$X.Patient.Identifier)
rownames(Pheno_tcga) <- Pheno_tcga$X.Patient.Identifier
CommonSamples <- intersect(colnames(Expr_tcga), rownames(Pheno_tcga))
Pheno_tcga <- Pheno_tcga[CommonSamples, ]

Pheno_tcga$Progression.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Progression.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Progression.Free.Status)

all(rownames(Pheno_tcga) == colnames(Expr_tcga))

Expr_tcga <- as.matrix(Expr_tcga)

############################################################
## Keep only the relevant information (Metastasis Event and Time)

#table(Pheno_tcga$Relapse.Free.Status)
#Pheno_tcga$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Relapse.Free.Status)
Pheno_tcga$Progression.Free.Status <- as.numeric(Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Progression.Free.Status <- as.numeric(Pheno_tcga$Progression.Free.Status)

table(Pheno_tcga$Progression.Free.Status)
table(Pheno_tcga$Progression.Free.Status)

Pheno_tcga$Progress.Free.Survival..Months. <- as.numeric(Pheno_tcga$Progress.Free.Survival..Months.)
Pheno_tcga$Progress.Free.Survival..Months. <- as.numeric(Pheno_tcga$Progress.Free.Survival..Months.)

Phenotype <- cbind(Pheno_tcga[, c("Progression.Free.Status", "Progress.Free.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.")])

# create a merged pdata and Z-scores object
SurvData_Genes_cbioportal <- data.frame(Phenotype, t(Expr_tcga))

#################################################################################
## Workon on the ET124 genes

ET124 <- rownames(Expr_tcga)

######

for (i in ET124){
  SurvData_Genes_cbioportal[, i] <- ifelse(SurvData_Genes_cbioportal[, i] >= -2 & SurvData_Genes_cbioportal[, i] <= 2, 'non-altered', 'altered')
}


## For the positive/Up genes
# CutPoint_Genes <- surv_cutpoint(data = CoxData, time = "Progress.Free.Survival..Months.", event = "Progression.Free.Status", variables = ET124)
# CutPoint_Genes
# 
# SurvData_Genes <- surv_categorize(CutPoint_Genes)


########################################################################  
## Fit genes: Kaplan Mier logrank
ET124_list <- as.list(ET124)

surv_func <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(surv_fit(f, data = SurvData_Genes_cbioportal))
}

ET124_fit_list <- lapply(ET124_list, surv_func)
names(ET124_fit_list) <- ET124

# res_list_10yrs <- lapply(ET124_fit_list, summary, times = 120)
# 
# res_list_10yrs_pval <- lapply(res_list_10yrs, function(x){
#   diffSE <- sqrt(x$std.err[2]^2 + x$std.err[1]^2)
#   #a z-test test statistic is
#   zStat <- (x$surv[1] - x$surv[2])/diffSE
#   p <- 2*pnorm(abs(zStat), lower.tail=FALSE)
#   p
# })
# Pval_10yrs_df <- as.data.frame(do.call(rbind, res_list_10yrs_pval))
# Pval_10yrs_df$variable <- rownames(Pval_10yrs_df)
# colnames(Pval_10yrs_df) <- c('pval', 'variable')
# Pval_10yrs_df_fil <- Pval_10yrs_df[Pval_10yrs_df$pval < 0.05, ] 
# write.csv(Pval_10yrs_df_fil, 'objs/MetaBric_10yrsOS_sig.csv')
# CommGenes_10yrs_ET60 <- intersect(Pval_10yrs_df_fil$variable, ET60)
# CommGenes_10yrs_ET9 <- intersect(Pval_10yrs_df_fil$variable, All$`ET-9`)



Pval_list <- surv_pvalue(ET124_fit_list)

Pval_df <- do.call(rbind.data.frame, Pval_list)

Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

# save the results
write.csv(Pval_df_fil, 'objs/TCGA_PFS_byGene.csv')

CommGenes_ET60 <- intersect(Pval_df_fil$variable, ET60)
CommGenes_ET9 <- intersect(Pval_df_fil$variable, All$`ET-9`)

##################################################################
# fit genes: coxph model
ET124_list <- as.list(ET124)

surv_func_coxph <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(coxph(f, data = SurvData_Genes_cbioportal))
}

ET124_fit_list_coxph <- lapply(ET124_list, surv_func_coxph)
names(ET124_fit_list_coxph) <- ET124

ET124_summary_coxph <- lapply(ET124_fit_list_coxph, summary)

# get the HR
ET124_HR_coxph <- lapply(ET124_summary_coxph, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})




HR_df_coxph <- as.data.frame(do.call(rbind, ET124_HR_coxph))
HR_df_coxph$variable <- rownames(HR_df_coxph)
HR_df_coxph <- HR_df_coxph[order(HR_df_coxph$HR, decreasing = T), ]

HR_df_coxph_fil <- HR_df_coxph[HR_df_coxph$Pvalue_logrank_test < 0.05, ] 

# save the results
write.csv(HR_df_coxph, 'objs/TCGA_PFS_coxph.csv')

CommGenes_ET60_coxph <- HR_df_coxph[intersect(HR_df_coxph$variable, ET60), ]
CommGenes_ET9_coxph <- HR_df_coxph[intersect(HR_df_coxph$variable, All$`ET-9`), ]



