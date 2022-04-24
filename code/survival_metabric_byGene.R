## Clean working space
rm(list = ls())

library(survminer)
library(survival)
library(patchwork)
library(readxl)

#######################
# Goal: refine the initial 125 genes to 60 genes (ET60) based on survival in the Metabric dataset

################
# Load the METABRIC expression and phenotype data
Expr <- read.delim("./Data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt")
Pheno <- read.delim("./Data/brca_metabric_cbioportal/brca_metabric_clinical_data.tsv")

#################
## get the ET125 genes
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
ET125 <- All$`ET-125`




## get the ET60 genes
ET60 <- All$`ET-60`
ET60 <- ET60[!is.na(ET60)]

###########################
## Annotation
head(rownames(Expr))
Expr <- Expr[!duplicated(Expr$Hugo_Symbol), ]
rownames(Expr) <- Expr$Hugo_Symbol
Expr$Hugo_Symbol <- NULL
Expr$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr)))
sel <- which(apply(Expr, 1, function(x) all(is.finite(x)) ))
Expr <- Expr[sel, ]
Expr <- Expr[!is.na(rownames(Expr)),]
dim(Expr)

range(Expr)
x <- intersect(ET125, rownames(Expr))
#Expr <- t(scale(t(Expr), center = TRUE, scale = TRUE))

## subset the expression to ET125 genes
Expr <- Expr[rownames(Expr) %in% ET125, ]

`%!in%` <- Negate(`%in%`)
missingET125<- ET125[ET125 %!in% rownames(Expr)]


# GPR56
####################################
### Modify the phenotype

Pheno$Sample.ID <- gsub("\\-", "\\.", Pheno$Sample.ID)
rownames(Pheno) <- Pheno$Sample.ID
CommonSamples <- intersect(colnames(Expr), rownames(Pheno))
Pheno <- Pheno[CommonSamples, ]

Pheno$Relapse.Free.Status <- gsub("\\:.+", "", Pheno$Relapse.Free.Status)
Pheno$Overall.Survival.Status <- gsub("\\:.+", "", Pheno$Overall.Survival.Status)

all(rownames(Pheno) == colnames(Expr))

Expr <- as.matrix(Expr)


############################################################
## Keep only the relevant information (Metastasis Event and Time)

#table(Pheno$Relapse.Free.Status)
#Pheno$Relapse.Free.Status <- gsub("\\:.+", "", Pheno$Relapse.Free.Status)
Pheno$Relapse.Free.Status <- as.numeric(Pheno$Relapse.Free.Status)
Pheno$Overall.Survival.Status <- as.numeric(Pheno$Overall.Survival.Status)

table(Pheno$Relapse.Free.Status)
table(Pheno$Overall.Survival.Status)

Pheno$Relapse.Free.Status..Months. <- as.numeric(Pheno$Relapse.Free.Status..Months.)
Pheno$Overall.Survival..Months. <- as.numeric(Pheno$Overall.Survival..Months.)

Phenotype <- cbind(Pheno[, c("Relapse.Free.Status", "Relapse.Free.Status..Months.", "Overall.Survival.Status", "Overall.Survival..Months.")])

# create a merged pdata and Z-scores object
CoxData <- data.frame(Phenotype, t(Expr))


#################################################################################
## Define optimal cutpoints for each gene (converting the absolute expression into categorical low/high expression)

ET115 <- rownames(Expr)


######
## categorize the expression values into altered and non-altered
SurvData_Genes_cbioportal <- CoxData

for (i in ET115){
  SurvData_Genes_cbioportal[, i] <- ifelse(SurvData_Genes_cbioportal[, i] >= -2 & SurvData_Genes_cbioportal[, i] <= 2, 'non-altered', 'altered')
}


## For the positive/Up genes
CutPoint_Genes <- surv_cutpoint(data = CoxData, time = "Overall.Survival..Months.", event = "Overall.Survival.Status", variables = ET115)
CutPoint_Genes

SurvData_Genes <- surv_categorize(CutPoint_Genes)


########################################################################  
## Fit genes: Kaplan Mier logrank
ET115_list <- as.list(ET115)

surv_func <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = SurvData_Genes))
}

ET115_fit_list <- lapply(ET115_list, surv_func)
names(ET115_fit_list) <- ET115

res_list_10yrs <- lapply(ET115_fit_list, summary, times = 120)

res_list_10yrs_pval <- lapply(res_list_10yrs, function(x){
  diffSE <- sqrt(x$std.err[2]^2 + x$std.err[1]^2)
  #a z-test test statistic is
  zStat <- (x$surv[1] - x$surv[2])/diffSE
  p <- 2*pnorm(abs(zStat), lower.tail=FALSE)
  p
})
Pval_10yrs_df <- as.data.frame(do.call(rbind, res_list_10yrs_pval))
Pval_10yrs_df$variable <- rownames(Pval_10yrs_df)
colnames(Pval_10yrs_df) <- c('pval', 'variable')
Pval_10yrs_df_fil <- Pval_10yrs_df[Pval_10yrs_df$pval < 0.05, ] 
write.csv(Pval_10yrs_df_fil, 'objs/MetaBric_10yrsOS_sig.csv')
CommGenes_10yrs_ET60 <- intersect(Pval_10yrs_df_fil$variable, ET60)
CommGenes_10yrs_ET9 <- intersect(Pval_10yrs_df_fil$variable, All$`ET-9`)



Pval_list <- surv_pvalue(ET115_fit_list)

Pval_df <- do.call(rbind.data.frame, Pval_list)

Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

# save the results
write.csv(Pval_df_fil, 'objs/MetaBric_OS_sig.csv')

CommGenes_ET60 <- intersect(Pval_df_fil$variable, ET60)
CommGenes_ET9 <- intersect(Pval_df_fil$variable, All$`ET-9`)


##################################################################
# fit genes: coxph model
ET115_list <- as.list(ET115)

surv_func_coxph <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(coxph(f, data = SurvData_Genes))
}

ET115_fit_list_coxph <- lapply(ET115_list, surv_func_coxph)
names(ET115_fit_list_coxph) <- ET115

ET115_summary_coxph <- lapply(ET115_fit_list_coxph, summary)
ET115_pval_coxph <- lapply(ET115_summary_coxph, function(x){
  x <- x$logtest['pvalue']
})


Pval_df_coxph <- as.data.frame(do.call(rbind, ET115_pval_coxph))
Pval_df_coxph$variable <- rownames(Pval_df_coxph)

Pval_df_fil_coxph <- Pval_df_coxph[Pval_df_coxph$pval < 0.05, ] 

# save the results
write.csv(Pval_df_fil_coxph, 'objs/MetaBric_OS_sig_coxph.csv')

CommGenes_ET60_coxph <- intersect(Pval_df_fil_coxph$variable, ET60)
CommGenes_ET9_coxph <- intersect(Pval_df_fil_coxph$variable, All$`ET-9`)







