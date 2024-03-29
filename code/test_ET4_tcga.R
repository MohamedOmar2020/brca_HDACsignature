## Clean working space
rm(list = ls())

library(survminer)
library(survival)
library(patchwork)
library(readxl)

#######################
# Goal: test the ET4 signature on the tcga cohort

# load the ET signatures
load("./objs/ET_metabric.rda")

#######################3
`%!in%` <- Negate(`%in%`)
################
# Load the expression and Phenotype data
Expr_tcga <- read.delim("./Data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt")
Pheno_tcga <- read.delim("./Data/brca_tcga/data_clinical_patient.txt")
Pheno_tcga <- Pheno_tcga[-c(1:4), ]

#################
## get the ET125 genes
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
ET125 <- All$`ET-125`

ET125[ET125=='GPR56'] <- 'ADGRG1'
ET125[ET125=='FAM116B'] <- 'DENND6B'
#ET125[ET125=='FAM46B'] <- 'TENT5B'
ET125[ET125=='HSA011916'] <- 'CTDNEP1'

## get the ET60 genes
ET60 <- All$`ET-60`
ET60 <- ET60[!is.na(ET60)]

x <- All$`ET-65`[!(All$`ET-65` %in% All$`ET-60`)]
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

## subset the Expr_tcgaession to ET125 genes
#Expr_tcga <- Expr_tcga[rownames(Expr_tcga) %in% ET5, ]

# fix the column names
colnames(Expr_tcga) <- gsub('\\.01', '', colnames(Expr_tcga))

####################################
### Modify the Pheno_tcgatype

Pheno_tcga$X.Patient.Identifier <- gsub("\\-", "\\.", Pheno_tcga$X.Patient.Identifier)
rownames(Pheno_tcga) <- Pheno_tcga$X.Patient.Identifier
CommonSamples <- intersect(colnames(Expr_tcga), rownames(Pheno_tcga))
Pheno_tcga <- Pheno_tcga[CommonSamples, ]

Pheno_tcga$Progression.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_tcga$Overall.Survival.Status)

all(rownames(Pheno_tcga) == colnames(Expr_tcga))

Expr_tcga <- as.matrix(Expr_tcga)

############################################################
## Keep only the relevant information (Metastasis Event and Time)

#table(Pheno_tcga$Relapse.Free.Status)
#Pheno_tcga$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Relapse.Free.Status)
Pheno_tcga$Progression.Free.Status <- as.numeric(Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- as.numeric(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Progression.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Progress.Free.Survival..Months. <- as.numeric(Pheno_tcga$Progress.Free.Survival..Months.)
Pheno_tcga$Overall.Survival..Months. <- as.numeric(Pheno_tcga$Overall.Survival..Months.)

Phenotype <- cbind(Pheno_tcga[, c("Progression.Free.Status", "Progress.Free.Survival..Months.", "Overall.Survival.Status", "Overall.Survival..Months.")])

# create a merged pdata and Z-scores object
SurvData_Genes_cbioportal <- data.frame(Phenotype, t(Expr_tcga))

##########################################################################
## surv curves

# subset
temp_ET37 <- SurvData_Genes_cbioportal[, ET37]
temp_ET33 <- SurvData_Genes_cbioportal[, ET33]
temp_ET4 <- SurvData_Genes_cbioportal[, ET4]

# convert each gene value to altered vs non-altered based on the expression
for (g in seq_along(colnames(temp_ET37))){
  temp_ET37[, g] <- ifelse(temp_ET37[, g] > -2 & temp_ET37[, g] < 2, 'non-altered', 'altered')
}

for (g in seq_along(colnames(temp_ET33))){
  temp_ET33[, g] <- ifelse(temp_ET33[, g] > -2 & temp_ET33[, g] < 2, 'non-altered', 'altered')
}

for (g in seq_along(colnames(temp_ET4))){
  temp_ET4[, g] <- ifelse(temp_ET4[, g] > -2 & temp_ET4[, g] < 2, 'non-altered', 'altered')
}

temp_ET37$status <- rep(NA, length(rownames(temp_ET37)))
temp_ET33$status <- rep(NA, length(rownames(temp_ET33)))
temp_ET4$status <- rep(NA, length(rownames(temp_ET4)))


# calculate the sample alteration status
for (i in seq_along(rownames(temp_ET37))){
  temp_ET37$status[i] <- ifelse(any(temp_ET37[i, ] %in% 'altered'), 'altered', 'non-altered')
}

for (i in seq_along(rownames(temp_ET33))){
  temp_ET33$status[i] <- ifelse(any(temp_ET33[i, ] %in% 'altered'), 'altered', 'non-altered')
}

for (i in seq_along(rownames(temp_ET4))){
  temp_ET4$status[i] <- ifelse(any(temp_ET4[i, ] %in% 'altered'), 'altered', 'non-altered')
}

# re-add the survival info
SurvData_ET37 <- temp_ET37   
SurvData_ET37$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
SurvData_ET37$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status

SurvData_ET33 <- temp_ET33   
SurvData_ET33$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
SurvData_ET33$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status


SurvData_ET4 <- temp_ET4   
SurvData_ET4$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
SurvData_ET4$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status

# calculate survival
Fit_ET37 <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = SurvData_ET37)
Fit_ET33 <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = SurvData_ET33)
Fit_ET4 <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = SurvData_ET4)
# get the logrank p-value
#Pval <- surv_pvalue(Fit)
#Pval <- Pval$pval

plot_ET37 <- ggsurvplot(Fit_ET37,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "ET37")


plot_ET33 <- ggsurvplot(Fit_ET33,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "ET33")

plot_ET4 <- ggsurvplot(Fit_ET4,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "ET4")


PlotList <- list(plot_ET37, plot_ET33, plot_ET4)
names(PlotList) <- c('ET37', 'ET33', 'ET4')

Splot <- arrange_ggsurvplots(PlotList, title = "TCGA Progression Free Survival", ncol = 3, nrow = 1)
ggsave("./figs/ET_tcga_pfs.png", Splot, width = 40, height = 15, units = "cm")



length(ET37) <- length(All$`ET-125`)
length(ET33) <- length(All$`ET-125`)
length(ET4) <- length(All$`ET-125`)

ET4sel <- data.frame(ET124 = All$`ET-125`, ET37 = ET37, ET33 = ET33, ET4 = ET4,
                     check.rows = FALSE)

write.csv(ET4sel, file = './objs/ET4sel.csv', col.names = T)



