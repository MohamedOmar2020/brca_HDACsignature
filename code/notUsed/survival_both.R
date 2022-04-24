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

Expr_metabric <- read.delim("./Data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt")
Pheno_metabric <- read.delim("./Data/brca_metabric_cbioportal/brca_metabric_clinical_data.tsv")

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
## Annotation: TCGA
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

## subset the Expr_tcgaession to ET125 genes
Expr_tcga <- Expr_tcga[rownames(Expr_tcga) %in% ET125, ]

# fix the column names
colnames(Expr_tcga) <- gsub('\\.01', '', colnames(Expr_tcga))

##############
## Annotation: metabric 
head(rownames(Expr_metabric))
Expr_metabric <- Expr_metabric[!duplicated(Expr_metabric$Hugo_Symbol), ]
rownames(Expr_metabric) <- Expr_metabric$Hugo_Symbol
Expr_metabric$Hugo_Symbol <- NULL
Expr_metabric$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_metabric)))
sel <- which(apply(Expr_metabric, 1, function(x) all(is.finite(x)) ))
Expr_metabric <- Expr_metabric[sel, ]
Expr_metabric <- Expr_metabric[!is.na(rownames(Expr_metabric)),]
dim(Expr_metabric)

range(Expr_metabric)
x <- intersect(ET125, rownames(Expr_metabric))

## subset the expression to ET125 genes
Expr_metabric <- Expr_metabric[rownames(Expr_metabric) %in% ET125, ]

`%!in%` <- Negate(`%in%`)
missingET125<- ET125[ET125 %!in% rownames(Expr_metabric)]

####################################
### Modify the Pheno_tcga

Pheno_tcga$X.Patient.Identifier <- gsub("\\-", "\\.", Pheno_tcga$X.Patient.Identifier)
rownames(Pheno_tcga) <- Pheno_tcga$X.Patient.Identifier
CommonSamples <- intersect(colnames(Expr_tcga), rownames(Pheno_tcga))
Pheno_tcga <- Pheno_tcga[CommonSamples, ]

Pheno_tcga$Progression.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_tcga$Overall.Survival.Status)

all(rownames(Pheno_tcga) == colnames(Expr_tcga))

Expr_tcga <- as.matrix(Expr_tcga)

############
### Modify the metabric phenotype

Pheno_metabric$Sample.ID <- gsub("\\-", "\\.", Pheno_metabric$Sample.ID)
rownames(Pheno_metabric) <- Pheno_metabric$Sample.ID
CommonSamples <- intersect(colnames(Expr_metabric), rownames(Pheno_metabric))
Pheno_metabric <- Pheno_metabric[CommonSamples, ]

Pheno_metabric$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_metabric$Overall.Survival.Status)

all(rownames(Pheno_metabric) == colnames(Expr_metabric))

Expr_metabric <- as.matrix(Expr_metabric)

############################################################
## Keep only the relevant information (Metastasis Event and Time)

# TCGA
Pheno_tcga$Progression.Free.Status <- as.numeric(Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- as.numeric(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Progression.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Progress.Free.Survival..Months. <- as.numeric(Pheno_tcga$Progress.Free.Survival..Months.)
Pheno_tcga$Overall.Survival..Months. <- as.numeric(Pheno_tcga$Overall.Survival..Months.)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Progression.Free.Status", "Progress.Free.Survival..Months.", "Overall.Survival.Status", "Overall.Survival..Months.")])

# create a merged pdata and Z-scores object
SurvData_Genes_cbioportal_tcga <- data.frame(Phenotype_tcga, t(Expr_tcga))

######
# Metabric 
Pheno_metabric$Relapse.Free.Status <- as.numeric(Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- as.numeric(Pheno_metabric$Overall.Survival.Status)

table(Pheno_metabric$Relapse.Free.Status)
table(Pheno_metabric$Overall.Survival.Status)

Pheno_metabric$Relapse.Free.Status..Months. <- as.numeric(Pheno_metabric$Relapse.Free.Status..Months.)
Pheno_metabric$Overall.Survival..Months. <- as.numeric(Pheno_metabric$Overall.Survival..Months.)

Phenotype_metabric <- cbind(Pheno_metabric[, c("Relapse.Free.Status", "Relapse.Free.Status..Months.", "Overall.Survival.Status", "Overall.Survival..Months.")])

# create a merged pdata and Z-scores object
SurvData_Genes_cbioportal_metabric <- data.frame(Phenotype_metabric, t(Expr_metabric))

#################################################################################
## Workon on the ET124 genes

ET125 <- intersect(All$`ET-125`, rownames(Expr_tcga))
ET125 <- intersect(ET125, rownames(Expr_metabric))

#ET125 <- sort(ET125, decreasing = T)
ET60 <- intersect(All$`ET-60`, rownames(Expr_tcga))

# ======================
## categorize the Expression values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET125
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}


###############################################
# gene_list <- list()
# combo_list <- list()
# for (i in 1:length(ET60)){
#   missing_gene <- ET60[i]
#   missing_gene
#   gene_list[[i]] <- list()
#   combo_list[[i]] <- list()
#   names(gene_list)[i] <- missing_gene
#   names(combo_list)[i] <- missing_gene
#   for (j in 1:(length(ET60)-1)){
#     loop_length <- length(ET60)-j
#     combs <- combn(x = ET60, m = loop_length)
#     combos <- as.list(data.frame(combs))
#     combos
#     combos_missing_gene <- combos[unlist(lapply(combos, function(x) !any(missing_gene %in% x)))]
#     gene_list[[i]][j] <- list(combos_missing_gene)
#   }
#   combo_list[[i]] <- unlist(gene_list[[i]],recursive = F)
# }
#################################################

list_of_genes <- list_of_genes[-c(114, 115)]

surv_res_ET125_tcga <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal_tcga[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal_tcga$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal_tcga$Progression.Free.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_res_ET125_metabric <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal_metabric[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal_metabric$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal_metabric$Overall.Survival.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET125_tcga <- do.call(rbind, surv_res_ET125_tcga)
surv_df_ET125_metabric <- do.call(rbind, surv_res_ET125_metabric)

#ET76 <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1", "KRT7", "LIMK2", "LOXL1", "MYO1B", "PCDH1", "PDLIM1", "SCNN1A", "SUSD2", "BNIPL", "BOC", "CCND2", "CPA4", "FAM46B", "FGF1", "IL6", "LOXL2", "LRP1", "MB", "MGAT1", "VIPR1", "XPC",  "ABCC2", "ABCC3", "ACSF2", "ANK1", "AZIN2", "B3GAT3", "BAD", "CARD14", "CCDC106", "CDC42EP3", "CERS4", "CHCHD6", "CLIC1", "CLIP2", "COL1A1", "COL5A1", "COQ8B", "CTSH", "DBP")
ET40_tcga <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1", "KRT7", "LIMK2", "LOXL1", "MYO1B")
ET39_metabric <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1", "KRT7", "LIMK2", "LOXL1")

ET9int_tcga <- intersect(ET40_tcga, All$`ET-9`)
ET9int_metabric <- intersect(ET39_metabric, All$`ET-9`)
both_int <- intersect(ET40_tcga, ET39_metabric)


#################################################################################
## Work on on the intersection between ET40_tcga and ET39_metabric

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- both_int
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(38, 39)]

surv_res_ET39_tcga <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal_tcga[, both_int]
  temp <- temp[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal_tcga$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal_tcga$Progression.Free.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})


surv_res_ET39_metabric <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal_metabric[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal_metabric$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal_metabric$Overall.Survival.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET39_tcga <- do.call(rbind, surv_res_ET39_tcga)
surv_df_ET39_metabric <- do.call(rbind, surv_res_ET39_metabric)

ET38_tcga <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1", "KRT7", "LIMK2")
ET38_metabric <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1", "KRT7", "LIMK2")

ET9int_tcga <- intersect(ET38_tcga, All$`ET-9`)
ET9int_metabric <- intersect(ET38_metabric, All$`ET-9`)
both_int <- intersect(ET38_tcga, ET38_metabric)


#################################################################################
## Work on on the ET38

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- both_int
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(37, 38)]

surv_res_ET38_tcga <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal_tcga[, both_int]
  temp <- temp[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal_tcga$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal_tcga$Progression.Free.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})


surv_res_ET38_metabric <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal_metabric[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal_metabric$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal_metabric$Overall.Survival.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET38_tcga <- do.call(rbind, surv_res_ET38_tcga)
surv_df_ET38_metabric <- do.call(rbind, surv_res_ET38_metabric)

ET36_tcga <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1")
ET36_metabric <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1")

ET9int_tcga <- intersect(ET36_tcga, All$`ET-9`)
ET9int_metabric <- intersect(ET36_metabric, All$`ET-9`)
both_int <- intersect(ET36_tcga, ET36_metabric)


#################################################################################
## Work on on the ET36

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- both_int
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(35, 36)]

surv_res_ET36_tcga <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal_tcga[, both_int]
  temp <- temp[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal_tcga$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal_tcga$Progression.Free.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})


surv_res_ET36_metabric <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal_metabric[, x]
  
  # convert each gene value to altered vs non-altered based on the expression
  for (g in seq_along(colnames(temp))){
    temp[, g] <- ifelse(temp[, g] > -2 & temp[, g] < 2, 'non-altered', 'altered')
  }
  
  
  temp$status <- rep(NA, length(rownames(temp)))
  
  
  # calculate the sample alteration status
  for (i in seq_along(rownames(temp))){
    temp$status[i] <- ifelse(any(temp[i, ] %in% 'altered'), 'altered', 'non-altered')
  }
  
  
  # re-add the survival info
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal_metabric$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal_metabric$Overall.Survival.Status
  
  # calculate survival
  Fit <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = temp)
  # get the logrank p-value
  Pval <- surv_pvalue(Fit)
  Pval <- Pval$pval
  #temp$Pval <- Pval
  #temp
  
  # store the results in a dataframe
  temp_res <- list(
    num_genes = ncol(temp)-3,
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], sep = ', '),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET36_tcga <- do.call(rbind, surv_res_ET36_tcga)
surv_df_ET36_metabric <- do.call(rbind, surv_res_ET36_metabric)

ET33_tcga <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1")
ET36_metabric <- c("ABTB1", "BCAS4", "CCDC69", "CROCC", "DKK3", "EPHB3", "FCHO1", "GDPD5", "IGFBP5", "MAP6", "PLA2G6", "RABGGTA", "SHKBP1", "TMEM53", "CSK", "CYP2S1", "DEF6", "ENO3", "FZD2", "HDAC11", "ID3", "IER5L", "LAMA5", "NAB2", "PREX1", "PRPF40B", "S100A6", "SHC1", "SNPH", "THBS1", "WNT10A", "CACNG4", "CUX1", "CX3CL1", "ECH1", "FIBCD1")

ET9int_tcga <- intersect(ET33_tcga, All$`ET-9`)
ET9int_metabric <- intersect(ET36_metabric, All$`ET-9`)
both_int <- intersect(ET36_tcga, ET36_metabric)





