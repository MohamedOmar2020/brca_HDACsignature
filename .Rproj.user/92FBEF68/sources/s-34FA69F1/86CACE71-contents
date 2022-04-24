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


# subset
temp <- SurvData_Genes_cbioportal[, ET3]

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
SurvData <- temp   
SurvData$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
SurvData$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status

# calculate survival
Fit <- surv_fit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ status, data = SurvData)
# get the logrank p-value
Pval <- surv_pvalue(Fit)
Pval <- Pval$pval

ggsurvplot(Fit)


#################################################################################
## Workon on the ET124 genes

ET125 <- intersect(All$`ET-125`, rownames(Expr_tcga))
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

#################################################

list_of_genes <- list_of_genes[-c(120, 121)]

surv_res <- lapply(list_of_genes, function(x){
  
  # subset
  #temp <- SurvData_Genes_cbioportal[, c(5:length(colnames(SurvData_Genes_cbioportal)))]
  temp <- SurvData_Genes_cbioportal[, x]
  
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df <- do.call(rbind, surv_res)
surv_df <- as.data.frame(surv_df)

for (i in colnames(surv_df)){
  surv_df[,i] <- as.character(surv_df[,i])
  surv_df
}

surv_df$pvalue <- as.numeric(surv_df$pvalue)   
surv_df <- surv_df[order(surv_df$pvalue), ]

ET36 <- surv_df[1, 'gene_names']
ET36 <- unlist(strsplit(ET36, split = ", "))

x <- intersect(ET36, All$`ET-9`)


#################################################################################
## Workon on the ET36 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET36
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(35, 36)]

surv_res_ET36 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET36]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET36 <- do.call(rbind, surv_res_ET36)
surv_df_ET36 <- as.data.frame(surv_df_ET36)

for (i in colnames(surv_df_ET36)){
  surv_df_ET36[,i] <- as.character(surv_df_ET36[,i])
  surv_df_ET36
}

surv_df_ET36$pvalue <- as.numeric(surv_df_ET36$pvalue)   
surv_df_ET36 <- surv_df_ET36[order(surv_df_ET36$pvalue), ]

ET32 <- surv_df_ET36[1, 'gene_names']
ET32 <- unlist(strsplit(ET32, split = ", "))

x <- intersect(ET32, All$`ET-9`)




#################################################################################
## Workon on the ET32 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET32
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(31, 32)]

surv_res_ET32 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET32]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET32 <- do.call(rbind, surv_res_ET32)
surv_df_ET32 <- as.data.frame(surv_df_ET32)

for (i in colnames(surv_df_ET32)){
  surv_df_ET32[,i] <- as.character(surv_df_ET32[,i])
  surv_df_ET32
}

surv_df_ET32$pvalue <- as.numeric(surv_df_ET32$pvalue)   
surv_df_ET32 <- surv_df_ET32[order(surv_df_ET32$pvalue), ]

ET28 <- surv_df_ET32[1, 'gene_names']
ET28 <- unlist(strsplit(ET30, split = ", "))

x <- intersect(ET28, All$`ET-9`)


#################################################################################
## Workon on the ET28 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET28
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(27, 28)]

surv_res_ET28 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET28]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET28 <- do.call(rbind, surv_res_ET28)
surv_df_ET28 <- as.data.frame(surv_df_ET28)

for (i in colnames(surv_df_ET28)){
  surv_df_ET28[,i] <- as.character(surv_df_ET28[,i])
  surv_df_ET28
}

surv_df_ET28$pvalue <- as.numeric(surv_df_ET28$pvalue)   
surv_df_ET28 <- surv_df_ET28[order(surv_df_ET28$pvalue), ]

ET25 <- surv_df_ET28[2, 'gene_names']
ET25 <- unlist(strsplit(ET25, split = ", "))

x <- intersect(ET25, All$`ET-9`)



#################################################################################
## Workon on the ET25 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET25
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(24, 25)]

surv_res_ET25 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET25]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET25 <- do.call(rbind, surv_res_ET25)
surv_df_ET25 <- as.data.frame(surv_df_ET25)

for (i in colnames(surv_df_ET25)){
  surv_df_ET25[,i] <- as.character(surv_df_ET25[,i])
  surv_df_ET25
}

surv_df_ET25$pvalue <- as.numeric(surv_df_ET25$pvalue)   
surv_df_ET25 <- surv_df_ET25[order(surv_df_ET25$pvalue), ]

ET23 <- surv_df_ET25[2, 'gene_names']
ET23 <- unlist(strsplit(ET23, split = ", "))

x <- intersect(ET23, All$`ET-9`)


#################################################################################
## Workon on the ET23 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET23
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(22, 23)]

surv_res_ET23 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET23]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET23 <- do.call(rbind, surv_res_ET23)
surv_df_ET23 <- as.data.frame(surv_df_ET23)

for (i in colnames(surv_df_ET23)){
  surv_df_ET23[,i] <- as.character(surv_df_ET23[,i])
  surv_df_ET23
}

surv_df_ET23$pvalue <- as.numeric(surv_df_ET23$pvalue)   
surv_df_ET23 <- surv_df_ET23[order(surv_df_ET23$pvalue), ]

ET21 <- surv_df_ET23[2, 'gene_names']
ET21 <- unlist(strsplit(ET21, split = ", "))

x <- intersect(ET21, All$`ET-9`)


#################################################################################
## Workon on the ET21 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET21
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(20, 21)]

surv_res_ET21 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET21]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET21 <- do.call(rbind, surv_res_ET21)
surv_df_ET21 <- as.data.frame(surv_df_ET21)

for (i in colnames(surv_df_ET21)){
  surv_df_ET21[,i] <- as.character(surv_df_ET21[,i])
  surv_df_ET21
}

surv_df_ET21$pvalue <- as.numeric(surv_df_ET21$pvalue)   
surv_df_ET21 <- surv_df_ET21[order(surv_df_ET21$pvalue), ]

ET20 <- surv_df_ET21[2, 'gene_names']
ET20 <- unlist(strsplit(ET20, split = ", "))

x <- intersect(ET20, All$`ET-9`)


#################################################################################
## Workon on the ET20 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET20
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(19, 20)]

surv_res_ET20 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET20]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET20 <- do.call(rbind, surv_res_ET20)
surv_df_ET20 <- as.data.frame(surv_df_ET20)

for (i in colnames(surv_df_ET20)){
  surv_df_ET20[,i] <- as.character(surv_df_ET20[,i])
  surv_df_ET20
}

surv_df_ET20$pvalue <- as.numeric(surv_df_ET20$pvalue)   
surv_df_ET20 <- surv_df_ET20[order(surv_df_ET20$pvalue), ]

ET18 <- surv_df_ET20[1, 'gene_names']
ET18 <- unlist(strsplit(ET18, split = ", "))

x <- intersect(ET18, All$`ET-9`)


#################################################################################
## Workon on the ET18 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET18
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(17, 18)]

surv_res_ET18 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET18]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET18 <- do.call(rbind, surv_res_ET18)
surv_df_ET18 <- as.data.frame(surv_df_ET18)

for (i in colnames(surv_df_ET18)){
  surv_df_ET18[,i] <- as.character(surv_df_ET18[,i])
  surv_df_ET18
}

surv_df_ET18$pvalue <- as.numeric(surv_df_ET18$pvalue)   
surv_df_ET18 <- surv_df_ET18[order(surv_df_ET18$pvalue), ]

ET15 <- surv_df_ET18[1, 'gene_names']
ET15 <- unlist(strsplit(ET15, split = ", "))

x <- intersect(ET15, All$`ET-9`)



#################################################################################
## Workon on the ET15 genes

# ======================
## categorize the Expr_tcgaession values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET15
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

list_of_genes <- list_of_genes[-c(14, 15)]

surv_res_ET15 <- lapply(list_of_genes, function(x){
  
  # subset
  temp <- SurvData_Genes_cbioportal[, ET15]
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
  temp$Progress.Free.Survival..Months. <- SurvData_Genes_cbioportal$Progress.Free.Survival..Months.
  temp$Progression.Free.Status <- SurvData_Genes_cbioportal$Progression.Free.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET15 <- do.call(rbind, surv_res_ET15)
surv_df_ET15 <- as.data.frame(surv_df_ET15)

for (i in colnames(surv_df_ET15)){
  surv_df_ET15[,i] <- as.character(surv_df_ET15[,i])
  surv_df_ET15
}

surv_df_ET15$pvalue <- as.numeric(surv_df_ET15$pvalue)   
surv_df_ET15 <- surv_df_ET15[order(surv_df_ET15$pvalue), ]

ET7 <- surv_df_ET15[2, 'gene_names']
ET7 <- unlist(strsplit(ET7, split = ", "))

ET7_ET9_inters <- intersect(ET7, All$`ET-9`)


##############################
## save the surv dataframes
save(surv_df, surv_df_ET36, surv_df_ET32, surv_df_ET28, surv_df_ET25, surv_df_ET23, 
     surv_df_ET21, surv_df_ET20, surv_df_ET18, surv_df_ET15, file = './objs/surv_df.rda')

# save the ET signatures
save(ET36, ET32, ET28, ET25, ET23, ET21, ET20, ET18, ET15, ET7, file = './objs/ET.rda')





