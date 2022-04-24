## Clean working space
rm(list = ls())

library(survminer)
library(survival)
library(patchwork)
library(readxl)

#######################

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

################
# Load the METABRIC Expr_metabricession and Pheno_metabrictype data
Expr_metabric <- read.delim("./Data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt")
Pheno_metabric <- read.delim("./Data/brca_metabric_cbioportal/brca_metabric_clinical_data.tsv")

###########################
## Annotation
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
x <- intersect(ET7, rownames(Expr_metabric))
#Expr_metabric <- t(scale(t(Expr_metabric), center = TRUE, scale = TRUE))

## subset the Expr_metabricession to ET125 genes
Expr_metabric <- Expr_metabric[rownames(Expr_metabric) %in% ET125, ]

#`%!in%` <- Negate(`%in%`)
#missingET7<- ET7[ET7 %!in% rownames(Expr_metabric)]


# GPR56
####################################
### Modify the Pheno_metabric

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

#table(Pheno_metabric$Relapse.Free.Status)
#Pheno_metabric$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Relapse.Free.Status <- as.numeric(Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- as.numeric(Pheno_metabric$Overall.Survival.Status)

table(Pheno_metabric$Relapse.Free.Status)
table(Pheno_metabric$Overall.Survival.Status)

Pheno_metabric$Relapse.Free.Status..Months. <- as.numeric(Pheno_metabric$Relapse.Free.Status..Months.)
Pheno_metabric$Overall.Survival..Months. <- as.numeric(Pheno_metabric$Overall.Survival..Months.)

Phenotype<- cbind(Pheno_metabric[, c("Relapse.Free.Status", "Relapse.Free.Status..Months.", "Overall.Survival.Status", "Overall.Survival..Months.")])


# create a merged pdata and Z-scores object
SurvData_Genes_cbioportal <- data.frame(Phenotype, t(Expr_metabric))


#################################################################################
## Workon on the ET124 genes

ET125 <- intersect(All$`ET-125`, rownames(Expr_metabric))
#ET125 <- sort(ET125, decreasing = T)
ET60 <- intersect(All$`ET-60`, rownames(Expr_metabric))

# ======================
## categorize the Expression values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET125
for (x in seq_along(temp_ET)){
  #list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

##########

list_of_genes <- list_of_genes[-c(114, 115)]

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
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status
  
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

ET37 <- surv_df[6, 'gene_names']
ET37 <- unlist(strsplit(ET37, split = ", "))

x <- intersect(ET37, All$`ET-9`)


#################################################################################
## Workon on the ET37 genes

# ======================
## categorize the Expression values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET37
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

################

list_of_genes <- list_of_genes[-c(36, 37)]

surv_res_ET37 <- lapply(list_of_genes, function(x){
  
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
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET37 <- do.call(rbind, surv_res_ET37)
surv_df_ET37 <- as.data.frame(surv_df_ET37)

for (i in colnames(surv_df_ET37)){
  surv_df_ET37[,i] <- as.character(surv_df_ET37[,i])
  surv_df_ET37
}

surv_df_ET37$pvalue <- as.numeric(surv_df_ET37$pvalue)   
surv_df_ET37 <- surv_df_ET37[order(surv_df_ET37$pvalue), ]

ET33 <- surv_df_ET37[2, 'gene_names']
ET33 <- unlist(strsplit(ET33, split = ", "))

x <- intersect(ET33, All$`ET-9`)


#################################################################################
## Workon on the ET33 genes

# ======================
## categorize the Expression values into altered and non-altered
list_of_genes <- list()
temp_ET <- ET33
for (x in seq_along(temp_ET)){
  list_of_genes[[x]] <- sample(temp_ET, size = length(temp_ET)-1,replace = FALSE)
  #list_of_genes[[x]] <- temp_ET[-length(temp_ET)]
  temp_ET <- list_of_genes[[x]] 
}

################

list_of_genes <- list_of_genes[-c(32, 33)]

surv_res_ET33 <- lapply(list_of_genes, function(x){
  
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
  temp$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
  temp$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status
  
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
    gene_names = paste(colnames(temp)[1:(length(colnames(temp))-3)], collapse=", "),
    pvalue = Pval
  )
  
  temp_res
  
})

surv_df_ET33 <- do.call(rbind, surv_res_ET33)
surv_df_ET33 <- as.data.frame(surv_df_ET33)

for (i in colnames(surv_df_ET33)){
  surv_df_ET33[,i] <- as.character(surv_df_ET33[,i])
  surv_df_ET33
}

surv_df_ET33$pvalue <- as.numeric(surv_df_ET33$pvalue)   
surv_df_ET33 <- surv_df_ET33[order(surv_df_ET33$pvalue), ]

ET4 <- surv_df_ET33[3, 'gene_names']
ET4 <- unlist(strsplit(ET4, split = ", "))

ET3 <- intersect(ET4, All$`ET-9`)

###############################################################################
## save the results
save(surv_df, surv_df_ET37, surv_df_ET33, file = './objs/surv_df_metabric.rda')
save(ET37, ET33, ET3, ET4, file = './objs/ET_metabric.rda')

# load
load('./objs/ET_metabric.rda')
load('./objs/surv_df_metabric.rda')
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
SurvData_ET37$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
SurvData_ET37$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status

SurvData_ET33 <- temp_ET33   
SurvData_ET33$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
SurvData_ET33$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status


SurvData_ET4 <- temp_ET4   
SurvData_ET4$Overall.Survival..Months. <- SurvData_Genes_cbioportal$Overall.Survival..Months.
SurvData_ET4$Overall.Survival.Status <- SurvData_Genes_cbioportal$Overall.Survival.Status

# calculate survival
Fit_ET37 <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = SurvData_ET37)
Fit_ET33 <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = SurvData_ET33)
Fit_ET4 <- surv_fit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ status, data = SurvData_ET4)
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

Splot <- arrange_ggsurvplots(PlotList, title = "Metabric Overall Survival", ncol = 3, nrow = 1)
ggsave("./figs/ET_Metabric_os.png", Splot, width = 40, height = 15, units = "cm")





