###### 
# Clean Work space
rm(list = ls())

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(mltools)
library(xtable)
library(readxl)
library(survminer)
library(survival)

#################
load("./objs/ET_metabric.rda")

## get the ET125 genes
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
# ET60 <- All$`ET-60`
# 
# ET60[ET60=='GPR56'] <- 'ADGRG1'
# ET60[ET60=='FAM116B'] <- 'DENND6B'
# #ET125[ET125=='FAM46B'] <- 'TENT5B'
# ET60[ET60=='HSA011916'] <- 'CTDNEP1'

# ET60 <- ET60[!is.na(ET60)]

myTSPs <- t(combn(All$`ET-60`, 2))
colnames(myTSPs) <- c("gene1", "gene2")

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')
###################################

### Common genes
keepGns_datasets <- intersect(rownames(Expr_tcga), rownames(Expr_metabric))
keepGns_tsp <- intersect(as.vector(myTSPs), keepGns_datasets)

Expr_metabric <- Expr_metabric[keepGns_tsp, ]
Expr_tcga <- Expr_tcga[keepGns_tsp, ]

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns_tsp & myTSPs[,2] %in% keepGns_tsp , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################
### Set Feature number and max k
ktsp <- c(3:10) #7
featNo <- nrow(Expr_metabric)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktsp_tcga <- SWAP.Train.KTSP(
  Expr_tcga, group_tcga, krange=ktsp, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktsp_tcga


badGns <- c(ktsp_tcga$TSPs[,1])
goodGns <- c(ktsp_tcga$TSPs[,2])

save(badGns, goodGns, file = "./objs/KTSPforCmap.rda")
############################################################
# Get the genes in common with ET9
ClassifierGenes <- as.vector(ktsp_tcga$TSPs)

# intersection with ET9
TSPs_ET9_int <- intersect(All$`ET-9`, ClassifierGenes)
TSPs_ET9_int

# ET4 intersection with ET9
TSPs_ET4_int <- intersect(ClassifierGenes, ET4)
TSPs_ET4_int

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStats_tcga <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_tcga, CombineFunc = sum)
summary(ktspStats_tcga$statistics)

### Threshold
thr <- coords(roc(group_tcga, ktspStats_tcga$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(group_tcga, ktspStats_tcga$statistics, levels = c(0, 1), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROC_tcga <- roc(group_tcga, ktspStats_tcga$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga

### Get predictions based on best threshold from ROC curve
prediction_tcga <- SWAP.KTSP.Classify(Expr_tcga, ktsp_tcga, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(prediction_tcga, group_tcga, positive = '1', mode = "everything")

MCC_tcga <- mltools::mcc(pred = prediction_tcga, actuals = group_tcga)
MCC_tcga

########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStats_metabric <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_tcga, CombineFunc = sum)
summary(ktspStats_metabric$statistics)

## Plot curve
ROC_metabric <- roc(group_metabric, ktspStats_metabric$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_metabric

### Get predictions based on best threshold from ROC curve
prediction_metabric <- SWAP.KTSP.Classify(Expr_metabric, ktsp_tcga, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusion_metabric <- confusionMatrix(prediction_metabric, group_metabric, positive = "1", mode = "everything")
confusion_metabric

MCC_metabric <- mltools::mcc(pred = prediction_metabric, actuals = group_metabric)
MCC_metabric
