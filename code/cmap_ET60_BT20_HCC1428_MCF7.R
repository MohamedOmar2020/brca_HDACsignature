
# clean ws
rm(list = ls())

# load necessary packages
library(data.table)
library(pheatmap)
library(cmapR)
library(readxl)

#############################
## load the ET signatures
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
ET60 <- All$`ET-60`[!is.na(All$`ET-60`)]

###############################
ds_path <- "./data/cmap/level5_beta_trt_cp_n720216x12328.gctx"
siginfo_path <- "./data/cmap/siginfo_beta.txt"
geneinfo_path <- "./data/cmap/geneinfo_beta.txt"

# read signature annotations (corresponding to columns of level 5 matrix)
siginfo <- data.table::fread(siginfo_path)
str(siginfo)

# read the metadata: make sure we read the gene_ids as characters
geneinfo <- data.table::fread(geneinfo_path, colClasses = c("gene_id" = "character"))

# Look up the corresponding gene_id, since this will correspond to row ids in the data matrix.
Gns_id <- geneinfo[gene_symbol %in% ET60]$gene_id
Gns_symbol <- geneinfo[gene_symbol %in% ET60]$gene_symbol

# consider only those signatures marked as exemplars, 
# which will ensure we get only one signature for each compound / cell line combination. 
# These exemplar signatures were selected as those with the best replicate recall and are indicated by the field is_ncs_exemplar.
sigs_of_interest <- siginfo[pert_type == "trt_cp" &
                              cell_iname %in% c("BT20", "HCC1428", "MCF7")]

sig_ids <- sigs_of_interest$sig_id
length(sig_ids)


# slice out the corresponding rows and columns from the data matrix using parse_gctx, 
# which will return an object of class GCT. 
# This function takes arguments rid and cid which we can use to specify which rows/columns we want to extract.
ds_Et60Gns <- cmapR::parse_gctx(ds_path,
                                rid = Gns_id,
                                cid = sig_ids)
ds_Et60Gns


###############################################################
###############################################################
## Identify modulators
# extract the data matrix
zs_mat <- mat(ds_Et60Gns)

# look at the z-score distributions. add vertical lines for the thresholds for modulation 
hist(zs_mat, col="dodgerblue", border="white", breaks=30, main=ET60[2], xlab="z-score")
abline(v=2, lty=2, col=2, lwd=1.3)
abline(v=-2, lty=2, col=2, lwd=1.3)

######
## identify the signatures in which the genes are significantly modulated:

mod_idx <- c()
for (i in Gns_id){
  mod_idx <- which(abs(zs_mat[i, ]) > 2)
}

length(mod_idx)

######
## Compare to non-modulators
# For comparison, we'll include roughly the same number of signatures in which the gene was NOT modulated.

non_mod_idx <- c()
for (i in Gns_id){
  non_mod_idx <- which(abs(zs_mat[i, ]) < 0.01)
}

length(non_mod_idx)

#####
# Slice out the signatures from the matrix, restricting to landmark space (i.e. directly measured genes only).
lm_ids <- geneinfo[feature_space == "landmark"]$gene_id

#intersect(lm_ids, Gns_id)
#ds2_Et60Gns <- parse_gctx(ds_path, rid=lm_ids,
#                          cid=c(names(mod_idx), names(non_mod_idx)))

ds2_Et60Gns <- parse_gctx(ds_path, rid=Gns_id,
                          cid=c(names(mod_idx), names(non_mod_idx)))

######
# Compute all pairwise correlations between signatures.
corr <- cor(mat(ds2_Et60Gns), method="spearman")


# ########
# ## Cluster the correlations, overlaying a color bar indicating the degree of modulation of our gene of interest.
# 
# mod_df <- list()
# for (i in Gns_id){
#   mod_df[[i]] <- data.frame(abs_zscore=abs(zs_mat[i, c(mod_idx, non_mod_idx)]))
#   rownames(mod_df[[i]]) <- c(names(mod_idx), names(non_mod_idx))
# }
# 
# names(mod_df) <- Gns_symbol
# 
# 
# ##############
# # plot the heatmaps
# par(pty="s")
# 
# for (i in Gns_symbol){
#   pheatmap::pheatmap(corr, annotation_row=mod_df[[i]], annotation_col=mod_df[[i]],
#                      show_rownames=F, show_colnames=F, main = i, 
#                      filename = paste0('./figs/cmap/ET60Gns_corr/', i, '.png')) 
# }




# look at the number of genes modulated in each class
ss_list <- list(
  mod = siginfo[sig_id %in% names(mod_idx)]$ss_ngene,
  non_mod = siginfo[sig_id %in% names(non_mod_idx)]$ss_ngene
)

boxplot(ss_list, main="Signature Strength",  ylab="# of genes modulated")


#############
# look at the more specific compounds (for the genes of interest) and restrict to those that have a canonical name.

compounds_ET60 <- siginfo[sig_id %in% names(mod_idx) &
                               ss_ngene == 0 &
                               qc_pass == 1 &
                               !grepl("^BRD-", cmap_name),
                             .(sig_id, pert_id, cmap_name, ss_ngene)][
                               order(ss_ngene, decreasing=F)]

write.csv(compounds_ET60, file = './objs/cmap/compounds_ET60.csv')

############################################################
## correlation heatmaps

final_comp <- compounds_ET60$sig_id

# ds for just the final compounds
ds2_Et60Gns_mod <- parse_gctx(ds_path, rid=Gns_id,
                              cid=final_comp)

## annotation

# annotate compounds
siginfo_mod <- siginfo[siginfo$sig_id %in% final_comp, ]
ds2_Et60Gns_mod_annot <- annotate_gct(ds2_Et60Gns_mod, siginfo_mod, dim="col",
                                      keyfield="sig_id")
# annotate genes
#geneinfo_mod <- geneinfo[geneinfo$gene_id %in% Gns_id, ]
#ds2_Et60Gns_mod_annot <- annotate_gct(ds2_Et60Gns_mod, geneinfo_mod, dim="row",
#                                      keyfield="gene_id")


#rownames(mat(ds2_Et60Gns_mod_annot)) <- ds2_Et60Gns_mod_annot@rdesc$gene_symbol
colnames(mat(ds2_Et60Gns_mod_annot)) <- ds2_Et60Gns_mod_annot@cdesc$cmap_name

# get the z score matrix
zs_mat_mod <- mat(ds2_Et60Gns_mod_annot)

# corr for the final compounds only
corr_mod <- cor(zs_mat_mod, method="spearman")
                              
########
## Cluster the correlations, overlaying a color bar indicating the degree of modulation of our gene of interest.

mod_df <- list()
for (i in Gns_id){
  mod_df[[i]] <- data.frame(abs_zscore=abs(zs_mat_mod[i, ]))
  mod_df[[i]] <- as.matrix(mod_df[[i]])
  rownames(mod_df[[i]]) <- colnames(zs_mat_mod)
  mod_df[[i]] <- as.data.frame(mod_df[[i]])
}

names(mod_df) <- Gns_symbol


##############
# plot the heatmaps
par(pty="s")

for (i in Gns_symbol){
  pheatmap::pheatmap(corr_mod, annotation_row=mod_df[[i]], annotation_col=mod_df[[i]],
                     show_rownames=T, show_colnames=T, main = i,
                     fontsize_row = 5, fontsize_col = 5,
                     width = 15, height = 15,
                     filename = paste0('./figs/cmap/ET60Gns_corr/', i, '.png')) 
}









# use the sig_id of these compounds to look up its impact (z-score) on ET60
compounds_ET60$direction <- rep('NA', nrow(compounds_ET60))
compounds_ET60 <- as.matrix(compounds_ET60)
rownames(compounds_ET60) <- compounds_ET60[, 'sig_id']


compounds_ET60[i, 'direction'] <- round(zs_mat[, i], 2)

for (i in compounds_ET60[, 'sig_id']){
  compounds_ET60[i, 'direction'] <- round(zs_mat[, i], 2)
}

### change gene ids to symbols
rownames(zs_mat) <- Gns_symbol
round(zs_mat[, "RAD001_MCF7_24H:BRD-A23723433-001-01-2:10"], 2)

###########################################
# read the column metadata
#col_meta <- read_gctx_meta(ds_path, dim="col")
#cid <- ids(ds_Et60Gns, dimension = "column")

# in memory slice using the cid parameter
#bumetanide_ds <- subset_gct(ds2_Et60Gns,
#                              cid=which(col_meta$pert_iname=="bumetanide"))
#identical(vemurafenib_ds, vemurafenib_ds3)

# note we need to specifiy which dimension to annotate (dim)
# and which column in the annotations corresponds to the column
# ids in the matrix (keyfield)
#ds2_Et60Gns_annot <- annotate_gct(ds2_Et60Gns, col_meta, dim="col",
#                               keyfield="id")

ds_Et60Gns_melted <- melt_gct(ds_Et60Gns)

ds_Et60Gns_melted <- ds_Et60Gns_melted[compounds_ET60$sig_id, ]

#################
# figure out which signatures correspond to vorinostat by searching the 'pert_iname' column
bumetanide_idx <- which(siginfo$cmap_name=="bumetanide")

# read only those columns from the GCTX file by using the 'cid' parameter
bumetanide_ds <- parse_gctx(ds_path, cid=bumetanide_idx, rid = Gns_id)
bumetanide_ds_melted <- melt_gct(bumetanide_ds)

##########
# plot the matrix values grouped by gene
library(ggplot2)

ggplot(bumetanide_ds_melted, aes(x=id.x, y=value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))







