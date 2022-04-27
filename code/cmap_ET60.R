
# clean ws
rm(list = ls())

# load necessary packages
library(data.table)
library(pheatmap)
library(cmapR)

#############################
## load the ET signatures
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
ET60 <- All$`ET-60`[!is.na(All$`ET-60`)]

# load the ktsp signature based on ET60 in TCGA PFS
load("./objs/KTSPforCmap.rda")

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
goodGns_id <- geneinfo[gene_symbol %in% goodGns]$gene_id
goodGns_symbol <- geneinfo[gene_symbol %in% goodGns]$gene_symbol

badGns_id <- geneinfo[gene_symbol %in% badGns]$gene_id
badGns_symbol <- geneinfo[gene_symbol %in% badGns]$gene_symbol

# consider only those signatures marked as exemplars, 
# which will ensure we get only one signature for each compound / cell line combination. 
# These exemplar signatures were selected as those with the best replicate recall and are indicated by the field is_ncs_exemplar.
sigs_of_interest <- siginfo[pert_type == "trt_cp" &
                              cell_iname %in% c("BT20", "HCC1428", "MCF7") &
                              is_exemplar_sig == 1]

sig_ids <- sigs_of_interest$sig_id
length(sig_ids)


# slice out the corresponding rows and columns from the data matrix using parse_gctx, 
# which will return an object of class GCT. 
# This function takes arguments rid and cid which we can use to specify which rows/columns we want to extract.
ds_goodGns <- cmapR::parse_gctx(ds_path,
                        rid = goodGns_id,
                        cid = sig_ids)
ds_goodGns

ds_badGns <- cmapR::parse_gctx(ds_path,
                                rid = badGns_id,
                                cid = sig_ids)
ds_badGns

###############################################################
###############################################################
## Identify modulators
# extract the data matrix
zs_mat_goodGns <- mat(ds_goodGns)
zs_mat_badGns <- mat(ds_badGns)

# look at the z-score distributions. add vertical lines for the thresholds for modulation 
hist(zs_mat_goodGns, col="dodgerblue", border="white", breaks=30, main=badGns[2], xlab="z-score")
abline(v=1, lty=2, col=2, lwd=1.3)
abline(v=-1, lty=2, col=2, lwd=1.3)

######
## identify the signatures in which the genes are significantly modulated:

# good gns
mod_idx_goodGns <- c()
for (i in goodGns_id){
  mod_idx_goodGns <- which(abs(zs_mat_goodGns[i, ]) > 1)
}

length(mod_idx_goodGns)

# bad gns
mod_idx_badGns <- c()
for (i in badGns_id){
  mod_idx_badGns <- which(abs(zs_mat_badGns[i, ]) > 1)
}

length(mod_idx_badGns)

######
## Compare to non-modulators
# For comparison, we'll include roughly the same number of signatures in which the gene was NOT modulated.

# good gns
non_mod_idx_goodGns <- c()
for (i in goodGns_id){
  non_mod_idx_goodGns <- which(abs(zs_mat_goodGns[i, ]) < 0.15)
}

length(non_mod_idx_goodGns)

# bad gns
non_mod_idx_badGns <- c()
for (i in badGns_id){
  non_mod_idx_badGns <- which(abs(zs_mat_badGns[i, ]) < 0.15)
}

length(non_mod_idx_badGns)


#####
# Slice out the signatures from the matrix, restricting to landmark space (i.e. directly measured genes only).
lm_ids <- geneinfo[feature_space == "landmark"]$gene_id

# for good genes
ds2_goodGns <- parse_gctx(ds_path, rid=lm_ids,
                  cid=c(names(mod_idx_goodGns), names(non_mod_idx_goodGns)))

# for bad genes
ds2_badGns <- parse_gctx(ds_path, rid=lm_ids,
                          cid=c(names(mod_idx_badGns), names(non_mod_idx_badGns)))


######
# Compute all pairwise correlations between signatures.
corr_goodGns <- cor(mat(ds2_goodGns), method="spearman")
corr_badGns <- cor(mat(ds2_badGns), method="spearman")

########
## Cluster the correlations, overlaying a color bar indicating the degree of modulation of our gene of interest.

# for good genes
mod_df_goodGns <- list()
for (i in goodGns_id){
  mod_df_goodGns[[i]] <- data.frame(abs_zscore=abs(zs_mat_goodGns[i, c(mod_idx_goodGns, non_mod_idx_goodGns)]))
  rownames(mod_df_goodGns[[i]]) <- c(names(mod_idx_goodGns), names(non_mod_idx_goodGns))
}

names(mod_df_goodGns) <- goodGns_symbol

####
# for bad genes
mod_df_badGns <- list()
for (i in badGns_id){
  mod_df_badGns[[i]] <- data.frame(abs_zscore=abs(zs_mat_badGns[i, c(mod_idx_badGns, non_mod_idx_badGns)]))
  rownames(mod_df_badGns[[i]]) <- c(names(mod_idx_badGns), names(non_mod_idx_badGns))
}

names(mod_df_badGns) <- badGns_symbol

##############
# plot the heatmaps for the good genes
par(pty="s")

for (i in goodGns_symbol){
  pheatmap::pheatmap(corr_goodGns, annotation_row=mod_df_goodGns[[i]], annotation_col=mod_df_goodGns[[i]],
                     show_rownames=F, show_colnames=F, main = i, 
                     filename = paste0('./figs/cmap/goodGns/', i, '.png')) 
}

for (i in badGns_symbol){
  pheatmap::pheatmap(corr_badGns, annotation_row=mod_df_badGns[[i]], annotation_col=mod_df_badGns[[i]],
                     show_rownames=F, show_colnames=F, main = i, 
                     filename = paste0('./figs/cmap/badGns/', i, '.png')) 
}


# look at the number of genes modulated in each class
ss_list_goodGns <- list(
  mod = siginfo[sig_id %in% names(mod_idx_goodGns)]$ss_ngene,
  non_mod = siginfo[sig_id %in% names(non_mod_idx_goodGns)]$ss_ngene
)

boxplot(ss_list_goodGns, main="Signature Strength",  ylab="# of genes modulated")


ss_list_badGns <- list(
  mod = siginfo[sig_id %in% names(mod_idx_badGns)]$ss_ngene,
  non_mod = siginfo[sig_id %in% names(non_mod_idx_badGns)]$ss_ngene
)

boxplot(ss_list_badGns, main="Signature Strength",  ylab="# of genes modulated")


#############
# look at the more specific compounds (for the genes of interest) and restrict to those that have a canonical name.

compounds_goodGns <- siginfo[sig_id %in% names(mod_idx_goodGns) &
          ss_ngene <= 120 & 
            !grepl("^BRD-", cmap_name),
        .(sig_id, pert_id, cmap_name, ss_ngene)][
          order(ss_ngene, decreasing=F)]


compounds_badGns <- siginfo[sig_id %in% names(mod_idx_badGns) &
                               ss_ngene <= 120 & 
                               !grepl("^BRD-", cmap_name),
                             .(sig_id, pert_id, cmap_name, ss_ngene)][
                               order(ss_ngene, decreasing=F)]



















