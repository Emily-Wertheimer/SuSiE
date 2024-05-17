#--------------------
# title:"SuSiE on all sumstats"
# author: "Joy Lee"
# date: "10/29/2023"
# version: V3 use new 1000G.EUR
#--------------------
## install and load package
# install.packages("susieR")
library(susieR)
library(data.table)

### MS
#--------
## import MS sum stats data
# MS summary stats
MS_meta <- read.table("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Summary statistics/MS/discovery_metav3.0.meta", header = TRUE)
# remove SNP that does not have p-value
MS_meta_P <- MS_meta[complete.cases(MS_meta$P),]
# subset SNP reported in more than 7 datasets
MS_meta_P_N <- MS_meta_P[MS_meta_P$N > 7,]
# z-score calculation
MS_meta_P_N$z_scores <- qnorm(MS_meta_P_N$P/2)

## run SuSiE
# set SNP calculator into 0
nSNP = 0
# read in SNP in the clumping data and LD matrix
# run SuSiE and record the number of SNP
for (locus in 1:129) {
  print(locus)
  # import clumping and LD matrix
  plink_clumping_SNP <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_clumping/MS_clump/MS_locus_", locus, ".txt"), header = FALSE)
  colnames(plink_clumping_SNP) <- c("SNP")
  plink_clumping_matrix <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_calculation/MS_LD/LDmatrix_MS_", locus, ".ld"), header = FALSE)
  # prepare data for SuSiE run
  MS_SNP_list <- MS_meta_P_N[MS_meta_P_N$SNP %in% plink_clumping_SNP$SNP, ]
  matrix <- data.matrix(plink_clumping_matrix)
  # check LD matrix and SumStats correlation
  lambda <- estimate_s_rss(MS_SNP_list$z_scores, matrix, n=38000)
  print(lambda)
  condz_in <- kriging_rss(MS_SNP_list$z_scores, matrix, n=38000)
  #condz_in$plot
  index <- which(abs(condz_in$conditional_dist$z_std_diff)>3)
  # remove index that does not match between SumStats and LD matrix
  if (lambda > 0.01 & length(index) > 0 & nrow(MS_SNP_list ) > 10) {
    MS_SNP_list_remove <- MS_SNP_list[-index,]
    matrix_remove <- matrix[-index,-index]
    lambda_remove <- estimate_s_rss(MS_SNP_list_remove$z_scores, matrix_remove, n=38000)
    print(lambda_remove)
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(MS_SNP_list_remove$z_scores, matrix_remove, n = 38000, L = 5)
  } else {
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(MS_SNP_list$z_scores, matrix, n = 38000, L = 5)
  }
  # record number of SNP
  for (i in 1:5) {
    CS95 = paste0("fitted_rss$sets$cs$L",i)
    nSNP <- nSNP + print(length(eval(parse(text = CS95))))
  }
}

print(nSNP)

## result: 
# 2200 SNPs (5 out of 129 locus cannot be resolved, start with 30209 SNP with p < 10^(-6))
# locus 97, 117, 122, 127,128 canot be resolved

# before fine-map
plot(-log(MS_SNP_list$P), ylab = "-logP", xlab="variants", xaxt="n", pch=19, col="black")
axis(side=1, at=c(1:nrow(MS_SNP_list)), labels=FALSE)
abline(h = -log(5e-08), col = "red")
text(x=1:nrow(MS_SNP_list),y = par("usr")[3]-1.1,labels=as.graphicsAnnot(MS_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)

susie_plot(fitted_rss, y="PIP", b=log(MS_SNP_list$OR), add_legend = "topleft", xlab="variants", main="SuSiE credible sets 1000G.EUR", xaxt="n")
axis(side=1, at=c(1:nrow(MS_SNP_list)), labels=FALSE)
text(x=1:nrow(MS_SNP_list),y = par("usr")[3]-0.01,labels=as.graphicsAnnot(MS_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)
#text(x=1:nrow(MS_SNP_list),y = par("usr")[3]-0.05,labels=as.graphicsAnnot(MS_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)
#----





### AD
#--------
## read data in chunk
file_path <- "/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Summary statistics/AD/GCST90027158_buildGRCh38.tsv"
AD_meta <- fread(file_path, header = TRUE)
# z-score calculation
#AD_meta$z_scores <- AD_meta$beta / AD_meta$standard_error
AD_meta$z_scores <- qnorm(as.numeric(AD_meta$p_value)/2)

## run SuSiE
# set SNP calculator into 0
nSNP = 0
# skipping locus
skipping <- c(11, 37, 181)
skipping_nSNP <- 3 + 2 + 3

# read in SNP in the clumping data and LD matrix 1:223
# run SuSiE and record the number of SNP
for (locus in 1:128) {
  if (locus %in% skipping) {
    next  # Skips the current iteration when locus is 5
  }
  print(locus)
  # import clumping and LD matrix
  plink_clumping_SNP <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_clumping/AD_clump/AD_locus_", locus, ".txt"), header = FALSE)
  colnames(plink_clumping_SNP) <- c("SNP")
  plink_clumping_matrix <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_calculation/AD_LD/LDmatrix_AD_", locus, ".ld"), header = FALSE)
  # prepare data for SuSiE run
  AD_SNP_list <- AD_meta[AD_meta$variant_id %in% plink_clumping_SNP$SNP, ]
  matrix <- data.matrix(plink_clumping_matrix)
  # check LD matrix and SumStats correlation
  lambda <- estimate_s_rss(AD_SNP_list$z_scores, matrix, n=50000)
  print(lambda)
  condz_in <- kriging_rss(AD_SNP_list$z_scores, matrix, n=50000)
  condz_in$plot
  index <- which(abs(condz_in$conditional_dist$z_std_diff)>3)
  # remove index that does not match between SumStats and LD matrix
  if (lambda > 0.01 & length(index) > 0 & nrow(AD_SNP_list) > 10) {
    AD_SNP_list_remove <- AD_SNP_list[-index,]
    matrix_remove <- matrix[-index,-index]
    lambda_remove <- estimate_s_rss(AD_SNP_list_remove$z_scores, matrix_remove, n=38000)
    print(lambda_remove)
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(AD_SNP_list_remove$z_scores, matrix_remove, n = 38000, L = 5)
  } else {
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(AD_SNP_list$z_scores, matrix, n = 38000, L = 5)
  }
  # record number of SNP
  for (i in 1:5) {
    CS95 = paste0("fitted_rss$sets$cs$L",i)
    nSNP <- nSNP + print(length(eval(parse(text = CS95))))
  }
}

print(nSNP)
# 6 locus cannot resolve (3193)

# before fine-map
plot(-log(AD_SNP_list$p_value), ylab = "-logP", xlab="variants", xaxt="n", pch=19, col="black")
axis(side=1, at=c(1:nrow(AD_SNP_list)), labels=FALSE)
abline(h = -log(5e-08), col = "red")
text(x=1:nrow(AD_SNP_list),y = par("usr")[3]-1.1,labels=as.graphicsAnnot(AD_SNP_list$variant_id),xpd=NA,srt=45,cex=0.65,adj=1)
# SuSiE plot
susie_plot(fitted_rss, y="PIP", b=log(AD_SNP_list$odds_ratio), add_legend = "topleft", xlab="variants", main="SuSiE credible sets 1000G.EUR", xaxt="n")
axis(side=1, at=c(1:nrow(AD_SNP_list)), labels=FALSE)
text(x=1:nrow(AD_SNP_list),y = par("usr")[3]-0.01,labels=as.graphicsAnnot(AD_SNP_list$variant_id),xpd=NA,srt=45,cex=0.65,adj=1)
# plot LD
library(LDheatmap)
library(grid)
matrix_LDmap <- matrix*matrix
colnames(matrix_LDmap) <- AD_SNP_list$SNP
rownames(matrix_LDmap) <- AD_SNP_list$SNP
LDheatmap(matrix_LDmap)

#----





### PD
#--------
## import PD sum stats data
# PD summary stats
PD_meta <- fread("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Summary statistics/PD/PD_SumStats_rsID.tab", header = TRUE, sep = "\t")
# z-score calculation
PD_meta$z_scores <- PD_meta$b / PD_meta$se
PD_meta$z_scores_P <- qnorm(PD_meta$P/2)
## run SuSiE
# set SNP calculator into 0
nSNP = 0
# read in SNP in the clumping data and LD matrix
# run SuSiE and record the number of SNP
for (locus in 1:42) {
  print(locus)
  # import clumping and LD matrix
  plink_clumping_SNP <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_clumping/PD_clump/PD_locus_", locus, ".txt"), header = FALSE)
  colnames(plink_clumping_SNP) <- c("SNP")
  plink_clumping_matrix <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_calculation/PD_LD/LDmatrix_PD_", locus, ".ld"), header = FALSE)
  # prepare data for SuSiE run
  PD_SNP_list <- PD_meta[PD_meta$SNP %in% plink_clumping_SNP$SNP, ]
  matrix <- data.matrix(plink_clumping_matrix)
  # check LD matrix and SumStats correlation
  lambda <- estimate_s_rss(PD_SNP_list$z_scores, matrix, n=38000)
  print(lambda)
  condz_in <- kriging_rss(PD_SNP_list$z_scores, matrix, n=38000)
  condz_in$plot
  index <- which(abs(condz_in$conditional_dist$z_std_diff)>3)
  # remove index that does not match between SumStats and LD matrix
  if (lambda > 0.01 & length(index) > 0 & nrow(PD_SNP_list ) > 10) {
    PD_SNP_list_remove <- PD_SNP_list[-index,]
    matrix_remove <- matrix[-index,-index]
    lambda_remove <- estimate_s_rss(PD_SNP_list_remove$z_scores, matrix_remove, n=38000)
    print(lambda_remove)
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(PD_SNP_list_remove$z_scores, matrix_remove, n = 38000, L = 5)
  } else {
    # run SuSiE with sum stats
    fitted_rss <- susie_rss(PD_SNP_list$z_scores, matrix, n = 38000, L = 5)
  }
  # record number of SNP
  for (i in 1:5) {
    CS95 = paste0("fitted_rss$sets$cs$L",i)
    nSNP <- nSNP + print(length(eval(parse(text = CS95))))
  }
}

print(nSNP)

# 3 locus (1805)

#### TEST
locus = 1
# original
plink_clumping_SNP <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_clumping/PD_clump/PD_locus_", locus, ".txt"), header = FALSE)
colnames(plink_clumping_SNP) <- c("SNP")
plink_clumping_matrix <- read.table(paste0("/Users/cl2375/OneDrive - Yale University/PhD study/Labs/Huckins & Cotsapas Lab/Fine-mapping/LD_calculation/PD_LD/LDmatrix_PD_", locus, ".ld"), header = FALSE)
PD_SNP_list <- PD_meta[PD_meta$SNP %in% plink_clumping_SNP$SNP, ]
PD_SNP_list$index <- seq(1:80)
matrix <- data.matrix(plink_clumping_matrix)
lambda <- estimate_s_rss(PD_SNP_list$z_scores, matrix, n=38000)
print(lambda)
condz_in <- kriging_rss(PD_SNP_list$z_scores, matrix, n=38000)
condz_in$plot
# adjusted
# positive
matrix_test_pos <- matrix[PD_SNP_list$z_scores > 0, PD_SNP_list$z_scores > 0]
PD_SNP_list_test_pos <- PD_SNP_list[PD_SNP_list$z_scores > 0,]
lambda <- estimate_s_rss(PD_SNP_list_test_pos$z_scores, matrix_test_pos, n=38000)
print(lambda)
condz_in <- kriging_rss(PD_SNP_list_test_pos$z_scores, matrix_test_pos, n=38000)
condz_in$plot
index_pos_sub <- which(abs(condz_in$conditional_dist$z_std_diff)>3)
index_pos <- PD_SNP_list_test_pos[which(abs(condz_in$conditional_dist$z_std_diff)>3), index]
#PD_SNP_list_test_pos_remove <- PD_SNP_list_test_pos[index_pos,]
PD_SNP_list_test_pos_remove <- PD_SNP_list_test_pos[which(abs(condz_in$conditional_dist$z_std_diff)>3),]
fitted_rss <- susie_rss(PD_SNP_list_test_pos$z_scores, matrix_test_pos, n = 38000, L = 5)
susie_plot(fitted_rss, y ="PIP")
# negative
matrix_test_neg <- matrix[PD_SNP_list$z_scores < 0, PD_SNP_list$z_scores < 0]
PD_SNP_list_test_neg <- PD_SNP_list[PD_SNP_list$z_scores < 0,]
lambda <- estimate_s_rss(PD_SNP_list_test_neg$z_scores, matrix_test_neg, n=38000)
print(lambda)
condz_in <- kriging_rss(PD_SNP_list_test_neg$z_scores, matrix_test_neg, n=38000)
condz_in$plot
index_neg_sub <- which(abs(condz_in$conditional_dist$z_std_diff)>3)
index_neg <- PD_SNP_list_test_neg[which(abs(condz_in$conditional_dist$z_std_diff)>3), index]
fitted_rss <- susie_rss(PD_SNP_list_test_neg$z_scores, matrix_test_neg, n = 38000, L = 5)
susie_plot(fitted_rss, y ="PIP")
# combine
# check the index for removal
index_remove <- PD_SNP_list$SNP %in% PD_SNP_list_test_pos_remove$SNP
PD_SNP_list_adjust <- PD_SNP_list[!(index_remove),]
matrix_adjust <- matrix[!(index_remove),!(index_remove)]
fitted_rss <- susie_rss(PD_SNP_list_adjust$z_scores, matrix_adjust, n = 38000, L = 5)

# plot LD
library(LDheatmap)
library(grid)
matrix_LDmap <- matrix*matrix
colnames(matrix_LDmap) <- PD_SNP_list$SNP
rownames(matrix_LDmap) <- PD_SNP_list$SNP
LDheatmap(matrix_LDmap)


# before fine-map
plot(-log(PD_SNP_list$P), ylab = "-logP", xlab="variants", xaxt="n", pch=19, col="black")
axis(side=1, at=c(1:nrow(PD_SNP_list)), labels=FALSE)
abline(h = -log(5e-08), col = "red")
text(x=1:nrow(PD_SNP_list),y = par("usr")[3]-1.1,labels=as.graphicsAnnot(PD_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)

susie_plot(fitted_rss, y="PIP", b=log(PD_SNP_list$OR), add_legend = "topleft", xlab="variants", main="SuSiE credible sets 1000G.EUR", xaxt="n")
axis(side=1, at=c(1:nrow(PD_SNP_list)), labels=FALSE)
text(x=1:nrow(PD_SNP_list),y = par("usr")[3]-0.01,labels=as.graphicsAnnot(PD_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)
#text(x=1:nrow(PD_SNP_list),y = par("usr")[3]-0.05,labels=as.graphicsAnnot(PD_SNP_list$SNP),xpd=NA,srt=45,cex=0.65,adj=1)
#----
