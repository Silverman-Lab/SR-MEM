library(ALDEx3)
library(ANCOMBC)
library(LinDA)
library(NBZIMM)
library(ALDEx3)
library(lmerSeq)
library(Maaslin2)
library(phyloseq)
library(reader)
library(DESeq2)

run_ancom_bc2 <- function(phyloseq_obj, fix_formula, rand_formula) {
    output <- ancombc2(data=phyloseq_obj,
                       fix_formula=fix_formula, rand_formula=rand_formula,
                       p_adj_method="BH", n_cl=10)
    return(output)
}

run_linda <- function(Y, X, formula, n.cores=2) {
    fit.linda <- linda(Y, X, formula,
                       imputation=TRUE,
                       n.cores=n.cores)
    return(fit.linda)
}

run_NBZIMM <- function(Y, X, fixed_formula, random_formula) {
    nbzimm.fit <- mms(y=t(Y), fixed=fixed_formula,
                      data=X, random=random_formula,
                      min.p=0.2, method="zinb")
    return(data.frame(NBZIMM::fixed(nbzimm.fit)$dist))
}

run_lmerseq <- function(Y, X, formula) {
  dds <- DESeqDataSetFromMatrix(countData=Y, colData=X, design=~disease)
  dds <- estimateSizeFactors(dds, type="poscounts")
  vst_expr <- assay(vst(dds, nsub=sum(rowMeans(counts(dds, normalized=TRUE))>5)))
  fit.lmerSeq <- lmerSeq.fit(form=formula,
                          expr_mat=vst_expr,
                          sample_data=X,
                          REML=T)
  lmerseq.res <- lmerSeq.summary(lmerSeq_results=fit.lmerSeq,
                                 coefficient="disease",
                                 p_adj_method = 'BH',
                                 ddf = 'Satterthwaite',
                                 sort_results=F)
  return(lmerseq.res$summary_table)
}

run_maaslin2 <- function(Y, X, fixed_effects, random_effects, normalization,
                         n.cores=2) {
    fit.maaslin2 <- Maaslin2(input_data=Y,
                             input_metadata=X,
                             output="massline_output",
                             fixed_effects=fixed_effects,
                             random_effects=random_effects,
                             normalization=normalization,
                             cores=n.cores,
                             plot_scatter=F,
                             plot_heatmap=F)
    return(fit.maaslin2)
}

args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]
base_file <- tools::file_path_sans_ext(basename(file_name))
sim_data <- readRDS(file_name)
out_dir <- args[2]
tp <- 1.3
tp_s <- 0.25
n.cores <- 2

Y <- sim_data$Y
X <- sim_data$X
D <- nrow(Y)
N <- ncol(Y)
other <- colSums(Y[rowSums(Y==0)>(0.7*N),])
seq.trim <- Y[rowSums(Y==0)<=(0.7*N),]
seq.trim <- rbind(seq.trim, other)

## ALDEx3
tryCatch({
  aldex.res <- aldex(seq.trim, X=~disease, data=X, method="lm",
                     nsample=300, scale=coef.sm, c.mu=c(0, tp),
                     c.cor=rbind(c(0, 0), c(0, tp_s)))
  p <- aldex.res$p.val.adj["disease",]
  est <- apply(aldex.res$estimate, c(1,2), mean)
  est <- est["disease",]

  m <- cbind(p, est)
  row.names(m) <- colnames(aldex.res$p.val.adj)
  colnames(m) <- c("p.val.adj", "lfc")
  saveRDS(m, cat.path(out_dir, paste0(base_file, "_aldex_lin.RDS")))
}, error=function(x){print(x);})

## lmerSeq
tryCatch({
  lmerseq.res <- run_lmerseq(seq.trim, X, ~disease+(1|subject))
  saveRDS(lmerseq.res, cat.path(out_dir, paste0(base_file, "_lmerseq.RDS")))
}, error=function(x){print(x);})

## Ancombc2
tryCatch({
  df <- data.frame(sim_data$X[,c("disease", "subject")])
  df[,"disease"] <- factor(df[,"disease"], levels=c(0, 1))
  df[,"subject"] <- factor(df[,"subject"])
  row.names(df) <- colnames(seq.trim)
  otu_table <- otu_table(seq.trim, taxa_are_rows=TRUE)
  samples <- sample_data(df)
  phyloseq_obj <- phyloseq(otu_table, NULL, samples)
  ancombc2_res <- run_ancom_bc2(phyloseq_obj, "disease", "(1|subject)")
  ancombc2_res_df <- ancombc2_res$res
  saveRDS(ancombc2_res_df, cat.path(out_dir, paste0(base_file, "_ancombc2.RDS")))
}, error=function(x){print(x);})

## NBZIMM
tryCatch({
  df <- data.frame(sim_data$X[,c("disease", "subject")])
  df[,"disease"] <- factor(df[,"disease"], levels=c(0, 1))
  df[,"subject"] <- factor(df[,"subject"])
  row.names(df) <- colnames(seq.trim)
  Nd <- colSums(seq.trim)
  Xs <- cbind(df, Nd)
  fixed_formula <- ~disease+offset(log(Nd))
  random_formula <- ~1|subject
  nbzimm.res <- run_NBZIMM(seq.trim, Xs, fixed_formula, random_formula)
  saveRDS(nbzimm.res, cat.path(out_dir, paste0(base_file, "_NBZIMM.RDS")))
}, error=function(x){print(x);})

## Maaslin2
tryCatch({
  df <- data.frame(sim_data$X[,c("disease", "subject")])
  df[,"disease"] <- factor(df[,"disease"], levels=c(0, 1))
  df[,"subject"] <- factor(df[,"subject"])
  row.names(df) <- colnames(seq.trim)
  maaslin2.res <- run_maaslin2(seq.trim, df, fixed_effects=c("disease"),
                               random_effects="subject", normalization="TSS",
                               n.cores=n.cores) ## TSS Default
  saveRDS(maaslin2.res, cat.path(out_dir, paste0(base_file, "_maaslin2.RDS")))
}, error=function(x){print(x);})

## LINDA
tryCatch({
  df <- data.frame(sim_data$X[,c("disease", "subject")])
  df[,"disease"] <- factor(df[,"disease"], levels=c(0, 1))
  df[,"subject"] <- factor(df[,"subject"])
  row.names(df) <- colnames(seq.trim)
  formula <- "~disease+(1|subject)"
  linda.res <- run_linda(seq.trim, df, formula)
  saveRDS(linda.res, cat.path(out_dir, paste0(base_file, "_linda.RDS")))
}, error=function(x){print(x);})

## SR-MEM
tryCatch({
  aldex.res <- aldex(seq.trim, X=~disease+(1|subject), data=X, method="lme4",
                     nsample=300, scale=coef.sm, c.mu=c(0, tp),
                     c.cor=rbind(c(0, 0), c(0, tp_s)), n.cores=n.cores)
  p <- aldex.res$p.val.adj["disease",]
  est <- apply(aldex.res$estimate, c(1,2), mean)
  est <- est["disease",]

  m <- cbind(p, est)
  row.names(m) <- colnames(aldex.res$p.val.adj)
  colnames(m) <- c("p.val.adj", "lfc")
  saveRDS(m, cat.path(out_dir, paste0(base_file, "_aldex.RDS")))
}, error=function(x){print(x);})
