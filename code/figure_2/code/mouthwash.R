library(ALDEx3)
library(biomformat)
library(phyloseq)
library(nlme)
library(LinDA)
library(DESeq2)
library(ANCOMBC, lib.loc="~/local_R/")
library(Maaslin2)
library(lmerSeq)
library(NBZIMM)
library(dplyr)
library(tidyr)
library(lmerTest)
library(xlsx)

## Setup
set.seed(564527)
mc.samples <- 1000

##### Function Definitions #####

## Taken from Vandeputtee QMP https://github.com/raeslab/QMP
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) {
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE)
        stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table)
  cell_counts_table = t(cell_counts_table[row.names(cnv_corrected_abundance_table),])
  sample_sizes = rowSums(cnv_corrected_abundance_table)
  sampling_depths = sample_sizes / cell_counts_table
  minimum_sampling_depth = min(sampling_depths)
  rarefy_to = cell_counts_table * minimum_sampling_depth
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table,
                                                     taxa_are_rows = FALSE)
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq),
                         ncol = ncol(cnv_corrected_abundance_table_phyloseq),
                         dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq),
                                         colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,],
                           sample.size = rarefy_to[i], rngseed = 711,
                           replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}

run_ancombc2 <- function(phyloseq_obj, fix_formula, rand_formula) {
    output <- ancombc2(data=phyloseq_obj,
                       fix_formula=fix_formula, rand_formula=rand_formula,
                       p_adj_method="BH")
    return(output)
}

run_linda <- function(Y, X, formula, n.cores=5) {
    fit.linda <- linda(Y, X, formula,
                       n.cores=n.cores)
    return(fit.linda)
}

run_NBZIMM <- function(Y, X, fixed_formula, random_formula) {
    nbzimm.fit <- mms(y=t(Y), fixed=fixed_formula,
                      data=X, random=random_formula,
                      min.p=0.2, method="zinb")
    return(nbzimm.fit)
}

run_lmerseq <- function(Y, X, formula) {
  dds <- DESeqDataSetFromMatrix(countData=Y, colData=X, design=~treat)
  dds <- estimateSizeFactors(dds, type="poscounts")
  vst_expr <- assay(vst(dds, nsub=sum(rowMeans(counts(dds, normalized=TRUE))>5)))
  fit.lmerSeq <- lmerSeq.fit(form=formula,
                          expr_mat=vst_expr,
                          sample_data=X,
                          REML=T)
  lmerseq.res <- lmerSeq.summary(lmerSeq_results=fit.lmerSeq,
                                 coefficient="treatalcoholfree-mouthwash:timec15min-later",
                                 p_adj_method = 'BH',
                                 ddf = 'Satterthwaite',
                                 sort_results=F)
  return(lmerseq.res$summary_table)
}

###### Sequence Count Data Processing #####

## Get sequence count data
biom_data <- read_hdf5_biom("../data/T4_SRS_final_QF_table_wTaxa.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$metadata$Taxon})

## Collapse to Genus
Y_genus <- c()
genus <- unname(sapply(row.names(Y), function(x) {strsplit(x, split=";")[[1]][6]}))
unique_genus <- unique(genus[!is.na(genus)])
all_inds <- c()
for(genus_n in unique_genus) {
    genus_inds <- which(genus%in%genus_n)
    all_inds <- c(all_inds, genus_inds)
    if(length(genus_inds)>1) {
        new_row <- colSums(Y[genus_inds,])
    } else {
        new_row <- Y[genus_inds,]
    }
    Y_genus <- rbind(Y_genus, new_row)
}
row.names(Y_genus) <- unique_genus

##### Metadata Processing ######

## Read and process metadata, get only living cells (PMA)
metadata <- data.frame(read.csv("../data/T3_SRS_metadata_ms.txt", sep="\t"))
meta_trim <- metadata[metadata$processing=="PMA",]
meta_trim$timec <- factor(meta_trim$collection_description,
                          levels=c("before", "15min-later", "2hrs-later"))
meta_trim$treat <- factor(meta_trim$Treatment_description,
                          levels=c("water-mouthwash", "antiseptic-mouthwash",
                              "alcoholfree-mouthwash", "soda"))
meta_trim$log_scale <- log2(meta_trim$FC_avg_cells_per_ul)
scale_sd <- apply(log2(
  meta_trim[,c("FC_cells_per_ul_r1", "FC_cells_per_ul_r2")]), 1, sd)
meta_trim$log_scale_sd <- scale_sd

## Make sure sample names all line up
meta_trim_filt <- meta_trim[meta_trim$X.SampleID%in%colnames(Y_genus),]
row.names(meta_trim_filt) <- meta_trim_filt$X.SampleID
Y_genus <- Y_genus[,row.names(meta_trim_filt)]

## Amalgamate low abundance genera 70\% zeros
keep_names <- row.names(Y_genus[((rowSums(Y_genus==0))/ncol(Y_genus))<0.7,])
other <- colSums(Y_genus[((rowSums(Y_genus==0))/ncol(Y_genus))>=0.7,])
Y_genus <- Y_genus[keep_names,]
Y_genus <- rbind(Y_genus, other)

## Interaction Term Between Time and Treatment
meta_trim_filt$treat_timec <- interaction(meta_trim_filt$treat,
                                          meta_trim_filt$timec, drop = TRUE)

##### MEM Modeling #####

## SR-MEM
scale_regression <- lme(log_scale~treat*timec,
                        ~1|participant_id,
                        data=meta_trim_filt,
                        weights=varIdent(form=~1|treat_timec))
scale_coefs <- coef(summary(scale_regression))
scale_cov_matrix <- vcov(scale_regression)
aldex.res <- aldex(Y_genus, X=~treat*timec+(1|participant_id),
                   data=data.frame(meta_trim_filt), method="lme4",
                   nsample=2000, scale=coef.sm, c.mu=scale_coefs[,1],
                   c.cor=scale_cov_matrix, n.cores=10)

## Linda
linda.res <- run_linda(Y_genus, meta_trim_filt, "~treat*timec+(1|participant_id)")

## lmerSeq
lmerseq.res <- run_lmerseq(Y_genus, meta_trim_filt, ~treat*timec+(1|participant_id))

## ANCOM-BC2
otu_table <- otu_table(Y_genus, taxa_are_rows=TRUE)
samples <- sample_data(meta_trim_filt)
phyloseq_obj <- phyloseq(otu_table, NULL, samples)
ancom_bc2_res <- run_ancombc2(phyloseq_obj, "treat*timec", "(1|participant_id)")

## NBZIMM
Nd <- colSums(Y_genus)
Xs <- cbind(meta_trim_filt, Nd)
fixed_formula <- ~treat*timec+offset(log(Nd))
random_formula <- ~1|participant_id
nbzimm.res <- run_NBZIMM(Y_genus, Xs, fixed_formula, random_formula)
fixed.nbzimm.res <- data.frame(fixed(nbzimm.res)$dist)
fixed.nbzimm.res$id <- row.names(fixed.nbzimm.res)
fixed.nbzimm.trans <- fixed.nbzimm.res %>%
  separate(id, into = c("group", "replicate"), sep = "--")

## QMP w/ pseudo-count
fc_live_cells <- cbind(meta_trim_filt$FC_avg_cells_per_ul)
row.names(fc_live_cells) <- row.names(meta_trim_filt)
q <- t(rarefy_even_sampling_depth(t(Y_genus), fc_live_cells))
A <- log2(q+0.5)
qmp.res <- data.frame()
for(genus in 1:nrow(A)) {
  d <- meta_trim_filt
  d$A <- A[genus,]
  qmp_res <- lmer(A~treat*timec+(1|participant_id), data=d)
  qmp_df <- coef(summary(qmp_res))[,c("Estimate", "Pr(>|t|)")]
  colnames(qmp_df) <- c("Est", "p.val")
  qmp_df <- data.frame(row = rownames(qmp_df), qmp_df)
  qmp_df$genus <- row.names(A)[genus]
  qmp.res <- rbind(qmp.res, qmp_df)
}
qmp.res <- qmp.res %>%
  group_by(row) %>%
  mutate(pval_adj=p.adjust(p.val, method="BH")) %>%
  ungroup()
qmp.res <- data.frame(qmp.res)

fn <- "../../../supplement/Fig2_Results.xlsx"
## Comparison Data
if (file.exists(fn)) {
  #Delete file if it exists
  file.remove(fn)
}
cols <- c("treatantiseptic-mouthwash", "treatalcoholfree-mouthwash",
          "treatsoda", "timec15min-later",
          "timec2hrs-later", "treatantiseptic-mouthwash:timec15min-later",
          "treatalcoholfree-mouthwash:timec15min-later",
          "treatsoda:timec15min-later", "treatantiseptic-mouthwash:timec2hrs-later",
          "treatalcoholfree-mouthwash:timec2hrs-later",
          "treatsoda:timec2hrs-later")
df <- data.frame(cbind(paste0("fixed_effect_", 1:length(cols)), cols))
colnames(df) <- c("Sheet", "Fixed Effect Name")
write.xlsx(df, file=fn, append=TRUE, sheetName="legend")
counter <- 1
for(col in cols) {
  col2 <- gsub("-", ".", gsub(":", ".", col))
  col3 <- gsub("-", "_", gsub(":", "_", col))
  print(c(counter, col3))
  ## Create CSVs
  aldex.sig <- cbind(apply(aldex.res$estimate, c(1,2), mean)[col,],
                    aldex.res$p.val.adj[col,])
  colnames(aldex.sig) <- c(paste0(col, ":est"), paste0(col, ":pval.adj"))
  all_res <- data.frame(aldex.sig)
  all_res$method <- "aldex_mme"
  all_res$genus <- row.names(aldex.sig)
  colnames(all_res) <- c("est", "padj", "method", "genus")

  linda.df <- data.frame(linda.res$output[col][[1]])[,c("log2FoldChange", "padj")]
  linda.df$method <- "linda"
  linda.df$genus <- row.names(linda.df)
  colnames(linda.df) <- c("est", "padj", "method", "genus")
  all_res <- rbind(all_res, linda.df)

  lmer.df <- data.frame(lmerseq.res)
  row.names(lmer.df) <- lmer.df$gene
  lmer.df$method <- "lmerSeq"
  lmer.df$genus <- lmer.df$gene
  lmer.df <- lmer.df[,c("Estimate", "p_val_adj", "method", "genus")]
  colnames(lmer.df) <- c("est", "padj", "method", "genus")
  all_res <- rbind(all_res, lmer.df)

  ancom_bc2.df <- data.frame(ancom_bc2_res$res)
  ancom_bc2.df <- ancom_bc2.df[,c("taxon", paste0("lfc_", col2), paste0("q_", col2), paste0("passed_ss_", col2))]
  row.names(ancom_bc2.df) <- ancom_bc2.df$taxon
  ancom_bc2.df[ancom_bc2.df[,paste0("passed_ss_", col2)]==FALSE, paste0("q_", col2)] <- 1
  ancom_bc2.df$method <- "ancom_bc2"
  ancom_bc2.df$genus <- row.names(ancom_bc2.df)
  ancom_bc2.df <- ancom_bc2.df[,c(paste0("lfc_", col2), paste0("q_", col2), "method", "genus")]
  colnames(ancom_bc2.df) <- c("est", "padj", "method", "genus")
  all_res <- rbind(all_res, ancom_bc2.df)

  nbzimm.df <- data.frame(fixed.nbzimm.trans[fixed.nbzimm.trans[,"variables"]==col,])
  row.names(nbzimm.df) <- nbzimm.df$responses
  nbzimm.df <- nbzimm.df[,c("Estimate", "padj")]
  nbzimm.df$method <- "nbzimm"
  nbzimm.df$genus <- row.names(nbzimm.df)
  colnames(nbzimm.df) <- c("est", "padj", "method", "genus")
  all_res <- rbind(all_res, nbzimm.df)

  qmp.df <- qmp.res[qmp.res$row==col,]
  qmp.df <- data.frame(qmp.df[,c("Est", "pval_adj", "genus")])
  row.names(qmp.df) <- qmp.df$genus
  qmp.df <- qmp.df[,c(1,2)]
  qmp.df$method <- "qmp"
  qmp.df$genus <- row.names(qmp.df)
  colnames(qmp.df) <- c("est", "padj", "method", "genus")
  all_res <- rbind(all_res, qmp.df)
  row.names(all_res) <- NULL

  # Find all unique combinations of a and b
  all_combinations <- expand.grid(
    method = unique(all_res$method),
    genus = unique(all_res$genus)
  )
  ## Left join to find missing rows
  df_full <- all_combinations %>%
  left_join(all_res, by = c("method", "genus")) %>%
  mutate(
    est = ifelse(is.na(est), 0, est),
    padj = ifelse(is.na(padj), 1, padj)
  )
  df_full <- df_full %>% mutate(method=dplyr::recode(method, aldex_mme='ALDEx2-MEM', linda='LinDA',
                                              qmp='QMP', nbzimm="NBZIMM", ancom_bc2="ANCOM-BC2",
                                              lmerSeq='lmerSeq', aldex_mme_tp="ALDEx2-MEM-TP"))
  df_full$method <- factor(df_full$method, levels=c("ALDEx2-MEM", "QMP", "LinDA", "lmerSeq",
                                                    "NBZIMM", "ANCOM-BC2", "ALDEx2-MEM-TP"))

  df_full <- df_full[order(df_full$method, df_full$genus),]
  sparsity <- rowSums((Y_genus==0)/(ncol(Y_genus)))[df_full$genus]
  log2_mean <- log2(apply(Y_genus, 1, mean))[df_full$genus]
  df_full$sparsity <- round(sparsity*100)
  df_full$log2_mean <- round(log2_mean)
  if(col3=="treatalcoholfree_mouthwash_timec15min_later") {
    write.csv(df_full, "../results/treatalcoholfree_mouthwash_timec15min_later.csv")
  }
  write.xlsx(df_full, file=fn, append=TRUE,
             sheetName=paste0("fixed_effect_", counter))
  counter <- counter + 1
}
