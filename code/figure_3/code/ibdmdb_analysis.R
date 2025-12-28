library(ALDEx3)
library(lmerTest)
library(MASS)
library(dplyr)
library(ggplot2)
library(RTMBdist)

set.seed(54637)

## Vandeputtee Analysis
metadata <- read.csv("../data/41586_2017_BFnature24460_MOESM10_ESM.csv")
row.names(metadata) <- metadata$Sample

taxonomy <- read.csv("../data/taxa_assignments.csv", row.names=1, header=T)
Y <- apply(t(read.csv("../data/OTU_nochim.csv", row.names=1)), c(1,2), as.numeric)
Y <- Y[row.names(Y)%in%row.names(taxonomy),]
Y <- Y[,colnames(Y)%in%metadata$Sample]
genus <- taxonomy[row.names(Y), "Genus"]
genus[is.na(genus)] <- "unclassified"
row.names(Y) <- genus
Y_genus <- rowsum(Y, group=rownames(Y))
other <- colSums(Y_genus[(rowSums(Y_genus==0)/ncol(Y_genus))>0.75,])
Y_genus <- Y_genus[(rowSums(Y_genus==0)/ncol(Y_genus))<=0.75,]
Y_genus <- rbind(Y_genus, other)

metadata <- metadata[colnames(Y_genus),]
scales <- log2(metadata$Average.cell.count..per.gram.of.frozen.feces.)
X <- data.frame(disease=metadata$Health.status)
X$disease <- factor(X$disease, levels=c("Control", "CD"))
scale_m <- replicate(1000, rnorm(length(scales), scales, 0.7))
vand.res <- aldex(Y_genus, X=~disease, data=X, nsample=1000, scale=scale_m)
est <- apply(vand.res$estimate, c(1,2), mean)["diseaseCD",]
padj <- vand.res$p.val.adj["diseaseCD",]
vand.res <- cbind(est, padj)
write.csv(vand.res, "../results/vand.res.csv")

#### IBDMDB Analysis ####
Y_idb <- read.csv("../data/ibddmb.counts.csv", row.names=1)
dysbiosis <- read.table("../data/dysbiosis_scores.tsv", row.names=1)
ibd_meta <- read.csv("../data/hmp2_metadata_2018-08-20.csv")
ibd_meta <- ibd_meta[ibd_meta$data_type=="metagenomics",]
row.names(ibd_meta) <- ibd_meta$External.ID
ibd_meta <- ibd_meta[row.names(ibd_meta)%in%row.names(dysbiosis),]
ibd_meta <- ibd_meta[,c("diagnosis", "Age.at.diagnosis",
                        "Antibiotics", "Participant.ID",
                        "site_name")]
ibd_meta <- ibd_meta[row.names(ibd_meta)%in%colnames(Y_idb),]
Y_idb <- Y_idb[,row.names(ibd_meta)]
other <- colSums(Y_idb[(rowSums(Y_idb==0)/ncol(Y_idb))>0.75,])
Y_idb <- Y_idb[(rowSums(Y_idb==0)/ncol(Y_idb))<=0.75,]
Y_idb <- rbind(Y_idb, other)
ibd_meta$dysbiosis <- ifelse(dysbiosis[colnames(Y_idb),2], "dys", "non_dys")
ibd_meta$dys_disease <- paste(ibd_meta$diagnosis, ibd_meta$dysbiosis, sep="_")
ibd_meta$dys_disease <- factor(ibd_meta$dys_disease,
                               levels=c("nonIBD_non_dys", "CD_non_dys",
                                        "nonIBD_dys", "CD_dys",
                                        "UC_non_dys", "UC_dys"))
sparsity <- round(rowSums((Y_idb==0)/(ncol(Y_idb)))*100,1)
log2_mean <- round(apply(log2(Y_idb+0.5), 1, mean),1)
count_data <- data.frame(sparsity=sparsity, log2_mean=log2_mean)
write.csv(count_data, "../results/count.stats.csv")


##### MEM Modeling #####
sparsity <- replicate(1000, {
  pp <- t(rdirmult(length(colSums(Y_idb)), colSums(Y_idb), t(Y_idb+0.5)))
  colSums(pp==0)
})
sparsity_intervals <- apply(
  sparsity, 1,
  quantile, probs = c(0.025, 0.5, 0.975)
)
df_plot <- data.frame(
  sample_name = colnames(Y_idb==0),
  observed = colSums(Y_idb==0),
  q025 = sparsity_intervals[1, ],
  q50  = sparsity_intervals[2, ],
  q975 = sparsity_intervals[3, ]
)
df_plot <- df_plot[sample(1:nrow(df_plot), 50), ]
g <- ggplot(df_plot, aes(y = sample_name)) +
  # predictive interval
  geom_errorbar(
    aes(xmin = q025, xmax = q975),
    width = 0.2,
    color = "gray40"
  ) +
  # predictive median
  geom_point(
    aes(x = q50),
    size = 2,
    color = "black"
  ) +
  # observed value
  geom_point(
    aes(x = observed),
    color = "red",
    size = 2
  ) +
  labs(
    y = "Sample Name",
    x = "Number of zero counts",
    title = "Sparsity Posterior Predictive IBDMDB Dataset"
  ) +
  theme_bw()
ggsave("../../../supplement/SFigure_7.png", g, units="in", width=6, height=6)

## Original
P_idb <- apply(Y_idb, 2, function(col) col/sum(col))
ARC <- apply(P_idb, c(1,2), function(item) asin(sqrt(item)))
arc.res.CD <- c()
arc.res.UC <- c()
for(i in 1:nrow(ARC)) {
  row <- ARC[i,]
  A <- row
  data <- cbind(A, ibd_meta)
  one.arc <- coef(summary(lmer(
    A~dys_disease++Antibiotics+(1|Participant.ID)+(1|site_name),
    data=data)))
  arc.res.CD <- rbind(arc.res.CD, one.arc["dys_diseaseCD_dys",])
  arc.res.UC <- rbind(arc.res.UC, one.arc["dys_diseaseUC_dys",])
}
padj <- p.adjust(arc.res.CD[, 5], method="BH")
arc.res.CD <- cbind(arc.res.CD, padj)
row.names(arc.res.CD) <- row.names(ARC)
write.csv(arc.res.CD, "../results/original.res.CD.csv")

padj <- p.adjust(arc.res.UC[, 5], method="BH")
arc.res.UC <- cbind(arc.res.UC, padj)
row.names(arc.res.UC) <- row.names(ARC)
write.csv(arc.res.UC, "../results/original.res.UC.csv")

## CLR Estimate Scale
P_idb <- apply(Y_idb, 2, function(col) (col+0.5)/sum(col+0.5))
clr_scales <- apply(P_idb, 2, function(col) -mean(log2(col)))
data <- cbind(clr_scales, ibd_meta)
coef(summary(lmer(
  clr_scales~dys_disease++Antibiotics+(1|Participant.ID)+(1|site_name),
  data=data)))

## SR-MEM
control_sarr <- log2(read.csv("../data/control.csv", header=F)[,2])
cd_sarr <- log2(read.csv("../data/crohn.csv", header=F)[,2])
uc_sarr <- log2(read.csv("../data/ulcerative.csv", header=F)[,2])

cd_diff <- mean(cd_sarr)-mean(control_sarr)
uc_diff <- mean(uc_sarr)-mean(control_sarr)

scale.m <- c()
for(i in 1:1000) {
  scales <- rep(0, nrow(ibd_meta))
  scales[ibd_meta$dys_disease=="CD_dys"] <- rnorm(1, cd_diff, 0.5)
  scales[ibd_meta$dys_disease=="UC_dys"] <- rnorm(1, uc_diff, 0.5)
  scales[!(ibd_meta$dys_disease%in%c("CD_dys", "UC_dys"))] <- 0
  scale.m <- cbind(scale.m, scales)
}

sr.mem.res <- aldex(Y_idb, n.cores=6,
      X=~dys_disease+Antibiotics+(1|Participant.ID)+(1|site_name),
      data=ibd_meta, scale=scale.m, nsample=1000,
      method="lme4")
est <- apply(sr.mem.res$estimate, c(1,2), mean)["dys_diseaseCD_dys",]
padj <- sr.mem.res$p.val.adj["dys_diseaseCD_dys",]
sr.mem.res.f <- cbind(est, padj)
write.csv(sr.mem.res.f, "../results/sr.mem.res.CD.csv")

est <- apply(sr.mem.res$estimate, c(1,2), mean)["dys_diseaseUC_dys",]
padj <- sr.mem.res$p.val.adj["dys_diseaseUC_dys",]
sr.mem.res.f <- cbind(est, padj)
write.csv(sr.mem.res.f, "../results/sr.mem.res.UC.csv")

## SR-MEM Bias Sensitivity
for(b in seq(-5, 5, 1)) {
  scale.m <- c()
  for(i in 1:1000) {
    scales <- rep(0, nrow(ibd_meta))
    scales[ibd_meta$dys_disease=="CD_dys"] <- rnorm(1, cd_diff, 0.5) +
      rnorm(sum(ibd_meta$dys_disease=="CD_dys"), 0, 0.7) + b
    scales[ibd_meta$dys_disease=="UC_dys"] <- rnorm(1, uc_diff, 0.5) +
      rnorm(sum(ibd_meta$dys_disease=="UC_dys"), 0, 0.7) + b
    scales[!(ibd_meta$dys_disease%in%c("CD_dys", "UC_dys"))] <-
      rnorm(sum(!(ibd_meta$dys_disease%in%c("CD_dys", "UC_dys"))), 0, 0.7)
    scale.m <- cbind(scale.m, scales)
  }

  sr.mem.res <- aldex(Y_idb, n.cores=6,
                      X=~dys_disease+Antibiotics+(1|Participant.ID)+(1|site_name),
                      data=ibd_meta, scale=scale.m, nsample=1000,
                      method="lme4")
  est <- apply(sr.mem.res$estimate, c(1,2), mean)["dys_diseaseCD_dys",]
  padj <- sr.mem.res$p.val.adj["dys_diseaseCD_dys",]
  sr.mem.res <- cbind(est, padj)
  write.csv(sr.mem.res, paste0("../results/sr.mem.res.bias.", b, ".csv"))
}

## CLR-MEM
clr.mem.res <- aldex(Y_idb, n.cores=6,
      X=~dys_disease+Antibiotics+(1|Participant.ID)+(1|site_name),
      data=ibd_meta, scale=clr.sm, nsample=1000,
      method="lme4")
est <- apply(clr.mem.res$estimate, c(1,2), mean)["dys_diseaseCD_dys",]
padj <- clr.mem.res$p.val.adj["dys_diseaseCD_dys",]
clr.mem.res.f <- cbind(est, padj)
write.csv(clr.mem.res.f, "../results/clr.mem.res.CD.csv")

est <- apply(clr.mem.res$estimate, c(1,2), mean)["dys_diseaseUC_dys",]
padj <- clr.mem.res$p.val.adj["dys_diseaseUC_dys",]
clr.mem.res.f <- cbind(est, padj)
write.csv(clr.mem.res.f, "../results/clr.mem.res.UC.csv")
