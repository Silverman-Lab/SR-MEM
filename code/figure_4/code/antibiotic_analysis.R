library(nlme)
library(ALDEx3)
library(dplyr)
library(tidyr)
library(RTMBdist)
library(ggplot2)

set.seed(654)

### Load Data
counts <- read.csv("../data/tblcounts_asv_melt.csv")
antibiotics <- read.csv("../data/tbldrug.csv")
meta <- read.csv("../data/tblASVsamples.csv")
qpcr <- read.csv("../data/tblqpcr.csv")

## Pre-processing antibiotics, qpcr, meta
main_antibiotics <- c("carbapenems", "cephalosporins", "glycopeptide antibiotics",
                      "macrolide derivatives", "oxazolidinone antibiotics",
                      "penicillins", "quinolones", "sulfonamides", "aztreonam",
                      "metronidazole")
other_antibiotics <- c("aminoglycosides", "antituberculosis agents",
                       "glycylcyclines", "leprostatics", "lincomycin derivatives",
                       "miscellaneous antibiotics", "tetracyclines")
antibiotics[antibiotics$Factor=="aztreonam","Category"] <- "aztreonam"
antibiotics[antibiotics$Factor=="metronidazole","Category"] <- "metronidazole"
antibiotics <- antibiotics[
   (antibiotics$Category%in%main_antibiotics)|(antibiotics$Category%in%other_antibiotics),]
antibiotics[antibiotics$Category%in%other_antibiotics,"Category"] <- "other_antibiotics"
antibiotics <- antibiotics[antibiotics$Route%in%c("oral", "intravenous"),]
antibiotics$cat_route <- paste(antibiotics$Category, antibiotics$Route, sep="_")
qpcr <- qpcr[qpcr$qPCR16S>0,]
row.names(qpcr) <- qpcr$SampleID
meta <- meta[meta$SampleID%in%qpcr$SampleID,]
meta <- meta[!duplicated(meta$SampleID),]
meta <- meta %>%
  group_by(PatientID) %>%
  filter(n() >= 10) %>%
  ungroup()

## Build Read Count Matrix
asvs <- c("ASV_8", "ASV_2", "ASV_51", "ASV_32", "ASV_16", "ASV_19", "ASV_6",
          "ASV_239", "ASV_128", "ASV_500", "ASV_34", "ASV_20", "ASV_81", "ASV_150",
          "ASV_39", "ASV_258", "ASV_29", "ASV_310", "ASV_635", "ASV_107")
Y_df <- counts %>% pivot_wider(
    names_from = SampleID,
    values_from = Count,
    values_fill = 0
  )
Y <- as.matrix(Y_df[,-1])
rownames(Y) <- Y_df[[1]]
other <- colSums(Y[!(row.names(Y)%in%asvs),])
Y <- Y[row.names(Y)%in%asvs,]
Y <- rbind(Y, other)

## Build Regression Metadata
all_antib <- unique(antibiotics$cat_route)
X <- c()
for(row_i in 1:nrow(meta)) {
  row <- meta[row_i,]
  antib <- antibiotics[(antibiotics$PatientID==row["PatientID"][[1]])
                       &(antibiotics$StartTimepoint<=as.numeric(row["Timepoint"][[1]]))
                       &(antibiotics$StopTimepoint>=as.numeric(row["Timepoint"][[1]])),]
  X <- rbind(X, all_antib%in%antib$cat_route)
}
X <- apply(X, c(1,2), as.numeric)
colnames(X) <- all_antib
row.names(X) <- meta$SampleID

## Filter Within Window
filt <- c()
for (rn in row.names(X)) {
  row <- X[rn,]
  if(sum(row)==0) {
    patient_id <- unique(meta[meta$SampleID==rn,"PatientID"])
    if(nrow(antibiotics[antibiotics$PatientID==patient_id,])==0) {
      filt <- c(filt, 0)
    } else {
      days <- unname(unlist(apply(antibiotics[antibiotics$PatientID==patient_id,],
                                  1, function(arow) {
            sort(unique(as.numeric(arow["StartTimepoint"][[1]]):
                  (as.numeric(arow["StopTimepoint"][[1]])+21)))})))
      time_point <- as.numeric(meta[meta$SampleID==rn,"Timepoint"])
      if(any(time_point%in%days)) {
        filt <- c(filt, 1)
      } else {
        filt <- c(filt, 0)
      }
    }
  } else {
    filt <- c(filt, 0)
  }
}
X <- X[filt==0,]
colnames(X) <- gsub(" ", "_", colnames(X))
Y <- Y[,row.names(X)]
qpcr <- qpcr[colnames(Y),]
data <- data.frame(X)
patient_ids <- unname(unlist(sapply(colnames(Y), function(item) {
  unique(meta[meta$SampleID==item,"PatientID"])
})))
time_points <- unname(unlist(sapply(colnames(Y), function(item) {
  as.numeric(unique(meta[meta$SampleID==item,"Timepoint"]))
})))
data$PatientID <- patient_ids
data$time_point <- time_points


##### MEM Modeling #####
sparsity <- replicate(1000, {
  pp <- t(rdirmult(length(colSums(Y)), colSums(Y), t(Y+0.5)))
  colSums(pp==0)
})
sparsity_intervals <- apply(
  sparsity, 1,
  quantile, probs = c(0.025, 0.5, 0.975)
)
df_plot <- data.frame(
  sample_name = colnames(Y==0),
  observed = colSums(Y==0),
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
    title = "Sparsity Posterior Predictive Antibiogram Dataset"
  ) +
  theme_bw()
ggsave("../../../supplement/SFigure_8.png", g, units="in", width=6, height=6)

## Run ALDEx3
formula <- as.formula(paste("~", paste(colnames(X), collapse=" + ")))
aldex.res <- aldex(Y, X=formula, scale=sample.sm, nsample=3000, n.cores=5, data=data,
      method="nlme", correlation=corAR1(form=~time_point|PatientID), random=~1|PatientID,
      s.mu=log2(qpcr$qPCR16S), s.var=rep(0.25, ncol(Y)))
p <- aldex.res$p.val.adj
est <- apply(aldex.res$estimate, c(1,2), mean)
write.csv(p, "../results/antibiotic_adj.p.vals.csv")
write.csv(est, "../results/antibiotic_est.csv")
