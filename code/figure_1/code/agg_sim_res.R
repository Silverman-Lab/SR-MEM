library(dplyr)

all_res <- c()
sim_mat <- rbind(c("oral", 0.1, 0.4, -2, 3, 1, 2), c("oral", 0.1, 0.7, -2, 3, 1, 2),
                 c("armpit", 0, 0.5, 0, 3, 1, 2), c("gut", 0.1, 0.7, -2, 3, 1, 2),
                 c("oral", 0, 0.2, 0, 3, 1, 2), c("gut", 0, 0.2, 0, 3, 1, 2),
                 c("gut", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5),
                 c("oral", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5))
for(n_s in c(6, 8, 12, 20, 40, 80, 120, 160, 200, 250, 300)) {
	print(n_s)
for(row in 1:nrow(sim_mat)) {
  sim_row <- sim_mat[row,]
  for(method in c("ancombc2", "linda", "maaslin2", "maaslin3", "NBZIMM", "lmerseq", "aldex", "aldex_lin")) {
    for(ind in 1:20) {
      e <- tryCatch({
        x <- readRDS(paste0("../raw_data/sim_inputs/", sim_row[1], "_", n_s, "_20_1_", sim_row[2], "_", sim_row[3], "_", sim_row[4], "_", sim_row[5], "_", ind, ".RDS"))
        ref_lfcs <- x$all_lfcs
        output <- readRDS(paste0("../raw_data/sim_outputs/", sim_row[1], "_", n_s, "_20_1_", sim_row[2], "_", sim_row[3], "_", sim_row[4], "_", sim_row[5], "_", ind, "_", method, ".RDS"))
      }, error=function(x){;})
      if(is.null(e)) {next}

      all_lfcs <- rep(NA, length(ref_lfcs))
      names(all_lfcs) <- names(ref_lfcs)
      all_p <- rep(1, length(ref_lfcs))
      names(all_p) <- names(ref_lfcs)

      if(method=="ancombc2") {
        p <- output$q_disease1
        p[!output$passed_ss_disease1] <- 1
        names(p) <- output$taxon
        lfcs <- output$lfc_disease1
        names(lfcs) <- output$taxon
      } else if (method=="aldex") {
        p <- output[,1]
        names(p) <- row.names(output)
        lfcs <- output[,2]
        names(lfcs) <- row.names(output)
      } else if (method=="aldex_lin") {
        p <- output[,1]
        names(p) <- row.names(output)
        lfcs <- output[,2]
        names(lfcs) <- row.names(output)
      } else if (method=="linda") {
        p <- output$output$disease$padj
        names(p) <- row.names(output$output$disease)
        lfcs <- output$output$disease$log2FoldChange
        names(lfcs) <- row.names(output$output$disease)
      } else if (method=="maaslin2") {
        p <- output$results$qval
        names(p) <- output$results$feature
        lfcs <- output$results$coef
        names(lfcs) <- output$results$feature
      } else if(method=="NBZIMM") {
        output <- output[output$variables=="disease1",]
        p <- output$padj
        names(p) <- output$responses
        lfcs <- output$Estimate
        names(lfcs) <- output$responses
      } else if(method=="lmerseq") {
        p <- output$p_val_adj
        names(p) <- output$gene
        lfcs <- output$Estimate
        names(lfcs) <- output$gene
      } else if(method=="maaslin3") {
        p <- output$pval_joint
        p <- p.adjust(p, method="BH")
        names(p) <- output$feature
        lfcs <- output$coef
        names(lfcs) <- output$feature
      }
      lfcs <- lfcs[names(lfcs)%in%names(ref_lfcs)]
      p <- p[names(p)%in%names(ref_lfcs)]
      all_p[names(p)] <- p
      all_lfcs[names(lfcs)] <- lfcs
      res <- cbind(all_p, all_lfcs, ref_lfcs)
      NTP <- sum((res[,1]<=0.05)&(res[,3]!=0)&(sign(res[,2])!=sign(res[,3])))
      TP <- sum((res[,1]<=0.05)&(res[,3]!=0)&(sign(res[,2])==sign(res[,3])))
      FP <- sum((res[,1]<=0.05)&(res[,3]==0))
      FN <- sum((res[,1]>0.05)&(res[,3]!=0))
      TN <- sum((res[,1]>0.05)&(res[,3]==0))
      if (is.na(TP)|is.na(FP)) {
        next
      } else if ((TP+FP)!=0) {
        FDR <- FP/(TP+FP)
      } else {
        FDR <- 0
      }
      all_res <- rbind(all_res, c(method, n_s, ind, sim_row[1], sim_row[2], sim_row[3],
                                  sim_row[4], sim_row[5],
                                  TP, NTP, FP, FN, TN, FDR, TP/(TP+FN+NTP)))
    }
  }
}
}

all_res <- data.frame(all_res)
colnames(all_res) <- c("method", "n_s", "i", "system", "pc_l", "pc_u", "lfc_l", "lfc_u",
                       "TP", "NTP", "FP", "FN", "TN", "FDR", "Power")
all_res$i <- as.numeric(all_res$i)
all_res$pc_l <- as.numeric(all_res$pc_l)
all_res$pc_u <- as.numeric(all_res$pc_u)
all_res$lfc_l <- as.numeric(all_res$lfc_l)
all_res$lfc_u <- as.numeric(all_res$lfc_u)
all_res$TP <- as.numeric(all_res$TP)
all_res$NTP <- as.numeric(all_res$NTP)
all_res$FP <- as.numeric(all_res$FP)
all_res$FN <- as.numeric(all_res$FN)
all_res$TN <- as.numeric(all_res$TN)
all_res$FDR <- as.numeric(all_res$FDR)
all_res$Power <- as.numeric(all_res$Power)
all_res <- all_res %>% group_by(n_s, method, system, pc_l, pc_u, lfc_l, lfc_u) %>% summarize(FDR_mean=mean(FDR),
                                                       Power_mean=mean(Power))
data.frame(all_res)
write.csv(all_res, "../results/simulation_results.csv")
