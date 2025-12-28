library(reader)
library(maaslin3)

run_maaslin3 <- function(Y, X, formula, normalization, n.cores=2) {
	fit.maaslin3 <- maaslin3(input_data=Y,
				 input_metadata=X,
	  			 formula=formula,
				 median_comparison_abundance=T,
				 median_comparison_prevalence=F,
				 standardize=F,
				 max_pngs=0,
				 output="masslin3",
				 plot_summary_plot=F,
				 plot_associations=F,
				 normalization=normalization,
				 cores=n.cores)
	return(fit.maaslin3$fit_data_abundance$results)
}

args <- commandArgs(trailingOnly=TRUE)
file_name = args[1]
base_file <- tools::file_path_sans_ext(basename(file_name))
sim_data <- readRDS(file_name)
out_dir <- args[2]
tp <- 1.3
tp_s <- 0.25
n.cores <- 5

Y <- sim_data$Y
X <- sim_data$X
D <- nrow(Y)
N <- ncol(Y)
other <- colSums(Y[rowSums(Y==0)>(0.7*N),])
seq.trim <- Y[rowSums(Y==0)<=(0.7*N),]
seq.trim <- rbind(seq.trim, other)

tryCatch({
	df <- data.frame(sim_data$X[,c("disease", "subject")])
	df[, "disease"] <- factor(df[,"disease"], levels=c(0,1))
	df[,"subject"] <- factor(df[,"subject"])
	row.names(df) <- colnames(seq.trim)
	maaslin3.res <- run_maaslin3(seq.trim, df, "disease+(1|subject)",
	                             normalization="TSS", n.cores=n.cores)
	saveRDS(maaslin3.res, cat.path(out_dir, paste0(base_file, "_maaslin3.RDS")))
}, error=function(x){print(x);})
