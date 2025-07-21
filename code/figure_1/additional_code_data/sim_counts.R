library(SparseDOSSA2)
library(lme4)
source("sparse_trim.R")

args <- commandArgs(trailingOnly=TRUE)

n_s <- as.numeric(args[1])
n_r <- as.numeric(args[2])
vs <- as.numeric(args[3])
ind <- as.numeric(args[4])
model <- args[5]
pc_l <- as.numeric(args[6])
pc_u <- as.numeric(args[7])
lb <- as.numeric(args[8])
ub <- as.numeric(args[9])
lbs <- as.numeric(args[10])
ubs <- as.numeric(args[11])
seq_depth <- 80000
out_dir <- args[12]

set.seed(round(110*pc_l)+round(225*pc_u)+n_s+n_r+ind*4517)

if(model=="oral") {
  oral_model <- readRDS("oral.RDS")
  d <- SparseDOSSA2_m(template=oral_model,
                      new_features=F,
                      n_sample=10000)
} else if(model=="armpit"){
  vagina_model <- readRDS("armpit.RDS")
  d <- SparseDOSSA2_m(template=vagina_model,
                    new_features=F,
                    n_sample=10000)
}

Y <- apply(d$simulated_matrices$rel, 2, function(col) rmultinom(1, seq_depth, col)) ## 100000
row.names(Y) <- row.names(d$simulated_matrices$rel)
ntaxa <- nrow(Y)
use_names <- row.names(Y[rowSums(Y==0)<=7000,])
print(length(use_names))

counter <- 1
while(TRUE) {

	taxa_names <- row.names(d$simulated_matrices$rel)
	ninc <- round(length(use_names)*pc_u)
	ndec <- round(length(use_names)*pc_l)

	all_lfcs <- rep(0, ntaxa)
	names(all_lfcs) <- taxa_names
	inc_names <- sample(use_names, ninc, replace=F)
	print(length(inc_names))
	dec_names <- sample(use_names[!(use_names%in%inc_names)], ndec, replace=F)
	print(length(dec_names))
	all_lfcs[inc_names] <- runif(length(inc_names), 0, ub)
	all_lfcs[dec_names] <- runif(length(dec_names), lb, 0)

	n <- n_s * n_r
	mat_metadata <- data.frame(disease=c(rep(0, n/2), rep(1, n/2)))
	##mat_metadata <- data.frame(disease=rep(c(rep(0, n_r/2), rep(1, n_r/2)), n_s))
	mat_metadata <- cbind(mat_metadata, matrix(0, nrow=n, ncol=length(all_lfcs)))
	colnames(mat_metadata) <- c("disease",
				    paste("microbe_eff", 1:length(all_lfcs), sep="_"))
	mat_metadata[,2:ncol(mat_metadata)] <- apply(mat_metadata[,2:ncol(mat_metadata)],
						     2, function(col) {
							 rep(rnorm(n=n_s, 0, vs), each=n_r)
						     })
	spike_in_df <- data.frame()
	for(i in 1:(length(all_lfcs)-1)) {
	    spike_in_df <- rbind(spike_in_df, data.frame(metadata_datum=1,
							 feature_spiked=taxa_names[i],
							 effect_size=all_lfcs[i],
							 associated_property="abundance"))
	    spike_in_df <- rbind(spike_in_df, data.frame(metadata_datum=1,
							 feature_spiked=taxa_names[i],
							 effect_size=all_lfcs[i]/4,
							 associated_property="prevalence"))
	    spike_in_df <- rbind(spike_in_df, data.frame(metadata_datum=i+1,
							 feature_spiked=taxa_names[i],
							 effect_size=1,
							 associated_property="abundance"))
	}

	if(model=="oral") {
	  oral_model <- readRDS("oral.RDS")
	  d <- SparseDOSSA2_m(template=oral_model,
			      new_features=F,
			      metadata_matrix=as.matrix(mat_metadata),
			      spike_metadata=spike_in_df,
			      n_sample=nrow(mat_metadata))
	} else if(model=="armpit") {
	  vagina_model <- readRDS("armpit.RDS")
	  d <- SparseDOSSA2_m(template=vagina_model,
			    new_features=F,
			    metadata_matrix=as.matrix(mat_metadata),
			    spike_metadata=spike_in_df,
			    n_sample=nrow(mat_metadata))
	}

	#X <- data.frame(disease=rep(c(rep(0, n_r/2), rep(1, n_r/2)), n_s),
	#		subject=paste0("subject_", rep(1:n_s, each=n_r)))
	X <- data.frame(disease=c(rep(0, n/2), rep(1, n/2)),
			subject=paste0("subject_", rep(1:n_s, each=n_r)))
	A <- d$simulated_matrices$a_spiked
	Al <- log2(colSums(A))
	sdata <- data.frame(Al=Al, disease=X$disease, subject=X$subject)
	theta_perp <- coef(summary(lmer(Al~disease+(1|subject), data=sdata)))["disease",1]
	print(c(theta_perp, counter))
	if((theta_perp >= lbs) & (theta_perp <= ubs)) {
		break
	} else {counter <- counter + 1}
}

Y <- apply(d$simulated_matrices$rel, 2, function(col) rmultinom(1, seq_depth, col))
print(names(d$simulated_matrices))
row.names(Y) <- taxa_names
colnames(Y) <- paste0("sample_", 1:ncol(Y))
row.names(X) <- paste0("sample_", 1:ncol(Y))
saveRDS(
  list(Y=Y, A=A, X=X, all_lfcs=all_lfcs),
  paste(paste0(out_dir, model), n_s, n_r, vs, pc_l, pc_u, lb, ub, paste0(ind, ".RDS"), sep="_")
)
