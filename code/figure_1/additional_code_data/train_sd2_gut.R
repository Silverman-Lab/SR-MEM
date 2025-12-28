library(HMP2Data)
library(biomformat)
library(microbiomeDataSets)
library(phyloseq)
library(mia)
library(SparseDOSSA2)
library(microbiome)

set.seed(4326)

gut_data <- T2D16S()
gut_data <- prune_taxa(taxa_sums(gut_data) > 0, gut_data)
ps_genus <- gut_data #tax_glom(gut_data, "Species")
depths <- sample_sums(ps_genus)
top100_names <- names(sort(depths, decreasing = TRUE))[1:100]
ps_top100 <- prune_samples(top100_names, ps_genus)
ps_top100 <- prune_taxa(taxa_sums(ps_top100) > 0, ps_top100)
Y_gut <- otu_table(ps_top100, taxa_are_rows=T)
row.names(Y_gut) <- paste0("taxon_", 1:nrow(Y_gut))
Y_gut <- Y_gut[rowSums(Y_gut==0)<=70,]
dim(Y_gut)
total_rows <- nrow(Y_gut)
sampled_indices <- sample(1:total_rows, size = 200, replace = FALSE)
Y_gut <- Y_gut[sampled_indices, ]
dim(Y_gut)

## Create Counts
plan(future::multisession(workers=4))
fitted_parallel <- fit_SparseDOSSA2(data=Y_gut, lambda=0.3)
saveRDS(fitted_parallel, "gut.RDS")
