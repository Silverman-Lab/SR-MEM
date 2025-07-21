library(HMP2Data)
library(biomformat)
library(microbiomeDataSets)
library(phyloseq)
library(mia)
library(SparseDOSSA2)
library(microbiome)

set.seed(1235)

## Oral
ps_mom <- momspi16S()
ps_oral <- subset_samples(ps_mom, sample_body_site=="buccal mucosa")
ps_oral <- prune_taxa(taxa_sums(ps_oral) > 0, ps_oral)
ps_genus <- ps_oral ## tax_glom(ps_oral, "Genus")
depths <- sample_sums(ps_genus)
top100_names <- names(sort(depths, decreasing = TRUE))[1:100]
ps_top100 <- prune_samples(top100_names, ps_genus)
ps_top100 <- prune_taxa(taxa_sums(ps_top100) > 0, ps_top100)
Y_oral <- otu_table(ps_top100, taxa_are_rows=T)
row.names(Y_oral) <- paste0("taxon_", 1:nrow(Y_oral))
Y_oral <- Y_oral[rowSums(Y_oral==0)<=70,]

## Create Counts
plan(future::multisession(workers=10))
fitted_parallel <- fit_SparseDOSSA2(data=Y_oral,
                                    lambda=0.3)
saveRDS(fitted_parallel, "oral.RDS")

## Read Data
biom_file <- "./46352_otu_table.biom"
biom_data <- read_hdf5_biom(biom_file)
counts <- do.call(rbind, biom_data[[12]])
genera <- sapply(biom_data[[9]], function(item) item$metadata$taxonomy[4])
row.names(counts) <- genera

## sample type feces, treatment control, diabetes 0
metadata <- read.table("10370_20230206-075001.txt", sep="\t", header=T)
meta <- metadata[metadata[,"sample_name"]%in%colnames(counts),]
meta <- meta[(meta[,"timepoint"]=="Time 0")&(meta[,"host_body_site"]=="uberon:arm pit"),]
Y <- counts[,meta[,"sample_name"]]
Y <- rowsum(Y, group = rownames(Y))
Y <- Y[row.names(Y)!="f__",]
Y <- Y[,colSums(Y)>=1000]
Y <- Y[rowSums(Y==0)<=(ncol(Y)*0.7),]
dim(Y)

## Create Counts
plan(future::multisession(workers=10))
fitted_parallel <- fit_SparseDOSSA2(data=Y,
                                    lambda=0.3)
saveRDS(fitted_parallel, "armpit.RDS")
