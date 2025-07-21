in_files <- list.files(pattern="*_profile.txt")

l <- list()
count <- 1
for(in_file in in_files) {
	df <- read.table(in_file, comment.char="#")
	row.names(df) <- df[,1]
	df <- df[,2:ncol(df)]
	taxa_len <- unname(sapply(row.names(df), function(i) length(strsplit(i, "\\|")[[1]])))
	df <- df[taxa_len==6,]
	row.names(df) <- unname(sapply(row.names(df), function(i) strsplit(i, "\\|")[[1]][6]))
	l[[count]] <- df
	count <- count + 1
}


col_names <- sub("_metaphlan_profile\\.txt$", "", in_files)
genus_names <- unique(unlist(lapply(l, function(df) row.names(df))))

df <- data.frame(matrix(0, dimnames=list(genus_names, col_names),
			nrow=length(genus_names), ncol=length(col_names)))
for(i in 1:length(l)) {
	sdf <- l[[i]]
	df[row.names(sdf),i] <- sdf$V5
}

##df <- apply(df, 2, function(col) col/sum(col))

write.csv(df, "ibddmb.counts.csv")
