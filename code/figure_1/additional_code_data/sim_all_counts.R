#sim_mat <- rbind(c("oral", 0.1, 0.4, -2, 3, 1, 2), c("oral", 0.1, 0.7, -2, 3, 1, 2),
#                 c("armpit", 0, 0.5, 0, 3, 1, 2), c("gut", 0.1, 0.7, -2, 3, 1, 2),
#		 c("oral", 0, 0.2, 0, 3, 1, 2), c("gut", 0, 0.2, 0, 3, 1, 2))
sim_mat <- rbind(c("oral", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5),
		 c("gut", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5))
for(n_s in c(6, 8, 12, 20, 40, 80, 120, 160, 200, 250, 300)) {
 for(i in 1:nrow(sim_mat)) {
    d <- sim_mat[i,]
##    for(j in 1:20) {
for(j in 1:20) {
      print(paste0("Rscript sim_counts.R ", n_s, " 20 1 ", j,
                    " ", d[1], " ", d[2], " ", d[3], " ", d[4], " ", d[5], " ", d[6], " ", d[7], " ./output/"))
      system(paste0("Rscript sim_counts.R ", n_s, " 20 1 ", j,
                    " ", d[1], " ", d[2], " ", d[3], " ", d[4], " ", d[5], " ", d[6], " ", d[7], " ./output/"))
    }
  }
}
