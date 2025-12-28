#sim_mat <- rbind(c("oral", 0.1, 0.4, -2, 3, 1, 2), c("oral", 0.1, 0.7, -2, 3, 1, 2),
#                 c("armpit", 0, 0.5, 0, 3, 1, 2), c("gut", 0.1, 0.7, -2, 3, 1, 2),
#                 c("oral", 0, 0.2, 0, 3, 1, 2), c("gut", 0, 0.2, 0, 3, 1, 2))
sim_mat <- rbind(c("oral", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5),
                 c("gut", 0.2, 0.2, -2.5, 2.5, -0.5, 0.5))
for(n_s in c(80, 120, 160, 200, 250, 300)) {
 for(i in 1:nrow(sim_mat)) {
    d <- sim_mat[i,]
##    for(j in 1:20) {
for(j in 1:20) {
     system(paste0("Rscript run_maaslin3.R ../raw_data/sim_inputs/", d[1], "_", n_s,"_20_1_", d[2], "_", d[3],
            "_", d[4], "_", d[5], "_", j, ".RDS ../raw_data/sim_outputs/"))
  }
 }
}
