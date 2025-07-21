sim_mat <- rbind(c("oral", 0.1, 0.4, -2, 3, 1, 2), c("oral", 0.1, 0.7, -2, 3, 1, 2),
                 c("armpit", 0, 0.5, 0, 3, 1, 2))
for(n_s in c(6, 8, 12, 20, 40, 80, 120, 160, 200, 250, 300)) {
 for(i in 1:3) {
    d <- sim_mat[i,]
##    for(j in 1:20) {
for(j in 1:20) {
     cat(paste0("\"sbatch run_one.sh ./output/", d[1], "_", n_s,"_20_1_", d[2], "_", d[3],
                  "_", d[4], "_", d[5], "_", j, ".RDS ./results/\" \n"))
  }
 }
}
