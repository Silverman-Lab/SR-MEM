set.seed(25334)

df <- read.csv("hmp2_metadata_2018-08-20.csv")
df <- df[df$data_type=="metagenomics",]
names_already <- list.files("~/work/figure_3/", "_metaphlan_profile.txt")
already <- sapply(strsplit(names_already, "_"), "[", 1)

for(eid in df$External.ID) {
	if(!(eid %in% already)){
		system(paste0("wget -P /storage/home/kvm6065/scratch https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/raw/HMP2/MGX/2018-05-04/", eid, ".tar"))
		cat(paste0('\"sbatch run_metaphlan.sh /storage/home/kvm6065/scratch/', eid, '_R1.fastq.gz /storage/home/kvm6065/scratch/', eid, '_R2.fastq.gz\"\n'))
	}
}

