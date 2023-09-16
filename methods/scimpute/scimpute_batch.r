library(scImpute)
library(R.matlab)


dataset_name = c("cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas")
for(i in c(5)){
filename = paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".rds")
scimpute(count_path = filename, infile = "rds", outfile = "rds", out_dir = paste0('/home/suyanchi/project/dab/results/batch/', dataset_name[i],'_scimpute.rds'), 
labeled = FALSE, drop_thre = 0.5, Kcluster = 5, ncores = 20)
}