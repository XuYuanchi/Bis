library(scImpute)

filename = paste0("/home/suyanchi/project/dab/data/ti/gse.csv")
scimpute(count_path = filename, infile = "csv", outfile = "rds", out_dir = '/home/suyanchi/project/dab/results/ti/scimpute_gse.rds', 
labeled = FALSE, drop_thre = 0.5, Kcluster = 5, ncores = 5)