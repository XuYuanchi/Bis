library(SAVER)

dat_sc = read.table(file='/home/suyanchi/project/dab/data/ti/gse.csv', sep=',', header=T)

rownames(dat_sc) = dat_sc[,1]
dat_sc = dat_sc[,-1]

saver_result <- saver(dat_sc, estimates.only = TRUE, ncores=30)

saveRDS(saver_result, file='/home/suyanchi/project/dab/results/ti/saver_gse.rds')