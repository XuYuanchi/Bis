# ALRA for ti datasets
source("/home/suyanchi/project/ALRA/alra.R")

#load data
dat_sc = read.table(file = '/home/suyanchi/project/dab/data/ti/gse.csv', sep = ',', header = T)
rownames(dat_sc) = dat_sc[,1]
dat_sc = dat_sc[,-1]


dat_sc <- normalize_data(t(dat_sc))
# Choose k.
k_choice <- choose_k(dat_sc)

A_norm_completed <- alra(dat_sc,k=k_choice$k)[[3]]

saveRDS(t(A_norm_completed), file='/home/suyanchi/project/dab/results/ti/alra_gse.rds')
