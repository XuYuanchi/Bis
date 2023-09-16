# ALRA for 1M datasets
source("/home/suyanchi/project/ALRA/alra.R")
library(rhdf5)
#load data
h5file=H5Fopen('/home/suyanchi/project/dab/data/1M/1M.h5')
h5dump(h5file, load = FALSE)
dat_sc = h5file$data
h5closeAll()

dat_sc <- normalize_data(t(dat_sc))
# Choose k.
k_choice <- choose_k(dat_sc)

A_norm_completed <- alra(dat_sc,k=k_choice$k)[[3]]

#saveRDS(t(A_norm_completed), file='/home/suyanchi/project/dab/results/1M/alra.rds')
h5createFile('/home/suyanchi/project/dab/results/1M/alra.h5')
h5write(A_norm_completed, '/home/suyanchi/project/dab/results/1M/alra.h5', 'data')
