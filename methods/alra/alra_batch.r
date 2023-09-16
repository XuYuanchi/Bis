source("/home/suyanchi/project/ALRA/alra.R")
library(rhdf5)
library(R.matlab)
library(anndata)

# load data
dataset_name = c("cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas")
for (i in c(5)) {
    # h5file = H5Fopen(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".h5"))
    # h5dump(h5file, load = FALSE)
    # dat_sc = h5file$data
    # h5closeAll()
    # print(dataset_name[i])
    # print(dim(dat_sc))
    # dat_sc = read.table(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".csv"), sep = ",", header = T)
    dat_sc = read_h5ad(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".h5ad"))
    # rownames(dat_sc) = dat_sc[, 1]
    # dat_sc = dat_sc[, -1]
    # dat_sc = readMat('/home/suyanchi/project/dab/data/batch/pancreas.mat')[[1]]

    dat_sc <- normalize_data(as.matrix(dat_sc$X))
    # Choose k.
    k_choice <- choose_k(dat_sc)

    A_norm_completed <- alra(dat_sc, k = k_choice$k)[[3]]

    # saveRDS(t(A_norm_completed), file='/home/suyanchi/project/dab/results/1M/alra.rds')
    # writeMat(paste0("/home/suyanchi/project/dab/results/batch/", dataset_name[i], "_alra.mat"), data = A_norm_completed)
    path = paste0("/home/suyanchi/project/dab/results/batch/", dataset_name[i], "_alra.h5")
    h5createFile(path)
    h5write(A_norm_completed, path, "data")
    h5closeAll()
}
