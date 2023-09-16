library(SAVER)
library(R.matlab)
library(rhdf5)
library(anndata)

dataset_name = c("cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas")
for (i in c(5)) {
    # dat_sc = read.table(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".csv"), sep = ",", header = T)
    # dat_sc = readMat(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".mat"))[[1]]
    dat_sc = read_h5ad(paste0("/home/suyanchi/project/dab/data/batch/", dataset_name[i], ".h5ad"))

    saver_result <- saver(t(as.matrix(dat_sc$X)), estimates.only = TRUE, ncores = 20)

    # writeMat(paste0("/home/suyanchi/project/dab/results/batch/", dataset_name[i], "_saver.mat"), data = saver_result)
    
    path = paste0("/home/suyanchi/project/dab/results/batch/", dataset_name[i], '_saver.h5')
    h5createFile(path)
    h5write(saver_result, path, "data")
    h5closeAll()
}