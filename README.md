# Optimal Transport based and External Prior Guided Autoencoder for Single-Cell Data Imputation
Single-cell RNA sequencing (scRNA-seq) is a powerful technique that allows researchers to study gene expression patterns in individual cells. However, scRNA-seq data is often noisy and incomplete due to the technical limitations of the sequencing process. To address this issue, researchers have developed various methods for the imputation of scRNA-seq data. In this study, we propose a model that does not assume a specific probability distribution for the scRNA-seq data and combines optimal transport and an autoencoder to infer missing expression values while preserving the global structure of the scRNA-seq data. In addition, our model leverages the geometry of the cellular embedding space, as well as the prior information of consistency between scRNA-seq and bulk RNA-seq data, to improve the imputation accuracy. Extensive experimental results demonstrate that the proposed model achieves competitive performance for simulated and real scRNA-seq data imputation while enhancing the performance of downstream clustering tasks and maintaining attractive computational efficiency.

## Overview
<div align=center>
<img src="https://github.com/XuYuanchi/Bis/blob/main/framework.png" height="600" width="400">
</div>
## Installation
The code is implemented with Python, pytorch, and pytorch-lightning...
Please refer to the configuration file in details.

## Run the demos
python train.py --dataset_name Spleen
