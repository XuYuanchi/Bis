# Distribution-agnostic Deep Learning Enables Accurate Single‐Cell Data Imputation and Transcriptional Regulation Interpretation

## Abstract
Single-cell RNA sequencing (scRNA-seq) offers a robust methodology for investigating gene expression at the single-cell level. However, accurate quantification of genetic material is often hindered by the limited capture of intracellular mRNA, resulting in a large number of missing expression values, which impedes downstream analysis. Existing imputation methods rely heavily on stringent data assumptions, such as specific probability distributions, that restrict their broader application. Moreover, the recovery process lacks reliable supervision, leading to bias in gene expression signal imputation. To address these challenges, we developed a distribution-agnostic deep learning model, called Bis, for the accurate imputation of scRNA-seq data from multiple platforms. Bis is an optimal transport-based autoencoder model that can capture the intricate distribution of scRNA-seq data while addressing the characteristic sparsity by regularizing the cellular embedding space. After that, we propose a transcriptional expression consistency module that leverages bulk RNA-seq data as external priors to guide the imputation process and constrain the model to ensure consistency of the average gene expression between the aggregated imputed and bulk RNA-seq data. Experimental results validated that Bis outperforms other state-of-the-art models across numerous simulated datasets and diverse real scRNA-seq data generated by different representative single-cell sequencing platforms. Moreover, we showcase that Bis consistently achieved accuracy and efficacy in varied types of downstream analyses encompassing batch effect removal, clustering analysis, differential expression analysis, and trajectory inference. In addition, we demonstrated that in a tumor-matched peripheral blood dataset, Bis successfully restored the gene expression levels of rare cell subsets to unveil the developmental characteristics of cytokine-induced NK cells within a head and neck squamous cell carcinoma microenvironment.

## Overview
<div align=center>
<img src="https://github.com/XuYuanchi/Bis/blob/main/framework.png" height="400" width="800">
</div>
## Installation
The code is implemented with Python, pytorch, and pytorch-lightning...
Please refer to the configuration file in details.

## Run the demos
python train.py --dataset_name Spleen
