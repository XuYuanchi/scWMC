# scWMC
Source code of "scWMC: weighted matrix completion-based imputation of scRNA-seq data via imperfect prior subspace information".

Single-cell RNA sequencing (scRNA-seq) can provide insight into the gene expression patterns at the resolution of individual cells, which offers new opportunities to study the behavior of different cell types. However, it is often plagued by dropout events, a phenomenon where the expression value of a gene tends to be measured as zero in the expression matrix due to various technical defects. In this paper, we argue that borrowing genes and cells information across column and row subspaces directly results in suboptimal solution due to noise contamination in imputing dropout values. Thus, to more precisely imputing the dropout events in scRNA-seq data, we develop a regularization for
leveraging those imperfect prior information to estimate underlying true prior subspace and embed it in a typical low-rank matrix completion-based framework, named scWMC. To evaluate the performance of the proposed method, we conduct comprehensive experiments on simulated and real scRNA-seq data. Extensive data analysis, including simulated analysis, cell clustering, differential expression analysis, and cell trajectory inference, demonstrate that our method can produce more accurate imputation results than competing methods and benefit subsequent downstream analysis.

## Overview
<img src="https://github.com/XuYuanchi/scWMC/blob/main/model.png" height="600" width="1000">

## Installation
### Requirements
Our algorithm is implemented using Matlab R2020b.
### Installation from GitHub
To clone the repository, run the following from a terminal:
```
git clone git://github.com/XuYuanchi/scWMC.git
```

## Quick start
There is a simulated datasets in data/demo_data.mat file. There are the true data set, Dropout data set, and the imputed dataset by scWMC.

1. Demo.m -- Shows how to run scWMC.

The result is as follows

<img src="https://github.com/XuYuanchi/scWMC/blob/main/result_sWMC.png" height="300" width="1200">

## Usage
### load the dropout data, such as:
```
data_dropout = load('demo_data.mat');
```
Or you can use function <b>readtable</b> to load CSV/TSV/TXT-formatted raw count matrix with genes in rows and cells in columns. Cell and gene labels aren't necessary.
### Prepare the parameters for scWMC
Set up the parameters used in example
```
Par.lam  = 0.8;
Par.rho  = 0.8;
Par.mu1  = 0.00001;
Par.mu2  = 0.00001;
Par.iter = 100;
```
lam, rho     -- denote the weights reflecting the uncertainty in the prior subspace information, set as <b>dropout rate</b>

mu1, mu2 -- denote the regularization parameters, generally do not need to change

iter         -- denote the number of iteration, generally do not need to change
### Run scWMC
```
## imputation
dataRecovered        = impute(data_dropout, Par);
## 
dataRecovered        = max(dataRecovered, 0);
index                = find(data_dropout);
dataRecovered(index) = data_dropout(index);
```
### Visualize the results
```
gcf = figure(1);
set(gcf, 'Position', [100, 500, 1200, 300])
subplot(1,3,1)
imagesc(log10(data_true+1))
title('True Data')
axis off
subplot(1,3,2)
imagesc(log10(data_dropout+1))
title('Drop-out Data')
axis off
subplot(1,3,3)
imagesc(log10(dataRecovered+1))
title('Imputed Data by scWMC')
axis off
```
## Main functions
[impute.m](https://github.com/XuYuanchi/scWMC/blob/main/utils/impute.m)

main sWMC algorithm, there are two parameters in this function

```
Ori_P -- gene expression matrix which genes in rows and cells in columns needs to imputation
Par   -- parameters for scWMC, including lam, rho, mu1, mu2 and iter. lam and rho need to adjust according to dropout rate and others don't.
```
## scWMC_Reproducibility
The data and codes for reproducing all Figures and Tables in the manuscript of scWMC can be download from:
```
https://pan.baidu.com/s/1__KdQuVDOgrGx5fDEe1CLQ
```
code: abcd
## Calling in R
### install R package matlabr
You can install the stable version on CRAN:
```
install.packages('matlabr', dependencies = TRUE)
```
or you can install matlabr from GitHub with:
```
# install.packages("remotes")
remotes::install_github("muschellij2/matlabr")
```
### run scWMC in R such as
```
library(matlabr)
run_matlab_script('Demo.m')
```
## Contact

* Please feel free to contact Yanchi Su (suyanchi@gmail.com) if you have any questions about the software.

## License

This project is licensed under the MIT License.
