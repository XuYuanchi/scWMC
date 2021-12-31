# scWMC
Source code of "scWMC: weighted matrix completion-based imputation of scRNA-seq data via imperfect prior subspace information".

Single-cell RNA sequencing (scRNA-seq) can provide insight into the gene expression patterns at the resolution of individual cells, which offers new opportunities to study the behavior of different cell types. However, it is often plagued by dropout events, a phenomenon where the expression value of a gene tends to be measured as zero in the expression matrix due to various technical defects. In this paper, we argue that borrowing genes and cells information across column and row subspaces directly results in suboptimal solution due to noise contamination in imputing dropout values. Thus, to more precisely imputing the dropout events in scRNA-seq data, we develop a regularization for
leveraging those imperfect prior information to estimate underlying true prior subspace and embed it in a typical low-rank matrix completion-based framework, named scWMC. To evaluate the performance of the proposed method, we conduct comprehensive experiments on simulated and real scRNA-seq data. Extensive data analysis, including simulated analysis, cell clustering, differential expression analysis, and cell trajectory inference, demonstrate that our method can produce more accurate imputation results than competing methods and benefit subsequent downstream analysis.

## Overview
<img src="https://github.com/XuYuanchi/scWMC/blob/main/model.png" height="600" width="1000">

run Demo.m

The result is as follows

<img src="https://github.com/XuYuanchi/scWMC/blob/main/result_sWMC.png" height="300" width="1200">
