# Machine_learning_multiomics_MASH

This repository contains scripts used for the integrative analysis of bulk and single-nucleus RNA-seq data and other downstream analyses currently included in the manuscript. 

All analyses were carried using R (v.4.4.1) and python (v.3.8)

### Packages and libraries used

The following packages are required to run the scripts:

*Bulk RNAseq data processing and integration:*
- DESeq2 (v.1.46.0)
- RobustRankAggreg (v.1.2.1)
- sva (v.3.54.0)

*machine learning:*
- glmnet (v.4.1-10)
- randomForest (v.4.7-1.2)
- e1071 (v.1.7-16)

*Single cell data processing and integration:*
- Seurat (v.4.4.0) 
- scDblFinder (v. 1.20.2)
- DecontX (v.1.4.1)

*Downstream analyses:*
- FactoMineR (v.2.12)
- factoextra (v.1.0.7)
- WGCNA (v.1.73)
- fgsea (v.1.32.4)
- GSVA (v.2.0.7)
- ggplot2 (v.3.5.2)
- pheatmap (v.1.0.13)
- tidyverse (v.2.0.0)
- data.table (v.1.18.6.1)
- reticulate (v.1.18)
- RColorBrewer (v.1.1-3)
- ggrepel(v.0.9.8)
- dplyr (v.1.2.1)
- stringr (v.1.6.0)
- pROC (v.1.19.0.1)
- scCODA (v.0.1.9)
- CellPhoneDB (v.2.1.7)
- SCENIC (v.1.3.1)

### GEO data
The data generated for this study are deposited in GEO database with the following accession numbers: GSE130970, GSE135251, GSE162694, GSE225740, GSE174478, GSE202379, GSE205846, GSE137449, GSE114261, GSE162863, GSE233767, GSE162869, and GSE119340.