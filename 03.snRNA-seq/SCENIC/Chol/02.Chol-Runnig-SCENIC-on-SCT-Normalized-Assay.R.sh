#!/bin/bash
# Args

Rscript=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_scripts/SCENIC/Chol/01.Chol-Runnig-SCENIC-on-SCT-Normalized-Assay.R
dir=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_scripts/SCENIC/Chol/

cd $dir
/share/sequence/github/R-4.4.1/lib/R/bin/R CMD BATCH $Rscript 01.Chol-Runnig-SCENIC-on-SCT-Normalized-Assay.R.sh.log

