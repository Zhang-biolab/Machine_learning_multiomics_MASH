#!/bin/bash
# Args

Rscript=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_scripts/SCENIC/Hep/01.Hep-Runnig-SCENIC-on-SCT-Normalized-Assay.R
dir=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_scripts/SCENIC/Hep/

cd $dir
/home/dell/github/R-4.4.1/bin/R CMD BATCH $Rscript 01.Hep-Runnig-SCENIC-on-SCT-Normalized-Assay.R.sh.log

