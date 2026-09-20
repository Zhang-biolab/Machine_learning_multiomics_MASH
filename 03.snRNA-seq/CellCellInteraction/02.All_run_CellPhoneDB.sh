#!/bin/bash

## Usage
# 1. Make sure there is the METADATA folder in your h5ad dir
# 2. Run: cpdb cpdb.sh META.txt OUTDIR
# 3. For multiple files: use make_cpdb_script.py
# 4. Then use custom scripts for downstream analysis: python3.6 cpdb_heatmap.py

# Args
set -e

source activate /home/dell/share/data/junjie/anaconda3/envs/cpdb

cellphonedb=~/anaconda3/envs/cpdb/bin/cellphonedb
inDir=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_output/Cell-Cell-Interaction/1.CellPhoneDB/AllCells
outDir=~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_output/Cell-Cell-Interaction/1.CellPhoneDB/AllCells/Mild_Advanced_out
cd $inDir


# Notes:
## method can be one of: 'analysis' , 'statistical_analysis'
## counts-data  #[ensembl | gene_name | hgnc_symbol] Type of gene identifiers in the counts data

## 1. Run method with statistical analysis
${cellphonedb} method statistical_analysis AllCells_Mild_Advanced.celltype.CellphoneDB.metadata.txt AllCells_Mild_Advanced.celltype.CellphoneDB.count.txt \
--counts-data hgnc_symbol \
--threads=8 \
--threshold=0.1 \
--output-path ${outDir}
