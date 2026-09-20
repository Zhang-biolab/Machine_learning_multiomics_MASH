library(tidyverse)
library(dplyr)
library(sva)

setwd("~/projects/NASH_Fibrosis_datamining/")

input_dir <- "data/HSA/"
output_dir <- "MergeCohort_construction/"

# merge sample information for four datasets
tmp_sampleinfo_list <- list()
tmp_file_path <- list.files(path = input_dir, "^phenotype_clean.txt",recursive = T,full.names = T)[-5]

for(i in 1:length(tmp_file_path)){
  dataset_name <- str_match(tmp_file_path[[i]],"GSE[0-9]+|PRJNA_merge") %>% 
    as.character()
  rt <- read.delim(tmp_file_path[[i]],check.names = F)
  tmp_sampleinfo_list[[i]] <- data.frame(rt,GSE_name=dataset_name,Batch_GSE=i) %>% 
    mutate(sample_id=paste0(GSE_name,'_',sample_name2))
}

sample_info <- do.call("rbind", tmp_sampleinfo_list) %>% 
  mutate(group=factor(group,group_levles)) %>%
  arrange(group,Batch_GSE)

write.table(
  sample_info,
  paste0(output_dir,'sample_info_with_fibrosisScore.txt'),
  sep = '\t',
  quote = F,
  row.names = F
)

# merge gene expression matrix for four datasets 
tmp_list <- list()
tmp_sampleinfo_list <- list()
tmp_file_path <- list.files(path = input_dir, '^vst_normalized_no_duplicate_remove_outliers.xls', recursive = T, full.names = T)

for(i in 1:length(tmp_file_path)){
  dataset_name <- str_match(tmp_file_path[[i]],'GSE[0-9]+|PRJNA_merge') %>% 
    as.character()
  rt <- read.delim(tmp_file_path[[i]],check.names = F)
  tmp_sampleifo_list[[i]] <- data.frame(sample_id=colnames(rt[,-1]),GSE_name=dataset_name,Batch_GSE=i) %>% 
    mutate(group=gsub('-[0-9]+$','',sample_id)) %>% 
    rowwise() %>% 
    mutate(sample_id=paste0(GSE_name,'_',sample_id)) %>% 
    dplyr::select(sample_id,GSE_name,Batch_GSE,group)
    colnames(rt) <- c('gene',paste0(dataset_name,'_',colnames(rt[,-1])))
    tmp_list[[i]] <- rt
}

group_levles <- c('Mild','Adv')
tmp_mtx <- Reduce(function(df1,df2) merge(df1,df2,by='gene'),tmp_list)
sample_info <- bind_rows(tmp_sampleinfo_list) %>% 
  mutate(group=factor(group,group_levles)) %>%
  arrange(group, Batch_GSE)
write.table(
    sample_info,
    paste0(output_dir,'sample_info_sva.txt'),
    sep = '\t',
    quote = F,
    row.names = F
)
exp_mtx <- tmp_mtx[,c('gene',sample_info$sample_id)] # the number of intersect gene is 13505
write.table(
  exp_mtx,
  paste0(output_dir,'vst_mtx_sva.txt'),
  sep = '\t',
  quote = F,
  row.names = F
)

exp <- exp_mtx %>% column_to_rownames('gene') %>% as.matrix()
sample_info$hsaAdv <- as.numeric(sample_info$group == "Adv")
mod <- model.matrix(~hsaAdv,data=sample_info)
batchtype <- as.numeric(sample_info$Batch_GSE)
outTab <- as.data.frame(ComBat(exp, batchtype, mod, par.prior=T))
batch_correct_mtx <- outTab %>% rownames_to_column('gene')
write.table(
  batch_correct_mtx,
  paste0(output_dir,'batch_correct_vst_mtx_sva.txt'),
  sep = '\t',
  quote = F,
  row.names = F
)

