library(ktplots)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(RColorBrewer)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

margin_spacer <- function(x) {
  # where x is the column in your dataset
  left_length <- nchar(levels(factor(x)))[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}

cellPhoneDBDir <- '20250509_output/Cell-Cell-Interaction/1.CellPhoneDB/AllCells'

meta <- read.table('20250509_output/Cell-Cell-Interaction/1.CellPhoneDB/AllCells/AllCells_Mild_Advanced.celltype.CellphoneDB.metadata.txt', header = T,check.names = F,sep = '\t',row.names = 1)
unique(meta$cell_type)
# specify the directory to use for in/out
expt.list=list('Mild_out'='Mild', 
               'Advanced_out'='Advanced',
               'Mild_Advanced_out'='Mild_Advanced')

res.list=vector('list', length = 3)
names(res.list) <- names(expt.list)

experimentDir <- names(expt.list)[3]
for(experimentDir in names(expt.list)){
  cellGroupOfInterest <- c("Hepatocytes|Cholangiocytes", "Cholangiocytes|Hepatocytes")
  currDisStates <- expt.list[[experimentDir]][1]
  
  # input directory and output directory
  inDir <- paste0(cellPhoneDBDir, experimentDir, '/')
  outDir <- paste0(cellPhoneDBDir, experimentDir, '/plot/')
  
  # read the data ----
  # the identity columns in common across the input
  id.cols <- c('id_cp_interaction', 'interacting_pair', 
               'partner_a', 'partner_b', 'gene_a', 'gene_b',
               'secreted', 'receptor_a', 'receptor_b', 
               'annotation_strategy', 'is_integrin')
  
  # significant means
  smeans <- fread(paste0(inDir, 'significant_means.txt'))

  smeans.m <- unique(melt(smeans, id.vars=c(id.cols, 'rank'),
                          variable.name='interacting_cells',
                          value.name='mean_expression'))

  # means <- fread(paste0(inDir, 'means.txt'))
  pvals <- fread(paste0(inDir, 'pvalues.txt'))
  
  pvals.m <- unique(melt(pvals, id.vars=id.cols,
                         variable.name='interacting_cells', 
                         value.name='p_value'))
  
  # combine the expression and p-vals
  results <- merge(smeans.m, pvals.m, 
                   by=c(id.cols, 'interacting_cells'), all.y=T)
  
  rm(pvals, pvals.m, smeans, smeans.m)
  
  
  # filter and process the data ----
  
  # filter for comparisons between cell group of interest and others
  # results.filt <- results[grepl(cellGroupOfInterest, results$interacting_cells),]
  results.filt <- results[results$interacting_cells %in% cellGroupOfInterest, ]
  
  if (nrow(results.filt)==0) {
    print(paste0('no  interactions for cluster : ', cellGroupOfInterest, ' in ', currDisStates))
    print('check this is a real cluster in the data!')
    exit()
  }
  
  # apply FDR multiple correction for p-val
  results.filt$p_value_corrected <- p.adjust(results.filt$p_value, 
                                             method='BH')
  
  # filter for significant
  # results.filt <- results.filt[results.filt$p_value_corrected<=0.05,]
  results.filt <- results.filt[results.filt$p_value<=0.05, ]
  # ranking according to the mean expression
  results.filt <- results.filt[order(results.filt$mean_expression, decreasing = T), ]
  results.filt$currDisStates <- currDisStates
  results.filt$interacting_cells_currDisStates <- paste0(results.filt$interacting_cells, '_', results.filt$currDisStates)
  res.list[[experimentDir]] <- results.filt
}

res <- do.call('rbind', res.list)
res$interacting_cells_currDisStates <- factor(res$interacting_cells_currDisStates, unique(res$interacting_cells_currDisStates))

p1 <- ggplot(res, aes(y=interacting_pair, x=interacting_cells_currDisStates))+
  geom_point(aes(color=mean_expression), size=3)+
  # geom_point(aes(color=mean_expression, size=mean_expression))+
  # geom_point(aes(color=mean_expression))+
  # scale_color_gradientn(colours = viridis::viridis(20)) +
  scale_colour_gradientn(colours = colorRampPalette(brewer.pal(9, "OrRd"))(100)) +
  xlab('')+ ylab('Gene pair (ligand–receptor)')+
  coord_flip() +
  # ggtitle(paste0('first element is expressed in ', cellGroupOfInterest))+
  theme_bw()+  
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.x=element_text(size=8, angle=30, hjust=1, vjust=1),
        axis.text.y=element_text(size=14),
        axis.title.x = element_text(size = 20),   
        axis.title.y = element_text(size = 10))

p1 

pdf(paste0(cellPhoneDBDir, strsplit(cellGroupOfInterest, "\\|")[[1]][1], '_', strsplit(cellGroupOfInterest, "\\|")[[1]][2], '_filtered_significance_plots_coord_flip.pdf'),
    width=20, height=6)
print(p1)
dev.off()

## fill with gradient color
p1 <- ggplot(res, aes(y=interacting_pair, x=interacting_cells_currDisStates, fill=mean_expression))+ 
  geom_tile(colour='white', size=0.5)+
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, "OrRd"))(100))+
  xlab('')+
  theme(#plot.background=element_blank(),
    panel.background=element_blank(),
    #panel.grid.major=element_line(colour=NA),
    #panel.grid.minor=element_line(colour=NA),
    #axis.ticks=element_blank(),
    axis.text.x=element_text(size = 14, angle=30, hjust=1, vjust=1),
    axis.text.y=element_text(size = 10)
  )

p1 
pdf(paste0(cellPhoneDBDir, strsplit(cellGroupOfInterest, "\\|")[[1]][1], '_', strsplit(cellGroupOfInterest, "\\|")[[1]][2], '_filtered_significance_plots.pdf'),
    width=10, height=15)
print(p1)
dev.off()
