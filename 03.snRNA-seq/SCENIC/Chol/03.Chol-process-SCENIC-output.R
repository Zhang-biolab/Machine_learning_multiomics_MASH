library(Seurat) # packageVersion("Seurat") V4.4.0
library(SingleCellExperiment)
library(tidyverse)
library(foreach)
library(RcisTarget)
library(doParallel)
library(SCopeLoomR)
library(AUCell)
library(doRNG)
library(GENIE3)
library(SCENIC)
library(dplyr) 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")
set.seed(123)
seuDir <- '~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/20250509_output/Subset_Chol/'
outDir <- 'plot/'
dir.create(outDir, recursive = T)

Chol.seu <- readRDS(paste0(seuDir, 'Aggr_Sep2024_filt_Chol_harmony_th=0.rds'))
m <- seu@meta.data
Idents(Chol.seu) <- 'SCT_snn_harmony_t.0.0.1'
cellInfo <- data.frame(CellType=Idents(Chol.seu))

scenicOptions <- readRDS('int/scenicOptions.Rds')
skipBinaryThresholds=FALSE
skipHeatmap=FALSE
skipTsne=FALSE

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC") # 408 regulon totally
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ] # 269 unique regulon
identical(rownames(cellInfo), colnames(regulonAUC))

selectedResolution <- "CellType"

# 1.RSS dotplot
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), selectedResolution], )
rss <- na.omit(rss) 
pdf(paste0(outDir, '01.RSS_plot.pdf'), width = 10, height = 20)
rssPlot <- plotRSS(rss,
                   labelsToDiscard = NULL, 
                   zThreshold = 1,
                   cluster_columns = FALSE, 
                   order_rows = T, 
                   thr = 0.01, 
                   varName = "cellType",
                   col.low = '#330066',  
                   col.mid = '#66CC66',  
                   col.high= '#FFCC33',
                   revCol = F,
                   verbose = TRUE)

rssPlot$plot
rssPlot$df
dev.off()

# 2.RSS rank set 
for (i in unique(cellInfo$CellType)) {
  pdf(paste0(outDir, '02.RSS_plot_',i,'.pdf'),width = 8, height = 8)
  print(plotRSS_oneSet(rss, setName = i) )
  dev.off()
}

# 3.Regulon
cellsPerGroup <- split(rownames(cellInfo), cellInfo[, selectedResolution]) 
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(x) 
                                    rowMeans(getAUC(regulonAUC)[,x]))
range(regulonActivity_byGroup)

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 


dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)
pdf(paste0(outDir,'03.heatmap.pdf'),width = 10, height = 20)
print(Heatmap(
  regulonActivity_byGroup_Scaled,
  name = "z-score",
  col = colorRamp2(seq(from=-2, to=2, length=11),
                          rev(brewer.pal(11, "Spectral"))),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 4),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows  = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE
))
dev.off()

# 4.top3 Regulon
rss <- regulonActivity_byGroup_Scaled 
head(rss)
df <- do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss), 
                 cluster = colnames(rss)[i], 
                 sd.1 = rss[,i], 
                 sd.2 = apply(rss[,-i], 1, median)
               )
             }))
df$fc = df$sd.1 - df$sd.2

top10 <- df %>% 
  group_by(cluster) %>% 
  top_n(10, fc)
rowcn <- data.frame(path = top10$cluster) 
rss_top10_percells <- rss[top10$path,] 

breaksList = seq(-1.5, 1.5, by = 0.1)
colors <- colorRampPalette(c("#336699", "white", "tomato"))(length(breaksList))
pdf(paste0(outDir, "04.Top10_regulon_scaledAUC_heatmap.pdf"), width = 12, height = 12)
print(pheatmap(rss_top10_percells,
               color = colors,
               cluster_rows = F,
               cluster_cols = FALSE,
               show_rownames = T,
               #gaps_col = cumsum(table(annCol$Type)),  
               #gaps_row = cumsum(table(annRow$Methods)),
               fontsize_row = 12,
               fontsize_col = 12,
               annotation_names_row = FALSE))
dev.off()


sce <- Chol.seu
sce@meta.data<- cbind(sce@meta.data, t(regulonAUC@assays@data@listData[["AUC"]])) 
Idents(sce) <- sce$SCT_snn_harmony_t.0.0.1
p <- DotPlot(sce, features = unique(top3$path))+   
#     theme_bw()+  
     theme(plot.background=element_blank(),
           panel.background=element_blank(),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
           axis.text.x=element_text(size=14, angle=30, hjust=1, vjust=1),
           axis.text.y=element_text(size=14),
           axis.title.x = element_text(size = 20), 
           axis.title.y = element_text(size = 10))
p
ggsave(paste0(outDir, '05.Top3_regulon_dotplot.pdf'), width = 15, height = 6)

### Dotplot_2
p <- DotPlot(sce,features = unique(top3$path))
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)
exp$id <- as.factor(exp$id)
exp$id <- fct_inorder(exp$id)
ggplot(exp, aes(x=id, y= features.plot))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp.scaled))+
  geom_point(aes(size=pct.exp, color=avg.exp.scaled),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black", size = 1),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=11,color="black",angle=30, hjust=1, vjust=1), 
        axis.text.y=element_text(size=11,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)+coord_flip()
ggsave(paste0(outDir, '06.Top3_regulon_dotplot2.pdf'), width = 24, height = 12)

# Clustered_DotPlot
library(scCustomize) 
library(viridis)
library(gridExtra)
colnames(sce@meta.data)[1:60]
length(colnames(sce@meta.data))
pdf(paste0(outDir, '07.All_regulon_clustered_dotplot.pdf'),height = 40, width = 12)
Clustered_DotPlot(seurat_object = sce, 
                  colors_use_exp = c('gray','red'),
                  print_exp_quantiles = F,
                  features = colnames(sce@meta.data)[59:327],
                  cluster_ident = F)
dev.off()

pdf(paste0(outDir, '08.Top3_regulon_clustered_dotplot.pdf'), height = 8, width = 6)
Clustered_DotPlot(seurat_object = sce, 
                  colors_use_exp = c('gray','red'),
                  print_exp_quantiles = F,
                  features = unique(top3$path),
                  cluster_ident = T)
dev.off()

# get the top10 TF in cluster 3
library(readr) 
regulonTargetsInfo <- read_tsv('output/Step2_regulonTargetsInfo.tsv')
cluter3_top10 <- data.frame(top10[top10$cluster==3, 'path'])
cluter3_top10$path_rm <- trimws(gsub("\\(.*\\)", "", cluter3_top10$path))
cluter3_top10_regulonTargetsInfo <- regulonTargetsInfo[which(regulonTargetsInfo$TF %in% cluter3_top10$path_rm), ]
write_tsv(cluter3_top10_regulonTargetsInfo, file = paste0(outDir, 'cluter3_top10_regulonTargetsInfo.tsv'))


