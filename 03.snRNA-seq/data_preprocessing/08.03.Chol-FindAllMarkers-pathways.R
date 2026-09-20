library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)
library(dplyr)
library(broom)
library(Matrix)
library(lme4)
library(cowplot)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

source('~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/scripts/figure_scripts/plot_functions.R')
source('~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/scripts/figure_scripts/do_DGE_functions.R')
EnrichR.gseapy <- '~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/scripts/figure_scripts/EnrichR_gseapy.py'

inFile <- "20250509_output/Subset_Chol/Aggr_Sep2024_filt_Chol_harmony_th=0.rds"
outDir <- '20250509_output/Subset_Chol/'
cellStr <- 'Chol' # 'All', 'Chol', 'Chol' 'Chol-Chol'
plotDir <- paste0(outDir, 'FindAllMarkers_plots_', cellStr, '/') # plots
dir.create(plotDir, recursive = T)

# see powerpoint figure 3 explanations for what Chris has requested for this
# "slides" refer to this powerpoint ---

# summary plots of the harmony reduction
theta <- 0
reduction.name <- paste0('harmony_t.', theta)

seu <- readRDS(inFile)
aggregate(seu@meta.data$SCT_snn_harmony_t.0.0.2, by=list(seu@meta.data$Fibrosis.stage), FUN=table)
aggregate(seu@meta.data$SCT_snn_harmony_t.0.0.2, by=list(seu@meta.data$Fibrosis.score), FUN=table)

p1 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Fibrosis.stage', cols = c("#4C8BBD", "#D7282C"))
p2 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Patient.ID')
p3 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='orig.ident')
p4 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Fibrosis.score..F0.4.')
clust.to.plot <- paste0('SCT_snn_harmony_t.', theta, '.0.2') # plot clustering for current theta, c.res=0.1
p5 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by=clust.to.plot, label = TRUE)
p6 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Disease.status')
p.all <- plot_grid(p1, p2, p3, p4, p5, p6, ncol=2)
ggsave(file=paste0(plotDir, 'summary_UMAPs_theta=', theta, '.pdf'), p.all,
       width=20, height=10)
ggsave(file=paste0(plotDir, 'fig1_groups_UMAPs_theta=', theta, '.pdf'), p1,
       width=6, height=4)
ggsave(file=paste0(plotDir, 'fig2_clusters_UMAPs_theta=', theta, '.pdf'), p5,
       width=6, height=4)
ggsave(file=paste0(plotDir, 'fig3_fibroscore_UMAPs_theta=', theta, '.pdf'), p4,
       width=6, height=4)

p7 <- DimPlot(seu,
              reduction=paste0('umap_', reduction.name),
              group.by=paste0('SCT_snn_harmony_t.', theta, '.0.2'), 
              split.by=c('Fibrosis.stage'),
              label=TRUE,
              label.size=6)
ggsave(file=paste0(plotDir, 'fig4_fibrosisStage_Split_UMAPs_theta=', theta, '.pdf'), p7,
       width=6, height=4)


p8 <- DimPlot(seu,
              reduction=paste0('umap_', reduction.name),
              group.by=paste0('SCT_snn_harmony_t.', theta, '.0.2'), 
              split.by=c('Fibrosis.score..F0.4.'),
              label=TRUE,
              label.size=6)
p8

p9 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), 
               group.by='Fibrosis.stage', 
               split.by='Patient.ID', 
               ncol=7)

p9

########################################################################################
# Cells split by fibrosis stage
seu <- readRDS(inFile)
m <- seu@meta.data
harmony_reduction <- 'umap_harmony_t.0'
aggregate(seu@meta.data$SCT_snn_harmony_t.0.0.2, by=list(seu@meta.data$Fibrosis.stage), FUN=table)
table(seu@meta.data$SCT_snn_harmony_t.0.0.2)
Idents(seu) <- m$SCT_snn_harmony_t.0.0.2

# sort out some factor levels
seu@meta.data$Fibrosis.score..F0.4. <- as.factor(seu@meta.data$Fibrosis.score..F0.4.)

disease.order <- c('Mild.fibrosis..F0.2.',
                   'Advanced.fibrosis..F3.4.')

seu@meta.data$Fibrosis.stage <- factor(seu@meta.data$Fibrosis.stage,
                                       levels=disease.order)

# The number of nuclei recovered per sample ------
m$barcode <- rownames(m)
m <- data.table(m)

nCells.group <- m[, .('nCells'=.N), by=.(Patient.ID, SCT_snn_harmony_t.0.0.2, Fibrosis.stage)]
aggregate(nCells ~ SCT_snn_harmony_t.0.0.2+Fibrosis.stage, nCells.group, sum) 

fibrosis.cols <- c("#4C8BBD", "#D7282C")
p <- ggplot(nCells.group, aes(x=SCT_snn_harmony_t.0.0.2, weight=nCells, fill=Fibrosis.stage))+
  geom_bar(position="fill")+
  scale_fill_manual(values=fibrosis.cols) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.text.x=element_text(vjust=1),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")
p

ggsave(paste0(plotDir, 'fig5_cluster_splitbyfibrosis.pdf'), device = 'pdf',  height=5, width=6)

################################################################################
## test %change for a cell_type 
seu <- readRDS(inFile)
table(seu@meta.data$SCT_snn_harmony_t.0.0.2)

# sort out some factor levels
seu@meta.data$Fibrosis.score..F0.4. <- as.factor(seu@meta.data$Fibrosis.score..F0.4.)

disease.order <- c('Mild.fibrosis..F0.2.',
                   'Advanced.fibrosis..F3.4.')

seu@meta.data$Fibrosis.stage <- factor(seu@meta.data$Fibrosis.stage,
                                       levels=disease.order)
m <- seu@meta.data

mdata = data.table(m, keep.rownames = T)

mdata[, 'N' := .N, by = list(Patient.ID, Fibrosis.stage)] # 按照Patient.ID分组，计算每个患者的总细胞数

mdata1 = mdata[Fibrosis.stage == 'Mild.fibrosis..F0.2.']
mdata2 = mdata[Fibrosis.stage == 'Advanced.fibrosis..F3.4.']
mdata1[, 'n' := .N, by = list(Patient.ID, SCT_snn_harmony_t.0.0.2)]
mdata2[, 'n' := .N, by = list(Patient.ID, SCT_snn_harmony_t.0.0.2)]

mdata = rbind(mdata1, mdata2)
mdata[, 'frac' := n/N]
mdata = subset(mdata, select = c(Fibrosis.score, Fibrosis.stage, SCT_snn_harmony_t.0.0.2, frac, n, N, Patient.ID)) %>% unique()
# mdata = mdata[mdata$N > 80, ]
# mdata = mdata[-which(mdata$Fibrosis.score %in% c(2)), ]
# PID1 = unique(mdata1$Patient.ID)
# PID2 = unique(mdata2$Patient.ID)
PID1 = unique(mdata[Fibrosis.stage=='Mild.fibrosis..F0.2.']$Patient.ID)
PID2 = unique(mdata[Fibrosis.stage=='Advanced.fibrosis..F3.4.']$Patient.ID)

proportions.list = list()
for(cl0 in unique(mdata$SCT_snn_harmony_t.0.0.2)){
  mdata0 = mdata[SCT_snn_harmony_t.0.0.2 == cl0]
  setkey(mdata0, Patient.ID)
  mdata0[, nSample := .N, by = Patient.ID]
  mdata0 = subset(mdata0, select = c(Fibrosis.stage, frac, Patient.ID))
  mdata01 = mdata0[ Fibrosis.stage == 'Mild.fibrosis..F0.2.']
  mdata02 = mdata0[ Fibrosis.stage == 'Advanced.fibrosis..F3.4.']
  y = mdata02[PID2]$frac
  x = mdata01[PID1]$frac
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  pdata = data.table('Patient.ID' = c(PID1, PID2), 
                     'SCT_snn_harmony_t.0.0.2' = cl0,  
                     'Fibrosis.stage' = rep(c('Mild.fibrosis..F0.2.', 'Advanced.fibrosis..F3.4.'), c(length(PID1), length(PID2))),
                     'proportion' = c(x, y)
  )
  proportions.list[[cl0]] = pdata
  
}

proportions <- do.call(rbind, proportions.list)

proportions$SCT_snn_harmony_t.0.0.2 <- factor(proportions$SCT_snn_harmony_t.0.0.2,
                                              levels(unique(mdata$SCT_snn_harmony_t.0.0.2)))

proportions$Fibrosis.stage <- factor(proportions$Fibrosis.stage,
                                     disease.order)

compute_p_values <- function(data) {
  data %>%
    group_by(SCT_snn_harmony_t.0.0.2) %>%
    do(tidy(wilcox.test(proportion ~ Fibrosis.stage, exact = FALSE, data = .))) # 运行tidy(t.test()or wilcox.test)需要加载broom包
}

compute_p_values <- function(data) {
  data %>%
    group_by(SCT_snn_harmony_t.0.0.2) %>%
    do(tidy(t.test(proportion ~ Fibrosis.stage, data = .))) # 运行tidy(t.test()or wilcox.test)需要加载broom包
}

# Get p-values
p_values <- compute_p_values(proportions)
p_values

format_pvalue <- function(pvalue) {
  if (pvalue < 0.001) {
    return(sprintf("p = %.1e", pvalue))
  } else {
    return(sprintf("p = %.3f", pvalue))
  }
}
p_values <- p_values %>%
  mutate(label = sapply(p.value, format_pvalue))

p <- ggplot(proportions, aes(x = SCT_snn_harmony_t.0.0.2, y = proportion, fill = Fibrosis.stage)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 0.5, alpha = 1) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) + 
  geom_vline(xintercept = seq(1.5, length(unique(proportions$SCT_snn_harmony_t.0.0.2)) - 0.5, by = 1), 
             linetype = "dotted", color = "gray") +
  scale_fill_manual(values = c("#4C8BBD", "#D7282C")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(colour = "black"),
        #legend.position = "none"
  ) +
  geom_text(data = p_values, aes(x = SCT_snn_harmony_t.0.0.2, y = max(proportions$proportion) + 0.1, label = label),
            vjust = 0, size = 3.5, inherit.aes = FALSE)
ggsave(paste0(plotDir, 'fig5_cellType_frac_byFibrosisStage_wilcox_test.pdf'), device = 'pdf', width = 10, height = 6)

######################################################################################################################
# do DGE of ste clusters vs all stes -------------------------
# DE parameters:
seu <- readRDS(inFile)
DefaultAssay(seu) <- "SCT"
harmony_reduction <- 'umap_harmony_t.0'

# 如果在 PrepSCTFindMarkers 之后你想提取子集然后寻找差异基因,需要加上 recorrect_umi = FALSE 参数
seu <- PrepSCTFindMarkers(seu, assay = "SCT", verbose = TRUE)
Idents(seu) <- 'SCT_snn_harmony_t.0.0.2'

all.markers <- FindAllMarkers(seu,
                              assay = "SCT", 
                              logfc.threshold=0.25,
                              test.use='wilcox',
                              min.pct=0.1,
                              only.pos = FALSE,
                              base=2,
                              recorrect_umi = FALSE) # base = which to use for logs
fwrite(file=paste0(plotDir, 'clusterSpecific_DE_table.csv'), 
       all.markers)

# Heatmap showing the overexpressed genes between each subcluster and all other subclusters.
top6Markers <- all.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC) %>% print(.,n=Inf)

p <- DoHeatmap(seu, features = top6Markers$gene, 
               group.by = "SCT_snn_harmony_t.0.0.2", 
               group.bar = T, size = 4) + NoLegend() # using scale.data slot for the SCT assay
ggsave(paste0(plotDir, 'fig6_top6_markers_heatmap.pdf'), height=12, width=15)

top2Markers <- all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% print(.,n=Inf)
p <- VlnPlot(seu, features = top2Markers$gene, pt.size = 0, ncol=2)
ggsave(paste0(plotDir, 'fig7_top2_markers_Vlnplot.pdf'), width = 4, height = 12)

p <- FeaturePlot(seu, features = top2Markers$gene, cols = c("grey", "red"), slot='data',
            reduction = harmony_reduction, min.cutoff = 'q2', max.cutoff = 'q98',ncol=2)
ggsave(paste0(plotDir, 'fig8_top2_markers_FeaturePlot.pdf'), width = 8, height = 12)

p <- DotPlot(seu, assay='SCT',
             features=top2Markers$gene, 
             group.by=clust.to.plot,
             scale=T)
p <- p+coord_flip()+
  xlab('')+
  ylab('')+
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p.nolab <- p + theme(legend.position = 'none')
ggsave(paste0(plotDir, 'fig9_top2_markers_Dotplot.pdf'), width = 4, height = 6)

###########################################################################################################
#1. Enrichment analysis for up regulated genes in cluster specific way
all.markers <- read.csv(file=paste0(wd, plotDir, 'clusterSpecific_DE_table.csv'), header = T, sep = ',') 
all.markers.up <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC > 0.25)

pathwaysDir <- paste0(wd, plotDir, 'ClusterSpecific_up_pathways')
dir.create(pathwaysDir, recursive = T)

for(i in unique(all.markers.up$cluster)){
  degs <- data.frame(gene = all.markers.up[all.markers.up$cluster == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}

#2. Enrichment analysis for down regulated genes in cluster specific way
all.markers <- read.csv(file=paste0(wd, plotDir, 'clusterSpecific_DE_table.csv'), header = T, sep = ',') 
all.markers.down <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC < -0.25)

pathwaysDir <- paste0(wd, plotDir, 'ClusterSpecific_down_pathways')
dir.create(pathwaysDir, recursive = T)

for(i in unique(all.markers.down$cluster)){
  degs <- data.frame(gene = all.markers.down[all.markers.down$cluster == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}
