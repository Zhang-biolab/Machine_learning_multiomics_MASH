library(data.table)
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

# Goal:
# use expression of marker genes per cluster to manually annotate the Sep2024
# dataset with celltypes

# Functions ----
mega_melt <- function(df, source) {
  df.m <- melt(df, measure.vars=names(df))
  names(df.m) = c('cellType', 'marker.gene')
  df.m$source = source
  
  return(df.m)
}

update_annotation <- function(annotation, meta, name, clusterNums, clusterCol) {
  update.rows <- meta[[clusterCol]] %in% clusterNums
  annotation[update.rows] <- name
  return(annotation)
}

# MAIN ----

# path to preprocessed seurat object, (produced by do_preprocessing.sh)
seuObjDir <- '20250509_output/'
seuObjPath <- paste0(seuObjDir, 'Aggr_31Samples_Sep_2024_pre-processed.rds')
outObjFile <- 'Aggr_31Samples_Sep_2024_annotated.rds'
seuObjOutPath <- paste0(seuObjDir, outObjFile)

# path to plot output directory
plotDir <- '20250509_output/cellType_annotation/'
dir.create(plotDir)

#1. load the Seurat object
seu <- readRDS(seuObjPath)
DefaultAssay(seu) <- 'SCT'

tmp.seu <- RunPCA(seu, npcs = 30)
correlation <- cor(tmp.seu$nCount_RNA, tmp.seu@reductions$pca@cell.embeddings[,1])

## 2. Load Marker Genes
markerGeneDir <- 'Marker_genes_sets/Vallier_NAFLD-Gribben_Analysis_Markers-expression/'

# now ONLY use the most discriminative markers Ruben previously identified, 
# and the B-cell ones
all_the_markers <- fread(paste0(markerGeneDir, "marker_genes_v7.txt"), sep = "\t")
am_m <- mega_melt(all_the_markers, 'AM')
am_m$cellType <- as.character(am_m$cellType)
am_m$cellType[am_m$cellType=='Stellate cells'] <- 'Stellate'

all.markers <- am_m
rm(all_the_markers, am_m)

# 3. tidy up all.markers:
# remove duplicates
all.markers <- all.markers[!(duplicated(all.markers[, c('cellType', 'marker.gene')])), ]
# remove blank genes (from the csv format)
all.markers <- all.markers[all.markers$marker.gene != '']
all.markers <- all.markers[-which(all.markers$marker.gene %in% c('RAYL')), ]

# remove any putative markers not detected in the data
cluster_counts <- table(all.markers$cellType)[rev(unique(all.markers$cellType))]
clusters <- names(cluster_counts)
n_clusters <- length(clusters)

breaks <- cumsum(cluster_counts) + 0.5
breaks <- head(breaks, -1)  

### plot expression of marker genes by the clusters of cells identified in pre-processing
resns = c(2.0)
r = 2.0
for (r in resns) {
  clustercol = paste0('SCT_snn_res.', r)
  
  print(clustercol)
  p.umap <- DimPlot(seu, reduction='umap', group.by=clustercol,
                    label = TRUE,
                    label.size = 4)
  p.umap.split <- DimPlot(seu, reduction='umap', 
                          group.by=clustercol,split.by = clustercol, ncol=7,
                          label = TRUE,
                          label.size = 4)
  p.heatmap <- DoHeatmap(seu, features=all.markers$marker.gene, group.by=clustercol) + 
    geom_hline(
      yintercept = breaks, 
      linetype = "dashed", 
      color = "red",        
      linewidth = 0.5        
    )
  
  ggsave(paste0(plotDir,
                'marker_umap_', clustercol, '.pdf'), 
         plot = p.umap, 
         width=20, height=20)
  ggsave(paste0(plotDir, 
                'marker_umap_split_', clustercol, '.pdf'),
         plot=p.umap.split,
         width=20, height=20)
  ggsave(paste0(plotDir, 'marker_heatmap_', clustercol, '.pdf'),
         plot=p.heatmap,
         width=30, height=30)
}

# UMAPs with expression of the markers from 1E, mainly look at the expression of B cell markers
marker.genes <- c('IGKC', 'MZB1', # B cell1
                  'BANK1', 'PAX5', 'FCRL5', 'MS4A1', 'EBF1', # B cells2
                  'TRAC', 'KLRB1', 'CD8A', 'CD3E', # T-cells
                  'CD247', 'CD2', 'COL3A1', 'DCN','MNDA', 
                  'FCN1', 'CD163', 'MARCO', 'CFTR', 'KRT7', 'STAB2', 'ASGR1', 'CYP3A4')

plist <- list()
mg <- marker.genes[3]
for (mg in marker.genes) {
  p <- FeaturePlot(seu, reduction='umap', features=mg, slot='data')
  # p <- p+scale_color_viridis_c()+
  p <- p+ xlab('UMAP1') +
    ylab('UMAP2')+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank())
  plist[[mg]] <- p
}

pdf(file=paste0(plotDir, 'Marker_UMAPs.pdf'), height=7, width=7)
for (mg in names(plist)) {
  print(plist[[mg]])
}
dev.off()


# based on the above plots, assign the clusters to cellTypes ----
# A FIRST PASS TO ASSIGN CELLS BASED ON CLUSTERS
# based on marker by cluster in "marker_heatmap_SCT_snn_res.2_annotated.pdf"
chol.clusters = c(24,26,37,40)
stellate.clusters = c(31,34,36)
endothelial.clusters = c(5,10,18,21,41)
lympho.clusters = c(17,23,33,39) # half of 30 is endothelial
macro.clusters = c(27,30) # half of 25 is neutrophils 
hep.clusters = sort(setdiff(c(0:43), 
                            c(chol.clusters, stellate.clusters, 
                              endothelial.clusters, lympho.clusters, macro.clusters))) # 16 half hep, ste, endo, lymp - will need refining
Bcell.clusters = c()
neutro.clusters = c()
unassigned.clusters = c()

# assign cellType annotations to seurat object
meta <- seu@meta.data
annotation = rep('notDone', nrow(meta))
annotation <- update_annotation(annotation, meta, 'Hepatocytes', hep.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Cholangiocytes', chol.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Stellate', stellate.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Endothelial', endothelial.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Lymphocytes', lympho.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Macrophages', macro.clusters, 'SCT_snn_res.2')
# annotation <- update_annotation(annotation, meta, 'B-cell', Bcell.clusters, 'SCT_snn_res.2')
# annotation <- update_annotation(annotation, meta, 'Neutrophil', neutro.clusters, 'SCT_snn_res.2')
# unannotated labelled with clusters
annotation[annotation=='notDone'] <- as.character(seu@meta.data$SCT_snn_res.2[annotation=='notDone'])
#annotation <- update_annotation(annotation, meta, 'unknown', unassigned.clusters, 'SCT_snn_res.2')

# add the annotation to the seu obj
seu <- AddMetaData(seu, 'Annot.Sep2024.tmp', col.name='cell.annotation.version')
seu <- AddMetaData(seu, annotation, col.name='cell.annotation')

unique(seu@meta.data[, c('cell.annotation', 'SCT_snn_res.2')])

# sanity plot how it looks currently
# sanity plot UMAP 
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation')
p2 <- DimPlot(seu, reduction='umap', group.by='SCT_snn_res.2', label = TRUE, label.size = 4)
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation', split.by='cell.annotation', ncol=3)
p.both <- DimPlot(seu, reduction='umap', group.by='SCT_snn_res.2', split.by='cell.annotation', ncol=3,
                  label = TRUE, label.size = 4)

ggsave(paste0(plotDir, 'cellType_UMAP.pdf'), 
       plot=p, 
       width=10, height=10)
ggsave(paste0(plotDir, 'cellType_UMAP_facetted.pdf'), 
       plot=psplit, 
       width=20, height=20)


# mainly refining the cluster according to following results
p.list <- list()
for (ct in unique(seu@meta.data$cell.annotation)) {
  Idents(seu) <- seu@meta.data$cell.annotation
  seu.sub <- subset(seu, idents=ct)
  p <- DimPlot(seu.sub, reduction='umap', group.by='SCT_snn_res.2', label = TRUE, label.size = 4)+
    ggtitle(ct)
  p.list[[ct]] <- p
}
p.all <- plot_grid(plotlist=p.list)
ggsave(paste0(plotDir, 'cellType_UMAPgroup&celltype.pdf'), 
       plot=p.all, 
       width=15, height=10)


# tidy the clusters that look like have multiple celltypes:
meta <- seu@meta.data
umap.embeddings <- seu@reductions$umap@cell.embeddings
expr <- seu@assays$SCT@scale.data

# fixing Cholangiocytes
meta$cell.annotation[meta$cell.annotation=='Cholangiocytes' &
                       (umap.embeddings[, 2] > -5)] <- 'Hepatocytes'

meta$cell.annotation[meta$cell.annotation=='Cholangiocytes' &
                       (umap.embeddings[, 2] < -10)] <- 'Lymphocytes'


# fixing Endothelial
meta$cell.annotation[meta$cell.annotation=='Endothelial' &
                       (umap.embeddings[, 1] < 5)] <- 'Hepatocytes'

meta$cell.annotation[meta$cell.annotation=='Endothelial' &
                       (umap.embeddings[, 1] > 5 &
                        umap.embeddings[, 2] < -5)] <- 'Hepatocytes'



# fixing lymphocytes
meta$cell.annotation[meta$cell.annotation=='Lymphocytes' &
                       (umap.embeddings[, 2] > -10 &
                          umap.embeddings[, 2] < -5)] <- 'Macrophages'

meta$cell.annotation[meta$cell.annotation=='Lymphocytes' &
                     umap.embeddings[, 2] > -5] <- 'Hepatocytes'


meta$cell.annotation[meta$cell.annotation=='Lymphocytes' &
                       umap.embeddings[, 1] < -6] <- 'B-cell 1'


meta$cell.annotation[meta$cell.annotation=='Lymphocytes' &
                       (umap.embeddings[, 1] > 2 &
                          umap.embeddings[, 2] < -10)] <- 'B-cell 2'

# fixing Macrophages
meta$cell.annotation[meta$cell.annotation=='Macrophages' &
                       (umap.embeddings[, 2] > -6)] <- 'Hepatocytes'

meta$cell.annotation[meta$cell.annotation=='Macrophages' &
                       (umap.embeddings[, 1] > 5)] <- 'Cholangiocytes'

meta$cell.annotation[meta$cell.annotation=='Macrophages' &
                       expr['FCN1',] > 0] <- 'Neutrophils'

# fixing Stellate
meta$cell.annotation[meta$cell.annotation=='Stellate' &
                       (umap.embeddings[, 1] < 5 &
                          umap.embeddings[, 2] > -2)] <- 'Hepatocytes'

meta$cell.annotation[meta$cell.annotation=='Stellate' &
                       (umap.embeddings[, 1] > 5 &
                          umap.embeddings[, 2] > -5)] <- 'Cholangiocytes'

meta$cell.annotation[meta$cell.annotation=='Stellate' &
                       (umap.embeddings[, 1] > 5 &
                          umap.embeddings[, 2] < -5)] <- 'Endothelial'

# fixing Hepatocytes 
meta$cell.annotation[meta$cell.annotation=='Hepatocytes' &
                       (umap.embeddings[, 1] > 5 &
                        umap.embeddings[, 2] < -5)] <- 'Cholangiocytes'

meta$cell.annotation[meta$cell.annotation=='Hepatocytes' &
                       (umap.embeddings[, 2] < -10)] <- 'Lymphocytes'

meta$cell.annotation[meta$cell.annotation=='Hepatocytes' &
                       (umap.embeddings[, 1] > 5 &
                          umap.embeddings[, 2] > -2.5 &
                          umap.embeddings[, 2] < 5)] <- 'Endothelial'


seu <- AddMetaData(seu, meta$cell.annotation, col.name='cell.annotation.refined')

# sanity plot how it looks currently
# sanity plot UMAP
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation.refined')
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation.refined', split.by='cell.annotation.refined', ncol=3)
p.both <- DimPlot(seu, reduction='umap', group.by='SCT_snn_res.2', split.by='cell.annotation.refined', ncol=3)

ggsave(paste0(plotDir, 'cellType_UMAP_refined.pdf'),
       plot=p,
       width=10, height=10)
ggsave(paste0(plotDir, 'cellType_UMAP_facetted_refined.pdf'),
       plot=psplit,
       width=20, height=20)

# check how many cells of each type identified
table(seu@meta.data$cell.annotation.refined)

# only keep the refined annotation
seu@meta.data$cell.annotation <- seu@meta.data$cell.annotation.refined
seu@meta.data$cell.annotation.refined <- NULL

# sanity plot UMAP 
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation')
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation', split.by='cell.annotation', ncol=3)

# 5. save the cellType annotated Seuarat object
saveRDS(seu,
        file=seuObjOutPath)
