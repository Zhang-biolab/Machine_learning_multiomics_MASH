# Quick pre-processing of a Seurat object. 
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(data.table)
library(stringr)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")


harmony_and_reprocess <- function(seu, theta, reduction.name, cluster.resns) {
  # Run Harmony --- 
  seu <- RunHarmony(seu, 
                    group.by.vars = 'Patient.ID',
                    reduction='pca',
                    theta=theta, 
                    reduction.save=reduction.name)
  
  seu <- FindNeighbors(seu, 
                       dims=1:30, 
                       k.param=20,
                       reduction = reduction.name) 
  
  
  # UMAP - runs on the reduced dims directly - not on neighborhood graph
  seu <- RunUMAP(seu, 
                 reduction=reduction.name,
                 dims = 1:30, 
                 n.neighbors = 20, 
                 min.dist = 0.3,
                 reduction.name=paste0('umap_', reduction.name),
  )
  
  # Clusters
  for (c.res in cluster.resns) {
    # uses the Neighbors last found
    seu <- FindClusters(seu, 
                        resolution = c.res, 
                        alg = 1 # Louvain algorithm
    )
    seu <- AddMetaData(seu, 
                       seu@meta.data$seurat_clusters, 
                       col.name = paste0('SCT_snn_', reduction.name, '.', c.res))
  }
  
  return(seu)
}

### INPUT FILES-----
inFile <- "20250509_output/Aggr_31Samples_Sep_2024_annotated.rds"
outDir <- '20250509_output/'
cellStr <- 'All' # 'All', 'Chol', 'Hep' 'Chol-Hep'
plotDir <- paste0(outDir, 'harmony_plots_', cellStr, '/') # plots
dir.create(plotDir, recursive = T)
theta <- 0 
outFile <- paste0('Aggr_Sep2024_filt_', cellStr, '_harmony_th=', theta, '.rds')

seu <- readRDS(inFile)

# ElbowPlot(seu, ndims = 50)

# clust.resns <- c(0.4)
clust.resns <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.6)

reduction.name <- paste0('harmony_t.', theta)

seu <- harmony_and_reprocess(seu, theta, reduction.name,
                             clust.resns)

disease.order <- c('Mild.fibrosis..F0.2.',
                   'Advanced.fibrosis..F3.4.')

seu@meta.data$Fibrosis.stage <- factor(seu@meta.data$Fibrosis.stage,
                                       levels=disease.order)

gc()
print('about to save SeuObj:')
print(seu)

saveRDS(seu, 
        file=paste0(outDir, outFile))