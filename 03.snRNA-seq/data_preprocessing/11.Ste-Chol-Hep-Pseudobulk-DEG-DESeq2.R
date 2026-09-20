library(Seurat)
library(data.table)
library(dplyr)
library(DESeq2)
library(tibble)
library(dplyr)
library(tidyr)
library(ggrepel)
library(stringr)

setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

inFile <- "20250509_output/Aggr_Sep2024_filt_All_harmony_th=0.rds"
outDir <- '20250509_output/Pseudobulk_DEG_DESeq2_with_sex_chromosome_genes/'
dir.create(outDir, recursive = T)

# Load the MASH atlas
seu <- readRDS(inFile)
DefaultAssay(seu)

# seu <- readRDS("Human_MASH_Atlas.rds")
metadata <- seu@meta.data

pseudobulk_DESeq2 <- function(subset, Entity, reference) {
  cts <- AggregateExpression(subset, group.by = Entity, slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  colData <- data.frame(samples = colnames(cts))
  
  # Extract the "condition" column
  colData <- colData %>% mutate(condition = str_split(samples, "_", simplify = TRUE)[, 1])
  colData$condition <- sub("_[^_]+$", "", colData$samples)
  rownames(colData) <- colData$samples
  colData$samples <- NULL

  dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = reference)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds)
  
  result_col <- resultsNames(dds)[2]
  print(result_col <- resultsNames(dds)[2])
  res <- results(dds, name = result_col)
  res <- data.frame(res)
  res$Gene <- rownames(res)
  
  return(res)
}

#Process the pseudobulk DEG for the Stellate in respect to the MASH fibrosis stage
celltypes <- c("Stellate", "Cholangiocytes", "Hepatocytes")
for (ct in celltypes) {
  subset <- subset(seu, cell.annotation == ct) # P71 without Stellate
  res_Disease <- pseudobulk_DESeq2(subset, c("Fibrosis.stage", "Patient.ID"), "Mild.fibrosis..F0.2.")
  write.table(res_Disease, file=paste0(outDir, ct, "_pseudobulk_DEG_DESeq2_disease_status.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  p_thresh <- 0.05
  fc_thresh <- log2(1.5)
  data <- res_Disease %>%
  mutate(threshold = case_when(log2FoldChange > fc_thresh & pvalue < p_thresh ~ "Upregulated",
                               abs(log2FoldChange) < fc_thresh | pvalue > p_thresh ~ "Not significant",
                               log2FoldChange < -fc_thresh & pvalue < p_thresh ~ "Downregulated")) %>%
  mutate(threshold_factor = factor(threshold, levels=c("Not significant","Upregulated","Downregulated"))) %>%
  dplyr::select(Gene, log2FoldChange, pvalue, padj, threshold, threshold_factor)
  data$neg_logP <- -log10(data$pvalue)
  data <- data[order(data$threshold_factor), ] 
  if (sum(data$neg_logP != Inf) == 0) {
    data$neg_logP <- 500
    } else {
      data$neg_logP[data$neg_logP == Inf] <- max(data$neg_logP[data$neg_logP != Inf]) + 10
      }

  topUp <- data %>%
  filter(threshold_factor=='Upregulated') %>%
  distinct(Gene, .keep_all = T) %>% 
  arrange(desc(neg_logP)) %>%
  slice_head(n = 10) 
  
  topDown <- data %>%
  filter(threshold_factor=='Downregulated') %>%
  distinct(Gene, .keep_all = T) %>% 
  arrange(desc(neg_logP)) %>%
  slice_head(n = 10)
  
  mycol <- c("grey", "red","blue")
  p1 <- ggplot(data=data, aes(x=log2FoldChange, y=neg_logP, colour=threshold_factor))+
  scale_color_manual(values=alpha(mycol, 0.7)) +
  geom_point(size=1.5, pch=19) +
  ylim(c(0, max(data$neg_logP)+10)) +
  geom_vline(xintercept=c(-fc_thresh, fc_thresh), lty="dashed", col="black", size=0.7)+ #FC的阈值是-2，2
  geom_hline(yintercept = -log10(0.05), lty="dashed", col="black", size=0.7) +
  theme_classic() +
  scale_y_sqrt()+
  labs(x=expression("log"["2"]*"FC"), y=expression("-log"["10"]*" (p value)")) 

  p2 <- p1 +
  geom_text_repel(data = topUp,
                  aes(x = log2FoldChange, y = neg_logP, label = Gene),
                  color = "red",
                  seed = 233,
                  size = 4,
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black', 
                  segment.alpha = 0.5, 
                  nudge_x = 1.5 - topUp$log2FoldChange, 
                  direction = "y", 
                  hjust = 0 
  )

  p3 <- p2 +
  geom_text_repel(data = topDown,
                  aes(x = log2FoldChange, y = neg_logP, label = Gene),
                  color = "blue",
                  seed = 233,
                  size = 4,
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.8,
                  nudge_x = -1.5 - topDown$log2FoldChange,
                  direction = "y",
                  hjust = 1 
  ) 

  size <- 20
  p4 <- p3 +
  theme(plot.title = element_text(),
        axis.text.x= element_text(color="black", size=size),
        axis.text.y = element_text(color="black", size=size),
        axis.title.x = element_text(color="black",size=size),
        axis.title.y = element_text(color="black",size=size),
        legend.position="none")
  ggsave(paste0(outDir, "pseudobulkAnalysis_", ct, "_Volcano.pdf"), p4, width = 7, height = 6.8)
}