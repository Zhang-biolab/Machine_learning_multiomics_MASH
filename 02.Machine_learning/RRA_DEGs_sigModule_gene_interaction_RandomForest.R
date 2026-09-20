library(randomForest)    
set.seed(12345)

#Read input file
datExpr <- read.table('RRA_SVA_WGCNA_datExpr.txt' , header=T, sep="\t", check.names=F, row.names=1)
datExpr <- t(datExpr)

#Model construct
x <- as.matrix(datExpr)
y <- rownames(datExpr)
y[grep('Mild', y)] <- 'Mild_Fibrosis'
y[grep('Adv',y)] <- 'Advanced_Fibrosis'
y <- factor(y, levels = c('Mild_Fibrosis', 'Advanced_Fibrosis'))

rf <- randomForest(y~., data=x, ntree=500)

pdf(file="RandomForest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

optionTrees <- which.min(rf$err.rate[, 1])
rf2 <- randomForest(y~., data=x, ntree=optionTrees)

pdf(file="RandomForest_geneImportance.pdf", width=6.2, height=7)
varImpPlot(rf2, main="")
dev.off()

importance_data <- importance(rf2) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Variable") %>% 
  arrange(desc(MeanDecreaseGini))

rfGenes <- importance_data[importance_data$MeanDecreaseGini > 2, "Variable"]

write.table(rfGenes, file="RandomForest_genes.txt", sep="\t", quote=F, col.names=F, row.names=F)

# visualization
plotDir <- "~/projects/NASH_Fibrosis_datamining/manuscript_figure/ML_GSVA_ROC_prediction/"

# errorRate point plot
err_rate <- rf$err.rate

err_df <- data.frame(
  Trees = 1:nrow(err_rate),
  err_rate
)

err_long <- melt(err_df, id.vars = "Trees", variable.name = "Error_Type", value.name = "Error_Rate")

size <- 18
theme <- theme(plot.title= element_text(size=18, color="black", hjust = 0.5),
               axis.title.x=element_text(color = "black", size=size),
               axis.title.y=element_text(color= "black", size=size),
               axis.text.x=element_text(color = "black", size=size),
               axis.text.y=element_text(color = "black", size=size),
               axis.line = element_line(colour = "black"),
               axis.ticks=element_line(colour="black", size=0.25, linetype=1, lineend=1),
               legend.title = element_text(colour="black", size=size),
               legend.text = element_text(color="black", size = size),
               legend.position='none')

ggplot(err_long, aes(x = Trees, y = Error_Rate, color = Error_Type)) +
  geom_point(size=0.3) + 
  geom_line(size = 0.3) +   
  scale_color_manual(
    values = c("OOB" = "#090905", "Mild_Fibrosis" = "#c85367", "Advanced_Fibrosis" = "#78bb51")) +
  labs(title = "Random forest", x = "trees", y = "Error", color = "Error type") +
  theme_classic() +
  theme
ggsave(paste0(plotDir, "RandomForest_trees_errorRate_point_plot.pdf"), height=5, width=5)


# 2. MeanDecreaseGini barplot
importance_data <- importance(rf2) %>% 
  as.data.frame() %>% 
  filter(MeanDecreaseGini > 2) %>%
  tibble::rownames_to_column("Gene") %>% 
  dplyr::select(Gene, MeanDecreaseGini) %>% 
  arrange(desc(MeanDecreaseGini)) %>%  
  mutate(Gene = fct_reorder(Gene, MeanDecreaseGini)) 

ggplot(importance_data, aes(x=MeanDecreaseGini, y=Gene)) +
  geom_segment(aes(yend=Gene), xend=0, colour='grey50') +
  geom_point(colour="#1E90FF", size=5) +
  labs(title = "", x = "Importance", y = "") +
  xlim(c(0, max(importance_data$MeanDecreaseGini))) +
  theme_classic() +
  theme
ggsave(paste0(plotDir, "RandomForest_importance_point_plot.pdf"), height=7.5, width=6)
