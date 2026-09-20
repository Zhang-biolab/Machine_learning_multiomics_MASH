library(tidyverse)
library(glmnet) 
set.seed(12345)

#Read input file
datExpr <- read.table('RRA_SVA_WGCNA_datExpr.txt', header=T, sep="\t", check.names=F, row.names=1)
datExpr <- t(datExpr)

#Model construct
x <- as.matrix(datExpr)
y <- rownames(datExpr)
y[grep('Mild', y)] <- 'Mild_Fibrosis'
y[grep('Adv',y)] <- 'Advanced_Fibrosis'
y <- factor(y, levels = c('Mild_Fibrosis', 'Advanced_Fibrosis'))

fit <- glmnet(x, y, family = "binomial", alpha=1)
cvfit <- cv.glmnet(x, y, family="binomial", alpha=1, type.measure='deviance',nfolds = 10)

#Draw a graph of Lasso regression
pdf(file='lasso.pdf', width=6, height=5.5)
par(mar=c(5, 6, 4, 4), mgp = c(3, 1, 0),
  cex.lab=1.8, cex.axis=1.8)
plot(fit, xvar = 'lambda')
dev.off()

#Draw cross-validated graphs
pdf(file='cvfit.pdf',width=6,height=5.5)
par(mar=c(5, 6, 4, 4), mgp = c(3, 1, 0),
    cex.lab=1.8, cex.axis=1.8)
plot(cvfit, ylab="Binomial deviance")
dev.off()

#Output the screened feature genes
coef <- coef(fit, s=cvfit$lambda.min)
index <- which(coef != 0)
lassoGene <- row.names(coef)[index]
lassoGene <- lassoGene[-1]
write.table(lassoGene, file='LASSO.gene.txt', sep="\t", quote=F, row.names=F, col.names=F)
