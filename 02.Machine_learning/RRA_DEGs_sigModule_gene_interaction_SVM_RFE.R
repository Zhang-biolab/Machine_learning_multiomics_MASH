library(e1071)
source("~/projects/NASH_Fibrosis_datamining/src/20231207/Function/msvmRFE_V2.R")
set.seed(12345)

#Read input file
datExpr <- read.table('RRA_SVA_WGCNA_datExpr.txt' , header=T, sep="\t", check.names=F, row.names=1)
datExpr <- as.data.frame(t(datExpr))

#Obtain sample grouping information
group <- rownames(datExpr)
group[grep('Mild', group)] <- 'Mild_Fibrosis'
group[grep('Adv',group)] <- 'Advanced_Fibrosis'

inputDat <- cbind(group, datExpr)
inputDat$group <- factor(inputDat$group, levels = c('Mild_Fibrosis', 'Advanced_Fibrosis'))

#Construct Machine Learning-Support Vector Machine Recursive Feature Elimination Algorithm (SVM-RFE)
svmRFE(inputDat, k=10, halve.above=50)
nfold <- 10
nrows <- nrow(inputDat)
folds <- rep(1:nfold, len=nrows)[sample(nrows)]
folds <- lapply(1:nfold, function(x) which(folds == x))

# perform the feature ranking for all 10 training sets
results <- lapply(folds, svmRFE.wrap, inputDat, k=10, halve.above=50)

#Rank the importance of the characteristic genes
top.features <- WriteFeatures(results, inputDat, save=F)
#Output the sorted result
write.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

# Estimate generalization error using a varying number of top features
# This process is repeated while varying the number of top features that are used as input, 
# and there will typically be a "sweet spot" where there are not too many nor too few features.
featsweep <- lapply(1:74, FeatSweep.wrap, results, inputDat) 
save(featsweep,file = "featsweep.RData")

#Get the error of cross validation
no.info <- min(prop.table(table(inputDat[,1])))
errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#Plot cross-validation errors
pdf(file="SVM-REF_errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

#Draw graphs that cross-verify accuracy
pdf(file="SVM-REF_accuracy.pdf", width=5, height=5)
PlotAccuracy(1-errors, no.info=no.info) # the PlotAccuracy function was not included in msvmRFE.R, we should add it additionally. 
dev.off()

#The characteristic genes of SVM were output
featureGenes <- top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

