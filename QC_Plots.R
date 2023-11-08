library(pheatmap)
library(ggplot2)
library(ggrepel)

args=commandArgs(trailingOnly=T)
counts.raw = as.character(args[1]) #Raw counts file from MAGeCK
counts.norm = as.character(args[2]) #Normalized counts file from MAGeCK

#POST-ALIGNMENT COUNT CORRELATION PLOTS
rawCounts = read.table(counts.raw)
normCounts = read.table(counts.norm)

rawCorrelation = cor(apply(rawCounts[,-1],2,as.numeric))
pheatmap(rawCorrelation,scale="none",cluster_cols=F,cluster_row=F,main="Correlation arranged by sample: raw counts")
normCorrelation = cor(apply(normCounts[,-1],2,as.numeric]))
pheatmap(normCorrelation,scale="none",cluster_cols=F,cluster_rows=F,main="Correlation arranged by sample: normalized counts")

#POST-ALIGNMENT PCA PLOTS
pca.raw = princomp(apply(rawCounts[,-1],2,as.numeric))
var.pct = pca.raw$sdev^2/(sum(pca.raw$sdev^2))
plot(pca$loadings,xlab=paste0("PC1 (",var.pct[1],"% variance)"), ylab=paste0("PC2 (",var.pct[],"% variance)"),main="PCA: raw counts")

pca.norm = princomp(apply(normCounts[,-1],2,as.numeric))
var.pct = pca.norm$sdev^2/(sum(pca.norm$sdev^2))
plot(pca$loadings,xlab=paste0("PC1 (",var.pct[1],"% variance)"), ylab=paste0("PC2 (",var.pct[],"% variance)"),main="PCA: normalized counts")

#THESE ARE ONLY USED WHEN STUDYING ESSENTIAL GENES WITH BAGEL
#Precision-Recall From BAGEL Bayes Factors ESSENTIAL GENES
# prTable = read.table()

#Fold-change plots of essential genes


#PLOTS SPECIFIC TO DRUGZ-CALCULATED DRUG SYNERGY
