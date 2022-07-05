#! /usr/bin/Rscript
library(DESeq2)
library(gplots)

file_groups<-commandArgs(TRUE)

case_control <- file(file_groups,"r")
first_line <- readLines(case_control,n=1)
close(case_control)
splitted_line = strsplit(first_line, '\t')
vec = as.vector(unlist(splitted_line))
groups = vec[2:length(vec)]

#MA*********************************************************************
counts <- read.delim(file_groups, sep="\t", header=T, row.names=1)
counts <- as.matrix(counts)
design <- data.frame( condition=factor( groups ) )
rownames(design) <- colnames(counts)
dataset <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)
dataset <- DESeq(dataset)
de_results <- results(dataset)
new_columns <- data.frame(GeneID=rownames(de_results))
de_results <- cbind( new_columns, de_results)

de_results <- de_results[ de_results$padj < 0.05 & complete.cases(de_results$padj), ]
de_results <- de_results[order(de_results$padj),]

write.table(de_results, file='deseq2.tsv', sep="\t", quote=F, row.names=F)

pdf("plots/differential/case_vs_control_MA_plot.pdf")
plotMA(dataset, main="MA plot", ylim=c(-2, 2))
dev.off()
#***********************************************************************


#HEATMAP****************************************************************
normalized_expression <- counts(dataset, normalized=T)
de_genes <- results(dataset)

genes <- normalized_expression[de_genes$padj < 0.05 & complete.cases(de_genes$padj),]
pdf("plots/differential/case_vs_control_heatmap.pdf")
heatmap.2( as.matrix(genes), labRow=F, col=redgreen(100), scale="row", dendrogram="column", trace="none", cexCol=1.2, hclust=function(x) hclust(x,method="centroid"), distfun=function(x) as.dist(1-cor(t(x))) )
dev.off()
#***********************************************************************


#VOLCANO****************************************************************
pdf("plots/differential/case_vs_control_VolcanoPlot.pdf")
plot(main = "case vs control", de_genes$log2FoldChange,-log10(de_genes$padj),pch=19,cex=0.5,xlab="Log2FoldChange",ylab="-log10(Adjusted P-value)",col=ifelse(de_genes$padj<0.05,"red","black"))
dev.off()
#***********************************************************************


#PCA********************************************************************

counts <- read.delim(file_groups, sep="\t", header=T, row.names=1)
counts <- as.matrix(counts)
design <- data.frame( condition=factor(groups ) )
rownames(design) <- colnames(counts)
dataset <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)
dds <- DESeq(dataset)
rld <- rlog(dds, blind=FALSE)
pdf("plots/differential/PCA.pdf")
plotPCA(rld, intgroup=c("condition"))
dev.off()

#***********************************************************************


#HISTOGRAM**************************************************************
counts <- read.delim(file_groups, sep="\t", header=T, row.names=1)
counts <- as.matrix(counts)
design <- data.frame( condition=factor( c( "c", "c", "c", "c", "t", "t", "t") ) )
rownames(design) <- colnames(counts)
dataset <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)
dds <- DESeq(dataset)
rld <- rlog(dds)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sep = " : "))
hc <- hclust(distsRL)

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("plots/differential/clustering.pdf")
heatmap.2(mat, Rowv = as.dendrogram(hc), symm = TRUE, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.off() 
#***********************************************************************
