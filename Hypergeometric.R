#!/usr/bin/env Rscript

#############################################################
# Author: Gabriela Chavez Calvillo
# Copyright (c) gatoo, `r paste(format(Sys.Date(), "%Y"))`
# Email:  gabvet.c@gmail.com
#
# Date: `r paste(Sys.Date())`
#
# Script Name: RNA-seq analysis
#
# Script Description: This script was created based on "Differential gene and transcript
# expression analysis of RNA-seq experiments with TopHat and Cufflinks" on Nature Communications
# the script can run after cufflinks analysis.

#Install and load packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cummeRbund")
 library(cummeRbund)
cuff_data <- readCufflinks('diff_out')

#Analysis visualization insitu
csDensity(genes(cuff_data))
csScatter(genes(cuff_data), 'T0', 'T1')
csVolcano(genes(cuff_data), 'T0', 'T1')

#A quick overview to inspect the number of genes and transcripts that are DE between samples:
cuff_data

gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sig_gene_data)

#Write a table with significant expressed Genes
gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
write.table(sig_gene_data, 'diff_genes.txt', sep='\t',row.names = F, col.names = T, quote = F)

#To export the graphs
jpeg(file="Distrubution_of_expression_levels.jpeg")
csDensity(genes(cuff_data))
dev.off()

jpeg(file="Scatter_plot.jpeg")
csScatter(genes(cuff_data), 'T0', 'T1')
dev.off()

jpeg(file="Volcano_plot.jpeg")
csVolcano(genes(cuff_data), 'T0', 'T1')
dev.off()

mygene <- getGene(cuff_data,'KLF4')
jpeg(file="KLF4.jpeg")
expressionBarplot(mygene)
dev.off()

#Plot individual isoform expression level of selected genes with barplot
jpeg(file="KLF4_isoforms.jpeg")
expressionBarplot(isoforms(mygene))
dev.off()

#Venn Diagram needs this analysis is called "Hypergeometric testing" it can run as: dhyper, phyper, qhyper or rhyper
m <- 10; n <- 7; k <- 8
x <- 0:(k+1)
rbind(phyper(x, m, n, k), dhyper(x, m, n, k))
all(phyper(x, m, n, k) == cumsum(dhyper(x, m, n, k)))  # FALSE
## but error is very small:
signif(phyper(x, m, n, k) - cumsum(dhyper(x, m, n, k)), digits = 3)


####
cuff_data <- readCufflinks('diff_out')

#Write a table with significant expressed Genes
gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
#write.table(sig_gene_data, 'diff_genes.txt', sep='\t',row.names = F, col.names = T, quote = F)
#Retrive significant gene IDs (XLOC) with a pre-specified alpha
diffGeneIDs <- getSig(cuff_data,level="genes",alpha=0.05)

#Use returned identifiers to create a CuffGeneSet object with all relevant info for given genes
diffGenes<-getGenes(cuff_data,diffGeneIDs)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
names<-featureNames(diffGenes)
row.names(names)=names$tracking_id
diffGenesNames<-as.matrix(names)
diffGenesNames<-diffGenesNames[,-1]

# get the data for the significant genes
diffGenesData<-diffData(diffGenes)
row.names(diffGenesData)=diffGenesData$gene_id
diffGenesData<-diffGenesData[,-1]

# merge the two matrices by row names
diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")
write.table(diffGenesOutput, 'diff_genes.txt', sep='\t',row.names = F, col.names = T, quote = F)
