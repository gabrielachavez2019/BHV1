#!/usr/bin/env Rscript

#############################################################
# Author: Gabriela Chavez Calvillo
# Copyright (c) gatoo, `r paste(format(Sys.Date(), "%Y"))`
# Email:  gabvet.c@gmail.com
# 
# Date: `r paste(Sys.Date())`
#
# Script Name: RNA-seq analysis post mapping with TOPHAT2 and CUFFLINKS  ##### TONSILS UNINFECTED VS LATENCY ####
#
# Script Description: This script was created based on "Differential gene and transcript 
# expression analysis of RNA-seq experiments with TopHat and Cufflinks" on Nature Communications
# the script can run after cufflinks analysis. 

#Install and load packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("cummeRbund")

#I use this script for Mouse transcriptomic analysis

library(cummeRbund)
setwd("/Users/joneslab/Documents/Analysis/Tonsils_UvsL")

cuff_data <- readCufflinks('diff_out_tonsil')
csDensity(genes(cuff_data))
###csVolcano(genes(cuff_data), 'Lat', 'Uni')

#Write a table with significant expressed Genes
gene_diff_data <- diffData(genes(cuff_data),'Lat', 'Uni') #Perform the group pairings that have biological significance
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sig_gene_data) #Count how many we have
#write.table(sig_gene_data, 'diff_genes_m_t0.txt', sep='\t',row.names = F, col.names = T, quote = F)

#Below command is equivalent to looking at the isoform_exp.diff file 

isoform_diff_data <-diffData(isoforms(cuff_data))
sig_isoform_data <- subset(isoform_diff_data, (significant == 'yes'))
nrow(sig_isoform_data) #Count how many we have

up_gene_data  <- subset(sig_gene_data, (log2_fold_change > 1))
#How many
nrow(up_gene_data)
down_gene_data <- subset(sig_gene_data, (log2_fold_change < -1))
nrow(down_gene_data)

geneids <- c(up_gene_data$gene_id)
myGenes <- getGenes(cuff_data, geneids)
pdf("heatmap.pdf")
csHeatmap(myGenes, cluster="both")
dev.off()

csBoxplot(genes(cuff_data))
csBoxplot(isoforms(cuff_data))


#Retrive significant gene IDs (XLOC) with a pre-specified alpha
#diffGeneIDs <- getSig(sig_gene_data,level="genes",alpha=0.05)
diffGeneIDs <- getSig(cuff_data,level='genes','Lat', 'Uni',alpha=0.05)

#Use returned identifiers to create a CuffGeneSet object with all relevant info for given genes
diffGenes<-getGenes(cuff_data,diffGeneIDs)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
#featureNames(diffGenes)
names<-featureNames(diffGenes)
row.names(names)=names$tracking_id
diffGenesNames<-as.matrix(names)
diffGenesNames<-diffGenesNames[,-1]
head(diffGenesNames)  ### Just to check is working
#To get a table with the gene names for GO analysis
write.table(diffGenesNames, 'diff_genes_names_t0-t1.txt', sep='\t',row.names = F, col.names = T, quote = F)

# get the data for the significant genes, yes AGAIN
diffGenesData<-diffData(diffGenes)
row.names(diffGenesData)=diffGenesData$gene_id
diffGenesData<-diffGenesData[,-1]

# merge the two matrices by row names
diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")

write.table(diffGenesOutput, 'diff_genes_Tonsils2.txt', sep='\t',row.names = F, col.names = T, quote = F)
write.table(names, 'names.txt', sep='\t',row.names = F, col.names = T, quote = F)

tail(diffGenesOutput)

########################

isoform_diff_data <-diffData(isoforms(cuff_data), 'Lat', 'Uni')
sig_isoform_data <- subset(isoform_diff_data, (significant == 'yes'))
nrow(sig_isoform_data)
tss_diff_data <-diffData(TSS(cuff_data), 'Lat', 'Uni')

#####################

disp<-dispersionPlot(genes(cuff_data))
disp 
genes.scv<-fpkmSCVPlot(genes(cuff_data))
isoforms.scv<-fpkmSCVPlot(isoforms(cuff_data))

densRep<-csDensity(genes(cuff_data),replicates=T)
densRep
b<-csBoxplot(genes(cuff_data))
b

brep<-csBoxplot(genes(cuff_data),replicates=T)
brep

s<-csScatterMatrix(genes(cuff_data))
s

#sp<-csScatter(genes(cuff_data),"Fem_dLAT","Fem_U1",smooth=T)
sp<-csScatter(genes(cuff_data),'Lat', 'Uni',smooth=T)
sp

v<-csVolcanoMatrix(genes(cuff_data))
v

replicates(cuff_data)

dend<-csDendro(genes(cuff_data))  #To build a dendogram, but does not work !!
dend

m<-MAplot(genes(cuff_data),'Lat', 'Uni')
m

mCount<-MAplot(genes(cuff_data),'Lat', 'Uni',useCount=T)
mCount

##To get detailed information for now I don't see the importance, I knew it
runInfo(cuff_data)
replicates(cuff_data)
gene.features<-annotation(genes(cuff_data))
head(gene.features)
gene.fpkm<-fpkm(genes(cuff_data))
head(gene.fpkm)
gene.matrix<-fpkmMatrix(genes(cuff_data))  ###To create a matrix with FPKM information
head(gene.matrix)
gene.count.matrix<-countMatrix(genes(cuff_data))  ## To normalize that data count
head(gene.count.matrix)

myGenes<-getGenes(cuff_data,(names$gene_short_name))
myGenes
head(fpkm(myGenes))

#h<-csHeatmap(myGenes,cluster='both')

h<-csHeatmap(myGenes,cluster='both',replicates=T)
h
b<-expressionBarplot(myGenes)
b

s<-csScatter(myGenes, 'Lat', 'Uni',smooth=T)
s

v<-csVolcano(myGenes, 'Lat', 'Uni')
v
############
# Analysis of specific gene candidates, I sued my awk script before but his one is more efficient
############
jones_list <- read.csv("diff_genes_names.csv", TRUE, ",")
myGeneIds <- (jones_list$few1)
myGenes2 <- getGenes(cuff_data, myGeneIds)

head(fpkm(myGenes2))

h<-csHeatmap(myGenes2,cluster='both')
#h<-csHeatmap(myGenes,cluster='both',replicates=T)
h
b<-expressionBarplot(myGenes2)
b

myGenes<-getGenes(cuff_data,(jones_list$Gene))
myGenes
head(fpkm(myGenes))
head(fpkm(isoforms(myGenes)))

gl.iso.rep<-expressionPlot(isoforms(myGenes2),replicates=F)
gl.iso.rep

gl.rep<-expressionPlot(myGenes2,replicates=TRUE)
gl.rep

#### Or individual gene
#myGeneId<-"LRP1"
myGene<-getGene(cuff_data,"RCAN1")
myGene
head(fpkm(isoforms(myGene)))

#h<-csHeatmap(myGenes,cluster='both')
h<-csHeatmap(myGene,cluster='both',replicates=T)
h
b<-expressionBarplot("myGenes", replicates=T)
b

myGene<-getGene(cuff_data,"XLOC_000046")
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb

#### Refine for some specific genes #####
stat3_path <- read.csv("../HSV1-mouse/jones_list2.csv", TRUE, ",")
myGenes<-getGenes(cuff_data,stat3_path$DiffPathway)
gl.iso.rep<-expressionPlot(isoforms(myGenes),replicates=T)
gl.iso.rep
igbs<-expressionBarplot(isoforms(myGenes),replicates=F)
igbs

head(fpkm(isoforms(myGenes)))


mySigMat<-sigMatrix(cuff_data,level='genes',alpha=0.05)
mySigMat

###Distance matrix
myDistHeat<-csDistHeat(genes(cuff_data), replicates=T)
myDistHeat

###Dimentional reduction
genes.PCA<-PCAplot(genes(cuff_data),"PC1","PC2")
genes.PCA

#genes.MDS<-MDSplot(genes(cuff_data))
#genes.MDS

genes.PCA.rep<-PCAplot(genes(cuff_data),"PC1","PC2",replicates=T)
genes.PCA.rep

genes.MDS.rep<-MDSplot(genes(cuff_data),replicates=T)
genes.MDS.rep

ic<-csCluster(myGenes,k=4)
head(ic$cluster)

icp <- csClusterPlot(ic)
icp

###Finding similar expression
mySimilar<-findSimilar(cuff_data,"RCAN1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
mySimilar.expression

write.table(most_t1, 'most_0-1.txt', sep='\t',row.names = F, col.names = T, quote = F)

myGene<-getGene(cuff_data,"XLOC_000743")
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb

