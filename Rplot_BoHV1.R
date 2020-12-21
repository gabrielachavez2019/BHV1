#Install and load packages
#.libPaths(.libPaths()[2])

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cummeRbund")

library(cummeRbund)

getwd()
##Set as working directory in order to read files
setwd("~/_OSU_Project/BoHV-1/2ndRun-2020/Bt-BoHV1")

cuff_data <- readCufflinks('diff_out')

#warnings()
#Analysis visualization insitu
csDensity(genes(cuff_data))
csScatter(genes(cuff_data), 'Lat', 'Uni')
csVolcano(genes(cuff_data), 'Lat', 'Uni')

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

head(diffGenesOutput)

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

sp<-csScatter(genes(cuff_data),"Uni","Lat",smooth=T)
sp

v<-csVolcanoMatrix(genes(cuff_data))
v

replicates(cuff_data)

h<-csHeatmap(diffGenesData,cluster='both', replicates=T)
h
