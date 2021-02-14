#BosTaurus-BHV1
#install.packages("VennDiagram","ggplot2")
library(ggplot2)
library(VennDiagram)
#install.packages("yarr")
library(yarrr)
library(RColorBrewer)
library(gplots)
myCol <- brewer.pal(3, "Pastel2")

#getwd()
##Set as working directory in order to read files
setwd("~/OSU_Project/R") #"C:/Users/gabve/OneDrive/Documents/OSU_Project/R"

T0 <- read.csv("Mock.csv", TRUE, ",")
T1 <- read.csv("30min.csv", TRUE, ",")
T2 <- read.csv("1.5h.csv", TRUE, ",")
T3 <- read.csv("3h.csv", TRUE, ",")

#All_FoldChange <- read.csv("Log2FoldChange.csv", TRUE, ",")
All_FoldChange <- read.csv("Log2FoldChange_modif.csv", TRUE, ",")

#scarna_data <- filter (x, grepl('SCARNA', Gene))
id <- All_FoldChange$ï..Gene
All_FoldChange <- data.matrix(All_FoldChange[,2:5])
rownames(All_FoldChange) <- id
#plot(All_FoldChange)
heatmap(All_FoldChange[1:800,], Colv=NA, col=brewer.pal(11,"RdYlGn"))
library("RColorBrewer")
heatmap.2(All_FoldChange[1:800,],col=brewer.pal(11,"RdBu"),scale="row", trace="none", Colv = NA)

#To manually check values
#The first row in the heatmap
data_matrix["Gene_08743",]
#T1a   T1b    T2    T3    N1    N2 
#10404  4626 32103 14288 19626 16079 
#summary(data_matrix["Gene_08743",])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4626   11380   15180   16190   18740   32100
test <- heatmap.2(data_matrix,scale="row")
#create function to calculate z-score
z_score <- function(x){
      (x-mean(x))/sd(x)
   }
z_score(data_matrix["Gene_08743",])
#T1a         T1b          T2          T3          N1          N2 
#-0.61945966 -1.23831240  1.70461190 -0.20346381  0.36826272 -0.01163874
#compare the manually calculating values with the ones used for the heatmap
#they are the same, except the ordering
test$carpet[,"Gene_08743"]
#N2         T1b         T1a          T2          T3          N1 
#-0.01163874 -1.23831240 -0.61945966  1.70461190 -0.20346381  0.36826272




#if (T0$Log2.Fold.Change/T1$Log2.Fold.Change > 4) 

test <- heatmap.2(All_FoldChange)
All_FoldChange[rev(test$rowInd), test$colInd]
## Row clustering (adjust here distance/linkage methods to what you need!)
hr <- hclust(as.dist(1-cor(t(All_FoldChange), method="pearson")), method="complete")

## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(All_FoldChange, method="spearman")), method="complete")

## Plot heatmap
heatmap.2(All_FoldChange, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row", density.info="none", trace="none")

## Return matrix with row/column sorting as in heatmap
All_FoldChange[rev(hr$labels[hr$order]), hc$labels[hc$order]]



pca <- prcomp(t(All_FoldChange), scale=F) 
## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])
## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="BHV1", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])
plot(gene_score_ranked)
plot(gene_score_ranked[1:800])
top_1000_genes <- names(gene_score_ranked[1:800])

top_10_genes ## show the names of the top 10 genes
pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

write.csv(top_1000_genes, "top_1000_genes.csv")
write.csv(names(gene_score_ranked[150:800]), "top_150-800_genes.csv")

#######


#Plot the volcano plot hightligting the proteins of intererst
#T0
plot(T0$Log2.Fold.Change, T0$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="Mock", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))
#points(results[45,1:2], pch=16, col="blue", cex=1.5) #RPUSD2

#T1
plot(T1$Log2.Fold.Change, T1$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="30 min", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))

#T2
plot(T2$Log2.Fold.Change, T2$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="1h 30 min", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))

#T3
plot(T3$Log2.Fold.Change, T3$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="3h", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))

#I selected only genes that have more than 4-fold compared Uninfected vs Latency

MvsT1 <- read.csv("MvsT1.csv", TRUE, ",")

#scarna_data <- filter (x, grepl('SCARNA', Gene))
id <- MvsT1$ï..Gene
rownames(MvsT1) <- id
MvsT1 <- data.matrix(MvsT1)
MvsT1 <- MvsT1[,-1]
plot(MvsT1, pch=16, col=brewer.pal(11,"RdYlGn"))
heatmap(MvsT1, col=brewer.pal(11,"RdYlGn"))

