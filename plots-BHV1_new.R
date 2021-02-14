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
setwd("~/_OSU_Project/R") #"C:/Users/gabve/OneDrive/Documents/OSU_Project/R"

T0 <- read.csv("Mock.csv", TRUE, ",")
T1 <- read.csv("30min.csv", TRUE, ",")
T2 <- read.csv("1.5h.csv", TRUE, ",")
T3 <- read.csv("3h.csv", TRUE, ",")

#All_FoldChange <- read.csv("Log2FoldChange.csv", TRUE, ",")
#All_FoldChange <- read.csv("Log2FoldChange_modif.csv", TRUE, ",")
All_FoldChange <- read.csv("Log2FoldChange_noRNA.csv", TRUE, ",")


#scarna_data <- filter (x, grepl('SCARNA', Gene))
id <- All_FoldChange$ï..Gene
All_FoldChange <- data.matrix(All_FoldChange[,2:5])
rownames(All_FoldChange) <- id
#plot(All_FoldChange)
heatmap(All_FoldChange[1:100,], Colv=NA, col=brewer.pal(11,"RdYlGn"))
library("RColorBrewer")
heatmap.2(All_FoldChange[50:80,],col=brewer.pal(11,"RdBu"),scale="row", trace="none", Colv = NA)

#To manually check values
#The first row in the heatmap
All_FoldChange["FAT2-2",]
#Mock   X30min   X90min  X180min 
#2.948320 2.213020 1.820950 0.898381 

summary(All_FoldChange["FAT2-2",])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4626   11380   15180   16190   18740   32100
test <- heatmap.2(All_FoldChange,scale="row")
#create function to calculate z-score
z_score <- function(x){
      (x-mean(x))/sd(x)
   }
z_score(All_FoldChange["FAT2-2",])
#T1a         T1b          T2          T3          N1          N2 
#-0.61945966 -1.23831240  1.70461190 -0.20346381  0.36826272 -0.01163874
#compare the manually calculating values with the ones used for the heatmap
#they are the same, except the ordering
barplot(z_score(All_FoldChange["FAT2-2",]))
test$carpet[,"FAT2-2"]
#N2         T1b         T1a          T2          T3          N1 
#-0.01163874 -1.23831240 -0.61945966  1.70461190 -0.20346381  0.36826272
barplot(All_FoldChange["FAT2-2",])
#boxplot(All_FoldChange[1:10,])


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
top_100_genes <- names(gene_score_ranked[1:100])

top_100_genes ## show the names of the top 10 genes
pca$rotation[top_100_genes,1] ## show the scores (and +/- sign)

write.csv(top_100_genes, "top_100_genes_noRNA.csv")
write.csv(names(gene_score_ranked[150:800]), "top_150-800_genes.csv")

#######


#Plot the volcano plot hightligting the proteins of intererst
#T0
plot(T0$Log2.Fold.Change, T0$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="Mock", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))
points(T0[5302,4], T0[5302,6], pch=16, col="blue", cex=1.5) #FAT2-2
points(T0[102,4], T0[102,6], pch=16, col="blue", cex=1.5) #ABL2
points(T0[1205,4], T0[1205,6], pch=16, col="blue", cex=1.5) #ATP7A
points(T0[27442,4], T0[27442,6], pch=16, col="blue", cex=1.5) #USP9Y
points(T0[21658,4], T0[21658,6], pch=16, col="blue", cex=1.5) #RAB39
points(T0[3061,4], T0[3061,6], pch=16, col="blue", cex=1.5) #CHRNB3
points(T0[4306,4], T0[4306,6], pch=16, col="blue", cex=1.5) #DOCK9-2
points(T0[4795,4], T0[4795,6], pch=16, col="blue", cex=1.5) #EPHA4
points(T0[6107,4], T0[6107,6], pch=16, col="blue", cex=1.5) #GNAQ
points(T0[8281,4], T0[8281,6], pch=16, col="blue", cex=1.5) #LNPEP
points(T0[19053,4], T0[19053,6], pch=16, col="blue", cex=1.5) #NBEAL1
points(T0[19762,4], T0[19762,6], pch=16, col="blue", cex=1.5) #OPRM1
points(T0[20314,4], T0[20314,6], pch=16, col="blue", cex=1.5) #PCGF5
points(T0[20821,4], T0[20821,6], pch=16, col="blue", cex=1.5) #PLEKHM3
points(T0[22640,4], T0[22640,6], pch=16, col="blue", cex=1.5) #SCN9A
points(T0[24161,4], T0[24161,6], pch=16, col="blue", cex=1.5) #SV2C
points(T0[25035,4], T0[25035,6], pch=16, col="blue", cex=1.5) #TMTC1
points(T0[27310,4], T0[27310,6], pch=16, col="blue", cex=1.5) #UHMK1

points(T0[27746,4], T0[27746,6], pch=16, col="red", cex=1.5) #WNT10A
points(T0[4146,4], T0[4146,6], pch=16, col="red", cex=1.5) #DLK1
points(T0[3682,4], T0[3682,6], pch=16, col="red", cex=1.5) #CTNNB1 b-cat



#T1
plot(T1$Log2.Fold.Change, T1$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="30 min", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))
points(T1[5302,4], T1[5302,6], pch=16, col="blue", cex=1.5) #FAT2-2
points(T1[102,4], T1[102,6], pch=16, col="blue", cex=1.5) #ABL2
points(T1[1205,4], T1[1205,6], pch=16, col="blue", cex=1.5) #ATP7A
points(T1[27442,4], T1[27442,6], pch=16, col="blue", cex=1.5) #USP9Y
points(T1[21658,4], T1[21658,6], pch=16, col="blue", cex=1.5) #RAB39
points(T1[3061,4], T1[3061,6], pch=16, col="blue", cex=1.5) #CHRNB3
points(T1[4306,4], T1[4306,6], pch=16, col="blue", cex=1.5) #DOCK9-2
points(T1[4795,4], T1[4795,6], pch=16, col="blue", cex=1.5) #EPHA4
points(T1[6107,4], T1[6107,6], pch=16, col="blue", cex=1.5) #GNAQ
points(T1[8281,4], T1[8281,6], pch=16, col="blue", cex=1.5) #LNPEP
points(T1[19053,4], T1[19053,6], pch=16, col="blue", cex=1.5) #NBEAL1
points(T1[19762,4], T1[19762,6], pch=16, col="blue", cex=1.5) #OPRM1
points(T1[20314,4], T1[20314,6], pch=16, col="blue", cex=1.5) #PCGF5
points(T1[20821,4], T1[20821,6], pch=16, col="blue", cex=1.5) #PLEKHM3
points(T1[22640,4], T1[22640,6], pch=16, col="blue", cex=1.5) #SCN9A
points(T1[24161,4], T1[24161,6], pch=16, col="blue", cex=1.5) #SV2C
points(T1[25035,4], T1[25035,6], pch=16, col="blue", cex=1.5) #TMTC1
points(T1[27310,4], T1[27310,6], pch=16, col="blue", cex=1.5) #UHMK1

points(T1[27746,4], T1[27746,6], pch=16, col="red", cex=1.5) #WNT10A
points(T1[4146,4], T1[4146,6], pch=16, col="red", cex=1.5) #DLK1
points(T1[3682,4], T1[3682,6], pch=16, col="red", cex=1.5) #CTNNB1 b-cat


#T2
plot(T2$Log2.Fold.Change, T2$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="1h 30 min", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))
points(T2[5302,4], T2[5302,6], pch=16, col="blue", cex=1.5) #FAT2-2
points(T2[102,4], T2[102,6], pch=16, col="blue", cex=1.5) #ABL2
points(T2[1205,4], T2[1205,6], pch=16, col="blue", cex=1.5) #ATP7A
points(T2[27442,4], T2[27442,6], pch=16, col="blue", cex=1.5) #USP9Y
points(T2[21658,4], T2[21658,6], pch=16, col="blue", cex=1.5) #RAB39
points(T2[3061,4], T2[3061,6], pch=16, col="blue", cex=1.5) #CHRNB3
points(T2[4306,4], T2[4306,6], pch=16, col="blue", cex=1.5) #DOCK9-2
points(T2[4795,4], T2[4795,6], pch=16, col="blue", cex=1.5) #EPHA4
points(T2[6107,4], T2[6107,6], pch=16, col="blue", cex=1.5) #GNAQ
points(T2[8281,4], T2[8281,6], pch=16, col="blue", cex=1.5) #LNPEP
points(T2[19053,4], T2[19053,6], pch=16, col="blue", cex=1.5) #NBEAL1
points(T2[19762,4], T2[19762,6], pch=16, col="blue", cex=1.5) #OPRM1
points(T2[20314,4], T2[20314,6], pch=16, col="blue", cex=1.5) #PCGF5
points(T2[20821,4], T2[20821,6], pch=16, col="blue", cex=1.5) #PLEKHM3
points(T2[22640,4], T2[22640,6], pch=16, col="blue", cex=1.5) #SCN9A
points(T2[24161,4], T2[24161,6], pch=16, col="blue", cex=1.5) #SV2C
points(T2[25035,4], T2[25035,6], pch=16, col="blue", cex=1.5) #TMTC1
points(T2[27310,4], T2[27310,6], pch=16, col="blue", cex=1.5) #UHMK1

points(T2[27746,4], T2[27746,6], pch=16, col="red", cex=1.5) #WNT10A
points(T2[4146,4], T2[4146,6], pch=16, col="red", cex=1.5) #DLK1
points(T2[3682,4], T2[3682,6], pch=16, col="red", cex=1.5) #CTNNB1 b-cat


#T3
plot(T3$Log2.Fold.Change, T3$Adjusted.P.value,
     pch=16, cex=0.7, col=8,
     main="3h", xlab ="Fold Change", ylab = "Adjusted p-value", cex.lab=1.5, 
     cex.axis=1.2, ylim=rev(range(c(1:0))))
points(T3[5302,4], T3[5302,6], pch=16, col="blue", cex=1.5) #FAT2-2
points(T3[102,4], T3[102,6], pch=16, col="blue", cex=1.5) #ABL2
points(T3[1205,4], T3[1205,6], pch=16, col="blue", cex=1.5) #ATP7A
points(T3[27442,4], T3[27442,6], pch=16, col="blue", cex=1.5) #USP9Y
points(T3[21658,4], T3[21658,6], pch=16, col="blue", cex=1.5) #RAB39
points(T3[3061,4], T3[3061,6], pch=16, col="blue", cex=1.5) #CHRNB3
points(T3[4306,4], T3[4306,6], pch=16, col="blue", cex=1.5) #DOCK9-2
points(T3[4795,4], T3[4795,6], pch=16, col="blue", cex=1.5) #EPHA4
points(T3[6107,4], T3[6107,6], pch=16, col="blue", cex=1.5) #GNAQ
points(T3[8281,4], T3[8281,6], pch=16, col="blue", cex=1.5) #LNPEP
points(T3[19053,4], T3[19053,6], pch=16, col="blue", cex=1.5) #NBEAL1
points(T3[19762,4], T3[19762,6], pch=16, col="blue", cex=1.5) #OPRM1
points(T3[20314,4], T3[20314,6], pch=16, col="blue", cex=1.5) #PCGF5
points(T3[20821,4], T3[20821,6], pch=16, col="blue", cex=1.5) #PLEKHM3
points(T3[22640,4], T3[22640,6], pch=16, col="blue", cex=1.5) #SCN9A
points(T3[24161,4], T3[24161,6], pch=16, col="blue", cex=1.5) #SV2C
points(T3[25035,4], T3[25035,6], pch=16, col="blue", cex=1.5) #TMTC1
points(T3[27310,4], T3[27310,6], pch=16, col="blue", cex=1.5) #UHMK1

points(T3[27746,4], T3[27746,6], pch=16, col="red", cex=1.5) #WNT10A
points(T3[4146,4], T3[4146,6], pch=16, col="red", cex=1.5) #DLK1
points(T3[3682,4], T3[3682,6], pch=16, col="red", cex=1.5) #CTNNB1 b-cat

#I selected only genes that have more than 4-fold compared Uninfected vs Latency
#Idividual barplots for: 

barplot(z_score(All_FoldChange["DLK1",]), main="DLK1")
barplot(z_score(All_FoldChange["WNT10A",]),main="WNT10A")
barplot(z_score(All_FoldChange["FAT2-2",]),main="FAT2-2")
barplot(z_score(All_FoldChange["ABL2",]), main="ABL2")
barplot(z_score(All_FoldChange["CTNNB1",]), main="beta-catenin")


MvsT1 <- read.csv("MvsT1.csv", TRUE, ",")

#scarna_data <- filter (x, grepl('SCARNA', Gene)) example how to filter
id <- MvsT1$ï..Gene
rownames(MvsT1) <- id
MvsT1 <- data.matrix(MvsT1)
MvsT1 <- MvsT1[,-1]
plot(MvsT1, pch=16, col=brewer.pal(11,"RdYlGn"))
heatmap(MvsT1, col=brewer.pal(11,"RdYlGn"))

###############################################VennDiagrams###############
# Prepare a palette of 3 colors with R colorbrewer:
#library(RColorBrewer)
myCol5 <- brewer.pal(5, "Pastel2")
myCol2 <- brewer.pal(2, "Pastel2")
myCol4 <- brewer.pal(4, "Pastel2")

#Venndiagram for MS
venn.diagram(
  x= list(All_FoldChange$Mock, All_FoldChange$X30min, All_FoldChange$X90min, All_FoldChange$X180min),
  category.names = c("Mock", "30min", "1.5h", "3h"),
  filename = "Shared_genes_BHV1infection_in_TG.png",
  output=TRUE,
  # Output features
  imagetype="png", 
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol4,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 10),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  #rotation = 1
)

######## Plot of all genes #########

plot(All_FoldChange$X180min, col=4, pch=16)
 points(All_FoldChange$X90min, col=3, pch=16)
 points(All_FoldChange$X30min, col=2, pch=16)
 legend(1,20, legend = c("3h","1.5h","30min"), col = 4:2, lty=1, cex = 1)
#points(All_FoldChange$Mock, col=5, pch=16)

#APO <- filter(All_FoldChange, grepl('CASP', ï..Gene))
#head(All_FoldChange$ï..Gene)

######################## ViSEAGO ##################
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("ViSEAGO")
 
 