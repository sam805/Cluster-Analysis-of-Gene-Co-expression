#Cluster Analysis of Gene Co-expression Data.
#Preprocess and filter the GSE29617 microarray gene expression data from the Emory flu study 
#to obtain a gene x subject dataset for log2(day7/day0) gene expression changes.

#computing the correlation or distance matrix as needed for greedycommunity (modularity) and WGCNA to get 
# biological interpretations of the clusters, and determining this by using MsigDb
#(http://software.broadinstitute.org/gsea/msigdb/annotate.jsp)
#to investigate gene sets with Reactome pathways. 
#Other databases that could be useful are KEGG, GO Molecular Function, and Immunologic signatures. 
#Based on similar MsigDB annotations, we can match WGCNA clusters with modularity clusters.
  
#Data used: GSE29617_series_matrix.txt.gz

# load full series data
setwd("Users/sara/Desktop/PhD Research/Biological Network")
library(GEOquery)
library(Biobase)
# download GSE to current directory
# 2008 TIV data
gsm <- getGEO("GSE29617", destdir=".")

# or load from current directory
# gsm <- getGEO(filename="GSE29617_series_matrix.txt", 
#              destdir=".",
#              GSEMatrix=TRUE, AnnotGPL=FALSE)

# load expression data into expression set (eset) object
exprData <- exprs(gsm$GSE29617_series_matrix.txt.gz)
dim(exprData) # 54715 x 80

exprFeatures <- featureData(gsm$GSE29617_series_matrix.txt.gz)
dim(exprFeatures)  # 16 "features" of genes
colnames(exprFeatures)
exprFeatures$ID # affy probe ids
geneNames <- exprFeatures$"Gene Symbol"  # gene names, factor, may want as.character
sum(geneNames=="")/length(geneNames) # 23% do not have gene name
# replace affy probe id rownames with gene symbols
rownames(exprData) <- as.character(geneNames)
# let's remove data with no name probes
hasName <- !(geneNames=="")
exprData.noNameFilter <- exprData[hasName,]
geneNames.noMissing <- geneNames[hasName]
dim(exprData.noNameFilter) # 41796
# there are repeated genes because multiple probes per gene
# add a gene name column to aggregate under
tmp.data <- data.frame(gene=geneNames.noMissing,exprData.noNameFilter)
# compute median for each repeated gene
# takes a couple minutes
exprData.noNameFilter.medianCollapse <- aggregate(. ~ gene, data = tmp.data, median)
# we lose the rownames in the process, so we would like to move the gene column there
rownames(exprData.noNameFilter.medianCollapse) <- exprData.noNameFilter.medianCollapse$gene
dim(exprData.noNameFilter.medianCollapse) # 20741 x 81  need to get rid of gene column
# remove gene column
exprData.noNameFilter.medianCollapse$gene <- NULL
dim(exprData.noNameFilter.medianCollapse) # 20741 x 80  

# We have mapped genes to rownames and collapsed repeated genes to median
# Finally quantile normalize and log2 transform
library(preprocessCore)
exprData.normalize <- log2(normalize.quantiles(as.matrix(exprData.noNameFilter.medianCollapse)))
rownames(exprData.normalize)<-rownames(exprData.noNameFilter.medianCollapse)
colnames(exprData.normalize)<-colnames(exprData.noNameFilter.medianCollapse)
sum(exprData.normalize[,2])

# let's filter
# not using first two
genelowvalfilter <- function(dataMatrix) {
  # Remove gene profiles with low absolute values in dataMatrix. Returns:
  # 1) a logical vector mask identifying gene expression profiles in dataMatrix
  #    that have absolute expression levels in the lowest 10% of the data set.
  # 2) a data matrix containing filtered expression profiles.
  threshold <- quantile(dataMatrix, c(0.1))
  mask <- apply(dataMatrix, 1, function(x) all(x < threshold))
  fdata <- dataMatrix[!mask, ]
  
  # return the row mask and filtered data
  list(mask=!mask, fdata=fdata)
}

genevarfilter <- function(dataMatrix, percentile) {
  # calculates the variance for each gene expression profile in Data and returns Mask, 
  # which identifies the gene expression profiles with a variance less than the 10th 
  # percentile. Mask is a logical vector with one element for each row in Data. 
  # The elements of Mask corresponding to rows with a variance greater than the threshold 
  # have a value of 1, and those with a variance less than the threshold are 0.
  probeVariances <- apply(dataMatrix, 1, var)
  threshold <- quantile(probeVariances, c(percentile))
  mask <- apply(dataMatrix, 1, function(x) var(x) > threshold)
  fdata <- dataMatrix[mask, ]
  
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

cv.filter <- function(dataMatrix, threshold) {
  # coefficient of variation filter
  mask <- apply(dataMatrix, 1, function(x) {(sd(x)/abs(mean(x))) < threshold})
  fdata <- dataMatrix[mask, ]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

# shooting for about 5000 genes
exprData.filter <- cv.filter(exprData.normalize,.025)$fdata  # 5002 genes
dim(exprData.filter)
which(rownames(exprData.filter)=="GAPDH") # just curious
plot(exprData.filter[3298,])
#varFilter <- genevarfilter(exprData.normalize,.75)
#dim(varFilter$fdata)

# Within all of these columns of data, are day0, 3 and day7 measurements on the same individuals
# We want to separate these data out by day

# get meta data on the expression data, like day of blood draw and titer phenotype info
metaData <- phenoData(gsm$GSE29617_series_matrix.txt.gz) # biobase function
names(pData(metaData))  # fields
metaData$data_processing  # 
metaData$description  # example: "2008-TIV-78-D7," field 4 contains D0, D3, D7 info for 80 subjects

# let's make this description the column header
colnames(exprData.filter) <- metaData$description
# grab columns by day
getDayFromDesc <- function(valueString) {
  # example: 2008-TIV-2-D0
  splitValues <- unlist(strsplit(valueString, "-"))
  # return the 4th token
  splitValues[4]
}

#----------- Create data set of log FC between day7 and day0
# grab post-vaccine-day of each sample from description text
sampleDays <- unlist(lapply(as.character(metaData$description), getDayFromDesc))
sampleDays
# separate data by columns that have the same blood draw day
day0.expr <- exprData.filter[,sampleDays=="D0"]
day3.expr <- exprData.filter[,sampleDays=="D3"]
day7.expr <- exprData.filter[,sampleDays=="D7"]
colnames(day0.expr)
colnames(day3.expr)
colnames(day7.expr)
# not all same, need to look at overlapping ids at some point (3rd token)
dim(day0.expr); dim(day3.expr); dim(day7.expr)

# grab columns by day
getIDFromDesc <- function(valueString) {
  # example: 2008-TIV-2-D0
  splitValues <- unlist(strsplit(valueString, "-"))
  # return the 4th token
  splitValues[3]
}

day0subjectIds<-unlist(lapply(as.character(colnames(day0.expr)), getIDFromDesc))
day0subjectIds
day7subjectIds<-unlist(lapply(as.character(colnames(day7.expr)), getIDFromDesc))
day7subjectIds
length(intersect(day0subjectIds, day7subjectIds))
intersect(day0subjectIds, day7subjectIds)
day0subjectIds %in% day7subjectIds

colnames(day0.expr) <- day0subjectIds
colnames(day7.expr) <- day7subjectIds
# day0 is larger than day7
logFCd7d0.expr <- day7.expr - day0.expr[,day7subjectIds]
length(colnames(logFCd7d0.expr))
#------------- sort the average log FC
lfc.mean.vec <- rowMeans(logFCd7d0.expr)
names(lfc.mean.vec) <- rownames(logFCd7d0.expr)
lfc.mean.vec[1]

# should we take absolute value??
index.ordered.by.lfc <- order(abs(as.numeric(as.character(lfc.mean.vec))), decreasing=TRUE)
# all genes sorted by mean log fold change (lfc)
lfc.mean.sorted <- lfc.mean.vec[index.ordered.by.lfc]
lfc.mean.sorted[1:10]  # top 10
# gene names listed for pasting into msigdb
write.table(names(lfc.mean.sorted[1:10]),row.names=FALSE,quote=FALSE)

# extract an lfc data set from genes that show top 10 mean lfc:
top200.lfc.data <- logFCd7d0.expr[index.ordered.by.lfc[1:200],] # "top 10" genes
dim(top200.lfc.data)

# compute the cor of log fold change day7/day0
top200.lfc.data <- t(top200.lfc.data)
mycorr <- cor(top200.lfc.data)
# WGCNA
library(WGCNA)
# selecting the best threshold
thresh <- .2
flu.absThresh <- (abs(mycorr) > thresh) + 0

dissTOM <- GTOMdist(flu.absThresh, degree = 2)

FluCorrTree <- hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
plot(FluCorrTree, xlab="", sub="", main = "hclust Flu correlation on GTOM-based dissimilarity (threshold .2)")
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = FluCorrTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 20)

init_mod_sizes <- table(dynamicMods)
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
# Plot the dendrogram and colors underneath
par(mar = c(2,2,2,2))
plotDendroAndColors(FluCorrTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHand = 0.05,
                    main = "Flu dynamic cut module colors (Degree 2, threshold .2)")


clus1_GeneNames <- as.factor(colnames(top200.lfc.data)[as.numeric(which(dynamicMods==1))])
clus2_GeneNames <- as.factor(colnames(top200.lfc.data)[as.numeric(which(dynamicMods==2))])
clus3_GeneNames <- as.factor(colnames(top200.lfc.data)[as.numeric(which(dynamicMods==3))])
clus4_GeneNames <- as.factor(colnames(top200.lfc.data)[as.numeric(which(dynamicMods==4))])

# --- partial correlation
n <- nrow(mycorr)
pcorr.pvals <- matrix(0,n,n)
for (i in seq(1, n-1)){
  for (j in seq(i+1, n)){
    rowi <- mycorr[i,-c(i,j)]
    rowj <- mycorr[j,-c(i,j)]
    # vector of partial correlations for genes other than i and j
    tmp <- (mycorr[i,j]-rowi*rowj)/sqrt((1-rowi^2)*(1-rowj^2))
    # vector of fisher z transforms of the partial corrs
    # normally distributed with mean 0 and variance 1/(n-4)
    tmp.zvals <- .5 * log((1+tmp)/(1-tmp))
    tmp.s.zvals <- sqrt(n-4)*tmp.zvals
    tmp.pvals <- 2*pnorm(abs(tmp.s.zvals), 0, 1, lower.tail = FALSE)
    #  r
    pcorr.pvals[i,j] <- max(tmp.pvals)
    
  }
}

library(igraph)
pcorr.pvals.vec <- pcorr.pvals[upper.tri(pcorr.pvals)]
pcorr.pvals.adjust <- p.adjust(pcorr.pvals.vec, "BH")
pcorr.edges <- pcorr.pvals.adjust < .05
# length(pcorr.pvals.adjust[pcorr.edges]) # another way to show the edges
sum(pcorr.edges)
pcorr.A <- matrix(0,n,n)
pcorr.A[upper.tri(pcorr.A)] <- as.numeric(pcorr.edges)
g.pcorr <- graph.adjacency(pcorr.A, "undirected")
plot(g.pcorr, vertex.size=3, vertex.label=NA)

# modularity
thresh <- .2
Flu.absThresh <- (abs(mycorr) > thresh) + 0
g.fluDB <- graph.adjacency(Flu.absThresh, "undirected")

deg.vec <- degree(g.fluDB)
low.deg.vertices <- which(deg.vec <=2)
g.filtered <-delete.vertices(g.fluDB, as.numeric(low.deg.vertices) )

# subgraph centrality for Emory flu data	
A <- get.adjacency(g.filtered)	
Mat_Decompos <- eigen(A)	
N <- vcount(g.filtered)	
sub_cent <- numeric(N)	
for (i in 1:N){
  for (j in 1:N){	
    eigen_value <- Mat_Decompos$values[j]	
    eigen_vec <- Mat_Decompos$vectors[i,j]	
    sub_cent[i] <- sub_cent[i]+((eigen_vec^2)*(exp(eigen_value)))	
  }	
}
sorted_sub_cent <- sort(sub_cent, decreasing=T)	
top_200_cent_index <- match(sorted_sub_cent[1:200], sub_cent)	
top_200_gene_sets <- colnames(top1000.lfc.data[,top_200_cent_index])	
# Sub_Cent_scores <- subgraph_centrality(g.filtered, diag = FALSE)	
h <- hist(sorted_sub_cent[1:200], breaks = 200, main= "Histogram of SC_scores for 200 genes")	
h$counts	
h$breaks
  
##################
Flu.modularity<-fastgreedy.community(g.filtered,merges=TRUE, modularity=TRUE, membership=TRUE)
Flu.modularity$membership
Flu.modularity$merges
membership.ids <- unique(Flu.modularity$membership)
membership.ids
cat(paste('Number of detected communities =',length(membership.ids)))
cat("community sizes: ")
sapply(membership.ids,function(x) {sum(x==Flu.modularity$membership)})
cat("modularity: ")
max(Flu.modularity$modularity)
#karate.modularity$modularity

V(g.filtered)$color[Flu.modularity$membership==1] <- "green"
V(g.filtered)$color[Flu.modularity$membership==2] <- "red"
V(g.filtered)$color[Flu.modularity$membership==3] <- "blue"

plot(g.filtered,vertex.size=10,vertex.label=V(g.fluDB)$label,vertex.color=V(g.fluDB)$color)

write.table(V(g.filtered)$name[Flu.modularity$membership==1], quote = F, row.names = F)
write.table(V(g.filtered)$name[Flu.modularity$membership==2], quote = F, row.names = F)
write.table(V(g.filtered)$name[Flu.modularity$membership==3], quote = F, row.names = F)



fc.fn <- function(exprdata, i){
      genei.expr <- all.nol.filt[i,]
      gene1.fit <- glm(pheno.factor.relevel~gene1.expr,family=binomial)
      gene1.tdy <- tidy(gene1.fit)
      coefvec <- gene1.tdy$estimate # intercept, gene
      pvec <- gene1.tdy$p.value     # intercept, gene
      cbind(rownames(all.nol.filt)[i], coefvec[2], pvec[2])
}

# sapply the function to all genes
num.genes<-nrow(all.nol.filt)
lr.results.df<-data.frame(t(sapply(1:num.genes, lr.fn)))
colnames(lr.results.df) <- c("gene", "slope_coef", "slope_p")
# sort results by slope coefficient p-value
lr.results.sorted <- lr.results.df[order(as.numeric(as.character(lr.results.df$slope_p))), ]
lr.results.sorted[1:10,]


# a function for stripping the label from the GEO data "characteristics" entries
# grab the value, strip the label
stripLabel <- function(valueString) {
  # example 1 metaData$characteristics_ch1.1 "subject id: 2" returns 2
  # example 2 metaData$characteristics_ch1.4 "hai titer (day 0) - a/uruguay (h3n2): 5" returns 5
  splitValues <- unlist(strsplit(valueString, ": "))
  # return the second part of the string
  splitValues[2]
}

# all sample ids
as.character(metaData$geo_accession)

cbind(subjectIds,sampleDays)
table(subjectIds) # subjects 2(D0/D7), 4(D0/D3), 51(D0/D3), 70(D0/D7) have incomplete days
length(unique(subjectIds))

# 3 titers measured for each sample at day0 and day28
# phenotypes: uraguay h3n2, brisbane h1n1, brisbane
# unpack all the values from strings
# Some of this phenotype data is redundant because the samples repeat subjects for pre and post-vac expression 
# So, the same titer characteristic might be repeated for 2 or 3 times for a subject
# hai titer (day 0) - a/uruguay (h3n2):
a00 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.4), stripLabel)))
# hai titer (day 28) - a/uruguay (h3n2):
a28 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.5), stripLabel)))
# hai titer (day 0) - a/brisbane/59/2007 (h1n1):
b00 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.6), stripLabel)))
# hai titer (day 28) - a/brisbane/59/2007 (h1n1):
b28 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.7), stripLabel)))
# hai titer (day 0) - b/brisbane/3/2007:
c00 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.8), stripLabel)))
# hai titer (day 28) - b/brisbane/3/2007:
c28 <- as.numeric(unlist(lapply(as.character(metaData$characteristics_ch1.9), stripLabel)))

# compute fold changes
afc <- a28 / a00
bfc <- b28 / b00
cfc <- c28 / c00
# find the max fold change
bestPheno <- pmax(afc, bfc, cfc)
# threshold the class column based on max fold change
classColumn <- ifelse(bestPheno < 4, 0, 1)
