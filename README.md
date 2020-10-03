# Cluster-Analysis-of-Gene-Co-expression

#Preprocess and filter the GSE29617 microarray gene expression data from the Emory flu study 
#to obtain a gene x subject dataset for log2(day7/day0) gene expression changes.

#computing the correlation or distance matrix as needed for greedycommunity (modularity) and WGCNA to get 
# biological interpretations of the clusters, and determining this by using MsigDb
#(http://software.broadinstitute.org/gsea/msigdb/annotate.jsp)
#to investigate gene sets with Reactome pathways. 
#Other databases that could be useful are KEGG, GO Molecular Function, and Immunologic signatures. 
#Based on similar MsigDB annotations, we can match WGCNA clusters with modularity clusters.
  
#Data used: GSE29617_series_matrix.txt.gz
