#load libraries

library(Seurat)
library(psych)
----------------------------------------------------
#calculate average expression levels of intended sets of genes and upregulated genes in one neutrophil cluster
#input:list of pathways and list of thier GENS(gene-sets of intended pathways),neutrophil_cluster and its upregulated genes
#output:ex_table(table of average expression levels of intended pathways and cluster's upregulated-genes in each cell)


neutrophil_cluster<-AddModuleScore(neutrophil_cluster,features = 'list of gene-set + upregulated genes in each cluster', name = 'vectors of pathway names+upregulated-gene')
ex_table<-neutrophil_cluster@meta.data[,c('vectors of pathway-columns and upregulated-genes')]

----------------------------------------------------
#correlation between these pathways and upregulated genes
#input:ex_table(table of average expression levels of intended pathways and cluster's upregulated-genes in each cell)
#output:cor[1](correlation table) and matpadj(its p.adjvalue table)

ex_table<-ex_table[,-1]
cor<-corr.test(ex_table)
padj<-cor$p.adj
matpadj<-diag(x=0,nrow = 'number of ex_table columns',ncol = 'number of ex_table columns')
matpadj[upper.tri(matpadj, diag = F)] <- padj
matpadj[lower.tri(matpadj, diag = F)] <- padj
