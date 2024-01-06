#load libraries


library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(MAST)
--------------------------------------------------
#convert h5 file to seuratobject
#input:dic_list(list of directory which contains count matrixes of each sample)
#output:so_list(list of seurat-objects for each sample)


dic_list<-list('list of directories')
m<-0
so_list<-list()
for (i in dic_list){
alldata<-list.files(i,full.names =T,pattern = '*.h5ad');
alldata<-lapply(alldata,zellkonverter::readH5AD);
z<-lapply(alldata,assay);
z<-lapply(z,CreateSeuratObject);
seurat-object<-merg(z);so_list[[m]]<-seurat-object;m<-m+1}

--------------------------------------------------
#calculate percentage of mitochondrial  gens and filter out low-quality  cells
#input:so_list(list of seurat-objects for each sample)
#output:filtered_list(filtered-seuratobject for each sample)


so_list <- lapply(X =so_list, FUN = function(x) {
   x[['percent.mt']]<-PercentageFeatureSet(object =x, pattern = "^MT-")})


filtered_list <- lapply(X =so_list, FUN = function(x) {
x<-subset(x=x, subset=nCount_RNA < quantile(nCount_RNA, .98) &  nFeature_RNA >quantile(nFeature_RNA, .02) & percent.mt < 40)})

----------------------------------------------------
#normalize data and find most variable gens
#input:filtered_list(filtered-seuratobject for each sample)
#output:norm_list(normalized-seuratobjects)

norm_list <- lapply(X =filtered_list, FUN = function(x) {
   x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
 })

------------------------------------------------------
#find anchors and integrate data to remove batch effect
#input:norm_list(normalized-seuratobjects)
#output:integrated_SO(integrated seurat object)


features <- SelectIntegrationFeatures(object.list = norm_list)
anchors <- FindIntegrationAnchors(object.list = norm_list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)
--------------------------------------------------------
#remove transposable elements temporarily from integrated data
input:integrated(integrated seurat object),ret(names of all transposable elements)
output:temp_integrated(integrated seurat object without retroelements)

DefaultAssay(integrated)<-'RNA'
temp_integrated<-integrated[- which(rownames(integrated) %in% ret$gene_id),]
DefaultAssay(temp_integrated)<-'integrated'
temp_integrated<-temp_integrated[- which(rownames(temp_integrated) %in% ret$gene_id),]
---------------------------------------------------------
#sacle data and perform Principal component analysis(pca) and find optimum number of principle component
#input:temp_integrated(integrated seurat object without retroelements)
#output:temp_integrated(seurat object without retroelements and with pca information),dim(optimum Dimensions of PCA reduction)


DefaultAssay(temp_integrated)<-'integrated'
temp_integrated<- ScaleData(object =temp_integrated)
temp_integrated <-RunPCA(object = temp_integrated)
pct <- temp_integrated[["pca"]]@stdev / sum(temp_integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
dim<-min(co1,co2)
----------------------------------------------------------
#find neighbors and cluster cells
#input:temp_integrated(seurat object without retroelements and with pca information),dim(optimum Dimensions of PCA reduction)
#output:temp_integrated(clustered seurat object)


temp_integrated<-FindNeighbors(temp_integrated, reduction = "pca",dims=1:dim)
temp_integrated<-FindClusters(temp_integrated,reduction='pca')
----------------------------------------------------------
#assign clustering to integrated data with retro elements
#input:temp_integrated(clustered seurat object)
#output:integrated(clustered data with retro elements)

DefaultAssay(temp_integrated)<-'RNA'
clusters<-temp_integrated@meta.data$seurat_clusters
integrated[['cluster']]<-clusters
-----------------------------------------------------------
#label clusters based on average expression and percent expressed(dot plot)
#input:integrated(clustered data with retro elements)
#output:dotplot 



DefaultAssay(integrated)<-'RNA'
Idents(integrated)<-'cluster'
integrated<-NormalizeData(integrated)
fe<-c("cell-markers")
DotPlot(integrated,features=fe)+RotatedAxis()
-------------------------------------------------------------
#Differential expression between all clusters
#input:integrated(clustered data with retro elements)
#output:all.markers(table of differentially expressed genes)


DefaultAssay(integrated)<-'RNA'
Idents(integrated)<-'integrated_snn_res.0.8'
integrated<-NormalizeData(integrated)
all.markers <- FindAllMarkers(object =CH.integrated,test.use='MAST'))


