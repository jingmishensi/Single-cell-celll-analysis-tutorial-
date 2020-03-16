library(ggplot2)

library(Seurat)

file = "C:\\Users\\xxj07\\Desktop\\R practice\\python scanpy\\Seurat pbmc\\hg19"

data = Seurat::Read10X(data.dir = file)

#Construct the Seurat object
data.Seurat = CreateSeuratObject(counts = data,project = "pbmc",min.cells = 3,min.features = 200)

data.Seurat

mito_genes = grep(pattern = "^MT",x = rownames(data.Seurat))

mito_genes

library(Matrix)

#Compute the total counts and mitochondrial genes counts
total_counts = Matrix::colSums(data.Seurat)
mito_genes_counts = Matrix::colSums(data.Seurat[mito_genes,])

mito.gene.percent = mito_genes_counts/total_counts

data.Seurat =  AddMetaData(data.Seurat,metadata = mito.gene.percent,col.name = "mito.gene.percent")

VlnPlot(data.Seurat, features = c("nFeature_RNA", "nCount_RNA", "mito.gene.percent"), ncol = 3)

FeatureScatter(object = data.Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

FeatureScatter(object = data.Seurat,feature1 = "nCount_RNA",feature2 = "mito.gene.percent")

data.Seurat = NormalizeData(data.Seurat,normalization.method =  "LogNormalize", scale.factor = 10000)

data.Seurat@assays$RNA

#Normalize data
data.Seurat[["RNA"]]@data 

#find the high variable genes 
data.Seurat = FindVariableFeatures(object = data.Seurat,selection.method = "vst",nfeatures = 2000)

VariableFeatures(data.Seurat)

VariableFeaturePlot(data.Seurat)

plot1 <- VariableFeaturePlot(data.Seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data 
data.Seurat = ScaleData(data.Seurat, vars.to.regress  = c("nFeature_RNA", "nCount_RNA", "mito.gene.percent"))

# perform the PCA inear dimensional reduction
data.Seurat <- RunPCA(data.Seurat, features = VariableFeatures(object = data.Seurat))

DimPlot(data.Seurat,reduction = "pca")

DimHeatmap(data.Seurat, dims = 1:15, cells = 500, balanced = TRUE)

data.Seurat = JackStraw(data.Seurat,num.replicate = 100)

data.Seurat = ScoreJackStraw(data.Seurat)

JackStrawPlot(data.Seurat)

data.Seurat = FindNeighbors(data.Seurat,dims = 1:10)

seq(0.5,1.2,0.1)

# When we set resolution at 0.5 , I got 8 communities,but the tutorial got 9 communities 
data.Seurat = FindClusters(data.Seurat,resolution = seq(0.5,1.2,0.1))

data.Seurat = FindClusters(data.Seurat,resolution = 0.6)

#umap dimension reduction 
data.Seurat = RunUMAP(data.Seurat,dims = 1:10)

DimPlot(data.Seurat,reduction = "umap")

data.Seurat = RunTSNE(data.Seurat,dims = 1:10)

RunTSNE(data.Seurat,reduction = "tsne")

DimPlot(data.Seurat,reduction = "tsne")

allmarkers = FindAllMarkers(data.Seurat,min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)

allmarkers %>%group_by(cluster)%>%top_n(n = 2, wt = avg_logFC)

 VlnPlot(data.Seurat, features = c("MS4A1", "CD79A"))

VlnPlot(data.Seurat, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(data.Seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",  "CD8A"))

top10 <-allmarkers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
 DoHeatmap(data.Seurat, features = top10$gene) + NoLegend()

top10

filter(top10,cluster==3)

filter(top10,cluster==4)

filter(top10,cluster==5)

new.cluster.ids  = c("Naive CD4 T", "CD14+ Mono","Memory CD4 T" ,"B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

names(new.cluster.ids) = levels(data.Seurat)
 

names(new.cluster.ids)

data.Seurat = RenameIdents(data.Seurat, new.cluster.ids)
 DimPlot(data.Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) 


