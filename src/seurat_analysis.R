library(dplyr)
library(Seurat)
library(Matrix)


# Load the  dataset
data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
ob <- CreateSeuratObject(counts = data, project = "DO", min.cells = 3, min.features = 200)


#percent MT
ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))



ob <- NormalizeData(ob, normalization.method = "LogNormalize", scale.factor = 10000)

#identify highly variable features (high cell-to-cell variation)
ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ob), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ob)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(ob)
ob <- ScaleData(ob, features = all.genes)

ob <- RunPCA(ob,features = VariableFeatures(object = ob))

# Examine and visualize PCA results a few different ways
print(ob[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(ob, dims = 1:2, reduction = "pca")

DimPlot(ob, reduction = "pca")

#Determine the ‘dimensionality’ of the dataset

#jackstraw: randomly permute 1% of data, and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
#We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

ob <- JackStraw(ob, num.replicate = 100)
ob <- ScoreJackStraw(ob, dims = 1:20)

JackStrawPlot(ob, dims = 1:20)

ElbowPlot(ob)


##
ob <- FindNeighbors(ob, dims = 1:10)
ob <- FindClusters(ob, resolution = 0.5)
##

ob <- RunUMAP(ob, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ob, reduction = "umap", label=TRUE,split.by = "seurat_clusters")

cluster1.markers <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25)

ob.markers <- FindAllMarkers(ob, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ob.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


ob@meta.data$seurat_clusters
Idents(ob)[5]

#bulk rna counts
Matrix::rowSums(ob[["RNA"]]@counts)

#rna counts per cluster
Idents(ob)

c0Counts = ob[["RNA"]]@counts[,names(Idents(ob)[which(Idents(ob) == 0)])]
c0Sums = Matrix::rowSums(c0Counts)


cluster_RNA_sums = matrix(nrow = 17401, ncol=11)
colnames(cluster_RNA_sums) = paste0("c_",c(0:10))

for(i in c(0:10)){
  counts = ob[["RNA"]]@counts[,names(Idents(ob)[which(Idents(ob) == i)])]
  sums = as.matrix(Matrix::rowSums(counts))
  
  cluster_RNA_sums[,i+1] = sums
  rownames(cluster_RNA_sums) = rownames(counts)
}
