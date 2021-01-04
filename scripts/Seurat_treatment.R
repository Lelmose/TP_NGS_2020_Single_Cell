# Installs:
remotes::install_github('satijalab/seurat-wrappers')
install_github("velocyto-team/velocyto.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")


library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(patchwork)
library('umap-learn')

# Read the loom file, the resulting object contains 4 tabs: spliced, unspliced, ambiguous and spanning
tooth <- ReadVelocity('/ifb/data/mydatalocal/loom_files/onefilepercell_SRR11201752_and_others_DFMGU.loom')

# Create the Seurat object from tooth
incisor <- Seurat::as.Seurat(x = tooth)

# In a first time we only consider spliced RNA because we look for  gene expression
incisor[["RNA"]]<-incisor[["spliced"]]

# Rename the percentag of mitochondrial RNA
incisor[["percent.mt"]] <- Seurat::PercentageFeatureSet(incisor, pattern = "^mt-")

# Visualize the distribution of the three parameters
Seurat::VlnPlot(incisor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Threshold of the parameters on arbitrary basis in order to select cells with enough counts and not too much mitochondrial RNA

# First selection of the number of counts

ncount <- incisor[["nCount_RNA"]]
dim(ncount)
ncount <- sort(ncount[["nCount_RNA"]])
length(ncount)

# Select the lowest 5 percent values of nCount
low_born <- as.integer(5*length(ncount)/100)
high_born <- as.integer(95*length(ncount)/100)
ncount <- ncount[low_born:high_born]
# Print the lenght to ensure that selection has been done
length(ncount)

# First selection of the number of different genes

nfeature <- incisor[["nFeature_RNA"]]
dim(nfeature)
nfeature <- sort(nfeature[["nFeature_RNA"]])
length(nfeature)

# Select the lowest 5 percent values of nCount
low_born <- as.integer(5*length(nfeature)/100)
high_born <- as.integer(95*length(nfeature)/100)
nfeature <- nfeature[low_born:high_born]

length(nfeature)

# Purification 
incisor <- subset(incisor, subset = nCount_RNA %in% ncount & nFeature_RNA %in% nfeature & percent.mt < 15)
Seurat::VlnPlot(incisor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Display the number of cells that passed through the filter
dim(incisor)
# Plot the different parameters in function of each other to search for correlation
plot1 <- Seurat::FeatureScatter(incisor, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- Seurat::FeatureScatter(incisor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalization step
incisor<- Seurat::NormalizeData(incisor)

# Set the threshold on variability at 10 percents.
variable <- as.integer(0.1 * dim(incisor)[1])
variable
incisor <- Seurat::FindVariableFeatures(incisor, selection.method = "vst", nfeatures = variable)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(incisor), 10)

# plot variable features with and without labels

Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)



# Scaling the data
all.genes <- rownames(incisor)
incisor <- Seurat::ScaleData(incisor, features = all.genes)


# Linear dimensional reduction

# First: PCA
incisor <- Seurat::RunPCA(incisor, features = Seurat::VariableFeatures(object = incisor))

# Print the five main components of the five main axis
print(incisor[["pca"]], dims = 1:5, nfeatures = 5)

# Display the genes involved in the first two axis
Seurat::VizDimLoadings(incisor, dims = 1:2, reduction = "pca")

# Display a two dimensions plot of the data with the two first PCA components as axis
Seurat::DimPlot(incisor, reduction = "pca")                       



# Attempts to determine the optimal dimension of the PCA

# First method: print heatmap on each axis to have a glance at the variability 
# Here only the first and last axis are displayed as a comparison.
Seurat::DimHeatmap(incisor, dims = 1, cells = 500, balanced = TRUE)
Seurat::DimHeatmap(incisor, dims = 20, cells = 500, balanced = TRUE)


# Second method: Jackstraw 

# This first line create the null hypothese plot in dashed line
incisor <- Seurat::JackStraw(incisor, num.replicate = 100)
# Display the 20 axis
incisor <- Seurat::ScoreJackStraw(incisor, dims = 1:20)


# Third method: Elbowplot

Seurat::JackStrawPlot(incisor, dims = 1:20)
Seurat::ElbowPlot(incisor)


# First attempts at clustering

# As the clustering algorithm is similar to the Louvain algorithm for graphs,
# the first step is the compute the detection of Neighbors for each cell. The data
# are embedded in the 20 dimensions space previously defined.
incisor <- Seurat::FindNeighbors(incisor, dims = 1:20)
incisor <- Seurat::FindClusters(incisor, resolution = 0.5)

# Display the clusters on the PCA plot 
Seurat::DimPlot(incisor, reduction = "pca")

# Visualization t-SNE

incisor <- Seurat::RunTSNE(incisor, dims =1:20)
Seurat::DimPlot(incisor, reduction = "tsne")

# Visualization UMAP

incisor <- Seurat::RunUMAP(incisor, dims =1:20,n.neighbors = 10)
Seurat::DimPlot(incisor, reduction = "umap")

# Differentially expressed features

# Cluster 1
cluster1.markers <- FindMarkers(incisor, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- Seurat::FindMarkers(incisor, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
row.names(cluster5.markers)[1:100]

# find markers for every cluster compared to all remaining cells, report only the positive ones
incisor.markers <- Seurat::FindAllMarkers(incisor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
incisor.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(incisor, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

head(cluster1.markers, n = 5)
VlnPlot(incisor, features = c("Selenop", "Ctsc"))

Seurat::FeaturePlot(incisor, features = c("Spry4"))
# On nomme les cluster en fonction des types qu'on a inférés
new.cluster.ids <- c("Macrophages 1","Distal Pulp", "Distal Pulp", "Endothelium","Dental Epithelium","Macrophages 2",
                     "Preodontoblasts", "Dental follicle","Innate Lymphocytes","Leukocytes","Oseoblasts","Glia","Perivascular")
names(new.cluster.ids) <- levels(incisor)
incisor <- RenameIdents(incisor, new.cluster.ids)
DimPlot(incisor, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

                     
# Velocity estimation

emat <- tooth$spliced
nmat <- tooth$unspliced
emat <- emat[,incisor@assays$RNA@data@Dimnames[[2]]]
nmat <- nmat[,incisor@assays$RNA@data@Dimnames[[2]]]

cluster.label <- incisor@active.ident
PCA <- incisor@reductions$pca@cell.embeddings
cell.dist <-as.dist(1-armaCor(t(PCA)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=13,cell.dist=cell.dist,fit.quantile=fit.quantile)
emb <- incisor@reductions[["umap"]]@cell.embeddings
cell.colors <- ident.colors[Idents(object = incisor)]
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
cell.colors

# Nouvelle tentative
incisor <- RunVelocity(object = incisor, deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 3, reduction = "pca")
t <- Tool(object = incisor, slot = "RunVelocity")
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = incisor)))
names(x = ident.colors) <- levels(x = incisor)
cell.colors <- ident.colors[Idents(object = incisor)]
names(x = cell.colors) <- colnames(x = incisor)
show.velocity.on.embedding.cor(emb = Embeddings(object = incisor, reduction = "umap"), vel = Tool(object = incisor, 
                                                                                                  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


test_incisor<-incisor
dental_pulp = c("Macrophages 1","Apical Pulp", "Distal Pulp", "Endothelium","Dental Epithelium","Dental follicle")
test_incisor <- subset(test_incisor, subset = seurat_clusters %in% dental_pulp)
Seurat::DimPlot(test_incisor, reduction = "umap")
test_incisor <- RunVelocity(object = test_incisor, deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 3, reduction = "pca")
test_t <- Tool(object = test_incisor, slot = "RunVelocity")
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = test_incisor)))
names(x = ident.colors) <- levels(x = incisor)
cell.colors <- ident.colors[Idents(object = test_incisor)]
names(x = cell.colors) <- colnames(x = test_incisor)
show.velocity.on.embedding.cor(emb = Embeddings(object = test_incisor, reduction = "umap"), vel = Tool(object = test_incisor, 
                                                                                                  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

# Tentative pour isoler les gènes qui sont à l'origine de la vélocité

velo <- test_incisor@tools$RunVelocity$gamma
velo <- sort(velo, decreasing = TRUE)
head(velo, n=10)
Seurat::FeaturePlot(test_incisor, features = c("Slc25a21"))
length(velo)
