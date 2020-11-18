# Create by Sumeet Gupta in August 2020 for analysis of subset of cells in the 10x library. Cells were subset based classification identified by cellranger.
# cellranger identified human and mouse cells. Human cells were used for downstream analysis.

library(Seurat)

data_dir <- "/lab/solexa_public/Weinberg/200718_WIGTC-HISEQ2B_CE57AANXX/cellranger/v4/L21_577_hg38/filtered_feature_bc_matrix/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
allcells.577 <- CreateSeuratObject(counts = data)
rm(data)
cellidentity <- read.csv(file.path(data_dir,"../../L21_577_Mixed/analysis/gem_classification.csv"), row.names = 1)
cellidentity$lib <- 1
allcells.577 <- AddMetaData(allcells.577, cellidentity)
allcells.577 <- subset(allcells.577, subset = call == "GRCh38")

data_dir <- "/lab/solexa_public/Weinberg/200718_WIGTC-HISEQ2B_CE57AANXX/cellranger/v4/L21_578_hg38/filtered_feature_bc_matrix/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
allcells.578 <- CreateSeuratObject(counts = data)
rm(data)
cellidentity <- read.csv(file.path(data_dir,"../../L21_578_Mixed/analysis/gem_classification.csv"), row.names = 1)
cellidentity$lib <- 2
allcells.578 <- AddMetaData(allcells.578, cellidentity)
allcells.578 <- subset(allcells.578, subset = call == "GRCh38")

allcells <- merge(x = allcells.577, y=allcells.578)

allcells[["percent.mt"]] <- PercentageFeatureSet(object = allcells, pattern = "^MT-*")

jpeg(paste("QC-AllCells.jpeg", sep=""), width=600, height=600)
VlnPlot(object = allcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(paste("QC-AllCells-2.jpeg", sep=""), width=600, height=600)
FeatureScatter(allcells, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

jpeg(paste("QC-AllCells-3.jpeg", sep=""), width=600, height=600)
FeatureScatter(allcells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

allcells <- subset(allcells, subset = nFeature_RNA > 500 & percent.mt < 25)

jpeg(paste("QC-AllCells-Filtered.jpeg", sep=""), width=600, height=600)
VlnPlot(object = allcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(paste("QC-AllCells-Filtered-2.jpeg", sep=""), width=600, height=600)
FeatureScatter(allcells, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

jpeg(paste("QC-AllCells-Filtered-3.jpeg", sep=""), width=600, height=600)
FeatureScatter(allcells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

allcells <- NormalizeData(allcells)

allcells <- FindVariableFeatures(allcells, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(allcells)
allcells <- ScaleData(allcells, features = all.genes)

allcells <- RunPCA(allcells, features = VariableFeatures(object = allcells))

ElbowPlot(object = allcells)

allcells <- FindNeighbors(allcells, dims = 1:10)
allcells <- FindClusters(allcells, resolution = 0.5)
allcells <- RunUMAP(object = allcells, dims = 1:10)
allcells <- RunTSNE(object = allcells, dims = 1:10)

current.cluster.ids = FetchData(object = allcells, vars = c("ident"))
new.cluster.ids = current.cluster.ids
new.cluster.ids$ident = as.numeric(as.matrix(current.cluster.ids$ident)) + 1
Idents(object = allcells) = new.cluster.ids
Idents(object = allcells) = factor(Idents(object = allcells), levels = sort(as.numeric(levels(Idents(object = allcells)))))

jpeg(paste("AllCells-PCA.jpeg", sep=""), width=600, height=600)
DimPlot(object = allcells, reduction = "pca")
dev.off()

jpeg(paste("AllCells-UMAP.jpeg", sep=""), width=600, height=600)
DimPlot(allcells, reduction = "umap")
dev.off()

jpeg(paste("AllCells-TSNE.jpeg", sep=""), width=600, height=600)
DimPlot(allcells, reduction = "tsne")
dev.off()

jpeg(paste("AllCells-UMAP-GroupByLibrary.jpeg", sep=""), width=600, height=600)
DimPlot(allcells, reduction = "umap", group.by = "lib")
dev.off()

DimHeatmap(object = allcells, dims = 1, balanced = TRUE)

allcells.markers <- FindAllMarkers(allcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allcells.markers.sortedByPval = allcells.markers[order(allcells.markers$p_val),]

suppressMessages(library(dplyr))
top10 <- allcells.markers.sortedByPval %>%  group_by(cluster)  %>% do(head(., n=10))

jpeg(paste("AllCells-HeatMap.jpeg", sep=""), width=1000, height=600)
DoHeatmap(object = allcells, features = top10$gene) + NoLegend()
dev.off()

genes_uniquely_DE <- allcells.markers.sortedByPval %>% dplyr::filter(avg_logFC >= 0.58) %>% group_by(gene) %>%  dplyr::summarize(n=n()) %>%  dplyr::filter(n==1)

genes_uniquely_DE.markers.sortedByPval <- allcells.markers.sortedByPval[allcells.markers.sortedByPval$gene %in% genes_uniquely_DE$gene & 
                                                                          allcells.markers.sortedByPval$avg_logFC >= 0.58,]

top_marker_each <- genes_uniquely_DE.markers.sortedByPval %>% dplyr::group_by(cluster) %>% do(head(., n=10))

write.table(top_marker_each, file = "top_markers_uniq.xls", quote = FALSE, sep = "\t")
write.table(top10, file = "top_10.markers.xls", quote = FALSE, sep = "\t")
write.table(allcells.markers.sortedByPval, file = "All.markers.xls", quote = FALSE, sep = "\t")
