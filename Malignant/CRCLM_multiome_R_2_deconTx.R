library('Seurat')
library('Signac')
library('scater')
library('SummarizedExperiment')
#library(loomR)
library('dplyr')
library(celda)
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 24000 * 1024^2)

colPal = DiscretePalette(24,palette="polychrome")

CRC <- readRDS('CRC_GEX_integrated.rds')
CRC
head(CRC@meta.data,2)

DimPlot(CRC, label=T)

FeaturePlot(CRC, features = 'CD45')

DefaultAssay(CRC) <- "RNA"
library(limma)
CRC.markers <- FindAllMarkers(CRC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(CRC.markers, file = 'CRC_markers_preDecon.csv')

CRC[["Clusters_all_cells_preDecon"]] <- Idents(object = CRC)
new.cluster.ids <- c("Epithelial", "Epithelial","Epithelial", "Epithelial","Myeloid", "Epithelial",
                     "T/NK/ILC", "Epithelial","Fibroblast", "Hepatocytes","Endothelial",
                     "Epithelial", "Epithelial","Epithelial","Epithelial","B",
                     "B")
names(new.cluster.ids) <- levels(CRC)
CRC <- RenameIdents(CRC, new.cluster.ids)
DimPlot(CRC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
CRC[["Cell_type_preDecon"]] <- Idents(object = CRC)

DefaultAssay(CRC) <- "RNA"
CRC.sce <- as.SingleCellExperiment(CRC)

###decontX
CRC.sce <- decontX(CRC.sce,
                   z = CRC.sce$Cell_type_preDecon,
                   batch = CRC.sce$Sample)

head(colData(CRC.sce)$decontX_contamination,2)
head(colData(CRC.sce),2)
tail(colnames(CRC.sce),2)
tail(Cells(CRC),2)

# plot percentage contamination from each cell:

CRC$decontX_contamination <- CRC.sce$decontX_contamination
CRC$decontX_clusters <- CRC.sce$decontX_clusters

pdf('Plots_Seurat/4_integrated_umap_deconX_contamination.pdf', height = 5, width = 5)
FeaturePlot(CRC, features = 'decontX_contamination')
dev.off()

markers <- list(Tcell_Markers = c("TRAC", "CD3D","CD3E"),
                Bcell_Markers = c("CD79A", "CD79B", "MS4A1"),
                Myeloid_Markers = c("SLC11A1", "CD68", "LYZ"),
                Epithelial_Markers = c("EPCAM","ELF3","HNF4A"),
                Endothelial_Markers = c('PECAM1','ARL15','VWF'),
                Fibroblast_Markers = c('DCN','ACTA2','COL3A1'),
                Hepatocyte_Markers = c('ALB','CP','DEFB1')
)
pdf('Plots_Seurat/5_markerPercentages.pdf', height = 9, width = 9)
plotDecontXMarkerPercentage(CRC.sce,
                            markers = markers,
                            z = 'Cell_type_preDecon',
                            assayName = "counts")
plotDecontXMarkerPercentage(CRC.sce,
                            markers = markers,
                            z = 'Cell_type_preDecon',
                            assayName = "decontXcounts")
dev.off()

saveRDS(CRC.sce, file = 'CRC_GEX_decon.sce.rds')
CRC.sce.decon.counts <- decontXcounts(CRC.sce)
saveRDS(CRC.sce.decon.counts, file = 'CRC_GEX_decon_counts.rds')
CRC.sce

CRC.decon.counts <- readRDS('CRC_GEX_decon_counts.rds')
CRC.sce <- readRDS('CRC_GEX_decon.sce.rds')

CRC.decon <- CreateSeuratObject(counts = CRC.decon.counts, meta.data = as.data.frame(colData(CRC.sce)))

### remove mito and ribo genes
CRC.decon
genesToKeep <- grep("^MT-|^RP[SL][ [:digit:]]|^RPLP[[:digit:]]|^RPSA", row.names(CRC.decon), value = T, invert=T)
length(genesToKeep)
CRC.decon <- CRC.decon[genesToKeep,]
CRC.decon

CRC.decon <- subset(CRC.decon, subset = decontX_contamination < 0.5)
CRC.decon

CRC.decon <- NormalizeData(CRC.decon)
CRC.decon <- FindVariableFeatures(CRC.decon, selection.method = "vst", nfeatures = 2000)
# vars.to.regress = c("Gex_nGenes")
CRC.decon <- ScaleData(CRC.decon,features = VariableFeatures(CRC.decon))
CRC.decon <- RunPCA(CRC.decon, features = VariableFeatures(object = CRC.decon))
CRC.decon <- FindNeighbors(CRC.decon, dims = 1:30)
CRC.decon <- FindClusters(CRC.decon, resolution = 0.5, algorithm = 3)
CRC.decon <- RunUMAP(CRC.decon, dims = 1:30)

### integrate, first split data by sample
CRC.decon.list <- SplitObject(CRC.decon, split.by = "Sample")
CRC.decon.list <- lapply(X = CRC.decon.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = CRC.decon.list)
CRC.decon.anchors <- FindIntegrationAnchors(object.list = CRC.decon.list, anchor.features = features, dims = 1:30)
CRC.decon.integrated <- IntegrateData(anchorset = CRC.decon.anchors, dims = 1:30, new.assay.name = 'integratedRNADecon')

### Perform integrated analysis on all cells
DefaultAssay(CRC.decon.integrated) <- "integratedRNADecon"
# Run the standard workflow for visualization and clustering
CRC.decon.integrated <- ScaleData(CRC.decon.integrated, verbose = FALSE,
                                  features = VariableFeatures(CRC.decon.integrated)
                                  )
CRC.decon.integrated <- RunPCA(CRC.decon.integrated, npcs = 30, verbose = FALSE)
CRC.decon.integrated <- RunUMAP(CRC.decon.integrated, reduction = "pca", dims = 1:30)
CRC.decon.integrated <- FindNeighbors(CRC.decon.integrated, reduction = "pca", dims = 1:30)
CRC.decon.integrated <- FindClusters(CRC.decon.integrated, resolution = 0.5, algorithm = 3)

saveRDS(CRC.decon.integrated, 'CRC_GEX_decon_integrated.rds')
CRC.decon.integrated <- readRDS('CRC_GEX_decon_integrated.rds')

### umap plots
pdf("Plots_Seurat/6_integrated_umap.pdf", height = 5, width = 5)
DimPlot(CRC.decon.integrated, reduction = "umap", label = T, cols=colPal)
DimPlot(CRC.decon.integrated, reduction = "umap", group.by = "Sample", shuffle = T, cols=colPal)
DimPlot(CRC.decon.integrated, reduction = "umap", group.by = "Cell_type_preDecon", cols=colPal)
dev.off()

pdf("Plots_Seurat/6_integrated_umap_QC.pdf", height = 10, width = 10)
FeaturePlot(CRC.decon.integrated, reduction = "umap",
            features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),
            max.cutoff = c(15000,6000,20,50)
)
FeaturePlot(CRC.decon.integrated, reduction = "umap",
            features = c("decontX_contamination")
)
dev.off()

### Use RNA assay for plotting gene expression
DefaultAssay(CRC.decon.integrated) <- "RNA"

pdf("Plots_Seurat/6_integrated_umap_HSP_CC.pdf", height = 8, width = 8)
FeaturePlot(CRC.decon.integrated, reduction = "umap", features = c("HSPA1A","HSPA1B",'FOS','JUN'))
FeaturePlot(CRC.decon.integrated, reduction = "umap", features = c("MKI67","TOP2A","PCNA"))
dev.off()

genes = c('EPCAM','ELF3','HNF4A',
          'SLC11A1','CD68','LYZ',
          'CD79A','BANK1','BLK',
          'CD3E','PTPRC','TRAC',
          'ALB','CP','DEFB1',
          'PECAM1','ARL15','VWF',
          'DCN','ACTA2','COL3A1',
          'PROX1')
pdf("Plots_Seurat/6_integrated_umap_markers.pdf", height = 24, width = 20)
FeaturePlot(CRC.decon.integrated, reduction = "umap", features = genes,
            order = T, min.cutoff = 1, cols = c('lightgrey','red3'))
dev.off()

pdf("Plots_Seurat/6-contamination_vln.pdf", height=5, width=9)
VlnPlot(CRC.decon.integrated, features="decontX_contamination", pt.size=0, group.by='Sample', y.max=0.75, cols=colPal)
dev.off()

### Save matrix
DefaultAssay(CRC.decon.integrated) <- "RNA"
writeMM(GetAssayData(object = CRC.decon.integrated, slot = "counts"), file='Matrices/GEX_decontaminated.mtx')
write.csv(rownames(CRC.decon.integrated), file='Matrices/GEX_decontaminated_genes.csv', quote = F, row.names=F)
write.csv(colnames(CRC.decon.integrated), file='Matrices/GEX_decontaminated_barcodes.csv', quote = F, row.names=F)
write.csv(CRC.decon.integrated@meta.data, file = 'Matrices/GEX_decontaminated_metadata.csv', quote=F, row.names=T)
write.csv(as.data.frame(Embeddings(CRC.decon.integrated[['umap']])), file = 'Matrices/GEX_decontaminated.CCA.csv', quote = F)
