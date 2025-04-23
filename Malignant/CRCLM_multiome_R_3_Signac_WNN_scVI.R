library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(dplyr)
library(ggplot2)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 24000 * 1024^2)
set.seed(1234)

annotation <- readRDS('/data/BCI-CRC/SO/genomes/UCSC_hg38_annotation.rds')

# CRC <- readRDS('CRC_unfiltered.rds')
# CRC
# head(CRC@meta.data, 2)

CRC_RNA <- readRDS('CRC_GEX_decon_integrated.rds')
CRC_RNA
head(CRC_RNA@meta.data, 2)

### add metadata and embeddings from scanpy

X_umap <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/CRCLM_decon_scvi_Xumap.csv', header = F)
X_scVI <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/CRCLM_decon_scvi_XscVI.csv', header = F)
meta <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/CRCLM_decon_scvi_metadata.csv', header = T, row.names=1)
X_umap[1:5,1:2]
X_scVI[1:5,1:5]
meta[1,]

row.names(X_umap) <- row.names(meta)
row.names(X_scVI) <- row.names(meta)
colnames(X_umap) <- c('UMAP_1','UMAP_2')
colnames(X_scVI) <- paste0("scVI_", 1:10)
head(X_umap)
head(X_scVI)
dim(X_umap)
dim(X_scVI)
dim(meta)
dim(CRC_RNA@meta.data)

CRC_RNA@meta.data[1,]
length(Cells(CRC_RNA))
CRC_RNA <- CRC_RNA[,row.names(meta)]

CRC_RNA[["scVI"]] <- CreateDimReducObject(
  embeddings = as.matrix(X_scVI),
  key = 'scVI_',
  assay='integratedRNADecon')
CRC_RNA[["X_umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(X_umap),
  key = 'X_umap_',
  assay='integratedRNADecon')
CRC_RNA@meta.data <- meta

newNames <- gsub('_LM_','_LM#',colnames(CRC))
newNames <- gsub('_T_','_T#',newNames)
head(newNames,2)
CRC <- RenameCells(CRC, new.names = newNames)
head(colnames(CRC),2)

CRC <- CRC[,Cells(CRC_RNA)]
length(Cells(CRC))
tail(Cells(CRC),2)
tail(Cells(CRC_RNA),2)
CRC@meta.data <- CRC_RNA@meta.data

##### Process and integrate ATAC data
### Process ATAC
DefaultAssay(CRC) <- "ATAC"
# call peaks using MACS2
peaks <- CallPeaks(CRC, group.by = 'Cell_type')
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(CRC),
  features = peaks,
  cells = colnames(CRC)
)
# create a new assay using the MACS2 peak set and add it to the Seurat object
CRC[["Peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(CRC),
  annotation = annotation
)

saveRDS(CRC, 'CRC_ATAC_scVI.rds')
CRC <- readRDS('CRC_ATAC_scVI.rds')

# ATAC processing, perform LSI
DefaultAssay(CRC) <- "Peaks"
CRC <- FindTopFeatures(CRC, min.cutoff = 5)
CRC <- RunTFIDF(CRC)
CRC <- RunSVD(CRC)
CRC <- RunUMAP(CRC, reduction = "lsi", dims = 2:30, reduction.name = "umap_ATAC",)
# ATAC umap plots
pdf("Plots_Seurat_scVI/7_ATAC_umap.pdf", height = 5, width = 5)
DimPlot(CRC, reduction = "umap_ATAC", label = T)
DimPlot(CRC, reduction = "umap_ATAC", group.by = "Sample")
dev.off()

#####
# integrate ATAC data
CRC.list <- SplitObject(CRC, split.by = "Sample")
CRC.list
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = CRC.list,
  anchor.features = rownames(CRC.list),
  reduction = "rlsi",
  dims = 2:30
)
# integrate LSI embeddings
CRC <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = CRC[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
# create a new UMAP using the integrated embeddings
CRC <- RunUMAP(CRC, reduction = "integrated_lsi", dims = 2:30, reduction.name = "umap_ATAC_integrated")

# ATAC umap plots
pdf("Plots_Seurat_scVI/8_ATAC_integrated_umap.pdf", height = 12, width = 12)
DimPlot(CRC, reduction = "umap_ATAC_integrated", label = T)
DimPlot(CRC, reduction = "umap_ATAC_integrated", group.by = "Sample")
DefaultAssay(CRC) <- 'RNA'
CRC <- NormalizeData(CRC)
FeaturePlot(CRC, reduction = "umap_ATAC_integrated",
            features = c('HNF4A','COL3A1','VWF','CP','PTPRC','CD79A','SLC11A1','TRAC'))
dev.off()
saveRDS(CRC, 'CRC_ATAC_integrated_scVI.rds')

### embeddings are lost after data integration, add decontaminated RNA embeddings 
CRC_ATAC <- readRDS('CRC_ATAC_integrated_scVI.rds')
DefaultAssay(CRC_ATAC)
CRC_ATAC
CRC_RNA

CRC_ATAC[['integratedRNADecon']] <- CreateAssayObject(counts = CRC_RNA$RNA@counts)
CRC_ATAC[["scVI"]] <- CRC_RNA[["scVI"]]
CRC_ATAC[["X_umap"]] <- CRC_RNA[["X_umap"]]
dim(Embeddings(CRC_ATAC, reduction = "scVI"))
dim(Embeddings(CRC_ATAC, reduction = "X_umap"))
dim(Embeddings(CRC_ATAC, reduction = "integrated_lsi"))

CRC <- FindMultiModalNeighbors(CRC_ATAC, reduction.list = list("scVI", "integrated_lsi"),
                               dims.list = list(1:10, 2:30), modality.weight.name = c("RNA.weight","ATAC.weight"))
CRC <- RunUMAP(CRC, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CRC <- FindClusters(CRC, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)

DefaultAssay(CRC) <- 'integratedRNADecon'
CRC <- NormalizeData(CRC)

pdf('Plots_Seurat_scVI/9_WNN.pdf')
DimPlot(CRC, reduction = "X_umap", group.by = "Cell_type", label = T, repel = T,  cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("RNA X_umap")
DimPlot(CRC, reduction = "X_umap", group.by = "seurat_clusters", label = T, repel = T,  cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("RNA X_umap")
DimPlot(CRC, reduction = "X_umap", group.by = "Sample", label = F,shuffle=T,  cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("RNA X_umap")
FeaturePlot(CRC, reduction = "X_umap", features = "percent.mt") + ggtitle("percent.mt (RNA) X_umap")
FeaturePlot(CRC, reduction = "X_umap", features = "percent.ribo") + ggtitle("percent.ribo (RNA) X_umap")
FeaturePlot(CRC, reduction = "X_umap", features = "percent.ribo", max.cutoff = 20) + ggtitle("percent.ribo (RNA) X_umap")
FeaturePlot(CRC, reduction = "X_umap", features = "decontX_contamination") + ggtitle("decontX_contamination (RNA) X_umap")
FeaturePlot(CRC, reduction = "X_umap", features = "decontX_contamination", max.cutoff = 0.4) + ggtitle("decontX_contamination (RNA) X_umap")
DimPlot(CRC, reduction = "umap_ATAC_integrated", group.by = "Cell_type", label = T, repel = T, cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("ATAC")
DimPlot(CRC, reduction = "umap_ATAC_integrated", group.by = "seurat_clusters", label = T, repel = T, cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("ATAC")
DimPlot(CRC, reduction = "umap_ATAC_integrated", group.by = "Sample", label=F,shuffle=T,cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("ATAC")
FeaturePlot(CRC, reduction = "umap_ATAC_integrated", features = "nFrags") + ggtitle("nFrags (ATAC)")
FeaturePlot(CRC, reduction = "umap_ATAC_integrated", features = "nFrags", max.cutoff = 25000) + ggtitle("nFrags (ATAC)")
FeaturePlot(CRC, reduction = "umap_ATAC_integrated", features = "TSSEnrichment") + ggtitle("TSSEnrichment (ATAC)")
FeaturePlot(CRC, reduction = "umap_ATAC_integrated", features = "TSSEnrichment", max.cutoff = 10) + ggtitle("TSSEnrichment (ATAC)")
DimPlot(CRC, reduction = "wnn.umap", group.by = "Cell_type", label=T,repel=T, cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("WNN")
DimPlot(CRC, reduction = "wnn.umap", group.by = "seurat_clusters", label=T,repel=T, cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("WNN")
DimPlot(CRC, reduction = "wnn.umap", group.by = "Sample", label=F,shuffle=T, cols = DiscretePalette(24, palette = "polychrome")) + ggtitle("WNN")
FeaturePlot(CRC, reduction = "wnn.umap", features = "nFrags") + ggtitle("nFrags (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "nFrags", max.cutoff = 25000) + ggtitle("nFrags (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "TSSEnrichment") + ggtitle("TSSEnrichment (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "TSSEnrichment", max.cutoff = 10) + ggtitle("TSSEnrichment (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "percent.mt") + ggtitle("percent.mt (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "percent.ribo") + ggtitle("percent.ribo (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "percent.ribo", max.cutoff = 20) + ggtitle("percent.ribo (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "decontX_contamination") + ggtitle("decontX_contamination (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "decontX_contamination", max.cutoff = 0.4) + ggtitle("decontX_contamination (WNN)")
VlnPlot(CRC, features = "RNA.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend()
VlnPlot(CRC, features = "ATAC.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend()
plot(density(CRC@meta.data$decontX_contamination))
dev.off()

pdf('Plots_Seurat_scVI/9_WNN_GEX.pdf', height=12, width=12)
FeaturePlot(CRC, reduction = "wnn.umap",
            features = c('HNF4A','COL3A1','VWF','CP','PTPRC','CD79A','SLC11A1','TRAC'))
FeaturePlot(CRC, reduction = "umap_ATAC_integrated",
            features = c('HNF4A','COL3A1','VWF','CP','PTPRC','CD79A','SLC11A1','TRAC'))
dev.off()

saveRDS(CRC, file='CRC_wnn_scVI.rds')

DimPlot(CRC, reduction = 'X_umap', label=T, raster=T)

pdf('Plots_Seurat_scVI/9_WNN_deconTx.pdf', height=5, width=5)
FeaturePlot(CRC, reduction = "wnn.umap", features = "decontX_contamination") + ggtitle("decontX_contamination (WNN)")
FeaturePlot(CRC, reduction = "wnn.umap", features = "decontX_contamination", max.cutoff = 0.4) + ggtitle("decontX_contamination (WNN)")
FeaturePlot(subset(CRC, subset = decontX_contamination < 0.5), reduction = "wnn.umap", features = "decontX_contamination") + ggtitle("decontX_contamination < 0.5 (WNN)")
FeaturePlot(subset(CRC, subset = decontX_contamination < 0.5), reduction = "wnn.umap", features = "decontX_contamination", max.cutoff = 0.4) + ggtitle("decontX_contamination < 0.5 (WNN)")
FeaturePlot(subset(CRC, subset = decontX_contamination < 0.25), reduction = "wnn.umap", features = "decontX_contamination") + ggtitle("decontX_contamination < 0.25 (WNN)")
FeaturePlot(subset(CRC, subset = decontX_contamination < 0.25), reduction = "wnn.umap", features = "decontX_contamination", max.cutoff = 0.4) + ggtitle("decontX_contamination < 0.25 (WNN)")
dev.off()

### Save coordinates and metadata to load in scanpy
write.csv(CRC@meta.data, file = 'Matrices/WNN_decontaminated_metadata.csv', quote=F, row.names=T)
write.csv(as.data.frame(Embeddings(CRC[['wnn.umap']])), file = 'Matrices/WNN_decontaminated.WNN.csv', quote = F)
write.csv(as.data.frame(Embeddings(CRC[['integrated_lsi']])), file = 'Matrices/WNN_decontaminated.lsi.csv', quote = F)
write.csv(as.data.frame(Embeddings(CRC[['umap_ATAC_integrated']])), file = 'Matrices/WNN_decontaminated.ATAC.csv', quote = F)
