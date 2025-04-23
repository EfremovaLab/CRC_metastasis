library('ArchR')
library('Seurat')
set.seed(1)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

input_CRC01_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC01_LM/outs/')[1]
names(input_CRC01_LM) <- "CRC01_LM"

input_CRC02_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC02_LM/outs/')[1]
names(input_CRC02_LM) <- "CRC02_LM"

input_CRC03_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC03_LM/outs/')[1]
names(input_CRC03_LM) <- "CRC03_LM"

input_CRC04_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC04_LM/outs/')[1]
names(input_CRC04_LM) <- "CRC04_LM"

input_CRC05_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC05_LM/outs/')[1]
names(input_CRC05_LM) <- "CRC05_LM"

input_CRC06_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC06_LM/outs/')[1]
names(input_CRC06_LM) <- "CRC06_LM"

input_CRC07_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC07_LM/outs/')[1]
names(input_CRC07_LM) <- "CRC07_LM"

input_CRC08_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC08_LM/outs/')[1]
names(input_CRC08_LM) <- "CRC08_LM"

input_CRC09_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC09_LM/outs/')[1]
names(input_CRC09_LM) <- "CRC09_LM"

input_CRC10_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC10_LM/outs/')[1]
names(input_CRC10_LM) <- "CRC10_LM"

input_CRC11_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC11_LM/outs/')[1]
names(input_CRC11_LM) <- "CRC11_LM"

input_CRC12_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC12_LM/outs/')[1]
names(input_CRC12_LM) <- "CRC12_LM"

input_CRC13_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC13_LM/outs/')[1]
names(input_CRC13_LM) <- "CRC13_LM"

input_CRC14_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC14_LM/outs/')[1]
names(input_CRC14_LM) <- "CRC14_LM"

input_CRC15_LM <- getInputFiles('/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC15_LM/outs/')[1]
names(input_CRC15_LM) <- "CRC15_LM"

ArrowFiles <- createArrowFiles(
 inputFiles = c(input_CRC01_LM, input_CRC02_LM, input_CRC03_LM, input_CRC04_LM, input_CRC05_LM,
                input_CRC06_LM, input_CRC07_LM, input_CRC08_LM, input_CRC09_LM, input_CRC10_LM,
                input_CRC11_LM, input_CRC12_LM, input_CRC13_LM, input_CRC14_LM, input_CRC15_LM),
 sampleNames = c(names(input_CRC01_LM), names(input_CRC02_LM), names(input_CRC03_LM), names(input_CRC04_LM), names(input_CRC05_LM),
                 names(input_CRC06_LM), names(input_CRC07_LM), names(input_CRC08_LM), names(input_CRC09_LM), names(input_CRC10_LM),
                 names(input_CRC11_LM), names(input_CRC12_LM), names(input_CRC13_LM), names(input_CRC14_LM), names(input_CRC15_LM)),
 minTSS = 4, #Dont set this too high because you can always increase later
 minFrags = 1000,
 addTileMat = TRUE,
 addGeneScoreMat = TRUE)

ArrowFiles

#ArchRProject
proj <- ArchRProject(ArrowFiles)

# DO NOT IMPORT GEX INTO ARCHR UNTIL DECONTAMINATED THE COUNTS
### Filter Cells based on ATAC metrics only!!!
proj

pdf("Plots_ArchR/1-1-QC_prefiltering.pdf")
plotGroups(ArchRProj = proj,colorBy = "cellColData",groupBy = "Sample",plotAs="violin",
           alpha=0.4, addBoxPlot=T, name = "TSSEnrichment", baseSize = 12)
plotGroups(ArchRProj = proj,colorBy = "cellColData",groupBy = "Sample",plotAs="violin",
           alpha=0.4, addBoxPlot=T, name = "nFrags", baseSize = 12)
dev.off()

head(proj@cellColData,2)
table(proj@cellColData$Sample)
proj <- proj[proj$TSSEnrichment > 4 & proj$nFrags > 1500]
table(proj@cellColData$Sample)
proj

### Add doublet scores based on ATAC
proj <- addDoubletScores(proj, force = T)

pdf("Plots_ArchR/1-2-QC.pdf")
plotGroups(ArchRProj = proj,colorBy = "cellColData",groupBy = "Sample",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE,
                 name = "TSSEnrichment", baseSize = 12)
plotGroups(ArchRProj = proj,colorBy = "cellColData",groupBy = "Sample",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE,
                 name = "nFrags", baseSize = 12)
dev.off()

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject_allCells/", load = FALSE)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)
library(future)
plan("multicore", workers = 8)
plan()
options(future.globals.maxSize = 24000 * 1024^2)
colPal = DiscretePalette(24,palette="polychrome")

### Run Seurat on GEX. Assign cell types, then use deconTx to create decontaminated count matrix
annotation <- readRDS('/data/BCI-CRC/SO/genomes/UCSC_hg38_annotation.rds')

### CRC01_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC01_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC01_LM/outs/atac_fragments.tsv.gz"
CRC01_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC01_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC01_LM@meta.data$Sample <- 'CRC01_LM'
CRC01_LM@meta.data$Patient <- 'CRC01'
CRC01_LM@meta.data$Therapy <- 'NAC'
CRC01_LM@meta.data$Tissue <- 'LM'
### CRC02_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC02_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC02_LM/outs/atac_fragments.tsv.gz"
CRC02_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC02_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC02_LM@meta.data$Sample <- 'CRC02_LM'
CRC02_LM@meta.data$Patient <- 'CRC02'
CRC02_LM@meta.data$Therapy <- 'NAC'
CRC02_LM@meta.data$Tissue <- 'LM'
### CRC03_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC03_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC03_LM/outs/atac_fragments.tsv.gz"
CRC03_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC03_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC03_LM@meta.data$Sample <- 'CRC03_LM'
CRC03_LM@meta.data$Patient <- 'CRC03'
CRC03_LM@meta.data$Therapy <- 'NAC'
CRC03_LM@meta.data$Tissue <- 'LM'
### CRC04_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC04_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC04_LM/outs/atac_fragments.tsv.gz"
CRC04_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC04_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC04_LM@meta.data$Sample <- 'CRC04_LM'
CRC04_LM@meta.data$Patient <- 'CRC04'
CRC04_LM@meta.data$Therapy <- 'NAC'
CRC04_LM@meta.data$Tissue <- 'LM'
### CRC05_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC05_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC05_LM/outs/atac_fragments.tsv.gz"
CRC05_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC05_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC05_LM@meta.data$Sample <- 'CRC05_LM'
CRC05_LM@meta.data$Patient <- 'CRC05'
CRC05_LM@meta.data$Therapy <- 'naive'
CRC05_LM@meta.data$Tissue <- 'LM'
### CRC06_LM
# load the RNA and ATAC data
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC06_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC06_LM/outs/atac_fragments.tsv.gz"
CRC06_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC06_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC06_LM@meta.data$Sample <- 'CRC06_LM'
CRC06_LM@meta.data$Patient <- 'CRC06'
CRC06_LM@meta.data$Therapy <- 'naive'
CRC06_LM@meta.data$Tissue <- 'LM'
### CRC07_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC07_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC07_LM/outs/atac_fragments.tsv.gz"
CRC07_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC07_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC07_LM@meta.data$Sample <- 'CRC07_LM'
CRC07_LM@meta.data$Patient <- 'CRC07'
CRC07_LM@meta.data$Therapy <- 'naive'
CRC07_LM@meta.data$Tissue <- 'LM'
### CRC08_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC08_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC08_LM/outs/atac_fragments.tsv.gz"
CRC08_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC08_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC08_LM@meta.data$Sample <- 'CRC08_LM'
CRC08_LM@meta.data$Patient <- 'CRC08'
CRC08_LM@meta.data$Therapy <- 'naive'
CRC08_LM@meta.data$Tissue <- 'LM'
### CRC09_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC09_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC09_LM/outs/atac_fragments.tsv.gz"
CRC09_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC09_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC09_LM@meta.data$Sample <- 'CRC09_LM'
CRC09_LM@meta.data$Patient <- 'CRC09'
CRC09_LM@meta.data$Therapy <- 'naive'
CRC09_LM@meta.data$Tissue <- 'LM'
### CRC10_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC10_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC10_LM/outs/atac_fragments.tsv.gz"
CRC10_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC10_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC10_LM@meta.data$Sample <- 'CRC10_LM'
CRC10_LM@meta.data$Patient <- 'CRC10'
CRC10_LM@meta.data$Therapy <- 'naive'
CRC10_LM@meta.data$Tissue <- 'LM'
### CRC11_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC11_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC11_LM/outs/atac_fragments.tsv.gz"
CRC11_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC11_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC11_LM@meta.data$Sample <- 'CRC11_LM'
CRC11_LM@meta.data$Patient <- 'CRC11'
CRC11_LM@meta.data$Therapy <- 'naive'
CRC11_LM@meta.data$Tissue <- 'LM'
### CRC12_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC12_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC12_LM/outs/atac_fragments.tsv.gz"
CRC12_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC12_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC12_LM@meta.data$Sample <- 'CRC12_LM'
CRC12_LM@meta.data$Patient <- 'CRC12'
CRC12_LM@meta.data$Therapy <- 'naive'
CRC12_LM@meta.data$Tissue <- 'LM'
### CRC13_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC13_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC13_LM/outs/atac_fragments.tsv.gz"
CRC13_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC13_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC13_LM@meta.data$Sample <- 'CRC13_LM'
CRC13_LM@meta.data$Patient <- 'CRC13'
CRC13_LM@meta.data$Therapy <- 'NAC'
CRC13_LM@meta.data$Tissue <- 'LM'
### CRC14_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC14_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC14_LM/outs/atac_fragments.tsv.gz"
CRC14_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC14_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC14_LM@meta.data$Sample <- 'CRC14_LM'
CRC14_LM@meta.data$Patient <- 'CRC14'
CRC14_LM@meta.data$Therapy <- 'NAC'
CRC14_LM@meta.data$Tissue <- 'LM'
### CRC15_LM
counts <- Read10X_h5("/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC15_LM/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/CRC15_LM/outs/atac_fragments.tsv.gz"
CRC15_LM <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA")
CRC15_LM[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath, annotation = annotation)
CRC15_LM@meta.data$Sample <- 'CRC15_LM'
CRC15_LM@meta.data$Patient <- 'CRC15'
CRC15_LM@meta.data$Therapy <- 'NAC'
CRC15_LM@meta.data$Tissue <- 'LM'

# merge
CRC_LM <- merge(CRC01_LM, y=c(CRC02_LM,CRC03_LM,CRC04_LM,CRC05_LM,CRC06_LM,CRC07_LM,CRC08_LM,CRC09_LM,
                              CRC10_LM,CRC11_LM,CRC12_LM,CRC13_LM,CRC14_LM,CRC15_LM),
                add.cell.ids = c("CRC01_LM","CRC02_LM","CRC03_LM","CRC04_LM","CRC05_LM",
                                 "CRC06_LM","CRC07_LM","CRC08_LM","CRC09_LM","CRC10_LM",
                                 "CRC11_LM","CRC12_LM","CRC13_LM","CRC14_LM","CRC15_LM"),
                project = "CRC_LM")
saveRDS(CRC_LM, 'CRC_unfiltered.rds')

CRC_LM

DefaultAssay(CRC_LM) <- "RNA"

### Rename cells from Seurat format to ArchR format
proj <- loadArchRProject(path = "ArchRSavedProject_allCells/", showLogo = F)
proj
length(proj$cellNames)
length(Cells(CRC_LM))
head(proj$cellNames,2)
head(colnames(CRC_LM),2)
newNames <- gsub('_LM_','_LM#',colnames(CRC_LM))
newNames <- gsub('_T_','_T#', newNames)
head(newNames,2)
CRC_LM <- RenameCells(CRC_LM, new.names = newNames)
head(colnames(CRC_LM),2)
common_cells <- intersect(Cells(CRC_LM),proj$cellNames)
length(common_cells)
head(setdiff(Cells(CRC_LM),proj$cellNames))
head(setdiff(proj$cellNames,Cells(CRC_LM)))
### filter Seurat object by ArchR project cell names to keep cells passing ATAC QC
proj <- proj[common_cells,]
proj #
CRC_LM <- CRC_LM[,proj$cellNames]
CRC_LM #
head(colnames(CRC_LM))
head(proj$cellNames)
head(CRC_LM@meta.data,2)
CRC_LM@meta.data$TSSEnrichment <- proj$TSSEnrichment
CRC_LM@meta.data$nFrags <- proj$nFrags

CRC_LM[["percent.mt"]] <- PercentageFeatureSet(CRC_LM, pattern = "^MT-")
CRC_LM[["percent.ribo"]] <- PercentageFeatureSet(CRC_LM, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
head(CRC_LM@meta.data,2)
saveRDS(CRC_LM, 'CRC_GEX_unfiltered.rds')

pdf("Plots_Seurat/1-1_QC_vln.pdf", height=5, width=9)
VlnPlot(CRC_LM, features="TSSEnrichment", pt.size=0, group.by='Sample', y.max=10, cols=colPal)
VlnPlot(CRC_LM, features="nFrags", pt.size=0, group.by='Sample', y.max=10000, cols=colPal)
VlnPlot(CRC_LM, features="nCount_RNA", pt.size=0, group.by='Sample', y.max=20000, cols=colPal)
VlnPlot(CRC_LM, features="nFeature_RNA", pt.size=0, group.by='Sample', y.max=10000, cols=colPal)
VlnPlot(CRC_LM, features="percent.mt", pt.size=0, group.by='Sample', y.max=50, cols=colPal)
VlnPlot(CRC_LM, features="percent.ribo", pt.size=0, group.by='Sample', y.max=50, cols=colPal)
dev.off()

### remove mito and ribo genes
CRC_LM
genesToKeep <- grep("^MT-|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", row.names(CRC_LM), value = T, invert=T)
CRC_LM <- CRC_LM[genesToKeep,]
CRC_LM #24404 cells

### subset by GEX QC metrics
CRC_LM <- subset(CRC_LM, subset = nFeature_RNA > 300 & percent.mt < 10)

pdf("Plots_Seurat/1-2_QC_vln.pdf", height=5, width=9)
VlnPlot(CRC_LM, features="TSSEnrichment", pt.size=0, group.by='Sample', y.max=10, cols=colPal)
VlnPlot(CRC_LM, features="nFrags", pt.size=0, group.by='Sample', y.max=10000, cols=colPal)
VlnPlot(CRC_LM, features="nCount_RNA", pt.size=0, group.by='Sample', y.max=20000, cols=colPal)
VlnPlot(CRC_LM, features="nFeature_RNA", pt.size=0, group.by='Sample', y.max=8000, cols=colPal)
VlnPlot(CRC_LM, features="percent.mt", pt.size=0, group.by='Sample', y.max=20, cols=colPal)
VlnPlot(CRC_LM, features="percent.ribo", pt.size=0, group.by='Sample', y.max=50, cols=colPal)
dev.off()

saveRDS(CRC_LM, 'CRC_GEX_filtered.rds')

### umap pre-integration
DefaultAssay(CRC_LM) <- "RNA"
CRC_LM <- NormalizeData(CRC_LM)
CRC_LM <- FindVariableFeatures(CRC_LM, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(CRC_LM)
# alternatively can use VariableFeatures(CRC_LM) to scale variable genes only
# vars.to.regress = c("Gex_nGenes")
CRC_LM <- ScaleData(CRC_LM, features = VariableFeatures(CRC_LM))
CRC_LM <- RunPCA(CRC_LM, features = VariableFeatures(object = CRC_LM))
CRC_LM <- FindNeighbors(CRC_LM, dims = 1:40)
CRC_LM <- FindClusters(CRC_LM, resolution = 0.5, algorithm = 3)
CRC_LM <- RunUMAP(CRC_LM, dims = 1:40)

pdf("Plots_Seurat/1-3_umap.pdf", height = 5, width = 5)
DimPlot(CRC_LM, reduction = "umap", label = T)
DimPlot(CRC_LM, reduction = "umap", group.by = "Sample", shuffle=T, cols=colPal)
dev.off()

### integrate, first split data by sample
CRC_LM.list <- SplitObject(CRC_LM, split.by = "Sample")
CRC_LM.list
# normalize and identify variable features for each dataset independently
CRC_LM.list <- lapply(X = CRC_LM.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = CRC_LM.list)
### Perform integration #We then identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
CRC_LM.anchors <- FindIntegrationAnchors(object.list = CRC_LM.list, anchor.features = features, dims = 1:30)
CRC_LM.integrated <- IntegrateData(anchorset = CRC_LM.anchors, dims = 1:30)
### Perform integrated analysis on all cells
DefaultAssay(CRC_LM.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
length(VariableFeatures(CRC_LM.integrated))
#vars.to.regress = c("Gex_nGenes"),
CRC_LM.integrated <- ScaleData(CRC_LM.integrated, verbose = T, features = VariableFeatures(CRC_LM.integrated))
CRC_LM.integrated <- RunPCA(CRC_LM.integrated, npcs = 30, verbose = FALSE)

pdf("Plots_Seurat/1-4_elbow.pdf", height = 5, width = 5)
ElbowPlot(CRC_LM.integrated, ndims = 30)
dev.off()

CRC_LM.integrated <- RunUMAP(CRC_LM.integrated, reduction = "pca", dims = 1:30)
CRC_LM.integrated <- FindNeighbors(CRC_LM.integrated, reduction = "pca", dims = 1:30)
CRC_LM.integrated <- FindClusters(CRC_LM.integrated, resolution = 0.5, algorithm = 3)
saveRDS(CRC_LM.integrated, 'CRC_GEX_integrated.rds')

### umap plots
pdf("Plots_Seurat/2_integrated_umap.pdf", height = 5, width = 5)
DimPlot(CRC_LM.integrated, reduction = "umap", label = T, cols=colPal)
DimPlot(CRC_LM.integrated, reduction = "umap", group.by = "Sample", shuffle=T, cols=colPal)
FeaturePlot(CRC_LM.integrated, reduction = "umap", features = "nCount_RNA", max.cutoff = 20000)
FeaturePlot(CRC_LM.integrated, reduction = "umap", features = "nFeature_RNA", max.cutoff = 6000)
FeaturePlot(CRC_LM.integrated, reduction = "umap", features = "percent.mt")
FeaturePlot(CRC_LM.integrated, reduction = "umap", features = "percent.ribo")
dev.off()

DefaultAssay(CRC_LM.integrated) <- "RNA"
genes = c('EPCAM','ELF3','HNF4A','ALB','CP','DEFB1',
          'SLC11A1','CD68','LYZ','CD79A','BANK1','BLK','CD3E','PTPRC','TRAC',
          'DCN','ACTA2','COL3A1','PECAM1','ARL15','VWF','PROX1')
pdf("Plots_Seurat/3_integrated_umap_markers.pdf", height = 24, width = 20)
FeaturePlot(CRC_LM.integrated, reduction = "umap", features = genes,
            order = T, min.cutoff = 1)
dev.off()

##### next run deconTx on GEX.
