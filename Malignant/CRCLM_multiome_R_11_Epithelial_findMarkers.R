setwd("/data/BCI-CRC/SO/data/CRC_multiome/ArchR_final_analysis/")
library('ArchR')
library('Seurat')
#library('Signac')
set.seed(1)
addArchRThreads(threads = 2)
addArchRGenome("hg38")
library(dplyr)
#library(stringr)
library(parallel)

### Find DEGs for mCRC

metadata <- read.csv("/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Matrices/Epithelial_scvi_annotations_metadata.csv", row.names=1)
features <- read.csv("/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Matrices/Epithelial_scvi_annotations_features.csv")

CRC.decon.counts <- readRDS('CRC_GEX_decon_counts.rds')
CRC.sce <- readRDS('CRC_GEX_decon.sce.rds')

CRC <- CreateSeuratObject(counts = CRC.decon.counts, meta.data = as.data.frame(colData(CRC.sce)))
CRC <- CRC[row.names(features),row.names(metadata)]
head(CRC@meta.data[,1:2])
head(metadata[,1:2])
CRC@meta.data <- metadata
CRC <- NormalizeData(CRC)

saveRDS(CRC, file="CRC_Epithelial.rds")

CRC <- readRDS("CRC_Epithelial.rds")

Idents(CRC) <- "Cell_subtype"
markers <- FindAllMarkers(CRC, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group

library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs_seurat_wilcoxon.xlsx")

DE_stem <- FindMarkers(CRC, ident.1 = c('Stem', 'Stem (NOTUM high)'),only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_stem <- DE_stem %>%
  arrange(desc(avg_log2FC))

DE_ipEMT <- FindMarkers(CRC, ident.1 = 'ipEMT',only.pos = TRUE,
                        min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_ipEMT <- DE_ipEMT %>%
  arrange(desc(avg_log2FC))

DE_pEMT <- FindMarkers(CRC, ident.1 = 'pEMT',only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_pEMT <- DE_pEMT %>%
  arrange(desc(avg_log2FC))

write.csv(DE_stem, file="DEGs/DEGs_seurat_wilcoxon_stem.csv")
write.csv(DE_ipEMT, file="DEGs/DEGs_seurat_wilcoxon_ipEMT.csv")
write.csv(DE_pEMT, file="DEGs/DEGs_seurat_wilcoxon_pEMT.csv")

### Find DEGs between mCRC ipEMT and pEMT
CRC <- readRDS("CRC_Epithelial.rds")
Idents(CRC) <- "Cell_subtype"
DE_ipEMT_pEMT <- FindMarkers(CRC, ident.1 = 'ipEMT',ident.2='pEMT', logfc.threshold = -Inf,
                        min.pct = 0.1,  min.diff.pct = -Inf, only.pos = FALSE)
dim(DE_ipEMT_pEMT)
write.csv(DE_ipEMT_pEMT, file="DEGs/DEGs_seurat_wilcoxon_pEMT_vs_ipEMT_volcano.csv")

DE_ipEMT <- FindMarkers(CRC, ident.1 = 'pEMT',ident.2='ipEMT', only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_ipEMT <- DE_ipEMT %>%
  arrange(desc(avg_log2FC))
write.csv(DE_ipEMT, file="DEGs/DEGs_seurat_wilcoxon_ipEMT_vs_pEMT.csv")

DE_pEMT <- FindMarkers(CRC, ident.1 = 'pEMT',ident.2='ipEMT', only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_pEMT <- DE_pEMT %>%
  arrange(desc(avg_log2FC))
write.csv(DE_pEMT, file="DEGs/DEGs_seurat_wilcoxon_pEMT_vs_ipEMT.csv")

### Find DEGs between mCRC ipEMT and pEMT
CRC <- readRDS("CRC_Epithelial.rds")
Idents(CRC) <- "Cell_subtype"
DE_ipEMT_StemNOTUM <- FindMarkers(CRC, ident.1 = 'ipEMT',ident.2='Stem (NOTUM high)', logfc.threshold = -Inf,
                             min.pct = 0.1,  min.diff.pct = -Inf, only.pos = FALSE)
dim(DE_ipEMT_StemNOTUM)
write.csv(DE_ipEMT_StemNOTUM, file="DEGs/DEGs_seurat_wilcoxon_ipEMT_vs_StemNOTUM_volcano.csv")






table(CRC$Sample, CRC$Therapy)
Idents(CRC) <- "Therapy"
markers <- FindAllMarkers(CRC, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group
library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs/DEGs_seurat_wilcoxon_therapy.xlsx")

CRC@meta.data$Therapy_Elise_3 <- plyr::mapvalues(CRC@meta.data$Sample,
                                                 c('CRC01_LM','CRC02_LM','CRC03_LM','CRC04_LM','CRC05_LM','CRC06_LM','CRC07_LM','CRC08_LM','CRC09_LM','CRC10_LM','CRC11_LM','CRC12_LM','CRC13_LM','CRC14_LM','CRC15_LM'),
                                                 c('NAC','NAC','NAC','Treated','Naive','Naive','Naive','Naive','Naive','Naive','Naive','Naive','Treated','Treated','Treated')
)
table(CRC$Sample, CRC$Therapy_Elise_3)
Idents(CRC) <- "Therapy_Elise_3"
markers <- FindAllMarkers(CRC, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group
library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs/DEGs_seurat_wilcoxon_therapy_Elise_3.xlsx")

### Find DEGs for primary CRC

counts <- Matrix::readMM("/data/BCI-CRC/SO/data/public/primaryCRC_20mt/Harmony/Matrices/SMC_KUL_Pelka_Che_Wu_CRC_epithelial_regressCC_harmonyPatient_20_filtered_harmony_inferCNV_harmony_RAW.mtx")
metadata <- read.csv("/data/BCI-CRC/SO/data/public/primaryCRC_20mt/Harmony/Matrices/SMC_KUL_Pelka_Che_Wu_CRC_epithelial_regressCC_harmonyPatient_20_filtered_harmony_inferCNV_harmony_metadata.csv", row.names=1)
features <- read.csv("/data/BCI-CRC/SO/data/public/primaryCRC_20mt/Harmony/Matrices/SMC_KUL_Pelka_Che_Wu_CRC_epithelial_regressCC_harmonyPatient_20_filtered_harmony_inferCNV_harmony_features.csv")
dim(features)
dim(metadata)
dim(counts)

dimnames(counts) = list(row.names(metadata),features$X)
counts

pCRC <- CreateSeuratObject(counts = t(counts), meta.data = metadata)
pCRC

saveRDS(pCRC, file="primaryCRC_Epithelial.rds")

pCRC <- readRDS("primaryCRC_Epithelial.rds")

Idents(pCRC) <- "Cell_subtype"

markers <- FindAllMarkers(pCRC, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group

library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs/primaryEpithelial_DEGs_seurat_wilcoxon.xlsx")

DE_stem <- FindMarkers(pCRC, ident.1 = c('Stem', 'Stem (NOTUM high)'),only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_stem <- DE_stem %>%
  arrange(desc(avg_log2FC))

DE_ipEMT <- FindMarkers(pCRC, ident.1 = 'ipEMT',only.pos = TRUE,
                        min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_ipEMT <- DE_ipEMT %>%
  arrange(desc(avg_log2FC))

DE_pEMT <- FindMarkers(pCRC, ident.1 = 'pEMT',only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_pEMT <- DE_pEMT %>%
  arrange(desc(avg_log2FC))

write.csv(DE_stem, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_stem.csv")
write.csv(DE_ipEMT, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_ipEMT.csv")
write.csv(DE_pEMT, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_pEMT.csv")

### Find DEGs between pCRC ipEMT and pEMT
CRC <- readRDS("primaryCRC_Epithelial.rds")
Idents(CRC) <- "Cell_subtype"
DE_ipEMT <- FindMarkers(CRC, ident.1 = 'ipEMT',ident.2='pEMT', only.pos = TRUE,
                        min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_ipEMT <- DE_ipEMT %>%
  arrange(desc(avg_log2FC))
write.csv(DE_ipEMT, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_ipEMT_vs_pEMT.csv")
DE_pEMT <- FindMarkers(CRC, ident.1 = 'pEMT',ident.2='ipEMT', only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01)
DE_pEMT <- DE_pEMT %>%
  arrange(desc(avg_log2FC))
write.csv(DE_pEMT, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_pEMT_vs_ipEMT.csv")

# iREC vs. REC for volcano primary
DE_ipEMT_pEMT <- FindMarkers(CRC, ident.1 = 'ipEMT',ident.2='pEMT', logfc.threshold = -Inf,
                             min.pct = 0.1,  min.diff.pct = -Inf, only.pos = FALSE)
dim(DE_ipEMT_pEMT)
DE_ipEMT_pEMT <- DE_ipEMT_pEMT %>%
  arrange(desc(avg_log2FC))
write.csv(DE_ipEMT_pEMT, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_iREC_vs_REC_volcano.csv")

# iREC vs. Stem NOTUM for volcano primary
DE_ipEMT_StemNOTUM <- FindMarkers(CRC, ident.1 = 'ipEMT',ident.2='Stem (NOTUM high)', logfc.threshold = -Inf,
                             min.pct = 0.1,  min.diff.pct = -Inf, only.pos = FALSE)
dim(DE_ipEMT_StemNOTUM)
DE_ipEMT_StemNOTUM <- DE_ipEMT_StemNOTUM %>%
  arrange(desc(avg_log2FC))
write.csv(DE_ipEMT_StemNOTUM, file="DEGs/primaryEpithelial_DEGs_seurat_wilcoxon_iREC_vs_StemNOTUM_volcano.csv")

### Find DEGs for normal colon

counts <- Matrix::readMM("/data/BCI-CRC/SO/data/public/normal_colonFinal/Matrices/normalColon_harmony_annotations_RAW.mtx")
metadata <- read.csv("/data/BCI-CRC/SO/data/public/normal_colonFinal/Matrices/normalColon_harmony_annotationsmetadata.csv", row.names=1)
features <- read.csv("/data/BCI-CRC/SO/data/public/normal_colonFinal/Matrices/normalColon_harmony_annotationsfeatures.csv", row.names=1)
dim(features)
dim(metadata)
dim(counts)

dimnames(counts) = list(row.names(features),row.names(metadata))
counts[1:5,1:5]

colon <- CreateSeuratObject(counts = counts, meta.data = metadata)
colon

saveRDS(colon, file="normalColon_Epithelial.rds")

Idents(colon) <- "Cell_subtype"

markers <- FindAllMarkers(colon, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group
library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs/normalColon_DEGs_seurat_wilcoxon.xlsx")

### Find DEGs for integrated with public CRC

CRC <- readRDS('/data/BCI-CRC/Mirjana/CRC_metastasis/multiome_wang_sathe_che_mt.rds')
CRC

library(Matrix)
writeMM(GetAssayData(object = CRC, slot = "counts"), file = 'Matrices/multiome_wang_sathe_che_mt.mtx')
write.csv(CRC@meta.data, file='Matrices/multiome_wang_sathe_che_mt_metadata.csv')
write.csv(row.names(GetAssayData(object = CRC, slot = "counts")), file='Matrices/multiome_wang_sathe_che_mt_features.csv')
write.csv(as.data.frame(Embeddings(CRC[['umap']])), file = 'Matrices/multiome_wang_sathe_che_mt_umap.csv', quote = F)
write.csv(as.data.frame(Embeddings(CRC[['pca']])), file = 'Matrices/multiome_wang_sathe_che_mt_pca.csv', quote = F)

CRC@meta.data[1:2,]

Idents(CRC) <- "Annotation"

markers <- FindAllMarkers(CRC, test.use = 'wilcox', only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.01
)
markers$pct.diff <- markers$pct.1 - markers$pct.2
markers$cluster <- gsub(pattern = 'pEMT/ipEMT', replacement = 'REC_iREC', x = markers$cluster)
markers.group <- markers %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group = T)
markers.group

library(writexl)
write_xlsx(split(markers.group, markers.group$cluster), "DEGs/metastasis_publicIntegration_DEGs_seurat_wilcoxon.xlsx")


### DARs between iREC and REC

proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj

markersPeaks_iREC <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                       groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                       useGroups='ipEMT', bgdGroups = 'pEMT')
df <- data.frame(log2FC = assays(markersPeaks_iREC)$Log2FC, FDR = assays(markersPeaks_iREC)$FDR)
row.names(df) <- row.names(markersPeaks_iREC)
write.csv(df, 'DEGs/DARs_iREC_vs_REC.csv', quote=F)

markersPeaks_iREC <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                       groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                       useGroups='ipEMT', bgdGroups = 'Stem (NOTUM high)')
df <- data.frame(log2FC = assays(markersPeaks_iREC)$Log2FC, FDR = assays(markersPeaks_iREC)$FDR)
row.names(df) <- row.names(markersPeaks_iREC)
write.csv(df, 'DEGs/DARs_iREC_vs_StemNOTUM.csv', quote=F)



