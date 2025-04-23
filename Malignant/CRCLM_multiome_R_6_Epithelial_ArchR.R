library('ArchR')
library('Seurat')
library('Signac')
set.seed(1)
addArchRThreads(threads = 4)
addArchRGenome("hg38")
library(dplyr)
library(stringr)
library(parallel)

proj <- loadArchRProject(path = "./ArchRSavedProject1_allCells/", force = FALSE, showLogo = FALSE)
proj

meta <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Epithelial_scvi_annotations_metadata.csv', header = T, row.names=1)
dim(meta)

proj <- subsetArchRProject(ArchRProj = proj,
                            cells = row.names(meta),
                            outputDirectory = "ArchRSavedProject2_Epithelial",
                            dropCells = TRUE, logFile = NULL,threads = getArchRThreads(), force = T)

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject2_Epithelial/", load = FALSE)
proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj

### add metadata to ArchR
proj$nFeature_RNA <- meta[['nFeature_RNA']]
proj$leiden <- meta[['leiden']]
proj$Cell_subtype <- meta[['Cell_subtype']]

### add embeddings
X_umap <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Epithelial_scvi_annotations_Xumap.csv', header = F)
X_scVI <- read.csv('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Epithelial_scvi_annotations_XscVI.csv', header = F)

df_umap <- as.data.frame(X_umap)
colnames(df_umap) <- c("#UMAP1", "#UMAP2")
row.names(df_umap) <- row.names(meta)
proj@embeddings$umap_Epi_scVI <- SimpleList(df = df_umap, params = list())

df_scVI <- as.data.frame(X_scVI)
row.names(df_scVI) <- row.names(meta)
latents = 20
colnames(df_scVI) <- paste0("scVI_", 1:latents)
proj@reducedDims$integrated_scVI_Epi <- SimpleList(df = df_scVI, params = list())
proj@reducedDims$integrated_scVI_Epi$corToDepth <- proj@reducedDims$LSI_Combined$corToDepth
proj@reducedDims$integrated_scVI_Epi$corToDepth$scaled <- proj@reducedDims$LSI_Combined$corToDepth$scaled[1:latents]
proj@reducedDims$integrated_scVI_Epi$corToDepth$none <- proj@reducedDims$LSI_Combined$corToDepth$none[1:latents]
proj@reducedDims$integrated_scVI_Epi$scaleDims <- proj@reducedDims$LSI_Combined$scaleDims

##### Recall peaks
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cell_subtype", force=T)
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(ArchRProj = proj,groupBy = "Cell_subtype", pathToMacs2 = pathToMacs2,
                               cutOff = 0.00001)

# ### Access peakset:
# temp <- proj@peakSet[(elementMetadata(proj@peakSet)[, "score"] > 10)]
# temp$Cell_subtype <- names(temp)
# temp_df <- as.data.frame(unname(temp))
# as.data.frame(table(temp_df$peakType, temp_df$Cell_subtype))
# ggplot(as.data.frame(table(temp_df$peakType, temp_df$Cell_subtype)), aes(fill=Var1, y=Freq, x=Var2)) +
#   geom_bar(position="stack", stat="identity")
# 
# mean(elementMetadata(proj@peakSet)[, "score"])
# quantile(elementMetadata(proj@peakSet)[, "score"], probs=seq(0, 1, 0.1))
# proj@peakSet <- proj@peakSet[(elementMetadata(proj@peakSet)[, "score"] > 10)]

proj <- addPeakMatrix(proj)

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject2_Epithelial/", load = FALSE)
proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj

### extract bed file from peaksset (granges object)
df <- data.frame(seqnames=seqnames(proj@peakSet),
                 starts=start(proj@peakSet)-1,
                 ends=end(proj@peakSet),
                 names=c(rep(".", length(proj@peakSet))),
                 scores=c(rep(".", length(proj@peakSet))),
                 strands=strand(proj@peakSet))
write.table(df, file="Epithelial_peaks/Epithelial_unionPeakset.bed", quote=F, sep="\t", row.names=F, col.names=F)

##### Peak heatmap
markersPeaks <- getMarkerFeatures(
 ArchRProj = proj, useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon",
 groupBy = "Cell_subtype")
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.75")
markerList$ipEMT[1:5,]

### write bed file from peaksset (granges object)
df <- markerList$ipEMT[,c(1,3,4)]
df$stand <- '+'
df
write.table(df, file="./ArchRSavedProject2_Epithelial/ipEMT_specificPeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

heatmapPeaks <- plotMarkerHeatmap(
 seMarker = markersPeaks,
 cutOff = "FDR < 0.1 & Log2FC > 0.5",
 transpose = F
)
plotPDF(heatmapPeaks, name = "Epithelial-1-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
### GEX marker genes
markerGenes = c('LGR5','ASCL2','SMOC2','PROM1','TFF3','MUC2','ATOH1','SLC26A3','KRT19','KRT20','FABP1','CA2','MMP7',
                'KRT7','KRT17','PLAUR','ANXA1','LAMC2','EMP1','LAMA3','L1CAM','IFI6','IFI44','IFIT1','IFIT3','MKI67',
                'TOP2A','PCNA','VEGFA','LDLR','CHGA','CHGB','NEUROD1','LRMP')
markersGEX <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneExpressionMatrix",
                              groupBy = "Cell_subtype", bias = "log10(nFeature_RNA)",
                              testMethod = "wilcoxon")
saveRDS(markersGEX, file='markers_Epithelial_GEX.rds')
markerList <- getMarkers(markersGEX, cutOff = "FDR < 0.01 & Log2FC > 1.0")
heatmapGEX <- plotMarkerHeatmap(
 seMarker = markersGEX,
 cutOff = "FDR < 0.01 & Log2FC > 1.0",
 labelMarkers = markerGenes,
 transpose = F,
 clusterCols = TRUE)
ComplexHeatmap::draw(heatmapGEX, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGEX, name = "Epithelial-1-GeneExpression-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

p <- plotEmbedding(proj, name = 'Cell_subtype', embedding = "umap_Epi_scVI",
                   size = 1.5, labelAsFactors=F, plotAs = 'points',
                   imputeWeights = NULL)
p
head(proj@cellColData,1)

proj$Cell_subtype <- as.character(proj$Cell_subtype)
p <- plotBrowserTrack(
 ArchRProj = proj,
 groupBy = "Cell_subtype",
 geneSymbol = c("LGR5",'EMP1','LAMC2','SMOC2',"HNF4A",'CDX2',"ASCL2","MUC2","SLC26A3","KRT19",'PLAUR','ANXA1'),
 #features =  getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0", returnGR = TRUE)["0"],
 upstream = 50000,
 downstream = 200000
)
grid::grid.newpage()
grid::grid.draw(p$LGR5)
plotPDF(p, name = "Epithelial-1-Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj
### motif enrichment in differential peaks
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = T)

########### Chrom-var
if("Motif" %ni% names(proj@peakAnnotation)){
 proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

proj <- addBgdPeaks(proj, force = T)

proj <- addDeviationsMatrix(
 ArchRProj = proj,
 peakAnnotation = "Motif",
 force = T)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Epithelial-3-chromVar-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
VarDevDF <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE, n=25)

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject2_Epithelial/", load = FALSE)
proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj

motifs <- c("HNF4A",'ASCL2','ATOH1','NFKB1','JUN')
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Cell_subtype", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = NULL)

plotPDF(p, name = "Epithelial-3-chromVar-Groups", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

# plot motif deviations on UMAP
p <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", 
  name = sort(markerMotifs), embedding = "umap_Epi_scVI",
  imputeWeights = getImputeWeights(proj), plotAs = 'points')

plotPDF(p, name = "Epithelial-3-chromVar-Deviation-UMAP", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

# check mRNA expression of these TFs

markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneExpressionMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("HNF4A-AS1","CEBPA-DT")]
markerRNA

proj <- addImputeWeights(proj, reducedDims = 'integrated_scVI_Epi', dimsToUse=1:20)

p1 <- plotEmbedding(
  ArchRProj = proj, colorBy = "GeneExpressionMatrix", 
  name = sort(markerRNA), embedding = "umap_Epi_scVI",
  imputeWeights = getImputeWeights(proj), plotAs = 'points')

p2 <- plotEmbedding(
  ArchRProj = proj, colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), embedding = "umap_Epi_scVI",
  imputeWeights = getImputeWeights(proj), plotAs = 'points')

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", 
                   name = sort(markerMotifs), embedding = "umap_Epi_scVI",
                   imputeWeights = getImputeWeights(proj), plotAs = 'points')
plotPDF(p1,p2,p3, name = "Epithelial-3-chromVar-UMAP-RNAandATACscore", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

### Identification of positive TF regulators
#get deviant motifs
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Cell_subtype")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix",
  reducedDims = "integrated_scVI_Epi", dimsToUse = 1:20, removeFromName1=NULL)

corGEM_MM <- correlateMatrices(
  ArchRProj = proj, useMatrix1 = "GeneExpressionMatrix", useMatrix2 = "MotifMatrix",
  reducedDims = "integrated_scVI_Epi", dimsToUse = 1:20, removeFromName1=NULL)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corCutOff <- 0.3
corGSM_MM$TFRegulator[which(corGSM_MM$cor > corCutOff & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "ACTIVATOR"
corGSM_MM$TFRegulator[which(corGSM_MM$cor < -1*corCutOff & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "REPRESSOR"
sort(corGSM_MM[corGSM_MM$TFRegulator=="ACTIVATOR",1])
sort(corGSM_MM[corGSM_MM$TFRegulator=="REPRESSOR",1])

write.csv(corGEM_MM, file='correlationMatrix_GEX_MM.csv', row.names = T)
library(ggrepel)

p1 <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point(size = 1) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dotted") +
  geom_vline(xintercept = corCutOff, lty = "dotted") + 
  geom_vline(xintercept = -corCutOff, lty = "dotted") + 
  geom_hline(yintercept = quantile(corGSM_MM$maxDelta, 0.75), lty = "dotted") + 
  scale_color_manual(values = c("REPRESSOR"="royalblue","ACTIVATOR"="red")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  geom_text_repel( 
    data = subset(data.frame(corGSM_MM), TFRegulator == 'ACTIVATOR'),
    aes(x = cor, y = maxDelta, label = GeneScoreMatrix_name),
    nudge_x = 0.01, nudge_y = 0.01, size = 1.25, max.overlaps = 20
  ) +
  geom_text_repel( 
    data = subset(data.frame(corGSM_MM), TFRegulator == 'REPRESSOR'),
    aes(x = cor, y = maxDelta, label = GeneScoreMatrix_name),
    nudge_x = 0.01, nudge_y = 0.01, size = 1.25, max.overlaps = 20
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p1

corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM[,"MotifMatrix_name"]))), ]
corGEM_MM$TFRegulator <- "NO"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > corCutOff & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Activator"
corGEM_MM$TFRegulator[which(corGEM_MM$cor < -corCutOff & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Repressor"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > -corCutOff & corGEM_MM$cor < corCutOff & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.97))] <- "No correlation"
sort(corGEM_MM[corGEM_MM$TFRegulator=="Activator",1])
sort(corGEM_MM[corGEM_MM$TFRegulator=="Repressor",1])
sort(corGEM_MM[corGEM_MM$TFRegulator=="No correlation",1])

p2 <- ggplot(data.frame(corGEM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point(size = 1) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dotted") +
  geom_vline(xintercept = corCutOff, lty = "dotted") + 
  geom_vline(xintercept = -corCutOff, lty = "dotted") +
  geom_hline(yintercept = quantile(corGEM_MM$maxDelta, 0.75), lty = "dotted") + 
  scale_color_manual(values = c("Repressor"="royalblue","Activator"="red", "No correlation"="black")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  geom_text_repel( 
    data = subset(data.frame(corGEM_MM), TFRegulator == 'Activator'),
    aes(x = cor, y = maxDelta, label = GeneExpressionMatrix_name),
    nudge_x = 0.01, nudge_y = 0.01, size = 1.25, max.overlaps = 20
  ) +
  geom_text_repel( 
    data = subset(data.frame(corGEM_MM), TFRegulator == 'Repressor'),
    aes(x = cor, y = maxDelta, label = GeneExpressionMatrix_name),
    nudge_x = 0.01, nudge_y = 0.01, size = 1.25, max.overlaps = 20
  ) +
  geom_text_repel( 
    data = subset(data.frame(corGEM_MM), TFRegulator == 'No correlation'),
    aes(x = cor, y = maxDelta, label = GeneExpressionMatrix_name),
    nudge_x = 0.01, nudge_y = 0.01, size = 1.25, max.overlaps = 20
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGEM_MM$maxDelta)*1.05)
  )
p2

plotPDF(p1,p2, name = "Epithelial-4-ChromVar-Correlation-MotifDelta-GeneExpressionOrScore", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

motifs_a <- sort(corGEM_MM[corGEM_MM$TFRegulator=="Activator",1])
motifs_a
#motifs_r <- sort(corGEM_MM[corGEM_MM$TFRegulator=="Repressor",1])
#motifs_r
motifs <- c(motifs_a)#, motifs_r)
motifs
markerMotifs_a <- getFeatures(proj, select = paste(motifs_a, collapse="|"), useMatrix = "MotifMatrix")
#markerMotifs_r <- getFeatures(proj, select = paste(motifs_r, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs_a <- grep("z:", markerMotifs_a, value = TRUE)
#markerMotifs_r <- grep("z:", markerMotifs_r, value = TRUE)
markerMotifs_a
#markerMotifs_r
#if need to remove a motif:
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs <- c(sort(markerMotifs_a))# #sort(markerMotifs_r))

# plot motif deviations on UMAP
p <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", 
                   name = markerMotifs, embedding = "umap_Epi_scVI",
                   imputeWeights = NULL, plotAs = 'points')

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", 
                   name = markerMotifs, embedding = "umap_Epi_scVI",
                   imputeWeights = getImputeWeights(proj),
                   plotAs = 'points')

plotPDF(p, p1, name = "Epithelial-4-Variable-Motif-Deviation-TFregulators", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# check mRNA expression of these TFs

markerRNA_a <- getFeatures(proj, select = paste(motifs_a, collapse="|"), useMatrix = "GeneExpressionMatrix")
#markerRNA_r <- getFeatures(proj, select = paste(motifs_a, collapse="|"), useMatrix = "GeneExpressionMatrix")
markerRNA <- c(sort(markerRNA_a))#,sort(markerRNA_r))

p <- plotEmbedding(
  ArchRProj = proj, colorBy = "GeneExpressionMatrix", 
  name = motifs, embedding = "umap_Epi_scVI", imputeWeights = getImputeWeights(proj),
  plotAs = 'points', pal = ArchRPalettes$greenBlue, log2Norm = T)
p1 <- plotEmbedding(
  ArchRProj = proj, colorBy = "GeneExpressionMatrix", 
  name = motifs, embedding = "umap_Epi_scVI", imputeWeights = NULL,
  plotAs = 'points', pal = ArchRPalettes$greenBlue, log2Norm = T)

plotPDF(p1,p, name = "Epithelial-4-Variable-Motif-Deviation-TFregulators-RNA", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

##### Motif heatmap
markerMotifs <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix",groupBy = "Cell_subtype",
  testMethod = "wilcoxon", useSeqnames = 'z')
markerList <- getMarkers(markerMotifs, cutOff = "FDR < 0.05 & MeanDiff > 0.75")
df <-  as.data.frame(do.call(rbind, markerList))
length(unique(df$name))
grep(pattern = 'ASCL', x = df$name, value = T)
write.csv(x=df, file = 'differentialMotifs_chromVARz.csv')

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerMotifs, transpose = F, log2Norm = F, labelMarkers = gsub(x=sort(markerMotifs_a), pattern = "z:", replacement = "" ),
  cutOff = "FDR < 0.05 & MeanDiff > 0.75")
p <- plotEmbedding(proj, name = 'Cell_subtype', embedding = "umap_Epi_scVI",
                   size = 1.5, labelAsFactors=F, plotAs = 'points',
                   imputeWeights = NULL)
plotPDF(p,heatmapPeaks, name = "Epithelial-4-Variable-Motif-Deviation-heatmap", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject2_Epithelial/", load = FALSE)

proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
### Heatmaps:

### GEM
markerGenes = c('EPCAM','ELF3','HNF4A','LGR5','ASCL2','PROM1','SMOC2',
               'PLAUR','ANXA1','LAMC2','LAMA3','MMP7','IFI6','IFIT3','IFIT1','CXCL8','IFI44','KRT20',
               'KRT7','KRT17','KRT19',
               'MUC2','TFF3','ATOH1','CHGA','CHGB','NEUROD1',
               'SLC26A3','CA1','CA2','FABP1',
               'MKI67','PCNA','TOP2A',
               'VEGFA','LDLR', 'DUOXA1','DUOXA2','EMP1','KRT7','KRT19')
markersGEX <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneExpressionMatrix",
                                groupBy = "Cell_subtype", bias = c("nFeature_RNA"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGEX, cutOff = "FDR < 0.01 & Log2FC > 1.0")
heatmapGEX <- plotMarkerHeatmap(seMarker = markersGEX, cutOff = "FDR < 0.01 & Log2FC > 1.0", labelMarkers = markerGenes, 
                                transpose = F, clusterCols = TRUE, pal = ArchRPalettes$solarExtra, plotLog2FC = T) #show_column_dend = FALSE)
ComplexHeatmap::draw(heatmapGEX, heatmap_legend_side = "bot", show_column_dend = FALSE)

### GSM
markersGS <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneScoreMatrix",
                               groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR < 0.05 & Log2FC > 1.0")
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR < 0.05 & Log2FC > 1.0", labelMarkers = markerGenes,
                               transpose = F, clusterCols = TRUE, pal = ArchRPalettes$blueYellow)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot")

### peakMatrix
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.75")
do.call(rbind, markerList)
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.75",
                                  transpose = F, pal = ArchRPalettes$purpleOrange, labelRows=FALSE)
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot")#, Cell_type_legend_side = "bot")

markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.5")
DARs <- do.call(rbind, markerList)
write.csv(do.call(rbind, markerList), file = "EpithelialPeaks_DARs_FDR05_log2fc05.csv")
length(unique(row.names(DARs))) #1861
length(row.names(DARs)) #1909

##### motifMatrix
markerMotifs <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix",groupBy = "Cell_subtype",
                                  testMethod = "wilcoxon", useSeqnames = 'z')
markerMotifs

rowData(markerMotifs)$name <- str_split_fixed(rowData(markerMotifs)$name, '_', 2)[,1]
rowData(markerMotifs)$name
markerList <- getMarkers(markerMotifs, cutOff = "FDR < 0.05 & MeanDiff > 0.75")
heatmapMotifs <- plotMarkerHeatmap(
  seMarker = markerMotifs, transpose = F, log2Norm = F, labelMarkers = c(motifs_a,'POU2F3','ASCL2','ASCL1','IRF4','ATOH1'),
  cutOff = "FDR < 0.05 & MeanDiff > 0.75", pal = ArchRPalettes$horizonExtra)
heatmapMotifs
plotPDF(heatmapGEX,heatmapGS,heatmapPeaks,heatmapMotifs, name = "Epithelial-5-heatmaps", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)
plotPDF(heatmapMotifs, name = "Epithelial-5a-heatmaps", width = 7, height = 12, ArchRProj = proj, addDOC = FALSE)
plotPDF(heatmapMotifs, name = "Epithelial-5b-heatmaps", width = 6, height = 12, ArchRProj = proj, addDOC = FALSE)
plotPDF(heatmapGEX,heatmapGS,heatmapPeaks,heatmapMotifs, name = "Epithelial-5d-heatmaps", width = 6, height = 9, ArchRProj = proj, addDOC = FALSE)
### Get key regulators:
activators_genes <- sort(corGEM_MM[corGEM_MM$TFRegulator=="Activator",1])
repressors_genes <- sort(corGEM_MM[corGEM_MM$TFRegulator=="Repressor",1])
repressors_genes
genes_plot <- c(activators_genes)#, repressors_genes)

seMotifs <- getGroupSE(ArchRProj = proj, useMatrix="MotifMatrix", groupBy = "Cell_subtype")
seMotifs <- seMotifs[rowData(seMotifs)$seqnames == "z", ]
library(stringr)
rownames(seMotifs) <- str_split_fixed(rowData(seMotifs)$name, '_', 2)[,1]
seMotifs <- seMotifs[genes_plot,] #[-12] to remove HOXA11-AS - want genes not antisense transcripts!
seMotifs

library(circlize)
col_BuY = colorRamp2(c(-2, 2), c("blue", "yellow"))
col_horizonExtra = colorRamp2(c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2), ArchRPalettes$horizonExtra)
col_BuRd = colorRamp2(c(-2,0,2), c("blue","white","red"))
library(ComplexHeatmap)
row_dend = as.dendrogram(hclust(dist(assay(seMotifs))))
col_dend = as.dendrogram(hclust(dist(t(assay(seMotifs)))))

heatmapMotifs <- Heatmap(assay(seMotifs), cluster_rows = row_dend, cluster_columns = col_dend,
                        col = col_horizonExtra, show_row_names = T,
                        heatmap_legend_param = list(title = 'Chrom var z-score'),
                        column_title = 'Motif enrichment')

##### seRNA
seRNA <- getGroupSE(ArchRProj = proj, useMatrix="GeneExpressionMatrix", groupBy = "Cell_subtype", scaleTo = 10000)
rownames(seRNA) <- rowData(seRNA)$name
seRNA <- seRNA[genes_plot,]
mat_seRNA <- t(scale(t(log2(assay(seRNA)+1))))
heatmapMotifsRNA <- Heatmap(t(scale(t(log2(assay(seRNA)+1)))), cluster_rows = row_dend, cluster_columns = col_dend,
                         col = col_BuRd, show_row_names = T,
                         heatmap_legend_param = list(title = 'z-score'),
                         column_title = 'RNA')
plotPDF(heatmapMotifsRNA+heatmapMotifs, name = "Epithelial-5-2-motifs", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)

###############################
GSM <- getMatrixFromProject(ArchRProj=proj,useMatrix="GeneScoreMatrix",verbose=T,binarize = F)
GSM_df <- as.data.frame(rowData(GSM))
assays(GSM)$GeneScoreMatrix[1:4,1:4]
length(rowData(GSM)$name)
length(unique(rowData(GSM)$name))
setequal(rowData(GSM)$name,unique(rowData(GSM)$name))
rownames(GSM) <- rowData(GSM)$name
writeMM(assays(GSM)$GeneScoreMatrix,file='Matrices/ArchRmatrices_Epithelial/gsmMatrix.mtx')
write.csv(colnames(GSM), file = 'Matrices/ArchRmatrices_Epithelial/gsmMatrix_barcodes.csv', quote = F, row.names = F)
write.csv(rownames(GSM), file = 'Matrices/ArchRmatrices_Epithelial/gsmMatrix_features.csv', quote = F, row.names = F)

# getAvailableMatrices(proj)
# GEM <- getMatrixFromProject(ArchRProj=proj,useMatrix="GeneExpressionMatrix",verbose=T,binarize = F)
# assays(GEM)$GeneExpressionMatrix[1:4,1:4]
# length(rowData(GEM)$name)
# length(unique(rowData(GEM)$name))
# setequal(rowData(GEM)$name,unique(rowData(GEM)$name))
# writeMM(assays(GEM)$GeneScoreMatrix,file='Matrices/ArchRmatrices_Epithelial/gemMatrix.mtx')
# write.csv(colnames(GEM), file = 'Matrices/ArchRmatrices_Epithelial/gemMatrix_barcodes.csv', quote = F, row.names = F)
# write.csv(rownames(GEM), file = 'Matrices/ArchRmatrices_Epithelial/gemMatrix_features.csv', quote = F, row.names = F)

getAvailableMatrices(proj)
MM <- getMatrixFromProject(ArchRProj=proj,useMatrix="MotifMatrix",verbose=T,binarize = F)
assays(MM)$z[1:4,1:4]
length(rowData(MM)$name)
length(unique(rowData(MM)$name))
setequal(rowData(MM)$name,unique(rowData(MM)$name)) 
writeMM(assays(MM)$z,file='Matrices/ArchRmatrices_Epithelial/motifMatrix.mtx')
write.csv(colnames(MM), file = 'Matrices/ArchRmatrices_Epithelial/motifMatrix_barcodes.csv', quote = F, row.names = F)
write.csv(rownames(MM), file = 'Matrices/ArchRmatrices_Epithelial/motifMatrix_features.csv', quote = F, row.names = F)

PM <- getMatrixFromProject(ArchRProj=proj,useMatrix="PeakMatrix",verbose=T,binarize = F)
peak_df <- as.data.frame(rowRanges(PM))
rownames(PM) <- paste(peak_df$seqnames,peak_df$start,peak_df$end,sep = '_')
writeMM(assays(PM)$PeakMatrix,file='Matrices/ArchRmatrices_Epithelial/peakMatrix.mtx')
write.csv(colnames(PM), file = 'Matrices/ArchRmatrices_Epithelial/peakMatrix_barcodes.csv', quote = F, row.names = F)
write.csv(rownames(PM), file = 'Matrices/ArchRmatrices_Epithelial/peakMatrix_features.csv', quote = F, row.names = F)

PM <- getMatrixFromProject(ArchRProj=proj,useMatrix="PeakMatrix",verbose=T,binarize = T)
writeMM(assays(PM)$PeakMatrix,file='Matrices/ArchRmatrices_Epithelial/peakMatrix_binary.mtx')

### save bigwigs

table(proj@cellColData$Cell_subtype)

getGroupBW(ArchRProj = proj,normMethod = "ReadsInTSS",
  groupBy = "Cell_subtype",
  tileSize = 50, #width of bins
  maxCells = 1000, # max no. of cells
  ceiling = 4, #Maximum contribution of accessibility per cell in each tile.
)

#### DE analysis stem vs. ipEMT
### peakMatrix
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                  useGroups='ipEMT', bgdGroups = 'Stem')
pv <- plotMarkers(seMarker = markersPeaks, name = "ipEMT", cutOff = "FDR <= 0.01 & Log2FC >= 1 | FDR <= 0.01 & Log2FC <= -1", plotAs = "Volcano")
pv
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 1")
df <- markerList$ipEMT[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_ipEMTvsStem_specificIpEMT.bed", quote=F, sep="\t", row.names=F, col.names=F)
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC < -1")
df <- markerList$ipEMT[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_ipEMTvsStem_specificStem.bed", quote=F, sep="\t", row.names=F, col.names=F)

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                  useGroups='ipEMT', bgdGroups = 'pEMT')
pv2 <- plotMarkers(seMarker = markersPeaks, name = "ipEMT", cutOff = "FDR <= 0.01 & Log2FC >= 1 | FDR <= 0.01 & Log2FC <= -1", plotAs = "Volcano")
pv2
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 1")
df <- markerList$ipEMT[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_ipEMTvsPEMT_specificIpEMT.bed", quote=F, sep="\t", row.names=F, col.names=F)
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC < -1")
df <- markerList$ipEMT[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_ipEMTvsPEMT_specificPEMT.bed", quote=F, sep="\t", row.names=F, col.names=F)

markersGEX <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneExpressionMatrix",
                                groupBy = "Cell_subtype", bias = c("nFeature_RNA"), testMethod = "wilcoxon",
                                useGroups = 'ipEMT', bgdGroups = 'Stem')
pv3 <- plotMarkers(seMarker = markersGEX, name = "ipEMT", cutOff = "FDR <= 0.01 & Log2FC >= 1 | FDR <= 0.01 & Log2FC <= -1", plotAs = "Volcano")

markersGEX <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneExpressionMatrix",
                                groupBy = "Cell_subtype", bias = c("nFeature_RNA"), testMethod = "wilcoxon",
                                useGroups = 'ipEMT', bgdGroups = 'pEMT')
pv4 <- plotMarkers(seMarker = markersGEX, name = "ipEMT", cutOff = "FDR <= 0.01 & Log2FC >= 1 | FDR <= 0.01 & Log2FC <= -1", plotAs = "Volcano")
pv4

plotPDF(pv,pv2,pv3,pv4, name = "Epithelial-6-volcanos", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)

#######
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                  useGroups= 'Colonocyte', bgdGroups = 'Stem')
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.5")
df <- markerList$Colonocyte[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_ColonocytevsStem_specificIpEMT.bed", quote=F, sep="\t", row.names=F, col.names=F)

#######
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                  useGroups= 'Goblet', bgdGroups = 'Stem')
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.5")
df <- markerList$Goblet[,c(1,3,4)]
df$strand <- '+'
write.table(df, file="Epithelial_peaks/DARs_GobletvsStem_specificIpEMT.bed", quote=F, sep="\t", row.names=F, col.names=F)

#####
markerMotifs <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Cell_subtype",
                                  testMethod = "wilcoxon", useSeqnames = 'z',
                                  useGroups = 'ipEMT', bgdGroups = 'Stem')

assays(markerMotifs)$Log2FC <- assays(markerMotifs)$MeanDiff
markersGEX
pma <- plotMarkers(seMarker = markerMotifs, name = 'ipEMT', plotAs = "Volcano")
pma


##### motifMatrix
markerMotifs <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix",groupBy = "Cell_subtype",
                                  testMethod = "wilcoxon", useSeqnames = 'z')
markerMotifs
rowData(markerMotifs)$name <- str_split_fixed(rowData(markerMotifs)$name, '_', 2)[,1]
rowData(markerMotifs)$name
markerList <- getMarkers(markerMotifs, cutOff = "FDR < 0.01 & MeanDiff > 1")
heatmapMotifs <- plotMarkerHeatmap(
  seMarker = markerMotifs, transpose = F, log2Norm = F, labelMarkers = motifs_a,
  cutOff = "FDR < 0.05 & MeanDiff > 0.75", pal = ArchRPalettes$horizonExtra)

##### for sckinetics
### peakMatrix
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "Cell_subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",
                                  useGroups = c("Stem", "Colonocyte", "ipEMT", "Intermediate", "TA1", "TA2", "Hypoxia", "pEMT"),
                                  bgdGroups = "Stem (NOTUM high)"
                                  )
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.75")
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05")
DARs <- do.call(rbind, markerList)
DARs_dedup <- DARs[c(1,2,3,4)]
DARs_dedup
length(unique(row.names(DARs)))

as_tibble(DARs_dedup)
DARs_dedup <- distinct(as_tibble(DARs_dedup))
write.table(DARs_dedup, file = "Epithelial_peaks/DARs_vsStemNOTUM.txt", sep="\t")
DARs_dedup


###
proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = FALSE)
proj

proj <- proj[row.names(subset(proj@cellColData, Cell_subtype %in% c('ipEMT','pEMT', 'Stem (NOTUM high)', 'Stem'))),]

addArchRThreads(threads = 1)
markerMotifs <- getMarkerFeatures(ArchRProj = proj, useMatrix = "MotifMatrix",groupBy = "Cell_subtype",
                                  testMethod = "wilcoxon", useSeqnames = 'z')
markerList <- getMarkers(markerMotifs, cutOff = "FDR < 0.05 & MeanDiff > 0.75")

heatmapMotifs <- plotMarkerHeatmap(
  seMarker = markerMotifs, transpose = F, log2Norm = F, scaleRows = F,
  #labelMarkers = gsub(x=sort(markerMotifs_a), pattern = "z:", replacement = "" ),
  cutOff = "FDR < 0.05 & MeanDiff > 0.75",)

heatmapMotifs

grep('ATF3', markerList$pEMT$name, value = TRUE)

df <-  as.data.frame(do.call(rbind, markerList))
length(unique(df$name))

grep('ASCL2_', markerList$`Stem`$name, value = TRUE)

