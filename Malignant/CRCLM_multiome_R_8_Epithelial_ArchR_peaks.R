library('ArchR')
library('Seurat')
library('Signac')
set.seed(1)
addArchRThreads(threads = 1)
#addArchRThreads(threads = as.integer(nslots))
addArchRGenome("hg38")
library(dplyr)
library(stringr)
library(ggplot2)
library(parallel)

proj <- loadArchRProject(path = "./ArchRSavedProject2_Epithelial/", force = FALSE, showLogo = F)
proj

#proj@cellColData$Cell_subtype <- gsub('Stem (NOTUM high)', 'Stem NOTUM', proj@cellColData$Cell_subtype)
#proj@cellColData$Cell_subtype <- gsub('iREC', 'ipEMT', proj@cellColData$Cell_subtype)
#proj@cellColData$Cell_subtype <- gsub('REC', 'pEMT', proj@cellColData$Cell_subtype)

table(proj$Cell_subtype)

### Pie chart of peak types

peaks <- getPeakSet(proj)
peaks$Cell_subtype <- row.names(peaks)
write.csv(as.data.frame(peaks, row.names = NULL), file = 'Epithelial_peaks/Epithelial_unionPeakset.csv')

peaks <- getPeakSet(proj)
#peaks$peakType
as.data.frame(table(peaks$peakType))

peakType <- as.data.frame(prop.table(table(peaks$peakType)))
colnames(peakType) <- c('Peak_type','Freq')

# ggplot(as.data.frame(table(peaks$peakType)), aes(x="", y=Freq, fill=Var1)) +
#   geom_bar(stat="identity", width=1) +
#   coord_polar("y", start=0) +
#   theme_void() # remove background, grid, numeric labels

# Create a basic bar
pie = ggplot(peakType, aes(x="", y=Freq, fill=Peak_type)) + geom_bar(stat="identity", width=1)
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(Freq*100), "%")), position = position_stack(vjust = 0.5)) + theme_void()
pie = pie + scale_fill_manual(values=c("#6a3d9a", "#a6cee3", "#cab2d6", "#fb9a99"))

# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "Peak type")

# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5))

pdf('Epithelial_peaks/Epithelial_peaks_pie.pdf', height=3.5, width=3.5)
pie
dev.off()

proj <- addPeak2GeneLinks(ArchRProj = proj,reducedDims = "integrated_scVI_Epi", useMatrix = 'GeneExpressionMatrix',
                          dimsToUse = 1:20, k = 100
)
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSavedProject2_Epithelial/", load = FALSE)

### peak2gene linkage heatmap - activators!!!!!
sePeaks <- getGroupSE(ArchRProj = proj, useMatrix="PeakMatrix", groupBy = "Cell_subtype", scaleTo = 10000)
seRNA <- getGroupSE(ArchRProj = proj, useMatrix="GeneExpressionMatrix", groupBy = "Cell_subtype", scaleTo = 10000)

p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF$linkage <- paste(p2geneDF$peakName, p2geneDF$geneName, sep='-')
dim(p2geneDF)
head(p2geneDF)
p2geneDF$distNearestGeneTSS <- proj@peakSet$distToTSS[p2geneDF$idxATAC]
p2geneDF$nearestGene <- proj@peakSet$nearestGene[p2geneDF$idxATAC]

p2g_tss <- as.data.frame(getGenes(proj))
head(p2g_tss)
dim(p2g_tss)
library(dplyr)
p2g_tss <- p2g_tss %>% 
  mutate(TSS = if_else(strand == '+', start, end))
dim(p2g_tss)
head(p2g_tss)
dim(na.omit(p2g_tss))
p2g_tss <- na.omit(p2g_tss)
library(stringr)
p2geneDF$peakChr <- str_split_fixed(p2geneDF$peakName, '_', n=3)[,1]
p2geneDF$peakStart <- str_split_fixed(p2geneDF$peakName, '_', n=3)[,2]
p2geneDF$peakEnd <-  str_split_fixed(p2geneDF$peakName, '_', n=3)[,3]
p2geneDF$peakLength <- as.numeric(p2geneDF$peakEnd) - as.numeric(p2geneDF$peakStart)
min(p2geneDF$peakLength)
max(p2geneDF$peakLength)
p2geneDF$peakCentre <- as.numeric(p2geneDF$peakEnd) - 250
head(p2geneDF,3)

p2geneDF_TSS <- merge(x = p2geneDF,y = p2g_tss, by.x = 'geneName', by.y = 'symbol', all.x = TRUE)
dim(p2geneDF_TSS)
p2geneDF_TSS <- na.omit(p2geneDF_TSS)
dim(p2geneDF_TSS)
p2geneDF_TSS$PeakDistanceToTSS <- abs(as.numeric(p2geneDF_TSS$TSS)-p2geneDF_TSS$peakCentre)
plot(density(na.omit(p2geneDF_TSS$PeakDistanceToTSS)))
p2geneDF_TSS <- as.data.frame(p2geneDF_TSS) %>% 
  mutate(Linkage = if_else(Correlation > 0.40 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000 & VarQATAC > 0.25 & VarQRNA > 0.25 & FDR < 1e-04, 'Putative enhancer', 'No linkage'))

dim(subset(p2geneDF_TSS, PeakDistanceToTSS > 0 & PeakDistanceToTSS < 250000)) #613160
dim(subset(p2geneDF_TSS, PeakDistanceToTSS > 1000 & PeakDistanceToTSS < 250000)) #590687
dim(subset(p2geneDF_TSS, PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000)) #584132
#p2geneDF_TSS <- subset(p2geneDF_TSS, PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000)

p2geneDF_TSS$PeakDistanceToTSSKb <- p2geneDF_TSS$PeakDistanceToTSS/1000
PeakDistanceToTSS_activators <- subset(p2geneDF_TSS, Correlation > 0.40 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000 & VarQATAC > 0.25 & VarQRNA > 0.25 & FDR < 1e-04)
PeakDistanceToTSS_noLinkage <- subset(p2geneDF_TSS, Correlation <= 0.40 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000)
PeakDistanceToTSS_random <- sample_n(na.omit(subset(p2geneDF_TSS, PeakDistanceToTSS > 2000)), 5000)

dim(subset(p2geneDF_TSS, Correlation > 0.4 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000 & VarQATAC > 0.25 & VarQRNA > 0.25 & FDR < 1e-04))

write.csv(p2geneDF_TSS, 'Epithelial_peaks/peak2genes.csv', quote=F)
p2geneDF_TSS <- read.csv("Epithelial_peaks/peak2genes.csv")
p2geneDF_TSS

# test normality, p-value less than 0.05 indicates data is not normally distributed
shapiro.test(PeakDistanceToTSS_activators$PeakDistanceToTSS) # p-value < 2.2e-16
plot(density(PeakDistanceToTSS_activators$PeakDistanceToTSS))
wilcox.test(PeakDistanceToTSS_activators$PeakDistanceToTSS, PeakDistanceToTSS_noLinkage$PeakDistanceToTSS, alternative = "two.sided")
# W = 331633315, p-value < 2.2e-16
wilcox.test(PeakDistanceToTSS_activators$PeakDistanceToTSS, PeakDistanceToTSS_random$PeakDistanceToTSS, alternative = "two.sided")
# W = 2505551, p-value < 2.2e-16
wilcox.test(PeakDistanceToTSS_random$PeakDistanceToTSS, PeakDistanceToTSS_noLinkage$PeakDistanceToTSS, alternative = "two.sided")
# W = 1596234042, p-value = 6.247e-06

pdf('Epithelial_peaks/cumulativeDistribution.pdf', height=5, width=5)
ggplot(na.omit(p2geneDF_TSS), aes(x=Linkage, y=PeakDistanceToTSSKb)) + 
  geom_boxplot(notch=F, fill = c("#fb9a99","#6a3d9a")) +
  ylim(0,250) + 
  theme_bw() + ylab('Distance to TSS (kb)') + xlab('') + 
  theme(text = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

plot(ecdf(PeakDistanceToTSS_activators$PeakDistanceToTSSKb),do.points=FALSE, xlim = c(0,250), col='#6a3d9a', 
     axes = T, yaxt = "n", xlab='Distance (kb)', ylab='Cumulative distribution', main = "Distance to TSS")
plot(ecdf(PeakDistanceToTSS_noLinkage$PeakDistanceToTSSKb),do.points=FALSE, add=TRUE, col='#fb9a99', xlim = c(0,250))
plot(ecdf(PeakDistanceToTSS_random$PeakDistanceToTSSKb),do.points=FALSE, add=TRUE, col='black', xlim = c(0,250))
legend(x = 'right', legend = c('Putative enhancer','No linkage','random'), box.lty = 0, 
       text.col = c('#6a3d9a','#fb9a99','black'),
       bg = rgb(0,0,0, alpha = 0))
axis(side = 2, at = c(0,0.25,0.5,0.75,1.0))
dev.off()

p2geneDF_TSS_noLinkage <- subset(p2geneDF_TSS, Correlation <= 0.40 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000)
p2geneDF_TSS <- subset(p2geneDF_TSS, Correlation > 0.40 & PeakDistanceToTSS > 2000 & PeakDistanceToTSS < 250000 & VarQATAC > 0.25 & VarQRNA > 0.25 & FDR < 1e-04)

write.table(p2geneDF_TSS_noLinkage[,c('peakChr','peakStart','peakEnd','seqnames','TSS','TSS','Correlation','linkage')],
            file = 'Epithelial_peaks/Epithelial_noLinkage_p2g.bedpe', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(p2geneDF_TSS[,c('peakChr','peakStart','peakEnd','seqnames','TSS','TSS','Correlation','linkage')],
            file = 'Epithelial_peaks/Epithelial_activators_p2g.bedpe', quote = F, sep = '\t', row.names = F, col.names = F)

PeakDistanceToTSS_random_2 <- sample_n(PeakDistanceToTSS_random, dim(p2geneDF_TSS)[1])
write.table(PeakDistanceToTSS_random_2[,c('peakChr','peakStart','peakEnd','seqnames','start','start','Correlation','linkage')],
            file = 'Epithelial_peaks/Epithelial_random_p2g.bedpe', quote = F, sep = '\t', row.names = F, col.names = F)


### Use bedtools intersect on cluster: /data/BCI-CRC/SO/data/CRC_multiome/ArchR_final_analysis/Epithelial_peaks/intersectBed.sh
#https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
#assume a <= b
#sum(dhyper(t:b, a, n - a, b))
a <- dim(p2geneDF_TSS)[1] # no. of putative enhancers (p2g links) and no. of random peaks
b <- 33130 # no. of enhancers (Della Chiara)
n <- length(getPeakSet(proj)) # total no. peaks
t1 <- 802 # intersection of a and b
sum(dhyper(t1:b, a, n-a, b)) #P = 1.033799e-32

library(eulerr)

# intersection of random:
t2 <- 447
sum(dhyper(t2:b, a, n - a, b)) #P=1

# intersection of union peakset:
t <- 31171 # overlap of Dell Chiara enhancers and union peakset
a <- 33130 # no. of Della Chiara enhancers
b <- length(getPeakSet(proj)) # no. of peaks in union peakset, identical to n
n <- length(getPeakSet(proj)) # no. of peaks in union peakset
sum(dhyper(t:b, a, n-a, b))

### Pie chart overlaps
library(ggplot2)
df <- data.frame(c(t1,t2,a-t1,a-t2), c('Putative_enhancers','Random','Putative_enhancers','Random'),c('Yes','Yes', 'No','No'))
colnames(df) <- c('Value','Peakset','ChromHMM_enhancer')
df

library(ggplot2)
#df <- data.frame(c(t1,t2,t3,a-t1,a-t2,b-t3), c('Putative_enhancers','Random',"All_peaks",'Putative_enhancers','Random',"All_peaks"),c('Yes','Yes',"Yes","No",'No','No'))
#colnames(df) <- c('Value','Peakset','ChromHMM_enhancer')
#df

# Stacked + percent
pdf('Epithelial_peaks/intersect_p2g_DellaChiara.pdf', height = 5, width = 5)
ggplot(df, aes(fill=ChromHMM_enhancer, y=Value, x=Peakset)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#a6cee3","#6a3d9a", "#cab2d6", "#fb9a99")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black")
    )
ggplot(df, aes(fill=ChromHMM_enhancer, y=Value, x=Peakset)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#a6cee3","#6a3d9a", "#cab2d6", "#fb9a99")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black")
        )
dev.off()

########## 
########## 
########## 

peaks <- read.table(file = "Epithelial_peaks/Epithelial_unionPeakset.bed", header = F)
peaks$Peak_size <- peaks$V3 - peaks$V2
accessible_genome_size <- sum(peaks$Peak_size)

########## 
########## 
########## sePeaks

colnames(sePeaks) <- c('Colonocyte','Enteroendocrine','Goblet','Hypoxia','Intermediate','iREC','REC','Stem','Stem NOTUM',
                       'TA1','TA2','Tuft','UPR')
rowData(sePeaks) <- as.data.frame(rowData(sePeaks)) %>%
  mutate(peak = paste0(seqnames, "_", start, "_", end))
dim(rowData(sePeaks))
rownames(sePeaks) <- rowData(sePeaks)$peak
sePeaks <- sePeaks[p2geneDF_TSS$peakName,]
dim(rowData(sePeaks))
#rownames(p2geneDF_TSS) <- p2geneDF_TSS$peakName
#p2geneDF_TSS <- p2geneDF_TSS[rowData(sePeaks)$peak,]
rowData(sePeaks)$linkage <- p2geneDF_TSS$linkage
rownames(sePeaks) <- p2geneDF_TSS$linkage

########## seRNA
colnames(seRNA) <- c('Colonocyte','Enteroendocrine','Goblet','Hypoxia','Intermediate','iREC','REC','Stem','Stem NOTUM',
                       'TA1','TA2','Tuft','UPR')
rowData(seRNA)$gene <- rowData(seRNA)$name
rownames(seRNA) <- rowData(seRNA)$name
seRNA <- seRNA[p2geneDF_TSS$geneName]
rowData(seRNA)$linkage <- p2geneDF_TSS$linkage
rownames(seRNA) <- p2geneDF_TSS$linkage

#####

mat_sePeaks <- log2(assay(sePeaks)+1)
mat_seRNA <- log2(assay(seRNA)+1)
z_score_sePeaks <- t(scale(t(mat_sePeaks)))
z_score_seRNA <- t(scale(t(mat_seRNA)))
head(z_score_sePeaks)
saveRDS(z_score_sePeaks,file='Epithelial_peaks/Epithelial_activators_p2g_Peaks_zscore.rds')
saveRDS(z_score_seRNA,file='Epithelial_peaks/Epithelial_activators_p2g_RNA_zscore.rds')

library(circlize)
col_BuY = colorRamp2(c(-2, 2), c("blue", "yellow"))
col_BuRd = colorRamp2(c(-2,0,2), c("blue","white","red"))
library(ComplexHeatmap)
row_dends = as.dendrogram(hclust(dist(z_score_sePeaks)))
col_dends = as.dendrogram(hclust(dist(t(z_score_sePeaks))))

patternToGrep <- "-LGR5$|-ASCL2$|-SMOC2$|-PROM1$|-TFF3$|-MUC2$|-ATOH1$|-SLC26A3$|-KRT19$|-KRT20$|-FABP1$|-CA2$|-MMP7$|-KRT7$|-KRT17$|-PLAUR$|-ANXA1$|-LAMC2$|-EMP1$|-LAMA3$|-IFI6$|-IFI44$|-IFIT1$|-IFIT3$|-MKI67$|-TOP2A$|-PCNA$|-VEGFA$|-LDLR$|-CHGA$|-CHGB$|-NEUROD1$|-LRMP$|-FOSB$|-KRT19$|-KRT20$|-FOSB$|-FOS"
ha = rowAnnotation(z_score_sePeaks = anno_mark(at = grep(pattern = patternToGrep, rownames(z_score_sePeaks)), 
                                               labels = grep(patternToGrep, rownames(z_score_sePeaks), value = T),
                                               labels_gp = gpar(fontsize = 8)
                                               )
                   )
ha_test = rowAnnotation(z_score_sePeaks = anno_mark(at = grep(pattern = patternToGrep, rownames(z_score_sePeaks)), 
                                               labels = grep(patternToGrep, rownames(z_score_sePeaks), value = T),
                                               labels_gp = gpar(fontsize = 8)),
                   #annotation_name_gp= gpar(fontsize = 12),
                   gp = gpar(fontsize = 1))#which = 'row')
heatmapPeaks <- draw(Heatmap(z_score_sePeaks, cluster_rows = row_dends, cluster_columns = col_dends,
                             col = col_BuY, show_row_names = F, right_annotation = ha,
                             heatmap_legend_param = list(title = 'ATAC z-score'),
                             column_title = 'Accessibility'
                             )
                     )
heatmapRNA <- draw(Heatmap(z_score_seRNA, cluster_rows = row_dends, cluster_columns = col_dends,
                           col = col_BuRd, show_row_names = F, right_annotation = ha,
                           heatmap_legend_param = list(title = 'RNA z-score'),
                           column_title = 'Expression'))
### k means
# determine number of clusters to use
library(cluster)
hcluster = clusGap(z_score_sePeaks, FUN = kmeans, K.max = 9, B = 250)
dat <- data.table(hcluster$Tab)
dat[, k := .I]
p <- ggplot(dat, aes(k, gap)) + geom_line() + geom_point(size = 3) +
  geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), width = 0.25) +
  ggtitle("Clustering Results") +
  labs(x = "Number of Clusters", y = "Gap Statistic") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))
p

k = kmeans(z_score_sePeaks, centers = 6)
saveRDS(k, 'kmeans.rds')
cl = k$cluster
cl_df <- as.data.frame(cl)
cl_df$linkage <- row.names(cl_df)
library(stringr)
cl_df$chr <- str_split_fixed(cl_df$linkage, '_', n = 3)[,1]
cl_df$st <- str_split_fixed(cl_df$linkage, '_', n = 3)[,2]
cl_df$ed <- str_split_fixed(str_split_fixed(cl_df$linkage, '_', n = 3)[,3], '-', n=2)[,1]
dim(subset(cl_df, subset = cl == 1)[,c(3,4,5,2)])

for (i in unique(cl_df$cl)){
  write.table(subset(cl_df, subset = cl == i)[,c(3,4,5,2)],
              file = paste('Epithelial_peaks/kmeans_',i,'.bed', sep = ""),
              row.names = F, quote = F, sep='\t', col.names = F)
}

# classes from k-means are always put as the first column in `row_split`
heatmapPeaks2 <- draw(Heatmap(z_score_sePeaks, row_split = cl, cluster_columns = col_dends, cluster_row_slices = F,
                              col = col_BuY, show_row_names = F, right_annotation = ha,
                              heatmap_legend_param = list(title = 'ATAC z-score'),
                              column_title = 'Accessibility'))
heatmapRNA2 <- draw(Heatmap(z_score_seRNA, row_split = cl, cluster_columns = col_dends, cluster_row_slices = F,
                            col = col_BuRd, show_row_names = F, right_annotation = ha,
                            heatmap_legend_param = list(title = 'RNA z-score'),
                            column_title = 'Expression'))

split <- paste0("Cluster\n", cl)
split <- factor(paste0("Cluster\n", cl), levels=c("Cluster\n5","Cluster\n3","Cluster\n2","Cluster\n1","Cluster\n6","Cluster\n4"))
heatmapPeaks3 <- draw(Heatmap(z_score_sePeaks, split = split, cluster_columns = col_dends, cluster_row_slices = F,
                              col = col_BuY, show_row_names = F, right_annotation = ha,
                              heatmap_legend_param = list(title = 'ATAC z-score'),
                              column_title = 'Accessibility'))
heatmapRNA3 <- draw(Heatmap(z_score_seRNA, split = split, cluster_columns = col_dends, cluster_row_slices = F,
                            col = col_BuRd, show_row_names = F, right_annotation = ha,
                            heatmap_legend_param = list(title = 'RNA z-score'),
                            column_title = 'Expression'))

patternToGrep <- "-LGR5$|-ASCL2$|-SMOC2$|-PROM1$|-TFF3$|-MUC2$|-ATOH1$|-KRT19$|-KRT20$|-CA2$|-MMP7$|-KRT7$|-KRT17$|-PLAUR$|-LAMC2$|-EMP1$|-LAMA3$|-IFI6$|-IFI44$|-IFIT1$|-IFIT3$|-MKI67$|-TOP2A$|-PCNA$|-LDLR$|-CHGA$|-CHGB$|-NEUROD1$|-LRMP$|-FOSB$|-KRT19$|-KRT20$|-FOSB$|-FOS$|chr7_107843123_107843623-SLC26A3|chr7_107860920|chr6_44005151|chr6_43787376"
patternToGrep <- "chr18_23872608_23873108???|chr7_107843123_107843623???|chr7_107860920_107861420???|-LGR5$|-ASCL2$|-EMP1$|-PLAUR$|-FOS$|-LAMC2$|-KRT19$|chr7_107843123_107843623-SLC26A3|chr7_107860920|chr6_44005151|chr6_43787376|chr6_168192138_168192638???"
ha2 = rowAnnotation(z_score_sePeaks = anno_mark(at = grep(pattern = patternToGrep, rownames(z_score_sePeaks)), 
                                               labels = grep(patternToGrep, rownames(z_score_sePeaks), value = T),
                                               labels_gp = gpar(fontsize = 8)
                                               )
                    )

heatmapPeaks4 <- draw(Heatmap(z_score_sePeaks, split = split, cluster_columns = col_dends, cluster_row_slices = F,
                              col = col_BuY, show_row_names = F, right_annotation = ha2, show_row_dend = FALSE,
                              heatmap_legend_param = list(title = 'ATAC z-score'),
                              column_title = 'Accessibility'))
heatmapRNA4 <- draw(Heatmap(z_score_seRNA, split = split, cluster_columns = col_dends, row_order = unlist(row_order(heatmapPeaks3, name = NULL)),
                            col = col_BuRd, show_row_names = F, right_annotation = ha2,  
                            heatmap_legend_param = list(title = 'RNA z-score'), 
                            column_title = 'Expression'))

cl2 <- as.data.frame(cl)
cl2$peaks <- labels(cl)
cl2$genes <- str_split_fixed(cl2$peaks, '-', 2)[,2]
write.csv(cl2, file = 'Epithelial_peaks/Epithelial_kmeans_final.csv')

pdf('Epithelial_peaks/Epithelial_activators_peak2gene_heatmap_final.pdf')
heatmapPeaks
heatmapRNA
heatmapPeaks2
heatmapRNA2
heatmapPeaks3
heatmapRNA3
heatmapPeaks4
heatmapRNA4
dev.off()

