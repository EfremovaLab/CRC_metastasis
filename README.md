# Code for the analyses in Ogden, Metic et al. "Phenotypic heterogeneity and plasticity in colorectal cancer metastasis"

Phenotypic heterogeneity and plasticity in colorectal cancer (CRC) is increasingly recognised as a driver of tumour progression, metastasis and therapy resistance. However, the regulatory factors and the extrinsic signals from the microenvironment driving phenotypic heterogeneity in primary and metastatic CRC remain unknown. 

Using a combination of single-cell multiomics and spatial transcriptomics data from primary and metastatic CRC patients, we reveal cancer cell states with regenerative and inflammatory phenotypes that closely resemble metastasis-initiating cells in mouse models. We identify an intermediate population with a hybrid regenerative and stem phenotype, indicating phenotypic transitions between stem and pro-metastatic cells. We show that the regenerative cells have elevated activity of transcription factors AP-1 and NF-κB and validate the role of AP-1 using patient-derived organoids. 

Our spatial analyses show localisation of the regenerative inflammatory states in an immunosuppressive niche both at the invasive edge in primary CRC and in liver metastasis. We uncover ligand-receptor interactions driven by cancer-associated fibroblasts (CAFs), immunosuppressive and inflammatory macrophages and exhausted CD8 T cells predicted to activate the regenerative and inflammatory phenotype in cancer cells. 

Together, our findings reveal regulatory and signalling factors that mediate distinct cancer cell states and can serve as potential targets to restrict transition into invasive states to impair metastasis.

## Highlights

* Profiled liver metastatic samples from colorectal cancer patients using single-cell multiomics and integrated them with public primary and metastatic data.
* Identified regenerative and inflammatory cancer cell states and showed that the transcription factors AP-1 and NF-κB are their key regulators.
* The regenerative cells are localised at the invasive edge and reside in an immunosuppressive niche.

## Overview
This project includes workflows covering data preprocessing, quality control, scRNA-seq/snRNA-seq integration (via scVI), dimensionality reduction, clustering, spatial deconvolution (cell2location), spatial niche inference and cell-cell communication (CellPhoneDB, NicheNet). Code for our [manuscript](https://doi.org/10.21203/rs.3.rs-3846377/v1) Ogden, Metic et al. "Phenotypic heterogeneity and plasticity in colorectal cancer metastasis".

<figure>
<img src="workflow_overall.png"
     alt=""
     width="100%" />
<figcaption>
    Overview of the computational workflow. Schematic describing the computational analysis done in this project. To examine the cellular heterogeneity and composition of malignant cell states and the tumour microenvironment in primary CRC and see how this landscape changes in liver metastasis, we integrated published pCRC scRNA-seq data and generated single nucleus Multiome RNA+ATAC data of liver metastasis. To predict transcriptional regulators of cancer cell states, we investigated the accessible chromatin landscape of cancer cells in our mCRC dataset. To infer the spatial arrangement of cancer cell states and TME cells defined by scRNA-seq data, we coupled scRNA-seq data to spatial transcriptomics data (10X Visium), revealing cellular niches and spatial interactions associated with specific spatial niches.
  </figcaption>
</figure>

<br>  <!-- Force a line break -->

## Data Requirements

The pipeline relies on multiple datasets:

- Public scRNA-seq data for primary CRC (pCRC) and metastatic CRC (mCRC)
- In-house 10X Multiome data (snRNA-seq + snATAC-seq) generated by colleague Dr Sam Ogden
- 10X Visium spatial transcriptomics data (in-house and/or publicly available)
- RNA-seq data from TCGA

## Dependencies

`.yml` files are provided in `vens/` folder. Key libraries include (but are not limited to):
- Python 3.8+
- Scanpy
- scVI-tools
- Scrublet
- CellPhoneDB
- NicheNet (R package)
- cell2location
- SpatialDE2
- Squidpy
- ArchR (for snATAC-seq in R)
- SCENIC+ (Python + R)
- GSEApy
- Other standard scientific libraries (NumPy, SciPy, pandas, matplotlib, seaborn, statsmodels etc.)
- vistools which can be found here https://github.com/Nasrine26/vistools

Make sure that you have GPU support if running scVI or large scale cell2location analyses. All pipelines were run on high performance clustering (HPC) clusters to efficiently handle large datasets (for scVI and cell2location with GPU support).

## Repository Structure

- `ccc/`: Directory for cell-cell communication analysis for lmCRC
- `scrna-seq/`: Directory for single cell RNA-seq analysis, single-nuclei RNA-seq analysis, RNA-seq analayis for pCRC and lmCRC
- `spatial/`: Directory for spatial analysis for pCRC and lmCRC
- `venvs/`: yml files for envs used

## Analysis Workflow Overview

**Raw Data Ingestion**
- Load `.h5ad` or `.mtx` files for scRNA-seq/spatial transcripotmics
- Load TCGA RNA-seq data, normalise and preprocess
  
**Quality Control**
- Filter low-quality cells/nuclei (mito < 20%, genes > 300, etc.)
- Identify doublets
  
**Integration & Dimensionality Reduction**
- scVI: Variational autoencoder for batch correction and latent embedding. I align multiple scRNA-seq/snRNA-seq datasets while adjusting for patient and technical covariates. 
  
**Visualisation & Clustering**
- KNN graph (k varies) → UMAP → Leiden (r varies)
- Identify major cell types, iterative sub-clustering for TME compartments. I leverage Leiden clustering, differential expression with the Wilcoxon rank-sum test, and canonical marker genes to annotate major cell types and sub-clusters.

<figure>
<img src="workflow_scrnaseq.png"
     alt=""
     width="70%" />
<figcaption>
    Overview of scRNA-seq analysis workflow. Each count matrix was preprocessed individually. Doublet detection method Scrublet was applied. Individual count matrices were filtered to remove low quality cells. The combined count matrix was then analysed via the Scanpy pipeline, reducing the dimensionality of the data and batch-correcting the data with scVI. Data are visualised with UMAP in a two-dimensional embedding. Unsupervised clus- tering with Leiden community detection method was performed on KNN graph to identify major cell populations and subsequently, the same integration and clustering analysis was ap- plied iteratively to the cells of each major cell type separately to identify and annotate TME cell states.
  </figcaption>
</figure>

<br>  <!-- Force a line break -->
  
**Spatial Transcriptomics**
- Map scRNA-seq–derived cell states with cell2location
- Segment tissue with SpatialDE2 and cell type mRNA count estimates to identify spatial niches

<figure>
<img src="workflow_spatial.png"
     alt=""
     width="70%" />
<figcaption>
    Overview of spatial niche identification workflow. Schematic describing the analysis of the spatial organisation of pCRC and lmCRC samples and identification of local niches surrounding cancer cell states. To elucidate the influence of TME on cancer cell states in pCRC and lmCRC samples, we require better understanding of the tumour architecture. Cell2location was applied to assign cell types to spots by integrating scRNA-seq and ST transcriptomes, yielding cell type abundance estimates per spot. Subsequently, we deciphered local niches surrounding the cancer cell states and identified spatially co-occurring cell types by spatially clustering on the cell type abundance estimates per spot across samples.
  </figcaption>
</figure>

<br>  <!-- Force a line break -->

**Cell-Cell Interaction Analysis**
- CellPhoneDB, NicheNet for ligand-receptor and downstream TF network
- TF regulons, co-accessibility (snATAC + snRNA) were estimated with SCENIC+ by colleague Dr Sam Ogden

<figure>
<img src="workflow_ccc.png"
     alt=""
     width="80%" />
<figcaption>
    Overview of cell communication analysis workflow. Schematic describing the cell communication analysis in lmCRC samples after identification of a spatial niche surrounding iRECs. Combining CellPhoneDB and NicheNet methods, the analysis identified potential ligands from TME subpopulations in the iREC spatial niche predicted to induce expression of AP-1 and NF-κB regulons.
  </figcaption>
</figure>

<br>  <!-- Force a line break -->

**Statistical Tests & Survival**
- Differential expression (Wilcoxon)
- Gene/Pathway enrichment analysis with GSEApy, and Molecular Signatures Database (MSigDB) hallmark gene set (2020) and KEGG pathway gene sets (2021 Human)
- Spearman/Pearson correlation of cell state proportions
- Kruskal-Wallis for cell type enrichment analysis across histopathological regions and spatial niches
- Spatial ligand enrichment analysis
- Kaplan-Meier analysis for clinical outcome and niche signatures
- Clinical staging and cancer cell state signatures association.
