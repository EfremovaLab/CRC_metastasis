import os
import pandas as pd
import scanpy as sc
import dill
import pickle

sc.settings.set_figure_params(dpi=60, frameon=True, figsize=(5, 5), facecolor='white')
pd.set_option('display.max_columns', 10)

import scenicplus
print (scenicplus.__version__)

#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

work_dir = '/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/scenicplus_unionPeakset'
tmp_dir = '/data/scratch/union/'

n_cpu = int(os.environ["NSLOTS"])

### scATAC-seq preprocessing using pysistopic

import pycisTopic
#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scATAC')):
    os.makedirs(os.path.join(work_dir, 'scATAC'))
if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))

#path to fragments files
fragpath = '/data/BCI-CRC/SO/data/CRC_multiome/cellranger-ARC/count/'
fragpath2 = 'outs/atac_fragments.tsv.gz'
fragments_dict = {
    'CRC01_LM': os.path.join(fragpath, 'CRC01_LM', fragpath2),
    'CRC02_LM': os.path.join(fragpath, 'CRC02_LM', fragpath2),
    'CRC03_LM': os.path.join(fragpath, 'CRC03_LM', fragpath2),
    'CRC04_LM': os.path.join(fragpath, 'CRC04_LM', fragpath2),
    'CRC05_LM': os.path.join(fragpath, 'CRC05_LM', fragpath2),
    'CRC06_LM': os.path.join(fragpath, 'CRC06_LM', fragpath2),
    'CRC07_LM': os.path.join(fragpath, 'CRC07_LM', fragpath2),
    'CRC08_LM': os.path.join(fragpath, 'CRC08_LM', fragpath2),
    'CRC09_LM': os.path.join(fragpath, 'CRC09_LM', fragpath2),
    'CRC10_LM': os.path.join(fragpath, 'CRC10_LM', fragpath2),
    'CRC11_LM': os.path.join(fragpath, 'CRC11_LM', fragpath2),
    'CRC12_LM': os.path.join(fragpath, 'CRC12_LM', fragpath2),
    'CRC13_LM': os.path.join(fragpath, 'CRC13_LM', fragpath2),
    'CRC14_LM': os.path.join(fragpath, 'CRC14_LM', fragpath2),
    'CRC15_LM': os.path.join(fragpath, 'CRC15_LM', fragpath2)
}

# load the cell type annotation we generated in the scRNA-seq analysis above
adata = sc.read('/data/BCI-CRC/SO/data/CRC_multiome/scanpy/CRCLM_finalAnalysis/Epithelial_scvi_annotations.h5ad')
adata.obs['Cell_subtype'] = adata.obs['Cell_subtype'].replace({'Stem (NOTUM) high':'Stem_NOTUM'})
adata.write(os.path.join(work_dir, 'scRNA/adata.h5ad'))
cols = ['nCount_RNA','nFeature_RNA','nCount_ATAC','nFeature_ATAC']
cell_data = adata.obs[cols]
cell_data['sample_id'] = adata.obs['Sample'].astype(str)
cell_data['celltype'] = adata.obs['Cell_subtype'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
del(adata)
cell_data.index = [i.split('#')[1] for i in cell_data.index]
cell_data['barcode'] = cell_data.index
#save result
cell_data.to_csv(os.path.join(work_dir, 'scATAC/cell_data.tsv'), sep = '\t')
print(cell_data.head(2))

### Creat cisTopic object and topic modelling using union peakset from ArchR

path_to_regions=dict()
for key in fragments_dict:
    path_to_regions[key] = '/data/BCI-CRC/SO/data/CRC_multiome/ArchR_final_analysis/Epithelial_peaks/Epithelial_unionPeakset.bed'

from pycisTopic.cistopic_class import *

path_to_blacklist = '/data/BCI-CRC/SO/data/CRC_multiome/ArchR_WNN_analysis/blacklist.bed'

from pycisTopic.cistopic_class import *
cistopic_obj_list=[create_cistopic_object_from_fragments(path_to_fragments=fragments_dict[key],
                                                         path_to_regions=path_to_regions[key],
                                                         path_to_blacklist=path_to_blacklist,
                                                         metrics=cell_data[cell_data['sample_id']==key],
                                                         valid_bc=list(cell_data[cell_data['sample_id']==key].index),
                                                         n_cpu=n_cpu,
                                                         project=key) for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)
print(cistopic_obj)

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj_backup.pkl'), 'wb'))

cell_data = pd.read_csv(os.path.join(work_dir, 'scATAC/cell_data.tsv'), sep = '\t', index_col=0)
#cistopic makes use of the sample_id to match the correct cell barcodes to the metadata, let's add the sample_id as a suffix to the cell barcodes
cell_data['barcode'] = cell_data['barcode'] +'___'+ cell_data['sample_id']
print(cell_data['barcode'][0:5])
cell_data = cell_data.set_index('barcode')
cistopic_obj.add_cell_data(cell_data[['sample_id']])

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

# Run topic modelling:
## 1. finds set of co-accessible regions (topics), which are used as candidate enhancers with DARs
## 2. impute dropouts
# Can be very computationally intense - may want to submit job to cluster

from pycisTopic.cistopic_class import *
models=run_cgs_models(cistopic_obj,
                    n_topics=[12,16,24,32],
                    n_cpu=n_cpu,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + '/ray_spill'))

if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/CRCLM_models_500_iter_LDA.pkl'), 'wb'))


###
# Select model and identify candidate enhancers
###

model_number = 24

cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
print(cistopic_obj)

models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/CRCLM_models_500_iter_LDA.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
from pycisTopic.lda_models import *
model = evaluate_models(models, save='Figures/models.pdf',
                       select_model=model_number,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)

# Add selected models as metrics stabilise
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)

# Infer candidate enhancer regions
# 1. binarise region-topic probabilities
# 2. calculate DARs per cell type
# Use above in pycistarget - look for motif enrichment

# binarise with otsu method:
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

# calculate DARs per cell state
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions, split_pattern = '-')

if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))

pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))

###
### pycistarget for motif enrichment in topics and DARs
###

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

db_fpath = "create_cistarget_databases/"
motif_annot_fpath = "/data/BCI-CRC/SO/genomes/scenicplus"

rankings_db = os.path.join(db_fpath, 'unionPeakset.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'unionPeakset.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

### Run pyCistarget

if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

### may get ray error if parallelize
from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = n_cpu,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    )

###
# Scenicplus
###

adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/adata.h5ad'))
adata.obs_names = [i.split('#')[1]+'___'+i.split('#')[0] for i in adata.obs_names]
adata = adata.raw.to_adata()
sc.pp.filter_genes(adata, min_counts=None, min_cells=60, inplace=True, copy=False)
print(adata.X.shape)
cistopic_obj = dill.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

## create scenicplus object
# Cell metadata comming from the cistopic_obj will be prefixed with the string ACC_ and
# metadata comming from the adata object will be prefixed with the string GEX_.

from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
print(scplus_obj)

dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)

### Add adjacenices (generated using pyscenic) to scenicplus object
#arboreto_with_multiprocessing.py \
#    pyscenic.loom \
#    /data/BCI-CRC/SO/genomes/scenicplus/allTFs_hg38.txt \
#    --method grnboost2 \
#    --output adj.tsv \
#    --num_workers $NSLOTS \
#    --seed 777

scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

from scenicplus.TF_to_gene import load_TF2G_adj_from_file
load_TF2G_adj_from_file(scplus_obj,
                        f_adj = 'pyscenic/adj.tsv',
                        inplace = True,
                        key= 'TF2G_adj')

# set biomart host
biomart_host = "http://sep2019.archive.ensembl.org/"

### Run scenicplus

from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_Cell_subtype'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '/data/BCI-CRC/SO/genomes/scenicplus/allTFs_hg38.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = False,
        export_to_loom_file = True,
        export_to_UCSC_file = False,
        path_bedToBigBed = '/data/BCI-CRC/SO/packages',
        n_cpu = n_cpu,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill2'))

###        object_store_memory=500000000000) #500 gb
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj_failure.pkl'), 'wb'), protocol=-1)
    raise(e)

dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj_finished.pkl'), 'wb'), protocol=-1)
