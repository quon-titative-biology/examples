# data prep
## load the single cell data from Tasic et al. to make pseudo bulk samples 

(please contact hrhu@ucdavis.edu to request processed data)

```
tasic
# An object of class Seurat
# 28163 features across 13586 samples within 1 assay
# Active assay: RNA (28163 features, 2000 variable features)
table(tasic@meta.data$subclass)
#      Astro         CR       Endo    L2/3 IT         L4      L5 IT      L5 PT
#        361          7         76        969       1350        812        536
#      L6 CT      L6 IT        L6b      Lamp5 Macrophage      Meis2         NP
#        953       1836        329       1070         51         44        338
#      Oligo       Peri      Pvalb   Serpinf1        SMC       Sncg        Sst
#         91         28       1211         26         55        125       1567
#        Vip       VLMC
#       1690         61

scRNA <- tasic[,tasic@meta.data$subclass%in%c("L2/3 IT", "L4", "L5 IT", "L6 IT", "Lamp5", "Pvalb", "Sst", "Vip")]
table(scRNA@meta.data$subclass)
# L2/3 IT      L4   L5 IT   L6 IT   Lamp5   Pvalb     Sst     Vip
#     969    1350     812    1836    1070    1211    1567    1690
```

## generate pseudo bulk samples
```
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
from random import randint
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg
import scipy.io as sio

input_path = '/home/qlab/scProjection/pseudobulk/'
output_path2 = '/home/qlab/scProjection/pseudobulk/ct2/'
output_path4 = '/home/qlab/scProjection/pseudobulk/ct4/'
output_path6 = '/home/qlab/scProjection/pseudobulk/ct6/'
output_path8 = '/home/qlab/scProjection/pseudobulk/ct8/'

ct2 = ['L23', 'Lamp5']
ct4 = ['L23','L4', 'Lamp5','Pvalb']
ct6 = ['L23','L4','L5', 'Lamp5','Pvalb','Sst']
ct8 = ['L23','L4','L5','L6', 'Lamp5','Pvalb','Sst','Vip']

dir = list([output_path2, output_path4, output_path6, output_path8])
celltypes = list([ct2, ct4, ct6, ct8])

########################
#### load scRNA dat ####
########################
adata_sc = sc.read_h5ad(input_path + "allen_all8ct.h5ad")
adata_sc
# AnnData object with n_obs × n_vars = 10505 × 25478
#     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample_name', 'sample_id', 'donor', 'sex', 'age_days', 'eye_condition', 'genotype', 'driver_lines', 'reporter_lines', 'brain_hemisphere', 'brain_region', 'brain_subregion', 'injection_label_direction', 'injection_primary', 'injection_secondary', 'injection_tract', 'injection_material', 'injection_exclusion_criterion', 'facs_date', 'facs_container', 'facs_sort_criteria', 'rna_amplification_set', 'rna_amplification_pcr_cycles', 'library_prep_set', 'library_prep_avg_size_bp', 'seq_name', 'seq_tube', 'seq_batch', 'total_reads', 'percent_exon_reads', 'percent_intron_reads', 'percent_intergenic_reads', 'percent_rrna_reads', 'percent_mt_exon_reads', 'percent_reads_unique', 'percent_synth_reads', 'percent_ecoli_reads', 'percent_aligned_reads_total', 'complexity_cg', 'genes_detected_cpm_criterion', 'genes_detected_fpkm_criterion', 'tdt_cpm', 'gfp_cpm', 'class', 'subclass', 'cluster', 'confusion_score', 'cluster_correlation', 'core_intermediate_call'
#     var: 'features'
meta = adata_sc.obs.copy()
pd.unique(meta.subclass)
# ['Pvalb', 'L4', 'Vip', 'L23', 'Lamp5', 'Sst', 'L5', 'L6'],
meta.groupby("subclass").size()
# subclass
# L23       969
# L4       1350
# L5        812
# L6       1836
# Lamp5    1070
# Pvalb    1211
# Sst      1567
# Vip      1690

for i in range(len(dir)):
    meta_subset = meta[meta.subclass.isin(celltypes[i])]
    meta_subset.to_csv(dir[i] + "meta.csv")
    add_list = list()
    for j in range(5000):
        meta_cells = meta_subset.groupby("subclass").sample(n=1, random_state=j).index
        sample_mixture = list(meta_cells)
        print(j, " pair with ",  meta_subset.loc[meta_cells].subclass)
        add_list.append(sample_mixture)
    index_list = list()
    anno_list = list()
    for k in range(len(add_list)):
        print("make bulk cell", k)
        anno_list.append('&'.join((str(n) for n in meta_subset.loc[add_list[k]].subclass.values)))
        index_list.append('&'.join((str(n) for n in add_list[k])))
    samples = pd.DataFrame([index_list, anno_list], index=["pseudobulk_index","pseudobulk_celltype"]).T
    samples.to_csv(dir[i] + "samples_pseudobulk_allen.csv")
```

## using 8-cell-type mixtures as an exmaple
```
# cd /home/qlab/scProjection/pseudobulk/
import pathlib
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from sklearn.decomposition import PCA
from anndata import AnnData
import scanpy as sc
import squidpy as sq
import tangram as tg
import umap
import pickle 
import scProjection # or deconv for old version

celltype_mixture = "ct8/"
############################################
path = "/home/qlab/scProjection/pseudobulk/"

adata_sc = sc.read_h5ad(path + celltype_mixture + "allen.h5ad")
adata_sc

adata_mer = sc.read_h5ad(path + celltype_mixture + "allen_bulk.h5ad")
adata_mer

hvgs = adata_sc.var[adata_sc.var['vst.variable'] == 1].index

component_data  = np.array(adata_sc.raw.X.todense())
component_label = adata_sc.obs.subclass
mixture_data  = np.array(adata_mer.raw.X.todense())

deconvModel = deconv.deconvModel(component_data   = component_data,
                                 component_label  = component_label,
                                 mixture_data     = mixture_data) 

marker_gene_mask = list(np.ones(hvgs.size,))
deconvModel.deconvolve(marker_gene_mask         = marker_gene_mask,
                       ## Hyper params
                       max_steps_component      = 2500,
                       max_steps_proportion     = 2500, 
                       max_steps_combined       = 100,
                       component_learning_rate  = 1e-3,
                       proportion_learning_rate = 1e-2,
                       combined_learning_rate   = 1e-4,
                       ## Early stopping
                       early_stopping           = 'True',
                       early_min_step           = 500,
                       max_patience             = 50,
                       ## Output
                       log_step                 = 100,
                       print_step               = 100,
                       ## Seed
                       seed                     = 12,
                       ## VAE params
                       KL_weight                = 1.0,
                       heldout_percent          = 0.0,
                       num_latent_dims          = 200,     
                       num_layers               = 3,
                       hidden_unit_2power       = 9,
                       decoder_var              = "none",
                       ## Batch size
                       batch_size_component     = 150,
                       batch_size_mixture       = 1000,
                       ## Training method
                       training_method          = 'train',
                       dropout_rate             = 0.0,     
                       tpm                      = "True",  
                       deconvolution            = 'True',
                       batch_correction         = 'False',
                       decay_lr                 = 'False', 
                       batch_norm_layers        = 'False', 
                       corr_check               = 'True',
                       log_results              = 'True',
                       log_samples              = 'False',
                       save_to_disk             = 'False',
                       save_to_object           = 'True',
                       cuda_device              = 0)




deconvModel.deconvResults.deconv_data.keys()
# dict_keys(['component', 'purified'])
```


## save and load models
```
file_pickle = open(path + celltype_mixture + 'deconv_overfitting.obj', 'wb') 
pickle.dump(deconvModel, file_pickle)
file_pickle.close()
```
```
fileObj = open(path + celltype_mixture + 'deconv_overfitting.obj', 'rb')
exampleObj = pickle.load(fileObj)
fileObj.close()
print(exampleObj)
deconvModel = exampleObj
```
