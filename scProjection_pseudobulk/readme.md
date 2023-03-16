# data prep
## load the single cell data from Tasic et al. to make pseudo bulk samples 

(please contact hrhu@ucdavis.edu to request processed data)

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

input_path = '/group/gquongrp/workspaces/hongruhu/scProjection/sc_atlas/pseudobulk_data/'
output_path2 = '2_two/'
output_path4 = '4_four/'
output_path6 = '6_six/'
output_path8 = '8_eight/'

ct2 = ['L23', 'Lamp5']
ct4 = ['L23','L5', 'Lamp5','Sst']
ct6 = ['L23','L5','L6', 'Lamp5','Sst','Vip']
ct8 = ['L23','L5','L6','L45', 'Lamp5','Sst','Vip','Pvalb']

dir = list([input_path + output_path2, 
            input_path + output_path4, 
            input_path + output_path6, 
            input_path + output_path8])
celltypes = list([ct2, ct4, ct6, ct8])
########################
#### load scRNA dat ####
########################
adata_sc = sc.read_h5ad(input_path + "mop.h5ad")
adata_sc
# AnnData object with n_obs × n_vars = 4502 × 28140

meta = adata_sc.obs.copy()
pd.unique(meta.subclass_label)
# array(['Sst', 'L5', 'Vip', 'Lamp5', 'L45', 'Pvalb', 'L23', 'L6'],
meta.groupby("subclass_label").size()
# subclass_label
# L23       555
# L45      1229
# L5        465
# L6        345
# Lamp5     375
# Pvalb     528
# Sst       410
# Vip       595


for i in range(len(dir)):
    meta_subset = meta[meta.subclass_label.isin(celltypes[i])]
    meta_subset.to_csv(dir[i] + "meta.csv")
    add_list = list()
    for j in range(5000):
        meta_cells = meta_subset.groupby("subclass_label").sample(n=1, random_state=j).index
        sample_mixture = list(meta_cells)
        print(j, " pair with ",  meta_subset.loc[meta_cells].subclass_label)
        add_list.append(sample_mixture)
    index_list = list()
    anno_list = list()
    for k in range(len(add_list)):
        print("make bulk cell", k)
        anno_list.append('&'.join((str(n) for n in meta_subset.loc[add_list[k]].subclass_label.values)))
        index_list.append('&'.join((str(n) for n in add_list[k])))
    samples = pd.DataFrame([index_list, anno_list], index=["pseudobulk_index","pseudobulk_celltype"]).T
    samples.to_csv(dir[i] + "samples_pseudobulk_allen.csv")
```

## using 8-cell-type mixtures as an exmaple
```
import pathlib
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import MaxAbsScaler
from anndata import AnnData
import scanpy as sc
import squidpy as sq
import tangram as tg
import pickle 
import scProjection
import deconv

celltype_mixture = "8_eight/"
method = "scProjection/"
################################################################################
# data loading
path = "/group/gquongrp/workspaces/hongruhu/scProjection/sc_atlas/pseudobulk_data/"

adata_sc = sc.read_h5ad(path + celltype_mixture + "scRNA.h5ad")
adata_sc
adata_mer = sc.read_h5ad(path + celltype_mixture + "bulk.h5ad")
adata_mer
adata_ref = sc.read_h5ad(path + celltype_mixture + "scRNA_ref.h5ad")
adata_ref
################################################################################
# data set up
component_label = adata_ref.obs.subclass_label
component_data  = np.array(adata_ref.X.todense())
mixture_data  = np.array(adata_mer.X.todense()/8)
hvgs = adata_ref.var.index
marker_gene_mask = list(np.ones(hvgs.size,))
# model set up
deconvModel = deconv.deconvModel(component_data   = component_data,
                                 component_label  = component_label,
                                 mixture_data     = mixture_data) 
deconvModel.deconvolve(marker_gene_mask         = marker_gene_mask,
                       ## Hyper params
                       max_steps_component      = 2500,
                       max_steps_proportion     = 100, 
                       max_steps_combined       = 100,
                       component_learning_rate  = 1e-3,
                       proportion_learning_rate = 1e-2,
                       combined_learning_rate   = 1e-4,
                       ## Early stopping
                       early_stopping           = 'False', #
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
                       num_latent_dims          = 200,     #
                       num_layers               = 3,
                       hidden_unit_2power       = 9,
                       decoder_var              = "none",
                       ## Batch size
                       batch_size_component     = 150,
                       batch_size_mixture       = 1000,
                       ## Training method
                       training_method          = 'train',
                       dropout_rate             = 0.0,     #
                       tpm                      = "True",  #
                       deconvolution            = 'True',
                       batch_correction         = 'False',
                       decay_lr                 = 'False', #
                       batch_norm_layers        = 'False', #
                       corr_check               = 'True',
                       log_results              = 'True',
                       log_samples              = 'False',
                       save_to_disk             = 'False',
                       save_to_object           = 'True',
                       cuda_device              = 0)


deconvModel.deconvResults.deconv_data.keys()
# dict_keys(['component', 'purified'])
################################################################################
# projection evaluation
meta = adata_mer.obs[['pseudobulk_index','pseudobulk_celltype']] 
expr = pd.DataFrame(adata_sc.X.todense(), index=list(adata_sc.obs.index))
meta.pseudobulk_celltype.str.split("&",expand=True)
ct_name = list(meta.pseudobulk_celltype.str.split("&",expand=True).iloc[0,:])
for i in range(len(ct_name)):
    idx = meta.pseudobulk_index.str.split("&",expand=True).iloc[:,i]
    component_ct = expr.loc[idx].to_numpy()
    purified_ct = deconvModel.deconvResults.deconv_data['purified']['train'][ct_name[i]]
    purified_ct_ = purified_ct
    measured = component_ct
    deconv = purified_ct_
    cell_wise = list()
    for j in range(deconv.shape[0]):
        cell_wise.append(np.corrcoef(deconv[j,:],measured[j,:])[0,1])
    print(ct_name[i]+":", np.mean(cell_wise),np.median(cell_wise))
    pd.DataFrame(deconv).to_csv(path + celltype_mixture+ method + "scProjection_deconv_res_" + ct_name[i] +".csv")
    pd.DataFrame(measured).to_csv(path + celltype_mixture+ method + "scProjection_measured_res_" + ct_name[i] +".csv")
    idx.to_csv(path + celltype_mixture + method + "sample_id_" + ct_name[i] +".csv")


```
the training will take about less than 30 min with one GPU

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
