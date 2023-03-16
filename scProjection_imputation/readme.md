## Tutorial for gene imputation using osmFISH dataset


#### data loading: osmFISH data and reference singel cell RNA-seq data
```
import deconv
import pickle 
import scanpy
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from scipy.stats import spearmanr
from scvi.data import smfish, cortex
from scvi.external import GIMVI
path = "/group/gquongrp/workspaces/hongruhu/scProjection/sc_atlas/imputation/scProjection/"

spatial_data = smfish()
seq_data = cortex()

spatial_data
# AnnData object with n_obs × n_vars = 4530 × 33
#     obs: 'x_coord', 'y_coord', 'labels', 'str_labels', 'batch'
#     uns: 'cell_types'

seq_data
# AnnData object with n_obs × n_vars = 3005 × 19972
#     obs: 'labels', 'precise_labels', 'cell_type'
```

#### data processing
```
all_genes = list(spatial_data.var.index.sort_values())
spatial_data = spatial_data[:, all_genes].copy()
#seq_data = seq_data[:, all_genes].copy()
spatial_data
seq_data

scanpy.pp.filter_cells(spatial_data, min_counts= 1)
scanpy.pp.filter_cells(seq_data, min_counts = 1)
print(spatial_data, '\n', seq_data)

scanpy.pp.normalize_total(spatial_data, target_sum=1e4)
scanpy.pp.normalize_total(seq_data, target_sum=1e4)
spatial_data.X.sum(1)
seq_data.X.sum(1)

scanpy.pp.log1p(spatial_data)
scanpy.pp.log1p(seq_data)
spatial_data.X.sum(1)
seq_data.X.sum(1)

scanpy.pp.scale(spatial_data, max_value=10)
scanpy.pp.scale(seq_data, max_value=10)
spatial_data.X.sum(1)
seq_data.X.sum(1)

spatial_data.X
seq_data.X

spatial_data.X.shape
seq_data.X.shape
```

#### leave-one-gene-out per training
```
seq_data = seq_data[:, spatial_data.var_names].copy()
seq_data
seq_gene_names = seq_data.var_names
n_genes = seq_data.n_vars
n_train_genes = n_genes - 1 #int(n_genes*train_size)

# number of genes of the reference single cell dataset has been reduced to 33 which matches the spatial data
seq_data
# AnnData object with n_obs × n_vars = 3005 × 33
#     obs: 'labels', 'precise_labels', 'cell_type'
seq_gene_names
# Index(['Acta2', 'Aldoc', 'Anln', 'Apln', 'Bmp4', 'Cnr1', 'Cpne5', 'Crh',
#        'Crhbp', 'Ctps', 'Flt1', 'Foxj1', 'Gad2', 'Gfap', 'Hexb', 'Itpr2',
#        'Kcnip2', 'Lamp5', 'Mfge8', 'Mrc1', 'Pdgfra', 'Plp1', 'Pthlh', 'Rorb',
#        'Serpinf1', 'Slc32a1', 'Sox10', 'Syt6', 'Tbr1', 'Tmem2', 'Ttr', 'Vip',
#        'Vtn'],

print(n_genes, n_train_genes) # 33 32
seq_data.obs.cell_type[seq_data.obs.cell_type == "pyramidal SS"] = "pyramidal_SS"
seq_data.obs.cell_type[seq_data.obs.cell_type == "pyramidal CA1"] = "pyramidal_CA1"
seq_data.obs.cell_type[seq_data.obs.cell_type == "endothelial-mural"] = "endothelial_mural"
```

#### training the models (33 models, in each model, using the 32 gene expression patterns to impute the held-out one)
```
corr_list = []
df_original = pd.DataFrame(index=range(spatial_data.shape[0]), columns=all_genes)
df_imputed = pd.DataFrame(index=range(spatial_data.shape[0]), columns=all_genes)
#randomly select training_genes
for rand_test_gene_idx in range(n_genes):
    print("impute " + str(rand_test_gene_idx) + "th gene")
    rand_train_gene_idx =  sorted(set(range(n_genes)) - set([rand_test_gene_idx]))
    print("using " + str(rand_train_gene_idx) + "th gene")
    rand_train_genes = seq_gene_names[rand_train_gene_idx]
    rand_test_genes = seq_gene_names[rand_test_gene_idx]
    print(rand_train_genes, rand_test_genes)
    #spatial_data_partial has a subset of the genes to train on 
    spatial_data_partial = spatial_data[:,rand_train_genes].copy()
    component_data  = np.array(seq_data.X)
    component_label = (seq_data.obs.cell_type)
    mixture_data  = np.array(spatial_data_partial.X)
    deconvModel = deconv.deconvModel(component_data   = component_data,
                                     component_label  = component_label,
                                     mixture_data     = mixture_data) 
    marker_gene_mask = list(np.zeros(n_genes,))
    for i in rand_train_gene_idx:
        marker_gene_mask[i] = 1
    marker_gene_mask = [int(i) for i in marker_gene_mask]
    print(np.sum(marker_gene_mask))
    final = 2500 #2500
    deconvModel.deconvolve(marker_gene_mask         = marker_gene_mask,
                           ## Hyper params
                           max_steps_component      = 240, #1000,
                           max_steps_proportion     = final, 
                           max_steps_combined       = 100,
                           component_learning_rate  = 1e-3,
                           proportion_learning_rate = 1e-2,
                           combined_learning_rate   = 1e-4,
                           ## Early stopping
                           early_stopping           = 'True', #
                           early_min_step           = 120,
                           max_patience             = 24,
                           ## Output
                           log_step                 = 100,
                           print_step               = 100,
                           ## Seed
                           seed                     = 12,
                           ## VAE params
                           KL_weight                = 1.0,
                           heldout_percent          = 0.1,
                           num_latent_dims          = 128,     #
                           num_layers               = 2,
                           hidden_unit_2power       = 7,
                           decoder_var              = "none",
                           ## Batch size
                           batch_size_component     = 240,
                           batch_size_mixture       = 240,
                           ## Training method
                           training_method          = 'train',
                           infer_expr               = 'True',
                           dropout_rate             = 0.1,     #
                           tpm                      = "False",  #
                           deconvolution            = 'True',
                           batch_correction         = 'False',
                           decay_lr                 = 'False', #
                           batch_norm_layers        = 'False', #
                           corr_check               = 'True', # 'False', # 
                           log_results              = 'True',
                           log_samples              = 'False',
                           save_to_disk             = 'False',
                           save_to_object           = 'True',
                           cuda_device              = 0)
    ################################################################################
    file_pickle = open(path + rand_test_genes + '_deconv_impute_demo.obj', 'wb') 
    pickle.dump(deconvModel, file_pickle)
    file_pickle.close()
    fileObj = open(path + rand_test_genes + '_deconv_impute_demo.obj', 'rb')
    exampleObj = pickle.load(fileObj)
    fileObj.close()
    print(exampleObj)
    deconvModel = exampleObj
    deconvModel.deconvResults.deconv_data['purified']['train']["interneurons"].shape
    spatial_data_partial # 4530 × 26
    df_sum_purified = np.zeros(deconvModel.deconvResults.deconv_data['purified']['train']["interneurons"].shape)
    deconvModel.deconvResults.proportions[str(final)]
    deconvModel.deconvResults.deconv_data['purified']['train'].keys()
    cell_labels = list(deconvModel.deconvResults.deconv_data['purified']['train'].keys())
    prop = pd.DataFrame(deconvModel.deconvResults.proportions[str(final)], columns=cell_labels)
    for i in cell_labels:
        df_sum_purified = df_sum_purified + deconvModel.deconvResults.deconv_data['purified']['train'][i] * np.array(prop[i]).reshape(-1,1)
    df_sum_purified = pd.DataFrame(df_sum_purified, columns=list(seq_data.var.index)) 
    df_sum_purified_impute = df_sum_purified[rand_test_genes]
    spatial_data_test = spatial_data[:,rand_test_genes].copy()
    df_sum_purified_measured = pd.Series(spatial_data_test.X.reshape(-1,), name=rand_test_genes)
    fig=plt.figure(figsize=(12,7))
    plt.subplot(1,2,1)
    plt.scatter(x=spatial_data_test.obs.x_coord, y=spatial_data_test.obs.y_coord, c=df_sum_purified_measured, marker=".", cmap="viridis")
    #plt.axis('off')
    plt.subplot(1,2,2)
    plt.scatter(x=spatial_data_test.obs.x_coord, y=spatial_data_test.obs.y_coord, c=df_sum_purified_impute, marker=".", cmap="viridis")
    #plt.axis('off')
    fig.savefig(path + rand_test_genes + "rna_measured_imputed.jpg", bbox_inches='tight', dpi=300)
    corr = spearmanr(df_sum_purified_measured, df_sum_purified_impute)[0]
    print(corr)
    corr_list.append(corr)
    pd.Series(corr, name=rand_test_genes).to_csv(path + rand_test_genes + "rna_measured_imputed.csv")
    
    
    df_original[rand_test_genes] = df_sum_purified_measured
    df_imputed[rand_test_genes] = df_sum_purified_impute
    ```
