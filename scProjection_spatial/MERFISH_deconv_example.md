## Load pkgs

```
library(reticulate)
library(caret)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(gplots)
library(RColorBrewer)
library(caret)
library(class)
library(data.table)
scProjection = import("scProjection")
options(stringsAsFactors=F)
source('./nnDeconvClass.R')
source('./proportionHeatmap.R')
```
## load data
```
MERFISH <- readRDS("MERFISH_deconv.rds") # MERFISH data used for deconv task
scRNA.seurat <- readRDS("scRNA.rds")     # single-cell dataset from the same study
```
## determine feature space
```
hvgs.sc <- VariableFeatures(scRNA.seurat) # 2000
markers.mer <- rownames(MERFISH)          # 154
common.genes <- intersect(hvgs.sc, markers.mer)
```
## set current cell class level
```
Idents(scRNA.seurat) <- scRNA.seurat@meta.data$mapped_celltype
Idents(MERFISH)      <- MERFISH@meta.data$mapped_celltype
# ground truth
merfish.label <- MERFISH@meta.data$mapped_celltype
merfish.bregma <- factor(MERFISH@meta.data$Bregma)
```
## prepate for deconv task
```
component_data  <- scRNA.seurat@assays$RNA@scale.data[common.genes,]
component_label <- scRNA.seurat@meta.data$mapped_celltype
mixture_data    <- MERFISH@assays$RNA@scale.data[common.genes,]
dim(component_data) # 115 13904
dim(mixture_data)   # 115 33036
deconvModel = NULL
deconvModel = scProjection$deconvModel(component_data   = as.matrix(t(component_data)),
                                       component_label  = as.array(as.character(component_label)),
                                       mixture_data     = as.matrix(t(mixture_data))) 
```
## Now we define the network architecture
```
deconvModel$deconvolve(## Masks
                       #marker_gene_mask = marker_gene_mask,
                       ## Training steps
                       max_steps_component  = as.integer(1000),
                       max_steps_proportion = 10000L, ## 30 times over whole dataset of 1M MEFISH
                       max_steps_combined   = 100L,   ## was 100L
                       ## Learning rate
                       component_learning_rate  = 1e-3,
                       proportion_learning_rate = 1e-2, # was 1e-3
                       combined_learning_rate   = 1e-4,
                       ## Early stopping
                       early_stopping       = 'True',
                       early_min_step       = 500L,
                       max_patience         = 50L,
                       ## Output
                       log_step             = 100L,
                       print_step           = 100L,
                       ## Seed
                       seed                 = as.integer(1234),
                       ## VAE params
                       KL_weight            = 1.0,
                       num_latent_dims      = as.integer(32),
                       num_layers           = as.integer(3),
                       hidden_unit_2power   = as.integer(9),
                       decoder_var          = "per_sample",
                       ## Batch size
                       batch_size_component = 150L,
                       batch_size_mixture   = 1000L,
                       ## Training method
                       training_method     = 'train', ## train: use all scRNA data, "valid": split scRNA into training and testing sets
                       #infer_expr          = 'True',
                       ## Deconvolution flags
                       early_stopping      = 'True',
                       deconvolution       = 'True',
                       batch_correction    = 'False',
                       decay_lr            = 'False',
                       batch_norm_layers   = 'True',
                       corr_check          = 'True',
                       log_results         = 'True',
                       log_samples         = 'False',
                       save_to_disk        = 'False',
                       save_to_object      = 'True',
                       logdir              = paste0("/home/qlab/scProjection/models_ckpt/_long_training_",2500,"_",32, "_layer_",3,"_",9),
                       cuda_device         = 0L)
```
## get the results
```
deconvResults = convertDeconv(deconvModel,               ## Deconv model object
                              rownames(component_data),  ## Features used during deconvolution
                              colnames(component_data),  ## scRNA cell ids
                              colnames(mixture_data))    ## mixture sample ids
```

## Get the final proportion estimates
```
proportions = deconvResults@proportions$final
dim(proportions) #  33036     8
colnames(proportions) = deconvResults@celltypes
head(proportions)
MERFISH@meta.data$mapped_celltype
# Call cell type from proportions
predicted.labels = deconvModel$celltypes[apply(proportions, 1, which.max)]
predicted.label = deconvResults@celltypes[apply(deconvResults@proportions$final,1,which.max)]

prop = as.data.frame(cbind(proportions, MERFISH@meta.data$mapped_celltype))
prop_ = as.data.frame(predicted.label, MERFISH@meta.data$mapped_celltype)

fwrite(prop_, file = "prop_labels_match_MERFISH_common.csv", quote = F, row.names=T, col.names=T)
fwrite(prop, file = "prop_match_MERFISH_common.csv", quote = F, row.names=T, col.names=T)

saveRDS(deconvResults, "deconvResults_match_MERFISH_common.rds")
#saveRDS(deconvModel, "deconvModel_match_MERFISH_common.rds")
saveRDS(rownames(component_data), "rownames_component_data_match_MERFISH_common.rds")
saveRDS(colnames(component_data), "colnames_component_data_match_MERFISH_common.rds")
saveRDS(colnames(mixture_data), "colnames_mixture_data_match_MERFISH_common.rds")
saveRDS(proportions, "proportions_match_MERFISH_common.rds")
```
