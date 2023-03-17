# Tutorial: Deconvolution of CellBench mixtures

This tutorial provides a guided deconvolution of the benchmark dataset CellBench which contains simulated bulk RNA-seq. Deconvolution requires a reference atlas to define the gene expression patterns of individual cell types within the mixtures.

## Steps for tutorial
The following is a walkthrough of a standard deconvolution problem with `scProjection` and has been designed to provide an overview of setting up the deconvolution model and sanity checking the results. Here, our primary goals include:

1. Setup an environment for scProjection
2. Define and train `scProjection` on a reference atlas and performing deconvolution bulk RNA-seq.
3. Perform standard sanity checks on estimated cell type abundances.

## Setup python/system
Launch a virtual environment that contains the dependencies for scProjection, as detailed here: https://github.com/quon-titative-biology/scProjection#setup-of-virtualenv
```shell
source ./pyvTf2/bin/activate
```

## Setup R workspace
Now we launch R where the remaining steps will be performed. Please download the R scripts located here: https://github.com/quon-titative-biology/examples/tree/master/scProjection_cellbench/R

```R
library(reticulate) ## Let's us call Python from R!

##
scProjection = import("scProjection")

## Load in helper functions/scripts
source('nnDeconvClass.R')
source('proportionHeatmap.R')

```

## Load data
In this tutorial we will be using a matched scRNA-seq atlas also from CellBench, the data can be found here: https://github.com/quon-titative-biology/examples/blob/master/scProjection_cellbench/data/deconv_tutorial_CELBENCH.rda.
```R
load("deconv_tutorial_CELBENCH.rda")
```

## Deconv. model setup
Now that we have our cell type labels and marker gene set finalized it's time to build the deconvolution model.
`scProjection` requires 3 inputs:
```text
scRNA Atlas                   cell x gene
scRNA Atlas cell type labels  cell x 1
mixture                       sample x gene
```
```R
## Load scProjection with the component and mixture rna-seq data
deconvModel = scProjection$deconvModel(component_data   = as.matrix(component_data),
                                       component_label  = as.array(as.character(component_label)),
                                       mixture_data     = as.matrix(mixture_data))


```

Furthermore, we have to determine the marker genes to focus on for deconvolution.
```R
## Lets define a mask for the marker genes (1=marker, 0=non_marker). In the case, all the genes are marker genes!
marker_gene_mask = rep(0, ncol(component_data))
marker_gene_mask[which(colnames(component_data) %in% marker.genes)] = 1
```

## Deconvolution of CellBench
Now, we take our defined `deconvModel` and run the deconvolution task!

```R
## Now lets run deconvolution, I have left the hyperparameters exposed here so you can see the different ways that the model could be tuned.
deconvModel$deconvolve(## Masks
                       marker_gene_mask = marker_gene_mask,
                       ## Training steps
                       max_steps_component  = 2500L,
                       max_steps_proportion = 500L,
                       max_steps_combined   = 250L,
                       ## Learning rate
                       component_learning_rate  = 1e-3,
                       proportion_learning_rate = 1e-3,
                       combined_learning_rate   = 1e-3,
                       ## Early stopping
                       early_stopping       = 'True',
                       early_min_step       = 500L,
                       max_patience         = 100L,
                       ## Output
                       log_step             = 25L,
                       print_step           = 100L,
                       ## Seed
                       seed                 = as.integer(1234),
                       ## VAE params
                       KL_weight            = 1.0,
                       num_latent_dims      = as.integer(3),
                       num_layers           = as.integer(32),
                       hidden_unit_2power   = as.integer(9),
                       decoder_var          = "per_sample",
                       ## Batch size
                       batch_size_component = 100L,
                       batch_size_mixture   = 1000L,
                       ## GPU info
                       cuda_device         = 0L)
```

## Getting the deconv results
While `deconvModel` holds all the input data and results, it is a python class which R can't interpret. So we have a helper function to convert `deconvModel` to an R S4 class which will contain all the results from deconvolution plus data annotations.
```R
## Convert the python class to an R S4 class. check: str(deconvResults)
deconvResults = convertDeconv(deconvModel,               ## Deconv model object
                              rownames(component_data),  ## Features used during deconvolution
                              colnames(component_data),  ## scRNA cell ids
                              colnames(mixture_data))    ## mixture sample ids

## Get the final proportion estimates
proportions = deconvResults@proportions$final

## Call cell type from proportions
predicted.labels = deconvModel$celltypes[apply(proportions, 1, which.max)]
```

## Plot the estimated proportions
```R
png(paste0("cellbench_deconv_summary.png"), width=14, height=12, units="in", res=600)
propHeatmap(proportions, mixture.labels=ground.truth.label)
dev.off()
```

![]("https://github.com/quon-titative-biology/examples/blob/master/scProjection_cellbench/figures/cellbench_deconv_summary.png")

<!-- ## Proportion sanity check
Let's take a quick look at the proportions and make sure the marker gene expression supports the estimated cell type proportions.
```R
proportionSanityCheck(deconvResults,
                      mixture_data,
                      data = "scale.data", ## "data" or "scale.data" from Seurat object
                      filename="./mouse_patchseq_deconv_sanity.png",
                      marker.anno=marker.anno.df)
```
![proportion_check](https://github.com/ucdavis/quonlab/blob/master/development/deconvAllen/deconvTutorial/mouse_patchseq_example_run.png) -->

## Projection Tutorial
Finally, `scProjection` also outputs all the mixture samples purified to each cell type in the scRNAseq atlas. Let's take a look at the purified and scRNA data.
```R
## Let's first check all the types of purified data we have!
print(names(deconvResults@deconv_data$purified$train))

## We extract the purified data for the `Sst` cell type as follows:
pure.data.sst = t(deconvResults@deconv_data$purified$train$Sst) ## Make gene x cell for now
rownames(pure.data.sst) = deconvResults@genes
colnames(pure.data.sst) = deconvResults@mixtureCellNames

## Now lets plot the scRNA and purified data together to see how our model performed.
## First step, convert the purified data to a Seurat object (in a somewhat hacky way)
pure.seurat = CreateSeuratObject(pure.data.sst)
pure.seurat@assays$RNA@data = pure.seurat@assays$RNA@scale.data = pure.data.sst

## Merge the pure and scRNA Data
merged.seurat = merge(mouse, pure.seurat)
merged.seurat@assays$RNA@scale.data = cbind(GetAssayData(mouse, "scale.data")[rownames(pure.seurat),],
                                            GetAssayData(pure.seurat, "scale.data")[rownames(pure.seurat),])

## Update the celltype.deconv for the Sst purified data
merged.seurat@meta.data$celltype.deconv = c(mouse@meta.data$celltype.deconv, rep("Sst", ncol(pure.seurat)))
merged.seurat@meta.data$datatype = c(rep("scRNA", ncol(mouse)), rep("purifiedSst", ncol(pure.seurat)))

## Run dimensionality reduction on just the marker genes
merged.seurat <- RunPCA(merged.seurat, npcs = 30, features=marker.genes, verbose = FALSE)
merged.seurat <- RunUMAP(merged.seurat, reduction = "pca", dims = 1:30)

## Plot
celltype.plot <- DimPlot(merged.seurat, reduction = "umap", group.by = "celltype.deconv", label=TRUE, repel=TRUE)
datatype.plot <- DimPlot(merged.seurat, reduction = "umap", group.by = "datatype", label=TRUE, repel=TRUE)

## Save our plot
library(cowplot)
png("./mouse_pure_scRNA_atlas.png", width=24, height=12, res=150,  units="in")
plot_grid(celltype.plot, datatype.plot)
dev.off()
```
![purified_check]()

```R
sessionInfo()
```
