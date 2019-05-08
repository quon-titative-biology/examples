# Tutorial: Alignment of pancreatic islet cells sequenced using four different protocols.

This tutorial provides a guided alignment of pancreatic islet cells sequenced across four prococols: CEL-Seq (GSE81076), CEL-Seq2 (GSE85241), Fluidigm C1 (GSE86469), and Smart-Seq2 (E-MTAB-5061) as described by [Seurat](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html). In this tutorial we demonstrate the unsupervised alignment strategy of `scAlign` described in [Johansen et al, 2018](https://www.biorxiv.org/content/10.1101/504944v2) extended to alignment of multiple datasets (conditions) along with typical analysis utilizing the aligned dataset.

## Alignment goals
The following is a walkthrough of a multi-way alignment problem with `scAlign` and has been designed to provide an overview of data preprocessing, alignment and finally analysis on the alignment space. Here, our primary goals include:

1. Train scAlign on a large dataset by using an initial step of dimensionality reduction.
2. Learning a low-dimensional cell state space in which cells group by function and type, regardless of condition (sequencing protocol).

## Data setup
The gene count matrices used for this tutorial are hosted [here](https://www.dropbox.com/s/1zxbn92y5du9pu0/pancreas_v3_files.tar.gz?dl=1).

First, we perform a typical scRNA preprocessing step using the `Seurat` package. Then, reduce to the top 3,000 highly variable genes and perform an initial dimensionality reduction step using CCA improve convergence and reduce computation time.

```R
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)

## User paths
working.dir = "." #where the data files are located
results.dir = "." #where the output should be stored

## Load in data
pancreas.data = readRDS(file = paste0(working.dir, "pancreas_expression_matrix.rds"))
metadata      = readRDS(file = paste0(working.dir, "pancreas_metadata.rds"))

## Subset combined object into a list containing the four individual datasets
pancreas      = CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list = SplitObject(pancreas, attribute.1 = "tech")

## Independently preprocess each datasets
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]])
    pancreas.list[[i]] <- ScaleData(pancreas.list[[i]], do.scale=T, do.center=T, display.progress=T)
    pancreas.list[[i]] <- FindVariableGenes(pancreas.list[[i]], top.genes = 3000)
}

## Extract common set of genes across all datasets
genes.use = Reduce(union, lapply(pancreas.list, function(seurat.obj) seurat.obj@var.genes))

## Run CCA for input to scAlign
pancreas.multi.cca = RunMultiCCA(pancreas.list, genes.use, num.ccs = 10)
```

## scAlign setup
The general design of `scAlign` makes it agnostic to the input RNA-seq data representation. Thus, the input data can either be
gene-level counts, transformations of those gene level counts or a preliminary step of dimensionality reduction such
as canonical correlates or principal component scores. Here we convert the previously defined
`Seurat` objects to `SingleCellExperiment` objects in order to create the combined `scAlign` object.

```R
## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
pancreas.sce.list = lapply(pancreas.list,
                           function(seurat.obj){
                             SingleCellExperiment(assays = list(logcounts = seurat.obj@data[genes.use,],
                                                                scale.data = seurat.obj@scale.data[genes.use,]),
                                                  colData = seurat.obj@meta.data)
                            })
```
We now build the scAlign combined SCE object and append the results of RunMultiCCA. The output matrices of scAlign will always be ordered based on the sce.objects list order.

```R
## Create the combined scAlign object from list of SCE(s).
scAlignPancreas = scAlignCreateObject(sce.objects = pancreas.sce.list,
                                      project.name = "scAlign_Pancreatic_Islet")

## Append the results of MultiCCA to our scAlign object
reducedDim(scAlignPancreas, "MultiCCA") = pancreas.multi.cca@dr$cca@cell.embeddings
```

## Alignment of sequencing protocols
Now we align the young and old cpopulations for multiple input types which are specified by `encoder.data`. `scAlign` returns a
low-dimensional joint embedding space where the effect of age is removed allowing us to use the complete dataset for downstream analyses such as clustering or differential expression. For the gene level input we also run the decoder procedure which projects each cell into logcount space for both conditions to perform paired single cell differential expressional.

```R
## Run scAlign with CCA results as input to the encoder (alignment).
scAlignPancreas = scAlignMulti(scAlignPancreas,
                        options=scAlignOptions(steps=15000,
                                               log.every=5000,
                                               architecture="large",  ## 3 layer neural network
                                               num.dim=64),            ## Number of latent dimensions
                        encoder.data="MultiCCA",
                        supervised='none',
                        run.encoder=TRUE,
                        run.decoder=FALSE,
                        log.results=TRUE,
                        log.dir=file.path('./tmp'),
                        device="GPU")
```
