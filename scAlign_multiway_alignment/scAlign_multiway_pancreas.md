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

## PMA depends on 'impute' from bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("impute")
library(PMA)

## User paths
working.dir = "." #where the data files are located
results.dir = "." #where the output should be stored

## Load in data
pancreas.data = readRDS(file = file.path(working.dir, "pancreas_expression_matrix.rds"))
metadata      = readRDS(file = file.path(working.dir, "pancreas_metadata.rds"))

## Subset combined object into a list containing the four individual datasets
pancreas      = CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas.list = SplitObject(pancreas, split.by = "tech")

## Independently preprocess each datasets
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]])
    pancreas.list[[i]] <- ScaleData(pancreas.list[[i]], do.scale=T, do.center=T, display.progress=T)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], nFeatures = 3000)
}

## Extract common set of genes across all datasets
genes.use = Reduce(intersect, lapply(pancreas.list, function(seurat.obj) VariableFeatures(seurat.obj)))
```

## scAlign setup
The general design of `scAlign` makes it agnostic to the input RNA-seq data representation. Thus, the input data can either be
gene-level counts, transformations of those gene level counts or a preliminary step of dimensionality reduction such
as canonical correlates or principal component scores. Here we convert the previously defined
`Seurat` objects to `SingleCellExperiment` objects in order to create the combined `scAlign` object to which we append the results of Seurat's `RunMultiCCA`.

```R
## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
pancreas.sce.list = lapply(pancreas.list,
                           function(seurat.obj){
                             SingleCellExperiment(assays = list(counts = seurat.obj@assays$RNA@counts[genes.use,],
                                                                logcounts = seurat.obj@assays$RNA@data[genes.use,],
                                                                scale.data = seurat.obj@assays$RNA@scale.data[genes.use,]),
                                                  colData = seurat.obj@meta.data)
                            })

## Create the combined scAlign object from list of SCE(s).
scAlignPancreas = scAlignCreateObject(sce.objects = pancreas.sce.list,
                                      genes.use = genes.use,
                                      cca.reduce=T,
                                      ccs.compute=10,
                                      project.name = "scAlign_Pancreatic_Islet")
```

## Alignment of sequencing protocols
Now we align the pancreas islets sequenced across different protocols using the results of `RunMultiCCA` as input to `scAlign`.

```R
## Run scAlign with CCA results as input to the encoder (alignment).
scAlignPancreas = scAlignMulti(scAlignPancreas,
                        options=scAlignOptions(steps=15000,
                                               log.every=5000,
                                               batch.size=300,
                                               perplexity=30,
                                               norm=TRUE,
                                               batch.norm.layer=FALSE,
                                               architecture="large",  ## 3 layer neural network
                                               num.dim=64),            ## Number of latent dimensions
                        encoder.data="MultiCCA",
                        decoder.data="scale.data",
                        supervised='none',
                        run.encoder=TRUE,
                        run.decoder=TRUE,
                        log.results=TRUE,
                        log.dir=file.path('./tmp'),
                        device="GPU")
```
After alignment, `scAlign` returns a
low-dimensional joint embedding space where the effect of sequencing protocol is removed allowing us to use the complete dataset for downstream analyses such as clustering
or differential expression.
![Results](https://github.com/quon-titative-biology/examples/blob/master/scAlign_multiway_alignment/figures/pancreas_result.png)

## Differential expression using scAlign embeddings for clustering
After alignment we can use the resulting embeddings to cluster the cells without dataset specific biases then perform differential expression analysis. Further
information about using Seurat to find clusters and differentially expressed genes can be found [here](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html ).

```R
library(dplyr)
scAlignSeuratObj = as.Seurat(scAlignPancreas, counts="counts", scale.data="scale.data")

## Cluster the data based on scAlign embeddings, be sure to use all the embedding dimensions!
scAlignSeuratObj <- FindNeighbors(scAlignSeuratObj, dims = 1:64, reduction="ALIGNED-MultiCCA")
scAlignSeuratObj <- FindClusters(scAlignSeuratObj, resolution = 0.4)
table(Idents(scAlignSeuratObj))

## Find markers for each scAlign based cluster compared against the remaining cells.
pancreas.markers <- FindAllMarkers(scAlignSeuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## Alternatively, find markers which are conserved when performing DE on each dataset individually.
pancreas.markers <- FindConservedMarkers(scAlignSeuratObj, ident.1 = 0, ident.2 = 1, grouping.var = "group.by")
```R
