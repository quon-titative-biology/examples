# Tutorial: Supervised and Semi-supervised alignment of C57BL/6 hemoteopietic stem cells (HSCs) from young and old mice from Kowalcyzk et al.

In this tutorial we demonstrate the supervised and semi-supervised alignment strategy of `scAlign` described in [Johansen et al, 2018](https://www.biorxiv.org/content/10.1101/504944v2).

## Alignment goals
The following is a walkthrough of a typical alignment problem for `scAlign` and has been designed to provide an overview of data preprocessing, alignment and finally analysis in our joint embedding space. Here, our primary goals include:

1. Learning a low-dimensional cell state space in which cells group by function and type, regardless of condition (age).
2. Computing a single cell paired differential expression map from paired cell projections.

## Data setup
The gene count matrices used for this tutorial are hosted [here](https://github.com/quon-titative-biology/examples/blob/master/scAlign_kowalcyzk_et_al/kowalcyzk_gene_counts.rda).

First, we perform a typical scRNA preprocessing step using the `Seurat` package. Then, reduce to the top 3,000 highly variable genes from both datasets to improve convergence and reduce computation time.


```R
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)

## User paths
working.dir = "." #where our data file, kowalcyzk_gene_counts.rda is located
results.dir = "." #where the output should be stored

## Load in data
load(file.path(working.dir, 'kowalcyzk_gene_counts.rda'))

## Extract age and cell type labels
cell_age = unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 1))
cell_type = gsub('HSC', '', unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 2)))

## Separate young and old data
young_data = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "young")]
old_data   = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "old")]

## Set up young mouse Seurat object
youngMouseSeuratObj <- CreateSeuratObject(raw.data = young_data, project = "MOUSE_AGE", min.cells = 0)
youngMouseSeuratObj <- FilterCells(youngMouseSeuratObj, subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
youngMouseSeuratObj <- NormalizeData(youngMouseSeuratObj)
youngMouseSeuratObj <- ScaleData(youngMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)
youngMouseSeuratObj@meta.data$type = cell_type[which(cell_age == "young")]
youngMouseSeuratObj@meta.data$age  = "young"

## Set up old mouse Seurat object
oldMouseSeuratObj <- CreateSeuratObject(raw.data = old_data, project = "MOUSE_AGE", min.cells = 0)
oldMouseSeuratObj <- FilterCells(oldMouseSeuratObj, subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
oldMouseSeuratObj <- NormalizeData(oldMouseSeuratObj)
oldMouseSeuratObj <- ScaleData(oldMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)
oldMouseSeuratObj@meta.data$type = cell_type[which(cell_age == "old")]
oldMouseSeuratObj@meta.data$age  = "old"

## Gene selection
youngMouseSeuratObj <- FindVariableGenes(youngMouseSeuratObj, do.plot = F)
oldMouseSeuratObj <- FindVariableGenes(oldMouseSeuratObj, do.plot = F)
g.1 <- head(rownames(youngMouseSeuratObj@hvg.info), 3000)
g.2 <- head(rownames(oldMouseSeuratObj@hvg.info), 3000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(youngMouseSeuratObj@scale.data))
genes.use <- intersect(genes.use, rownames(oldMouseSeuratObj@scale.data))
```

## scAlign setup
The general design of `scAlign` makes it agnostic to the input RNA-seq data representation. Thus, the input data can either be
gene-level counts, transformations of those gene level counts or a preliminary step of dimensionality reduction such
as canonical correlates or principal component scores. Here we create the scAlign object from the previously defined
`Seurat` objects and perform both PCA and CCA on the unaligned data.

```R
## Create paired dataset SCE objects to pass into scAlignCreateObject
youngMouseSCE <- SingleCellExperiment(
    assays = list(counts = youngMouseSeuratObj@raw.data[genes.use,],
                  logcounts  = youngMouseSeuratObj@data[genes.use,],
                  scale.data = youngMouseSeuratObj@scale.data[genes.use,]),
    colData = youngMouseSeuratObj@meta.data
)

oldMouseSCE <- SingleCellExperiment(
  assays = list(counts = oldMouseSeuratObj@raw.data[genes.use,],
                logcounts  = oldMouseSeuratObj@data[genes.use,],
                scale.data = oldMouseSeuratObj@scale.data[genes.use,]),
  colData = oldMouseSeuratObj@meta.data
)

```

## Alignment of young and old HSCs using complete set type labels for both young and old mice.
We now build the scAlign SCE object and compute PCs and/or CCs using Seurat for the assay defined by `data.use`. It is assumed that `data.use`, which is being used for the initial step of dimensionality reduction, is properly normalized and scaled. Resulting combined matrices will always be ordered based on the sce.objects list order.

```R
scAlignHSC = scAlignCreateObject(sce.objects = list("YOUNG"=youngMouseSCE, "OLD"=oldMouseSCE),
                                 labels = list(youngMouseSCE@colData$type, oldMouseSCE@colData$type),
                                 data.use="scale.data",
                                 project.name = "scAlign_Kowalcyzk_HSC")

## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, early.stop=FALSE, architecture="small"),
                    encoder.data="scale.data",
                    supervised='both',
                    run.encoder=TRUE,
                    log.dir=file.path(results.dir, 'models','gene_input'),
                    device="CPU")
```

## Alignment of young and old HSCs using partial labels for the old (>2 years) mouse cells.
We now build the scAlign SCE object for semi-supervised alignment by remove some labels from the old mouse dataset.

```R
## Remove labels from two cell types for the old mouse cells
oldMouseSCE@colData$type[which(oldMouseSCE@colData$type == "LT")] = NA
oldMouseSCE@colData$type[which(oldMouseSCE@colData$type == "ST")] = NA

scAlignHSC = scAlignCreateObject(sce.objects = list("YOUNG"=youngMouseSCE, "OLD"=oldMouseSCE),
                                 labels = list(youngMouseSCE@colData$type, oldMouseSCE@colData$type),
                                 data.use="scale.data",
                                 project.name = "scAlign_Kowalcyzk_HSC")

## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlignHSC = scAlign(scAlignHSC,
                    options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, early.stop=FALSE, architecture="small"),
                    encoder.data="scale.data",
                    supervised='both',
                    run.encoder=TRUE,
                    log.dir=file.path(results.dir, 'models','gene_input'),
                    device="CPU")
```
